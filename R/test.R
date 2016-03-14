
library(rgdal)
library(GSIF)
library(gdalUtils)
library(raster)
library(plyr)
library(aqp)
library(psych)
library(mda)
library(classInt)
library(caret)
library(corrplot)
library(MASS)
library(splines)
library(glmnet)
library(glinternet)
library(hierNet)
library(magrittr)
library(dplyr)

gk_7 <- "+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m"
utm <- "+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m"

load(paste(getwd(),"inst","extdata","gridmaps.RDa",sep="/"))
grids$CD<-exp(-grids$DirDif)
grids$DD<-as.numeric(exp(-scale(grids$Dist,center=FALSE)))

bor <- join(read.csv(paste(getwd(),"inst","extdata","Profili_sredjeno_csv.csv",sep="/")), read.csv(paste(getwd(),"inst","extdata","Koordinate_csv.csv",sep="/")), type="inner")
bor$hdepth<-bor$Bottom-bor$Top
bor$altitude <- - (bor$Top / 100 + (( bor$Bottom - bor$Top ) / 2) / 100)
bor <- bor[, - c( 7:12,14,15,16,17 )]
names(bor) <- c("Soil.Type","ID","Horizont","Top" , "Bottom","pH","Humus","As","Co","x","y","hdepth","altitude")

bor.xy <- bor[complete.cases(bor[,c("ID","x","y","altitude","Humus")]),c("ID","x","y","hdepth","altitude","Humus")] 
bor.xy <- plyr::rename(bor.xy, replace=c("x" = "longitude", "y" = "latitude"))
coordinates(bor.xy) <- ~ longitude + latitude
proj4string(bor.xy) <- CRS(utm)
bor.xy <- spTransform(bor.xy, CRS(gk_7))

sm2D.lst<-names(gridmaps.sm2D)

#=========== Overlay with covariates (gridmaps.sm2D) ================================
ov.bor <- over(bor.xy, gridmaps.sm2D)
bor.xy@data[,c(sm2D.lst)] <- ov.bor
str(bor.xy)

bor.df <- data.frame(bor.xy)
bor.df <- bor.df[complete.cases(bor.df[,names(bor.df)]),]
str(bor.df)

#============== Removing correlated covariates (previously determined) ==============
names(bor.df)
bor.df<-bor.df[,-which(names(bor.df) %in% c("DirDif","Dist","ES","AnalyticalHills","LSFactor","ChannelNetBaseLevel","RelSlopePosition","CatchArea","WindEffectNorthWest","VelleyDepth","optional"))]
str(bor.df)

#bor.Humus.df<-bor.Humus.df[-which(bor.Humus.df$ID %in% c(141)),] #153,137,141,145,150

#============== Sorting by altitude =================================================
bor.df<-ddply(bor.df,.(altitude))
bor.df <- bor.df[order(-bor.df$altitude),] 
head(bor.df)
names(bor.df)

sm2D.lst<-names(bor.df)[c(5:17)]

#====================== formulas ===================================================
fm.H <- as.formula(paste("Co ~", paste(c(sm2D.lst,"altitude"), collapse="+")))
fm.int.H<- as.formula(paste("Humus ~",paste(sm2D.lst,"altitude",sep="*", collapse="+"),sep=""))
fm.GSIF.H <- as.formula(paste("Humus ~", paste(c(sm2D.lst,"altitude","ns(altitude,df=4)"), collapse="+")))
fm.GSIF.int.H <- as.formula(paste("Humus ~",paste(paste(sm2D.lst,"altitude",sep="*" ,collapse="+"),"ns(altitude,df=4)",sep="+")))
fm.GSIF.int.H1 <- as.formula(paste("Humus ~",paste(paste(sm2D.lst,"poly(altitude,3)",sep="*",collapse="+"))))

fm.poly.H <- as.formula(paste("Humus ~", paste(c(sm2D.lst,"altitude","poly(altitude,3)"), collapse="+")))

#===============================================================================================

#================== test for stratfold3d and penint3D ============================================
source(paste(getwd(),"R","stratFold3D.R",sep="/"))
source(paste(getwd(),"R","penint3D.R",sep="/"))
source(paste(getwd(),"R","plotfolds.R",sep="/"))


fun<-as.formula(paste("Humus ~", paste(c(sm2D.lst,"altitude"), collapse="+")))
regdat<-bor.df#[-which(bor.df$ID %in% c(66,129,14,51,69,130,164,165,166)),]
contVar<-sm2D.lst[1:11]
targetVar<-"Humus"

#regdat$altitude<-2+regdat$altitude
#regdat[1:20,1:4]

#set.seed(432)
detach("package:dplyr", unload=TRUE)
strat<-stratfold3d(targetVar="Humus",seed=123,regdat=regdat,folds=10,cent=3,dimensions="2D",IDs=FALSE,sum=TRUE)
flist<-strat$folds
#plotfolds(strat,"pH")
library(dplyr)

rez<-penint3D(fun=fun,regdat=regdat,hier = FALSE,lambda=10^seq(-5,5,length=110),contVar=contVar,int=FALSE,flist=flist,depth.fun="nspline")
rez$measure
summary(rez$measure[1:10,])

rezint<-penint3D(fun=fun,regdat=regdat,hier = FALSE,lambda=10^seq(-5,5,length=110),contVar=contVar,int=TRUE,flist=flist,depth.fun="nspline")
rezint$measure
summary(rezint$measure[1:10,])

rezinthier<-penint3D(fun=fun,regdat=regdat,hier=TRUE,lambda=10^seq(-5,5,length=110),contVar=contVar,int=TRUE,flist=flist,depth.fun="nspline")
rezinthier$measure
summary(rezinthier$measure[1:10,])




# Verovatno treba izbaciti dubinu iz skaliranja
# i sta sa onda sa dubinom...kako ona utice onda...
#











#================= edgeroi ===================================================================

data(edgeroi)

edgeroi$sites[edgeroi$sites$SOURCEID=="399_EDGEROI_ed095_1",]
edgeroi$horizons[edgeroi$horizons$SOURCEID=="399_EDGEROI_ed095_1",]


sites <- edgeroi$sites
coordinates(sites) <- ~ LONGDA94 + LATGDA94
proj4string(sites) <- CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs")
sites <- spTransform(sites, CRS("+init=epsg:28355"))

pnts <- list("sp.points", sites, pch="+", col="black")

## load the 250 m grids:
con <- url("http://gsif.isric.org/lib/exe/fetch.php?media=edgeroi.grids.rda")
load(con)
str(edgeroi.grids)
gridded(edgeroi.grids) <- ~x+y
proj4string(edgeroi.grids) <- CRS("+init=epsg:28355")
spplot(edgeroi.grids[1], sp.layout=pnts)


## load the 100 m grids:
con2 <- url("http://gsif.isric.org/lib/exe/fetch.php?media=edgeroi.grids100.rda")
load(con2)
str(edgeroi.grids100)
gridded(edgeroi.grids100) <- ~x+y
proj4string(edgeroi.grids100) <- CRS("+init=epsg:28355")
spplot(edgeroi.grids100["TI1LAN6"], sp.layout=pnts)
