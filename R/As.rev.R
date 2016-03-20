
library(rgdal)
library(GSIF)
library(gdalUtils)
library(raster)
library(dplyr)
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

gk_7 <- "+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m"
utm <- "+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m"

load(paste(getwd(),"inst","extdata","gridmaps.RDa",sep="/"))
gridmaps.sm2D$CD <- exp(-gridmaps.sm2D$DirDif)
gridmaps.sm2D$DD <- as.numeric(exp(-scale(gridmaps.sm2D$Dist,center=FALSE)))

bor <- join(read.csv(paste(getwd(),"inst","extdata","Profili_sredjeno_csv.csv",sep="/")), read.csv(paste(getwd(),"inst","extdata","Koordinate_csv.csv",sep="/")), type="inner")
bor$hdepth<-bor$Bottom-bor$Top
bor$altitude <- - (bor$Top / 100 + (( bor$Bottom - bor$Top ) / 2) / 100)
bor <- bor[, - c(7:12,14,15,16,17 )]
names(bor) <- c("Soil.Type","ID","Horizont","Top" , "Bottom","pH","Humus","As","Co","x","y","hdepth","altitude")

bor <- bor[-which(bor$ID %in% c(66,129,14,51,69,130,164,165,166)),]

## Predict soil properties of interest:
bor.profs <- bor[,c("ID","x","y","Soil.Type","Top","Bottom","Humus","pH","Co","As")]
depths(bor.profs) <- ID ~ Top + Bottom
site(bor.profs) <- ~ Soil.Type + x + y
coordinates(bor.profs) <- ~ x+y
proj4string(bor.profs) <- CRS(utm)
bor.geo<-as.geosamples(bor.profs)

bor.xy <- bor[complete.cases(bor[,c("ID","x","y","altitude","As")]),c("ID","x","y","hdepth","altitude","As")] 
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
bor.df <- bor.df[complete.cases(bor.df[,names(bor.df)]),] # Uvek razmisiti kod ove funkcije
str(bor.df)

#============== Removing correlated covariates (previously determined) ==============
names(bor.df)
bor.df<-bor.df[,-which(names(bor.df) %in% c("DirDif","Dist","AnalyticalHills","LSFactor","RelSlopePosition","VelleyDepth","optional"))]
str(bor.df)



#============== Sorting by altitude =================================================
bor.df<-ddply(bor.df,.(altitude))
bor.df <- bor.df[order(-bor.df$altitude),] 
head(bor.df)
names(bor.df)

sm2D.lst<-names(bor.df)[c(5:23)]


#====================== formulas ===================================================
fun <- as.formula(paste("As ~", paste(c(sm2D.lst,"altitude"), collapse="+")))
fm.int.H<- as.formula(paste("As ~",paste(sm2D.lst,"altitude",sep="*", collapse="+"),sep=""))
fm.GSIF.H <- as.formula(paste("As ~", paste(c(sm2D.lst,"altitude","ns(altitude,df=4)"), collapse="+")))
fm.GSIF.int.H <- as.formula(paste("As ~",paste(paste(sm2D.lst,"altitude",sep="*" ,collapse="+"),"ns(altitude,df=4)",sep="+")))
fm.GSIF.int.H1 <- as.formula(paste("As ~",paste(paste(sm2D.lst,"poly(altitude,3)",sep="*",collapse="+"))))

fm.poly.H <- as.formula(paste("As ~", paste(c(sm2D.lst,"altitude","poly(altitude,3)"), collapse="+")))

#===============================================================================================

#================== test for stratfold3d and penint3D ============================================
source(paste(getwd(),"R","stratFold3D.R",sep="/"))
source(paste(getwd(),"R","penint3D.R",sep="/"))
source(paste(getwd(),"R","plotfolds.R",sep="/"))


fun<-as.formula(paste("As ~", paste(c(sm2D.lst,"altitude"), collapse="+")))
regdat<-bor.df
contVar<-sm2D.lst[-which(sm2D.lst %in% c("clc","SoilType"))]
targetVar<-"As"

#regdat$altitude<-2+regdat$altitude
#regdat[1:20,1:4]

#set.seed(432)
detach("package:dplyr", unload=TRUE)
strat<-stratfold3d(targetVar="As",seed=123,regdat=regdat,folds=5,cent=3,dimensions="2D",IDs=TRUE,sum=TRUE)
flist<-strat$folds
#plotfolds(strat,"pH")
library(dplyr)

rez<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="linear")
rez
summary(rez$measure[1:5,])

rezint<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="linear")
rezint
summary(rezint$measure[1:5,])

rezinthier<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = TRUE, int=FALSE, depth.fun="linear")
rezinthier
summary(rezinthier$measure[1:5,])

#============= Poly ========================
rez.poly<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="poly")
rez.poly
summary(rez.poly$measure[1:5,])

rezint.poly<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
rezint.poly
summary(rezint.poly$measure[1:5,])

rezinthier.poly<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = TRUE, int=FALSE, depth.fun="poly")
rezinthier.poly
summary(rezinthier.poly$measure[1:5,])


#================= Prediction Final Model ========================================

rez.pred.l<-penint3Dpred(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="linear")
rez.pred.l$cv.par
rez.pred.l$results

rezint.pred.l<-penint3Dpred(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="linear")
rezint.pred.l$cv.par
rez.l.pred.l$results

rezinthier.pred.l<-penint3Dpred(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="linear")
rezinthier.pred.l$cv.par
rezinthier.pred.l$results

rez.pred.l<-penint3Dpred(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
rez.pred.l$cv.par
rez.pred.l$results

rezint.pred.l<-penint3Dpred(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
rezint.pred.l$cv.par
rez.l.pred.l$results

rezinthier.pred.p<-penint3Dpred(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
rezinthier.pred.p$cv.par
rezinthier.pred.p$results

penint3Dpred












stdepths <- c(-.1,-.3,-.5)
new3D <- sp3D(gridmaps.sm2D, stdepths=stdepths)
str(new3D)


gridmaps.sm2D$As<-rep(1,dim(gridmaps.sm2D@data)[1])
gridmaps.sm2D$altitude<-rep(-0.05,dim(gridmaps.sm2D@data)[1])


modmat <- model.matrix(fun ,gridmaps.sm2D@data)[,-1]
# removing nzv varaibles 
nzv <- nearZeroVar(modmat)
if(sum(nzv)!=0){modmat <- modmat[, -nzv]}else{modmat<-modmat}











