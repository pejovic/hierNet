
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
library(doParallel)
library(foreach)

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
source(paste(getwd(),"R","penint3D_def.R",sep="/"))
source(paste(getwd(),"R","plotfolds.R",sep="/"))
source(paste(getwd(),"R","penint3Dpred.R",sep="/"))


fun<-as.formula(paste("Humus ~", paste(c(head(sm2D.lst,(length(sm2D.lst)-2)),"altitude"), collapse="+")))
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
summary(Humus.result.lst$rez$measure[1:5,])

rezint<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="linear")
rezint
summary(rezint$measure[1:5,])
summary(Humus.result.lst$rezint$measure[1:5,])

rezinthier<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = TRUE, int=FALSE, depth.fun="linear")
rezinthier
summary(rezinthier$measure[1:5,])
summary(Humus.result.lst$rezinthier$measure[1:5,])

#============= Poly ========================
rez.poly<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="poly")
rez.poly
summary(rez.poly$measure[1:5,])
summary(Humus.result.lst$rez.poly$measure[1:5,])

rezint.poly<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
rezint.poly
summary(rezint.poly$measure[1:5,])
summary(Humus.result.lst$rezint.poly$measure[1:5,])

rezinthier.poly<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="poly")
rezinthier.poly
summary(rezinthier.poly$measure[1:5,])


#================= Prediction Final Model ========================================

rez.pred.l<-penint3Dpred(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="linear")
rez.pred.l$cv.par
rez.pred.l$results


rezint.pred.l<-penint3Dpred(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="linear")
rezint.pred.l$cv.par
rezint.pred.l$results
Humus.result.lst$rezint.pred.l$results

rezinthier.pred.l<-penint3Dpred(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="linear")
rezinthier.pred.l$cv.par
rezinthier.pred.l$results

rez.pred.p<-penint3Dpred(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="poly")
rez.pred.p$cv.par
rez.pred.p$results

rezint.pred.p<-penint3Dpred(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
rezint.pred.p$cv.par
rezint.pred.p$results

rezinthier.pred.p<-penint3Dpred(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="poly")
rezinthier.pred.p$cv.par
rezinthier.pred.p$results

Humus.result.lst<-list(rez=rez,rezint=rezint,rezinthier=rezinthier, rez.poly=rez.poly, rezint.poly=rezint.poly, rezinthier.poly=rezinthier.poly, rez.pred.l=rez.pred.l,rez.pred.p=rez.pred.p,rezint.pred.l=rezint.pred.l,rezint.pred.p=rezint.pred.p,rezinthier.pred.l=rezinthier.pred.l,rezinthier.pred.p=rezinthier.pred.p)

str(Humus.result.lst)
names(Humus.result.lst) # fali rezinthier.pred.p jer je javljao gresku kod kreiranja "results" dimenzije nisu dobre...

save(Humus.result.lst,file="HumusResults.Rda")

load("HumusResults.Rda")



#=================== predint3D ====================================================

fun<-as.formula(paste("As ~", paste(c(head(sm2D.lst,(length(sm2D.lst))),"altitude"), collapse="+")))
cogrids <- gridmaps.sm2D
profs<-bor.profs
int=TRUE; hier=TRUE; depth.fun="linear";pred=TRUE;depths=c(-.1,-.3); lambda=seq(0,5,0.1); deg=3;depth.fun="poly"; preProc=TRUE; cent=3; fold=5; seed=321;cores=8;chunk=20000;l=c(370000:450000)
source(paste(getwd(),"R","stratFold3D.R",sep="/"))


#predint3D<-function(fun, profs, cogrids, hier=FALSE,pred=TRUE,lambda=seq(0,5,0.1),deg=3,fold=5,cent=3,int=TRUE,depth.fun=list("linear","poly"),depths=c(-.1,-.3),chunk=20000,preProc=TRUE,cores=2,seed=321,l=c(1:621426)){
source(paste(getwd(),"R","predint3D.R",sep="/"))
predint3D.As.g<-predint3D(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=TRUE ,hier = TRUE, int=TRUE, depth.fun="linear",cores=8)

str(predint3D.As.g)
names(predint3D.As.g)

summary(predint3D.As.g$prediction[[1]]$pred)

(predint3D.As.g$summary$results)

pp<-raster(predint3D.As.g[[2]],"pred")

plot(pp)

spplot(pp,"pred")



#================== Visualization ===================================
p = get("cellsize", envir = GSIF.opts)[2]
s = c(-.1,-.3)#get("stdepths", envir = GSIF.opts)

sd.ll <- sapply(1:length(predint3D.pred), FUN=function(x){make.3Dgrid(predint3D.pred[[x]][c("pred")], pixelsize=p, stdepths=s[x])})

sd.bb<-raster(sd.ll[[1]])
plot(sd.bb)

spplot(sd.ll[[1]],"pred")

z0 = mean(gridmaps.sm2D$DEM, na.rm=TRUE)

for(j in 1:length(sd.ll)){
  kml(slot(SNDMHT.gsm, paste("sd", j, sep="")),
              folder.name = paste("eberg_sd", j, sep=""),
              file = paste("SNDMHT_sd", j, ".kml", sep=""),
              colour = M, z.lim=c(10,85),
              raster_name = paste("SNDMHT_sd", j, ".png", sep=""),
              altitude = z0+5000+(s[j]*2500))
   }


