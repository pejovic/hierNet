
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
fun <- as.formula(paste("pH ~", paste(c(sm2D.lst[-which(sm2D.lst %in% c("ES","CD","DD"))],"altitude"), collapse="+")))
#===============================================================================================

#================== test for stratfold3d and penint3D ============================================
source(paste(getwd(),"R","stratFold3D.R",sep="/"))
source(paste(getwd(),"R","penint3D_def.R",sep="/"))
source(paste(getwd(),"R","plotfolds.R",sep="/"))
#source(paste(getwd(),"R","penint3Dpred.R",sep="/"))


fun<-as.formula(paste("Humus ~", paste(c(head(sm2D.lst,(length(sm2D.lst))),"altitude"), collapse="+")))
regdat<-bor.df
contVar<-sm2D.lst[-which(sm2D.lst %in% c("clc","SoilType"))]
targetVar<-"As"

#regdat$altitude<-2+regdat$altitude
#regdat[1:20,1:4]

#set.seed(432)
detach("package:dplyr", unload=TRUE)
strat<-stratfold3d(targetVar="As",seed=123,regdat=regdat,folds=5,cent=3,dimensions="2D",IDs=TRUE,sum=TRUE)
flist<-strat$folds
plotfolds(strat,"As")
#library(dplyr)

rez<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="linear")
rez$measure
summary(rez$measure[1:5,])


rezint<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="linear")
rezint$measure
summary(rezint$measure[1:5,])


rezinthier<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = TRUE, int=FALSE, depth.fun="linear")
rezinthier$measure
summary(rezinthier$measure[1:5,])


#============= Poly ========================
rez.poly<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="poly")
rez.poly$measure
summary(rez.poly$measure[1:5,])


rezint.poly<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
rezint.poly$measure
summary(rezint.poly$measure[1:5,])


rezinthier.poly<-penint3Drev(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="poly")
rezinthier.poly$measure
summary(rezinthier.poly$measure[1:5,])


#=================== predint3D ====================================================

#fun<-as.formula(paste("As ~", paste(c(head(sm2D.lst,(length(sm2D.lst))),"altitude"), collapse="+")))
#cogrids <- gridmaps.sm2D
#profs<-bor.profs
#int=TRUE; hier=FALSE; depth.fun="poly";pred=TRUE;depths=c(-.1,-.3); lambda=seq(0,5,0.1); deg=3;depth.fun="poly"; preProc=TRUE; cent=3; fold=5; seed=321;cores=8;chunk=20000;l=c(370000:450000)
#source(paste(getwd(),"R","stratFold3D.R",sep="/"))

#fun=fun; profs = bor.profs; cogrids = gridmaps.sm2D; pred=TRUE ;hier = TRUE; int=TRUE; depth.fun="linear";cores=8
#predint3D<-function(fun, profs, cogrids, hier=FALSE,pred=TRUE,lambda=seq(0,5,0.1),deg=3,fold=5,cent=3,int=TRUE,depth.fun=list("linear","poly"),depths=c(-.1,-.3),chunk=20000,preProc=TRUE,cores=2,seed=321,l=c(1:621426)){
source(paste(getwd(),"R","predint3D.R",sep="/"))


# ====================== int=FALSE, fun=linear =========================================

FFL<-predint3D(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=TRUE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
FFP<-predint3D(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=TRUE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

TFL<-predint3D(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=TRUE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
TFP<-predint3D(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=TRUE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

TTL<-predint3D(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=TRUE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
TTP<-predint3D(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=TRUE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

#============================ Prediction Summary =============================================

summary.cv<-list(FFL=FFL$summary$results,FFP=FFP$summary$results,TFL=TFL$summary$results,TFP=TFP$summary$results,TTL=TTL$summary$results,TTP=TTP$summary$results)
summary.coef<-list(FFL=FFL$summary$coefficients,FFP=FFP$summary$coefficients,TFL=TFL$summary$coefficients,TFP=TFP$summary$coefficients,TTL=TTL$summary$coefficients,TTP=TTP$summary$coefficients)
summary.pred<-list(FFL=FFL$summary$pred,FFP=FFP$summary$pred,TFL=TFL$summary$pred,TFP=TFP$summary$pred,TTL=TTL$summary$pred,TTP=TTP$summary$pred)


summary.cv
summary.coef
summary.pred

#============================ Prediction 0.2 =================================================

prediction<-FFL$prediction[[1]][,"pred"]
names(prediction)<-"FFL0.1"

prediction$FFP0.1<-FFP$prediction[[1]]$pred

fun<-as.formula(paste("As ~", paste(c(head(sm2D.lst,(length(sm2D.lst))),"altitude"), collapse="+")))
cogrids <- gridmaps.sm2D
profs<-bor.profs
int=TRUE; hier=FALSE; depth.fun="poly";pred=TRUE;depths=c(-.1,-.3); lambda=seq(0,5,0.1); deg=3;depth.fun="poly"; preProc=TRUE; cent=3; fold=5; seed=321;cores=8;chunk=20000;l=c(370000:450000)
source(paste(getwd(),"R","stratFold3D.R",sep="/"))

prediction$TFL0.1 <- TFL$prediction[[1]]$pred
prediction$TFP0.1 <- TFP$prediction[[1]]$pred


prediction$TTL0.1 <- TTL$prediction[[1]]$pred
prediction$TTP0.1 <- TTP$prediction[[1]]$pred


prediction$FFL0.2<-FFL$prediction[[2]]$pred
prediction$FFP0.2<-FFP$prediction[[2]]$pred

prediction$TFL0.2 <- TFL$prediction[[2]]$pred
prediction$TFP0.2 <- TFP$prediction[[2]]$pred

prediction$TTL0.2 <- TTL$prediction[[2]]$pred
prediction$TTP0.2 <- TTP$prediction[[2]]$pred


prediction$FFL0.3<-FFL$prediction[[3]]$pred
prediction$FFP0.3<-FFP$prediction[[3]]$pred

prediction$TFL0.3 <- TFL$prediction[[3]]$pred
prediction$TFP0.3 <- TFP$prediction[[3]]$pred

prediction$TTL0.3 <- TTL$prediction[[3]]$pred
prediction$TTP0.3 <- TTP$prediction[[3]]$pred

str(prediction)

sort(grep("L", names(prediction), value=TRUE))


pred.stack<-stack(prediction)

minR<-min(minValue(subset(pred.stack, subset=sort(grep("L", names(prediction), value=TRUE)))))
maxR<-max(maxValue(subset(pred.stack, subset=sort(grep("L", names(prediction), value=TRUE)))))

color.pal <- colorRampPalette(c("dark red","orange","light Yellow"),space="rgb")

plot(subset(pred.stack, subset=sort(grep("L", names(prediction), value=TRUE)) , col=color.pal, breaks=seq(minR,maxR,length.out=30)))


#=========================== Prediction 0.2 ===================================================



prediction0.1<-FFL$prediction[[1]][,"pred"]
names(prediction0.1)<-"FFL0.1"

prediction0.1$FFP0.1<-FFP$prediction[[1]]$pred

prediction0.1$TFL0.1<-TFL$prediction[[1]]$pred
prediction0.1$TFP0.1<-TFP$prediction[[1]]$pred

prediction0.1$TTL0.1<-TTL$prediction[[1]]$pred
prediction0.1$TTP0.1<-TTP$prediction[[1]]$pred

minv<-min(prediction0.1@data)
maxv<-max(prediction0.1@data)


pred.stack0.1<-stack(prediction0.1)

plot(pred.stack0.1, col=rev( heat.colors( 10 ) ),breaks=seq(minv,maxv,length.out=10))



save.image("pHPrediction2632016.RData")





str(TFL)
names(predint3D.As.g)

summary(predint3D.As.g$prediction[[1]]$pred)

(predint3D.As.g$summary$results)


pp<-raster(new3D[[2]],"pred")

plot(pp)

spplot(pp,"pred")


load("Results2432016.RData")




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


