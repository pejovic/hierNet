
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

#======================== CRS ====================================
gk_7 <- "+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m"
utm <- "+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m"

#================== Loading covariates ===========================
load(paste(getwd(),"inst","extdata","gridmaps.RDa",sep="/"))
gridmaps.sm2D$CD <- exp(-gridmaps.sm2D$DirDif)
gridmaps.sm2D$DD <- as.numeric(exp(-scale(gridmaps.sm2D$Dist,center=FALSE)))
#=================================================================

#================= Loading observations ==========================
bor <- join(read.csv(paste(getwd(),"inst","extdata","Profili_sredjeno_csv.csv",sep="/")), read.csv(paste(getwd(),"inst","extdata","Koordinate_csv.csv",sep="/")), type="inner")
bor$hdepth<-bor$Bottom-bor$Top
bor$altitude <- - (bor$Top / 100 + (( bor$Bottom - bor$Top ) / 2) / 100)
bor <- bor[, - c( 7:12,14,15,16,17 )]
names(bor) <- c("Soil.Type","ID","Horizont","Top" , "Bottom","pH","Humus","As","Co","x","y","hdepth","altitude")


#======================= SPDF ====================================
bor.As.xy <- bor[complete.cases(bor[,c("ID","x","y","altitude","As")]),c("ID","x","y","hdepth","altitude","As")] 
bor.As.xy <- plyr::rename(bor.As.xy, replace=c("x" = "longitude", "y" = "latitude"))
bor.As.xy<-ddply(bor.As.xy,.(ID),here(summarize),As=As[1],longitude=longitude[1],latitude=latitude[1])
bor.As.xy<-bor.As.xy[-which(bor.As.xy$ID %in% c(66,129,14,51,69,130,164,165,166)),]
coordinates(bor.As.xy) <- ~ longitude + latitude
proj4string(bor.As.xy) <- CRS(utm)
bor.As.xy <- spTransform(bor.As.xy, CRS(gk_7))
#=================================================================

#================== formula ======================================
fun <- as.formula(paste("As ~", paste(c("DEM","ES","CD","DD","VertDistChannelNet","ChannelNetBaseLevel","NegOpenness","Slope","TWI","clc","SoilType"), collapse="+")))

#================ all columns of interest: =======================
tv <- all.vars(fun)[1]
seln <- names(gridmaps.sm2D) %in% all.vars(fun)[-1]

covsgrids<-gridmaps.sm2D[,seln]
sm2D.lst<-names(covsgrids@data)


#======================= Overlay =================================
ov <- over(bor.As.xy,covsgrids)

bor.As.xy@data<-cbind(bor.As.xy@data,ov)
bor.As.df<-bor.As.xy@data
bor.As.df<-as.data.frame(bor.As.xy)

#=================================================================

#============= removing nzv and factor coding ====================

dummy.par <- dummyVars(as.formula(paste("~", paste(all.vars(fun)[-1], collapse="+"))),bor.As.df,levelsOnly=FALSE)
bor.As.df <- predict(dummy.par,newdata=bor.As.df)

nzv.par<-preProcess(bor.As.df,method=c("nzv"))
regdat<-as.data.frame(predict(nzv.par,bor.As.df))
head(bor.As.df)



## all columns of interest:
tv <- all.vars(fun)[1]
seln <- names(gridmaps.sm2D) %in% all.vars(fun)[-1]

## prepare regression matrix:
ov <- over(x=cogrids, y=geo, method=methodid, var.type = "numeric")
if(nrow(ov)==0|is.null(ov$observedValue)) {
  warning("The 'over' operations resulted in an empty set. Check 'methodid' column.")
}
## geostats only possible with numeric variables:
ov[,methodid] = as.numeric(ov$observedValue)
ov$observedValue = NULL

regdat<-ov[,all.vars(fun)[-1]]



################################################3
cont.par <- preProcess(regdat[,sapply(regdat,is.numeric)],method=c("center", "scale"))
regdat <- predict(cont.par,regdat)

dummy.par <- dummyVars(as.formula(paste("~", paste(all.vars(fun)[-1], collapse="+"))),regdat,levelsOnly=FALSE)
regdat <- predict(dummy.par,newdata=regdat)

nzv.par<-preProcess(regdat,method=c("nzv"))
regdat<-as.data.frame(predict(nzv.par,regdat))
head(regdat)

names(regdat)<-gsub(" ",".",names(regdat))



#aaa<-lapply(new3D, function(x) as.data.frame(x)) %>% lapply(.,function(x) predict(cont.par,newdata=x)) %>% lapply(.,function(x) predict(dummy.par,newdata=x)) %>% lapply(., function(x) as.data.frame(predict(nzv.par,newdata=x)))
#colnames(aaa)<-gsub(" ",".",colnames(aaa))
#XX <- hierNet::compute.interactions.c(as.matrix(aaa[[1]]),diagonal=FALSE)







fun <- as.formula(paste("As ~", paste(c(sm2D.lst), collapse="+")))

omm2 <- fit.gstatModel(bor.As.xy, fun, gridmaps.sm2D, method="randomForest")
plot(omm2)

As.reg.m<-omm2@regModel
varImpPlot(omm2@regModel)

As.reg.m<-randomForest(fun,data=(bor.As.df),ntree=500,nodesize=6,mtry=6,importance = TRUE,keep.inbag=TRUE)#,xtest=AsTest[,paste(c(sm2D.lst, "SoilType","altitude"))],ytest=AsTest[,"As"]

AsFF<-forestFloor(As.reg.m,bor.As.df)
plotmo(As.reg.m,all2=TRUE)
plotres(As.reg.m)
varImpPlot(As.reg.m)
pp<-partialPlot(As.reg.m, As.df,DEM )
#ice(As.reg.m,AsTrain[,paste(c(sm2D.lst, "SoilType","altitude"))],AsTrain[,"As"],predictor="altitude",)


Col=fcol(AsFF,1,orderByImportance=FALSE)
plot(AsFF,plot_seq=c(1:11),col=Col,plot_GOF=TRUE)
#ggPlotForestFloor(AsFF, 1:9)
show3d(AsFF,c(5,1),col=Col,plot_GOF=TRUE,orderByImportance=FALSE)



ov.H <- over(bor.As.xy, gridmaps.sm2D)
bor.As.xy@data[,c(sm2D.lst)] <- ov.H
str(bor.As.xy)

bor.As.df <- data.frame(bor.As.xy)
bor.As.df <- bor.As.df[complete.cases(bor.As.df[,names(bor.As.df)]),]
str(bor.As.df)

names(bor.As.df)
bor.As.df<-bor.As.df[,-which(names(bor.As.df) %in% c("AnalyticalHills","LSFactor","ChannelNetBaseLevel","RelSlopePosition","CatchArea","WindEffectNorthWest","VelleyDepth","optional"))]
str(bor.As.df)



#bor.As.df<-bor.As.df[-which(bor.As.df$ID %in% c(141)),] #153,137,141,145,150



