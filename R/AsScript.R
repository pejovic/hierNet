
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
library(MASS)
library(splines)
library(glmnet)
library(hierNet)
library(magrittr)
library(doParallel)
library(foreach)
library(stargazer)
library(gstat)

#path <- "C:/Users/User/Dropbox/Extensions of soil 3D trend models/Data and Scripts"

#============ Functions============================================================
source(paste(getwd(),"R","stratFold3D.R",sep="/"))
#source(paste(getwd(),"R","penint3D_def.R",sep="/"))
source(paste(getwd(),"R","plotfolds.R",sep="/"))
#source(paste(getwd(),"R","predint3D.R",sep="/"))
source(paste(getwd(),"R","penint3D_defP.R",sep="/"))
source(paste(getwd(),"R","predint3DP.R",sep="/"))
source(paste(getwd(),"R","3Dkrigencv.R",sep="/"))
source(paste(getwd(),"R","multiplot.R",sep = "/"))
source(paste(getwd(),"R","krige3Dpred.R",sep = "/"))
#==================================================================================

load("C:/Users/Milutin/Dropbox/Extensions of soil 3D trend models/Data and Scripts/BorData.rda")
load("C:/Users/Milutin/Dropbox/Extensions of soil 3D trend models/Data and Scripts/covmaps.rda")

gk_7 <- "+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m"
utm <- "+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m"

#================== Covariates table =================================================

CovNames <- c("Digital Elevation Model", "Aspect", "Slope","Topographic Wetness Index", "Convergence Index" ,"Cross Sectional Curvature", "Longitudinal Curvature", "Channel Network Base Level" ,"Vertical Distance to Channel Network", "Negative Openness","Positive Openness", "Wind Effect (East)","Wind Effect (North-West)","Down-wind Dilution", "Cross-wind Dilution" ,"Corine Land Cover 2006", "Soil Type")
CovAbb <- c("DEM","Aspect","Slope","TWI","ConvInd","CrSectCurv","LongCurv","ChNetBLevel","VDistChNet","NegOp", "PosOp","WEeast","WEnw","DD","CD","clc","SoilType")

#covs<-data.frame(Name=CovNames,
#                 Abbrevation=CovAbb,
#                Type=c("C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","F","F")
#                 )

#stargazer(covs,summary=FALSE,type="latex",digits=2)
#=====================================================================================

#=================== DATA ============================================================
bor.profs <- bor[,c("ID","x","y","Soil.Type","Top","Bottom","SOM","pH","Co","As")]
depths(bor.profs) <- ID ~ Top + Bottom
site(bor.profs) <- ~ Soil.Type + x + y
coordinates(bor.profs) <- ~ x+y
proj4string(bor.profs) <- CRS(utm)
#=====================================================================================

#====================== formulas ================================================================
As.fun <- as.formula(paste("As ~", paste(c(CovAbb,"altitude"), collapse="+")))
SOM.fun <- as.formula(paste("SOM ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"altitude"), collapse="+")))
pH.fun <- as.formula(paste("pH ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"altitude"), collapse="+")))

#================================================================================================

#=================================== plot stratified fold - Figure 2 ============================================
rdat <- bor
rdat <- plyr::rename(rdat, replace=c("x" = "longitude", "y" = "latitude"))
rdat <- rdat[complete.cases(rdat[,c("ID","longitude","latitude","altitude","SOM")]),c("ID","longitude","latitude","hdepth","altitude","As")] 

coordinates(rdat)<-~longitude+latitude
proj4string(rdat) <- CRS(utm)
rdat <- spTransform(rdat, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
rdat <- as.data.frame(rdat)
head(rdat)

rdat <- rdat[complete.cases(rdat$As),]


rdat.folds <- stratfold3d(targetVar="As",regdat=rdat,folds=5,cent=3,seed=666,dimensions="3D",IDs=TRUE,sum=TRUE)
plotfolds(rdat.folds,targetvar="As")

stargazer(do.call(rbind, rdat.folds$`As summary`), summary=FALSE, digits=2, type="latex")
stargazer(do.call(rbind,as.list(rdat.folds$`altitude summary`)[[1]]), summary=FALSE, digits=2, type="text")
#================================================
bor.profs$depth<-profileApply(bor.profs, function(x) estimateSoilDepth(x,name ="Horizont", top = "Top", bottom = "Bottom"))
spBor<-as(bor.profs,'SpatialPointsDataFrame')
spBor<-spTransform(spBor,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
str(spBor)


Bor <- as.data.frame(spBor)

q <- ggplot(Bor,aes(x = x, y = y))
r <- q +geom_point(aes(size = sqrt(depth/pi)), pch = 21, show.legend = FALSE) + scale_size_continuous(range=c(1,10))
#r <- r + facet_wrap(~ Soil.Type)
r <- r + aes(fill = Soil.Type)
r

#========== Plot to file ========================
pdf("SOMfolds.pdf",width=10,height=12)
plotfolds(rdat.folds,targetvar="As")
dev.off()
#================================================

#============== Table 3 ===================================================
stargazer(do.call(rbind, rdat.folds$`SOM summary`), summary=FALSE, digits=2, type="latex")
stargazer(do.call(rbind,as.list(rdat.folds$`altitude summary`)[[1]]), summary=FALSE, digits=2, type="text")

#============= plot to file - Figure 2 ====================================
library(ggplot2)
postscript("SOMfolds.pdf",horizontal=FALSE, paper="special",height=11,width=8)
plotfolds(rdat.folds,targetvar="SOM")
ggsave("SOMfolds.pdf", height=11,width=8)

dev.off()

#=================================== As =======================================================================

#==============================================================================================================
#===================== Nested Cross-validation of trend =======================================================
#==============================================================================================================

#================ Without interactions ========================================================================
As.BaseL.ncv <- penint3DP(fun=As.fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="linear")
As.BaseP.ncv <- penint3DP(fun=As.fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="poly")

#================ With interactions and without hierarchy =====================================================
As.IntL.ncv <- penint3DP(fun=As.fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="linear")
As.IntP.ncv <- penint3DP(fun=As.fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")

#================ With interactions and hierarchy =============================================================
As.IntHL.ncv <- penint3DP(fun=As.fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="linear")
As.IntHP.ncv <- penint3DP(fun=As.fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="poly")

#=================== Nested 5 fold cross-validation results for trend - Table 5 ===================================================

summary.n5cv<-rbind(As.BaseL.ncv = As.BaseL.ncv$measure[6,], As.BaseP.ncv = As.BaseP.ncv$measure[6,], As.IntL.ncv = As.IntL.ncv$measure[6,], As.IntP.ncv = As.IntP.ncv$measure[6,], As.IntHL.ncv = As.IntHL.ncv$measure[6,], As.IntHP.ncv = As.IntHP.ncv$measure[6,])
summary.n5cv[,1] <- rownames(summary.n5cv)
rownames(summary.n5cv)<-NULL
names(summary.n5cv)<-c("Model","RMSE","R squared")
stargazer(summary.n5cv,summary=FALSE,digits=2,type = "text")

latex.table.by(summary.n5cv, num.by.vars = 2)

#==================================================================================================================================
#============================================= krige3D.ncv ========================================================================


As.krige.res <- krige3D.ncv(As.fun,reg.ncv = As.IntP.ncv,profs=bor.profs,cogrids = gridmaps.sm2D, model=TRUE,krige=TRUE,v.cutoff = 60,v.width = 3,sp.cutoff = 4000,sp.width = 480, sp.vgm = vgm(1200, "Sph", 2000, 300),v.vgm = vgm(250, "Gau", 30, 5))

#============== Plot Variograms ================
multiplot(plotlist = As.krige.res$var1D,cols=2)
multiplot(plotlist = As.krige.res$var2D,cols=2)
multiplot(plotlist = As.krige.res$var3D,cols=2)

As.final.results <- do.call(rbind,As.krige.res$test.results)

caret::R2(pred=As.final.results$predicted,obs=As.final.results$observed, formula="traditional") # Rsquared for RK
RMSE(pred=As.final.results$predicted,obs=As.final.results$observed)

caret::R2(pred=As.final.results$final.predicted,obs=As.final.results$observed, formula="traditional") # Rsquared for RK
RMSE(pred=As.final.results$final.predicted,obs=As.final.results$observed)

#============================================================================================================================================

#============================================================================================================================================
#============== PREDICTION ==================================================================================================================
#============================================================================================================================================


#=================== As =====================================================================================================================

BaseL.As <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
BaseP.As <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntL.As <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntP.As <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntHL.As <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntHP.As <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = TRUE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

As.pred <- list(BaseL.As=BaseL.As$summary$pred, BaseP.As=BaseP.As$summary$pred, IntL.As=IntL.As$summary$pred, IntP.As=IntP.As$summary$pred, IntHL.As=IntHL.As$summary$pred , IntHP.As=IntHP.As$summary$pred)
As.coef <- list(BaseL.As=BaseL.As$summary$coefficients, BaseP.As=BaseP.As$summary$coefficients, IntL.As=IntL.As$summary$coefficients, IntP.As=IntP.As$summary$coefficients, IntHL.As=IntHL.As$summary$coefficients ,IntHP.As=IntHP.As$summary$coefficients)

As.pred <- lapply(As.pred, function(x) {names(x)<-c("ID", "longitude", "latitude","hdepth", "altitude","observed", "predicted","residual" ); return(x)})

#============================= Coefficients tables for As models ===========================================================================================
ll <- length(IntL.As$summary$coefficients)
pp <- length(IntHL.As$summary$coefficients[,1])+1


#============================= Coefficients for models with linear depth function ==============================================================
cmL <- data.frame(variable=IntHL.As$summary$coefficients[,1], BaseL.As.me=BaseL.As$summary$coefficients[2:pp], IntL.As.me=IntL.As$summary$coefficients[2:pp],IntL.As.ie=c(IntL.As$summary$coefficients[(pp+1):ll],0),IntHL.As.me=IntHL.As$summary$coefficients[,2],IntHL.As.ie=IntHL.As$summary$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.As$summary$coefficients)
p <- length(IntHP.As$summary$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.As <- data.frame(variable=IntHP.As$summary$coefficients[,1], BaseP.As.me=BaseP.As$summary$coefficients[2:p], IntP.As.me=IntP.As$summary$coefficients[2:p],IntP.As.ie1=c(IntP.As$summary$coefficients[(p+1):l][i1],0,0,0),IntP.As.ie2=c(IntP.As$summary$coefficients[(p+1):l][i2],0,0,0),IntP.As.ie3=c(IntP.As$summary$coefficients[(p+1):l][i3],0,0,0),IntHP.As.me=IntHP.As$summary$coefficients[,2],IntHP.As.ie1=IntHP.As$summary$coefficients[,3],IntHP.As.ie2=IntHP.As$summary$coefficients[,4],IntHP.As.ie3=IntHP.As$summary$coefficients[,5] )
cmAs <- cmP.As[,c(1,7:10)]
#cmAs <- cmP.As[,c(1,3:6)]
#==============================================================================================================================================
stargazer(cmP.As,summary = FALSE ,digits=3,type='latex')
#=========================== 3D Kriging prediction ============================================================================================
#=========================== The best model is IntP ===========================================================================================
As.kriging.pred <- krige3Dpred(fun=As.fun, reg.pred=IntP.As, profs=bor.profs, model = TRUE, krige = FALSE, v.cutoff=40, v.width=3, sp.vgm = vgm(1200, "Sph", 2000, 300),v.vgm = vgm(250, "Gau", 30, 5), sp.cutoff=4000, sp.width=480)

plot(As.kriging.pred$var1D) # Vertical variogram
plot(As.kriging.pred$var2D) # Horizontal variogram
plot(As.kriging.pred$var3D) # 3D variogram

pdf("AsIntPVariogram1D.pdf",width=5,height=4)
plot(As.kriging.pred$var1D)
dev.off()

pdf("AsIntPVariogram2D.pdf",width=5,height=4)
plot(As.kriging.pred$var2D)
dev.off()

pdf("AsIntPVariogram3D.pdf",width=5,height=4)
plot(As.kriging.pred$var3D)
dev.off()

#=== Anistotropy parameter ==============
As.kriging.pred$anis

#===============================================================================================================================================


#===================== Prediction plots ========================================================================================================


str(As.kriging.pred$prediction)

prediction<-As.kriging.pred$prediction[[1]][,"final.prediction"]
names(prediction)<-"As0.1"

prediction$As0.2 <- As.kriging.pred$prediction[[2]]$final.prediction
prediction$As0.3 <- As.kriging.pred$prediction[[3]]$final.prediction

prediction$As0.1 <- pmax(prediction$As0.1,min(na.omit(bor$As))/3)
prediction$As0.2 <- pmax(prediction$As0.2,min(na.omit(bor$As))/3)
prediction$As0.3 <- pmax(prediction$As0.3,min(na.omit(bor$As))/3)

pred.stack<-stack(prediction)

gridpix <- as(pred.stack,"SpatialPixelsDataFrame")

zmax <- max(gridpix$As0.1, gridpix$As0.2,gridpix$As0.3)
zmin <- min(gridpix$As0.1, gridpix$As0.2,gridpix$As0.3)
zmax <- round(zmax) 
zmin <- round(zmin) 

ramp <- seq(from = zmin, to = zmax, by = 3)
color.As <- colorRampPalette(c("light blue","yellow","orange","dark red"))

ckey <- list(labels=list(cex=2))

p1 <- spplot(gridpix, "As0.1", asp = 1, at = ramp, col.regions = color.As,colorkey=FALSE)
p2 <- spplot(gridpix, "As0.2", asp = 1, at = ramp,  col.regions = color.As,colorkey=FALSE)
p3 <- spplot(gridpix, "As0.3", asp = 1, at = ramp, col.regions = color.As,colorkey=ckey) # last plot contains the legend

#================== plot to files ======================================================
pdf("Asplot0.1.pdf",width=6,height=12)
p1
dev.off()

pdf("Asplot0.2.pdf",width=6,height=12)
p2
dev.off()


pdf("Asplot0.3.pdf",width=7,height=12)
p3
dev.off()
#========================================================================================


save.image("D:/_Bor/Drugi rad/AsScript.RData")

save("bor",file="C:/Users/User/Dropbox/Extensions of soil 3D trend models/Data and Scripts/BorData.rda")
save("gridmaps.sm2D",file="C:/Users/User/Dropbox/Extensions of soil 3D trend models/Data and Scripts/covmaps.rda")

C:\Users\User\Dropbox\Extensions of soil 3D trend models\Data and Scripts













#============================ Prediction Summary =============================================

summary.cv<-list(BaseL=BaseL$summary$results[6,],BaseP=BaseP$summary$results[6,],IntL=IntL$summary$results[6,],IntP=IntP$summary$results[6,],IntHL=IntHL$summary$results[6,],IntHP=IntHP$summary$results[6,])
summary.coef<-list(BaseL=BaseL$summary$coefficients,BaseP=BaseP$summary$coefficients,IntL=IntL$summary$coefficients,IntP=IntP$summary$coefficients,IntHL=IntHL$summary$coefficients,IntHP=IntHP$summary$coefficients)
summary.pred<-list(BaseL=BaseL$summary$pred,BaseP=BaseP$summary$pred,IntL=IntL$summary$pred,IntP=IntP$summary$pred,IntHL=IntHL$summary$pred,IntHP=IntHP$summary$pred)


summary.cv
summary.coef
summary.pred

str(summary.pred)
lapply(summary.pred, function(x) {names(x)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")}, return(summary.pred))


#============================ Prediction 0.2 =================================================
fun <- As.fun

BaseL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
BaseP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntP <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=TRUE ,hier = FALSE, int=TRUE, depths=c(-0.1) ,depth.fun="poly",cores=8)

fun=As.fun; profs = bor.profs; cogrids = gridmaps.sm2D; pred=TRUE ;hier = FALSE; int=TRUE; depths=c(-0.1) ;depth.fun="poly";cores=8
lambda=seq(0,5,0.1);deg=3;fold=5;cent=3;chunk=20000;preProc=TRUE;seed=321


IntHL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
system.time(IntHP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1) ,depth.fun="poly",cores=8))

prediction<-FFL$prediction[[1]][,"pred"]
names(prediction)<-"FFL0.1"

prediction$FFP0.1<-FFP$prediction[[1]]$pred

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

plot(subset(pred.stack, subset=sort(grep("IntL", names(prediction), value=TRUE)) , col=color.pal, breaks=seq(minR,maxR,length.out=30)),axes=FALSE)

gridpix <- as(pred.stack,"SpatialPixelsDataFrame")

color.pal <- colorRampPalette(c("dark red","orange","light Yellow"))
spplot(gridpix,"IntL0.1",col.regions = color.pal)

str(prediction)
#### ====================== Prediction =============================================
zmax <- max(gridpix$IntHP0.1, gridpix$IntHP0.2,gridpix$IntHP0.3)

zmin <- min(gridpix$IntHP0.1, gridpix$IntHP0.2,gridpix$IntHP0.3)

zmax <- round(zmax) 

zmin <- round(zmin) 

ramp <- seq(from = zmin, to = zmax, by = 3)

color.pH <- colorRampPalette(c("dark red","red","orange","yellow","light green","green","light blue","blue","dark blue"))

color.SOM <- colorRampPalette(c("dark red","orange","light Yellow"))
color.As <- colorRampPalette(c("light blue","yellow","orange","dark red"))
#color.pal <- brewer.pal(8,name="YlOrRd")

ckey <- list(labels=list(cex=2))

p1 <- spplot(gridpix, "IntHP0.1", asp = 1, at = ramp, col.regions = color.As,colorkey=FALSE)
p2 <- spplot(gridpix, "IntHP0.2", asp = 1, at = ramp,  col.regions = color.As,colorkey=FALSE)
p3 <- spplot(gridpix, "IntHP0.3", asp = 1, at = ramp, col.regions = color.As,colorkey=ckey)
p1

pdf("Asplot0.1.pdf",width=6,height=12)
p1
dev.off()

pdf("Asplot0.2.pdf",width=6,height=12)
p2
dev.off()


pdf("Asplot0.3.pdf",width=7,height=12)
p3
dev.off()














































#=================== PREDICTION ==============================================================================================================

#=================== predint3D ===============================================================================================================

#=================== SOM =====================================================================================================================

BaseL.SOM <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
BaseP.SOM <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntL.SOM <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntP.SOM <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntHL.SOM <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntHP.SOM <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

SOM.pred <- list(BaseL.SOM=BaseL.SOM$summary$pred, BaseP.SOM=BaseP.SOM$summary$pred, IntL.SOM=IntL.SOM$summary$pred, IntP.SOM=IntP.SOM$summary$pred, IntHL.SOM=IntHL.SOM$summary$pred , IntHP.SOM=IntHP.SOM$summary$pred)
SOM.coef <- list(BaseL.SOM=BaseL.SOM$summary$coefficients, BaseP.SOM=BaseP.SOM$summary$coefficients, IntL.SOM=IntL.SOM$summary$coefficients, IntP.SOM=IntP.SOM$summary$coefficients, IntHL.SOM=IntHL.SOM$summary$coefficients ,IntHP.SOM=IntHP.SOM$summary$coefficients)

SOM.pred <- lapply(SOM.pred, function(x) {names(x)<-c("ID", "longitude", "latitude","hdepth", "altitude","observed", "predicted","residual" ); return(x)})
str(SOM.pred)

#============================= Coefficients tables for all SOM models ==========================================================================
ll <- length(IntL.SOM$summary$coefficients)
pp <- length(IntHL.SOM$summary$coefficients[,1])+1

#============================= Coefficients for models with linear depth function ==============================================================
cmL.SOM <- data.frame(variable=IntHL.SOM$summary$coefficients[,1], BaseL.SOM.me=BaseL.SOM$summary$coefficients[2:pp], IntL.SOM.me=IntL.SOM$summary$coefficients[2:pp],IntL.SOM.ie=c(IntL.SOM$summary$coefficients[(pp+1):ll],0),IntHL.SOM.me=IntHL.SOM$summary$coefficients[,2],IntHL.SOM.ie=IntHL.SOM$summary$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.SOM$summary$coefficients)
p <- length(IntHP.SOM$summary$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.SOM <- data.frame(variable=IntHP.SOM$summary$coefficients[,1], BaseP.SOM.me=BaseP.SOM$summary$coefficients[2:p], IntP.SOM.me=IntP.SOM$summary$coefficients[2:p],IntP.SOM.ie1=c(IntP.SOM$summary$coefficients[(p+1):l][i1],0,0,0),IntP.SOM.ie2=c(IntP.SOM$summary$coefficients[(p+1):l][i2],0,0,0),IntP.SOM.ie3=c(IntP.SOM$summary$coefficients[(p+1):l][i3],0,0,0),IntHP.SOM.me=IntHP.SOM$summary$coefficients[,2],IntHP.SOM.ie1=IntHP.SOM$summary$coefficients[,3],IntHP.SOM.ie2=IntHP.SOM$summary$coefficients[,4],IntHP.SOM.ie3=IntHP.SOM$summary$coefficients[,5] )
cmSOM <- cmP.SOM[,c(1,7:10)] # extraction only for SOM IntPH final models coefficients
#cmSOM <- cmP.SOM[,c(1,3:6)] # extraction only for SOM IntP.SOM final models coefficients

stargazer(cmP.SOM, summary=FALSE, digits=3,type='latex')

SOM.kriging.pred <- krige3Dpred(fun=SOM.fun, reg.pred=IntP.SOM, profs=bor.profs, model = TRUE, krige = FALSE, v.cutoff=40, v.width=3, sp.vgm = vgm(15, "Sph", 2000, 5),v.vgm = vgm(1.5, "Gau", 10, 0), sp.cutoff=4000, sp.width=420)

plot(SOM.kriging.pred$var1D) # Vertical variogram
plot(SOM.kriging.pred$var2D) # Horizontal variogram
plot(SOM.kriging.pred$var3D) # 3D variogram

pdf("SOMIntPVariogram1D.pdf",width=5,height=4)
plot(SOM.kriging.pred$var1D)
dev.off()

pdf("SOMIntPVariogram2D.pdf",width=5,height=4)
plot(SOM.kriging.pred$var2D)
dev.off()

pdf("SOMIntPVariogram3D.pdf",width=5,height=4)
plot(SOM.kriging.pred$var3D)
dev.off()


#=================== pH =========================================================================================================================

BaseL.pH <- predint3DP(fun=pH.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
BaseP.pH <- predint3DP(fun=pH.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntL.pH <- predint3DP(fun=pH.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntP.pH <- predint3DP(fun=pH.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntHL.pH <- predint3DP(fun=pH.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntHP.pH <- predint3DP(fun=pH.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

pH.pred <- list(BaseL.pH=BaseL.pH$summary$pred, BaseP.pH=BaseP.pH$summary$pred, IntL.pH=IntL.pH$summary$pred, IntP.pH=IntP.pH$summary$pred, IntHL.pH=IntHL.pH$summary$pred , IntHP.pH=IntHP.pH$summary$pred)
pH.coef <- list(BaseL.pH=BaseL.pH$summary$coefficients, BaseP.pH=BaseP.pH$summary$coefficients, IntL.pH=IntL.pH$summary$coefficients, IntP.pH=IntP.pH$summary$coefficients, IntHL.pH=IntHL.pH$summary$coefficients ,IntHP.pH=IntHP.pH$summary$coefficients)

pH.pred <- lapply(pH.pred, function(x) {names(x)<-c("ID", "longitude", "latitude","hdepth", "altitude","observed", "predicted","residual" ); return(x)})

#============================= Coefficients tables for pH models =============================================================================================

ll <- length(IntL.pH$summary$coefficients)
pp <- length(IntHL.pH$summary$coefficients[,1])+1


#============================= Coefficients for models with linear depth function ==============================================================
cmL.pH <- data.frame(variable=IntHL.pH$summary$coefficients[,1], BaseL.pH.me=BaseL.pH$summary$coefficients[2:pp], IntL.pH.me=IntL.pH$summary$coefficients[2:pp],IntL.pH.ie=c(IntL.pH$summary$coefficients[(pp+1):ll],0),IntHL.pH.me=IntHL.pH$summary$coefficients[,2],IntHL.pH.ie=IntHL.pH$summary$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.pH$summary$coefficients)
p <- length(IntHP.pH$summary$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.pH <- data.frame(variable=IntHP.pH$summary$coefficients[,1], BaseP.pH.me=BaseP.pH$summary$coefficients[2:p], IntP.pH.me=IntP.pH$summary$coefficients[2:p],IntP.pH.ie1=c(IntP.pH$summary$coefficients[(p+1):l][i1],0,0,0),IntP.pH.ie2=c(IntP.pH$summary$coefficients[(p+1):l][i2],0,0,0),IntP.pH.ie3=c(IntP.pH$summary$coefficients[(p+1):l][i3],0,0,0),IntHP.pH.me=IntHP.pH$summary$coefficients[,2],IntHP.pH.ie1=IntHP.pH$summary$coefficients[,3],IntHP.pH.ie2=IntHP.pH$summary$coefficients[,4],IntHP.pH.ie3=IntHP.pH$summary$coefficients[,5] )
cmpH <- cmP.pH[,c(1,7:10)]
#cmpH <- cmP.pH[,c(1,3:6)]

stargazer(cmP.pH,summary=FALSE,digits=3,type='latex')

#======================== Kriging ==============================================

pH.kriging.pred <- krige3Dpred(fun=pH.fun, reg.pred=IntP.pH, profs=bor.profs, model = TRUE, krige = FALSE, v.cutoff=50, v.width=3, sp.vgm = vgm(0.3, "Sph", 2000, 0.1),v.vgm = vgm(0.05, "Gau", 20, 0), sp.cutoff=4000, sp.width=585)

plot(pH.kriging.pred$var1D) # Vertical variogram
plot(pH.kriging.pred$var2D) # Horizontal variogram
plot(pH.kriging.pred$var3D) # 3D variogram

pdf("pHIntPVariogram1D.pdf",width=5,height=4)
plot(pH.kriging.pred$var1D)
dev.off()

pdf("pHIntPVariogram2D.pdf",width=5,height=4)
plot(pH.kriging.pred$var2D)
dev.off()

pdf("pHIntPVariogram3D.pdf",width=5,height=4)
plot(pH.kriging.pred$var3D)
dev.off()


#=================== As =====================================================================================================================

BaseL.As <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
BaseP.As <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntL.As <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntP.As <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntHL.As <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntHP.As <- predint3DP(fun=As.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

As.pred <- list(BaseL.As=BaseL.As$summary$pred, BaseP.As=BaseP.As$summary$pred, IntL.As=IntL.As$summary$pred, IntP.As=IntP.As$summary$pred, IntHL.As=IntHL.As$summary$pred , IntHP.As=IntHP.As$summary$pred)
As.coef <- list(BaseL.As=BaseL.As$summary$coefficients, BaseP.As=BaseP.As$summary$coefficients, IntL.As=IntL.As$summary$coefficients, IntP.As=IntP.As$summary$coefficients, IntHL.As=IntHL.As$summary$coefficients ,IntHP.As=IntHP.As$summary$coefficients)

As.pred <- lapply(As.pred, function(x) {names(x)<-c("ID", "longitude", "latitude","hdepth", "altitude","observed", "predicted","residual" ); return(x)})

#============================= Coefficients tables for As models ===========================================================================================
ll <- length(IntL.As$summary$coefficients)
pp <- length(IntHL.As$summary$coefficients[,1])+1


#============================= Coefficients for models with linear depth function ==============================================================
cmL <- data.frame(variable=IntHL.As$summary$coefficients[,1], BaseL.As.me=BaseL.As$summary$coefficients[2:pp], IntL.As.me=IntL.As$summary$coefficients[2:pp],IntL.As.ie=c(IntL.As$summary$coefficients[(pp+1):ll],0),IntHL.As.me=IntHL.As$summary$coefficients[,2],IntHL.As.ie=IntHL.As$summary$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.As$summary$coefficients)
p <- length(IntHP.As$summary$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.As <- data.frame(variable=IntHP.As$summary$coefficients[,1], BaseP.As.me=BaseP.As$summary$coefficients[2:p], IntP.As.me=IntP.As$summary$coefficients[2:p],IntP.As.ie1=c(IntP.As$summary$coefficients[(p+1):l][i1],0,0,0),IntP.As.ie2=c(IntP.As$summary$coefficients[(p+1):l][i2],0,0,0),IntP.As.ie3=c(IntP.As$summary$coefficients[(p+1):l][i3],0,0,0),IntHP.As.me=IntHP.As$summary$coefficients[,2],IntHP.As.ie1=IntHP.As$summary$coefficients[,3],IntHP.As.ie2=IntHP.As$summary$coefficients[,4],IntHP.As.ie3=IntHP.As$summary$coefficients[,5] )
cmAs <- cmP.As[,c(1,7:10)]
#cmAs <- cmP.As[,c(1,3:6)]
#==============================================================================================================================================

#============== Combining IntPH model coefficients for all variables in one table =============================================================

cm <- list(cmAs = cmAs, cmpH = cmpH, cmSOM = cmSOM)
save(cm, file="cm.rda")

cm1 <- cbind(cm$cmSOM,cm$cmpH[,-1])
cm1[,1] <- as.character(cm1[,1])
cm1 <- rbind(cm1[1:26,-1],rep(0,8),rep(0,8), cm1[27:29,-1])
cm1 <- cbind(cm$cmAs$variable,cm1)
cm1 <- cbind(cm$cmAs,cm1[,-1])
names(cm1) <- c("variable","As.me","As.ie1","As.ie2","As.ie3","SOM.me","SOM.ie1","SOM.ie2","SOM.ie3","pH.me","pH.ie1","pH.ie2","pH.ie3")

stargazer(cm1,summary=FALSE, digits=3, type="latex")

intCoefs <- cm1
intHCoefs <- cm1

stargazer(intHCoefs,summary=FALSE, digits=3, type="latex")

#===============================================================================================================================================

IntP <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=TRUE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntP.pred <- krige3Dpred(fun=SOM.fun, reg.pred=IntP, profs=bor.profs, model = TRUE, krige=TRUE, v.cutoff=60, v.width=3, v.vgm = vgm(1.5, "Gau", 10, 0), sp.cutoff=4000, sp.width=420, sp.vgm = vgm(15, "Sph", 2000, 5))

plot(IntP.pred$var1D)
plot(IntP.pred$var2D)
plot(IntP.pred$var3D)

str(IntP.pred$)

summary(IntP.pred$prediction[[3]]$res.pred)


#================================================================================================================================================


