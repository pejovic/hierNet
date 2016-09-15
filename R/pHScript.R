# Demonstration of penint3D functions for evaluation and prediction of 3D interaction models 
# Target variable : Arsenic
# Data: Bor

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

#path <- "~/Dropbox/Extensions of soil 3D trend models/Data and Scripts"
path <- "C:/Users/user/Dropbox/Extensions of soil 3D trend models/Data and Scripts"
path1 <- "D:/R_projects/int3D/R"
#============ Functions============================================================
source(paste(path1,"stratFold3D.R",sep="/"))
source(paste(path1,"plotfolds.R",sep="/"))
#source(paste(path1,"penint3D.R",sep="/"))
source(paste(path1,"predint3D.R",sep="/"))
source(paste(path1,"3Dkrigencv.R",sep="/"))
source(paste(path1,"multiplot.R",sep = "/"))
source(paste(path1,"krige3Dpred.R",sep = "/"))
source(paste(path1,"penint3D_defP.R",sep="/"))
#==================================================================================

load(paste(path,"BorData.rda",sep = "/"))
load(paste(path,"covmaps.rda",sep = "/"))

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
bor.profs <- bor[,c("ID","x","y","Soil.Type","Top","Bottom","SOM","pH","Co","SOM")]
bor.profs$logpH <- log(bor.profs$pH)
bor.profs$mid <- (bor.profs$Top+bor.profs$Bottom)/2
depths(bor.profs) <- ID ~ Top + Bottom
site(bor.profs) <- ~ Soil.Type + x + y
coordinates(bor.profs) <- ~ x+y
proj4string(bor.profs) <- CRS(utm)
#=====================================================================================

#====================== formulas ================================================================
As.fun <- as.formula(paste("SOM ~", paste(c(CovAbb,"altitude"), collapse="+")))
SOM.fun <- as.formula(paste("SOM ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"altitude"), collapse="+")))
pH.fun <- as.formula(paste("pH ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"altitude"), collapse="+")))

#================================================================================================

#=================================== plot stratified fold - Figure 2 ============================================
rdat <- bor
rdat <- plyr::rename(rdat, replace=c("x" = "longitude", "y" = "latitude"))
rdat <- rdat[complete.cases(rdat[,c("ID","longitude","latitude","altitude","pH")]),c("ID","longitude","latitude","hdepth","altitude","pH")] 

coordinates(rdat)<-~longitude+latitude
proj4string(rdat) <- CRS(utm)
rdat <- spTransform(rdat, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
rdat <- as.data.frame(rdat)
head(rdat)

rdat <- rdat[complete.cases(rdat$pH),]

rdat.folds <- stratfold3d(targetVar="pH",regdat=rdat,folds=5,cent=3,seed=666,dimensions="3D",IDs=TRUE,sum=TRUE)
plotfolds(rdat.folds,targetvar="pH")

#========== Plot to file ========================
pdf("pHfolds.pdf",width=10,height=12)
plotfolds(rdat.folds,targetvar="pH")
dev.off()
#================================================

#============== Table 3 ===================================================
stargazer(do.call(rbind, rdat.folds$`pH summary`), summary=FALSE, digits=2, type="latex")
stargazer(do.call(rbind,as.list(rdat.folds$`altitude summary`)[[1]]), summary=FALSE, digits=2, type="text")

#============= plot to file - Figure 2 ====================================
library(ggplot2)
postscript("pHfolds.pdf",horizontal=FALSE, paper="special",height=11,width=8)
plotfolds(rdat.folds,targetvar="pH")
ggsave("pHfolds.pdf", height=11,width=8)

dev.off()

#=================================== pH =======================================================================

#==============================================================================================================
#===================== Nested Cross-validation of trend =======================================================
#==============================================================================================================

#================ Without interactions ========================================================================
pH.BaseL.ncv <- penint3DP(fun=pH.fun, profs = bor.profs, seed=321, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="linear")
pH.BaseP.ncv <- penint3DP(fun=pH.fun, profs = bor.profs, seed=321, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="poly")

#================ With interactions and without hierarchy =====================================================
pH.IntL.ncv <- penint3DP(fun=pH.fun, profs = bor.profs, seed=321, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="linear")
pH.IntP.ncv <- penint3DP(fun=pH.fun, profs = bor.profs, seed=111, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")

#================ With interactions and hierarchy (takes more than 40 min!!! ) =============================================================
pH.IntHL.ncv <- penint3DP(fun=pH.fun, profs = bor.profs, seed=321, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="linear")
pH.IntHP.ncv <- penint3DP(fun=pH.fun, profs = bor.profs, seed=321, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="poly")
#================================
#================================ Checking the variance of residuals along the depth for pH ans Log(pH) models ==================
pH.fun <- as.formula(paste("log(pH) ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"altitude"), collapse="+")))
logpH.fun <- as.formula(paste("logpH ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"altitude"), collapse="+")))

pH.IntP.ncv <- penint3DP(fun=pH.fun, profs = bor.profs, seed=444, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
logpH.IntP.ncv <- penint3DP(fun=logpH.fun, profs = bor.profs, seed=444, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")


Int<-pH.IntP.ncv$test.results
logInt<-logpH.IntP.ncv$test.results


Int<-do.call(rbind,Int)
logInt<-do.call(rbind,logInt)



Int$residual <- Int$observed-Int$predicted
logInt$residual <- logInt$observed-logInt$predicted

d=0
#Int results
caret::R2(pred=Int[Int$altitude < d, "predicted"] ,obs=Int[Int$altitude < d, "observed"]) # Rsquared for RK
RMSE(pred=Int[Int$altitude < d, "predicted"] ,obs=Int[Int$altitude < d, "observed"])
#logInt results
caret::R2(pred=logInt[logInt$altitude < d, "predicted"] ,obs=logInt[logInt$altitude < d, "observed"]) # Rsquared for RK
RMSE(pred=logInt[logInt$altitude < d, "predicted"] ,obs=logInt[logInt$altitude < d, "observed"])


dseq <- seq(-0.6,0,0.1)
dseq[1] <- -1
dseq[2] <- -0.6
#R2.res <- rep(NA,length(dseq)-1)
var.res <- rep(NA,length(dseq)-1)
log.var.res <- rep(NA,length(dseq)-1)
ni <- rep(NA,length(dseq)-1)
log.ni <- rep(NA,length(dseq)-1)

for(i in 1:length(dseq)-1){
  d <- dseq[i+1]
  #R2.res[i] <- caret::R2(pred=Base[Base$altitude < d, "predicted"] ,obs=Base[Base$altitude < d, "observed"])
  var.res[i] <- var(x=(Int[Int$altitude <= d & Int$altitude > dseq[i] , "observed"]-Int[Int$altitude <= d & Int$altitude > dseq[i], "predicted"]))
  ni[i] <- length((Int[Int$altitude <= d & Int$altitude > dseq[i], "observed"]-Int[Int$altitude <= d & Int$altitude > dseq[i], "predicted"]))
}

for(i in 1:length(dseq)-1){
  d <- dseq[i+1]
  #R2.res[i] <- caret::R2(pred=Int[Int$altitude < d, "predicted"] ,obs=Int[Int$altitude < d, "observed"])
  log.var.res[i] <- var(x=(logInt[logInt$altitude <= d & Int$altitude > dseq[i], "observed"]-logInt[logInt$altitude <= d & Int$altitude > dseq[i], "predicted"]))
  log.ni[i] <- length((logInt[logInt$altitude <= d & logInt$altitude > dseq[i], "observed"]-logInt[logInt$altitude <= d & logInt$altitude > dseq[i], "predicted"]))
}

int.var.res <- var.res
log.int.var.res <- log.var.res

res.var <- data.frame(depth=c("0-10 cm","10-20 cm","20-30 cm","30-40 cm","40-60 cm","60- 100 cm"),
                      var.res=rev(int.var.res),
                      var.log.res=rev(log.int.var.res),
                      n.obs=rev(ni))

stargazer(res.var,summary = FALSE,digits=2,type='text')
#==================================================================================================================================

#================================
#=================== Nested 5 fold cross-validation results for trend - Table 5 ===================================================

summary.n5cv<-rbind(pH.BaseL.ncv = pH.BaseL.ncv$measure[6,], pH.BaseP.ncv = pH.BaseP.ncv$measure[6,], pH.IntL.ncv = pH.IntL.ncv$measure[6,], pH.IntP.ncv = pH.IntP.ncv$measure[6,], pH.IntHL.ncv = pH.IntHL.ncv$measure[6,], pH.IntHP.ncv = pH.IntHP.ncv$measure[6,])
summary.n5cv[,1] <- rownames(summary.n5cv)
rownames(summary.n5cv)<-NULL
names(summary.n5cv)<-c("Model","RMSE","R squared")
stargazer(summary.n5cv,summary=FALSE,digits=2,type = "text")

latex.table.by(summary.n5cv, num.by.vars = 2)

#==================================================================================================================================
#============================================= krige3D.ncv ========================================================================

# Try with model=FALSE and krige = FALSE to examine the experimental variograms. If the default model parameters are good enough for initial model parameters, run with mode=TRUE and krige=TRUE. 

pH.krige.res <- krige3D.ncv(pH.fun,reg.ncv = pH.IntP.ncv,profs=bor.profs,cogrids = gridmaps.sm2D, model=TRUE,krige=TRUE,v.cutoff = 40,v.width = 3,sp.cutoff = 3000,sp.width = 650, v.vgm = vgm(0.05, "Gau", 20, 0),sp.vgm = vgm(0.3, "Sph", 2000, 0.1))

#============== Plot Variograms (experimental with model=FALSE and krige=FALSE. If model=TRUE and krige=TRUE the modeled variograms will be ploted) ================
multiplot(plotlist = pH.krige.res$var1D,cols=2)
multiplot(plotlist = pH.krige.res$var2D,cols=2)
multiplot(plotlist = pH.krige.res$var3D,cols=2)

pH.final.results <- do.call(rbind,pH.krige.res$test.results)

#============== Results of nested cv for the best trend model ==================================================
caret::R2(pred=pH.final.results$predicted,obs=pH.final.results$observed, formula="traditional") # Rsquared for trend
RMSE(pred=pH.final.results$predicted,obs=pH.final.results$observed) # RMSE for trend
#===============================================================================================================

#============== Results of nested cv for the 3D regression kriging ==================================================
caret::R2(pred=pH.final.results$final.predicted,obs=pH.final.results$observed, formula="traditional") # Rsquared for RK
RMSE(pred=pH.final.results$final.predicted,obs=pH.final.results$observed) # RMSE for RK
#=================================================================================================================
#============================================================================================================================================

#============================================================================================================================================
#============== PREDICTION ==================================================================================================================
#============================================================================================================================================


#=================== pH =====================================================================================================================
#====== Without interactions ================================================================================================================
BaseL.pH <- predint3DP(fun=pH.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8, seed = 111)
BaseP.pH <- predint3DP(fun=pH.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8, seed = 111)

#====== With interactions but without hierarchy ================================================================================================================
IntL.pH <- predint3DP(fun=pH.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8, seed = 111)
IntP.pH <- predint3DP(fun=pH.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8, seed = 111)
# IntP.pH model is the best. This is the reason why the pred=TRUE. If pred=TRUE the prediction over grid will be performed. If pred=FALSE, only the prediction on observed points will be performed

#====== With interactions honoring the hierarchy ================================================================================================================
#==================================== Require free memory and takes more than 40 min!!! =====================================================
IntHL.pH <- predint3DP(fun=pH.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8, seed = 111)
IntHP.pH <- predint3DP(fun=pH.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8, seed = 111)
#============================================================================================================================================

#============ Plot glmnet model ======================
pdf("LassoCoef.pdf",width=8,height=6)
plot(IntL.pH$summary$model$glmnet.fit,xvar = "lambda", label = TRUE)
abline(v=log(c(IntL.pH$summary$model$lambda.min)),lty=2)
dev.off()

plot(IntL.pH$summary$model)
#====================================================


pH.pred <- list(BaseL.pH=BaseL.pH$summary$pred, BaseP.pH=BaseP.pH$summary$pred, IntL.pH=IntL.pH$summary$pred, IntP.pH=IntP.pH$summary$pred, IntHL.pH=IntHL.pH$summary$pred , IntHP.pH=IntHP.pH$summary$pred)
pH.coef <- list(BaseL.pH=BaseL.pH$summary$coefficients, BaseP.pH=BaseP.pH$summary$coefficients, IntL.pH=IntL.pH$summary$coefficients, IntP.pH=IntP.pH$summary$coefficients, IntHL.pH=IntHL.pH$summary$coefficients ,IntHP.pH=IntHP.pH$summary$coefficients)

pH.pred <- lapply(pH.pred, function(x) {names(x)<-c("ID", "longitude", "latitude","hdepth", "altitude","observed", "predicted","residual" ); return(x)})

#============================= Coefficients tables for pH models ===========================================================================================
ll <- length(IntL.pH$summary$coefficients)
pp <- length(IntHL.pH$summary$coefficients[,1])+1


#============================= Coefficients for models with linear depth function ==============================================================
cmL <- data.frame(variable=IntHL.pH$summary$coefficients[,1], BaseL.pH.me=BaseL.pH$summary$coefficients[2:pp], IntL.pH.me=IntL.pH$summary$coefficients[2:pp],IntL.pH.ie=c(IntL.pH$summary$coefficients[(pp+1):ll],0),IntHL.pH.me=IntHL.pH$summary$coefficients[,2],IntHL.pH.ie=IntHL.pH$summary$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.pH$summary$coefficients)
p <- length(IntHP.pH$summary$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.pH <- data.frame(variable=IntHP.pH$summary$coefficients[,1], BaseP.pH.me=BaseP.pH$summary$coefficients[2:p], IntP.pH.me=IntP.pH$summary$coefficients[2:p],IntP.pH.ie1=c(IntP.pH$summary$coefficients[(p+1):l][i1],0,0,0),IntP.pH.ie2=c(IntP.pH$summary$coefficients[(p+1):l][i2],0,0,0),IntP.pH.ie3=c(IntP.pH$summary$coefficients[(p+1):l][i3],0,0,0),IntHP.pH.me=IntHP.pH$summary$coefficients[,2],IntHP.pH.ie1=IntHP.pH$summary$coefficients[,3],IntHP.pH.ie2=IntHP.pH$summary$coefficients[,4],IntHP.pH.ie3=IntHP.pH$summary$coefficients[,5] )
cmpH <- cmP.pH[,c(1,7:10)]
#cmpH <- cmP.pH[,c(1,3:6)]

stargazer(cmP.pH, summary=FALSE, digits=3,type='text')
#==============================================================================================================================================

#============================== pH coef path plot ==================================================
altitude <- data.frame(altitude=seq(-0.0,-0.40,-0.01))

altitude.s <- as.numeric(predict(IntHP.pH$summary$preProc$alt.par , newdata = altitude)[,1])
variables <- c("DEM","Slope","TWI","LongCurv","ChNetBLevel","WEnw")
variables <- variables[order(variables)]

coefs.IntHP <- (cmP.pH[which(as.character(cmP.pH$variable) %in% variables),c(1,7:10)])
coefs.IntP <- (cmP.pH[which(as.character(cmP.pH$variable) %in% variables),c(1,3:6)])

coefs.IntHP <- coefs.IntHP[order(coefs.IntHP$variable),]
coefs.IntP <- coefs.IntP[order(coefs.IntP$variable),]

coefs.IntHP <- coefs.IntHP[,-1]
coefs.IntP <- coefs.IntP[,-1]


effects.IntHP <- as.numeric()
for(i in 1:dim(coefs.IntHP)[1]){
  effects <- coefs.IntHP[i,1] + coefs.IntHP[i,2]*altitude.s + coefs.IntHP[i,3]*altitude.s^2 + coefs.IntHP[i,4]*altitude.s^3
  effects <- data.frame(depth = altitude[,1], var = effects)
  effects.IntHP <- rbind(effects.IntHP, effects)
}  

effects.IntHP$Variables <- factor(rep(variables, each = length(altitude.s)))
intHPlot <- qplot(var, depth, data = effects.IntHP, geom = c("point", "line"), color = Variables)+theme_bw()+theme(legend.text=element_text(size=12))+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+labs(x="Coefficients",y="Depth")


effects.IntP <- as.numeric()
for(i in 1:dim(coefs.IntP)[1]){
  effects <- coefs.IntP[i,1] + coefs.IntP[i,2]*altitude.s + coefs.IntP[i,3]*altitude.s^2 + coefs.IntP[i,4]*altitude.s^3
  effects <- data.frame(depth = altitude[,1], var = effects)
  effects.IntP <- rbind(effects.IntP, effects)
}  

effects.IntP$Variables <- factor(rep(variables, each = length(altitude.s)))
intPlot <- qplot(var, depth, data = effects.IntP, geom = c("point", "line"), color = Variables)+theme_bw()+theme(legend.text=element_text(size=12))+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+labs(x="Coefficients",y="Depth")

pdf("pHIntHPlot.pdf",width=8,height=10)
intHPlot
dev.off()

pdf("pHIntPlot.pdf",width=8,height=10)
intPlot
dev.off()

#=========================== 3D Kriging prediction ============================================================================================
#=========================== The best model is IntP, so we use IntP.pH results! ===========================================================================================

# Try with model=FALSE and krige = FALSE to examine the experimental variograms. If the default model parameters are good enough for initial model parameters, run with mode=TRUE and krige=TRUE. 
pH.kriging.pred <- krige3Dpred(fun=pH.fun, reg.pred=IntP.pH, profs=bor.profs, model = TRUE, krige = FALSE, v.cutoff=40, v.width=3, v.vgm = vgm(0.05, "Gau", 20, 0),sp.vgm = vgm(0.3, "Sph", 1000, 0.1), sp.cutoff=3000, sp.width=640)

#============== Plot Variograms (experimental with model=FALSE and krige=FALSE. If model=TRUE and krige=TRUE the modeled variograms will be ploted) ================
plot(pH.kriging.pred$var1D) # Vertical variogram
plot(pH.kriging.pred$var2D) # Horizontal variogram
plot(pH.kriging.pred$var3D) # 3D variogram

#=== Anistotropy parameter ==============
pH.kriging.pred$anis

#===============================================================================================================================================
pdf("pHIntPVariogram1D.pdf",width=5,height=4)
plot(pH.kriging.pred$var1D)
dev.off()

pdf("pHIntPVariogram2D.pdf",width=5,height=4)
plot(pH.kriging.pred$var2D)
dev.off()

pdf("pHIntPVariogram3D.pdf",width=5,height=4)
plot(pH.kriging.pred$var3D)
dev.off()

#===================== Prediction plots ========================================================================================================


str(pH.kriging.pred$prediction)

prediction<-pH.kriging.pred$prediction[[1]][,"final.prediction"]
names(prediction)<-"pH0.1"

prediction$pH0.2 <- pH.kriging.pred$prediction[[2]]$final.prediction
prediction$pH0.3 <- pH.kriging.pred$prediction[[3]]$final.prediction

prediction$pH0.1 <- pmax(prediction$pH0.1,min(na.omit(bor$pH))/3)
prediction$pH0.2 <- pmax(prediction$pH0.2,min(na.omit(bor$pH))/3)
prediction$pH0.3 <- pmax(prediction$pH0.3,min(na.omit(bor$pH))/3)

pred.stack<-stack(prediction)

gridpix <- as(pred.stack,"SpatialPixelsDataFrame")

zmax <- max(gridpix$pH0.1, gridpix$pH0.2,gridpix$pH0.3)
zmin <- min(gridpix$pH0.1, gridpix$pH0.2,gridpix$pH0.3)
zmax <- round(zmax) 
zmin <- round(zmin) 

ramp <- seq(from = zmin, to = zmax, by = 3)
color.pH <- colorRampPalette(c("dark red","orange","light Yellow"))

ckey <- list(labels=list(cex=2))

p1 <- spplot(gridpix, "pH0.1", asp = 1, at = ramp, col.regions = color.pH,colorkey=FALSE)
p2 <- spplot(gridpix, "pH0.2", asp = 1, at = ramp,  col.regions = color.pH,colorkey=FALSE)
p3 <- spplot(gridpix, "pH0.3", asp = 1, at = ramp, col.regions = color.pH,colorkey=ckey) # last plot contains the legend

#================== plot to files ======================================================
pdf("pHplot0.1.pdf",width=6,height=12)
p1
dev.off()

pdf("pHplot0.2.pdf",width=6,height=12)
p2
dev.off()


pdf("pHplot0.3.pdf",width=7,height=12)
p3
dev.off()
#========================================================================================

