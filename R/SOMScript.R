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
bor.profs <- bor[,c("ID","x","y","Soil.Type","Top","Bottom","As","pH","Co","SOM")]
bor.profs$logSOM <- log(bor.profs$SOM)
bor.profs$logAs <- log(bor.profs$As)
bor.profs$mid <- (bor.profs$Top+bor.profs$Bottom)/2
depths(bor.profs) <- ID ~ Top + Bottom
site(bor.profs) <- ~ Soil.Type + x + y
coordinates(bor.profs) <- ~ x+y
proj4string(bor.profs) <- CRS(utm)
#=====================================================================================

#====================== formulas ================================================================
As.fun <- as.formula(paste("SOM ~", paste(c(CovAbb,"altitude"), collapse="+")))
SOM.fun <- as.formula(paste("log(SOM) ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"altitude"), collapse="+")))
pH.fun <- as.formula(paste("pH ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"altitude"), collapse="+")))

#================================================================================================

#=================================== plot stratified fold - Figure 2 ============================================
rdat <- bor
rdat <- plyr::rename(rdat, replace=c("x" = "longitude", "y" = "latitude"))
rdat <- rdat[complete.cases(rdat[,c("ID","longitude","latitude","altitude","SOM")]),c("ID","longitude","latitude","hdepth","altitude","SOM")] 

coordinates(rdat)<-~longitude+latitude
proj4string(rdat) <- CRS(utm)
rdat <- spTransform(rdat, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
rdat <- as.data.frame(rdat)
head(rdat)

rdat <- rdat[complete.cases(rdat$SOM),]

rdat.folds <- stratfold3d(targetVar="SOM",regdat=rdat,folds=5,cent=3,seed=666,dimensions="3D",IDs=TRUE,sum=TRUE)
plotfolds(rdat.folds,targetvar="SOM")

#========== Plot to file ========================
pdf("SOMfolds.pdf",width=10,height=12)
plotfolds(rdat.folds,targetvar="SOM")
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

#=================================== SOM =======================================================================

#==============================================================================================================
#===================== Nested Cross-validation of trend =======================================================
#==============================================================================================================

#================ Without interactions ========================================================================
SOM.BaseL.ncv <- penint3DP(fun=SOM.fun, profs = bor.profs, seed=444, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="linear")
SOM.BaseP.ncv <- penint3DP(fun=SOM.fun, profs = bor.profs, seed=444, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="poly")

#================ With interactions and without hierarchy =====================================================
SOM.IntL.ncv <- penint3DP(fun=SOM.fun, profs = bor.profs, seed=444, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="linear")
SOM.IntP.ncv <- penint3DP(fun=SOM.fun, profs = bor.profs, seed=444, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")

#================ With interactions and hierarchy (takes more than 40 min!!! ) =============================================================
SOM.IntHL.ncv <- penint3DP(fun=SOM.fun, profs = bor.profs, seed=444, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="linear")
SOM.IntHP.ncv <- penint3DP(fun=SOM.fun, profs = bor.profs, seed=444, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="poly")
#================================
SOM.fun <- as.formula(paste("log(SOM) ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"altitude"), collapse="+")))
logSOM.fun <- as.formula(paste("logSOM ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"altitude"), collapse="+")))

SOM.IntP.ncv <- penint3DP(fun=SOM.fun, profs = bor.profs, seed=444, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
logSOM.IntP.ncv <- penint3DP(fun=logSOM.fun, profs = bor.profs, seed=444, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")



Int<-SOM.IntP.ncv$test.results
logInt<-logSOM.IntP.ncv$test.results



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
sd.res <- rep(NA,length(dseq)-1)
log.sd.res <- rep(NA,length(dseq)-1)
ni <- rep(NA,length(dseq)-1)
log.ni <- rep(NA,length(dseq)-1)

for(i in 1:length(dseq)-1){
  d <- dseq[i+1]
  #R2.res[i] <- caret::R2(pred=Base[Base$altitude < d, "predicted"] ,obs=Base[Base$altitude < d, "observed"])
  sd.res[i] <- sd(x=(Int[Int$altitude <= d & Int$altitude > dseq[i] , "observed"]-Int[Int$altitude <= d & Int$altitude > dseq[i], "predicted"]))
  ni[i] <- length((Int[Int$altitude <= d & Int$altitude > dseq[i], "observed"]-Int[Int$altitude <= d & Int$altitude > dseq[i], "predicted"]))
}

for(i in 1:length(dseq)-1){
  d <- dseq[i+1]
  #R2.res[i] <- caret::R2(pred=Int[Int$altitude < d, "predicted"] ,obs=Int[Int$altitude < d, "observed"])
  log.sd.res[i] <- sd(x=(logInt[logInt$altitude <= d & Int$altitude > dseq[i], "observed"]-logInt[logInt$altitude <= d & Int$altitude > dseq[i], "predicted"]))
  log.ni[i] <- length((logInt[logInt$altitude <= d & logInt$altitude > dseq[i], "observed"]-logInt[logInt$altitude <= d & logInt$altitude > dseq[i], "predicted"]))
}

int.sd.res <- sd.res
log.int.sd.res <- log.sd.res

res.sd <- data.frame(depth=c("0-10 cm","10-20 cm","20-30 cm","30-40 cm","40-60 cm","60- 100 cm"),
                      sd.res=rev(int.sd.res),
                      sd.log.res=rev(log.int.sd.res),
                      n.obs=rev(ni))

stargazer(res.sd,summary = FALSE,digits=2,type='text')

for(i in 1:length(dseq)-1){
  d <- dseq[i+1]
  R2.res[i] <- caret::R2(pred=IntH[IntH$altitude < d, "predicted"] ,obs=IntH[IntH$altitude < d, "observed"])
  RMSE.res[i] <- RMSE(pred=IntH[IntH$altitude < d, "predicted"] ,obs=IntH[IntH$altitude < d, "observed"])
}


#Base results
caret::R2(pred=Base[Base$altitude < d, "predicted"] ,obs=Base[Base$altitude < d, "observed"]) # Rsquared for RK
RMSE(pred=Base[Base$altitude < d, "predicted"] ,obs=Base[Base$altitude < d, "observed"])

#Int results
caret::R2(pred=Int[Int$altitude < d, "predicted"] ,obs=Int[Int$altitude < d, "observed"]) # Rsquared for RK
RMSE(pred=Int[Int$altitude < d, "predicted"] ,obs=Int[Int$altitude < d, "observed"])

caret::R2(pred=logInt[logInt$altitude < d, "predicted"] ,obs=logInt[logInt$altitude < d, "observed"]) # Rsquared for RK
RMSE(pred=logInt[logInt$altitude < d, "predicted"] ,obs=logInt[logInt$altitude < d, "observed"])

#IntH results
caret::R2(pred=IntH[IntH$altitude < d, "predicted"] ,obs=IntH[IntH$altitude < d, "observed"]) # Rsquared for RK
RMSE(pred=IntH[IntH$altitude < d, "predicted"] ,obs=IntH[IntH$altitude < d, "observed"])

#================================
#=================== Nested 5 fold cross-validation results for trend - Table 5 ===================================================

summary.n5cv<-rbind(SOM.BaseL.ncv = SOM.BaseL.ncv$measure[6,], SOM.BaseP.ncv = SOM.BaseP.ncv$measure[6,], SOM.IntL.ncv = SOM.IntL.ncv$measure[6,], SOM.IntP.ncv = SOM.IntP.ncv$measure[6,], SOM.IntHL.ncv = SOM.IntHL.ncv$measure[6,], SOM.IntHP.ncv = SOM.IntHP.ncv$measure[6,])
summary.n5cv[,1] <- rownames(summary.n5cv)
rownames(summary.n5cv)<-NULL
names(summary.n5cv)<-c("Model","RMSE","R squared")
stargazer(summary.n5cv,summary=FALSE,digits=2,type = "text")

latex.table.by(summary.n5cv, num.by.vars = 2)

#==================================================================================================================================
#============================================= krige3D.ncv ========================================================================

# Try with model=FALSE and krige = FALSE to examine the experimental variograms. If the default model parameters are good enough for initial model parameters, run with mode=TRUE and krige=TRUE. 

SOM.krige.res <- krige3D.ncv(SOM.fun,reg.ncv = SOM.IntP.ncv,profs=bor.profs,cogrids = gridmaps.sm2D, model=TRUE,krige=TRUE,v.cutoff = 40,v.width = 3,sp.cutoff = 3000,sp.width = 420, v.vgm = vgm(1.5, "Gau", 10, 0),sp.vgm = vgm(15, "Sph", 2000, 5))

#============== Plot Variograms (experimental with model=FALSE and krige=FALSE. If model=TRUE and krige=TRUE the modeled variograms will be ploted) ================
multiplot(plotlist = SOM.krige.res$var1D,cols=2)
multiplot(plotlist = SOM.krige.res$var2D,cols=2)
multiplot(plotlist = SOM.krige.res$var3D,cols=2)

SOM.final.results <- do.call(rbind,SOM.krige.res$test.results)

#============== Results of nested cv for the best trend model ==================================================
caret::R2(pred=SOM.final.results$predicted,obs=SOM.final.results$observed, formula="traditional") # Rsquared for trend
RMSE(pred=SOM.final.results$predicted,obs=SOM.final.results$observed) # RMSE for trend
#===============================================================================================================

#============== Results of nested cv for the 3D regression kriging ==================================================
caret::R2(pred=SOM.final.results$final.predicted,obs=SOM.final.results$observed, formula="traditional") # Rsquared for RK
RMSE(pred=SOM.final.results$final.predicted,obs=SOM.final.results$observed) # RMSE for RK
#=================================================================================================================
#============================================================================================================================================

#============================================================================================================================================
#============== PREDICTION ==================================================================================================================
#============================================================================================================================================


#=================== SOM =====================================================================================================================
#====== Without interactions ================================================================================================================
BaseL.SOM <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8, seed = 444)
BaseP.SOM <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8, seed = 444)

#====== With interactions but without hierarchy ================================================================================================================
IntL.SOM <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8, seed = 444)
IntP.SOM <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8, seed = 444)
# IntP.SOM model is the best. This is the reason why the pred=TRUE. If pred=TRUE the prediction over grid will be performed. If pred=FALSE, only the prediction on observed points will be performed

#====== With interactions honoring the hierarchy ================================================================================================================
#==================================== Require free memory and takes more than 40 min!!! =====================================================
IntHL.SOM <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8, seed = 444)
IntHP.SOM <- predint3DP(fun=SOM.fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred = FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8, seed = 444)
#============================================================================================================================================

#============ Plot glmnet model ======================
pdf("LassoCoef.pdf",width=8,height=6)
plot(IntL.SOM$summary$model$glmnet.fit,xvar = "lambda", label = TRUE)
abline(v=log(c(IntL.SOM$summary$model$lambda.min)),lty=2)
dev.off()

plot(IntL.SOM$summary$model)
#====================================================


SOM.pred <- list(BaseL.SOM=BaseL.SOM$summary$pred, BaseP.SOM=BaseP.SOM$summary$pred, IntL.SOM=IntL.SOM$summary$pred, IntP.SOM=IntP.SOM$summary$pred, IntHL.SOM=IntHL.SOM$summary$pred , IntHP.SOM=IntHP.SOM$summary$pred)
SOM.coef <- list(BaseL.SOM=BaseL.SOM$summary$coefficients, BaseP.SOM=BaseP.SOM$summary$coefficients, IntL.SOM=IntL.SOM$summary$coefficients, IntP.SOM=IntP.SOM$summary$coefficients, IntHL.SOM=IntHL.SOM$summary$coefficients ,IntHP.SOM=IntHP.SOM$summary$coefficients)

SOM.pred <- lapply(SOM.pred, function(x) {names(x)<-c("ID", "longitude", "latitude","hdepth", "altitude","observed", "predicted","residual" ); return(x)})

#============================= Coefficients tables for SOM models ===========================================================================================
ll <- length(IntL.SOM$summary$coefficients)
pp <- length(IntHL.SOM$summary$coefficients[,1])+1


#============================= Coefficients for models with linear depth function ==============================================================
cmL <- data.frame(variable=IntHL.SOM$summary$coefficients[,1], BaseL.SOM.me=BaseL.SOM$summary$coefficients[2:pp], IntL.SOM.me=IntL.SOM$summary$coefficients[2:pp],IntL.SOM.ie=c(IntL.SOM$summary$coefficients[(pp+1):ll],0),IntHL.SOM.me=IntHL.SOM$summary$coefficients[,2],IntHL.SOM.ie=IntHL.SOM$summary$coefficients[,3] )


#============================= Coefficients for models with polynomial depth function ===========================================================
l <- length(IntP.SOM$summary$coefficients)
p <- length(IntHP.SOM$summary$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.SOM <- data.frame(variable=IntHP.SOM$summary$coefficients[,1], BaseP.SOM.me=BaseP.SOM$summary$coefficients[2:p], IntP.SOM.me=IntP.SOM$summary$coefficients[2:p],IntP.SOM.ie1=c(IntP.SOM$summary$coefficients[(p+1):l][i1],0,0,0),IntP.SOM.ie2=c(IntP.SOM$summary$coefficients[(p+1):l][i2],0,0,0),IntP.SOM.ie3=c(IntP.SOM$summary$coefficients[(p+1):l][i3],0,0,0),IntHP.SOM.me=IntHP.SOM$summary$coefficients[,2],IntHP.SOM.ie1=IntHP.SOM$summary$coefficients[,3],IntHP.SOM.ie2=IntHP.SOM$summary$coefficients[,4],IntHP.SOM.ie3=IntHP.SOM$summary$coefficients[,5] )
cmSOM <- cmP.SOM[,c(1,7:10)]
#cmSOM <- cmP.SOM[,c(1,3:6)]

stargazer(cmP.SOM, summary=FALSE, digits=3,type='text')
#==============================================================================================================================================

#============================== SOM coef path plot ==================================================
altitude <- data.frame(altitude=seq(-0.0,-0.40,-0.01))

altitude.s <- as.numeric(predict(IntHP.SOM$summary$preProc$alt.par , newdata = altitude)[,1])
variables <- c("DEM","Slope","TWI","LongCurv","ChNetBLevel","WEnw")
variables <- variables[order(variables)]

coefs.IntHP <- (cmP.SOM[which(as.character(cmP.SOM$variable) %in% variables),c(1,7:10)])
coefs.IntP <- (cmP.SOM[which(as.character(cmP.SOM$variable) %in% variables),c(1,3:6)])

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

pdf("SOMIntHPlot.pdf",width=8,height=10)
intHPlot
dev.off()

pdf("SOMIntPlot.pdf",width=8,height=10)
intPlot
dev.off()

#=========================== 3D Kriging prediction ============================================================================================
#=========================== The best model is IntP, so we use IntP.SOM results! ===========================================================================================

# Try with model=FALSE and krige = FALSE to examine the experimental variograms. If the default model parameters are good enough for initial model parameters, run with mode=TRUE and krige=TRUE. 
SOM.kriging.pred <- krige3Dpred(fun=SOM.fun, reg.pred=IntP.SOM, profs=bor.profs, model = TRUE, krige = TRUE, v.cutoff=60, v.width=3, v.vgm = vgm(1.5, "Gau", 10, 0),sp.vgm = vgm(15, "Sph", 2000, 5), sp.cutoff=4000, sp.width=480)

#============== Plot Variograms (experimental with model=FALSE and krige=FALSE. If model=TRUE and krige=TRUE the modeled variograms will be ploted) ================
plot(SOM.kriging.pred$var1D) # Vertical variogram
plot(SOM.kriging.pred$var2D) # Horizontal variogram
plot(SOM.kriging.pred$var3D) # 3D variogram

#=== Anistotropy parameter ==============
SOM.kriging.pred$anis

#===============================================================================================================================================


#===================== Prediction plots ========================================================================================================


str(SOM.kriging.pred$prediction)

prediction<-SOM.kriging.pred$prediction[[1]][,"final.prediction"]
names(prediction)<-"SOM0.1"

prediction$SOM0.2 <- SOM.kriging.pred$prediction[[2]]$final.prediction
prediction$SOM0.3 <- SOM.kriging.pred$prediction[[3]]$final.prediction

prediction$SOM0.1 <- pmax(prediction$SOM0.1,min(na.omit(bor$SOM))/3)
prediction$SOM0.2 <- pmax(prediction$SOM0.2,min(na.omit(bor$SOM))/3)
prediction$SOM0.3 <- pmax(prediction$SOM0.3,min(na.omit(bor$SOM))/3)

pred.stack<-stack(prediction)

gridpix <- as(pred.stack,"SpatialPixelsDataFrame")

zmax <- max(gridpix$SOM0.1, gridpix$SOM0.2,gridpix$SOM0.3)
zmin <- min(gridpix$SOM0.1, gridpix$SOM0.2,gridpix$SOM0.3)
zmax <- round(zmax) 
zmin <- round(zmin) 

ramp <- seq(from = zmin, to = zmax, by = 3)
color.SOM <- colorRampPalette(c("dark red","orange","light Yellow"))

ckey <- list(labels=list(cex=2))

p1 <- spplot(gridpix, "SOM0.1", asp = 1, at = ramp, col.regions = color.SOM,colorkey=FALSE)
p2 <- spplot(gridpix, "SOM0.2", asp = 1, at = ramp,  col.regions = color.SOM,colorkey=FALSE)
p3 <- spplot(gridpix, "SOM0.3", asp = 1, at = ramp, col.regions = color.SOM,colorkey=ckey) # last plot contains the legend

#================== plot to files ======================================================
pdf("SOMplot0.1.pdf",width=6,height=12)
p1
dev.off()

pdf("SOMplot0.2.pdf",width=6,height=12)
p2
dev.off()


pdf("SOMplot0.3.pdf",width=7,height=12)
p3
dev.off()
#========================================================================================

