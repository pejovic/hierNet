
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

#load("PredTables1452016.RData")

gk_7 <- "+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m"
utm <- "+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m"

load(paste(getwd(),"inst","extdata","gridmaps.RDa",sep="/"))
gridmaps.sm2D$CD <- exp(-gridmaps.sm2D$DirDif)
gridmaps.sm2D$DD <- as.numeric(exp(-scale(gridmaps.sm2D$Dist,center=FALSE)))

sm2D.lst<-names(gridmaps.sm2D)
sm2D.lst <- sm2D.lst[ -which(sm2D.lst %in% c("DirDif","Dist","AnalyticalHills","LSFactor","RelSlopePosition","VelleyDepth","optional")) ]

names(gridmaps.sm2D) <- c("AHils","Aspect","CatchArea","ChNetBLevel","ConvInd","CrSectCurv","DEM","DirDif","Dist","ES","LongCurv","LSFactor","NegOp","PosOp","RelSlopePosition","Slope","TWI","VelleyDepth","VDistChNet","WEeast","WEnw","clc","SoilType","CD","DD")

#================== Covariates table =================================================

CovNames <- c("Digital Elevation Model", "Aspect", "Slope","Topographic Wetness Index", "Convergence Index" ,"Cross Sectional Curvature", "Longitudinal Curvature", "Channel Network Base Level" ,"Vertical Distance to Channel Network", "Negative Openness","Positive Openness", "Wind Effect (East)","Wind Effect (North-West)","Down-wind Dilution", "Cross-wind Dilution" ,"Corine Land Cover 2006", "Soil Type")
CovAbb <- c("DEM","Aspect","Slope","TWI","ConvInd","CrSectCurv","LongCurv","ChNetBLevel","VDistChNet","NegOp", "PosOp","WEeast","WEnw","DD","CD","clc","SoilType")

#covs<-data.frame(Name=CovNames,
#                 Abbrevation=CovAbb,
#                Type=c("C","C","C","C","C","C","C","C","C","C","C","C","C","C","C","F","F")
#                 )

#stargazer(covs,summary=FALSE,type="latex",digits=2)


#gridmaps <- gridmaps.sm2D[,sm2D.lst]
#save(gridmaps,file="C:/Users/Milutin/Dropbox/ES3D/Data and Scripts/gridmaps.rda")


bor <- join(read.csv(paste(getwd(),"inst","extdata","Profili_sredjeno_csv.csv",sep="/")), read.csv(paste(getwd(),"inst","extdata","Koordinate_csv.csv",sep="/")), type="inner")

levels(bor$Soil.Type) <- c("Dystric Cambisol","Dystric Leptosol","Dystric Regosol","Eutric Cambisol","Eutric Leptosol","Eutric Regosol","Calcaric Cambisol","Mollic Leptosol","Colluvium","Luvisol","Planosol","Vertisol")


levels(bor$Soil.Type) <- c("Dystric Cambisol","Dystric Leptosol","Dystric Regosol","Eutric Cambisol","Eutric Leptosol","Eutric Regosol","Calcaric Cambisol","Mollic Leptosol","Colluvium","Luvisol","Planosol","Vertisol")

bor$hdepth<-bor$Bottom-bor$Top
bor$altitude <- - (bor$Top / 100 + (( bor$Bottom - bor$Top ) / 2) / 100)
bor <- bor[, - c(7:12,14,15,16,17 )]
names(bor) <- c("Soil.Type","ID","Horizont","Top" , "Bottom","pH","SOM","As","Co","x","y","hdepth","altitude")

bor[which(bor$ID %in% c(66,129,14,51,69,130,164,165,166)),"As"] <- NA
#bor <- bor[-which(bor$ID %in% c(66,129,14,51,69,130,164,165,166)),]

#========================= Creating Soil Profile Collections ====================================
bor.profs <- bor[,c("ID","x","y","Soil.Type","Top","Bottom","SOM","pH","Co","As")]
depths(bor.profs) <- ID ~ Top + Bottom
site(bor.profs) <- ~ Soil.Type + x + y
coordinates(bor.profs) <- ~ x+y
proj4string(bor.profs) <- CRS(utm)
bor.geo<-as.geosamples(bor.profs)


a<-mpspline(bor.profs,"SOM",vlow=1)
str(a)
ab<-ddply(bor.profs@horizons,.(ID), summarize, variation=var(SOM))
singleValueProfile<-ab[is.na(ab$variation),"ID"]

std<-cbind(a$idcol,a$var.std)
names(std)[1]<-"ID"
str(std)
br.meas <- apply(std[,c(2:7)],2,function(x) sum(!is.na(x))) # zakljucak...nema dovoljno merenja za dubine preko 60cm

stargazer(br.meas, summary=FALSE, type="latex")

#====================== formulas ================================================================
As.fun <- as.formula(paste("As ~", paste(c(CovAbb,"altitude"), collapse="+")))
SOM.fun <- as.formula(paste("SOM ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"altitude"), collapse="+")))
pH.fun <- as.formula(paste("pH ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"altitude"), collapse="+")))

#================================================================================================

#================== test for stratfold3d and penint3D ============================================

source(paste(getwd(),"R","stratFold3D.R",sep="/"))
#source(paste(getwd(),"R","penint3D_def.R",sep="/"))
source(paste(getwd(),"R","plotfolds.R",sep="/"))
#source(paste(getwd(),"R","predint3D.R",sep="/"))
source(paste(getwd(),"R","penint3D_defP.R",sep="/"))
source(paste(getwd(),"R","predint3DP.R",sep="/"))

fun <- As.fun

fun=As.fun; profs = bor.profs; seed=443; cogrids = gridmaps.sm2D; hier = FALSE; int=TRUE; depth.fun="poly"
lambda=seq(0,5,0.1);deg=3;fold=5;cent=3;preProc=TRUE

#=================================== plot stratified fold ============================================
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
plotfolds(rdat.folds,targetvar="SOM")

pdf("SOMfolds.pdf",width=10,height=12)
plotfolds(rdat.folds,targetvar="SOM")
dev.off()



stargazer(do.call(rbind, rdat.folds$`SOM summary`), summary=FALSE, digits=2, type="latex")
stargazer(do.call(rbind,as.list(rdat.folds$`altitude summary`)[[1]]), summary=FALSE, digits=2, type="text")

#============= eps file ====================================

library(ggplot2)
postscript("SOMfolds.pdf",horizontal=FALSE, paper="special",height=11,width=8)
plotfolds(rdat.folds,targetvar="SOM")
ggsave("SOMfolds.pdf", height=11,width=8)

dev.off()

#============================================================

#==============================================================================================================

BaseLcv <- penint3DP(fun=fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="linear")
BaseLcv$measure
summary(BaseLcv$measure[1:5,])


IntLcv <- penint3DP(fun=fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="linear")
IntLcv$measure
summary(IntLcv$measure[1:5,])


IntHLcv <- penint3DP(fun=fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = TRUE, int=FALSE, depth.fun="linear")
IntHLcv$measure
summary(IntHLcv$measure[1:5,])


#============= Poly ========================
BasePcv <- penint3DP(fun=fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="poly")
BasePcv$measure
summary(BasePcv$measure[1:5,])


IntPcv <- penint3DP(fun=fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
IntPcv$measure
summary(IntPcv$measure[1:5,])


IntPHcv <- penint3DP(fun=fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="poly")
IntPHcv$measure
summary(IntPHcv$measure[1:5,])

#=================== Nested 5 fold cross-validation results =======================

summary.n5cv<-rbind(BaseLcv = BaseLcv$measure[6,], BasePcv = BasePcv$measure[6,], IntLcv = IntLcv$measure[6,], IntPcv = IntPcv$measure[6,], IntHLcv = IntHLcv$measure[6,], IntHPcv = IntHPcv$measure[6,])
summary.n5cv[,1] <- rownames(summary.n5cv)
rownames(summary.n5cv)<-NULL
names(summary.n5cv)<-c("Model","RMSE","R squared")
stargazer(summary.n5cv,summary=FALSE,digits=2)

latex.table.by(summary.n5cv, num.by.vars = 2)

#=================== predint3D ====================================================


fun <- SOM.fun

BaseL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
BaseP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntHL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntHP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

SOM.pred <- list(BaseL=BaseL$summary$pred, BaseP=BaseP$summary$pred, IntL=IntL$summary$pred, IntP=IntP$summary$pred, IntHL=IntHL$summary$pred , IntHP=IntHP$summary$pred)
SOM.coef <- list(BaseL=BaseL$summary$coefficients, BaseP=BaseP$summary$coefficients, IntL=IntL$summary$coefficients, IntP=IntP$summary$coefficients, IntHL=IntHL$summary$coefficients ,IntHP=IntHP$summary$coefficients)

SOM.pred <- lapply(SOM.pred, function(x) {names(x)<-c("ID", "longitude", "latitude","hdepth", "altitude","observed", "predicted","residual" ); return(x)})
str(SOM.pred)

#============================= Coefficients matrix ==========================================
ll <- length(IntL$summary$coefficients)
pp <- length(IntHL$summary$coefficients[,1])+1


cmL <- data.frame(variable=IntHL$summary$coefficients[,1], BaseL.me=BaseL$summary$coefficients[2:pp], IntL.me=IntL$summary$coefficients[2:pp],IntL.ie=c(IntL$summary$coefficients[(pp+1):ll],0),IntHL.me=IntHL$summary$coefficients[,2],IntHL.ie=IntHL$summary$coefficients[,3] )


l <- length(IntP$summary$coefficients)
p <- length(IntHP$summary$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)


cmP.SOM <- data.frame(variable=IntHP$summary$coefficients[,1], BaseP.me=BaseP$summary$coefficients[2:p], IntP.me=IntP$summary$coefficients[2:p],IntP.ie1=c(IntP$summary$coefficients[(p+1):l][i1],0,0,0),IntP.ie2=c(IntP$summary$coefficients[(p+1):l][i2],0,0,0),IntP.ie3=c(IntP$summary$coefficients[(p+1):l][i3],0,0,0),IntHP.me=IntHP$summary$coefficients[,2],IntHP.ie1=IntHP$summary$coefficients[,3],IntHP.ie2=IntHP$summary$coefficients[,4],IntHP.ie3=IntHP$summary$coefficients[,5] )
cmSOM <- cmP.SOM[,c(1,7:10)]
#cmSOM <- cmP.SOM[,c(1,3:6)]

#=================== predint3D ====================================================
fun <- pH.fun

BaseL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
BaseP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntHL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntHP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

pH.pred <- list(BaseL=BaseL$summary$pred, BaseP=BaseP$summary$pred, IntL=IntL$summary$pred, IntP=IntP$summary$pred, IntHL=IntHL$summary$pred , IntHP=IntHP$summary$pred)
pH.coef <- list(BaseL=BaseL$summary$coefficients, BaseP=BaseP$summary$coefficients, IntL=IntL$summary$coefficients, IntP=IntP$summary$coefficients, IntHL=IntHL$summary$coefficients ,IntHP=IntHP$summary$coefficients)

pH.pred <- lapply(pH.pred, function(x) {names(x)<-c("ID", "longitude", "latitude","hdepth", "altitude","observed", "predicted","residual" ); return(x)})
#============================= Coefficients matrix ==========================================
ll <- length(IntL$summary$coefficients)
pp <- length(IntHL$summary$coefficients[,1])+1


cmL <- data.frame(variable=IntHL$summary$coefficients[,1], BaseL.me=BaseL$summary$coefficients[2:pp], IntL.me=IntL$summary$coefficients[2:pp],IntL.ie=c(IntL$summary$coefficients[(pp+1):ll],0),IntHL.me=IntHL$summary$coefficients[,2],IntHL.ie=IntHL$summary$coefficients[,3] )


l <- length(IntP$summary$coefficients)
p <- length(IntHP$summary$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.pH <- data.frame(variable=IntHP$summary$coefficients[,1], BaseP.me=BaseP$summary$coefficients[2:p], IntP.me=IntP$summary$coefficients[2:p],IntP.ie1=c(IntP$summary$coefficients[(p+1):l][i1],0,0,0),IntP.ie2=c(IntP$summary$coefficients[(p+1):l][i2],0,0,0),IntP.ie3=c(IntP$summary$coefficients[(p+1):l][i3],0,0,0),IntHP.me=IntHP$summary$coefficients[,2],IntHP.ie1=IntHP$summary$coefficients[,3],IntHP.ie2=IntHP$summary$coefficients[,4],IntHP.ie3=IntHP$summary$coefficients[,5] )
cmpH <- cmP.pH[,c(1,7:10)]
#cmpH <- cmP.pH[,c(1,3:6)]
#============================ As ====================

#bor <- bor[-which(bor$ID %in% c(66,129,14,51,69,130,164,165,166)),]

#========================= Creating Soil Profile Collections ====================================
#bor.profs <- bor[,c("ID","x","y","Soil.Type","Top","Bottom","SOM","pH","Co","As")]
#depths(bor.profs) <- ID ~ Top + Bottom
#site(bor.profs) <- ~ Soil.Type + x + y
#coordinates(bor.profs) <- ~ x+y
#proj4string(bor.profs) <- CRS(utm)
#bor.geo<-as.geosamples(bor.profs)

#=================== predint3D ====================================================
fun <- As.fun


BaseL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
BaseP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntHL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntHP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

As.pred <- list(BaseL=BaseL$summary$pred, BaseP=BaseP$summary$pred, IntL=IntL$summary$pred, IntP=IntP$summary$pred, IntHL=IntHL$summary$pred , IntHP=IntHP$summary$pred)
As.coef <- list(BaseL=BaseL$summary$coefficients, BaseP=BaseP$summary$coefficients, IntL=IntL$summary$coefficients, IntP=IntP$summary$coefficients, IntHL=IntHL$summary$coefficients ,IntHP=IntHP$summary$coefficients)

As.pred <- lapply(As.pred, function(x) {names(x)<-c("ID", "longitude", "latitude","hdepth", "altitude","observed", "predicted","residual" ); return(x)})

#============================= Coefficients matrix ==========================================
ll <- length(IntL$summary$coefficients)
pp <- length(IntHL$summary$coefficients[,1])+1


cmL <- data.frame(variable=IntHL$summary$coefficients[,1], BaseL.me=BaseL$summary$coefficients[2:pp], IntL.me=IntL$summary$coefficients[2:pp],IntL.ie=c(IntL$summary$coefficients[(pp+1):ll],0),IntHL.me=IntHL$summary$coefficients[,2],IntHL.ie=IntHL$summary$coefficients[,3] )


l <- length(IntP$summary$coefficients)
p <- length(IntHP$summary$coefficients[,1])+1
i1 <- seq(1,l-p,3)
i2 <- seq(2,l-p,3)
i3 <- seq(3,l-p,3)

cmP.As <- data.frame(variable=IntHP$summary$coefficients[,1], BaseP.me=BaseP$summary$coefficients[2:p], IntP.me=IntP$summary$coefficients[2:p],IntP.ie1=c(IntP$summary$coefficients[(p+1):l][i1],0,0,0),IntP.ie2=c(IntP$summary$coefficients[(p+1):l][i2],0,0,0),IntP.ie3=c(IntP$summary$coefficients[(p+1):l][i3],0,0,0),IntHP.me=IntHP$summary$coefficients[,2],IntHP.ie1=IntHP$summary$coefficients[,3],IntHP.ie2=IntHP$summary$coefficients[,4],IntHP.ie3=IntHP$summary$coefficients[,5] )
cmAs <- cmP.As[,c(1,7:10)]
#cmAs <- cmP.As[,c(1,3:6)]
#======================



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

#============================ Changing coefficients =================================================

altitude <- data.frame(altitude=seq(-0.0,-0.40,-0.01))

#altitude <- data.frame(altitude=c(-0.1,-0.3,-0.4))
altitude.s <- as.numeric(predict(IntHP$summary$preProc$alt.par , newdata = altitude)[,1])
variables <- c("CD","DD","DEM","Slope","VDistChNet","WEnw")#,"TWI")
variables <- variables[order(variables)]

#10.4895187  5.7506602 -0.4258973000
#cmP.As[cmP.As$variable=="VDistChNet","IntHP.ie1"] <- 5.7506602

coefs.IntHP <- (cmP.As[which(as.character(cmP.As$variable) %in% variables),c(1,7:10)])
coefs.IntP <- (cmP.As[which(as.character(cmP.As$variable) %in% variables),c(1,3:6)])

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

pdf("AsIntHPlot.pdf",width=8,height=10)
intHPlot
dev.off()

pdf("AsIntPlot.pdf",width=8,height=10)
intPlot
dev.off()

#============================== SOM coef path plot ==================================================
altitude <- data.frame(altitude=seq(-0.0,-0.40,-0.01))

altitude.s <- as.numeric(predict(IntHP$summary$preProc$alt.par , newdata = altitude)[,1])
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


#======================------ Changing coefficients table ============================================
var.names <- as.character(intHCoefs$variable)
var.names <- c(var.names[1:20],"CMca","CMdy","LPdy","RGdy","CMeu","LPeu","LPmo","VR","d","d2","d3")


As01 <- intHCoefs$As.me+intHCoefs$As.ie1*(-0.1)+intHCoefs$As.ie2*(-0.1)^2+intHCoefs$As.ie3*(-0.1)^3
As02 <- intHCoefs$As.me+intHCoefs$As.ie1*(-0.2)+intHCoefs$As.ie2*(-0.2)^2+intHCoefs$As.ie3*(-0.2)^3
As03 <- intHCoefs$As.me+intHCoefs$As.ie1*(-0.3)+intHCoefs$As.ie2*(-0.3)^2+intHCoefs$As.ie3*(-0.3)^3

As.depth.coef <- data.frame(cbind(As01,As02,As03))
names(As.depth.coef) <- var.names

SOM01 <- intHCoefs$SOM.me+intHCoefs$SOM.ie1*(-0.1)+intHCoefs$SOM.ie2*(-0.1)^2+intHCoefs$SOM.ie3*(-0.1)^3
SOM02 <- intHCoefs$SOM.me+intHCoefs$SOM.ie1*(-0.2)+intHCoefs$SOM.ie2*(-0.2)^2+intHCoefs$SOM.ie3*(-0.2)^3
SOM03 <- intHCoefs$SOM.me+intHCoefs$SOM.ie1*(-0.3)+intHCoefs$SOM.ie2*(-0.3)^2+intHCoefs$SOM.ie3*(-0.3)^3

SOM.depth.coef <- data.frame(cbind(SOM01,SOM02,SOM03))
names(SOM.depth.coef) <- var.names


pH01 <- intHCoefs$pH.me+intHCoefs$pH.ie1*(-0.1)+intHCoefs$pH.ie2*(-0.1)^2+intHCoefs$pH.ie3*(-0.1)^3
pH02 <- intHCoefs$pH.me+intHCoefs$pH.ie1*(-0.2)+intHCoefs$pH.ie2*(-0.2)^2+intHCoefs$pH.ie3*(-0.2)^3
pH03 <- intHCoefs$pH.me+intHCoefs$pH.ie1*(-0.3)+intHCoefs$pH.ie2*(-0.3)^2+intHCoefs$pH.ie3*(-0.3)^3

pH.depth.coef <- data.frame(cbind(pH01,pH02,pH03))
names(pH.depth.coef) <- var.names

depth.coef <- cbind(var.names,As.depth.coef,SOM.depth.coef,pH.depth.coef)

stargazer(depth.coef[1:28,],summary=FALSE, digits=3, type="latex")


#=====================================================================================================

intHCoefs$As.ie2*(-0.2)

#====================================================================================================
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
IntP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

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


#=========================== Profile plots ===================================================

FFL.p <- SOM.pred$BaseL
FFP.p <- SOM.pred$BaseP
TFL.p <- SOM.pred$IntL
TFP.p <- SOM.pred$IntP
TTL.p <- SOM.pred$IntHP
TTP.p <- SOM.pred$IntHP

res.profs<-data.frame(ID=TFL.p$ID,BaseL=FFL.p$residual,BaseP=FFP.p$residual,IntL=TFL.p$residual,IntP=TFP.p$residual,IntHL=TTL.p$residual,IntHP=TTP.p$residual)

borind<-ddply(bor,.(ID))

ind<-which(res.profs$ID %in% borind[,"ID"])

res.profs<-cbind(res.profs,borind[ind,c("Top","Bottom","Soil.Type","x","y","altitude")])

res.list<-list(sites=ddply(res.profs,.(ID),summarize,y=y[1],x=x[1],Soil.Type=Soil.Type[1]),horizons=res.profs[,c("ID","Top","Bottom","BaseL","BaseP","IntL","IntP","IntHL","IntHP")])
res.profs<-join(res.list$sites,res.list$horizons, type="inner")

res.profs.data<-res.profs

depths(res.profs) <- ID ~ Top + Bottom
site(res.profs) <- ~x+y
coordinates(res.profs)<-~x+y
proj4string(res.profs)<- CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m")

# again, this time along 1-cm slices, computing quantiles
agg <- slab(res.profs, fm= ~ BaseL + BaseP + IntL + IntP + IntHL + IntHP,slab.structure=seq(0,70,5)) #c(0,5,15,30,40,60,80)

## see ?slab for details on the default aggregate function
head(agg)


#agg$contributing_fraction <- do.call(c,c(alt.summary[,2:7]))

res.plot<-xyplot(top ~ p.q50 | variable, data=agg, ylab='Depth',
       xlab='median bounded by 25th and 75th percentiles',
       lower=agg$p.q25, upper=agg$p.q75, ylim=c(60,-2),
       panel=panel.depth_function,
       alpha=0.25, sync.colors=TRUE,
       par.settings=list(superpose.line=list(col='RoyalBlue', lwd=2)),
       prepanel=prepanel.depth_function,
       cf=agg$contributing_fraction, cf.col='black', cf.interval=5, 
       layout=c(6,1), strip=strip.custom(bg=grey(0.6)),
       scales=list(x=list(tick.number=4, alternating=3, relation='free'))
)

pdf("AsresPlots.pdf",width=12,height=8)
plot(res.plot) # Make plot
dev.off()


#========================= Residual comparison ===============================================
res.profs.data$altitude <- - (res.profs.data$Top / 100 + (( res.profs.data$Bottom - res.profs.data$Top ) / 2) / 100)

res.profs.data<-res.profs.data[order(-res.profs.data$altitude),]

Alt_quant<-classIntervals(res.profs.data$altitude,style="fixed",fixedBreaks=c(0,-0.05,-0.15,-0.30,-0.40,-0.60,min(res.profs.data$altitude)))
#Kreiranje faktora prema clasifikaciji koja je uradjena....
Alt_qclass<- cut(Alt_quant$var, breaks=Alt_quant$brks, #labels =c("-0.5","-0.5","3","4",),
                 include.lowest = TRUE, right = TRUE, dig.lab = 4,
                 ordered_result = FALSE)

res.profs.data$Alt.class<-Alt_qclass

sumfun <- function(x, ...){sqrt(sum(x^2)/(length(x)))}


alt.summary <- summaryBy(BaseL + BaseP + IntL + IntP + IntHL + IntHP~Alt.class, data=res.profs.data, FUN=mean)
alt.summary <- alt.summary[order(-as.numeric(alt.summary$Alt.class)),]
names(alt.summary) <- c("Depth","BaseL","BaseP","IntL","IntP","IntHL","IntHP")
row.names(alt.summary)<-NULL
alt.summary[,1]<-c("0-5 cm","5-15 cm","15-30 cm","30-40 cm","40-60 cm","60-125 cm")
stargazer(alt.summary,summary=FALSE,digits=2,type="latex")

#do.call(c,c(alt.summary[,2:7]))


ddply(res.profs.data[,c(1,7:12,14)],.(ID),summarize, r=max(BaseL)-min(BaseL))


################## Figure 2 ######################################################################################
# again, this time along 1-cm slices, computing quantiles
mean.and.sd <- function(values) {
  narm.values<-values[!is.na(values)]
  m <- mean(values, na.rm=TRUE)
  s <- sd(values, na.rm=TRUE)
  s.mean<-s/sqrt(length(narm.values))
  mini <- min(values, na.rm=TRUE)
  maxi <- max(values, na.rm=TRUE)
  q25<-quantile(values, probs = 0.25, na.rm=TRUE)
  q75<-quantile(values, probs = 0.75, na.rm=TRUE)
  cvar<-s/m
  med<-median(values,na.rm=TRUE)
  inqr <- IQR(values,na.rm=TRUE)
  n<-length(narm.values)
  res <- c(min=mini,q1=q25,mean=m, mean.sd=s.mean ,median=med, q3=q75 ,max=maxi, IQR=inqr,sd=s ,CV=cvar,obs=n)
  return(res)
}


mean.and.sd2 <- function(values) {
  narm.values<-values[!is.na(values)]
  m <- mean(values, na.rm=TRUE)
  s <- sd(values, na.rm=TRUE)
  s.mean<-s/sqrt(length(narm.values))
  n<-length(narm.values)
  res <- c(mean=m, mean.sd=s.mean ,sd=s,obs=n)
  return(res)
}


agg <- slab(bor.profs, fm= ~ As+SOM+pH, slab.structure=seq(0,70,5),slab.fun=mean.and.sd2)

## see ?slab for details on the default aggregate function
head(agg)

Figure2<-xyplot(top ~ p.q50 | variable, data=agg, ylab='Depth',
       xlab='median bounded by 25th and 75th percentiles',
       lower=agg$p.q25, upper=agg$p.q75, ylim=c(70,-2),
       panel=panel.depth_function,
       alpha=0.25, sync.colors=TRUE,
       par.settings=list(superpose.line=list(col='RoyalBlue', lwd=2)),
       prepanel=prepanel.depth_function,
       cf=agg$contributing_fraction, cf.col='black', cf.interval=5, 
       layout=c(3,1), strip=strip.custom(bg=grey(0.8)),
       scales=list(x=list(tick.number=4, alternating=3, relation='free'))
)

class(Figure2)

pdf("AsSOMpH.pdf",width=10,height=8)
plot(Figure2) # Make plot
dev.off()



#======================= Nested CV (for regression part and for RK) ================================================

#======================= As nested kros validacija ================================================================
As.IntP.cv <- penint3DP(fun=As.fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
As.IntP.cv$measure
summary(As.IntP.cv$measure[1:5,])


As.IntPH.cv <- penint3DP(fun=As.fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="poly")
As.IntPH.cv$measure
summary(As.IntPH.cv$measure[1:5,])

for(i in 1:5){

TTP.p <- As.IntP.cv$train.results[[i]]
TTP.p$residual <- TTP.p$observed-TTP.p$predicted

res.profs<-data.frame(ID=TTP.p$ID,IntHP=TTP.p$residual)

borind<-ddply(bor,.(ID))

borind<-borind[complete.cases(borind$As),] # Ovde treba zameniti...

"%ni%" <- Negate("%in%")

ind <- which(borind[,"ID"] %in% unique(res.profs$ID))

res.profs <- cbind(res.profs, borind[ind,c("Top","Bottom","Soil.Type","x","y","altitude")])

res.list<-list(sites=ddply(res.profs,.(ID),summarize,y=y[1],x=x[1],Soil.Type=Soil.Type[1]),horizons=res.profs[,c("ID","Top","Bottom","IntHP")])
res.profs<-join(res.list$sites,res.list$horizons, type="inner")

res.profs.data<-res.profs

depths(res.profs) <- ID ~ Top + Bottom
site(res.profs) <- ~x+y
coordinates(res.profs)<-~x+y
proj4string(res.profs)<- CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m")


res.data <- TTP.p

spline.profs <- mpspline(res.profs, "IntHP",d = t(c(0,5,15,30,60,80)))
spline.res <- spline.profs$var.1cm

cm <- rep(1:80,dim(spline.res)[2])

a <-ddply(res.data,.(ID),summarize ,longitude=longitude[1],latitude=latitude[1])

b <- a[rep(seq_len(nrow(a)), each=80),]

bm <- cbind(b,residual=as.numeric(spline.res),altitude=cm)

bm <- na.omit(bm)

vc <- variogram(residual~1, data=bm, locations=~longitude+latitude+altitude,cutoff=50,width=3)
vc
plot(vc)

vc.fit = fit.variogram(vc, vgm(250, "Gau", 30, 5)) #As
plot(vc, vc.fit)

res.pp <- ddply(res.data, .(ID),summarize,x=longitude[1],y=latitude[1],altitude=altitude[1],residual=residual[1],observed=observed[1])
coordinates(res.pp) <- ~x+y+altitude
proj4string(res.pp)<- CRS(gk_7)

bubble(res.pp,"residual")
sv <- variogram(residual~1, data=res.pp,cutoff=4000,width=480) # 480 za As
sv
plot(sv)

sv.fit =  fit.variogram(sv, vgm(1200, "Sph", 2000, 300)) #As

plot(sv, sv.fit)

##
coordinates(res.data) <- ~ longitude + latitude + altitude

proj4string(res.data) <- CRS(gk_7)

vr.def <- variogram(residual~1,res.data, cutoff=4000,width=480)  # 420 za As
plot(vr.def)

vgr.def <-fit.variogram(vr.def, vgm(1200, "Sph", 2000, 100,anis=c(0,0,0,1,(vc.fit$range/100)[2]/sv.fit$range[2])))  #za As
As.vgr <- vgr.def
vr.line <- (variogramLine(vgr.def, maxdist = max(vr.def$dist)))

plot(vr.def, vgr.def)

test.data <- As.IntP.cv$test.results[[i]]
test.data$residual <- test.data$observed - test.data$predicted
coordinates(test.data) <- ~ longitude + latitude + altitude
proj4string(test.data) <- CRS(gk_7)

test.data@data$krige.pred <- krige.cv(residual ~ 1, test.data, model = vgr.def)@data$var1.pred #,nfold=foldid

test.data$final.predicted <- test.data$predicted+test.data$krige.pred

As.IntP.cv$test.results[[i]] <- test.data@data

}

final.results <- do.call(rbind,As.IntP.cv$test.results)


caret::R2(pred=final.results$predicted,obs=final.results$observed, formula="traditional") # Rsquared for RK
RMSE(pred=final.results$predicted,obs=final.results$observed)


caret::R2(pred=final.results$final.predicted,obs=final.results$observed, formula="traditional") # Rsquared for RK
RMSE(pred=final.results$final.predicted,obs=final.results$observed)


#===================================================================================================================
#======================= SOM nested kros validacija ================================================================
#===================================================================================================================

SOM.IntP.cv <- penint3DP(fun=SOM.fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
SOM.IntP.cv$measure
summary(SOM.IntP.cv$measure[1:5,])


As.IntPH.cv <- penint3DP(fun=As.fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="poly")
As.IntPH.cv$measure
summary(As.IntPH.cv$measure[1:5,])

for(i in 1:5){
  
  TTP.p <- SOM.IntP.cv$train.results[[i]]
  TTP.p$residual <- TTP.p$observed-TTP.p$predicted
  
  res.profs<-data.frame(ID=TTP.p$ID,IntHP=TTP.p$residual)
  
  borind<-ddply(bor,.(ID))
  
  borind<-borind[complete.cases(borind$SOM),] # Ovde treba zameniti...
  
  "%ni%" <- Negate("%in%")
  
  ind <- which(borind[,"ID"] %in% unique(res.profs$ID))
  
  res.profs <- cbind(res.profs, borind[ind,c("Top","Bottom","Soil.Type","x","y","altitude")])
  
  res.list<-list(sites=ddply(res.profs,.(ID),summarize,y=y[1],x=x[1],Soil.Type=Soil.Type[1]),horizons=res.profs[,c("ID","Top","Bottom","IntHP")])
  res.profs<-join(res.list$sites,res.list$horizons, type="inner")
  
  res.profs.data<-res.profs
  
  depths(res.profs) <- ID ~ Top + Bottom
  site(res.profs) <- ~x+y
  coordinates(res.profs)<-~x+y
  proj4string(res.profs)<- CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m")
  
  
  res.data <- TTP.p
  
  spline.profs <- mpspline(res.profs, "IntHP",d = t(c(0,5,15,30,60,80)))
  spline.res <- spline.profs$var.1cm
  
  cm <- rep(1:80,dim(spline.res)[2])
  
  a <-ddply(res.data,.(ID),summarize ,longitude=longitude[1],latitude=latitude[1])
  
  b <- a[rep(seq_len(nrow(a)), each=80),]
  
  bm <- cbind(b,residual=as.numeric(spline.res),altitude=cm)
  
  bm <- na.omit(bm)
  
  vc <- variogram(residual~1, data=bm, locations=~longitude+latitude+altitude,cutoff=50,width=3)
  vc
  plot(vc)
  
  vc.fit = fit.variogram(vc, vgm(1.5, "Gau", 10, 0)) #SOM
  plot(vc, vc.fit)
  
  res.pp <- ddply(res.data, .(ID),summarize,x=longitude[1],y=latitude[1],altitude=altitude[1],residual=residual[1],observed=observed[1])
  coordinates(res.pp) <- ~x+y+altitude
  proj4string(res.pp)<- CRS(gk_7)
  
  bubble(res.pp,"residual")
  sv <- variogram(residual~1, data=res.pp,cutoff=4000,width= 420) # 480 za As
  sv
  plot(sv)
  
  sv.fit =  fit.variogram(sv, vgm(15, "Sph", 2000, 5)) #SOM
  
  plot(sv, sv.fit)
  
  ##
  coordinates(res.data) <- ~ longitude + latitude + altitude
  
  proj4string(res.data) <- CRS(gk_7)
  
  vr.def <- variogram(residual~1,res.data, cutoff=4000,width=550)  # 420 za SOM
  plot(vr.def)
  
  vgr.def <-fit.variogram(vr.def, vgm(15, "Sph", 2000, 5,anis=c(0,0,0,1,(vc.fit$range/100)[2]/sv.fit$range[2])))  #za SOM
  As.vgr <- vgr.def
  vr.line <- (variogramLine(vgr.def, maxdist = max(vr.def$dist)))
  
  plot(vr.def, vgr.def)
  
  test.data <- SOM.IntP.cv$test.results[[i]]
  test.data$residual <- test.data$observed - test.data$predicted
  coordinates(test.data) <- ~ longitude + latitude + altitude
  proj4string(test.data) <- CRS(gk_7)
  
  test.data@data$krige.pred <- krige.cv(residual ~ 1, test.data, model = vgr.def)@data$var1.pred #,nfold=foldid
  
  test.data$final.predicted <- test.data$predicted+test.data$krige.pred
  
  SOM.IntP.cv$test.results[[i]] <- test.data@data
  
}

final.results <- do.call(rbind,SOM.IntP.cv$test.results)


caret::R2(pred=final.results$predicted,obs=final.results$observed, formula="traditional") # Rsquared for RK
RMSE(pred=final.results$predicted,obs=final.results$observed)


caret::R2(pred=final.results$final.predicted,obs=final.results$observed, formula="traditional") # Rsquared for RK
RMSE(pred=final.results$final.predicted,obs=final.results$observed)



#===================================================================================================================
#======================= pH nested kros validacija ================================================================
#===================================================================================================================

pH.IntP.cv <- penint3DP(fun=pH.fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
pH.IntP.cv$measure
summary(pH.IntP.cv$measure[1:5,])


As.IntPH.cv <- penint3DP(fun=As.fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="poly")
As.IntPH.cv$measure
summary(As.IntPH.cv$measure[1:5,])

for(i in 1:5){
  
  TTP.p <- pH.IntP.cv$train.results[[i]]
  TTP.p$residual <- TTP.p$observed-TTP.p$predicted
  
  res.profs<-data.frame(ID=TTP.p$ID,IntHP=TTP.p$residual)
  
  borind<-ddply(bor,.(ID))
  
  borind<-borind[complete.cases(borind$pH),] # Ovde treba zameniti...
  
  "%ni%" <- Negate("%in%")
  
  ind <- which(borind[,"ID"] %in% unique(res.profs$ID))
  
  res.profs <- cbind(res.profs, borind[ind,c("Top","Bottom","Soil.Type","x","y","altitude")])
  
  res.list<-list(sites=ddply(res.profs,.(ID),summarize,y=y[1],x=x[1],Soil.Type=Soil.Type[1]),horizons=res.profs[,c("ID","Top","Bottom","IntHP")])
  res.profs<-join(res.list$sites,res.list$horizons, type="inner")
  
  res.profs.data<-res.profs
  
  depths(res.profs) <- ID ~ Top + Bottom
  site(res.profs) <- ~x+y
  coordinates(res.profs)<-~x+y
  proj4string(res.profs)<- CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m")
  
  
  res.data <- TTP.p
  
  spline.profs <- mpspline(res.profs, "IntHP",d = t(c(0,5,15,30,60,80)))
  spline.res <- spline.profs$var.1cm
  
  cm <- rep(1:80,dim(spline.res)[2])
  
  a <-ddply(res.data,.(ID),summarize ,longitude=longitude[1],latitude=latitude[1])
  
  b <- a[rep(seq_len(nrow(a)), each=80),]
  
  bm <- cbind(b,residual=as.numeric(spline.res),altitude=cm)
  
  bm <- na.omit(bm)
  
  vc <- variogram(residual~1, data=bm, locations=~longitude+latitude+altitude,cutoff=50,width=3)
  vc
  plot(vc)
  
  vc.fit = fit.variogram(vc, vgm(0.05, "Gau", 20, 0)) #pH
  
  plot(vc, vc.fit)
  
  res.pp <- ddply(res.data, .(ID),summarize,x=longitude[1],y=latitude[1],altitude=altitude[1],residual=residual[1],observed=observed[1])
  coordinates(res.pp) <- ~x+y
  proj4string(res.pp)<- CRS(gk_7)
  
  w <- c(650)
  
  #bubble(res.pp,"residual")
  sv <- variogram(residual~1, data=res.pp)#,cutoff=4000,width=650)
  sv
  plot(sv)
  
  sv.fit =  fit.variogram(sv, vgm(0.3, "Sph", 2000, 0.1)) #pH
  
  #sv.fit= autofitVariogram(residual~1, res.pp, miscFitOptions = list(merge.small.bins = TRUE))
  
  plot(sv, sv.fit)
  
  ##
  coordinates(res.data) <- ~ longitude + latitude + altitude
  
  proj4string(res.data) <- CRS(gk_7)
  
  w.def <- c(550)
  
  vr.def <- variogram(residual~1,res.data)#, cutoff=4000,width=550)
  plot(vr.def)
  
  vgr.def <-fit.variogram(vr.def, vgm(0.3, "Sph", 2000, 0.1,anis=c(0,0,0,1,(vc.fit$range/100)[2]/sv.fit$range[2])))  #za pH
  As.vgr <- vgr.def
  vr.line <- (variogramLine(vgr.def, maxdist = max(vr.def$dist)))
  
  plot(vr.def, vgr.def)
  
  test.data <- pH.IntP.cv$test.results[[i]]
  test.data$residual <- test.data$observed - test.data$predicted
  coordinates(test.data) <- ~ longitude + latitude + altitude
  proj4string(test.data) <- CRS(gk_7)
  
  test.data@data$krige.pred <- krige.cv(residual ~ 1, test.data, model = vgr.def)@data$var1.pred #,nfold=foldid
  
  test.data$final.predicted <- test.data$predicted+test.data$krige.pred
  
  pH.IntP.cv$test.results[[i]] <- test.data@data
  
}

final.results <- do.call(rbind,pH.IntP.cv$test.results)


caret::R2(pred=final.results$predicted,obs=final.results$observed, formula="traditional") # Rsquared for RK
RMSE(pred=final.results$predicted,obs=final.results$observed)


caret::R2(pred=final.results$final.predicted,obs=final.results$observed, formula="traditional") # Rsquared for RK
RMSE(pred=final.results$final.predicted,obs=final.results$observed)

############################################## Variogrami za krajnju predikciju #############################################################################


coordnames(gridmaps.sm2D) <- c("longitude","latitude")

pred <- As.pred

FFL.p <- pred$BaseL
FFP.p <- pred$BaseP
TFL.p <- pred$IntL
TFP.p <- pred$IntP
TTL.p <- pred$IntHP
TTP.p <- pred$IntHP

res.profs<-data.frame(ID=TFL.p$ID,BaseL=FFL.p$residual,BaseP=FFP.p$residual,IntL=TFL.p$residual,IntP=TFP.p$residual,IntHL=TTL.p$residual,IntHP=TTP.p$residual)

borind<-ddply(bor,.(ID))

borind<-borind[complete.cases(borind$As),] # Ovde treba zameniti...

"%ni%" <- Negate("%in%")

ind <- which(borind[,"ID"] %in% unique(res.profs$ID))

res.profs <- cbind(res.profs, borind[ind,c("Top","Bottom","Soil.Type","x","y","altitude")])

res.list<-list(sites=ddply(res.profs,.(ID),summarize,y=y[1],x=x[1],Soil.Type=Soil.Type[1]),horizons=res.profs[,c("ID","Top","Bottom","BaseL","BaseP","IntL","IntP","IntHL","IntHP")])
res.profs<-join(res.list$sites,res.list$horizons, type="inner")

res.profs.data<-res.profs

depths(res.profs) <- ID ~ Top + Bottom
site(res.profs) <- ~x+y
coordinates(res.profs)<-~x+y
proj4string(res.profs)<- CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m")


res.data <- pred$IntHP

spline.profs <- mpspline(res.profs, "IntHP",d = t(c(0,5,15,30,60,80)))
spline.res <- spline.profs$var.1cm

cm <- rep(1:80,dim(spline.res)[2])

a <-ddply(res.data,.(ID),summarize ,longitude=longitude[1],latitude=latitude[1])

b <- a[rep(seq_len(nrow(a)), each=80),]

bm <- cbind(b,residual=as.numeric(spline.res),altitude=cm)

bm <- na.omit(bm)

vc <- variogram(residual~1, data=bm, locations=~longitude+latitude+altitude,cutoff=50,width=3)
vc
plot(vc)

vc.fit = fit.variogram(vc, vgm(250, "Gau", 30, 5)) #As
plot(vc, vc.fit)

res.pp <- ddply(res.data, .(ID),summarize,x=longitude[1],y=latitude[1],altitude=altitude[1],residual=residual[1],observed=observed[1])
coordinates(res.pp) <- ~x+y+altitude
proj4string(res.pp)<- CRS(gk_7)

bubble(res.pp,"residual")
sv <- variogram(residual~1, data=res.pp,cutoff=4000,width=480) # 480 za As
sv
plot(sv)

sv.fit =  fit.variogram(sv, vgm(1200, "Sph", 2000, 300)) #As

plot(sv, sv.fit)

#=================================== plot stratified fold =======================================================
rdat <- bor
rdat <- plyr::rename(rdat, replace=c("x" = "longitude", "y" = "latitude"))
rdat <- rdat[complete.cases(rdat[,c("ID","longitude","latitude","altitude","As")]),c("ID","longitude","latitude","hdepth","altitude","As")] 

coordinates(rdat)<-~longitude+latitude
proj4string(rdat) <- CRS(utm)
rdat <- spTransform(rdat, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
rdat <- as.data.frame(rdat)
head(rdat)

#rdat <- rdat[complete.cases(rdat$As),]


rdat.folds <- stratfold3d(targetVar="As",regdat=rdat,folds=5,cent=3,seed=321,dimensions="3D",IDs=TRUE,sum=TRUE)
plotfolds(rdat.folds,targetvar="As")

rdat.folds$folds
flist<-rdat.folds$folds

folds.in.list <- as.list(rep(NA,length(flist)))
names(folds.in.list) <- paste("fold",c(1:length(flist)),sep = "")
foldid <- rep(NA,dim(res.data)[1])

for(j in 1:length(flist)){
  folds.in.list[[j]]<-which(res.data$ID %in% flist[[j]])
  foldid[folds.in.list[[j]]]<-j
}

foldid
#===============================================================================================================


#================================= 3D Kriging ==================================================================
coordinates(res.data) <- ~ longitude + latitude + altitude

proj4string(res.data) <- CRS(gk_7)

vr.def <- variogram(residual~1,res.data, cutoff=4000,width=480)  # 420 za As
plot(vr.def)

vgr.def <-fit.variogram(vr.def, vgm(1200, "Sph", 2000, 100,anis=c(0,0,0,1,(vc.fit$range/100)[2]/sv.fit$range[2])))  #za As
As.vgr <- vgr.def
vr.line <- (variogramLine(vgr.def, maxdist = max(vr.def$dist)))

plot(vr.def, vgr.def)
res.pred <- krige.cv(residual ~ 1, res.data, model = vgr.def) #,nfold=foldid


res.data@data$krige.res <- res.pred$var1.pred
res.data@data$krige.pred <- res.data@data$predicted + res.data@data$krige.res

#res.data <- res.data[complete.cases(res.data$krige.pred),]

caret::R2(pred=res.data@data$krige.pred,obs=res.data@data$observed, formula="traditional") # Rsquared for RK
RMSE(pred=res.data@data$krige.pred,obs=res.data@data$observed)

#====================================================================================================================
#================================== Raw 3D kriging ==================================================================

vr.obs <- variogram(observed~1,res.data, cutoff=4000,width=480)


vgr.obs <-fit.variogram(vr.obs, vgm(1200, "Exp", 2000, 100,anis=c(0,0,0,1,(vc.fit$range/100)[2]/sv.fit$range[2])))  #za As

vr.line.obs <- (variogramLine(vgr.obs, maxdist = max(vr.obs$dist)))
plot(vr.obs, vgr.obs)

obs.pred <- krige.cv(observed ~ 1, res.data, model = vgr.obs) #,nfold=foldid


caret::R2(pred=obs.pred$var1.pred,obs=res.data@data$observed, formula="traditional") # Rsquared for RK
RMSE(pred=obs.pred$var1.pred,obs=res.data@data$observed)

#================================= As results ======================================================================

As.reg.results <- defaultSummary(data.frame(pred=res.data@data$predicted,obs=res.data@data$observed))
As.raw.results <- defaultSummary(data.frame(pred=obs.pred$var1.pred,obs=res.data@data$observed))
As.rk.results <- defaultSummary(data.frame(pred=res.data@data$krige.pred,obs=res.data@data$observed))


#================================== Variogram plot ==================================================================

addlinetoplot <- function(dataset, varx, vary,color="#000000",sajz=0.7) { 
  list(
    geom_line(data=dataset, aes_string(x=varx, y=vary), colour=color,size=sajz)
  )
}

addpointtoplot <- function(dataset, varx, vary, size=3, color="#CC0000") { 
  list(
    geom_point(data=dataset, aes_string(x=varx, y=vary))
  )
}

vr <- rbind(vr.obs,vr.def)
vr$id <- c(rep("var1",dim(vr.obs)[1]),rep("var2",dim(vr.def)[1]))


var.plot <- qplot(dist, gamma, data = vr, geom = c("point"),color=id)+theme_bw() + 
  theme(legend.position="none")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  labs(x="Distance",y="Semivariance") + addlinetoplot(vr.line, varx = "dist", vary = "gamma",color="#56B4E9",sajz=0.65) + 
  addlinetoplot(vr.line.obs, varx = "dist", vary = "gamma",color="#CC0000",sajz=0.65) + addpointtoplot(vr.obs, varx = "dist", vary = "gamma")

var.plot


pdf("AsVariograms.pdf",width=8,height=8)
var.plot # Make plot
dev.off()

#====================================================================================================================
#====================================================================================================================
#=============================== 3D KRIGING =========================================================================
#====================================================================================================================
#================================== SOM =============================================================================

coordnames(gridmaps.sm2D) <- c("longitude","latitude")

pred <- SOM.pred

FFL.p <- pred$BaseL
FFP.p <- pred$BaseP
TFL.p <- pred$IntL
TFP.p <- pred$IntP
TTL.p <- pred$IntHP
TTP.p <- pred$IntHP

res.profs<-data.frame(ID=TFL.p$ID,BaseL=FFL.p$residual,BaseP=FFP.p$residual,IntL=TFL.p$residual,IntP=TFP.p$residual,IntHL=TTL.p$residual,IntHP=TTP.p$residual)

borind<-ddply(bor,.(ID))

borind<-borind[complete.cases(borind$SOM),] # Ovde treba zameniti...

"%ni%" <- Negate("%in%")

ind <- which(borind[,"ID"] %in% unique(res.profs$ID))

res.profs <- cbind(res.profs, borind[ind,c("Top","Bottom","Soil.Type","x","y","altitude")])

res.list<-list(sites=ddply(res.profs,.(ID),summarize,y=y[1],x=x[1],Soil.Type=Soil.Type[1]),horizons=res.profs[,c("ID","Top","Bottom","BaseL","BaseP","IntL","IntP","IntHL","IntHP")])
res.profs<-join(res.list$sites,res.list$horizons, type="inner")

res.profs.data<-res.profs

depths(res.profs) <- ID ~ Top + Bottom
site(res.profs) <- ~x+y
coordinates(res.profs)<-~x+y
proj4string(res.profs)<- CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m")


res.data <- pred$IntHP

spline.profs <- mpspline(res.profs, "IntHP",d = t(c(0,5,15,30,60,80)))
spline.res <- spline.profs$var.1cm

cm <- rep(1:80,dim(spline.res)[2])

a <-ddply(res.data,.(ID),summarize ,longitude=longitude[1],latitude=latitude[1])

b <- a[rep(seq_len(nrow(a)), each=80),]

bm <- cbind(b,residual=as.numeric(spline.res),altitude=cm)

bm <- na.omit(bm)

vc <- variogram(residual~1, data=bm, locations=~longitude+latitude+altitude,cutoff=50,width=3)
vc
plot(vc)


vc.fit = fit.variogram(vc, vgm(1.5, "Gau", 10, 0)) #SOM

plot(vc, vc.fit)

res.pp <- ddply(res.data, .(ID),summarize,x=longitude[1],y=latitude[1],altitude=altitude[1],residual=residual[1],observed=observed[1])
coordinates(res.pp) <- ~x+y+altitude
proj4string(res.pp)<- CRS(gk_7)

bubble(res.pp,"residual")
sv <- variogram(residual~1, data=res.pp,cutoff=4000,width=420) 
sv
plot(sv)

sv.fit =  fit.variogram(sv, vgm(15, "Sph", 2000, 5)) #SOM

plot(sv, sv.fit)

#=================================== plot stratified fold =======================================================
rdat <- bor
rdat <- plyr::rename(rdat, replace=c("x" = "longitude", "y" = "latitude"))
rdat <- rdat[complete.cases(rdat[,c("ID","longitude","latitude","altitude","SOM")]),c("ID","longitude","latitude","hdepth","altitude","SOM")] 

coordinates(rdat)<-~longitude+latitude
proj4string(rdat) <- CRS(utm)
rdat <- spTransform(rdat, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
rdat <- as.data.frame(rdat)
head(rdat)

#rdat <- rdat[complete.cases(rdat$SOM),]


rdat.folds <- stratfold3d(targetVar="SOM",regdat=rdat,folds=5,cent=3,seed=321,dimensions="3D",IDs=TRUE,sum=TRUE)
plotfolds(rdat.folds,targetvar="SOM")

rdat.folds$folds
flist<-rdat.folds$folds

folds.in.list <- as.list(rep(NA,length(flist)))
names(folds.in.list) <- paste("fold",c(1:length(flist)),sep = "")
foldid <- rep(NA,dim(res.data)[1])

for(j in 1:length(flist)){
  folds.in.list[[j]]<-which(res.data$ID %in% flist[[j]])
  foldid[folds.in.list[[j]]]<-j
}

foldid
#===============================================================================================================


#================================= 3D Kriging ==================================================================
coordinates(res.data) <- ~ longitude + latitude + altitude

proj4string(res.data) <- CRS(gk_7)

vr.def <- variogram(residual~1,res.data, cutoff=4000,width=420)
plot(vr.def)
vgr.def <-fit.variogram(vr.def, vgm(15, "Sph", 2000, 5,anis=c(0,0,0,1,(vc.fit$range/100)[2]/sv.fit$range[2])))
SOM.vgr <- vgr.def
vr.line <- (variogramLine(vgr.def, maxdist = max(vr.def$dist)))

plot(vr.def, vgr.def)

res.pred <- krige.cv(residual ~ 1, res.data, model = vgr.def) #,nfold=foldid


res.data@data$krige.res <- res.pred$var1.pred
res.data@data$krige.pred <- res.data@data$predicted + res.data@data$krige.res

#res.data <- res.data[complete.cases(res.data$krige.pred),]

caret::R2(pred=res.data@data$krige.pred,obs=res.data@data$observed, formula="traditional") # Rsquared for RK
RMSE(pred=res.data@data$krige.pred,obs=res.data@data$observed)

#====================================================================================================================
#================================== Raw 3D kriging ==================================================================

vr.obs <- variogram(observed~1,res.data, cutoff=4000,width=420)

vgr.obs <-fit.variogram(vr.obs, vgm(15, "Gau", 2000, 5,anis=c(0,0,0,1,(vc.fit$range/100)[2]/sv.fit$range[2])))

vr.line.obs <- (variogramLine(vgr.obs, maxdist = max(vr.obs$dist)))
plot(vr.obs, vgr.obs)

obs.pred <- krige.cv(observed ~ 1, res.data, model = vgr.obs) #,nfold=foldid


caret::R2(pred=obs.pred$var1.pred,obs=res.data@data$observed, formula="traditional") # Rsquared for RK
RMSE(pred=obs.pred$var1.pred,obs=res.data@data$observed)

#================================= pH results ======================================================================

SOM.reg.results <- defaultSummary(data.frame(pred=res.data@data$predicted,obs=res.data@data$observed))
SOM.raw.results <- defaultSummary(data.frame(pred=obs.pred$var1.pred,obs=res.data@data$observed))
SOM.rk.results <- defaultSummary(data.frame(pred=res.data@data$krige.pred,obs=res.data@data$observed))

#================================== Variogram plot ==================================================================

addlinetoplot <- function(dataset, varx, vary,color="#000000",sajz=0.7) { 
  list(
    geom_line(data=dataset, aes_string(x=varx, y=vary), colour=color,size=sajz)
  )
}

addpointtoplot <- function(dataset, varx, vary, size=3, color="#CC0000") { 
  list(
    geom_point(data=dataset, aes_string(x=varx, y=vary))
  )
}

vr <- rbind(vr.obs,vr.def)
vr$id <- c(rep("var1",dim(vr.obs)[1]),rep("var2",dim(vr.def)[1]))


var.plot <- qplot(dist, gamma, data = vr, geom = c("point"),color=id)+theme_bw() + 
  theme(legend.position="none")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  labs(x="Distance",y="Semivariance") + addlinetoplot(vr.line, varx = "dist", vary = "gamma",color="#56B4E9",sajz=0.65) + 
  addlinetoplot(vr.line.obs, varx = "dist", vary = "gamma",color="#CC0000",sajz=0.65) + addpointtoplot(vr.obs, varx = "dist", vary = "gamma")

var.plot


pdf("SOMVariograms.pdf",width=8,height=8)
var.plot # Make plot
dev.off()

#====================================================================================================================



#====================================================================================================================
#====================================================================================================================
#=============================== 3D KRIGING =========================================================================
#====================================================================================================================
#================================== pH =============================================================================

coordnames(gridmaps.sm2D) <- c("longitude","latitude")

pred <- pH.pred

FFL.p <- pred$BaseL
FFP.p <- pred$BaseP
TFL.p <- pred$IntL
TFP.p <- pred$IntP
TTL.p <- pred$IntHP
TTP.p <- pred$IntHP

res.profs<-data.frame(ID=TFL.p$ID,BaseL=FFL.p$residual,BaseP=FFP.p$residual,IntL=TFL.p$residual,IntP=TFP.p$residual,IntHL=TTL.p$residual,IntHP=TTP.p$residual)

borind<-ddply(bor,.(ID))

borind<-borind[complete.cases(borind$pH),] # Ovde treba zameniti...

"%ni%" <- Negate("%in%")

ind <- which(borind[,"ID"] %in% unique(res.profs$ID))

res.profs <- cbind(res.profs, borind[ind,c("Top","Bottom","Soil.Type","x","y","altitude")])

res.list<-list(sites=ddply(res.profs,.(ID),summarize,y=y[1],x=x[1],Soil.Type=Soil.Type[1]),horizons=res.profs[,c("ID","Top","Bottom","BaseL","BaseP","IntL","IntP","IntHL","IntHP")])
res.profs<-join(res.list$sites,res.list$horizons, type="inner")

res.profs.data<-res.profs

depths(res.profs) <- ID ~ Top + Bottom
site(res.profs) <- ~x+y
coordinates(res.profs)<-~x+y
proj4string(res.profs)<- CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m")


res.data <- pred$IntHP

spline.profs <- mpspline(res.profs, "IntHP",d = t(c(0,5,15,30,60,80)))
spline.res <- spline.profs$var.1cm

cm <- rep(1:80,dim(spline.res)[2])

a <-ddply(res.data,.(ID),summarize ,longitude=longitude[1],latitude=latitude[1])

b <- a[rep(seq_len(nrow(a)), each=80),]

bm <- cbind(b,residual=as.numeric(spline.res),altitude=cm)

bm <- na.omit(bm)

vc <- variogram(residual~1, data=bm, locations=~longitude+latitude+altitude,cutoff=50,width=3)
vc
plot(vc)

vc.fit = fit.variogram(vc, vgm(0.05, "Gau", 20, 0)) #pH

plot(vc, vc.fit)

res.pp <- ddply(res.data, .(ID),summarize,x=longitude[1],y=latitude[1],altitude=altitude[1],residual=residual[1],observed=observed[1])
coordinates(res.pp) <- ~x+y
proj4string(res.pp)<- CRS(gk_7)

#bubble(res.pp,"residual")
sv <- variogram(residual~1, data=res.pp,cutoff=4000,width=620)
sv
plot(sv)

sv.fit =  fit.variogram(sv, vgm(0.3, "Sph", 1000, 0.1)) #pH

plot(sv, sv.fit)


#=================================== plot stratified fold =======================================================
rdat <- bor
rdat <- plyr::rename(rdat, replace=c("x" = "longitude", "y" = "latitude"))
rdat <- rdat[complete.cases(rdat[,c("ID","longitude","latitude","altitude","pH")]),c("ID","longitude","latitude","hdepth","altitude","pH")] 

coordinates(rdat)<-~longitude+latitude
proj4string(rdat) <- CRS(utm)
rdat <- spTransform(rdat, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
rdat <- as.data.frame(rdat)
head(rdat)

#rdat <- rdat[complete.cases(rdat$pH),]


rdat.folds <- stratfold3d(targetVar="pH",regdat=rdat,folds=5,cent=3,seed=321,dimensions="3D",IDs=TRUE,sum=TRUE)
plotfolds(rdat.folds,targetvar="pH")

rdat.folds$folds
flist<-rdat.folds$folds

folds.in.list <- as.list(rep(NA,length(flist)))
names(folds.in.list) <- paste("fold",c(1:length(flist)),sep = "")
foldid <- rep(NA,dim(res.data)[1])

for(j in 1:length(flist)){
  folds.in.list[[j]]<-which(res.data$ID %in% flist[[j]])
  foldid[folds.in.list[[j]]]<-j
}

foldid
#===============================================================================================================


#================================= 3D Kriging ==================================================================
coordinates(res.data) <- ~ longitude + latitude + altitude

proj4string(res.data) <- CRS(gk_7)

vr.def <- variogram(residual~1,res.data, cutoff=4000,width=650)
plot(vr.def)

vgr.def <-fit.variogram(vr.def, vgm(0.3, "Sph", 2000, 0.05,anis=c(0,0,0,1,(vc.fit$range/100)[2]/sv.fit$range[2])))
pH.vgr <- vgr.def
vr.line <- (variogramLine(vgr.def, maxdist = max(vr.def$dist)))

plot(vr.def, vgr.def)
res.pred <- krige.cv(residual ~ 1, res.data, model = vgr.def) #,nfold=foldid


res.data@data$krige.res <- res.pred$var1.pred
res.data@data$krige.pred <- res.data@data$predicted + res.data@data$krige.res

#res.data <- res.data[complete.cases(res.data$krige.pred),]

caret::R2(pred=res.data@data$krige.pred,obs=res.data@data$observed, formula="traditional") # Rsquared for RK
RMSE(pred=res.data@data$krige.pred,obs=res.data@data$observed)

#====================================================================================================================
#================================== Raw 3D kriging ==================================================================

vr.obs <- variogram(observed~1,res.data, cutoff=4000,width=650)

vgr.obs <-fit.variogram(vr.obs, vgm(0.3, "Sph", 2000, 0.05,anis=c(0,0,0,1,(vc.fit$range/100)[2]/sv.fit$range[2])))

vr.line.obs <- (variogramLine(vgr.obs, maxdist = max(vr.obs$dist)))
plot(vr.obs, vgr.obs)

obs.pred <- krige.cv(observed ~ 1, res.data, model = vgr.obs) #,nfold=foldid

caret::R2(pred=obs.pred$var1.pred,obs=res.data@data$observed, formula="traditional") # Rsquared for RK
RMSE(pred=obs.pred$var1.pred,obs=res.data@data$observed)

#================================= pH results ======================================================================

pH.reg.results <- defaultSummary(data.frame(pred=res.data@data$predicted,obs=res.data@data$observed))
pH.raw.results <- defaultSummary(data.frame(pred=obs.pred$var1.pred,obs=res.data@data$observed))
pH.rk.results <- defaultSummary(data.frame(pred=res.data@data$krige.pred,obs=res.data@data$observed))


#================================== Variogram plot ==================================================================

addlinetoplot <- function(dataset, varx, vary,color="#000000",sajz=0.7) { 
  list(
    geom_line(data=dataset, aes_string(x=varx, y=vary), colour=color,size=sajz)
  )
}

addpointtoplot <- function(dataset, varx, vary, size=3, color="#CC0000") { 
  list(
    geom_point(data=dataset, aes_string(x=varx, y=vary))
  )
}

vr <- rbind(vr.obs,vr.def)
vr$id <- c(rep("var1",dim(vr.obs)[1]),rep("var2",dim(vr.def)[1]))


var.plot <- qplot(dist, gamma, data = vr, geom = c("point"),color=id)+theme_bw() + 
  theme(legend.position="none")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
  labs(x="Distance",y="Semivariance") + addlinetoplot(vr.line, varx = "dist", vary = "gamma",color="#56B4E9",sajz=0.65) + 
  addlinetoplot(vr.line.obs, varx = "dist", vary = "gamma",color="#CC0000",sajz=0.65) + addpointtoplot(vr.obs, varx = "dist", vary = "gamma")

var.plot


pdf("pHVariograms.pdf",width=8,height=8)
var.plot # Make plot
dev.off()

#====================================================================================================================

As.results <- rbind(As.reg.results,As.raw.results,As.rk.results)
SOM.results <- rbind(SOM.reg.results,SOM.raw.results,SOM.rk.results)
pH.results <- rbind(pH.reg.results,pH.raw.results,pH.rk.results)

results <- cbind(As.results,SOM.results,pH.results)
rownames(results) <- c("lasso","3D OK", "3D RK")

stargazer(results,summary=FALSE,digits=2,type="latex")

#========================= Variogram table =========================================================================

As.vgr
SOM.vgr
pH.vgr

As.vgr.data <- c(As.vgr$psill[1], As.vgr$psill[1]+As.vgr$psill[2], As.vgr$range[2], 0.05/As.vgr$anis2[2], As.vgr$psill[1]/(As.vgr$psill[1]+As.vgr$psill[2]))
SOM.vgr.data <- c(SOM.vgr$psill[1], SOM.vgr$psill[1]+SOM.vgr$psill[2], SOM.vgr$range[2], 0.05/SOM.vgr$anis2[2], SOM.vgr$psill[1]/(SOM.vgr$psill[1]+SOM.vgr$psill[2]))
pH.vgr.data <- c(pH.vgr$psill[1], pH.vgr$psill[1]+pH.vgr$psill[2], pH.vgr$range[2], 0.05/pH.vgr$anis2[2], pH.vgr$psill[1]/(pH.vgr$psill[1]+pH.vgr$psill[2]))

vgr.data <-rbind(As.vgr.data,SOM.vgr.data,pH.vgr.data)
rownames(vgr.data) <- c("As","SOM","pH")
colnames(vgr.data) <- c("Nugget","Sill","Range","Anisotropy (5cm depth=)","Range/Sill")

stargazer(vgr.data,summary=FALSE,digits=2,type="latex")

#====================================================================================================================
As.gls.res.sk3@data

vgr.gsif <- fit.gstatModel(hm.geo,residual~1,gridmaps.sm2D, model=vgr$vgm, vgmFun = "Exp")



res.geo<-as.geosamples(res.profs)



vgr.gsif <- fit.gstatModel(bor.geo,residual~1,gridmaps.sm2D, model=vgr$vgm, vgmFun = "Exp")
plot(vgr.gsif)

plot(variogram(residual~1, rmatrix=TFL), vgr$vgm)




#======================================================================================================


#s <- slice(res.profs, seq(10,30,10) ~ .,just.the.data=TRUE)
s <- slice(res.profs, 15 ~ .)

#s$s <- rep(1:3,dim(s)[1]/3)
s<-s[-which(is.na(s$.pctMissing)),]
sv <- variogram(IntHP~1, data=s,cutoff=4000,width=500)

plot(sv)#,type="o",add=TRUE)
sv.fit = fit.variogram(sv, vgm(1500, "Sph", 2000, 50)) #,fit.range=FALSE
plot(sv, sv.fit)

sv10<-sv
sv10$id <- rep("10 cm",dim(sv10)[1])
sv20<-sv
sv20$id <- rep("20 cm",dim(sv20)[1])
sv30<-sv
sv30$id <- rep("30 cm",dim(sv30)[1])


sv.sve<-rbind(sv10[,c(1:3,6)],sv20[,c(1:3,6)],sv30[,c(1:3,6)])
sv.sve$depth <- factor(sv.sve$id)

As.sve <- sv.sve
SOM.sve <- sv.sve
pH.sve <- sv.sve

As.vars <- qplot(dist, gamma, data = As.sve, geom = c("point", "line"), color = depth)+theme_bw()+theme(legend.text=element_text(size=12))+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+labs(x="Distance",y="Semivariance")
SOM.vars <- qplot(dist, gamma, data = SOM.sve, geom = c("point", "line"), color = depth)+theme_bw()+theme(legend.text=element_text(size=12))+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+labs(x="Distance",y="Semivariance")
pH.vars <- qplot(dist, gamma, data = pH.sve, geom = c("point", "line"), color = depth)+theme_bw()+theme(legend.text=element_text(size=12))+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+labs(x="Distance",y="Semivariance")


pdf("AsVariograms.pdf",width=5,height=5)
As.vars # Make plot
dev.off()

pdf("SOMVariograms.pdf",width=5,height=5)
SOM.vars # Make plot
dev.off()

pdf("pHVariograms.pdf",width=5,height=5)
pH.vars # Make plot
dev.off()





########################################################################################################
TFL.pp <- TFL.pp[complete.cases(TFL.pp@data[,c("residual")]),]
Bor.std<-Bor.std[,c("As1std","dmt","dwe","dwec","raspect","rsaspect","protection_index","we285R10000")]


vc<-variogram(residual~1,TFL.pp,cloud=TRUE)
#quantile(vc$gamma,0.7)
vc<-vc[vc$dist<300 & vc$gamma>quantile(vc$gamma,0.80),]

vci<-(sort(table(c(data.frame(vc)[,6],data.frame(vc)[,7])),decreasing=TRUE))
vcit<-as.numeric(dimnames(vci)[[1]])

#std1$gls.resc<-abs(std1$gls.res)+10

dgmplot<-function(spdf, ind, attrdis, maxrad=100){
  spd<-TFL.pp#spdf
  spd$sfactor<-as.factor(NA)
  spd$rn<-row(spd@data)[,1]
  spd@data<-spd@data[,c(attrdis,"rn","sfactor")]
  levels(spd$sfactor)<-c(1,2)
  spd@data[c(which(row(spd@data)[,1] %in% ind)),"sfactor"]<- levels(spd$sfactor)[1]
  spd@data[-c(which(row(spd@data)[,1] %in% ind)),"sfactor"]<- levels(spd$sfactor)[2]
  obs.bubble<-bubbleSP(spd, attrdis[1],max.radius=maxrad,do.sqrt=FALSE, scale_e=1)
  obs.bubble$sfactor<-as.factor(obs.bubble$sfactor)
  plotGoogleMaps(obs.bubble, zcol="sfactor")
}

dgmplot(TFL.pp,ind=vcit,attrdis="residual")

dgmplot(std3,vcit,"gls.res")

#############################################################################################################

v.fit = fit.variogram(v, vgm(0.05, "Sph", 0.15, 0.001))
plot(v, v.fit)

res.geo <- as.geosamples(res.profs)

vgr <- fit.vgmModel(residual~1, rmatrix=res.data, gridmaps.sm2D, anis=c(0,0,0,1,0.000112964),vgmFun = "Sph", dimensions="3D") #0.000112964

vgr.gsif <- fit.gstatModel(bor.geo,As.fun,gridmaps.sm2D, model=vgr$vgm, vgmFun = "Sph")
plot(vgr.gsif)

plot(variogram(residual~1, rmatrix=TFL), vgr$vgm)

plot(vgr$svgm,vgr$vgm)
