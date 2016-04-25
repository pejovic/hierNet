
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

#load("SOMPrediction2732016.RData")

gk_7 <- "+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m"
utm <- "+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m"

load(paste(getwd(),"inst","extdata","gridmaps.RDa",sep="/"))
gridmaps.sm2D$CD <- exp(-gridmaps.sm2D$DirDif)
gridmaps.sm2D$DD <- as.numeric(exp(-scale(gridmaps.sm2D$Dist,center=FALSE)))

sm2D.lst<-names(gridmaps.sm2D)
sm2D.lst <- sm2D.lst[ -which(sm2D.lst %in% c("DirDif","Dist","AnalyticalHills","LSFactor","RelSlopePosition","VelleyDepth","optional")) ] 

#gridmaps <- gridmaps.sm2D[,sm2D.lst]
#save(gridmaps,file="C:/Users/Milutin/Dropbox/ES3D/Data and Scripts/gridmaps.rda")


bor <- join(read.csv(paste(getwd(),"inst","extdata","Profili_sredjeno_csv.csv",sep="/")), read.csv(paste(getwd(),"inst","extdata","Koordinate_csv.csv",sep="/")), type="inner")
#save(bor,file="C:/Users/Milutin/Dropbox/ES3D/Data and Scripts/BorData.rda")

bor$hdepth<-bor$Bottom-bor$Top
bor$altitude <- - (bor$Top / 100 + (( bor$Bottom - bor$Top ) / 2) / 100)
bor <- bor[, - c(7:12,14,15,16,17 )]
names(bor) <- c("Soil.Type","ID","Horizont","Top" , "Bottom","pH","SOM","As","Co","x","y","hdepth","altitude")

bor <- bor[-which(bor$ID %in% c(66,129,14,51,69,130,164,165,166)),]

#========================= Creating Soil Profile Collections ====================================
bor.profs <- bor[,c("ID","x","y","Soil.Type","Top","Bottom","SOM","pH","Co","As")]
depths(bor.profs) <- ID ~ Top + Bottom
site(bor.profs) <- ~ Soil.Type + x + y
coordinates(bor.profs) <- ~ x+y
proj4string(bor.profs) <- CRS(utm)
bor.geo<-as.geosamples(bor.profs)


#====================== formulas ================================================================
As.fun <- as.formula(paste("As ~", paste(c(sm2D.lst,"altitude"), collapse="+")))
SOM.fun <- as.formula(paste("SOM ~", paste(c(sm2D.lst[-which(sm2D.lst %in% c("ES","CD","DD"))],"altitude"), collapse="+")))
pH.fun <- as.formula(paste("pH ~", paste(c(sm2D.lst[-which(sm2D.lst %in% c("ES","CD","DD"))],"altitude"), collapse="+")))

#================================================================================================

#================== test for stratfold3d and penint3D ============================================

source(paste(getwd(),"R","stratFold3D.R",sep="/"))
source(paste(getwd(),"R","penint3D_def.R",sep="/"))
source(paste(getwd(),"R","plotfolds.R",sep="/"))
source(paste(getwd(),"R","predint3D.R",sep="/"))
source(paste(getwd(),"R","penint3D_defP.R",sep="/"))
source(paste(getwd(),"R","predint3DP.R",sep="/"))

fun <- As.fun


#=================================== plot stratified fold ============================================
rdat <- bor
rdat <- plyr::rename(rdat, replace=c("x" = "longitude", "y" = "latitude"))
rdat <- rdat[complete.cases(rdat[,c("ID","longitude","latitude","altitude","SOM")]),c("ID","longitude","latitude","hdepth","altitude","SOM")] 

coordinates(rdat)<-~longitude+latitude
proj4string(rdat) <- CRS(utm)
rdat <- spTransform(rdat, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
rdat <- as.data.frame(rdat)
head(rdat)


rdat.folds <- stratfold3d(targetVar="SOM",regdat=rdat,folds=5,cent=3,seed=666,dimensions="3D",IDs=TRUE,sum=TRUE)
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

BaseLcv <- penint3D(fun=fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="linear")
BaseLcv$measure
summary(BaseLcv$measure[1:5,])


IntLcv <- penint3D(fun=fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="linear")
IntLcv$measure
summary(IntLcv$measure[1:5,])


IntHLcv <- penint3D(fun=fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = TRUE, int=FALSE, depth.fun="linear")
IntHLcv$measure
summary(IntHLcv$measure[1:5,])


#============= Poly ========================
BasePcv <- penint3D(fun=fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=FALSE, depth.fun="poly")
BasePcv$measure
summary(BasePcv$measure[1:5,])


IntPcv <- penint3D(fun=fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = FALSE, int=TRUE, depth.fun="poly")
IntPcv$measure
summary(IntPcv$measure[1:5,])


IntPHcv <- penint3D(fun=fun, profs = bor.profs, seed=443, cogrids = gridmaps.sm2D, hier = TRUE, int=TRUE, depth.fun="poly")
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

BaseL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
BaseP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=FALSE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = FALSE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

IntHL <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="linear",cores=8)
IntHP <- predint3DP(fun=fun, profs = bor.profs, cogrids = gridmaps.sm2D, pred=FALSE ,hier = TRUE, int=TRUE, depths=c(-0.1,-0.2,-0.3) ,depth.fun="poly",cores=8)

#============================= Coefficients matrix ==========================================
cmL <- data.frame(variable=IntHL$summary$coefficients[,1], BaseL.me=BaseL$summary$coefficients[2:27], IntL.me=IntL$summary$coefficients[2:27],IntL.ie=c(IntL$summary$coefficients[28:52],0),IntHL.me=IntHL$summary$coefficients[,2],IntHL.ie=IntHL$summary$coefficients[,3] )

l <- length(IntP$summary$coefficients)
i1 <- seq(1,l-29,3)
i2 <- seq(2,l-29,3)
i3 <- seq(3,l-29,3)

cmP <- data.frame(variable=IntHP$summary$coefficients[,1], BaseP.me=BaseP$summary$coefficients[2:29], IntP.me=IntP$summary$coefficients[2:29],IntP.ie1=c(IntP$summary$coefficients[30:l][i1],0,0,0),IntP.ie2=c(IntP$summary$coefficients[30:l][i2],0,0,0),IntP.ie3=c(IntP$summary$coefficients[30:l][i3],0,0,0),IntHP.me=IntHP$summary$coefficients[,2],IntHP.ie1=IntHP$summary$coefficients[,4],IntHP.ie2=IntHP$summary$coefficients[,5],IntHP.ie3=IntHP$summary$coefficients[,6] )
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

plot(subset(pred.stack, subset=sort(grep("TF", names(prediction), value=TRUE)) , col=color.pal, breaks=seq(minR,maxR,length.out=30)))


#=========================== Profile plots ===================================================

FFL.p <- BaseL$summary$pred
FFP.p <- BaseP$summary$pred
TFL.p <- IntL$summary$pred
TFP.p <- IntP$summary$pred
TTL.p <- IntHL$summary$pred
TTP.p <- IntHP$summary$pred

names(FFL.p)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")
names(FFP.p)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")
names(TFL.p)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")
names(TFP.p)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")
names(TTL.p)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")
names(TTP.p)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")

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
agg <- slab(bor.profs, fm= ~ As+SOM+pH, slab.structure=seq(0,70,5))

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




#======================= Computin variogram ================================================

v <- variogram(residual~1, data=TTL.p, locations=~x+y+altitude,cutoff=0.7,width=0.05)
v
plot(v)

v.fit = fit.variogram(v, vgm(0.05, "Sph", 0.15, 0.001))
plot(v, v.fit)

TFL.pp <- ddply(TFL.p, .(ID),summarize,x=x[1],y=y[1],altitude=altitude[1],residual=residual[1])
coordinates(TFL.pp) <- ~x+y
proj4string(TFL.pp)<- CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m")

bubble(TFL.pp,"residual")
sv <- variogram(residual~1, data=TFL.pp,cutoff=4000,width=470)
sv
plot(sv)
sv.fit = fit.variogram(sv, vgm(0.2, "Sph", 700, 0.1)) #,fit.range=FALSE
plot(sv, sv.fit)


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



v.fit = fit.variogram(v, vgm(0.05, "Sph", 0.15, 0.001))
plot(v, v.fit)



vgr <- fit.vgmModel(residual~1, rmatrix=TFL.p, gridmaps.sm2D, anis=c(0,0,0,1,0.0007635)  , dimensions="3D")

plot(variogram(residual~1, rmatrix=TFL), vgr$vgm)

plot(vgr$svgm,vgr$vgm)