
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

gk_7 <- "+proj=tmerc +lat_0=0 +lon_0=21 +k=0.9999 +x_0=7500000 +y_0=0 +ellps=bessel +towgs84=574.027,170.175,401.545,4.88786,-0.66524,-13.24673,0.99999311067 +units=m"
utm <- "+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m"

load(paste(getwd(),"inst","extdata","gridmaps.RDa",sep="/"))
bor <- join(read.csv(paste(getwd(),"inst","extdata","Profili_sredjeno_csv.csv",sep="/")), read.csv(paste(getwd(),"inst","extdata","Koordinate_csv.csv",sep="/")), type="inner")
bor$altitude <- - (bor$Top / 100 + (( bor$Bottom - bor$Top ) / 2) / 100)
bor <- bor[, - c( 7:12,14,15,16,17 )]
names(bor) <- c("Soil.Type","ID","Horizont","Top" , "Bottom","pH","Humus","As","Co","x","y","altitude")

bor.Humus.xy <- bor[complete.cases(bor[,c("ID","x","y","altitude","Humus")]),c("ID","x","y","altitude","Humus")] 
bor.Humus.xy <- plyr::rename(bor.Humus.xy, replace=c("x" = "longitude", "y" = "latitude"))
coordinates(bor.Humus.xy) <- ~ longitude + latitude
proj4string(bor.Humus.xy) <- CRS(utm)
bor.Humus.xy <- spTransform(bor.Humus.xy, CRS(gk_7))

sm2D.lst<-names(gridmaps.sm2D)

#=========== Overlay with covariates (gridmaps.sm2D) ================================
ov.H <- over(bor.Humus.xy, gridmaps.sm2D)
bor.Humus.xy@data[,c(sm2D.lst)] <- ov.H
str(bor.Humus.xy)

bor.Humus.df <- data.frame(bor.Humus.xy)
bor.Humus.df <- bor.Humus.df[complete.cases(bor.Humus.df[,names(bor.Humus.df)]),]
str(bor.Humus.df)

#============== Removing correlated covariates (previously determined) ==============
names(bor.Humus.df)
bor.Humus.df<-bor.Humus.df[,-which(names(bor.Humus.df) %in% c("AnalyticalHills","LSFactor","ChannelNetBaseLevel","RelSlopePosition","CatchArea","WindEffectNorthWest","VelleyDepth","optional"))]
str(bor.Humus.df)

#bor.Humus.df<-bor.Humus.df[-which(bor.Humus.df$ID %in% c(141)),] #153,137,141,145,150

#============== Sorting by altitude =================================================
H.df<-ddply(bor.Humus.df,.(altitude))
H.df <- H.df[order(-H.df$altitude),] 
head(H.df)
names(H.df)

sm2D.lst<-names(H.df)[c(4:16)]

#====================== formulas ===================================================
fm.H <- as.formula(paste("Humus ~", paste(c(sm2D.lst,"altitude"), collapse="+")))
fm.int.H<- as.formula(paste("Humus ~",paste(sm2D.lst,"altitude",sep="*", collapse="+"),sep=""))
fm.GSIF.H <- as.formula(paste("Humus ~", paste(c(sm2D.lst,"altitude","ns(altitude,df=4)"), collapse="+")))
fm.GSIF.int.H <- as.formula(paste("Humus ~",paste(paste(sm2D.lst,"altitude",sep="*" ,collapse="+"),"ns(altitude,df=4)",sep="+")))
fm.GSIF.int.H1 <- as.formula(paste("Humus ~",paste(paste(sm2D.lst,"poly(altitude,3)",sep="*",collapse="+"))))
#===============================================================================================

#================== test for stratfold3d and penint3D ============================================
source(paste(getwd(),"R","stratFold3D.R",sep="/"))
source(paste(getwd(),"R","penint3D.R",sep="/"))
source(paste(getwd(),"R","plotfolds.R",sep="/"))


fun<-fm.H
regdat<-H.df
contVar<-sm2D.lst[1:11]
targetVar<-"Humus"


strat<-stratfold3d(targetVar="Humus",regdat=H.df,folds=6,cent=3,dimensions="2D",IDs=TRUE,sum=TRUE)
flist<-strat$folds
plotfolds(strat,"Humus")
penint3D(fun=fm.H,regdat=H.df,lambda=10^seq(10,-2,length=100),contVar=contVar,int=TRUE,flist=flist)

#================= test for hierNet modification ===================================================================
source(paste(getwd(),"R","funcs.R",sep="/"))

H.df[,c("altitude",sm2D.lst[1:11])]<-apply(H.df[,c("altitude",sm2D.lst[1:11])],2,scale)

mm <- model.matrix(fm.GSIF.H ,H.df)[,-1] #fm.int.lm.As, # fm.GSIF.int.lm.As
nzv <- nearZeroVar(mm)
names(data.frame(mm)[, nzv])
if(sum(nzv) != 0){mm <- mm[, -nzv]}else{mm <- mm}

colnames(mm)<-gsub("\\(altitude,.df.=.4\\)","",colnames(mm))

H <- H.df$Humus

X <- hierNet::compute.interactions.c(mm,diagonal=FALSE)
head(X)


(columns.to.keep <- c(paste0(names(data.frame(mm)), ":",names(data.frame(mm))), X %>%
                        colnames() %>%
                        grep("altitude" %>% as.character(),., value = T)))

X[,colnames(X) %ni% columns.to.keep ] <- 0
fit = hierNet.path(mm,H, zz = X,diagonal=FALSE)
fitcv = hierNet.cv(fit,mm,H)

print(fitcv)
plot(fitcv)
# matrix of estimated interaction coefficients
(d <- fit$th[,,which(fitcv$lamhat.1se==fit$lamlist)])
# main effect estimates
fit$bp[,which(fitcv$lamhat.1se==fit$lamlist), drop = F] - fit$bn[,which(fitcv$lamhat.1se==fit$lamlist), drop = F]

# lamlist
# kvadratni clanovi
# izlaz iz hiernet
# 



hnet<-hierNet(x=mm[,c(1:19)],y=H,lam=0.3)#,diagonal=FALSE)#,zz=mm[,20:34])

nn<-names(data.frame(compute.interactions.c(mm[,c(1:19)],diagonal = FALSE)))
intind<-which(nn %in% grep("altitude", nn, value = TRUE))
head(compute.interactions.c(mm[,c(1:19)],diagonal = FALSE)[,intind])

interactions<-(compute.interactions.c(mm[,c(1:19)],diagonal = FALSE)[,intind])

hierNet.varimp(hnet,x=HTrans.df[,1:13],y=H)

hnet<-hierNet(x=mm[,c(1:19)],y=H,lam=0.3,diagonal=FALSE,zz=interactions)











