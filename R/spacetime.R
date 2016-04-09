#########################################
#                                       #
# 3D+T demo of the cookfarm data        #
#                                       #
# The demo script has been prepared by: #
# Benedikt Graeler                      #
# Institute for Geoinformatics          #
# ben.graeler@uni-muenster.de           #
# ifgi.de/graeler                       #
#                                       #
##########################################

load("HumusPrediction2732016.RData")

# time <- proc.time() # the full script takes 10-15 minutes
## krigeST
library(spacetime)
library(gstat)
library(xts)
library(sp)

FFL.p <- FFL$summary$pred
FFP.p <- FFP$summary$pred
TFL.p <- TFL$summary$pred
TFP.p <- TFP$summary$pred
TTL.p <- TTL$summary$pred
TTP.p <- TTP$summary$pred




names(FFL.p)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")
names(FFP.p)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")
names(TFL.p)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")
names(TFP.p)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")
names(TTL.p)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")
names(TTP.p)<-c("ID","x","y","hdepth","altitude","observed","predicted","residual")


res.profs<-data.frame(ID=TFL.p$ID,FFL.res=FFL.p$residual,FFP.res=FFP.p$residual,TFL.res=TFL.p$residual,TFP.res=TFP.p$residual,TTL.res=TTL.p$residual,TTP.res=TTP.p$residual)

borind<-ddply(bor,.(ID))

ind<-which(res.profs$ID %in% borind[,"ID"])

res.profs<-cbind(res.profs,borind[ind,c("Top","Bottom","Soil.Type","x","y","altitude")])

res.list<-list(sites=ddply(res.profs,.(ID),summarize,y=y[1],x=x[1],Soil.Type=Soil.Type[1]),horizons=res.profs[,c("ID","Top","Bottom","FFL.res","FFP.res","TFL.res","TFP.res","TTL.res","TTP.res")])
res.profs<-join(res.list$sites,res.list$horizons, type="inner")


depths(res.profs) <- ID ~ Top + Bottom
site(res.profs) <- ~x+y
coordinates(res.profs)<-~x+y
proj4string(res.profs)<- CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m")

# again, this time along 1-cm slices, computing quantiles
agg <- slab(res.profs, fm= ~ FFL.res + FFP.res + TFL.res + TFP.res + TTL.res + TTP.res,slab.structure=seq(0,80,5))

## see ?slab for details on the default aggregate function
head(agg)

xyplot(top ~ p.q50 | variable, data=agg, ylab='Depth',
       xlab='median bounded by 25th and 75th percentiles',
       lower=agg$p.q25, upper=agg$p.q75, ylim=c(80,-2),
       panel=panel.depth_function,
       alpha=0.25, sync.colors=TRUE,
       par.settings=list(superpose.line=list(col='RoyalBlue', lwd=2)),
       prepanel=prepanel.depth_function,
       cf=agg$contributing_fraction, cf.col='black', cf.interval=5, 
       layout=c(6,1), strip=strip.custom(bg=grey(0.8)),
       scales=list(x=list(tick.number=4, alternating=3, relation='free'))
)

head(FFL.p)


#============== Compare profiles =========================================

FFL.p$x <- FFL.p$x+5000
FFP.p$x <- FFP.p$x+10000
TFL.p$x <- TFL.p$x+15000
TFP.p$x <- TFP.p$x+20000
TTL.p$x <- TTL.p$x+25000
TTP.p$x <- TTP.p$x+30000

FFL.p$ID <- FFL.p$ID
FFP.p$ID <- FFP.p$ID+500
TFL.p$ID <- TFL.p$ID+1000
TFP.p$ID <- TFP.p$ID+1500
TTL.p$ID <- TTL.p$ID+2000
TTP.p$ID <- TTP.p$ID+2500

c.p <- rbind(FFL.p,FFP.p,TFL.p,TFP.p,TTL.p,TTP.p)

str(c.p)

#c.profs<-data.frame(ID=TFL.p$ID,FFL.res=FFL.p$residual,FFP.res=FFP.p$residual,TFL.res=TFL.p$residual,TFP.res=TFP.p$residual,TTL.res=TTL.p$residual,TTP.res=TTP.p$residual)

c.profs<-cbind(c.p,borind[ind,c("Top","Bottom")][rep(ind,6),])
c.profs$model <- factor(c(rep("FFL",455),rep("FFP",455),rep("TFL",455),rep("TFP",455),rep("TTL",455),rep("TTP",455)))

str(c.profs)

c.list<-list(sites=ddply(c.profs,.(ID),summarize,y=y[1],x=x[1],model=model[1]),horizons=c.profs[,c("ID","Top","Bottom","hdepth","altitude","observed","predicted","residual")])
c.profs<-join(c.list$sites,c.list$horizons, type="inner")


depths(c.profs) <- ID ~ Top + Bottom
site(c.profs) <- ~ model + x + y
coordinates(c.profs)<-~x+y
proj4string(c.profs)<- CRS("+proj=utm +zone=34 +ellps=GRS80 +towgs84=0.26901,0.18246,0.06872,-0.01017,0.00893,-0.01172,0.99999996031 +units=m")

agg <- slab(c.profs[which(c.profs@site$model %in% c('FFL')), ], fm= ~ residual, slab.structure=seq(0,50,5))


idx <- which(c.profs@site$model %in% c('FFL', 'TFL'))

FFL.and.TFL <- slab(c.profs[which(c.profs@site$model %in% c('FFL', 'TFL')), ], fm= ~ residual, slab.structure=seq(0,50,5))
FFL.and.TTL <- slab(c.profs[which(c.profs@site$model %in% c('TTL')), ], fm= ~ residual, slab.structure=seq(0,50,5))


g <- make.groups(FFL=agg, TTL=FFL.and.TTL)

xyplot(top ~ p.q50 | variable, groups=which, data=g, ylab='Depth',
       xlab='median bounded by 25th and 75th percentiles',
       lower=g$p.q25, upper=g$p.q75, ylim=c(42,-2),
       panel=panel.depth_function,
       alpha=0.25, sync.colors=TRUE, cf=g$contributing_fraction, cf.interval=10, 
       par.settings=list(superpose.line=list(col=c('RoyalBlue', 'Red4'), lwd=2, lty=c(1,2))),
       prepanel=prepanel.depth_function,
       layout=c(1,1), strip=strip.custom(bg=grey(0.8)),
       scales=list(x=list(tick.number=4, alternating=3, relation='free')),
       auto.key=list(columns=2, lines=TRUE, points=FALSE)
)







FFL$dd <- -FFL$altitude/5*1000*24*3600 # 0,5cm -> days in secs
TFL$dd <- -TFL$altitude/5*1000*24*3600 # 0,5cm -> days in secs


# create a spacetime data structure
FFLSt <- STIDF(SpatialPoints((FFL[,c("Easting","Northing")])),
                 as.POSIXct(FFL$dd, tz = "GMT", origin = "2005-01-01"),
                 FFL[,c("ID","residual")])

str(FFLSt)

FFLSt <- as(FFLSt,"STFDF")

# calculate the horizointal-vertical sample variogram
svgmCDepthAsTime <- variogramST(residual~1, FFLSt, tunit="day", tlags=seq(0,40,20) ,width=700, cutoff=2400)
                                tlags=0:3,  assumeRegular = T,
                                width=300, cutoff=1600, na.omit=F)

svgmCDepthAsTime$timelag <- svgmCDepthAsTime$timelag*0.005
plot(svgmCDepthAsTime, wireframe=T)


#=============================== TFL ==================================
# create a spacetime data structure
TFLSt <- STIDF(SpatialPoints((TFL[,c("Easting","Northing")])),
               as.POSIXct(TFL$dd, tz = "GMT", origin = "2005-01-01"),
               TFL[,c("ID","residual")])

str(TFLSt)
summary(TFLSt)

TFLSt <- as(TFLSt,"STFDF")

# calculate the horizointal-vertical sample variogram
svgmCDepthAsTime <- variogramST(residual~1, TFLSt,tuni, tlags=seq(10,60,20), assumeRegular = T, na.omit=TRUE ,width=600, cutoff=2400)
plot(svgmCDepthAsTime)#, wireframe=T)
plot(svgmCDepthAsTime, wireframe=T)

svgmCDepthAsTime$timelag <- svgmCDepthAsTime$timelag*0.005
plot(svgmCDepthAsTime)#, wireframe=T)

TFL.p$a<-rep(1,dim(TFL.p)[1])

ddply(TFL.p,.(ID),summarize, res.range=residual[2]-residual[1],alt.range=altitude[2]-altitude[1])

TFL.p <- TFL.p[-(which(TFL.p$ID %in% c(103,108))),]

v <- variogram(residual~1, data=TFL.p, locations=~x+y+altitude,cutoff=0.7,width=0.1)
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




# fit a variogram surface, assuming a metric model
fvgmCDepthAsTime <- fit.StVariogram(svgmCDepthAsTime, 
                                    vgmST("metric",
                                          joint = vgm(0.02, "Sph", 150, 0.02),
                                          stAni=100))
attr(fvgmCDepthAsTime, "optim.output")$value

plot(svgmCDepthAsTime, fvgmCDepthAsTime, wireframe=T, all=T,
     ylab="depth diff [m]", scales=list(arrows=F))




# add coordinates to the readings data.frame
matchIDs <- match(cookfarm$readings$SOURCEID, cookfarm$profiles$SOURCEID)
measC <- cookfarm$readings[,c("Date","Port1C","Port2C","Port3C","Port4C","Port5C")]
measC$Easting <-  cookfarm$profiles$Easting[matchIDs]
measC$Northing <-  cookfarm$profiles$Northing[matchIDs]
measC$DoY <- as.numeric(format(as.POSIXct(measC$Date), format="%j"))

# deseasonalize per depth
deseason <- function(data, var, par) {
  data[,var] - (par[1]+par[2]*sin(((par[3]+data[,"DoY"])/365)*2*pi))
}

varVec <- c("Port1C","Port2C","Port3C","Port4C","Port5C")
opt <- matrix(NA,length(varVec), 3)
for (dVar in 1:length(varVec)) {
  opt[dVar,] <- optim(c(0,1,0), 
                      function(par) {
                        sqrt(mean(deseason(measC, varVec[dVar], par)^2, na.rm = T))
                      })$par
  measC[,paste("resid",varVec[dVar],sep="")] <- deseason(measC, varVec[dVar], opt[dVar,])  
}

# re-order data.frame in a long format
measC <- data.frame(Easting=rep(measC$Easting, 5),
                    Northing=rep(measC$Northing, 5),
                    altitude=rep(1:5*-0.3, each=nrow(measC)),
                    Date=rep(measC$Date,5),
                    DoY=rep(measC$DoY, 5),
                    C=c(measC$Port1C, 
                        measC$Port2C,
                        measC$Port3C,
                        measC$Port4C,
                        measC$Port5C),
                    resid=c(measC$residPort1C, 
                            measC$residPort2C,
                            measC$residPort3C,
                            measC$residPort4C,
                            measC$residPort5C))

# abuse time as depth axis and encode depth levels as days
measC$dd <- -measC$altitude/3*10*24*3600 # depth level -> days in secs

# move the 3D data set by 2 km every day avoiding any overlap between days
measC$EastingShift <- measC$Easting - (as.numeric(measC$Date)-14976)*2000 

# create a spacetime data structure
measCSt <- STIDF(SpatialPoints(measC[,c("EastingShift","Northing")]),
                 as.POSIXct(measC$dd, tz = "GMT", origin = "2005-01-01"),
                 measC[,c("resid","C")])
measCSt <- as(measCSt,"STFDF")

# calculate the horizointal-vertical sample variogram
svgmCDepthAsTime <- variogramST(resid~1, measCSt, 
                                tlags=0:4, assumeRegular = T, 
                                width=75, cutoff=375, na.omit=F)
svgmCDepthAsTime$timelag <- svgmCDepthAsTime$timelag*0.3
plot(svgmCDepthAsTime, wireframe=T)

# fit a variogram surface, assuming a metric model
fvgmCDepthAsTime <- fit.StVariogram(svgmCDepthAsTime, 
                                    vgmST("metric",
                                          joint = vgm(0.02, "Sph", 150, 0.02),
                                          stAni=100))
attr(fvgmCDepthAsTime, "optim.output")$value

plot(svgmCDepthAsTime, fvgmCDepthAsTime, wireframe=T, all=T,
     ylab="depth diff [m]", scales=list(arrows=F))


fvgmCDepthAsTime$stAni # 1 m in depth corresponds to ~486 m in the plain

# build a re-scaled spatio-temporal data set with 3D spatial component
measC$scaleAlt <- measC$altitude * fvgmCDepthAsTime$stAni

measCSt <- STIDF(SpatialPoints(measC[,c("Easting","Northing","scaleAlt")]),
                 as.POSIXct(measC$Date, tz = "GMT", origin = "2005-01-01"),
                 measC[,c("resid", "C")])
measCSt <- as(measCSt, "STFDF")

# 3D space-time sample variogram
svgmC3DT <- variogramST(resid~1, measCSt, 
                        tlags=0:9, assumeRegular=T, 
                        width=75, cutoff=375)

plot(svgmC3DT, wireframe=T, all=T, zlim=c(0,1.5),
     xlab = "3D distance [m]", scales=list(arrows=F),
     col.regions=bpy.colors())

# fit a variogram
fvgmC3DT <- fit.StVariogram(svgmC3DT, 
                            vgmST("sumMetric",
                                  space = vgm(1, "Exp", 50, 0),
                                  time = vgm(1, "Exp", 9, 0),
                                  joint = vgm(1, "Sph", 20, 0),
                                  stAni=0.3),
                            lower=c(0,0,0,
                                    0,0,0,
                                    0,0,0,
                                    0),
                            method="L-BFGS-B")
attr(fvgmC3DT, "optim.output")$value # 0.031

plot(svgmC3DT, fvgmC3DT, wireframe=T, all=T, 
     scales=list(arrows=F), xlab = "3D distance [m]")

# predict depth level wise for July:
dLevel <- 1 # 1..5

## gridded data:
grid10m <- cookfarm$grids
gridded(grid10m) <- ~x+y
proj4string(grid10m) <- CRS(cookfarm$proj4string)

predSTF <- as.data.frame(cbind(coordinates(grid10m),
                               c(1:5*(-0.3)*fvgmCDepthAsTime$stAni)[dLevel]))
colnames(predSTF) <- c("Easting", "Northing", "altitude")
coordinates(predSTF) <- ~Easting+Northing+altitude
predSTF <- STF(predSTF, measCSt@time[547:577])

# drop NAs
measCSt <- as(measCSt, "STSDF")
noNA <- !is.na(measCSt@data$resid)
measCSt@index <- measCSt@index[noNA,]
measCSt@data <- measCSt@data[noNA,]

# do the prediction for the first 7 days of July 2012
preds <- NULL
pars <- opt[dLevel,]
for(day in 1:7) { # day <- 1
  cat(paste("day:", day), "\n")
  pred <- krigeST(resid~1, measCSt[,(541+day):(551+day)], predSTF[,day,drop=F],
                  fvgmC3DT)$var1.pred
  doy <- as.numeric(format(index(predSTF[,day,drop=F]@time), format="%j"))
  pred <- data.frame(residPred=pred, 
                     pred=pred + (pars[1]+pars[2]*sin(((pars[3]+doy)/365)*2*pi)))
  preds <- rbind(preds, pred)
}

stCpred <- STFDF(SpatialPixels(SpatialPoints(coordinates(predSTF@sp)[,1:2])),
                 predSTF@time[1:7], preds)

stplot(stCpred[,,"pred"])
# proc.time() - time