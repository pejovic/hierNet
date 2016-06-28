# This function perform the 3D kriging leave-one-out cross-validation. It was perform for the each fold in the outer loop in the nested cross-validation. 
# Input data: 

# fun - function, the function which relates the response variable and covariates. 
# reg.ncv - output from the penint3Dncv. This output contains slots, among others, with the train and test data/results. 
# The train data/results consists the n-fold lists with the observed and predicted values from the final trend model on the train data. These folds serve for 3D variogram modeling for each      

#fun=SOM.fun; reg.pred=IntP; profs=bor.profs; model = FALSE; krige=FALSE ;v.cutoff=60; v.width=3; v.vgm = vgm(1.5, "Gau", 10, 0); sp.cutoff=4000; sp.width=420; sp.vgm = vgm(15, "Sph", 2000, 5)

krige3Dpred <- function(fun, reg.pred, profs, model = FALSE, krige=FALSE ,v.cutoff, v.width, v.vgm = vgm(1.5, "Gau", 10, 0), sp.cutoff, sp.width, sp.vgm = vgm(15, "Sph", 2000, 5) ) {
  
  "%ni%" <- Negate("%in%")
  #n.cv <- length(reg.ncv$train.results)
  
  #var1D.list <- as.list(rep(NA,n.cv))
  #var2D.list <- as.list(rep(NA,n.cv))
  #var3D.list <- as.list(rep(NA,n.cv))
  
  #profs <- bor.profs; fun <- SOM.fun; cogrids <- gridmaps.sm2D
  
  tv <- all.vars(fun)[1]
  #seln <- names(cogrids) %in% all.vars(fun)[-1]
  #xyn <- attr(cogrids@bbox, "dimnames")[[1]]
  #===============================================================
  tv.data <- join(profs@horizons,data.frame(data.frame(profs@site),data.frame(profs@sp)),type="inner")
  
  #============== Names ==========================================
  sp.names <- tail(names(tv.data),2); pr.names <- head(names(tv.data),3); methods <- names(tv.data)[names(tv.data) %ni% c(sp.names,pr.names)]
  
  #============== Adding altitude and hdepth =====================
  tv.data$altitude <- - (tv.data$Top / 100 + (( tv.data$Bottom - tv.data$Top ) / 2) / 100)
  tv.data$hdepth<-tv.data$Bottom - tv.data$Top
  
  prof.sp <- tv.data[complete.cases(tv.data[,c("ID",sp.names,"altitude",tv)]),c("ID",sp.names,"hdepth","Top","Bottom","altitude",tv)]
  prof.sp<-ddply(prof.sp,.(ID))
  
    train.results <- reg.pred$summary$pred # results obtained on training data i.e. final model fitted on training data in order to obtain the residuals for variogram fitting
    #train.results$residual <- train.results$observed-train.results$predicted # residuals from final model fitted on training data
    
    #======================= Preaparing profiles with residuals for vertical variogram fitting ======================================
    res.profs<-data.frame(ID=train.results$ID,residuals=train.results$residual) # 
    
    ind <- which(prof.sp[,"ID"] %in% unique(res.profs$ID)) # extracting indices for ith fold
    
    res.profs <- cbind(res.profs, prof.sp[ind,c("Top","Bottom","x","y","altitude")])
    
    res.list<-list(sites=ddply(res.profs,.(ID),summarize,y=y[1],x=x[1]),horizons=res.profs[,c("ID","Top","Bottom","residuals")])
    res.profs<-join(res.list$sites,res.list$horizons, type="inner")
    
    res.profs.data<-res.profs
    
    depths(res.profs) <- ID ~ Top + Bottom
    site(res.profs) <- ~x+y
    coordinates(res.profs)<-~x+y
    proj4string(res.profs)<- proj4string(profs)
    
    res.data <- train.results
    
    spline.profs <- mpspline(res.profs, "residuals",d = t(c(0,5,15,30,60,80)))
    spline.res <- spline.profs$var.1cm
    
    cm <- rep(1:80,dim(spline.res)[2])
    a <-ddply(res.data,.(ID),summarize ,longitude=longitude[1],latitude=latitude[1])
    b <- a[rep(seq_len(nrow(a)), each=80),]
    bm <- cbind(b,residual=as.numeric(spline.res),altitude=cm)
    bm <- na.omit(bm)
    
    #====================== Fitting 1D variogram =============================================================
    
    if(!model){
      
      vc <- variogram(residual~1, data=bm, locations=~longitude+latitude+altitude,cutoff=v.cutoff,width=v.width)
      
      var1D.list <- qplot(dist, gamma, data = vc, geom = c("point"),color=id)+theme_bw() + 
        theme(legend.position="none")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
        labs(x="Distance",y="Semivariance")
      #names(var1D.list[[i]])<-paste("P",i,sep="")
      
    }else{
      vc <- variogram(residual~1, data=bm, locations=~longitude+latitude+altitude,cutoff=v.cutoff,width=v.width)
      vc.fit = fit.variogram(vc, v.vgm)
      vr.line <- (variogramLine(vc.fit, maxdist = max(vc$dist)))
      
      var1D.list <- qplot(dist, gamma, data = vc, geom = c("point"),color=id)+theme_bw() + 
        theme(legend.position="none")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
        labs(x="Distance",y="Semivariance") + addlinetoplot(vr.line, varx = "dist", vary = "gamma",color="#56B4E9",sajz=0.65)  
      #names(var1D.list[[i]])<-paste("P",i,sep="")
    }
    
    #================= fitting 2D variogram ============================================================
    
    
    res.pp <- ddply(res.data, .(ID),summarize,x=longitude[1],y=latitude[1],altitude=altitude[1],residual=residual[1],observed=observed[1])
    coordinates(res.pp) <- ~x+y+altitude
    proj4string(res.pp)<- proj4string(reg.pred$prediction[[1]])
    
    #bubble(res.pp,"residual")
    
    if(!model){
      sv <- variogram(residual~1, data=res.pp,cutoff=sp.cutoff,width= sp.width) # 480 za As
      var2D.list <- qplot(dist, gamma, data = sv, geom = c("point"),color=id)+theme_bw() + 
        theme(legend.position="none")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
        labs(x="Distance",y="Semivariance")
      #names(var2D.list[[i]])<-paste("P",i,sep="")
    }else{
      sv <- variogram(residual~1, data=res.pp,cutoff=sp.cutoff,width= sp.width) # 480 za As
      sv.fit =  fit.variogram(sv, sp.vgm) #SOM
      sv.line <- (variogramLine(sv.fit, maxdist = max(sv$dist)))
      
      var2D.list <- qplot(dist, gamma, data = sv, geom = c("point"),color=id)+theme_bw() + 
        theme(legend.position="none")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
        labs(x="Distance",y="Semivariance") + addlinetoplot(sv.line, varx = "dist", vary = "gamma",color="#56B4E9",sajz=0.65)  
      #names(var2D.list[[i]])<-paste("P",i,sep="")
    }
    
    
    #================ Fitting 3D variogram =============================================================
    coordinates(res.data) <- ~ longitude + latitude + altitude
    proj4string(res.data) <- CRS(gk_7)
    
    if(!model){
      
      vr.def <- variogram(residual~1,res.data, cutoff=sp.cutoff,width=sp.width)  # 420 za SOM
      var3D.list <- qplot(dist, gamma, data = vr.def, geom = c("point"),color=id)+theme_bw() + 
        theme(legend.position="none")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
        labs(x="Distance",y="Semivariance")
      #names(var3D.list[[i]])<-paste("P",i,sep="")
      
    }else {
      
      vr.def <- variogram(residual~1,res.data, cutoff=sp.cutoff,width=sp.width)  # 420 za SOM
      vgr.def <- fit.variogram(vr.def, vgm(sv.fit$psill[2], as.character(sv.fit$model[2]), sv.fit$range[2],sv.fit$psill[1] ,anis=c(0,0,0,1,(vc.fit$range/100)[2]/sv.fit$range[2])))
      vr.def.line <- (variogramLine(vgr.def, maxdist = max(vr.def$dist)))
      
      var3D.list <- qplot(dist, gamma, data = vr.def, geom = c("point"),color=id)+theme_bw() + 
        theme(legend.position="none")+theme(axis.text=element_text(size=12),axis.title=element_text(size=14,face="bold"))+
        labs(x="Distance",y="Semivariance") + addlinetoplot(vr.def.line, varx = "dist", vary = "gamma",color="#56B4E9",sajz=0.65)
      #names(var3D.list[[i]])<-paste("P",i,sep="")
    }
    
    #================================= 3D kriging =============================================================================
    if(krige){
      
      
      for( i in 1:length(reg.pred$prediction)){
        reg.pred$prediction[[i]]$res.pred <- krige(residual ~ 1, res.pp, model = vgr.def, newdata=reg.pred$prediction[[i]])$var1.pred
        reg.pred$prediction[[i]]$final.prediction <- reg.pred$prediction[[i]]$pred+reg.pred$prediction[[i]]$res.pred
      }
      
      
      }  

  if(!model){
  out <- list(prediction = reg.pred$prediction, var1D = var1D.list, var2D = var2D.list, var3D = var3D.list)}else{
    out <- list(prediction = reg.pred$prediction, var1D = var1D.list, var2D = var2D.list, var3D = var3D.list, anis=(vc.fit$range/100)[2]/sv.fit$range[2])
  }
  return(out)
}