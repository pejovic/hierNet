geo <- bor.geo
cogrids <- gridmaps.sm2D

penint3Dpred<-function(fun, geo, cogrids, flist,hier=TRUE,lambda=seq(0,5,0.1),int=TRUE,depth.fun=list("linear","nspline","poly"),df=4,deg=3,preProc=TRUE){
  
  "%ni%" <- Negate("%in%")
  
  if(hier==TRUE){int=TRUE}
  
  if(int == FALSE){hier = FALSE}
  
  
  #======= all columns of interest:
  methodid <- all.vars(fun)[1]
  seln <- names(cogrids) %in% all.vars(fun)[-1]
  xyn <- attr(cogrids@bbox, "dimnames")[[1]]
  ## check if all covariates are available:
  if(sum(!is.na(seln))==0){
    stop("None of the covariates in the 'formulaString' matches the names in the 'covariates' object")
  }
  
  #======== prepare regression matrix:
  ov <- over(x=cogrids, y=geo, method=methodid, var.type = "numeric")
  if(nrow(ov)==0|is.null(ov$observedValue)) {
    warning("The 'over' operations resulted in an empty set. Check 'methodid' column.")
  }
  
  ## geostats only possible with numeric variables:
  ov[,methodid] = as.numeric(ov$observedValue)
  ov$observedValue = NULL
  
  ov<-ov[-which(is.na(ov[,methodid])),]
  
  modmat<-ov[,c(all.vars(fun)[-1])]
  
  #================== depth.fun query ======================================
  
  if(depth.fun!="linear"){
    modmat<-cbind(modmat,poly(modmat$altitude,deg,raw=TRUE,simple=TRUE)[,-1])
    names(modmat)<-c(names(modmat)[1:(length(names(modmat))-(deg-1))],(paste("poly",c(2:deg),sep="")))
  }
  
  #==========================================================================
  
  if(preProc==TRUE){
    
    cont.par <- ov[,c("sampleid",all.vars(fun)[-1])] %>% ddply(.,.(sampleid),function(x) head(x,1)) %>% subset(.,select=-c(sampleid,altitude),drop=FALSE) %>% subset(., select=which(sapply(., is.numeric))) %>% preProcess(.,method=c("center", "scale"))
    
    if(depth.fun =="linear"){ alt.par <-  modmat %>% subset(.,select=altitude) %>% preProcess(.,method=c("center", "scale"))}
    
    
    if(depth.fun!="linear"){ poly.par <- modmat[,-which(names(modmat) %in% c(all.vars(fun)[-1]))] %>% preProcess(.,method=c("center", "scale"))}
    
    if(depth.fun!="linear"){ dummy.par <- dummyVars(as.formula(paste("~", paste(names(modmat), collapse="+"))),modmat,levelsOnly=FALSE)} else {
      dummy.par <- dummyVars(as.formula(paste("~", paste(names(modmat), collapse="+"))),modmat,levelsOnly=FALSE)
    }
    
    if(depth.fun!="linear"){ modmat <- predict(cont.par,newdata = modmat) %>% predict(alt.par,newdata = .) %>% predict(poly.par,newdata = .) %>% predict(dummy.par,newdata = .)} else {
      modmat <- predict(cont.par,newdata = modmat) %>% predict(alt.par,newdata = .) %>% predict(dummy.par,newdata = .)
    }
    
    nzv.par<-preProcess(modmat,method=c("nzv"))
    modmat<-as.data.frame(predict(nzv.par,modmat))
    
    names(modmat)<-gsub( "\\_|/|\\-|\"|\\s" , "." , names(modmat) )
    
  }
  
  #====================== Columns to keep and modmat ==================================================== 
  
  if (int==TRUE){
    if (hier!=TRUE){
      X <- hierNet::compute.interactions.c(as.matrix(modmat),diagonal=FALSE)
      
      if(depth.fun!="linear"){ columns.to.keep <- (X %>% as.data.frame() %>% subset(., select=grep("altitude", names(.), value=TRUE)) %>% subset(., select=(names(.) %ni% c(grep("poly", names(.), value=TRUE)))) %>% colnames())} else {
        columns.to.keep <- (X %>% as.data.frame() %>% subset(., select=grep("altitude", names(.), value=TRUE)) %>% colnames())
      }
      
      if(depth.fun!="linear"){ modmat <- cbind(modmat,X[,colnames(X) %in% columns.to.keep])} else {
        modmat <- cbind(modmat,X[,colnames(X) %in% columns.to.keep])
      }
      
    }  else { X <- hierNet::compute.interactions.c(as.matrix(modmat),diagonal=FALSE)
    
    if(depth.fun!="linear"){ columns.to.keep <- (X %>% as.data.frame() %>% subset(., select=grep("altitude", names(.), value=TRUE)) %>% subset(., select=(names(.) %ni% c(grep("poly", names(.), value=TRUE)))) %>% colnames())} else {
      columns.to.keep <- (X %>% as.data.frame() %>% subset(., select=grep("altitude", names(.), value=TRUE)) %>% colnames())
    }
    
    X[,colnames(X) %ni% columns.to.keep ] <- 0
    
    modmat <- modmat
    
    } 
    
  } else {
    
    modmat <- modmat
    
  }
  
  #=====================================================================================================
  
  regdat<-cbind(As=ov$As,modmat)
  
  #=====================================================================================================
  
  if(!hier){
    
    allData <- cbind(regdat,ov[,c("sampleid","longitude","latitude")])
    names(allData)[names(allData) == 'sampleid'] <- 'ID'
    allData$ID <- as.integer(allData$ID)
    allData$hdepth <- rep(1,dim(allData)[1])
    
    flist<-stratfold3d(methodid,allData, folds=5,sum=TRUE,IDs=TRUE,preProc=FALSE)$folds
    
    #results <- data.frame(lambda = rep(NA,length(flist)+1),RMSE=rep(NA,length(flist)+1),Rsquared=rep(NA,length(flist)+1))
    #coef.list = as.list(rep(NA,length(strat)))
    #pred <- data.frame()
    

      folds.in.list <- as.list(rep(NA,length(flist)))
      names(folds.in.list) <- paste("fold",c(1:length(flist)),sep = "")
      foldid <- rep(NA,dim(allData)[1])
      
      for(j in 1:length(flist)){
        folds.in.list[[j]]<-which(allData$ID %in% flist[[j]])
        foldid[folds.in.list[[j]]]<-j
      }
      
      #==========================================================================================
      
      allData <- allData[,1:(dim(allData)[2]-4)] # ovo mora da se menja...ne sme da stoji tri vec moraju da se uklone poslednje tri kolone lepse

      #==================== Lasso ===============================================================
      
      lasso.mod <- cv.glmnet(as.matrix(allData[,-1]),allData[,1],alpha=1,lambda=lambda,foldid=foldid,type.measure="mse")
      lasso.pred <- predict(lasso.mod,s=lasso.mod$lambda.min,newx=as.matrix(allData[,-1]))
      cvm <- sqrt(lasso.mod$cvm[which(lasso.mod$lambda==lasso.mod$lambda.min)])
      obs.pred <- data.frame(obs=allData[,1],pred=as.numeric(lasso.pred))
      coef.list <- predict(lasso.mod,type="coefficients",s=lasso.mod$lambda.min)
      dfresults <- data.frame(lambda=lasso.mod$lambda.min,RMSE=defaultSummary(obs.pred)[1],Rsquared=defaultSummary(obs.pred)[2])
      pred<-rbind(pred,obs.pred)


      regrez<-list(preProc=list(cont.par=cont.par,alt.par=alt.par,dummy.par=dummy.par, nzv.par=nzv.par),rezults=dfresults,coefficients=coef.list,pred=data.frame(obs=allData[,1],pred=lasso.pred,rez=allData[,1]-lasso.pred))
      return(regrez)
    
  } else {
    
    allData <- cbind(regdat,X,ov[,c("sampleid","longitude","latitude")])
    names(allData)[names(allData) == 'sampleid'] <- 'ID'
    allData$ID <- as.numeric(allData$ID)
    allData$hdepth <- rep(1,dim(allData)[1])
    
    flist<-stratfold3d(methodid,allData, folds=5,sum=TRUE,IDs=TRUE,preProc=FALSE)$folds
    
    #results <- data.frame(lambda=rep(NA,length(flist)+1),RMSE=rep(NA,length(flist)+1),Rsquared=rep(NA,length(flist)+1))
    #coef.list = as.list(rep(NA,length(flist)))
    #pred <- data.frame()
    
    folds.in.list <- as.list(rep(NA,length(flist)))
    names(folds.in.list) <- paste("fold",c(1:length(flist)),sep = "")
    foldid <- rep(NA,dim(allData)[1])
    
    for(j in 1:length(flist)){
      folds.in.list[[j]]<-which(allData$ID %in% flist[[j]])
      foldid[folds.in.list[[j]]]<-j
    }
    
    #==========================================================================================
    
    allData <- allData[,1:(dim(allData)[2]-4)] # ovo mora da se menja...ne sme da stoji tri vec moraju da se uklone poslednje tri kolone lepse
    
      
      allx <- as.matrix(allData[,colnames(modmat)]) 
      allzz <- as.matrix(allData[,colnames(X)])
      ally <- (allData[,methodid])

      fit = hierNet.path(allx,ally, zz = allzz,diagonal=FALSE,strong=TRUE,trace=0)
      fitcv = hierNet.cv(fit, allx, ally, folds=folds.in.list, trace=0)
      fit.def <- hierNet(allx,ally, zz = allzz,diagonal=FALSE,strong=TRUE,lam=fit$lamlist[which(fitcv$lamhat==fit$lamlist)])
      fit.pred <- predict(fit.def,newx=allx,newzz = allzz)
      
      ie <- as.matrix(fit$th[,,which(fitcv$lamhat==fit$lamlist)][,length(colnames(modmat))])
      me<-fit$bp[,which(fitcv$lamhat==fit$lamlist), drop = F] - fit$bn[,which(fitcv$lamhat==fit$lamlist), drop = F]
      rbind(ie,me)
      
      obs.pred <- data.frame(obs=ally,pred=fit.pred)
      coef.list<-data.frame(cov.name=colnames(allx),me,ie)
      dfresults<-data.frame(lambda=fitcv$lamhat,RMSE=defaultSummary(obs.pred)[1],Rsquared=defaultSummary(obs.pred)[2])

    }
    
      regrez<-list(preProc=list(cont.par=cont.par,alt.par=alt.par,dummy.par=dummy.par, nzv.par=nzv.par),rezults=dfresults,coefficients=coef.list,pred=data.frame(obs=ally,pred=fit.pred,rez=ally-fit.pred))
      return(regrez)
      
      
  }

  

