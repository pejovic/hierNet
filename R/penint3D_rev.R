geo <- bor.geo
cogrids <- gridmaps.sm2D

penint3Drev<-function(fun, geo, cogrids, flist,hier=TRUE,lambda=seq(0,5,0.1),deg=3,int=TRUE,depth.fun=list("linear","nspline","poly"),df=4,deg=3,preProc=TRUE){
  
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

    allData<-cbind(regdat,ov[,c("sampleid","longitude","latitude")])
    names(allData)[names(allData) == 'sampleid'] <- 'ID'
    allData$ID <- as.numeric(allData$ID)

    results <- data.frame(lambda = rep(NA,length(flist)+1),RMSE=rep(NA,length(flist)+1),Rsquared=rep(NA,length(flist)+1))
    coef.list = as.list(rep(NA,length(strat)))
    pred <- data.frame()
    
    for(i in 1:length(flist)){
      ind <- which(allData$ID %in% do.call(c,flist[-i]))
      TrainData <- allData[ind,]
      Train.ID.index <- flist[-i]
      TestData <- allData[-ind,]
      Test.ID.Index <- flist[i]
      
      #===========================================================================
      #unique.df<-ddply(allData,.(ID),here(summarize),tv=mean(eval(parse(text=methodid))),longitude=longitude[1],latitude=latitude[1])

      k.list <- as.list(rep(NA,length(Train.ID.index)))
      names(k.list)<-paste("k",c(1:length(k.list)),sep="")
      
      TrainData$hdepth <- rep(1,dim(TrainData)[1])
      
      ff<-stratfold3d(methodid,TrainData, folds=length(flist),sum=TRUE,IDs=TRUE,preProc=FALSE)$folds
      
      folds.in.list <- as.list(rep(NA,length(ff)))
      names(folds.in.list) <- paste("fold",c(1:length(ff)),sep = "")
      foldid <- rep(NA,dim(TrainData)[1])
      for(j in 1:length(ff)){
        folds.in.list[[j]]<-which(TrainData$ID %in% ff[[j]])
        foldid[folds.in.list[[j]]]<-j
      }
      
      #==========================================================================================
      
      TrainData <- TrainData[,1:(dim(TrainData)[2]-4)] # ovo mora da se menja...ne sme da stoji tri vec moraju da se uklone poslednje tri kolone lepse
      TestData <- TestData[,1:(dim(TestData)[2]-3)]
      #==================== Lasso ===============================================================
      
      lasso.mod <- cv.glmnet(as.matrix(TrainData[,-1]),TrainData[,1],alpha=1,lambda=lambda,foldid=foldid,type.measure="mse")
      lasso.pred <- predict(lasso.mod,s=lasso.mod$lambda.min,newx=as.matrix(TestData[,-1]))
      obs.pred <- data.frame(obs=TestData[,1],pred=as.numeric(lasso.pred))
      coef.list[[i]] <- predict(lasso.mod,type="coefficients",s=lasso.mod$lambda.min)
      dfresults <- data.frame(lambda=lasso.mod$lambda.min,RMSE=defaultSummary(obs.pred)[1],Rsquared=defaultSummary(obs.pred)[2])
      results[i,] <- dfresults
      pred<-rbind(pred,obs.pred)
    }
    
    results[length(flist)+1,] <- c(NA,RMSE=defaultSummary(pred)[1],Rsquared=defaultSummary(pred)[2])
    coef.mat <- do.call(cbind,coef.list)
    out <- list(measure=results,coef=coef.mat)
    return(out)
    
  } else {
    

    allData <- cbind(regdat,X,ov[,c("sampleid","longitude","latitude")])
    names(allData)[names(allData) == 'sampleid'] <- 'ID'
    allData$ID <- as.numeric(allData$ID)
    
    results <- data.frame(lambda=rep(NA,length(flist)+1),RMSE=rep(NA,length(flist)+1),Rsquared=rep(NA,length(flist)+1))
    coef.list = as.list(rep(NA,length(flist)))
    pred <- data.frame()
    
    for(i in 1:length(flist)){
      ind<-which(allData$ID %in% do.call(c,flist[-i])) # izdvajaje indeksa instanci koje pribadaju foldovima za trening
      TrainData<-as.data.frame(do.call(cbind,allData[ind,])) # izdvajanje trening podacima
      Train.ID.index<-flist[-i] # Izdvajanje dela liste foldova sa  trening podacima
      TestData<-as.data.frame(do.call(cbind,allData[-ind,])) #Izdvajanje test podataka
      Test.ID.Index<-flist[i] # Izdvajanje dela liste foldova sa test podacima
      
      #=======================================
      #unique.df<-ddply(allData,.(ID),here(summarize),tv=mean(eval(parse(text=methodid))),longitude=longitude[1],latitude=latitude[1])
      
      k.list <- as.list(rep(NA,length(Train.ID.index)))
      names(k.list)<-paste("k",c(1:length(k.list)),sep="")
      
      TrainData$hdepth <- rep(1,dim(TrainData)[1])
      
      ff<-stratfold3d(methodid,TrainData, folds=length(flist),sum=TRUE,IDs=TRUE,preProc=FALSE)$folds
      
      folds.in.list <- as.list(rep(NA,length(ff)))
      names(folds.in.list) <- paste("fold",c(1:length(ff)),sep = "")
      foldid <- rep(NA,dim(TrainData)[1])
      for(j in 1:length(ff)){
        folds.in.list[[j]]<-which(TrainData$ID %in% ff[[j]])
        foldid[folds.in.list[[j]]]<-j
      }
      

      trainx <- as.matrix(TrainData[,colnames(modmat)]) 
      testx <- as.matrix(TestData[,colnames(modmat)])
      trainzz <- as.matrix(TrainData[,colnames(X)])
      testzz <- as.matrix(TestData[,colnames(X)])
      trainy <- (TrainData[,methodid])
      testy <- (TestData[,methodid])
      
      fit = hierNet.path(trainx,trainy, zz = trainzz,diagonal=FALSE,strong=TRUE,trace=0)
      fitcv = hierNet.cv(fit, trainx, trainy, folds=folds.in.list, trace=0)
      fit.def <- hierNet(trainx,trainy, zz = trainzz,diagonal=FALSE,strong=TRUE,lam=fit$lamlist[which(fitcv$lamhat==fit$lamlist)])
      fit.pred <- predict(fit.def,newx=testx,newzz = testzz)
      
      ie <- as.matrix(fit$th[,,which(fitcv$lamhat==fit$lamlist)][,length(colnames(modmat))])
      me<-fit$bp[,which(fitcv$lamhat==fit$lamlist), drop = F] - fit$bn[,which(fitcv$lamhat==fit$lamlist), drop = F]
      rbind(ie,me)
      
      obs.pred <- data.frame(obs=testy,pred=fit.pred)
      coef.list[[i]]<-rbind(ie,me)
      dfresults<-data.frame(lambda=fitcv$lamhat,RMSE=defaultSummary(obs.pred)[1],Rsquared=defaultSummary(obs.pred)[2])
      results[i,]<-dfresults
      pred<-rbind(pred,obs.pred)
    }
    results[length(flist)+1,]<-c(NA,RMSE=defaultSummary(pred)[1],Rsquared=defaultSummary(pred)[2])
    coef.mat<-do.call(cbind,coef.list)
    out<-list(measure=results,coef=coef.mat)
    return(out)
  }
  
}
