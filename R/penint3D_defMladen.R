
# solution for including all polynomial terms in interactions.

# predint3D is function for making prediction based on interaction approach.
# fun - function
# profs - Soil Profile Collection
# cogrids - SpatialPixelsDataFrame with covariates.
# hier - logical, does the hierarchial principle need to be honor
# pred - logocal, if TRUE , prediction on grids will be made
# lambda - grid of lambda (regularization parameter for lasso)
# deg - degree of polynomial depth function
# fold - number of folds
# cent - number of centers for k-means clustering
# int - logical, if TRUE , interaction will be included in model
# preProc - logical, if TRUE, covariates will be scaled and centered
# chunk - number of rows in prediction matrix (subselection needed to speedup the prediction computation)
# cores - number of cores

# sparsereg3D
# fun-base.model
# profs-profiles
# cov.grids
# poly.deg
# cv.folds
# cent-num.means
# interactions
# depth.fun-remove
# preproc

As.fun <- as.formula(paste("SOM ~", paste(c(CovAbb,"depth"), collapse="+")))
SOM.fun <- as.formula(paste("log(SOM) ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"depth"), collapse="+")))
pH.fun <- as.formula(paste("pH ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"depth"), collapse="+")))

base.model=SOM.fun; profiles=bor.profs; cov.grids=gridmaps.sm2D;poly.deg=3; num.folds=5;num.means=5;interactions=TRUE

sparsereg3Dncv<-function(base.model, profiles, cov.grids, hier=FALSE, lambda=seq(0,5,0.1), poly.deg=3, num.folds=5,num.means=3,interactions=TRUE,preproc=TRUE,seed=321){
  
  "%ni%" <- Negate("%in%")
  
  if(hier==TRUE){interactions=TRUE}


  #======= all columns of interest ===============================
  target <- all.vars(base.model)[1] # target.variable
  ## check if all covariates are available:
  if(sum(names(cov.grids) %in% all.vars(base.model)[-1])==0){
    stop("None of the covariates in the 'formulaString' matches the names in the 'covariates' object")
  }
  #===============================================================
  p4s=proj4string(profiles)
  profiles <- join(profiles@horizons,data.frame(data.frame(profiles@site),data.frame(profiles@sp)),type="inner")
  
  #============== Names ==========================================
  coord.names <- tail(names(profiles),2); 
  pr.names <- head(names(profiles),3);
  
  #============== Adding altitude and hdepth =====================
  profiles$depth <- - (profiles$Top / 100 + (( profiles$Bottom - profiles$Top ) / 2) / 100)
  profiles$hdepth<-profiles$Bottom - profiles$Top
  
  #============== Overlay ========================================
  profiles <- profiles[complete.cases(profiles[,c("ID",coord.names,"hdepth","depth",target)]),c("ID",coord.names,"hdepth","depth",target)]
  coordinates(profiles) <- ~ x + y
  proj4string(profiles) <- p4s
  profiles <- spTransform(profiles, proj4string(cov.grids))
  
  ov <- over(profiles, cov.grids)
  
  factor.names <- ov %>% subset(., select=which(sapply(., is.factor))) %>% names()
  
  for(i in factor.names){
    ov[,i] <- factor(ov[,i])
  }
  
  #======== prepare regression matrix: ===========================
  profiles <- cbind(as.data.frame(profiles), ov[,which(names(ov) %in% c(all.vars(base.model)))])
  profiles <- profiles[complete.cases(profiles[,all.vars(base.model)]),] 
  
  #cov.data <- profiles[,c(all.vars(base.model)[-1])]
  
#================== polynomial ======================================
  
  if(poly.deg > 1){
    profiles<-cbind(profiles,poly(profiles$depth,poly.deg,raw=TRUE,simple=TRUE)[,-1])
    names(profiles)<-c(names(profiles)[1:(length(names(profiles))-(poly.deg-1))],(paste("depth",c(2:poly.deg),sep="")))
  }

#===================== Dummy coding =======================================
  
  dummy.par <- dummyVars(as.formula(paste("~", paste(names(profiles), collapse="+"))),profiles,levelsOnly=FALSE)
  profiles <- predict(dummy.par, newdata = profiles)
  colnames(profiles)<-gsub( "\\_|/|\\-|\"|\\s" , "." , colnames(profiles) )
#==========================================================================
  
  if (interactions == TRUE){
    
    X <- hierNet::compute.interactions.c(as.matrix(profiles[,-c(1:6)]),diagonal=FALSE)
    
    if(poly.deg > 1){ columns.to.keep <- colnames(X[,do.call(c,lapply(strsplit(colnames(X),":"), function(x) x[2] %in% c("depth",paste("depth",c(2:poly.deg),sep="")) & x[1] %ni% c("depth",paste("depth",c(2:poly.deg),sep=""))))]) } else {
                      columns.to.keep <- (X %>% as.data.frame() %>% subset(., select=grep("depth", names(.), value=TRUE)) %>% colnames())
                      }
      
    if(hier == TRUE) { X[,colnames(X) %ni% columns.to.keep ] <- 0 
                       profiles <- cbind(profiles,X)
                     }else{
                            profiles <- cbind(profiles,X[,colnames(X) %in% columns.to.keep])
                          }

    } 

  
  
#===================== Dummy coding =======================================

  if(preproc==TRUE){

     st.par <- as.data.frame(profiles[,-c(1:6)]) %>% preProcess(.,method=c("center", "scale"))
     profiles[,-c(1:6)] <- predict(st.par,newdata = profiles[,-c(1:6)])
     
     nzv.par<-preProcess(profiles[,-c(1:6)],method=c("nzv"))
     profiles[,-c(1:6)]<-as.data.frame(predict(nzv.par,profiles[,-c(1:6)]))

  }


#====================== Columns to keep and cov.data ==================================================== 
  
  if (interactions==TRUE){
    if (hier!=TRUE){
          X <- hierNet::compute.interactions.c(as.matrix(cov.data),diagonal=FALSE)
      
          if(depth.fun!="linear"){ columns.to.keep <- colnames(X[,do.call(c,lapply(strsplit(colnames(X),":"), function(x) x[2] %in% c("altitude",paste("poly",c(2:poly.deg),sep="")) & x[1] %ni% c("altitude",paste("poly",c(2:poly.deg),sep=""))))]) } else {
                                   columns.to.keep <- (X %>% as.data.frame() %>% subset(., select=grep("altitude", names(.), value=TRUE)) %>% colnames())
          }
          
          
          #colnames(X[,do.call(c,lapply(strsplit(colnames(X),":"), function(x) x[2] %in% c("altitude",paste("poly",c(2:poly.deg),sep="")) & x[1] %ni% c("altitude",paste("poly",c(2:poly.deg),sep=""))))])
            
          if(depth.fun!="linear"){ cov.data <- cbind(cov.data,X[,colnames(X) %in% columns.to.keep])} else {
                                   cov.data <- cbind(cov.data,X[,colnames(X) %in% columns.to.keep])
                                 }
      
                   }  else { X <- hierNet::compute.interactions.c(as.matrix(cov.data),diagonal=FALSE)
      
                   if(depth.fun!="linear"){ columns.to.keep <- colnames(X[,do.call(c,lapply(strsplit(colnames(X),":"), function(x) x[2] %in% c("altitude",paste("poly",c(2:poly.deg),sep="")) & x[1] %ni% c("altitude",paste("poly",c(2:poly.deg),sep=""))))]) } else {
                     columns.to.keep <- (X %>% as.data.frame() %>% subset(., select=grep("altitude", names(.), value=TRUE)) %>% colnames())
                   }
      
                            X[,colnames(X) %ni% columns.to.keep ] <- 0
                            
                            cov.data <- cov.data
      
                           } 
    
                   } else {
      
                           cov.data <- cov.data
    
             }
    
#=====================================================================================================
  
  regmat.def <- as.data.frame(cbind(regmat[, target],cov.data))
  names(regmat.def) <- c(target, names(regmat.def[-1]))
  target.min <- min(regmat.def[,1])
#=====================================================================================================

  if(!hier){

    allData <- cbind(regmat.def, regmat[,c(pr.names[1], coord.names,"hdepth")])
    allData <- plyr::rename(allData, replace=c("x" = "longitude", "y" = "latitude"))
    
    strat<-stratfold3d(targetVar=target,seed=seed,regdat=allData,folds=num.folds,cent=num.means,preProc=FALSE,dimensions="2D",IDs=TRUE,sum=TRUE)
    flist<-strat$folds
    cv.folds<-stratfold3d(targetVar=target,seed=seed,regdat=allData,folds=num.folds,cent=num.means,preProc=FALSE,dimensions="2D",IDs=FALSE,sum=TRUE)$folds
    
    results <- data.frame(lambda = rep(NA,length(flist)+1),RMSE=rep(NA,length(flist)+1),Rsquared=rep(NA,length(flist)+1))
    coef.list = as.list(rep(NA,length(flist)))
    pred <- data.frame()
    train.results <- as.list(rep(NA,length(flist)))
    test.results <- as.list(rep(NA,length(flist)))
    
    for(i in 1:length(flist)){
      ind <- which(allData$ID %in% do.call(c,flist[-i]))
      TrainData <- allData[ind,]
      Train.ID.index <- flist[-i]
      TestData <- allData[-ind,]
      Test.ID.Index <- flist[i]
      
      #===========================================================================

      ff<-stratfold3d(target, TrainData, folds=num.folds, seed=seed ,cent=num.means, dimensions="2D", sum=TRUE, IDs=TRUE, preProc=FALSE)$folds
      
      folds.in.list <- as.list(rep(NA,length(ff)))
      names(folds.in.list) <- paste("fold",c(1:length(ff)),sep = "")
      foldid <- rep(NA,dim(TrainData)[1])
      for(j in 1:length(ff)){
        folds.in.list[[j]]<-which(TrainData$ID %in% ff[[j]])
        foldid[folds.in.list[[j]]]<-j
      }
      
      #==========================================================================================
      
      TrainData <- TrainData[,1:(dim(TrainData)[2]-4)] # ovo mora da se menja...ne sme da stoji tri vec moraju da se uklone poslednje tri kolone lepse
      TestData <- TestData[,1:(dim(TestData)[2]-4)]
      #==================== Lasso ===============================================================
      
      lasso.mod <- cv.glmnet(as.matrix(TrainData[,-1]),TrainData[,1],alpha=1,lambda=lambda,foldid=foldid,type.measure="mse")
      lasso.pred <- predict(lasso.mod,s=lasso.mod$lambda.min,newx=as.matrix(TestData[,-1]))
      lasso.pred <- pmax(lasso.pred,target.min/3)
      
      train.pred <- predict(lasso.mod,s=lasso.mod$lambda.min,newx=as.matrix(TrainData[,-1]))
      train.pred <- pmax(train.pred,target.min/3)
      train.pred <- data.frame(predicted=as.numeric(train.pred))
      train.res <- data.frame(ID=allData[ind,"ID"], observed=TrainData[,1], predicted=train.pred,longitude=allData[ind,"longitude"],latitude=allData[ind,"latitude"],altitude=regmat[ind,"altitude"])
      train.results[[i]] <- train.res
      
      test.res <- data.frame(ID=allData[-ind,"ID"], observed=TestData[,1],predicted=as.numeric(lasso.pred),longitude=allData[-ind,"longitude"],latitude=allData[-ind,"latitude"],altitude=regmat[-ind,"altitude"])
      test.results[[i]] <- test.res
      
      obs.pred <- data.frame(obs=TestData[,1],pred=as.numeric(lasso.pred))
      coef.list[[i]] <- predict(lasso.mod,type="coefficients",s=lasso.mod$lambda.min)
      dfresults <- data.frame(lambda=lasso.mod$lambda.min,RMSE=defaultSummary(obs.pred)[1],Rsquared=defaultSummary(obs.pred)[2])
      results[i,] <- dfresults
      pred<-rbind(pred,obs.pred)
    }
    
    results[length(flist)+1,] <- c(NA,RMSE=defaultSummary(pred)[1],Rsquared=defaultSummary(pred)[2])
    coef.mat <- do.call(cbind,coef.list)
    out <- list(measure=results,coef=coef.mat,obs.pred.df = pred,folds=cv.folds,train.results=train.results,test.results = test.results)
    return(out)
    
  } else {
    

    allData <- cbind(regmat.def,X, regmat[,c(pr.names[1], coord.names,"hdepth")])
    allData <- plyr::rename(allData, replace=c("x" = "longitude", "y" = "latitude"))
    
    strat<-stratfold3d(targetVar=target,seed=123,regdat=allData,folds=num.folds,cent=3,dimensions="2D",IDs=TRUE,sum=TRUE)
    flist<-strat$folds
    cv.folds<-stratfold3d(targetVar=target,seed=seed,regdat=allData,folds=num.folds,cent=num.means,preProc=FALSE,dimensions="2D",IDs=FALSE,sum=TRUE)$folds
    
    results <- data.frame(lambda = rep(NA,length(flist)+1),RMSE=rep(NA,length(flist)+1),Rsquared=rep(NA,length(flist)+1))
    coef.list = as.list(rep(NA,length(flist)))
    pred <- data.frame()
    train.results <- as.list(rep(NA,length(flist)))
    test.results <- as.list(rep(NA,length(flist)))
    
    for(i in 1:length(flist)){
      ind<-which(allData$ID %in% do.call(c,flist[-i])) # izdvajaje indeksa instanci koje pribadaju foldovima za trening
      TrainData<-as.data.frame(do.call(cbind,allData[ind,])) # izdvajanje trening podacima
      Train.ID.index<-flist[-i] # Izdvajanje dela liste foldova sa  trening podacima
      TestData<-as.data.frame(do.call(cbind,allData[-ind,])) #Izdvajanje test podataka
      Test.ID.Index<-flist[i] # Izdvajanje dela liste foldova sa test podacima
      
      #=======================================
      
      ff<-stratfold3d(target,TrainData, folds=num.folds,sum=TRUE,IDs=TRUE,preProc=FALSE)$folds
      
      folds.in.list <- as.list(rep(NA,length(ff)))
      names(folds.in.list) <- paste("fold",c(1:length(ff)),sep = "")
      foldid <- rep(NA,dim(TrainData)[1])
      for(j in 1:length(ff)){
        folds.in.list[[j]]<-which(TrainData$ID %in% ff[[j]])
        foldid[folds.in.list[[j]]]<-j
      }
      
      
      trainx <- as.matrix(TrainData[,colnames(cov.data)]) 
      testx <- as.matrix(TestData[,colnames(cov.data)])
      trainzz <- as.matrix(TrainData[,colnames(X)])
      testzz <- as.matrix(TestData[,colnames(X)])
      trainy <- (TrainData[,target])
      testy <- (TestData[,target])
      
      fit = hierNet.path(trainx,trainy, zz = trainzz,diagonal=FALSE,strong=TRUE,trace=0)
      fitcv = hierNet.cv(fit, trainx, trainy, folds=folds.in.list, trace=0)
      fit.def <- hierNet(trainx,trainy, zz = trainzz,diagonal=FALSE,strong=TRUE,lam=fit$lamlist[which(fitcv$lamhat==fit$lamlist)])
      fit.pred <- predict(fit.def,newx=testx,newzz = testzz)
      fit.pred <- pmax(fit.pred,target.min/3)
      ie <- as.matrix(fit$th[,,which(fitcv$lamhat==fit$lamlist)][,length(colnames(cov.data))])
      me<-fit$bp[,which(fitcv$lamhat==fit$lamlist), drop = F] - fit$bn[,which(fitcv$lamhat==fit$lamlist), drop = F]
      rbind(ie,me)
      
      train.pred <- predict(fit.def,newx=trainx,newzz = trainzz)
      train.pred <- pmax(train.pred,target.min/3)
      train.pred <- data.frame(predicted=as.numeric(train.pred))
      train.res <- data.frame(ID=TrainData[,"ID"],observed=TrainData[,1], predicted=train.pred,longitude=TrainData[,"longitude"], latitude=TrainData[,"latitude"], altitude=regmat[ind,"altitude"])
      train.results[[i]] <- train.res
      
      test.res <- data.frame(ID=TestData[,"ID"], observed=TestData[,1], predicted=as.numeric(fit.pred), longitude=TestData[,"longitude"], latitude=TestData[,"latitude"], altitude=regmat[-ind,"altitude"])
      test.results[[i]] <- test.res
      
      obs.pred <- data.frame(obs=testy,pred=fit.pred)
      coef.list[[i]]<-rbind(ie,me)
      dfresults<-data.frame(lambda=fitcv$lamhat,RMSE=defaultSummary(obs.pred)[1],Rsquared=defaultSummary(obs.pred)[2])
      results[i,]<-dfresults
      pred<-rbind(pred,obs.pred)
    }
    results[length(flist)+1,]<-c(NA,RMSE=defaultSummary(pred)[1],Rsquared=defaultSummary(pred)[2])
    coef.mat<-do.call(cbind,coef.list)
    out<-list(measure=results,coef=coef.mat,obs.pred.df = pred,folds=cv.folds,train.results=train.results,test.results=test.results)
    return(out)
  }
  
}
