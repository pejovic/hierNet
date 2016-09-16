
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

As.fun <- as.formula(paste("As ~", paste(c(CovAbb,"depth"), collapse="+")))
SOM.fun <- as.formula(paste("SOM ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"depth"), collapse="+")))
pH.fun <- as.formula(paste("pH ~", paste(c(CovAbb[-which(CovAbb %in% c("CD","DD"))],"depth"), collapse="+")))

base.model=SOM.fun; profiles=bor.profs; cov.grids=gridmaps.sm2D; poly.deg=3; num.folds=5;num.means=5;interactions=TRUE; hier=TRUE;preproc=TRUE;seed=321

pre.sparsereg3D <- function(base.model, profiles, cov.grids, hier=FALSE, lambda=seq(0,5,0.1), poly.deg=3, num.folds=5,num.means=3,interactions=TRUE,preproc=TRUE,seed=321){
  
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
  profiles <- profiles[complete.cases(profiles[,c("ID",coord.names,"hdepth","depth",target)]),c("ID",target,"hdepth",coord.names,"depth")]
  coordinates(profiles) <- ~ x + y
  proj4string(profiles) <- p4s
  profiles <- spTransform(profiles, proj4string(cov.grids))
  
  ov <- over(profiles, cov.grids)
  
  factor.names <- ov %>% subset(., select=which(sapply(., is.factor))) %>% names()
  
  for(i in factor.names){
    ov[,i] <- factor(ov[,i])
  }
  
  #======== prepare regression matrix: ===========================
  sp.cov.names <- names(ov[,which(names(ov) %in% c(all.vars(base.model)))])
  profiles <- cbind(as.data.frame(profiles), ov[,sp.cov.names])
  profiles <- profiles[complete.cases(profiles[,all.vars(base.model)]),c("ID",target,"hdepth",coord.names, sp.cov.names, "depth")] 
  
#================== polynomial ======================================
  
  if(poly.deg > 1){
    profiles <- cbind(profiles,poly(profiles$depth,poly.deg,raw=TRUE,simple=TRUE)[,-1])
    names(profiles) <- c(names(profiles)[1:(length(names(profiles))-(poly.deg-1))],(paste("depth",c(2:poly.deg),sep="")))
  }

#===================== Dummy coding =======================================
  
  dummy.par <- dummyVars(as.formula(paste("~", paste(names(profiles), collapse="+"))),profiles,levelsOnly=FALSE)
  profiles <- predict(dummy.par, newdata = profiles)
  colnames(profiles) <- gsub( "\\_|/|\\-|\"|\\s" , "." , colnames(profiles) )
  main.effects <- colnames(profiles)[-which(colnames(profiles) %in% c("ID","hdepth",target,coord.names))]
#==========================================================================
  
  if (interactions == TRUE){
    
    X <- hierNet::compute.interactions.c(as.matrix(profiles[,-c(which(colnames(profiles) %in% c("ID",target,"hdepth",coord.names)))]),diagonal=FALSE)
    
    if(poly.deg > 1){ columns.to.keep <- colnames(X[,do.call(c,lapply(strsplit(colnames(X),":"), function(x) x[2] %in% c("depth",paste("depth",c(2:poly.deg),sep="")) & x[1] %ni% c("depth",paste("depth",c(2:poly.deg),sep=""))))]) } else {
                      columns.to.keep <- (X %>% as.data.frame() %>% subset(., select=grep("depth", names(.), value=TRUE)) %>% colnames())
                      }
      
    if(hier == TRUE) { X[,colnames(X) %ni% columns.to.keep ] <- 0 
                       profiles <- cbind(profiles,X)
                     }else{
                            profiles <- cbind(profiles,X[,colnames(X) %in% columns.to.keep])
                          }

    } 

  
  
#===================== preprocessing =======================================

  if(preproc==TRUE){

     st.par <- as.data.frame(profiles) %>% subset(., select = c(main.effects,columns.to.keep)) %>% preProcess(.,method=c("center", "scale"))
     profiles <- predict(st.par,newdata = profiles)
     
     #nzv.par<-preProcess(profiles[,-c(1:6)],method=c("nzv"))
     #profiles[,-c(1:6)]<-as.data.frame(predict(nzv.par,profiles[,-c(1:6)]))

  }


#======================== Stratification ====================================================
  
  profiles <- as.data.frame(profiles)
  profiles <- plyr::rename(profiles, replace=c("x" = "longitude", "y" = "latitude"))

  strat<-stratfold3d(targetVar=target,seed=seed,regdat=profiles,folds=num.folds,cent=num.means,preProc=FALSE,dimensions="2D",IDs=TRUE,sum=TRUE)
  profile.folds.list <- strat$folds
  obs.folds.list <- stratfold3d(targetVar=target, seed=seed,regdat=profiles, folds=num.folds,cent=num.means, preProc=FALSE, dimensions="2D", IDs=FALSE,sum=TRUE)$folds
#=====================================================================================================
  
  out <- list(profiles = profiles, model = list(target = target, main.effects = main.effects, poly.deg = poly.deg, hier = hier, depth.interactions = columns.to.keep, all.interactions = colnames(X)), num.folds = num.folds, num.means = num.means, preproc = list(dummy.par = dummy.par, st.par = st.par, columns.to.keep = columns.to.keep), folds = list(profile.folds.list = profile.folds.list, obs.folds.list = obs.folds.list, depth.summary = strat[[3]], target.summmary = strat[[4]]))
  return(out)
  
}  
  

pre.som <- pre.sparsereg3D(base.model = SOM.fun, hier = TRUE, profiles=bor.profs, cov.grids = gridmaps.sm2D)
  


sparsereg3D.ncv <- function(sparsereg = pre.som, lambda = seq(0,5,0.1), seed = 321){
  
  flist <- sparsereg$folds$profile.folds.list
  cv.folds <- sparsereg$folds$obs.folds.list
  profiles <- sparsereg$profiles
  num.folds <- sparsereg$num.folds
  num.means <- sparsereg$num.means
  
  target <- sparsereg$model$target
  target.min <- min(profiles[,target])
  hier = sparsereg$model$hier
  main.effects = sparsereg$model$main.effects
  poly.deg = sparsereg$model$poly.deg
  columns.to.keep = sparsereg$model$interactions
  all.interactions = sparsereg$model$all.interactions
  
  results <- data.frame(lambda = rep(NA,length(flist)+1), RMSE = rep(NA,length(flist)+1), Rsquared = rep(NA,length(flist)+1))
  coef.list = as.list(rep(NA,length(flist)))
  test.prediction <- data.frame()
  train.results <- as.list(rep(NA,length(flist)))
  test.results <- as.list(rep(NA,length(flist)))
  
  for(i in 1:length(flist)){
    train.obs.ind <- which(profiles$ID %in% do.call(c, flist[-i]))
    train.data <- profiles[train.obs.ind,]
    train.profiles.ids <- flist[-i]
    test.data <- profiles[-train.obs.ind,]
    test.profiles.ids <- flist[i]
    
    # =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = 
    
    inner.partitioning <- stratfold3d(target, train.data, folds = num.folds, seed = seed ,cent = num.means, dimensions = "2D", sum = TRUE, IDs = TRUE, preProc = FALSE)$folds
    
    new.folds <- as.list(rep(NA,length(inner.partitioning)))
    names(new.folds) <- paste("fold",c(1:length(inner.partitioning)),sep = "")
    
    which.new.fold <- rep(NA,dim(train.data)[1])
    
    for(j in 1:length(inner.partitioning)){
      new.folds[[j]]<-which(train.data$ID %in% inner.partitioning[[j]])
      which.new.fold[new.folds[[j]]]<-j
    }
    
    # =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = 
    if(!hier){
      
      train.data <- train.data %>% subset(., select = c(target, main.effects, columns.to.keep)) 
      test.data <- test.data %>% subset(., select = c(target, main.effects, columns.to.keep))
      # =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  Lasso  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = 
      
      lasso.cv <- cv.glmnet(as.matrix(train.data[,-1]), train.data[,1], alpha = 1,lambda = lambda, foldid = which.new.fold, type.measure = "mse")
      test.pred <- predict(lasso.cv, s = lasso.cv$lambda.min, newx = as.matrix(test.data[,-1]))
      test.pred <- pmax(test.pred, target.min/3)
      
      train.pred <- predict(lasso.cv,s = lasso.cv$lambda.min, newx = as.matrix(train.data[,-1]))
      train.pred <- pmax(train.pred,target.min/3)
      train.pred <- data.frame(predicted = as.numeric(train.pred))
      train.per.fold <- data.frame(ID = profiles[train.obs.ind,"ID"], observed = train.data[,1], predicted = train.pred,longitude = profiles[train.obs.ind,"longitude"], latitude = profiles[train.obs.ind,"latitude"], depth = profiles[train.obs.ind,"depth"])
      train.results[[i]] <- train.per.fold
      
      test.per.fold <- data.frame(ID = profiles[-train.obs.ind,"ID"], observed = test.data[,1], predicted = as.numeric(test.pred), longitude = profiles[-train.obs.ind,"longitude"], latitude = profiles[-train.obs.ind,"latitude"], depth = profiles[-train.obs.ind,"depth"])
      test.results[[i]] <- test.per.fold
      
      test.obs.pred <- data.frame(obs = test.data[,1], pred = as.numeric(test.pred))
      coef.list[[i]] <- predict(lasso.cv, type = "coefficients", s = lasso.cv$lambda.min)
      fold.measures <- data.frame(lambda = lasso.cv$lambda.min, RMSE = defaultSummary(test.obs.pred)[1], Rsquared = defaultSummary(test.obs.pred)[2])
      results[i,] <- fold.measures
      test.prediction <- rbind(test.prediction,test.obs.pred)
      
    }else{
      
      trainx <- as.matrix(train.data[,main.effects]) 
      testx <- as.matrix(test.data[,main.effects])
      trainzz <- as.matrix(train.data[,all.interactions])
      testzz <- as.matrix(test.data[,all.interactions])
      trainy <- train.data[,target]
      testy <- test.data[,target]
      
      fit = hierNet.path(trainx, trainy, zz = trainzz, diagonal=FALSE, strong=TRUE, trace=0, stand.main = FALSE)
      fitcv = hierNet.cv(fit, trainx, trainy, folds = new.folds, trace=0)
      fit.def <- hierNet(trainx, trainy, zz = trainzz, diagonal=FALSE, strong=TRUE, lam=fit$lamlist[which(fitcv$lamhat==fit$lamlist)])
      test.pred <- predict(fit.def,newx=testx,newzz = testzz)
      test.pred <- pmax(test.pred,target.min/3)
      ie <- as.matrix(fit$th[,,which(fitcv$lamhat==fit$lamlist)][,length(main.effects)])
      me<-fit$bp[,which(fitcv$lamhat==fit$lamlist), drop = F] - fit$bn[,which(fitcv$lamhat==fit$lamlist), drop = F]

      train.pred <- predict(fit.def,newx=trainx,newzz = trainzz)
      train.pred <- pmax(train.pred,target.min/3)
      train.pred <- data.frame(predicted=as.numeric(train.pred))
      train.per.fold <- data.frame(ID = train.data[,"ID"], observed = train.data[,1], predicted = train.pred, longitude = train.data[,"longitude"], latitude = train.data[,"latitude"], depth = profiles[train.obs.ind,"depth"])
      train.results[[i]] <- train.per.fold
      
      test.per.fold <- data.frame(ID = test.data[,"ID"], observed = test.data[,1], predicted = as.numeric(test.pred), longitude = test.data[,"longitude"], latitude = test.data[,"latitude"], depth = profiles[-train.obs.ind,"depth"])
      test.results[[i]] <- test.per.fold
      
      test.obs.pred <- data.frame(obs=testy,pred=test.pred)
      coef.list[[i]] <- rbind(ie,me)
      fold.measures <- data.frame(lambda=fitcv$lamhat,RMSE=defaultSummary(test.obs.pred)[1],Rsquared=defaultSummary(test.obs.pred)[2])
      results[i,] <- fold.measures
      test.prediction<-rbind(test.prediction,test.obs.pred)
      }
    }
  
  results[length(flist)+1,] <- c(NA, RMSE = defaultSummary(test.prediction)[1], Rsquared = defaultSummary(test.prediction)[2])
  coef.mat <- do.call(cbind, coef.list)
  out <- list(accuracy.measures = results, coefficients = coef.mat, test.prediction = test.prediction, folds = cv.folds, train.results = train.results, test.results = test.results)
  return(out)
  
}

pre.som <- pre.sparsereg3D(base.model = SOM.fun, hier = TRUE, profiles=bor.profs, cov.grids = gridmaps.sm2D, poly.deg = 3)    
sp.reg.som <- sparsereg3D.ncv(sparsereg = pre.som, lambda = seq(0,5,0.1), seed = 321)
    
sp.reg.som$accuracy.measures
    
    
    
    
    
    
    
    
    
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
