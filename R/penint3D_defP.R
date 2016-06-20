
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

penint3DP<-function(fun, profs, cogrids,hier=FALSE,lambda=seq(0,5,0.1),deg=3,fold=5,cent=3,int=TRUE,depth.fun=list("linear","poly"),preProc=TRUE,seed=321){
  
  "%ni%" <- Negate("%in%")
  
  if(hier==TRUE){int=TRUE}
  
  if(int == FALSE){hier = FALSE}
  

  #======= all columns of interest ===============================
  tv <- all.vars(fun)[1]
  seln <- names(cogrids) %in% all.vars(fun)[-1]
  xyn <- attr(cogrids@bbox, "dimnames")[[1]]
  ## check if all covariates are available:
  if(sum(!is.na(seln))==0){
    stop("None of the covariates in the 'formulaString' matches the names in the 'covariates' object")
  }
  #===============================================================
  tv.data <- join(profs@horizons,data.frame(data.frame(profs@site),data.frame(profs@sp)),type="inner")
  
  #============== Names ==========================================
  sp.names <- tail(names(tv.data),2); pr.names <- head(names(tv.data),3); methods <- names(tv.data)[names(tv.data) %ni% c(sp.names,pr.names)]
  
  #============== Adding altitude and hdepth =====================
  tv.data$altitude <- - (tv.data$Top / 100 + (( tv.data$Bottom - tv.data$Top ) / 2) / 100)
  tv.data$hdepth<-tv.data$Bottom - tv.data$Top
  
  #============== Overlay ========================================
  prof.sp <- tv.data[complete.cases(tv.data[,c("ID",sp.names,"altitude",tv)]),c("ID",sp.names,"hdepth","altitude",tv)]
  # prod.sp <- plyr::rename(prof.sp, replace=c(paste(sp.names[1]) = "x", paste(sp.names[2]) = "y"))
  coordinates(prof.sp) <- ~ x + y
  proj4string(prof.sp) <- proj4string(profs)
  prof.sp <- spTransform(prof.sp, proj4string(cogrids))
  
  ov <- over(prof.sp, cogrids)# %>% cbind(ID=prof.sp@data[,"ID"],.)
  
  fnames <- ov %>% subset(., select=which(sapply(., is.factor))) %>% names()
  
  for(i in fnames){
    ov[,i] <- factor(ov[,i])
  }
  
  #======== prepare regression matrix: ===========================
  regmat<-cbind(as.data.frame(prof.sp), ov)
  regmat <- regmat[complete.cases(regmat[,all.vars(fun)[-1]]),] 
  
  if (sum(is.na(regmat[,tv])) != 0) {regmat<-regmat[-which(is.na(regmat[,tv])),]}
  
  modmat<-regmat[,c(all.vars(fun)[-1])]
  
#================== depth.fun query ======================================
  
  if(depth.fun!="linear"){
    modmat<-cbind(modmat,poly(modmat$altitude,deg,raw=TRUE,simple=TRUE)[,-1])
    names(modmat)<-c(names(modmat)[1:(length(names(modmat))-(deg-1))],(paste("poly",c(2:deg),sep="")))
  }

#==========================================================================
  
  if(preProc==TRUE){

    cont.par <- regmat[,c("ID",all.vars(fun)[-1])] %>% ddply(.,.(ID),function(x) head(x,1)) %>% subset(.,select=-c(ID,altitude),drop=FALSE) %>% subset(., select=which(sapply(., is.numeric))) %>% preProcess(.,method=c("center", "scale"))
    
                           {alt.par <-  modmat %>% subset(.,select=altitude) %>% preProcess(.,method=c("center", "scale"))}
                           #if(depth.fun =="linear")
    
    if(depth.fun!="linear"){ poly.par <- modmat %>% subset(.,select=paste("poly",c(2:deg),sep="")) %>% preProcess(.,method=c("center", "scale"))}
    
    if(depth.fun!="linear"){ dummy.par <- dummyVars(as.formula(paste("~", paste(names(modmat), collapse="+"))),modmat,levelsOnly=FALSE)} else {
                             dummy.par <- dummyVars(as.formula(paste("~", paste(names(modmat), collapse="+"))),modmat,levelsOnly=FALSE)
                           }
                              #%>% predict(alt.par,newdata = .) %>% predict(poly.par,newdata = .)
    #if(depth.fun!="linear"){ mm <- predict(cont.par,newdata = modmat) ; modmat <- model.matrix(as.formula(paste("~", paste(c(colnames(mm)), collapse="+"))),mm)[,-1]} else {
                             #mm <- predict(cont.par,newdata = modmat) %>% predict(alt.par,newdata = .) ; modmat <- model.matrix(as.formula(paste("~", paste(c(colnames(mm)), collapse="+"))),mm)[,-1]
                           #} 
    
    
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
      
          if(depth.fun!="linear"){ columns.to.keep <- colnames(X[,do.call(c,lapply(strsplit(colnames(X),":"), function(x) x[2] %in% c("altitude",paste("poly",c(2:deg),sep="")) & x[1] %ni% c("altitude",paste("poly",c(2:deg),sep=""))))]) } else {
                                   columns.to.keep <- (X %>% as.data.frame() %>% subset(., select=grep("altitude", names(.), value=TRUE)) %>% colnames())
          }
          
          
          #colnames(X[,do.call(c,lapply(strsplit(colnames(X),":"), function(x) x[2] %in% c("altitude",paste("poly",c(2:deg),sep="")) & x[1] %ni% c("altitude",paste("poly",c(2:deg),sep=""))))])
            
          if(depth.fun!="linear"){ modmat <- cbind(modmat,X[,colnames(X) %in% columns.to.keep])} else {
                                   modmat <- cbind(modmat,X[,colnames(X) %in% columns.to.keep])
                                 }
      
                   }  else { X <- hierNet::compute.interactions.c(as.matrix(modmat),diagonal=FALSE)
      
                   if(depth.fun!="linear"){ columns.to.keep <- colnames(X[,do.call(c,lapply(strsplit(colnames(X),":"), function(x) x[2] %in% c("altitude",paste("poly",c(2:deg),sep="")) & x[1] %ni% c("altitude",paste("poly",c(2:deg),sep=""))))]) } else {
                     columns.to.keep <- (X %>% as.data.frame() %>% subset(., select=grep("altitude", names(.), value=TRUE)) %>% colnames())
                   }
      
                            X[,colnames(X) %ni% columns.to.keep ] <- 0
                            
                            modmat <- modmat
      
                           } 
    
                   } else {
      
                           modmat <- modmat
    
             }
    
#=====================================================================================================
  
  regmat.def <- as.data.frame(cbind(regmat[, tv],modmat))
  names(regmat.def) <- c(tv, names(regmat.def[-1]))
  tv.min <- min(regmat.def[,1])
#=====================================================================================================

  if(!hier){

    allData <- cbind(regmat.def, regmat[,c(pr.names[1], sp.names,"hdepth")])
    allData <- plyr::rename(allData, replace=c("x" = "longitude", "y" = "latitude"))
    
    strat<-stratfold3d(targetVar=tv,seed=seed,regdat=allData,folds=fold,cent=cent,preProc=FALSE,dimensions="2D",IDs=TRUE,sum=TRUE)
    flist<-strat$folds
    
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

      ff<-stratfold3d(tv, TrainData, folds=fold, seed=seed ,cent=cent, dimensions="2D", sum=TRUE, IDs=TRUE, preProc=FALSE)$folds
      
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
      lasso.pred <- pmax(lasso.pred,tv.min/3)
      obs.pred <- data.frame(obs=TestData[,1],pred=as.numeric(lasso.pred))
      coef.list[[i]] <- predict(lasso.mod,type="coefficients",s=lasso.mod$lambda.min)
      dfresults <- data.frame(lambda=lasso.mod$lambda.min,RMSE=defaultSummary(obs.pred)[1],Rsquared=defaultSummary(obs.pred)[2])
      results[i,] <- dfresults
      pred<-rbind(pred,obs.pred)
    }
    
    results[length(flist)+1,] <- c(NA,RMSE=defaultSummary(pred)[1],Rsquared=defaultSummary(pred)[2])
    coef.mat <- do.call(cbind,coef.list)
    out <- list(measure=results,coef=coef.mat,obs.pred.df = pred)
    return(out)
    
  } else {
    

    allData <- cbind(regmat.def,X, regmat[,c(pr.names[1], sp.names,"hdepth")])
    allData <- plyr::rename(allData, replace=c("x" = "longitude", "y" = "latitude"))
    
    strat<-stratfold3d(targetVar=tv,seed=123,regdat=allData,folds=fold,cent=3,dimensions="2D",IDs=TRUE,sum=TRUE)
    flist<-strat$folds

    results <- data.frame(lambda = rep(NA,length(flist)+1),RMSE=rep(NA,length(flist)+1),Rsquared=rep(NA,length(flist)+1))
    coef.list = as.list(rep(NA,length(strat)))
    pred <- data.frame()
    

    for(i in 1:length(flist)){
      ind<-which(allData$ID %in% do.call(c,flist[-i])) # izdvajaje indeksa instanci koje pribadaju foldovima za trening
      TrainData<-as.data.frame(do.call(cbind,allData[ind,])) # izdvajanje trening podacima
      Train.ID.index<-flist[-i] # Izdvajanje dela liste foldova sa  trening podacima
      TestData<-as.data.frame(do.call(cbind,allData[-ind,])) #Izdvajanje test podataka
      Test.ID.Index<-flist[i] # Izdvajanje dela liste foldova sa test podacima
      
      #=======================================

      ff<-stratfold3d(tv,TrainData, folds=fold,sum=TRUE,IDs=TRUE,preProc=FALSE)$folds
      
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
      trainy <- (TrainData[,tv])
      testy <- (TestData[,tv])
      
      fit = hierNet.path(trainx,trainy, zz = trainzz,diagonal=FALSE,strong=TRUE,trace=0)
      fitcv = hierNet.cv(fit, trainx, trainy, folds=folds.in.list, trace=0)
      fit.def <- hierNet(trainx,trainy, zz = trainzz,diagonal=FALSE,strong=TRUE,lam=fit$lamlist[which(fitcv$lamhat==fit$lamlist)])
      fit.pred <- predict(fit.def,newx=testx,newzz = testzz)
      fit.pred <- pmax(fit.pred,tv.min/3)
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
    out<-list(measure=results,coef=coef.mat,obs.pred.df = pred)
    return(out)
  }
  
}
