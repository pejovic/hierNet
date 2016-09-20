
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

base.model=SOM.fun; profiles=bor.profs; cov.grids=gridmaps.sm2D; poly.deg=1; num.folds=10; num.means=5; interactions=TRUE; hier=FALSE; preproc=TRUE; seed=321

pre.sparsereg3D <- function(base.model, profiles, cov.grids, hier=FALSE, lambda=seq(0,5,0.1), poly.deg=3, num.folds=5,num.means=3,interactions=TRUE,preproc=TRUE,seed=321){
  
  "%ni%" <- Negate("%in%")
  
  if(hier == TRUE){interactions = TRUE}
  if(interactions == FALSE){hier = FALSE}

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
  
  dummy.par <- dummyVars(as.formula(paste("~", paste(c(all.vars(base.model))[-1], collapse="+"))),profiles,levelsOnly=FALSE)
  profiles <- cbind(profiles[,which(colnames(profiles) %in% c("ID","hdepth",target,coord.names))], predict(dummy.par, newdata = profiles))
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
    if(interactions == TRUE) {
      st.par <- as.data.frame(profiles) %>% subset(., select = c(main.effects,columns.to.keep)) %>% preProcess(.,method=c("center", "scale"))
      profiles <- predict(st.par,newdata = profiles)
    }else{
      st.par <- as.data.frame(profiles) %>% subset(., select = c(main.effects)) %>% preProcess(.,method=c("center", "scale"))
      profiles <- predict(st.par,newdata = profiles)
    }
  }


#======================== Stratification ====================================================
  
  profiles <- as.data.frame(profiles)
  profiles <- plyr::rename(profiles, replace=c("x" = "longitude", "y" = "latitude"))

  strat<-stratfold3d(targetVar=target,seed=seed,regdat=profiles,folds=num.folds,cent=num.means,preProc=FALSE,dimensions="2D",IDs=TRUE,sum=TRUE)
  profile.folds.list <- strat$folds
  obs.folds.list <- stratfold3d(targetVar=target, seed=seed,regdat=profiles, folds=num.folds,cent=num.means, preProc=FALSE, dimensions="2D", IDs=FALSE,sum=TRUE)$folds
#=====================================================================================================
  if(interactions == TRUE){
    out <- list(profiles = profiles, cov.grids = cov.grids, model = list(base.model = base.model, target = target, main.effects = main.effects, poly.deg = poly.deg, interactions = interactions ,hier = hier, depth.interactions = columns.to.keep, all.interactions = colnames(X)), num.folds = num.folds, num.means = num.means, preproc = list(dummy.par = dummy.par, st.par = st.par, columns.to.keep = columns.to.keep), folds = list(profile.folds.list = profile.folds.list, obs.folds.list = obs.folds.list, depth.summary = strat[[3]], target.summmary = strat[[4]]))
  } else {
    out <- list(profiles = profiles, cov.grids = cov.grids, model = list(base.model = base.model, target = target, main.effects = main.effects, poly.deg = poly.deg, interactions = interactions ,hier = hier), num.folds = num.folds, num.means = num.means, preproc = list(dummy.par = dummy.par, st.par = st.par), folds = list(profile.folds.list = profile.folds.list, obs.folds.list = obs.folds.list, depth.summary = strat[[3]], target.summmary = strat[[4]]))
    
  }
  return(out)
  
}  
  

sparsereg3D.ncv <- function(sparsereg = pre.som, lambda = seq(0,5,0.1), seed = 321){
  
  flist <- sparsereg$folds$profile.folds.list
  cv.folds <- sparsereg$folds$obs.folds.list
  profiles <- sparsereg$profiles
  num.folds <- sparsereg$num.folds
  num.means <- sparsereg$num.means
  
  target <- sparsereg$model$target
  target.min <- min(profiles[,target])
  interactions = sparsereg$model$interactions
  hier = sparsereg$model$hier
  main.effects = sparsereg$model$main.effects
  poly.deg = sparsereg$model$poly.deg
  
  if(interactions){
    columns.to.keep = sparsereg$model$depth.interactions
    all.interactions = sparsereg$model$all.interactions
  }

  
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
      if(interactions){
        train.data <- train.data %>% subset(., select = c(target, main.effects, columns.to.keep)) 
        test.data <- test.data %>% subset(., select = c(target, main.effects, columns.to.keep))
      } else {
        train.data <- train.data %>% subset(., select = c(target, main.effects)) 
        test.data <- test.data %>% subset(., select = c(target, main.effects))
      }

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
      
      train.main.effects <- as.matrix(train.data[,main.effects]) 
      test.main.effects <- as.matrix(test.data[,main.effects])
      train.int.effects <- as.matrix(train.data[,all.interactions])
      test.int.effects <- as.matrix(test.data[,all.interactions])
      train.target <- train.data[,target]
      test.target <- test.data[,target]
      
      hier.path = hierNet.path(train.main.effects, train.target, zz = train.int.effects, diagonal=FALSE, strong=TRUE, trace=0, stand.main = FALSE, stand.int = FALSE)
      hier.lasso.cv = hierNet.cv(hier.path, train.main.effects, train.target, folds = new.folds, trace=0)
      final.hier.model <- hierNet(train.main.effects, train.target, zz = train.int.effects, diagonal=FALSE, strong=TRUE, lam=hier.path$lamlist[which(hier.lasso.cv$lamhat==hier.path$lamlist)])
      test.pred <- predict(final.hier.model,newx=test.main.effects,newzz = test.int.effects)
      test.pred <- pmax(test.pred,target.min/3)
      ie <- as.matrix(hier.path$th[,,which(hier.lasso.cv$lamhat==hier.path$lamlist)][,length(main.effects)])
      me<-hier.path$bp[,which(hier.lasso.cv$lamhat==hier.path$lamlist), drop = F] - hier.path$bn[,which(hier.lasso.cv$lamhat==hier.path$lamlist), drop = F]

      train.pred <- predict(final.hier.model,newx=train.main.effects,newzz = train.int.effects)
      train.pred <- pmax(train.pred,target.min/3)
      train.pred <- data.frame(predicted=as.numeric(train.pred))
      train.per.fold <- data.frame(ID = train.data[,"ID"], observed = train.data[,1], predicted = train.pred, longitude = train.data[,"longitude"], latitude = train.data[,"latitude"], depth = profiles[train.obs.ind,"depth"])
      train.results[[i]] <- train.per.fold
      
      test.per.fold <- data.frame(ID = test.data[,"ID"], observed = test.data[,1], predicted = as.numeric(test.pred), longitude = test.data[,"longitude"], latitude = test.data[,"latitude"], depth = profiles[-train.obs.ind,"depth"])
      test.results[[i]] <- test.per.fold
      
      test.obs.pred <- data.frame(obs=test.target,pred=test.pred)
      coef.list[[i]] <- rbind(ie,me)
      fold.measures <- data.frame(lambda=hier.lasso.cv$lamhat,RMSE=defaultSummary(test.obs.pred)[1],Rsquared=defaultSummary(test.obs.pred)[2])
      results[i,] <- fold.measures
      test.prediction<-rbind(test.prediction,test.obs.pred)
      }
    }
  
  results[length(flist)+1,] <- c(NA, RMSE = defaultSummary(test.prediction)[1], Rsquared = defaultSummary(test.prediction)[2])
  coef.mat <- do.call(cbind, coef.list)
  out <- list(accuracy.measures = results, coefficients = coef.mat, test.prediction = test.prediction, folds = cv.folds, train.results = train.results, test.results = test.results)
  return(out)
  
}


sparsereg3D.pred <- function(sparsereg = out, lambda = seq(0,5,0.1), prediction = TRUE, seed = 321, depths = c(-.1,-.3), chunk = 20000){
  
  flist = sparsereg$folds$profile.folds.list
  cv.folds = sparsereg$folds$obs.folds.list
  preproc.par = sparsereg$preproc
  profiles = sparsereg$profiles
  num.folds = sparsereg$num.folds
  num.means = sparsereg$num.means
  cov.grids = cov.grids
  
  target = sparsereg$model$target
  target.min = min(profiles[,target])
  interactions = sparsereg$model$interactions
  base.model = sparsereg$model$base.model
  hier = sparsereg$model$hier
  main.effects = sparsereg$model$main.effects
  poly.deg = sparsereg$model$poly.deg
  main.effects = sparsereg$model$main.effects
  
  if(interactions){
    columns.to.keep = sparsereg$model$depth.interactions
    all.interactions = sparsereg$model$all.interactions
  }
  
  which.folds.list <- as.list(rep(NA,length(flist)))
  names(which.folds.list) <- paste("fold",c(1:length(flist)),sep = "")
  
  which.fold <- rep(NA,dim(profiles)[1])
  
  for(j in 1:length(flist)){
    which.folds.list[[j]] <- which(profiles$ID %in% flist[[j]])
    which.fold[which.folds.list[[j]]]<-j
  }
 
  if(!hier){
    
    if(interactions){
      input.data <- profiles %>% subset(., select = c(target, main.effects, columns.to.keep)) 
    }else{
      input.data <- profiles %>% subset(., select = c(target, main.effects)) 
    }
    
    lasso.cv <- cv.glmnet(as.matrix(input.data[,-1]), input.data[,1], alpha = 1,lambda = lambda, foldid = which.fold, type.measure = "mse")
    lasso.pred <- predict(lasso.cv,s=lasso.cv$lambda.min,newx=as.matrix(input.data[,-1]))
    lasso.pred <- pmax(lasso.pred,target.min/3)
    obs.pred <- data.frame(obs=input.data[,1],pred=as.numeric(lasso.pred))
    coef.list <- predict(lasso.cv,type="coefficients",s=lasso.cv$lambda.min)
    results <- data.frame(lambda=lasso.cv$lambda.min,RMSE=defaultSummary(obs.pred)[1],Rsquared=defaultSummary(obs.pred)[2])
    
    
  }else{
    
    input.data.x <- as.matrix(profiles[,main.effects]) 
    input.data.z <- as.matrix(profiles[,all.interactions])
    input.data.y <- (profiles[,target])
    
    hier.path = hierNet.path(input.data.x,input.data.y, zz = input.data.z, diagonal=FALSE, strong=TRUE, trace=0, stand.main = FALSE, stand.int = FALSE)
    hier.lasso.cv = hierNet.cv(hier.path, input.data.x, input.data.y, folds = which.folds.list, trace=0)
    hier.lasso.final <- hierNet(input.data.x,input.data.y, zz = input.data.z, diagonal=FALSE, strong=TRUE, lam = hier.path$lamlist[which(hier.lasso.cv$lamhat == hier.path$lamlist)], center = TRUE, stand.main = FALSE, stand.int = FALSE)
    hier.lasso.pred <- predict(hier.lasso.final, newx = input.data.x, newzz = input.data.z)
    hier.lasso.pred <- pmax(hier.lasso.pred,target.min/3)
    
    if(poly.deg == 1){ ie <- as.matrix(hier.path$th[,,which(hier.lasso.cv$lamhat==hier.path$lamlist)][,length(columns.to.keep)]) 
              } else {
                       ie <- as.matrix(hier.path$th[,,which(hier.lasso.cv$lamhat==hier.path$lamlist)][,length(columns.to.keep)-poly.deg+1:length(columns.to.keep)]) 
                      }

    me <- hier.path$bp[,which(hier.lasso.cv$lamhat==hier.path$lamlist), drop = F] - hier.path$bn[,which(hier.lasso.cv$lamhat==hier.path$lamlist), drop = F]
    
    
    obs.pred <- data.frame(obs=input.data.y, pred=hier.lasso.pred)
    coef.list <- data.frame(cov.name=colnames(input.data.x),me,ie)
    dfresults <- data.frame(lambda = hier.lasso.cv$lamhat, RMSE = defaultSummary(obs.pred)[1], Rsquared = defaultSummary(obs.pred)[2])
   
  }
 
  if(!hier){
    regression.summary <- list(model = list(model = lasso.cv, lambda = lasso.cv$lambda.min), preproc.par = preproc.par, coefficients = coef.list, obs.pred = data.frame(profiles[,1:5], obs.pred))
    }else{
    regression.summary <- list(model = list(model = hier.lasso.final, lambda = hier.path$lamlist[which(hier.lasso.cv$lamhat == hier.path$lamlist)]), preproc.par = preproc.par, coefficients = coef.list, obs.pred = data.frame(profiles[,1:5], obs.pred))
    }

  if(prediction){
    
    grids.3D <- sp3D(cov.grids, stdepths = depths)
    cores = 1
    grids.3D <- mclapply(grids.3D, function(x) as.data.frame(x), mc.cores = cores)
    grids.3D <- lapply(grids.3D, function(x) {names(x) <- c(names(x)[1:(length(names(x))-(poly.deg))], c("depth"));return(x)})
    
    grids.3D <- mclapply(grids.3D, function(x) subset(x,select=c(all.vars(base.model))[-1],drop=FALSE),mc.cores=cores) %>% mclapply(.,function(x) predict(preproc.par$dummy.par,newdata=x),mc.cores=cores) %>% mclapply(., function(x) {colnames(x) <- gsub( "\\_|/|\\-|\"|\\s" , "." , colnames(x) );return(x)},mc.cores=cores)

    if(poly.deg > 1) { grids.3D <- lapply(grids.3D, function(x) x <- cbind(x, poly(x$depth, poly.deg, raw=TRUE, simple=TRUE)[,-1]))
    grids.3D <- lapply(grids.3D, function(x) {names(x) <- c(names(x)[1:(length(names(x))-(poly.deg))], c("depth",paste("depth",c(2:poly.deg),sep="")));return(x)})}
    
    #==================== Compute Interactions ==============================================================================
    
    if(interactions){
      chunk <- chunk # chunk = 20000
      n <- nrow(grids.3D[[1]])
      r  <- rep(1:ceiling(n/chunk),each = chunk)[1:n]
      grids.int.3D <- lapply(grids.3D, function(x) as.data.frame(x)) %>% lapply(., function(x) split(x,r))
      
      m.cores <- detectCores()
      registerDoParallel(cores=m.cores)
      
      for(i in 1:length(grids.int.3D)){ 
        grids.int.3D[[i]] <- foreach(j = grids.int.3D[[i]],.combine="rbind") %dopar% {hierNet::compute.interactions.c(as.matrix(j),diagonal = FALSE)}
      }
    }
    
    #======================= Standardization ============================
    for( i in 1:length(grids.int.3D)) {
      grids.3D[[i]] <- cbind(grids.3D[[i]],grids.int.3D[[i]])
    }
    
    grids.3D <- mclapply(grids.3D,function(x) predict(preproc.par$st.par,newdata=x),mc.cores=cores) %>% mclapply(.,function(x) subset(x[,which(colnames(x) %in% c(main.effects,columns.to.keep))]))
    
    #======================= Columns to keep ============================
    
    if(interactions == TRUE){
      if(hier == TRUE){
        
        for(i in 1:length(grids.int.3D)) {
          grids.3D[[i]] <- subset(grids.3D[[i]][,which(colnames(grids.3D[[i]]) %in% main.effects)])
          grids.int.3D[[i]][,colnames(grids.int.3D[[i]]) %ni% columns.to.keep ] <- 0
        } 
      }  
    }
    #====================== Prediction ==================================
    
    if(!hier){

      for(i in 1:length(grids.3D)){
        grids.3D[[i]]$pred <- as.numeric(lapply(grids.3D, function(x) predict(lasso.cv,s=lasso.cv$lambda.min,newx = as.matrix(x)))[[i]])
        grids.3D[[i]]$pred <- pmax(grids.3D[[i]]$pred,target.min/3)
        grids.3D[[i]] <- grids.3D[[i]][,"pred"]
      }
      
    } else {
      
      for( i in 1:length(grids.3D)){
        grids.3D[[i]]$pred <- as.numeric(predict(hier.lasso.final, newx=grids.3D[[i]], newzz = grids.int.3D[[i]]))
        grids.3D[[i]]$pred <- pmax(grids.3D[[i]]$pred,min.obs/3)
        grids.3D[[i]] <- grids.3D[[i]][,"pred"]
      }
      
    }
    return(list(prediction=grids.3D, regression.summary = regression.summary))
  } 
  
    else {return(regression.summary = regression.summary)}
  
}
  

pre.som <- pre.sparsereg3D(base.model = SOM.fun, hier = TRUE, profiles=bor.profs, cov.grids = gridmaps.sm2D)

sparsereg = out


pre.som <- pre.sparsereg3D(base.model = SOM.fun, hier = FALSE, profiles=bor.profs, interactions = TRUE, num.folds = 10, cov.grids = gridmaps.sm2D, poly.deg = 3)    
sp.reg.som <- sparsereg3D.ncv(sparsereg = pre.som, lambda = seq(0,5,0.1), seed = 321)

sp.reg.som <- sparsereg3D.pred(sparsereg = pre.som, prediction =  TRUE ,lambda = seq(0,5,0.1), seed = 321)

sp.reg.som$accuracy.measures  

  
  
  
