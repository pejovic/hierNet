# Create stratified folds, taking into account 3D position of observation.
#   1. In the first step, observations are clustered according to k-means clustering (k=cent)
#   2. then, in each cluster:
#       a. Each profiles are averaged by weigthted mean
#       b. profiles were sampled randomly according to number of folds
#   3. Merge each fold with corresponding folds in other clusters


# regdat - input data matrix with "ID", "hdept" and coordinates columns.
# targetVar - name of target variable
# folds - number of folds
# cent - number of centers for k-means clustering
# dimensions - "2D" or "3D"; if dimensions = "3D", k-means clustering takes depth into computation
# IDs (binary) - identify if the list of folds with profiles' IDs must be returned
# sum (binary) - identify if the summary statistics for each fold should be computed (to inspect the folds)


stratfold3d <- function(targetVar, regdat, folds = 5, cent = 3, seed = 666, dimensions = list("2D","3D"), IDs = TRUE, sum = FALSE){
  
  dimensions <- dimensions[[1]]
  if(dimensions == "2D"){
    
    unique.df <- ddply(regdat, .(ID), here(summarize), target=weighted.mean(eval(parse(text=targetVar)), hdepth), longitude=longitude[1], latitude=latitude[1])
    km <- kmeans(cbind(unique.df$longitude, unique.df$latitude),  centers = cent)
    unique.df$km <- as.factor(km$cluster)} else {
      
      unique.df <- ddply(regdat, .(ID), here(summarize), target=weighted.mean(eval(parse(text=targetVar)), hdepth), longitude=longitude[1], latitude=latitude[1], altitude=min(altitude))
      km <- kmeans(cbind(scale(unique.df$longitude), scale(unique.df$latitude), scale(unique.df$altitude)),  centers = cent)
      unique.df$km <- as.factor(km$cluster)
    }
  
  #============== Creating a list of cluster regions =======================================
  
  k.list <- as.list(rep(NA,length(unique(unique.df$km))))
  names(k.list) <- paste("k",c(1:length(k.list)),sep="")
  
  #============== Creating list of n folds per cluster region with profile indexes =========
  
  for(i in 1:length(k.list)){
    set.seed(seed)
    k.list[[i]] <- createFolds(unique.df[which(unique.df$km==levels(unique.df$km)[i]),"target"],k=folds)
    for(j in 1:folds){
      k.list[[i]][[j]] <- unique.df[which(unique.df$km == levels(unique.df$km)[i]),"ID"][k.list[[i]][[j]]]
    }
  }
  #=========================================================================================
  
  #============== Creating list of folds with indexes of profiles per each fold ============
  ID.list <- as.list(rep(NA,folds))
  names(ID.list) <- paste("fold",c(1:folds),sep = "")
  
  for(i in 1:folds){
    ID.list[[i]] <- do.call(c,lapply(k.list,function(x) x[[i]]))
    names(ID.list[[i]]) <- NULL
  }
  
  #=========================================================================================
  
  #================= Creating list of folds with observations' indexes =====================
  
  folds.list <- as.list(rep(NA,folds))
  names(folds.list) <- paste("fold",c(1:folds),sep = "")
  for(i in 1:length(ID.list)){
    folds.list[[i]] <- which(regdat$ID %in% ID.list[[i]])
  }
  #=========================================================================================
  
  pom <- data.frame()
  for(i in 1:length(folds.list)){
    allData1 <- regdat[folds.list[[i]],]
    allData1$fold <- paste("fold",i,sep="")
    allData <- rbind(allData1,pom)
    pom <- allData
  }
  
  allData$fold<-factor(allData$fold)
  sum.list=list(by(allData$depth,allData$fold,summary))
  if(IDs==TRUE){index.list=(ID.list)}else{index.list=(folds.list)}
  sum.list<-list(allData,index.list,sum.list,by(allData[,paste(targetVar)],allData$fold,summary))
  names(sum.list)<-c("Data","folds","depth summary",paste(targetVar,"summary", sep=" "))
  
  if(sum==TRUE){return(sum.list)}else{return(sum.list[[1]])}
}




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
  target <- all.vars(base.model)[1] # extract name of target variable
  ## check if all covariates are available:
  if(sum(names(cov.grids) %in% all.vars(base.model)[-1])==0){
    stop("None of the covariates in the 'formulaString' matches the names in the 'covariates' object")
  }
  #===============================================================
  p4s = proj4string(profiles) # extract spatial reference from profiles
  profiles <- join(profiles@horizons,data.frame(data.frame(profiles@site),data.frame(profiles@sp)),type="inner") #creating "profiles" data frame with all variables: "ID","TOp","Bottom","observed variables","x","y".
  
  #============== Names ==========================================
  coord.names <- tail(names(profiles),2); #extract names of coordiantes axes
  pr.names <- head(names(profiles),3); # extract names of profiles specific variables: "ID","Top","Bottom".
  
  #============== Adding altitude and hdepth =====================
  profiles$depth <- - (profiles$Top / 100 + (( profiles$Bottom - profiles$Top ) / 2) / 100) #calculating depth as depth of horizon mid point in "cm"
  profiles$hdepth <- profiles$Bottom - profiles$Top # Calculating depth of each horizon in each profile (necessary for data stratification)
  
  #============== Overlay ========================================
  profiles <- profiles[complete.cases(profiles[,c("ID",coord.names,"hdepth","depth",target)]),c("ID",target,"hdepth",coord.names,"depth")] #Removing NA values according to: "ID",coord.names,"hdepth","depth",target
  coordinates(profiles) <- ~ x + y #creating spatial class (sp) from profiles by defining the spatial coordinates in order to make spatial overlay possible
  proj4string(profiles) <- p4s # defining spatial reference
  profiles <- spTransform(profiles, proj4string(cov.grids)) # Spatial transform of cov.grids in order to provide the same spatial reference
  ov <- over(profiles, cov.grids) # spatial overlay with cov.grids in order to extract the covariate values at sample locations
  
  factor.names <- ov %>% subset(., select=which(sapply(., is.factor))) %>% names() # Extracting the names of categorical variagles
  
  for(i in factor.names){
    ov[,i] <- factor(ov[,i]) # removing empty classes from categorical variables
  }
  
  #======== Preparing data input matrix with following columns: "ID",target,"hdepth",coord.names, sp.cov.names, "depth" 
  
  sp.cov.names <- names(ov[,which(names(ov) %in% c(all.vars(base.model)))]) # extracting the names of all spatial covarites ( which participate in model formula )
  profiles <- cbind(as.data.frame(profiles), ov[,sp.cov.names]) # combining profiles with covariates data
  profiles <- profiles[complete.cases(profiles[,all.vars(base.model)]),c("ID",target,"hdepth",coord.names, sp.cov.names, "depth")] #removing NA values from input data matrix
  
  #======== Adding polynomial depth terms in input data matrix, only if poly.deg > 1  
  
  if(poly.deg > 1){
    profiles <- cbind(profiles,poly(profiles$depth,poly.deg,raw=TRUE,simple=TRUE)[,-1])
    names(profiles) <- c(names(profiles)[1:(length(names(profiles))-(poly.deg-1))],(paste("depth",c(2:poly.deg),sep="")))
  }

  #===================== Dummy coding =======================================
  
  dummy.par <- dummyVars(as.formula(paste("~", paste(c(all.vars(base.model))[-1], collapse="+"))),profiles,levelsOnly=FALSE) # extracting and storing dummy coding scheme for transformation of categorical variables
  profiles <- cbind(profiles[,which(colnames(profiles) %in% c("ID","hdepth",target,coord.names))], predict(dummy.par, newdata = profiles)) # applying dummy coding scheme on input data matrix (profiles)
  colnames(profiles) <- gsub( "\\_|/|\\-|\"|\\s" , "." , colnames(profiles) ) # Replacing special characters from variable names by "."
  main.effects <- colnames(profiles)[-which(colnames(profiles) %in% c("ID","hdepth",target,coord.names))] # Extracting names of all variables, including the new variables created from dummy classes.
  #==========================================================================

  #======= Computing interactions in input data matrix, only if interations == TRUE. =============================================    
  if (interactions == TRUE){
    
    X <- hierNet::compute.interactions.c(as.matrix(profiles[,-c(which(colnames(profiles) %in% c("ID",target,"hdepth",coord.names)))]),diagonal=FALSE) # Computing all two-pairs interactions between main effects; ("ID",target,"hdepth",coord.names) are excluded from computation
    
    #======= Defining the names of columns of final interaction model, i.e. names of main effects as well as intreactions effects (interactions between spatial covariates and depth terms)  
    if(poly.deg > 1){ columns.to.keep <- colnames(X[,do.call(c,lapply(strsplit(colnames(X),":"), function(x) x[2] %in% c("depth",paste("depth",c(2:poly.deg),sep="")) & x[1] %ni% c("depth",paste("depth",c(2:poly.deg),sep=""))))]) } else {
                      columns.to.keep <- (X %>% as.data.frame() %>% subset(., select=grep("depth", names(.), value=TRUE)) %>% colnames())
                      }
    
    #======= IF hier=TRUE, interactions columns other than interctions with depth become zero columns, while in non-interactions setting these columns are excluded.
    if(hier == TRUE) { X[,colnames(X) %ni% columns.to.keep ] <- 0 
                       profiles <- cbind(profiles,X)
                     }else{
                            profiles <- cbind(profiles,X[,colnames(X) %in% columns.to.keep])
                          }
    
  }

  
  
#===================== Preprocessing (standardization) of input data  =======================================

  if(preproc==TRUE){
    if(interactions == TRUE) {
                  # input data matrix   %>%    selection of columns to be transformed           %>%     computin preprocessing parameters
      st.par <- as.data.frame(profiles) %>% subset(., select = c(main.effects,columns.to.keep)) %>% preProcess(.,method=c("center", "scale"))
      profiles <- predict(st.par,newdata = profiles) # Applying preprocessing parameters on input data matrix
    }else{
      # if interactions = FALSE, only the main effects are transformed
      st.par <- as.data.frame(profiles) %>% subset(., select = c(main.effects)) %>% preProcess(.,method=c("center", "scale"))
      profiles <- predict(st.par,newdata = profiles)
    }
  }


#======================== Data stratification ====================================================
  
  profiles <- as.data.frame(profiles) #
  profiles <- plyr::rename(profiles, replace=c("x" = "longitude", "y" = "latitude")) #renaming coordinate axes to "longitude" "latitude" (requirement of stratfold3D function)

  strat <- stratfold3d(targetVar = target, seed = seed, regdat = profiles, folds = num.folds, cent = num.means, preProc = FALSE, dimensions = "2D", IDs = TRUE, sum = TRUE) # Applying "stratfold3D" function in order to make stratified folds
  profile.folds.list <- strat$folds # Extracting list of folds with profiles indexes (with indexes identifying what fold each profile contains)
  obs.folds.list <- stratfold3d(targetVar=target, seed=seed,regdat=profiles, folds=num.folds,cent=num.means, preProc=FALSE, dimensions="2D", IDs=FALSE,sum=TRUE)$folds ## Extracting list of folds of observation indexes (with indexes identifying what observation each profile contains)
#=====================================================================================================
  
  # Creating of object for output containing all necessary data for model selection, prediction and model assessment via nested cv.
  if(interactions == TRUE){
    out <- list(profiles = profiles, cov.grids = cov.grids, model = list(base.model = base.model, target = target, main.effects = main.effects, poly.deg = poly.deg, interactions = interactions ,hier = hier, depth.interactions = columns.to.keep, all.interactions = colnames(X)), num.folds = num.folds, num.means = num.means, preproc = list(dummy.par = dummy.par, st.par = st.par, columns.to.keep = columns.to.keep), folds = list(profile.folds.list = profile.folds.list, obs.folds.list = obs.folds.list, depth.summary = strat[[3]], target.summmary = strat[[4]]))
  } else {
    out <- list(profiles = profiles, cov.grids = cov.grids, model = list(base.model = base.model, target = target, main.effects = main.effects, poly.deg = poly.deg, interactions = interactions ,hier = hier), num.folds = num.folds, num.means = num.means, preproc = list(dummy.par = dummy.par, st.par = st.par), folds = list(profile.folds.list = profile.folds.list, obs.folds.list = obs.folds.list, depth.summary = strat[[3]], target.summmary = strat[[4]]))
    
  }
  return(out)
  
}  
  
#=========== sparsereg3D.ncv - function for model assessment via nested cv. =======
#==== Input data are :
# sparsereg - output from "pre.sparsereg3D" function
# lambda - vector of lambda values for lasso regression
# seed - random number generator (in order to make computation reproducible)

sparsereg3D.ncv <- function(pre.sparsereg, lambda = seq(0,5,0.1), seed = 321){
  
  #== Extracting data from "sparsereg" object:
  flist <- sparsereg$folds$profile.folds.list # list of folds with indexes of profiles
  cv.folds <- sparsereg$folds$obs.folds.list # list of folds with indexes of observations
  profiles <- sparsereg$profiles # input data
  num.folds <- sparsereg$num.folds # number of folds (necessary for data stratification in inner loop in ncv)
  num.means <- sparsereg$num.means # number of spatial cluster means (necessary for data stratification in inner loop in ncv)
  
  target <- sparsereg$model$target # name of target variable
  target.min <- min(profiles[,target]) # minimum of observed value of target variable (necessary for correcting negative predictions)
  interactions = sparsereg$model$interactions # binary parameter (identifying does interactions = TRUE or FALSE)
  hier = sparsereg$model$hier # binary parameter (identifying does hier = TRUE or FALSE)
  main.effects = sparsereg$model$main.effects # names of main effects
  poly.deg = sparsereg$model$poly.deg # degree of polynomial depth function
  
  if(interactions){
    columns.to.keep = sparsereg$model$depth.interactions # names of columns to keep (which variables constitute final model)
    all.interactions = sparsereg$model$all.interactions # names of all interactions terms
  }

  # Preparing empty data frames and lists which will contain the results of procedure.
  results <- data.frame(lambda = rep(NA,length(flist)+1), RMSE = rep(NA,length(flist)+1), Rsquared = rep(NA,length(flist)+1))
  coef.list = as.list(rep(NA,length(flist)))
  test.prediction <- data.frame()
  train.results <- as.list(rep(NA,length(flist)))
  test.results <- as.list(rep(NA,length(flist)))
  #==================================================================================
  
  #===== Inner data partitioning ====================================================
  
  for(i in 1:length(flist)){
    train.obs.ind <- which(profiles$ID %in% do.call(c, flist[-i]))
    train.data <- profiles[train.obs.ind,]
    train.profiles.ids <- flist[-i]
    test.data <- profiles[-train.obs.ind,]
    test.profiles.ids <- flist[i]
    
    # =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  
    
    inner.partitioning <- stratfold3d(target, train.data, folds = num.folds, seed = seed ,cent = num.means, dimensions = "2D", sum = TRUE, IDs = TRUE, preProc = FALSE)$folds
    
    new.folds <- as.list(rep(NA,length(inner.partitioning)))
    names(new.folds) <- paste("fold",c(1:length(inner.partitioning)),sep = "")
    
    which.new.fold <- rep(NA,dim(train.data)[1]) # vector of indexes identifying what inner fold each observation is in.
    
    for(j in 1:length(inner.partitioning)){
      new.folds[[j]]<-which(train.data$ID %in% inner.partitioning[[j]])
      which.new.fold[new.folds[[j]]]<-j
    }
    #======================================================================================
    
    #============== LASSO training ========================================================
    if(!hier){
      if(interactions){
        train.data <- train.data %>% subset(., select = c(target, main.effects, columns.to.keep)) 
        test.data <- test.data %>% subset(., select = c(target, main.effects, columns.to.keep))
      } else {
        train.data <- train.data %>% subset(., select = c(target, main.effects)) 
        test.data <- test.data %>% subset(., select = c(target, main.effects))
      }

      #   Lasso  
      
      lasso.cv <- cv.glmnet(as.matrix(train.data[,-1]), train.data[,1], alpha = 1,lambda = lambda, foldid = which.new.fold, type.measure = "mse") # lasso traning in inner loop
      test.pred <- predict(lasso.cv, s = lasso.cv$lambda.min, newx = as.matrix(test.data[,-1])) # prediction on test data (outer fold)
      test.pred <- pmax(test.pred, target.min/3) # correcting zero predictions to the 1/3 of min observation
      
      # Training prediction? ==============================================================================================================
      train.pred <- predict(lasso.cv,s = lasso.cv$lambda.min, newx = as.matrix(train.data[,-1])) # prediction on training data
      train.pred <- pmax(train.pred,target.min/3)  # correcting zero predictions to the 1/3 of min observation
      train.pred <- data.frame(predicted = as.numeric(train.pred)) # creating data frame from training prediction
      train.per.fold <- data.frame(ID = profiles[train.obs.ind,"ID"], observed = train.data[,1], predicted = train.pred,longitude = profiles[train.obs.ind,"longitude"], latitude = profiles[train.obs.ind,"latitude"], depth = profiles[train.obs.ind,"depth"]) # creating data frame of training predictions with additional columns that pertain to: longitude, latitude, depth etc. (In order to analyze training residuals for further vaiogram modeling) 
      train.results[[i]] <- train.per.fold # Storing training prediction
      #==========================================================================================================================================================
      
      # Test prediction =========================================================================================================================================
      ## creating data frame of test predictions with additional columns that pertain to: longitude, latitude, depth etc. (In order to analyze test residuals for further vaiogram modeling) 
      test.per.fold <- data.frame(ID = profiles[-train.obs.ind,"ID"], observed = test.data[,1], predicted = as.numeric(test.pred), longitude = profiles[-train.obs.ind,"longitude"], latitude = profiles[-train.obs.ind,"latitude"], depth = profiles[-train.obs.ind,"depth"])
      test.results[[i]] <- test.per.fold
      
      test.obs.pred <- data.frame(obs = test.data[,1], pred = as.numeric(test.pred)) # test obs/pred data frame for each ncv model
      coef.list[[i]] <- predict(lasso.cv, type = "coefficients", s = lasso.cv$lambda.min) # list of coefficients for each model in nested cv
      fold.measures <- data.frame(lambda = lasso.cv$lambda.min, RMSE = defaultSummary(test.obs.pred)[1], Rsquared = defaultSummary(test.obs.pred)[2]) # Computation of accuracy measures for each step in ncv
      results[i,] <- fold.measures # Storing computed accuracy measures for each step in ncv
      test.prediction <- rbind(test.prediction,test.obs.pred) # combining all test prediction in one data .frame
      #==========================================================================================================================================================
      
    }else{
      # Hierarchy setting requires the separation of main effects and interaction effects and also the training and test data.
      train.main.effects <- as.matrix(train.data[,main.effects]) 
      test.main.effects <- as.matrix(test.data[,main.effects])
      train.int.effects <- as.matrix(train.data[,all.interactions])
      test.int.effects <- as.matrix(test.data[,all.interactions])
      train.target <- train.data[,target]
      test.target <- test.data[,target]
      #==================================================================================================================================
      
      hier.path = hierNet.path(train.main.effects, train.target, zz = train.int.effects, diagonal=FALSE, strong=TRUE, trace=0, stand.main = FALSE, stand.int = FALSE) # Fit path of lambda on training data
      hier.lasso.cv = hierNet.cv(hier.path, train.main.effects, train.target, folds = new.folds, trace=0) # perform cross-validation on whole path
      final.hier.model <- hierNet(train.main.effects, train.target, zz = train.int.effects, diagonal=FALSE, strong=TRUE, lam=hier.path$lamlist[which(hier.lasso.cv$lamhat==hier.path$lamlist)]) # Fit the final model on training data
      test.pred <- predict(final.hier.model,newx=test.main.effects,newzz = test.int.effects) # Perform the prediction on test data
      test.pred <- pmax(test.pred,target.min/3) # Correcting negative predictions
      # Extracting the coefficients of each model in ncv; (ie-interaction effects; me- main effects)
      ie <- as.matrix(hier.path$th[,,which(hier.lasso.cv$lamhat==hier.path$lamlist)][,length(main.effects)]) #Extracting the intercation coefficients of final model (in inner loop)
      me <- hier.path$bp[,which(hier.lasso.cv$lamhat==hier.path$lamlist), drop = F] - hier.path$bn[,which(hier.lasso.cv$lamhat==hier.path$lamlist), drop = F] #Extracting the main coefficients of final model (in inner loop)
      
      # Training prediction? ==============================================================================================================
      train.pred <- predict(final.hier.model,newx=train.main.effects,newzz = train.int.effects)
      train.pred <- pmax(train.pred,target.min/3)
      train.pred <- data.frame(predicted=as.numeric(train.pred))
      train.per.fold <- data.frame(ID = train.data[,"ID"], observed = train.data[,1], predicted = train.pred, longitude = train.data[,"longitude"], latitude = train.data[,"latitude"], depth = profiles[train.obs.ind,"depth"])
      train.results[[i]] <- train.per.fold
      #====================================================================================================================================
    
      # Test prediction ==============================================================================================================
      test.per.fold <- data.frame(ID = test.data[,"ID"], observed = test.data[,1], predicted = as.numeric(test.pred), longitude = test.data[,"longitude"], latitude = test.data[,"latitude"], depth = profiles[-train.obs.ind,"depth"])
      test.results[[i]] <- test.per.fold
      
      test.obs.pred <- data.frame(obs=test.target,pred=test.pred) # test obs/pred data frame for each ncv model
      coef.list[[i]] <- rbind(ie,me) # matrix of coefficinets for each model in ncv
      fold.measures <- data.frame(lambda=hier.lasso.cv$lamhat,RMSE=defaultSummary(test.obs.pred)[1],Rsquared=defaultSummary(test.obs.pred)[2]) # Computation of accuracy measures for each step in ncv
      results[i,] <- fold.measures # Storing computed accuracy measures for each step in ncv
      test.prediction<-rbind(test.prediction,test.obs.pred) # combining all test prediction in one data .frame
      #===============================================================================================================================
      
      }
    }
  
  # Storing results ==============================================================================================================================
  results[length(flist)+1,] <- c(NA, RMSE = defaultSummary(test.prediction)[1], Rsquared = defaultSummary(test.prediction)[2])
  coef.mat <- do.call(cbind, coef.list)
  out <- list(accuracy.measures = results, coefficients = coef.mat, test.prediction = test.prediction, folds = cv.folds, train.results = train.results, test.results = test.results)
  return(out)
  
}

#=========== sparsereg3D.pred - function for sparse model selection and prediction =======
#==== Input data are :
# sparsereg - output from "pre.sparsereg3D" function
# lambda - vector of lambda values for lasso regression
# seed - random number generator (in order to make computation reproducible)
# prediction (binary) - if prediction over grid should be performed
# depths - depths defining on which depth prediction has to be made (only if prediction = TRUE)
# chunk - number that defines the number of cells which will be process in one step (in order to speed-up the computation process)

sparsereg3D.pred <- function(pre.sparsereg, lambda = seq(0,5,0.1), prediction = TRUE, seed = 321, depths = c(-.1,-.3), chunk = 1000){

  #== Extracting data from "sparsereg" object:
  flist = sparsereg$folds$profile.folds.list # list of folds with indexes of profiles
  cv.folds = sparsereg$folds$obs.folds.list # list of folds with indexes of observations
  preproc.par = sparsereg$preproc # preprocessing paramters (necessary for preparing grids for prediction)
  profiles = sparsereg$profiles # input data
  cov.grids = cov.grids # prediction grids with covariates
  
  target = sparsereg$model$target # name of target variable
  target.min = min(profiles[,target]) # minimum of observed value of target variable (necessary for correcting negative predictions)
  interactions = sparsereg$model$interactions # binary parameter (identifying does interactions = TRUE or FALSE)
  base.model = sparsereg$model$base.model # base model
  hier = sparsereg$model$hier # binary parameter (identifying does hier = TRUE or FALSE)
  main.effects = sparsereg$model$main.effects # names of main effects
  poly.deg = sparsereg$model$poly.deg # degree of polynomial depth function
  
  if(interactions){
    columns.to.keep = sparsereg$model$depth.interactions # names of columns to keep (which variables constitute final model)
    all.interactions = sparsereg$model$all.interactions # names of all interactions terms (necessary for "hier" setting)
  }
  
  which.folds.list <- as.list(rep(NA,length(flist)))
  names(which.folds.list) <- paste("fold",c(1:length(flist)),sep = "")
  
  which.fold <- rep(NA,dim(profiles)[1]) # vector of indexes identifying what fold each observation is in.
  
  for(j in 1:length(flist)){
    which.folds.list[[j]] <- which(profiles$ID %in% flist[[j]])
    which.fold[which.folds.list[[j]]]<-j
  }
  
  #============== LASSO training ========================================================
  
  if(!hier){
    
    if(interactions){
      input.data <- profiles %>% subset(., select = c(target, main.effects, columns.to.keep)) 
    }else{
      input.data <- profiles %>% subset(., select = c(target, main.effects)) 
    }
    
    #   Lasso  ===============================================================================================
    
    lasso.cv <- cv.glmnet(as.matrix(input.data[,-1]), input.data[,1], alpha = 1,lambda = lambda, foldid = which.fold, type.measure = "mse") # lasso cv based on data partitioning defined in "which.fold"
    lasso.pred <- predict(lasso.cv,s=lasso.cv$lambda.min,newx=as.matrix(input.data[,-1])) # Prediction of final model on whole data set (maybe not necessary)
    lasso.pred <- pmax(lasso.pred,target.min/3) # correcting of negative predictions
    obs.pred <- data.frame(obs=input.data[,1],pred=as.numeric(lasso.pred)) # obs/pred data frame
    coef.list <- predict(lasso.cv,type="coefficients",s=lasso.cv$lambda.min) # coefficients of final model
    results <- data.frame(lambda=lasso.cv$lambda.min,RMSE=defaultSummary(obs.pred)[1],Rsquared=defaultSummary(obs.pred)[2]) # Accyracy measures derived from prediction of final model on whole data (maybe not necessary)
    
    
  }else{
    
    # Hierarchy setting requires the separation of main effects and interaction effects and also the training and test data.
    input.data.x <- as.matrix(profiles[,main.effects]) 
    input.data.z <- as.matrix(profiles[,all.interactions])
    input.data.y <- (profiles[,target])
    
    hier.path = hierNet.path(input.data.x,input.data.y, zz = input.data.z, diagonal=FALSE, strong=TRUE, trace=0, stand.main = FALSE, stand.int = FALSE) # Fit path of lambda on whole data
    hier.lasso.cv = hierNet.cv(hier.path, input.data.x, input.data.y, folds = which.folds.list, trace=0) # perform cross-validation on whole path
    hier.lasso.final <- hierNet(input.data.x,input.data.y, zz = input.data.z, diagonal=FALSE, strong=TRUE, lam = hier.path$lamlist[which(hier.lasso.cv$lamhat == hier.path$lamlist)], center = TRUE, stand.main = FALSE, stand.int = FALSE) # Fit the final model on whole data
    hier.lasso.pred <- predict(hier.lasso.final, newx = input.data.x, newzz = input.data.z) # Perform the prediction on whole data (not necessary)
    hier.lasso.pred <- pmax(hier.lasso.pred,target.min/3) # correcting of negative predictions
    
    # Extracting the coefficients of final model (ie-interaction effects; me- main effects)
    if(poly.deg == 1){ ie <- as.matrix(hier.path$th[,,which(hier.lasso.cv$lamhat==hier.path$lamlist)][,length(columns.to.keep)]) 
              } else {
                       ie <- as.matrix(hier.path$th[,,which(hier.lasso.cv$lamhat==hier.path$lamlist)][,length(columns.to.keep)-poly.deg+1:length(columns.to.keep)]) 
                      }

    me <- hier.path$bp[,which(hier.lasso.cv$lamhat==hier.path$lamlist), drop = F] - hier.path$bn[,which(hier.lasso.cv$lamhat==hier.path$lamlist), drop = F]
    #============================================================================================================================================================
    
    obs.pred <- data.frame(obs=input.data.y, pred=hier.lasso.pred) # obs/pred data frame (not necessary)
    coef.list <- data.frame(cov.name=colnames(input.data.x),me,ie) # coefficients
    dfresults <- data.frame(lambda = hier.lasso.cv$lamhat, RMSE = defaultSummary(obs.pred)[1], Rsquared = defaultSummary(obs.pred)[2]) # Accyracy measures derived from prediction of final model on whole data (maybe not necessary)
   
  }
 
  # Regression summary list containing the final model, final lambda, coefficients, obs/pred data frames)
  if(!hier){
    regression.summary <- list(model = list(model = lasso.cv, lambda = lasso.cv$lambda.min), preproc.par = preproc.par, coefficients = coef.list, obs.pred = data.frame(profiles[,1:5], obs.pred))
    }else{
    regression.summary <- list(model = list(model = hier.lasso.final, lambda = hier.path$lamlist[which(hier.lasso.cv$lamhat == hier.path$lamlist)]), preproc.par = preproc.par, coefficients = coef.list, obs.pred = data.frame(profiles[,1:5], obs.pred))
    }

  
  #====== Prediction on grids =============================================================================================================================
  if(prediction){
    
    grids.3D <- sp3D(cov.grids, stdepths = depths) # Creating 3D grids (list of grids, each corresponds to different depth)
    cores = 1 # defining number of cores (necessary to be 1, because function bellow doesn't work on many cores)
    grids.3D <- mclapply(grids.3D, function(x) as.data.frame(x), mc.cores = cores) # Converting grids into data frames
    grids.3D <- lapply(grids.3D, function(x) {names(x) <- c(names(x)[1:(length(names(x))-(poly.deg))], c("depth"));return(x)}) # renameing the name of "altitude" column to "depth" column in each grid. (Besause function above uses word "altitude" to denote "depth") )
    
    grids.3D <- mclapply(grids.3D, function(x) subset(x,select=c(all.vars(base.model))[-1],drop=FALSE),mc.cores=cores) %>% mclapply(.,function(x) predict(preproc.par$dummy.par,newdata=x),mc.cores=cores) %>% mclapply(., function(x) {colnames(x) <- gsub( "\\_|/|\\-|\"|\\s" , "." , colnames(x) );return(x)},mc.cores=cores) # Applying dummy.par to transform categorical variables in each grid.
    
    # Adding polynomial depth variables in each grid (only if poly.deg > 1)
    if(poly.deg > 1) { grids.3D <- lapply(grids.3D, function(x) x <- cbind(x, poly(x$depth, poly.deg, raw=TRUE, simple=TRUE)[,-1]))
    grids.3D <- lapply(grids.3D, function(x) {names(x) <- c(names(x)[1:(length(names(x))-(poly.deg))], c("depth",paste("depth",c(2:poly.deg),sep="")));return(x)})}
    
    #==================== Compute Interactions ==============================================================================
    
    if(interactions){
      n <- nrow(grids.3D[[1]]) # defining how many rows each grid contains
      r  <- rep(1:ceiling(n/chunk),each = chunk)[1:n] # split each grid according to chunk
      grids.int.3D <- lapply(grids.3D, function(x) as.data.frame(x)) %>% lapply(., function(x) split(x,r)) # Split each grid in several smaller data frame, according to chunk values
      
      m.cores <- detectCores() # detect the number of cores (in order to split computation of interactions per cores)
      registerDoParallel(cores=m.cores) # defining parallel computation process
      # Computing interactions for each grid, chunk after chunk, and then combine.
      for(i in 1:length(grids.int.3D)){ 
        grids.int.3D[[i]] <- foreach(j = grids.int.3D[[i]],.combine="rbind") %dopar% {hierNet::compute.interactions.c(as.matrix(j),diagonal = FALSE)}
      }
    }
    
    #======================= Standardization of each grid according to the parameters stored previously (st.par) ============================
    for( i in 1:length(grids.int.3D)) {
      grids.3D[[i]] <- cbind(grids.3D[[i]],grids.int.3D[[i]])
    }
    
    grids.3D <- mclapply(grids.3D,function(x) predict(preproc.par$st.par,newdata=x),mc.cores=cores) %>% mclapply(.,function(x) subset(x[,which(colnames(x) %in% c(main.effects,columns.to.keep))]))
    
    #======================= Columns to keep ============================
    # For hierarchical setting interaction effects other than interactions with depth must be setted to zero columns
    if(interactions == TRUE){
      if(hier == TRUE){
        
        for(i in 1:length(grids.int.3D)) {
          grids.3D[[i]] <- subset(grids.3D[[i]][,which(colnames(grids.3D[[i]]) %in% main.effects)])
          grids.int.3D[[i]][,colnames(grids.int.3D[[i]]) %ni% columns.to.keep ] <- 0 
        } 
      }  
    }
    #====================== FINAL Prediction ==================================
    
    if(!hier){

      for(i in 1:length(grids.3D)){
        grids.3D[[i]]$pred <- as.numeric(predict(lasso.cv,s=lasso.cv$lambda.min, newx = grids.3D[[i]])) # Applying final model to predict on whole each grid
        grids.3D[[i]]$pred <- pmax(grids.3D[[i]]$pred,target.min/3) # Correcting the negative predictions
        grids.3D[[i]] <- grids.3D[[i]][,"pred"] # Creating of grids of final predictions
      }
      
    } else {
      
      for( i in 1:length(grids.3D)){
        grids.3D[[i]]$pred <- as.numeric(predict(hier.lasso.final, newx=grids.3D[[i]], newzz = grids.int.3D[[i]])) # Applying final model to predict on whole each grid
        grids.3D[[i]]$pred <- pmax(grids.3D[[i]]$pred,min.obs/3) # Correcting the negative predictions
        grids.3D[[i]] <- grids.3D[[i]][,"pred"] # Creating of grids of final predictions
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

  
  
  
