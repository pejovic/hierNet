

penint3D<-function(fun, geo, cogrids, flist,hier=TRUE,lambda=seq(0,5,0.1),int=TRUE,depth.fun=list("linear","nspline","poly"),df=4,deg=3,preProc=TRUE){
  
  if(hier==TRUE){int=TRUE}
  if(int == FALSE){hier = FALSE}
  

  ## all columns of interest:
  methodid <- all.vars(fun)[1]
  seln <- names(cogrids) %in% all.vars(fun)[-1]
  xyn <- attr(cogrids@bbox, "dimnames")[[1]]
  ## check if all covariates are available:
  if(sum(!is.na(seln))==0){
    stop("None of the covariates in the 'formulaString' matches the names in the 'covariates' object")
  }

  ## prepare regression matrix:
  ov <- over(x=cogrids, y=geo, method=methodid, var.type = "numeric")
  if(nrow(ov)==0|is.null(ov$observedValue)) {
    warning("The 'over' operations resulted in an empty set. Check 'methodid' column.")
  }
  ## geostats only possible with numeric variables:
  ov[,methodid] = as.numeric(ov$observedValue)
  ov$observedValue = NULL
  
  regdat<-ov[,all.vars(fun)[-1]]
  
  
#[,-which(names(regdat) == methodid)]
  ########## Pre-process parameters #########################
  
  if(preProc==TRUE){
    
    cont.par <- preProcess(regdat[,sapply(regdat,is.numeric)],method=c("center", "scale"))
    regdat <- predict(cont.par,regdat)
    
    dummy.par <- dummyVars(as.formula(paste("~", paste(all.vars(fun)[-1], collapse="+"))),regdat,levelsOnly=FALSE)
    regdat <- predict(dummy.par,newdata=regdat)
    
    nzv.par<-preProcess(regdat,method=c("nzv"))
    regdat<-as.data.frame(predict(nzv.par,regdat))
    head(regdat)
    
    names(regdat)<-gsub(" ",".",names(regdat))
    

    
    #aaa<-lapply(new3D, function(x) as.data.frame(x)) %>% lapply(.,function(x) predict(cont.par,newdata=x)) %>% lapply(.,function(x) predict(dummy.par,newdata=x)) %>% lapply(., function(x) as.data.frame(predict(nzv.par,newdata=x)))
    #colnames(aaa)<-gsub(" ",".",colnames(aaa))
    #XX <- hierNet::compute.interactions.c(as.matrix(aaa[[1]]),diagonal=FALSE)
    }
  
  if(depth.fun=="nspline"){fun <- as.formula(paste(methodid, "~" ,paste(c(names(regdat)[-which(names(regdat)=="altitude")],"ns(altitude,df=4)"), collapse="+"),sep=""))
  } else {if(depth.fun =="poly") {fun <- as.formula(paste(methodid, "~" ,paste(c(names(regdat)[-which(names(regdat)=="altitude")],"poly(altitude,3,raw=TRUE)"), collapse="+"),sep=""))}
    else {fun <- fun} # OVde treba ubaciti izbacivanje altituda kao main effect ako postoji poly i ns
  }
  
  regdat<-cbind(As=ov$As,regdat)
  modmat <- na.omit(regdat)#model.matrix(fun ,regdat)[,-1]
  # removing nzv varaibles 
  #nzv <- nearZeroVar(modmat)
  #if(sum(nzv)!=0){modmat <- modmat[, -nzv]}else{modmat<-modmat}
  
  colnames(modmat)<-gsub(" ",".",colnames(modmat))
  colnames(modmat)<-gsub("\\(altitude,.df.=.4\\)","",colnames(modmat))
  colnames(modmat)<-gsub("\\(altitude,.3\\)","",colnames(modmat))
  
  X <- hierNet::compute.interactions.c(as.matrix(modmat[,-1]),diagonal=FALSE)
  columns.to.keep <- (X %>% as.data.frame() %>% select(contains("altitude"))%>% colnames())
  modmat <- cbind(modmat,X[,colnames(X) %in% columns.to.keep])
  head(modmat)
  
  #modmat<-cbind(modmat,poly(modmat$altitude,3))
  
  
  
  
  
  if(int & !hier){
    X <- hierNet::compute.interactions.c(as.matrix(modmat),diagonal=FALSE)

        if(depth.fun == "nspline"){columns.to.keep <- (X %>% as.data.frame() %>% select(contains("altitude")) %>% select(-contains("ns"))  %>% colnames())
            } else {if(depth.fun =="poly") {columns.to.keep <- (X %>% as.data.frame() %>% select(contains("altitude")) %>% select(-contains("poly"))  %>% colnames())}
                    else {columns.to.keep <- (X %>% as.data.frame() %>% select(contains("altitude"))  %>% colnames())}
        }
         
        if(depth.fun =="nspline"){modmat <- cbind((modmat %>% as.data.frame() %>% select(-altitude) %>% select(-contains("ns")) %>% as.matrix()),X[,colnames(X) %in% columns.to.keep],(modmat %>% as.data.frame() %>% select(contains("ns")) %>% as.matrix()))
          } else {if(depth.fun =="poly") {modmat <- cbind((modmat %>% as.data.frame() %>% select(-altitude) %>% select(-contains("poly")) %>% as.matrix()),X[,colnames(X) %in% columns.to.keep],(modmat %>% as.data.frame() %>% select(contains("poly")) %>% as.matrix()))}
                  else {modmat <- cbind(modmat,X[,colnames(X) %in% columns.to.keep])}
                  }
  } else { if(hier) {

                    X <- hierNet::compute.interactions.c(modmat,diagonal=FALSE)
                    
                    "%ni%" <- Negate("%in%")
           
                    if(depth.fun == "nspline"){columns.to.keep <- (X %>% as.data.frame() %>% select(contains("altitude")) %>% select(-contains("ns"))  %>% colnames())
              
                                              } else {if(depth.fun == "poly") {columns.to.keep <- (X %>% as.data.frame() %>% select(contains("altitude")) %>% select(-contains("poly"))  %>% colnames())}
                                                
                                                else {columns.to.keep <- (X %>% as.data.frame() %>% select(contains("altitude"))  %>% colnames())}
                                              }
           
                    X[,colnames(X) %ni% columns.to.keep ] <- 0
                    

                    } else {
                      
                     
                      if(depth.fun =="nspline"){modmat <- cbind((modmat %>% as.data.frame() %>% select(-altitude) %>% as.matrix()))
                      } else {if(depth.fun =="poly") {modmat <- cbind((modmat %>% as.data.frame() %>% select(-altitude) %>% as.matrix()))}
                        else {modmat <- modmat}
                      }

                    }
                    
              }
    
  if(!hier){
    
    allData<-cbind(tvar=regdat[,paste(tvar)],modmat,regdat[,c("ID","longitude","latitude")])
    #names(allData)<-gsub("\\(altitude,.df.=.4\\)","",names(allData))
    
    results<-data.frame(lambda=rep(NA,length(flist)+1),RMSE=rep(NA,length(flist)+1),Rsquared=rep(NA,length(flist)+1))
    coef.list=as.list(rep(NA,length(strat)))
    pred<-data.frame()
    for(i in 1:length(flist)){
      ind<-which(allData$ID %in% do.call(c,flist[-i]))
      TrainData<-as.data.frame(do.call(cbind,allData[ind,]))
      Train.ID.index<-flist[-i]
      TestData<-as.data.frame(do.call(cbind,allData[-ind,]))
      Test.ID.Index<-flist[i]
      
      folds1<-length(Train.ID.index)
      folds.list1<-as.list(rep(NA,folds1))
      names(folds.list1)<-paste("fold",c(1:folds1),sep = "")
      foldid<-rep(NA,dim(TrainData)[1])
      for(j in 1:length(Train.ID.index)){
        folds.list1[[j]]<-which(TrainData$ID %in% Train.ID.index[[j]])
        foldid[folds.list1[[j]]]<-j
      }
      
      TrainData<-TrainData[,1:(dim(TrainData)[2]-3)] # ovo mora da se menja...ne sme da stoji tri vec moraju da se uklone poslednje tri kolone lepse
      TestData<-TestData[,1:(dim(TestData)[2]-3)]
      #lambdaGrid <- seq(0,2.5,0.05)#10^seq(10,-2, length =100)
      lasso.mod=cv.glmnet(as.matrix(TrainData[,-1]),TrainData[,1],alpha=1,lambda=lambda,foldid=foldid,type.measure="mse")
      lasso.pred<-predict(lasso.mod,s=lasso.mod$lambda.min,newx=as.matrix(TestData[,-1]))
      obs.pred<-data.frame(obs=TestData[,1],pred=as.numeric(lasso.pred))
      coef.list[[i]]<-predict(lasso.mod,type="coefficients",s=lasso.mod$lambda.min)
      dfresults<-data.frame(lambda=lasso.mod$lambda.min,RMSE=defaultSummary(obs.pred)[1],Rsquared=defaultSummary(obs.pred)[2])
      results[i,]<-dfresults
      pred<-rbind(pred,obs.pred)
    }
    results[length(flist)+1,]<-c(NA,RMSE=defaultSummary(pred)[1],Rsquared=defaultSummary(pred)[2])
    coef.mat<-do.call(cbind,coef.list)
    out<-list(measure=results,coef=coef.mat)
    return(out)
    
  }else{
    

    allData<-cbind(tvar=regdat[,paste(tvar)],modmat,regdat[,c("ID","longitude","latitude")],X)
    #names(allData)<-gsub("\\(altitude,.df.=.4\\)","",names(allData))
    head(allData)
    
    results<-data.frame(lambda=rep(NA,length(flist)+1),RMSE=rep(NA,length(flist)+1),Rsquared=rep(NA,length(flist)+1))
    coef.list=as.list(rep(NA,length(flist)))
    pred<-data.frame()
    
    for(i in 1:length(flist)){
      ind<-which(allData$ID %in% do.call(c,flist[-i])) # izdvajaje indeksa instanci koje pribadaju foldovima za trening
      TrainData<-as.data.frame(do.call(cbind,allData[ind,])) # izdvajanje trening podacima
      Train.ID.index<-flist[-i] # Izdvajanje dela liste foldova sa  trening podacima
      TestData<-as.data.frame(do.call(cbind,allData[-ind,])) #Izdvajanje test podataka
      Test.ID.Index<-flist[i] # Izdvajanje dela liste foldova sa test podacima
      
      folds1<-length(Train.ID.index) # duzina liste foldova sa trening podacima (tj. koliko foldova ide u tu listu)
      folds.list1<-as.list(rep(NA,folds1)) # kreiranje prazne liste iste duzine
      names(folds.list1)<-paste("fold",c(1:folds1),sep = "") # Nazivanje svakog clana liste po foldu...
      foldid<-rep(NA,dim(TrainData)[1]) # kreiranje atributa koji ce nositi informaciju kom foldu pripada podatak 
      
      for(j in 1:length(Train.ID.index)){
        folds.list1[[j]]<-which(TrainData$ID %in% Train.ID.index[[j]])
        foldid[folds.list1[[j]]]<-j
      }
      
      trainx <- as.matrix(TrainData[,colnames(modmat)]) 
      testx <- as.matrix(TestData[,colnames(modmat)])
      trainzz <- as.matrix(TrainData[,colnames(X)])
      testzz <- as.matrix(TestData[,colnames(X)])
      trainy <- (TrainData[,"tvar"])
      testy <- (TestData[,"tvar"])
      
      fit = hierNet.path(trainx,trainy, zz = trainzz,diagonal=FALSE,strong=TRUE,trace=0)
      fitcv = hierNet.cv(fit, trainx, trainy, folds=folds.list1, trace=0)
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
