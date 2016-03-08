## 3D Penalized regression with interactions 
## Ubaciti hier==TRUE parametar pa ponoviti petlju, razdvojiti ova dva pristupa.


penint3D<-function(fun,regdat,contVar,flist,lambda=seq(0,5,0.1),int=TRUE,depth.fun=list("linear","nspline","poly"),df=4,deg=3,preProc=TRUE){
  
  tvar<-as.character(fun)[2]
  covs<-gsub("[[:space:]]", "", as.character(fun)[3])
  covs<-unlist(strsplit(covs, "+", fixed = TRUE))
  depth<-covs[length(covs)]
  covs<-covs[-length(covs)]
  
  
  if(int==TRUE){
    fun<-as.formula(paste(tvar, "~" ,paste(covs,depth,sep="*", collapse="+"),sep=""))
  } else {fun<-fun}
  
  
  
  if(preProc==TRUE){
    preProcValues <- preProcess(regdat[,contVar], method = c("center", "scale"))
    Trans.df<-cbind(predict(preProcValues, regdat[,contVar]),regdat[,-which(names(regdat) %in% contVar)])
  }
  
  modmat <- model.matrix(fun ,Trans.df)[,-1] #fm.int.lm.As, # fm.GSIF.int.lm.As
  # removing nzv varaibles 
  nzv <- nearZeroVar(modmat)
  if(sum(nzv)!=0){modmat <- modmat[, -nzv]}else{modmat<-modmat}
  
  # removing correlated variables 
  corr_mat <- cor(modmat)
  too_high <- findCorrelation(corr_mat, cutoff = .9)
  if(sum(too_high)!=0){modmat <- modmat[, -too_high]}else{modmat<-modmat}
  
  
  allData<-cbind(regdat[,paste(tvar)],modmat,regdat[,c("ID","longitude","latitude")])
  names(allData)<-gsub("\\(altitude,.df.=.4\\)","",names(allData))
  
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
    
    TrainData<-TrainData[,1:(dim(TrainData)[2]-3)]
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
}
