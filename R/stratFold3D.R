# Create stratified folds
# Treba uraditi prosek uzimajuci u obzir sirinu horizonta...

stratfold3d <- function(targetVar,regdat,folds=6,cent=3,preProc=TRUE,seed=666,dimensions=list("2D","3D"),IDs=TRUE,sum=FALSE){
  
  dimensions<-dimensions[[1]]
  if(dimensions=="2D"){
    
    unique.df<-ddply(regdat,.(ID),here(summarize),target=weighted.mean(eval(parse(text=targetVar)),hdepth),longitude=longitude[1],latitude=latitude[1])
    km <- kmeans(cbind(unique.df$longitude,unique.df$latitude), centers = cent)
    #plot(unique.df$longitude,unique.df$latitude, col = km$cluster, pch = 20)
    unique.df$km<-as.factor(km$cluster)} else {
      
      unique.df<-ddply(regdat,.(ID),here(summarize),target=weighted.mean(eval(parse(text=targetVar)),hdepth),longitude=longitude[1],latitude=latitude[1],altitude=min(altitude))
      km <- kmeans(cbind(scale(unique.df$longitude),scale(unique.df$latitude),scale(unique.df$altitude)), centers = cent)
      #plot(unique.df$longitude,unique.df$latitude, col = km$cluster, pch = 20)
      unique.df$km<-as.factor(km$cluster)
    }
  
  
  k.list<-as.list(rep(NA,length(unique(unique.df$km))))
  names(k.list)<-paste("k",c(1:length(k.list)),sep="")
  
  
  ########## Creating list of profile indices of 6 folds per cluster region ###############
  for(i in 1:length(k.list)){
    set.seed(seed)
    k.list[[i]]<-createFolds(unique.df[which(unique.df$km==levels(unique.df$km)[i]),"target"],k=folds)
    for(j in 1:folds){
      k.list[[i]][[j]]<-unique.df[which(unique.df$km==levels(unique.df$km)[i]),"ID"][k.list[[i]][[j]]]
    }
  }
  ##########################################################################################
  
  ################### Creating list of profile indices per fold ############################
  ID.list<-as.list(rep(NA,folds))
  names(ID.list)<-paste("fold",c(1:folds),sep = "")
  for(i in 1:folds){
    ID.list[[i]]<-do.call(c,lapply(k.list,function(x) x[[i]]))
    names(ID.list[[i]])<-NULL
    #folds.list[[i]]<-as.character(folds.list[[i]])
  }
  
  ###########################################################################################
  
  ##############  Creating list of 
  folds.list<-as.list(rep(NA,folds))
  names(folds.list)<-paste("fold",c(1:folds),sep = "")
  for(i in 1:length(ID.list)){
    folds.list[[i]]<-which(regdat$ID %in% ID.list[[i]])
  }
  
  
  pom<-data.frame()
  for(i in 1:length(folds.list)){
    allData1<-regdat[folds.list[[i]],]
    allData1$fold<-paste("fold",i,sep="")
    allData<-rbind(allData1,pom)
    pom<-allData
  }
  allData$fold<-factor(allData$fold)
  sum.list=list(by(allData$altitude,allData$fold,summary))
  if(IDs==TRUE){index.list=(ID.list)}else{index.list=(folds.list)}
  sum.list<-list(allData,index.list,sum.list,by(allData[,paste(targetVar)],allData$fold,summary))
  names(sum.list)<-c("Data","folds","altitude summary",paste(targetVar,"summary", sep=" "))

  if(sum==TRUE){return(sum.list)}else{return(sum.list[[1]])}
}


