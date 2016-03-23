
gridfun<-as.formula(paste("~", paste(c(sm2D.lst,"altitude"), collapse="+")))

spc<-function(grids,fun,contVar,preProc==TRUE){
  
  formulaString<-fun
  obj<-grids
  
  tvar<-as.character(formulaString)[2]
  covs<-gsub("[[:space:]]", "", as.character(formulaString)[3])
  covs<-unlist(strsplit(covs, "+", fixed = TRUE))
  depth<-covs[length(covs)]
  covs<-covs[-length(covs)]
  
  formulaString <- as.formula(paste("~", paste(c(sm2D.lst), collapse="+")))
  
  vars = all.vars(formulaString)
  
  obj@data <- obj@data[,vars]
  
  if(nrow(obj)>10e6){
    warning('Operation not recommended for large grids', immediate. = TRUE)
  }
  
  if(preProc==TRUE){
    obj@data[,c(contVar)]<-apply(obj@data[,c(contVar)],2,scale)
  }
  
  modmat <- model.matrix(fun ,obj@data)[,-1]
  # removing nzv varaibles 
  nzv <- nearZeroVar(modmat)
  if(sum(nzv)!=0){modmat <- modmat[, -nzv]}else{modmat<-modmat}
  
  }

formulaString=fm.H
obj<-gridmaps.sm2D

if(missing(formulaString)) {
  formulaString <- as.formula(paste("~", paste(c(sm2D.lst), collapse="+")))
}
vars = all.vars(formulaString)
if(length(vars)< 2){
  stop("At least two covarites required to run Principal Component Analysis")
}
obj@data <- obj@data[,vars]
## print warning:
if(silent==FALSE){
  if(nrow(obj)>10e6){
    warning('Operation not recommended for large grids', immediate. = TRUE)
  }}
## convert every factor to indicators:
for(j in 1:length(vars)){
  if(is.factor(obj@data[,vars[j]])){
    # remove classes without pixels:
    obj@data[,vars[j]] <- as.factor(paste(obj@data[,vars[j]]))
    ln <- levels(obj@data[,vars[j]])
    for(k in 1:length(ln)){
      vn <- paste(vars[j], k, sep="_")
      obj@data[,vn] <- ifelse(obj@data[,vars[j]]==ln[k], 1, 0)
    }
    message(paste("Converting", vars[j], "to indicators..."))
  }
}
varsn = names(obj)[which(!sapply(obj@data, is.factor))]
obj@data <- obj@data[,varsn]
## filter the missing values:
if(scale. == TRUE){
  x <- scale(obj@data)
  x[is.na(x)] <- 0
  x <- as.data.frame(x)
  sd.l <- lapply(x, FUN=sd)
  x0 <- sd.l==0
  if(any(x0)){
    message(paste("Columns with zero variance removed:", names(x)[which(x0)]), immediate. = TRUE)
    formulaString.f = as.formula(paste("~", paste(varsn[-which(x0)], collapse="+")))
    ## principal component analysis:
    pcs <- prcomp(formula=formulaString.f, x)
  } else {
    formulaString = as.formula(paste("~", paste(varsn, collapse="+")))
    pcs <- prcomp(formula=formulaString, x)
  }
} else {
  formulaString = as.formula(paste("~", paste(varsn, collapse="+")))
  pcs <- prcomp(formula=formulaString, obj@data)
}
## copy values:
obj@data <- as.data.frame(pcs$x)
proj4string(obj) <- obj@proj4string
if(silent==FALSE){
  message(paste("Converting covariates to principal components..."))
  summary(pcs)
}
pcs <- new("SpatialComponents", predicted = obj, pca = pcs[-which(names(pcs)=="x")])
return(pcs)
}







stdepths <- c(-.1,-.3,-.5)
new3D <- sp3D(gridmaps.sm2D, stdepths=stdepths)
str(new3D)

grid1<-as.data.frame(new3D[[1]])

cogrids3D <- lapply(new3D, function(x) as.data.frame(x)) %>% lapply(.,function(x) predict(cont.par,newdata=x)) %>% lapply(.,function(x) predict(alt.par,newdata=x)) %>% lapply(.,function(x) predict(dummy.par,newdata=x)) %>% lapply(., function(x) as.data.frame(predict(nzv.par,newdata=x))) %>% lapply(., function(x) {colnames(x) <- gsub( "\\_|/|\\-|\"|\\s" , "." , colnames(x) );return(x)})
XX <- cogrids3D[[1]]
XX <- apply(XX, 1, function(x) hierNet::compute.interactions.c(t(as.matrix(x)),diagonal=FALSE)