plotfolds<-function(folds,targetvar){
  allData<-folds$Data
  targetVar<-targetvar
  allData.unique<-ddply(allData,.(ID),here(summarize),target=mean(eval(parse(text=targetVar))),longitude=longitude[1],latitude=latitude[1],fold=fold[1])
  q <- ggplot(allData.unique,aes(x = longitude, y = latitude))
  r <- q +geom_point(aes(size = sqrt(target/pi)), pch = 21, show.legend = FALSE) + scale_size_continuous(range=c(1,10))
  r <- r + facet_wrap(~ fold)
  r <- r + aes(fill = fold)
  plot(r)

}