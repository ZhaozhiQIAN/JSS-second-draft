library(rankdist)
set.seed(1080)
gen1<- GenerateExample(ranking=TRUE)
tail(gen1$ranking)
dat1 <- new("RankData",ranking=gen1$ranking,count=gen1$count)

gen2 <- GenerateExample(ranking=FALSE)
tail(gen2$ordering)
dat2 <- new("RankData",ordering=gen2$ordering,count=gen2$count)

str1 <- MomentsEst(dat1,500)
init1 <- new("RankInit",param.init=list(str1),modal_ranking.init=list(1:5),clu=1L)
init1c <- new("RankInit",param.init=list(str1,str1),modal_ranking.init=list(1:5,5:1),clu=2L)

ctrl1 <- new("RankControlKendall",SearchPi0_show_message=FALSE)
model1 <- RankDistanceModel(dat1,init1,ctrl1)
ModelSummary(model1)
model1c <- RankDistanceModel(dat1,init1c,ctrl1)
ModelSummary(model1c)

genq <- GenerateExampleTopQ()
head(genq$ranking)
datq <- new("RankData", ranking=genq$ranking, count=genq$count,nobj=5, topq=3)

initq <- new("RankInit",param.init=list(rep(0.5,3)),modal_ranking.init=list(1:5),clu=1L)
ctrlq <- new("RankControlWeightedKendall",SearchPi0_show_message=FALSE)
modelq <- RankDistanceModel(datq,initq,ctrlq)
ModelSummary(modelq)






