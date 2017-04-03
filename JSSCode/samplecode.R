######################## Fit APA data #############################
## Weighted Kendall distance for complete votes
library(rankdist)
set.seed(666)
apa_wk_ctrl<-new("RankControlWeightedKendall",SearchPi0_fast_traversal=TRUE,SearchPi0_show_message=FALSE)
geoKendall_model<-list()
for (i in 1L:3L){
  init<-new("RankInit",param.init=replicate(i,rep(0.1,4),FALSE),modal_ranking.init=replicate(i,sample(5,5),FALSE),clu=i,p.init=rep(1,i)/i)
  geoKendall_model[[i]]<-RankDistanceModel(apa_obj,init,apa_wk_ctrl)
}
dummy<-sapply(geoKendall_model,function(x)cat('BIC:',x$BIC,'\tSSR',x$SSR,'\n'))
ModelSummary(geoKendall_model[[3]])
## visualization
d1 <- DistanceMatrix(apa_obj@ranking)
sc <- scale(apa_obj@count,center=FALSE)
sc <- sc*6
pc <- log(apa_wk_c3$expectation)
pc <- round((pc-min(pc))/(max(pc)-min(pc))*100)
hc=heat.colors(100, alpha=1)
hc <- rev(hc)
cols <- hc[pc]
# library(png)
# png(filename="MDS_APA.png",width=800,height=800)
plot(cmdscale(d1),lwd=sc,pch=19,col=cols,
     main="Multidimensional scaling plot for APA data",
     xlab="x",ylab="y")
points(x=cmdscale(d1)[89,1],y=cmdscale(d1)[89,2]-0.3,col="blue",pch='A',lwd=10)
points(x=cmdscale(d1)[56,1],y=cmdscale(d1)[56,2]-0.2,col="blue",pch='B',lwd=10)
points(x=cmdscale(d1)[37,1],y=cmdscale(d1)[37,2]-0.2,col="blue",pch='C',lwd=10)
# dev.off()

## Weighted Kendall distance for all votes
# tied-rank
apa_wk_ctrl@assumption <- "tied-rank"
geoKendall_model_q<-list()
for (i in 1L:3L){
    init<-new("RankInit",param.init=replicate(i,rep(0.1,4),FALSE),modal_ranking.init=replicate(i,sample(5,5),FALSE),clu=i,p.init=rep(1,i)/i)
    geoKendall_model_q[[i]]<-RankDistanceModel(apa_partial_obj,init,apa_wk_ctrl)
}
dummy<-sapply(geoKendall_model_q,function(x)cat('BIC:',x$BIC,'\tSSR',x$SSR,'\n'))
ModelSummary(geoKendall_model_q[[3]])
# equal-probability
apa_wk_ctrl@assumption <- "equal-probability"
geoKendall_model_eq<-list()
for (i in 1L:3L){
    init<-new("RankInit",param.init=replicate(i,rep(0.5,4),FALSE),modal_ranking.init=replicate(i,sample(5,5),FALSE),clu=i,p.init=rep(1,i)/i)
    geoKendall_model_eq[[i]]<-RankDistanceModel(apa_partial_obj,init,apa_wk_ctrl)
}
dummy<-sapply(geoKendall_model_eq,function(x)cat('BIC:',x$BIC,'\tSSR',x$SSR,'\n'))
ModelSummary(geoKendall_model_eq[[3]])

## Mallow's Phi (Kendall distance) for complete votes
apa_k_ctrl<-new("RankControlKendall",SearchPi0_fast_traversal=TRUE,SearchPi0_show_message=FALSE)
Kendall_model<-list()
for (i in 1L:8L){
  init<-new("RankInit",param.init=replicate(i,0.1,FALSE),modal_ranking.init=replicate(i,sample(5,5),FALSE),clu=i,p.init=rep(1,i)/i)
  Kendall_model[[i]]<-RankDistanceModel(apa_obj,init,apa_k_ctrl)
}
dummy<-sapply(Kendall_model,function(x)cat('BIC:',x$BIC,'\tSSR',x$SSR,'\n'))
ModelSummary(Kendall_model[[5]])

## Phi-component for complete votes
apa_pc_ctrl<-new("RankControlPhiComponent",SearchPi0_fast_traversal=TRUE,SearchPi0_show_message=FALSE)
pc_model<-list()
for (i in 1L:3L){
  init<-new("RankInit",param.init=replicate(i,rep(0.1,4),FALSE),modal_ranking.init=replicate(i,sample(5,5),FALSE),clu=i,p.init=rep(1,i)/i)
  pc_model[[i]]<-RankDistanceModel(apa_obj,init,apa_pc_ctrl)
}
dummy<-sapply(pc_model,function(x)cat('BIC:',x$BIC,'\tSSR',x$SSR,'\n'))
ModelSummary(pc_model[[3]])

######################### Fit Cognitive Data ###########################
load("CogSurvey.RData")
set.seed(777)
cog_wk_ctrl<-new("RankControlWeightedKendall",SearchPi0_fast_traversal=TRUE,SearchPi0_show_message=FALSE)
autoInit = function(dat){
  init<-new("RankInit",param.init=list(rep(0.1,dat@nobj-1)),modal_ranking.init=list(sample(dat@nobj,dat@nobj)),clu=1L)
  init
}
rank_init_list = lapply(rank_data_list,autoInit)
agg_model = mapply(RankDistanceModel,rank_data_list,rank_init_list,MoreArgs=list(cog_wk_ctrl),SIMPLIFY=FALSE)
dist_to_truth = function(dat,model){
  dist = list()
  truth = seq_len(dat@nobj)
  dist$obs_truth = DistanceBlock(dat@ranking,truth)
  borda = order(order(colSums(dat@ranking)))
  dist$borda_truth = DistancePair(truth,borda)
  dist$model_truth = DistancePair(model$modal_ranking.est[[1]],truth)
  dist
}
d <- mapply(dist_to_truth,rank_data_list,agg_model,SIMPLIFY = FALSE)
BetterObs = sapply(d,function(x) sum(x$obs_truth<x$model_truth)/length(x$obs_truth))
DistToTruth = sapply(d,function(x)x$model_truth)
BordaToTruth = sapply(d,function(x)x$borda_truth)
