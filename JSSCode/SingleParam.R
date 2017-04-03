apa_s_ctrl<-new("RankControlSpearman",SearchPi0_fast_traversal=TRUE,SearchPi0_show_message=FALSE)
Spearman_model <- list()
for (i in 1L:8L){
    init<-new("RankInit",param.init=replicate(i,0.1,FALSE),modal_ranking.init=replicate(i,sample(5,5),FALSE),clu=i,p.init=rep(1,i)/i)
    Spearman_model[[i]]<-RankDistanceModel(apa_obj,init,apa_s_ctrl)
}
dummy<-sapply(Spearman_model,function(x)cat('BIC:',x$BIC,'\tSSR',x$SSR,'\n'))
ModelSummary(Spearman_model[[3]])

apa_s_ctrl<-new("RankControlFootrule",SearchPi0_fast_traversal=TRUE,SearchPi0_show_message=FALSE)
Footrule_model <- list()
for (i in 1L:8L){
    init<-new("RankInit",param.init=replicate(i,0.1,FALSE),modal_ranking.init=replicate(i,sample(5,5),FALSE),clu=i,p.init=rep(1,i)/i)
    Footrule_model[[i]]<-RankDistanceModel(apa_obj,init,apa_s_ctrl)
}
dummy<-sapply(Footrule_model,function(x)cat('BIC:',x$BIC,'\tSSR',x$SSR,'\n'))
ModelSummary(Footrule_model[[4]])


apa_s_ctrl<-new("RankControlHamming",SearchPi0_fast_traversal=TRUE,SearchPi0_show_message=FALSE)
Hamming_model <- list()
for (i in 1L:8L){
    init<-new("RankInit",param.init=replicate(i,0.1,FALSE),modal_ranking.init=replicate(i,sample(5,5),FALSE),clu=i,p.init=rep(1,i)/i)
    Hamming_model[[i]]<-RankDistanceModel(apa_obj,init,apa_s_ctrl)
}
dummy<-sapply(Hamming_model,function(x)cat('BIC:',x$BIC,'\tSSR',x$SSR,'\n'))
ModelSummary(Hamming_model[[7]])

apa_s_ctrl<-new("RankControlCayley",SearchPi0_fast_traversal=TRUE,SearchPi0_show_message=FALSE)
Cayley_model <- list()
for (i in 1L:8L){
    init<-new("RankInit",param.init=replicate(i,0.1,FALSE),modal_ranking.init=replicate(i,sample(5,5),FALSE),clu=i,p.init=rep(1,i)/i)
    Cayley_model[[i]]<-RankDistanceModel(apa_obj,init,apa_s_ctrl)
}
dummy<-sapply(Cayley_model,function(x)cat('BIC:',x$BIC,'\tSSR',x$SSR,'\n'))
ModelSummary(Cayley_model[[6]])
