library(rankdist)
library(reshape2)
library(ggplot2)

#' s is a partial ranking: 0 = unassigned
sample_step = function(pi0, w, s){
    # preprocess partially filled ranking s
    loc_unassigned = which(s == 0)
    n_unassigned = length(loc_unassigned)
    n_assigned = length(s) - n_unassigned
    weight_cum = c(0, cumsum(w[(n_assigned + 1): length(w)]))
    weight_ordered = weight_cum[pi0[loc_unassigned] - n_assigned]
    # sample
    next_obj = sample(loc_unassigned, prob = exp(-weight_ordered), size = 1)
    s[next_obj] = n_assigned + 1
    # update pi0
    original_rank = pi0[next_obj]
    pi0[pi0 < original_rank & pi0 > n_assigned] = pi0[pi0 < original_rank & pi0 > n_assigned] + 1
    pi0[next_obj] = n_assigned + 1
    # recursive
    if(n_assigned + 1 == length(w)){
        s[s == 0] = length(s)
        s
    } else 
        sample_step(pi0, w, s)
}

debug(sample_step)
sample_step(pi0 = c(2, 1, 4, 3), w = c(3, 2, 1), s = c(0, 0, 0, 0))
sample_step(pi0 = c(2, 1, 4, 3), w = c(3, 2, 1) / 100, s = c(0, 0, 0, 0))
undebug(sample_step)

#' pi0 is a ranking
#' w is the weights
#' size is sample size

sample_wk = function(pi0, w, size){
    hash_table = table(replicate(size, RanktoHash(sample_step(pi0, w, rep(0, length(pi0)))), simplify = TRUE))
    list(rankmat = HashtoRank(names(hash_table)), freq = as.numeric(hash_table))
}

s1 = sample_wk(pi0 = c(2, 1, 4, 3), w = c(3, 2, 1) / 5, 50000)
s1

# use vector approach
wToparam = function(w.true){
    param.true = numeric(length(w.true))
    param.true[1:(length(w.true)-1)] = -diff(w.true)
    param.true[length(param.true)] = w.true[length(w.true)]
    param.true
}

dist_vec = as.numeric(matrix(rankdist:::CWeightGivenPi(s1$rankmat,c(2, 1, 4, 3)),nrow=nrow(s1$rankmat)) %*% wToparam(c(3, 2, 1) / 5))
prob_vector = exp(-dist_vec)
s2 = table(sample(x = length(prob_vector), prob = prob_vector, size = 50000, replace = TRUE))

plot(s1$freq)
points(s2, col = 'blue')


##### test function: SingleClusterModel (Given pi0, estimate w)

run_given_pi0_estimate_w = function(pi0, w0, sample_size){
    s1 = sample_wk(pi0 = pi0, w = w0, sample_size)
    dat1 = new('RankData', nobj = length(pi0), nobs=sample_size,ranking=s1$rankmat, count=s1$freq)
    param.init = MomentsEst(dat1, size=sample_size/10, pi0 = pi0)
    init1 = new('RankInit', param.init = list(param.init), modal_ranking.init = list(pi0), clu = 1L, p.init = 1)
    ctrl1 = new('RankControlWeightedKendall',optimx_control=list(dowarn=FALSE))
    m1 = rankdist:::SingleClusterModel(dat1,init1,ctrl1,modal_ranking=pi0)
    m1$w.est
}

set.seed(42)

## toy
# pi0 = 1:4
# w0 =  c(3, 2, 1) / 5
# sample_size = 5000

pi0 = 1:40
w0 = log(40:2) / 5
sample_size = 500

res = replicate(500, run_given_pi0_estimate_w(pi0, w0, sample_size))

res_long = melt(as.data.frame(t(res)))
res_long$variable = factor(substr(res_long$variable, 2, 10), levels = unique(as.numeric(substr(res_long$variable, 2, 10))))
true_long = data.frame(variable = 1:39, value=w0)

ggplot(data = res_long, aes(x = variable, y=value)) +
    geom_boxplot() +
    geom_point(data = true_long, aes(x = variable, y=value), col = 'red') +
    labs(x = 'w', y = '') +
    theme_classic()

ggsave(filename = 'Estimating weights.png')

##### test function: SearchPi0_fix_wt (Given weight, estimate pi0)
run_given_w_estimate_pi0 = function(pi0, w0, sample_size){
    s1 = sample_wk(pi0 = pi0, w = w0, sample_size)
    dat1 = new('RankData', nobj = length(pi0), nobs=sample_size,ranking=s1$rankmat, count=s1$freq)
    
    avg_rank <- dat1@count %*% dat1@ranking
    modal_ranking.init = OrderingToRanking(order(avg_rank))
    
    init1 = new('RankInit', param.init = list(wToparam(w0)), modal_ranking.init = list(modal_ranking.init), clu = 1L, p.init = 1)
    ctrl1 = new('RankControlWeightedKendall',optimx_control=list(dowarn=FALSE), SearchPi0_fast_traversal = FALSE, SearchPi0_neighbour = 'Kendall')
    m1 = rankdist:::SearchPi0GivenW(dat1,init1,ctrl1)
    c(m1$pi0.ranking, m1$SearchPi0_step-1, modal_ranking.init)
}

set.seed(43)

## small sample size
pi0 = 1:40
w0 = log(40:2) / 5
sample_size = 200

res2 = replicate(500, run_given_w_estimate_pi0(pi0, w0, sample_size))

res2_rank = t(res2[1:length(pi0), ])
res2_step = res2[length(pi0) + 1, ]
res2_borda = t(res2[(length(pi0) + 2):nrow(res2), ])

table(res2_step)
table(DistanceBlock(res2_rank, pi0))
table(DistanceBlock(res2_borda, pi0))

## large sample size
sample_size = 500
res3 = replicate(500, run_given_w_estimate_pi0(pi0, w0, sample_size))

res3_rank = t(res3[1:length(pi0), ])
res3_step = res3[length(pi0) + 1, ]
res3_borda = t(res3[(length(pi0) + 2):nrow(res3), ])

table(res3_step)
table(DistanceBlock(res3_rank, pi0))
table(DistanceBlock(res3_borda, pi0))


##### test function: SearchPi0 (estimate w and pi0)

run_estimate_both = function(pi0, w0, sample_size){
    s1 = sample_wk(pi0 = pi0, w = w0, sample_size)
    dat1 = new('RankData', nobj = length(pi0), nobs=sample_size,ranking=s1$rankmat, count=s1$freq)
    
    avg_rank <- dat1@count %*% dat1@ranking
    modal_ranking.init = OrderingToRanking(order(avg_rank))
    param.init = MomentsEst(dat1, size=sample_size/10, pi0 = modal_ranking.init)
    
    init1 = new('RankInit', param.init = list(param.init), modal_ranking.init = list(modal_ranking.init), clu = 1L, p.init = 1)
    ctrl1 = new('RankControlWeightedKendall',optimx_control=list(dowarn=FALSE), SearchPi0_fast_traversal = FALSE, SearchPi0_neighbour = 'Kendall')
    m1 = rankdist:::SearchPi0(dat1,init1,ctrl1)
    c(m1$w.est, m1$pi0.ranking, modal_ranking.init)
}


set.seed(44)

pi0 = 1:40
w0 = log(40:2) / 5
sample_size = 500

res4 = replicate(500, run_estimate_both(pi0, w0, sample_size))

res4_w = res4[1:(length(pi0) - 1), ]
res4_rank = t(res4[length(pi0):(2*length(pi0) - 1), ])
res4_borda = t(res4[(2*length(pi0)):nrow(res4), ])

table(DistanceBlock(res4_rank, pi0))
table(DistanceBlock(res4_borda, pi0))

res_long = melt(as.data.frame(t(res4_w)))
res_long$variable = factor(substr(res_long$variable, 2, 10), levels = unique(as.numeric(substr(res_long$variable, 2, 10))))
true_long = data.frame(variable = 1:39, value=w0)

ggplot(data = res_long, aes(x = variable, y=value)) +
    geom_boxplot() +
    geom_point(data = true_long, aes(x = variable, y=value), col = 'red') +
    labs(x = 'w', y = '') +
    theme_classic()
ggsave(filename = 'Estimating weights full small.png')

saveRDS(res4, file = 'res4.rds')

## timing
system.time( run_estimate_both(pi0, w0, sample_size))
system.time( run_given_pi0_estimate_w(pi0, w0, sample_size))
system.time( run_given_w_estimate_pi0(pi0, w0, sample_size))
