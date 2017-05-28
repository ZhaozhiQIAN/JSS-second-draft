library(Rankcluster)

data("APA")
str(APA)

set.seed(42)
model_1 = rankclust(data=APA$data, K=1, detail=TRUE, Qsem=1000, Bsem=100, Ql=500, Bl=50, maxTry=20, run=10)
summary(model_1)

model_2 = rankclust(data=APA$data, K=2, detail=TRUE, Qsem=1000, Bsem=100, Ql=500, Bl=50, maxTry=20, run=10)
summary(model_2)

model_3 = rankclust(data=APA$data, K=3, detail=TRUE, Qsem=1000, Bsem=100, Ql=500, Bl=50, maxTry=20, run=10)
summary(model_3)

set.seed(41)
model_4 = rankclust(data=APA$data, K=4, detail=TRUE, Qsem=1000, Bsem=100, Ql=500, Bl=50, maxTry=20, run=10)
summary(model_4)
# 54335.77

# (54335.77 - 2*27120.28) / log(nrow(APA$data))

# full rankings -----------------------------------------------------------

# re-format data set to Rankcluster
library(rankdist)

str(apa_partial_obj)

head(apa_partial_obj@ranking)

rank_mat_compact = matrix(ncol = 5, nrow = 0)

for (i in 1:4){
    tmp_mat = apa_partial_obj@ranking[apa_partial_obj@q_ind[i]:(apa_partial_obj@q_ind[i + 1] - 1), ]
    if (i < 4) tmp_mat[tmp_mat == max(tmp_mat)] = NA
    rank_mat_compact = rbind(rank_mat_compact, tmp_mat)
}

nrow(rank_mat_compact)
nrow(apa_partial_obj@ranking)

head(rank_mat_compact)
head(apa_partial_obj@ranking)

tail(rank_mat_compact)
tail(apa_partial_obj@ranking)

rank_mat_full = matrix(ncol = 5, nrow = 0)
for(i in 1:length(apa_partial_obj@count)){
    r = rank_mat_compact[i, ]
    tmp_mat = matrix(data = rep(r, apa_partial_obj@count[i]), ncol = 5, byrow = TRUE)
    rank_mat_full = rbind(rank_mat_full, tmp_mat)
}

nrow(rank_mat_full) == sum(apa_partial_obj@count)
head(rank_mat_full)

# fit models
set.seed(39)
model_1_f = rankclust(data=rank_mat_full, K=1, detail=TRUE)
summary(model_1_f)
# 147363.3

set.seed(38)
model_2_f = rankclust(data=rank_mat_full, K=2, detail=TRUE)
summary(model_2_f)
# 145556.4

set.seed(37)
model_3_f = rankclust(data=rank_mat_full, K=3, detail=TRUE)
summary(model_3_f)
# 144774.2
# (144774.2 - 2*72348.53) / log(nrow(rank_mat_full))

set.seed(36)
model_4_f = rankclust(data=rank_mat_full, K=4, detail=TRUE)
summary(model_4_f)
# 146227.9
