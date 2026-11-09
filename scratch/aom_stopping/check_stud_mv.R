# MV studentisation check: does the observed REML pin hold element-wise
# (sum_j q~_j q~_j' = sum_j V_j), and rejection rates at several alpha.
R <- read.csv("scratch/aom_stopping/results/stage1_stud/dh_t3_k1.csv")
O <- R[R$status == "ok", ]
cat(sprintf("MV us: max |S.obs - S.target| / max diag(S.target): median %.4f, 90%% %.4f, max %.4f\n",
            median(O$obs.Sdev), quantile(O$obs.Sdev, .9), max(O$obs.Sdev)))
for (a in c(.01, .05, .10, .20))
    cat(sprintf("alpha %.2f:  original %.3f   scalar-stud %.3f   matrix-stud %.3f\n",
                a, mean(O$p.fwe <= a), mean(O$p.fwe.stud <= a), mean(O$p.fwe.studM <= a)))
