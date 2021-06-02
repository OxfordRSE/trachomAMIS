Data = read.csv("./tests/test_data/geostatistical_prevalence_data.csv")

prev = matrix(NA, ncol = n.map.sampl, nrow = nrow(Data))
for (i in 1:nrow(Data))
{
  set.seed(36)
  L = rnorm(n.map.sampl, Data$Logit[i], sd = Data$Sds[i])
  prev[i, ] = exp(L)/(1+exp(L))
}
prev <- prev*100
write.table(prev, file = "./tests/test_data/prevalence_map.csv",
            row.names = F, col.names = F)
