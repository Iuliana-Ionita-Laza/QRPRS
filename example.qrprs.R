library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
options(scipen = 999)



## phenotype for training set
pheno.train = fread(file = paste0("example/example.train.phenotype.tsv"),
                    header = T, sep = "\t", data.table = F, stringsAsFactors = F)



## quantile levels
tau.pool = as.character(c(0.025, (1:9) / 10, 0.975))



## QR-PRS for training set
prs.train = data.frame()
for (tau in tau.pool) {
  tmp = fread(file = paste0("example/example.train.qrprs.tau", tau, ".sscore"),
              header = T, sep = "\t", data.table = F, stringsAsFactors = F)
  prs.train = rbind(prs.train, 
                    tmp %>% select(`#IID`, SCORE1_SUM) %>% mutate(Tau = tau, .before = 1))
}
prs.train = prs.train %>% rename(IID = `#IID`, QRPRS = SCORE1_SUM) 



## conformal prediction
df.conform = data.frame(Rho = 0, Tau = "0.5", Calibration = 0, stringsAsFactors = F)
rho.pool = c(0.2, 0.4, 0.6, 0.8, 0.95)
for (rho in rho.pool) {
  tau.lo = as.character((1 - rho) / 2)
  tau.hi = as.character((1 + rho) / 2)
  error.lo = pheno.train %>% left_join(y = prs.train %>% filter(Tau == tau.lo), by = c("IID" = "IID")) %>% pull(QRPRS)
  error.hi = pheno.train %>% left_join(y = prs.train %>% filter(Tau == tau.hi), by = c("IID" = "IID")) %>% pull(QRPRS)
  error.true = pheno.train %>% pull(Trait)
  error.lo = error.lo - error.true
  error.hi = error.true - error.hi
  error = pmax(error.lo, error.hi)
  set.seed(seed = 256)
  idx.clbr = sample(1:length(error), size = round(length(error) / 2), replace = F)
  error.clbr = sort(error[idx.clbr])[ceiling(rho * (1 + length(idx.clbr)))]
  df.conform = df.conform %>% 
    add_row(Rho = rho, Tau = tau.lo, Calibration = -error.clbr) %>% 
    add_row(Rho = rho, Tau = tau.hi, Calibration = error.clbr)
}
df.conform = df.conform %>% mutate(Tau = as.character(Tau))



## phenotype for test set
pheno.test = fread(file = paste0("example/example.test.phenotype.tsv"),
                   header = T, sep = "\t", data.table = F, stringsAsFactors = F)



## QR-PRS for test set
prs.test = data.frame()
for (tau in tau.pool) {
  tmp = fread(file = paste0("example/example.test.qrprs.tau", tau, ".sscore"),
              header = T, sep = "\t", data.table = F, stringsAsFactors = F)
  prs.test = rbind(prs.test, 
                   tmp %>% select(`#IID`, SCORE1_SUM) %>% mutate(Tau = tau, .before = 1))
}
prs.test = prs.test %>% rename(IID = `#IID`, QRPRS = SCORE1_SUM) 



## conformalized QR-PRS (CQR-PRS) for test set
prs.test = prs.test %>% left_join(y = df.conform %>% select(Tau, Calibration), by = c("Tau" = "Tau"))
prs.test = prs.test %>% mutate(CQRPRS = QRPRS + Calibration) %>% select(-Calibration)



## risk stratification for test set 
## at threshold t = 90th percentile (top 10%)
th = quantile(x = prs.test %>% filter(Tau == "0.5") %>% pull(QRPRS), probs = 0.9) %>% unname()
strat.test = prs.test %>% pivot_wider(id_cols = IID, names_from = Tau, names_prefix = "QRPRS", values_from = QRPRS)
strat.test.gt = strat.test %>% 
  filter(QRPRS0.5 > th) %>% 
  mutate(N = rowSums(select(., starts_with("QRPRS")) > th)) %>% 
  mutate(Category = if_else(N == length(tau.pool), true = "Certain-above t", false = "Uncertain-above t"))
strat.test.lt = strat.test %>% 
  filter(QRPRS0.5 <= th) %>% 
  mutate(N = rowSums(select(., starts_with("QRPRS")) <= th)) %>% 
  mutate(Category = if_else(N == length(tau.pool), true = "Certain-below t", false = "Uncertain-below t"))
strat.test = rbind(strat.test.gt, strat.test.lt) %>% select(-N)
strat.test = strat.test %>% mutate(Category = factor(x = Category, 
                                                     levels = c("Certain-above t", "Uncertain-above t", "Uncertain-below t", "Certain-below t")))
print(strat.test %>% count(Category) %>% complete(Category, fill = list(n = 0)))



## coverage of prediction interval for test set
df.prs.test = pheno.test %>% 
  left_join(y = prs.test %>% pivot_wider(id_cols = IID, names_from = Tau, names_sep = "", values_from = c(QRPRS, CQRPRS)),
            by = c("IID" = "IID"))
qrprs.coverage.test = (df.prs.test %>% filter((QRPRS0.025 <= Trait & Trait <= QRPRS0.975)) %>% nrow()) / (df.prs.test %>% nrow())
cqrprs.coverage.test = (df.prs.test %>% filter((CQRPRS0.025 <= Trait & Trait <= CQRPRS0.975)) %>% nrow()) / (df.prs.test %>% nrow())
print(paste0("Coverage of QR-PRS prediction intervals: ", sprintf("%.1f%%", qrprs.coverage.test * 100)))
print(paste0("Coverage of CQR-PRS prediction intervals: ", sprintf("%.1f%%", cqrprs.coverage.test * 100)))


