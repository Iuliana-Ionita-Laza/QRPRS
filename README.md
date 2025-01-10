# Quantifying uncertainty in polygenic risk score prediction using quantile regression (QR-PRS)



## Install dependent packages in R

data.table [https://CRAN.R-project.org/package=data.table](https://CRAN.R-project.org/package=data.table)

dplyr [https://CRAN.R-project.org/package=dplyr](https://CRAN.R-project.org/package=dplyr)

quantreg [https://cran.r-project.org/package=quantreg](https://cran.r-project.org/package=quantreg)

stringr [https://CRAN.R-project.org/package=stringr](https://CRAN.R-project.org/package=stringr)

tidyr [https://CRAN.R-project.org/package=tidyr](https://CRAN.R-project.org/package=tidyr)



## QR-PRS and prediction intervals



### Step 1

Perform QR-GWAS on the training dataset using the code from [https://github.com/Iuliana-Ionita-Laza/QRGWAS](https://github.com/Iuliana-Ionita-Laza/QRGWAS). Set ```is.effect.estimated = T``` to enable the estimation of quantile-specific effect size in the R script for QR-GWAS. The expected output includes the quantile-specific p-value and effect size (beta) at different quantile levels for each variant. For downstream QR-PRS and risk uncertainty analysis, we suggest considering 11 quantile levels: 0.025, 0.1, 0.2, ..., 0.9, and 0.975 when running GWAS by setting ```qntl = c(0.025, (1:9) / 10, 0.975)``` in the QR-GWAS R script.



### Step 2

QR-PRS construction with Pruning and Thresholding (P+T). PLINK2 is requried (https://www.cog-genomics.org/plink/2.0/).

```bash
# QR-PRS using pruned variants with quantile-specific p-value < 0.05
file_genotype_train="genotype.train" # PLINK binary genotype fileset prefix for the training set
file_genotype_test="genotype.test" # PLINK binary genotype fileset prefix for the test set
file_sumstat="sumstat.train.tsv.gz" # QR-GWAS summary statistics file computed using the training set
tau_pool=( 0.025 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.975 ) # quantile level
bidx_pool=( $(seq 23 33) ) # column index of quantile-specific QR-GWAS beta in the summary statistics file
for idx in "${!tau_pool[@]}"; do
  tau="${tau_pool[idx]}"
  bidx="${bidx_pool[idx]}"
  colname_id="ID" # column name of variant ID from the summary statistics file
  colname_pvalue="P_QR${tau}" # column name of quantile-specific QR-GWAS p-value (quantile level = ${tau}) in the summary statistics file
  file_variant="variant.pruned.tau${tau}" # list of pruned variants
  file_prs_train="example.train.qrprs.tau${tau}" # quantile-specific QR-PRS for the training set
  file_prs_test="example.test.qrprs.tau${tau}" # quantile-specific QR-PRS for the test set
  plink2 \
  --bfile ${file_genotype_train} \
  --clump ${file_sumstat} \
  --clump-id-field ${colname_id} \
  --clump-p-field ${colname_pvalue} \
  --clump-p1 0.05 \
  --clump-r2 0.1 \
  --clump-kb 250 \
  --out ${file_variant}
  awk 'NR > 1 {print $3}' ${file_variant}.clumps > ${file_variant}.snplist
  plink2 \
  --bfile ${file_genotype_train} \
  --extract ${file_variant}.snplist \
  --score ${file_sumstat} 3 5 ${bidx} header cols=scoresums \
  --out ${file_prs_train}
  plink2 \
  --bfile ${file_genotype_test} \
  --extract ${file_variant}.snplist \
  --score ${file_sumstat} 3 5 ${bidx} header cols=scoresums \
  --out ${file_prs_test}
done
```

The expected output includes the quantile-specific QR-PRS in PLINK2 sscore files.



### Step 3

Compute plug-in and conformal prediction intervals.

We provide a toy dataset in [example](/example). 
The input data includes 
1) the phenotype tables for [training set](example/example.train.phenotype.tsv) and [test set](example/example.test.phenotype.tsv). 
2) QR-PRS scores produced by PLINK2 from Step 2 for [training set](example/example.train.qrprs.tau0.025.sscore) and [test set](example/example.test.qrprs.tau0.025.sscore). 

Run the R script [example.qrprs.R](example.qrprs.R) to compute plug-in and conformal prediction intervals for the test set in [example](/example), perform risk stratification, and evaluate the empirical coverage. 

The table below shows the number of individuals in each category based on risk stratification using QR-PRS on the test set.
| Category          |   n |
|-------------------|-----|
| Certain-above t   |   0 |
| Uncertain-above t |  35 |
| Uncertain-below t | 183 |
| Certain-below t   | 128 |

The R script also calculate the empirical coverage of plug-in and conformal prediction intervals as follows.

Coverage of QR-PRS prediction intervals: 72.8%.

Coverage of CQR-PRS prediction intervals: 97.1%.


