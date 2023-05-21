# AIM: TEST CANCER->MEDIATORS
#####################
# Install and load libraries #
#####################
library(tidyverse)
library(ieugwasr)
library(TwoSampleMR)
library(tibble)

########################################################
# MR
########################################################
setwd("~/analysis/cancer_mediators")

#1. SELECT SNP-EXPOSURE SUMMARY DATA#
# SNP Extraction from MRBase
ao <- available_outcomes()
exposure_dat <- extract_instruments(c("ieu-a-1126","ieu-a-1127","ieu-a-1128"))
#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_dat$SNP,
                                      outcomes =c("ieu-a-1095", "ieu-b-4824", "ieu-a-974", "ebi-a-GCST90020092", "prot-a-991", "met-c-842", "met-c-843", "ieu-a-299", "ieu-a-300","ieu-a-301","ieu-a-302","ieu-a-73","ieu-a-61")
 , proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
dat <- harmonise_data(exposure_dat, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, dat)
pdf("cancer_mediators_scatter_plot.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# 0.05/39 = 0.001
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.001, "pass_threshold", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_mediators_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(dat)
res_single_out <- split( res_single, f = res_single$id.exposure )

# Leave one out analysis
res_loo <- mr_leaveoneout(dat)

# ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(dat)
write.csv(mr_heterogeneity(dat),"heterogeneity_results.csv")

#Pleiotropy test
mr_pleiotropy_test(dat)
write.csv(mr_pleiotropy_test(dat),"pleiotropy_results.csv")

# MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
#BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
#EP_BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
EN_BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

################################################
save.image("cancer_mediators_analysis.RData")
############END########################
