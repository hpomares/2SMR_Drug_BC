# AIM: TEST MEDIATORS->CANCER
#####################
# Install and load libraries #
#####################
library(tidyverse)
library(ieugwasr)
library(TwoSampleMR)
library(tibble)
########################################################
# MR: Age_at_Menarche
########################################################
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/Age_at_menarche")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for age_menarche
exposure_data <- extract_instruments("ieu-a-1095", p1 = 5e-08)

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("age_menarche_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_age_menarche_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_age_menarche_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_age_menarche_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_age_menarche_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_age_menarche.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_age_menarche.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
save.image("age_menarche.RData")

########################################################
# MR: Age_at_menopause
########################################################
rm(list = ls())
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/Age_at_menopause")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for age_menopause
exposure_data <- extract_instruments("ieu-b-4824", p1 = 5e-08)

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("age_menopause_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_age_menopause_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_age_menopause_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_age_menopause_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_age_menopause_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_age_menopause.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_age_menopause.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
save.image("age_menopause.RData")

########################################################
# MR: ApoA
########################################################
rm(list = ls())
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/ApoA")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for ApoA
exposure_data <- extract_instruments("met-c-842", p1 = 5e-08)

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("ApoA_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_ApoA_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_ApoA_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_ApoA_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_ApoA_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_ApoA.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_ApoA.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
save.image("ApoA.RData")

########################################################
# MR: ApoB
########################################################
rm(list = ls())
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/ApoB")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for ApoB
exposure_data <- extract_instruments("met-c-843", p1 = 5e-08)

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("ApoB_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_ApoB_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_ApoB_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_ApoB_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_ApoB_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_ApoB.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_ApoB.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
save.image("ApoB.RData")

########################################################
# MR: BMI
########################################################
rm(list = ls())
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/BMI")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for BMI
exposure_data <- extract_instruments("ieu-a-974", p1 = 5e-08)

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("BMI_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_BMI_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_BMI_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_BMI_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_BMI_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_BMI.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_BMI.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
 
write.csv(mr_presso_global, "mr_presso_global.csv")

########################################################
# MR: ER
########################################################
rm(list = ls())
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/ER")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for ER
exposure_data <- extract_instruments("prot-a-991", p1 = 5e-08)

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("ER_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_ER_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_ER_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_ER_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_ER_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_ER.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_ER.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
save.image("ER.RData")

########################################################
# MR: estradiol
########################################################
rm(list = ls())
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/estradiol")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for estradiol
exposure_data <- extract_instruments("ebi-a-GCST90020092", p1 = 5e-08)

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("estradiol_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_estradiol_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_estradiol_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_estradiol_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_estradiol_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_estradiol.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_estradiol.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
save.image("estradiol.RData")

########################################################
# MR: HDL
########################################################
rm(list = ls())
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/HDL")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for HDL
exposure_data <- extract_instruments("ieu-a-299", p1 = 5e-08)

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("HDL_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_HDL_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_HDL_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_HDL_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_HDL_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_HDL.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_HDL.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
save.image("HDL.RData")

########################################################
# MR: LDL
########################################################
rm(list = ls())
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/LDL")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for LDL
exposure_data <- extract_instruments("ieu-a-300", p1 = 5e-08)

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("LDL_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_LDL_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_LDL_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_LDL_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_LDL_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_LDL.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_LDL.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
save.image("LDL.RData")

########################################################
# MR: TC
########################################################
rm(list = ls())
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/TC")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for TC
exposure_data <- extract_instruments("ieu-a-301", p1 = 5e-08)

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("TC_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_TC_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_TC_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_TC_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_TC_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_TC.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_TC.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
save.image("TC.RData")

########################################################
# MR: TG
########################################################
rm(list = ls())
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/TG")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for TG
exposure_data <- extract_instruments("ieu-a-302", p1 = 5e-08)

#2. SNP-OUTGOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("TG_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_TG_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_TG_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_TG_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_TG_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_TG.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_TG.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
save.image("TG.RData")

########################################################
# MR: WC
########################################################
rm(list = ls())
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/WC")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for WC
exposure_data <- extract_instruments("ieu-a-61", p1 = 5e-08)

#2. SNP-OUWCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("WC_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_threshold", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_WC_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_WC_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_WC_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_WC_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_WC.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_WC.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
save.image("WC.RData")

########################################################
# MR: WHR
########################################################
rm(list = ls())
# Set working directory
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/mediators_cancer/WHR")

#1. SNP-EXPOSURE SUMMARY DATA
# SNP Extraction from MRBase for WHR
exposure_data <- extract_instruments("ieu-a-73", p1 = 5e-08)

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data( snps = exposure_data$SNP,
                                      outcomes = c(
                                        "ieu-a-1126","ieu-a-1127","ieu-a-1128"#,"ieu-a-1165","ieu-b-4810"
                                      ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
write.table(out_dat_comb$SNP, "snplist_allSNPs.txt", row.names=F, col.names=F)
write.table(out_dat_comb, "exposure_data.txt", row.names=F, col.names=F)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(single_dat)
mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("WHR_cancer_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/3 = 0.016
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.016, "pass_threshold", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"mr_WHR_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)

# plot 1 to many
res<-split_outcome(mr_results) # to keep the Y axis label clean we exclude the exposure ID labels from the exposure column 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
# separate for plot
res<-subset_on_method(res) #default is to subset on either the IVW method (>1 instrumental SNP) or Wald ratio method (1 instrumental SNP). 
res<-sort_1_to_many(res,b="b",sort_action=4) #this sorts results by decreasing effect size (largest effect at top of the plot)
plot1<-forest_plot_1_to_many(res,b="b",se="se",
                             exponentiate=T,trans="log2",ao_slc=F,lo=0.004,
                             up=461,col1_width=2,TraitM="outcome",
                             col_text_size=3,xlab="")

pdf("exp_WHR_forest.pdf",  width = 15, height = 20)
plot1
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("exp_WHR_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"mr_WHR_cancer_res_loo.csv")

# #7. Calculate I-squared, r-squared and F-statistics
# #I-squared function
# Isq <- function(y,s){
#   k          = length(y)
#   w          = 1/s^2; sum.w  = sum(w)
#   mu.hat     = sum(y*w)/sum.w  
#   Q          = sum(w*(y-mu.hat)^2)
#   Isq        = (Q - (k-1))/Q
#   Isq        = max(0,Isq)
#   return(Isq)
# }
# #Calculate Isq wieghted and unweighted
# I2<-c()
# single_dat <- harmonise_data(exposure_data, out_dat_comb, action = 1)
# str(single_dat)
# 
# #F-statistic
# single_dat$samplesize.exposure <- 188578
# single_dat$samplesize.outcome <- 3301
# single_dat <- steiger_filtering(single_dat) 
# 
# N = single_dat$samplesize.exposure[1] #sample size
# K = length(single_dat$SNP) #number of SNPs
# total_r2 <- sum(single_dat$rsq.exposure) 
# Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
# total_r2
# Fstat
# 
# #Rename required columns
# single_dat $BetaXG<-single_dat $beta.exposure
# single_dat $seBetaXG<-single_dat $se.exposure
# BetaXG   = single_dat $BetaXG
# seBetaXG = single_dat $seBetaXG 
# seBetaYG<-single_dat $se.outcome
# 
# BXG = abs(BetaXG)         
# # Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
# F = BXG^2/seBetaXG^2
# mF = mean(F)
# Isq_unweighted <- Isq(BXG,seBetaXG) #unweighted
# Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted
# 
# #Save results
# output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
# I2<-rbind(I2, output)
# colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
# write.csv(I2, file="isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"heterogeneity_exp_WHR.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"pleiotropy_exp_WHR.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
BC_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
BCERN_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
save.image("WHR.RData")

############END########################