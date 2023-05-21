# AIM: TEST DRUGS->MEDIATORS
#####################
# Install and load libraries #
#####################
library(tidyverse)
library(ieugwasr)
library(TwoSampleMR)
library(tibble)
########################################################
# MR: Age at Menarche
########################################################
# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1095", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./Age_at_menarche/drugs_age_at_menarche_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./Age_at_menarche/mr_drugs_age_at_menarche_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./Age_at_menarche/drug_age_at_menarche_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./Age_at_menarche/drug_age_at_menarche_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./Age_at_menarche/mr_drugs_age_at_menarche_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./Age_at_menarche/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./Age_at_menarche/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./Age_at_menarche/pleiotropy.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
#APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("Age at menarche || id:ieu-a-1095")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./Age_at_menarche/mr_presso_global.csv")

########################################################
# MR: Age at menopause
########################################################
rm(list = ls())

# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-b-4824", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./Age_at_menopause/drugs_age_at_menopause_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./Age_at_menopause/mr_drugs_age_at_menopause_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./Age_at_menopause/drug_age_at_menopause_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./Age_at_menopause/drug_age_at_menopause_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./Age_at_menopause/mr_drugs_age_at_menopause_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./Age_at_menopause/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./Age_at_menopause/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./Age_at_menopause/pleiotropy.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("Age at menopause || id:ieu-a-1095")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./Age_at_menopause/mr_presso_global.csv")

########################################################
# MR: Body Mass Index
########################################################
rm(list = ls())

# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-974", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./BMI/drugs_BMI_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./BMI/mr_drugs_BMI_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./BMI/drug_BMI_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./BMI/drug_BMI_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./BMI/mr_drugs_BMI_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./BMI/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./BMI/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./BMI/pleiotropy.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("Body Mass Index || id:ieu-a-974")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./BMI/mr_presso_global.csv")

########################################################
# MR: Estradiol
########################################################
rm(list = ls())

# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ebi-a-GCST90020092", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./estradiol/drugs_estradiol_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./estradiol/mr_drugs_estradiol_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./estradiol/drug_estradiol_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./estradiol/drug_estradiol_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./estradiol/mr_drugs_estradiol_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./estradiol/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./estradiol/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./estradiol/pleiotropy.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("Estradiol levels || id:ebi-a-GCST90020092")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./estradiol/mr_presso_global.csv")

########################################################
# MR: ER expression
########################################################
rm(list = ls())

# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "prot-a-991", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./ER/drugs_ApoA_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./ER/mr_drugs_ApoA_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./ER/drug_ApoA_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./ER/drug_ApoA_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./ER/mr_drugs_ApoA_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./ER/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./ER/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./ER/pleiotropy.csv")

# 10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("Estrogen receptor || id:prot-a-991")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./ER/mr_presso_global.csv")

########################################################
# MR:  Apolipoprotein A
########################################################
rm(list = ls())

# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "met-c-842", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./ApoA/drugs_ApoA_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./ApoA/mr_drugs_ApoA_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./ApoA/drug_ApoA_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./ApoA/drug_ApoA_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./ApoA/mr_drugs_ApoA_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./ApoA/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./ApoA/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./ApoA/pleiotropy.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("Apolipoprotein A-I || id:met-c-842")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./ApoA/mr_presso_global.csv")

########################################################
# MR:  Apolipoprotein B
########################################################
rm(list = ls())

# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "met-c-843", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./ApoB/drugs_ApoB_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./ApoB/mr_drugs_ApoB_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./ApoB/drug_ApoB_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./ApoB/drug_ApoB_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./ApoB/mr_drugs_ApoB_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./ApoB/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./ApoB/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./ApoB/pleiotropy.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("Apolipoprotein B || id:met-c-843")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./ApoB/mr_presso_global.csv")

########################################################
# MR:  HDL
########################################################
rm(list = ls())

# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-299", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./HDL/drugs_HDL_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./HDL/mr_drugs_HDL_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./HDL/drug_HDL_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./HDL/drug_HDL_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./HDL/mr_drugs_HDL_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./HDL/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./HDL/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./HDL/pleiotropy.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("HDL cholesterol || id:ieu-a-299")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./HDL/mr_presso_global.csv")

########################################################
# MR:  LDL
########################################################
rm(list = ls())

# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-300", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./LDL/drugs_LDL_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./LDL/mr_drugs_LDL_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./LDL/drug_LDL_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./LDL/drug_LDL_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./LDL/mr_drugs_LDL_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./LDL/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./LDL/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./LDL/pleiotropy.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("LDL cholesterol || id:ieu-a-300")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./LDL/mr_presso_global.csv")

########################################################
# MR:  TC
########################################################
rm(list = ls())

# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-301", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./TC/drugs_TC_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./TC/mr_drugs_TC_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./TC/drug_TC_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./TC/drug_TC_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./TC/mr_drugs_TC_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./TC/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./TC/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./TC/pleiotropy.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("Total cholesterol || id:ieu-a-301")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./TC/mr_presso_global.csv")

########################################################
# MR:  TG
########################################################
rm(list = ls())

# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-302", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./TG/drugs_TG_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./TG/mr_drugs_TG_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./TG/drug_TG_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./TG/drug_TG_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./TG/mr_drugs_TG_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./TG/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./TG/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./TG/pleiotropy.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("Triglycerides || id:ieu-a-302")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./TG/mr_presso_global.csv")
########################################################
# MR:  WHR
########################################################
rm(list = ls())

# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA: Waist-to-hip ratio
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-73", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./WHR/drugs_WHR_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
#i1 <- mr_results$b > 0
#mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./WHR/mr_drugs_WHR_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./WHR/drug_WHR_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./WHR/drug_WHR_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./WHR/mr_drugs_WHR_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./WHR/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./WHR/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./WHR/pleiotropy.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("Waist-to-hip ratio || id:ieu-a-73")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./WHR/mr_presso_global.csv")
########################################################
# MR:  WC
########################################################
rm(list = ls())

# Set working directory
setwd("~/Desktop/Skeleton_projects/oxysterol/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA: Waist-to-hip ratio
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-61", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
setwd("/Users/hugopomaresmillan/Desktop/Skeleton_projects/oxysterol/analysis/drug_mediators")
mr_results <- mr(single_dat)
mr_results
# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, single_dat)
pdf("./WC/drugs_WC_scatter.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

#5. Estimate odds ratio and 95% confidence interval
# Flip the results (to protective effect)
#i1 <- mr_results$b > 0
#mr_results$b[i1] <- -1*(mr_results$b[i1])
##
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/7 = 0.007
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.007, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)
# Export results
write.csv(results,"./WC/mr_drugs_WC_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(single_dat)
mr_forest <- mr_forest_plot(res_single)

pdf("./WC/drug_WC_forest.pdf",  width = 15, height = 20)
mr_forest
dev.off()

#Leave one out analysis
res_loo <- mr_leaveoneout(single_dat)
mr_loo <- mr_leaveoneout_plot(res_loo)
pdf("./WC/drug_WC_leaveoneout.pdf",  width = 15, height = 20)
mr_loo
dev.off()

write.csv(res_loo,"./WC/mr_drugs_WC_res_loo.csv")

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
# single_dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
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
# write.csv(I2, file="./WC/isq_un_and_weighted.csv", row.names = FALSE)

#9. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(single_dat)
write.csv(mr_heterogeneity(single_dat),"./WC/heterogeneity.csv")

#Pleiotropy test
mr_pleiotropy_test(single_dat)
write.csv(mr_pleiotropy_test(single_dat),"./WC/pleiotropy.csv")

#10. MR PRESSO
library(MRPRESSO)
# Run MR-PRESSO global method
HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="HMGCR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="NPC1L1"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="PCSK9"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LDLR"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="CETP"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="LPL"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = single_dat[single_dat$exposure=="APOC3"&single_dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=7))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c("Waist circumference || id:ieu-a-61")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- NPC1L1_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- PCSK9_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- LDLR_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- CETP_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- LPL_presso$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- APOC3_presso$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "./WC/mr_presso_global.csv")
#################END###################################