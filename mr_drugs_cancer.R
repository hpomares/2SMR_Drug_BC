# AIM: TEST DRUGS->CANCERS
#####################
# Install and load libraries #
#####################
# install.packages("devtools")
# devtools::install_github("MRCIEU/MRInstruments", force = TRUE)
# devtools::install_github('MRCIEU/TwoSampleMR')
# install.packages("MendelianRandomization", force = TRUE)
# install.packages("LDlinkR")
# install.packages("plyr")
# install.packages("ggplot2")
# install.packages("ggpubr")
# install.packages("simex")
# devtools::install_github("rondolab/MR-PRESSO")

# sample size calculator for MR: https://shiny.cnsgenomics.com/mRnd/

library(devtools)
library(TwoSampleMR)
library(MRInstruments)
library(dplyr)
library(plyr) 
library(ggplot2)
library(MendelianRandomization)
library(gridExtra)
library(grid)
library(lattice)
library(LDlinkR)
library(ggpubr)
library(simex)
library(MRPRESSO)
library(DiagrammeR)
library(ieugwasr)

# Set working directory
setwd("~/analysis/")

# check ID in repository
ao <- available_outcomes()
ao1 <- subset(ao, id=="prot-a-2087"|id=="ieu-a-1126"|id=="ieu-a-1127"|id=="ieu-a-1128"|id=="ieu-a-1165"
              |id=="prot-a-991"|id=="ieu-a-1095"|id=="ieu-a-974"|id=="ieu-a-299"|id=="ieu-a-300" | id=="ebi-a-GCST90020092"
              |id=="ieu-a-301"|id=="ieu-a-302"|id=="met-c-842"|id=="met-c-843")
library(rio);export(ao1,"~/analysis/gwas_charac.xlsx")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12

#2. SNP-OUTCOME SUMMARY DATA
out_dat_comb <- extract_outcome_data(snps = drug_exp_data$SNP,
                                     outcomes = c(
                                       "ieu-a-1126","ieu-a-1127","ieu-a-1128" 
                                       ), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
table(out_dat_comb$outcome)
library(rio);export(out_dat_comb,"~/analysis/drug_cancer/out_dat_comb.csv")
write.table(out_dat_comb$SNP, "~/analysis/drug_cancer/snplist_allSNPs.txt", row.names=F, col.names=T)

#3. HARMONIZE DATASETS                                 
# Harmonise the datasets so that the effect alleles are the same (and reflect the increasing allele) 
dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 2)

#4. ESTIMATE THE CAUSAL EFFECTS ON THE OUTCOME
mr_results <- mr(dat); mr_results

# Run scatter plot code before flipping 
mr_scatter <- mr_scatter_plot(mr_results, dat)
#### EDIT PLOTS####
plot_HMGCR_BC<-mr_scatter$`lPaV5D.ieu-a-1126` + theme_bw()
tiff("plot_HMGCR_BC.tiff", units="in", width=5, height=5, res=300)
plot_HMGCR_BC + labs(y= "SNP effect on Breast cancer", x = "SNP efefect on HMGCR gene")
dev.off()
#
plot_HMGCR_EP_BC<-mr_scatter$`lPaV5D.ieu-a-1127` + theme_bw()
tiff("plot_HMGCR_EPBC.tiff", units="in", width=5, height=5, res=300)
plot_HMGCR_EP_BC + labs(y= "SNP efefect on estrogen receptor-positive Breast cancer", x = "SNP efefect on HMGCR gene")
dev.off()
#
plot_HMGCR_EN_BC<-mr_scatter$`lPaV5D.ieu-a-1128` + theme_bw()
tiff("plot_HMGCR_ENBC.tiff", units="in", width=5, height=5, res=300)
plot_HMGCR_EN_BC + labs(y= "SNP efefect on estrogen receptor-negative Breast cancer", x = "SNP efefect on HMGCR gene")
dev.off()

# all plots
pdf("~/analysis/drug_cancer/drugs_cancer_scatter_plot.pdf",  width = 15, height = 20)
print(mr_scatter)
dev.off()

ggarrange(mr_scatter$`7uELty.ieu-a-1126`,mr_scatter$`7uELty.ieu-a-1127`,mr_scatter$`7uELty.ieu-a-1128`, mr_scatter$`7uELty.ieu-a-1165`)

#5. ESTIMATE AND PLOT THE CAUSAL EFFECTS OF THE TRAIT ON THE OUTCOME
# Flip the results (to protective effect of LDL-lowering)
i1 <- mr_results$b > 0
mr_results$b[i1] <- -1*(mr_results$b[i1])
# Estimate odds ratio and 95% confidence interval
mr_results$or <- exp(mr_results$b)
mr_results$cil <- exp(mr_results$b-1.96*mr_results$se)
mr_results$ciu <- exp(mr_results$b+1.96*mr_results$se)
# Pvalue corrected (Bonferroni)
# 0.05/21 = 0.0023
require(tidyverse);mr_results<- mr_results %>% mutate(pval_corrected = if_else(pval < 0.0023, "pass_correction", "Not_pass"))
results<-cbind.data.frame(mr_results$outcome,mr_results$exposure,mr_results$nsnp,mr_results$method,mr_results$b,mr_results$se,mr_results$pval,mr_results$or,mr_results$cil,mr_results$ciu, mr_results$pval_corrected)

# Export results
write.csv(results,"~/analysis/drug_cancer/mr_drugs_cancer_results.csv")

# 6. PLOT FOREST AND LEAVE-ONE-OUT
#Plot causal estimates for effect of lipid-lowering drugs on head and neck cancer by individual SNP Forest Plot
res_single <- mr_singlesnp(dat)
res_single_out <- split( res_single, f = res_single$id.exposure )

# Leave one out analysis
res_loo <- mr_leaveoneout(dat)
write.csv(res_loo,"~/analysis/drug_cancer/mr_drugs_cancer_res_loo.csv")

# 7. Calculate I-squared, r-squared and F-statistics
# I-squared function
Isq <- function(y,s){
  k          = length(y)
  w          = 1/s^2; sum.w  = sum(w)
  mu.hat     = sum(y*w)/sum.w  
  Q          = sum(w*(y-mu.hat)^2)
  Isq        = (Q - (k-1))/Q
  Isq        = max(0,Isq)
  return(Isq)
}

# Calculate Isq weighted and unweighted
I2<-c()
dat <- harmonise_data(drug_exp_data, out_dat_comb, action = 1)
str(dat)

# F-statistic 
dat$samplesize.exposure <- 188577 # from the Global Lipids Genetic Consortium (GLGC)  up to 188,577 individuals of European ancestry
dat$samplesize.outcome <- 228951
dat <- steiger_filtering(dat)
N = dat$samplesize.exposure[1] #sample size
K = length(dat$SNP) #number of SNPs
total_r2 <- sum(dat$rsq.exposure)
Fstat <- (N-K-1)/K * total_r2 / (1-total_r2)
total_r2
Fstat

#Rename required columns
dat$BetaXG<-dat $beta.exposure
dat$seBetaXG<-dat $se.exposure
BetaXG   = dat $BetaXG
seBetaXG = dat $seBetaXG 
seBetaYG<-dat $se.outcome

BXG = abs(BetaXG)         

# Calculate F statistics and I-squared statistics to measure Instrument strength for MR-Egger
F = BXG^2/seBetaXG^2
mF = mean(F)
Isq_unweighted <- Isq(BXG,seBetaXG) # unweighted
Isq_weighted <- Isq((BXG/seBetaYG),(seBetaXG/seBetaYG)) #weighted

#Save results
output<-cbind (F, mF, Isq_unweighted, Isq_weighted)
I2<-rbind(I2, output)
colnames(I2) <- c("F", "mF", "Isq_unweighted", "Isq_weighted")
write.csv(I2, file="~/analysis/drug_cancer/regression_dilution_isq_weighted.csv", row.names = FALSE)

#8. ASSESS HETEROGENEITY AND PLEIOTROPY 
#Heterogeneity test
mr_heterogeneity(dat)
write.csv(mr_heterogeneity(dat),"~/analysis/drug_cancer/heterogeneity_results.csv")

#Pleiotropy test
mr_pleiotropy_test(dat)
write.csv(mr_pleiotropy_test(dat),"~/analysis/drug_cancer/pleiotropy_results.csv")

# MR PRESSO
# Run MR-PRESSO global method
# HMGCR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="HMGCR"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
HMGCR_presso1 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="HMGCR"&dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
HMGCR_presso2 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="HMGCR"&dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
HMGCR_presso3 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="HMGCR"&dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
HMGCR_presso4 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="HMGCR"&dat$outcome=="Breast cancer (Survival) || id:ieu-a-1165"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

# NPC1L1_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="NPC1L1"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso1 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="NPC1L1"&dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso2 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="NPC1L1"&dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso3 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="NPC1L1"&dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
NPC1L1_presso4 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="NPC1L1"&dat$outcome=="Breast cancer (Survival) || id:ieu-a-1165"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

# PCSK9_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="PCSK9"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso1 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="PCSK9"&dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso2 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="PCSK9"&dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso3 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="PCSK9"&dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
PCSK9_presso4 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="PCSK9"&dat$outcome=="Breast cancer (Survival) || id:ieu-a-1165"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

# LDLR_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="LDLR"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso1 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="LDLR"&dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso2 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="LDLR"&dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso3 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="LDLR"&dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LDLR_presso4 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="LDLR"&dat$outcome=="Breast cancer (Survival) || id:ieu-a-1165"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

# CETP_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="CETP"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso1 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="CETP"&dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso2 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="CETP"&dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso3 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="CETP"&dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
CETP_presso4 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="CETP"&dat$outcome=="Breast cancer (Survival) || id:ieu-a-1165"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

# LPL_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="LPL"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso1 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="LPL"&dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso2 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="LPL"&dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso3 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="LPL"&dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
LPL_presso4 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="LPL"&dat$outcome=="Breast cancer (Survival) || id:ieu-a-1165"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

# APOC3_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="APOC3"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso1 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="APOC3"&dat$outcome=="Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso2 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="APOC3"&dat$outcome=="ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso3 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="APOC3"&dat$outcome=="ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)
APOC3_presso4 <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat[dat$exposure=="APOC3"&dat$outcome=="Breast cancer (Survival) || id:ieu-a-1165"&dat$mr_keep=="TRUE",], NbDistribution = 1000,  SignifThreshold = 0.05)

headers<-c("outcome", "exposure", "RSSobs","Pvalue")
mr_presso_global <- as.data.frame(matrix(,ncol=4,nrow=28))
names(mr_presso_global)<-headers  
mr_presso_global$outcome <- c(
                              "Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1126", 
                              "ER+ Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1127", 
                              "ER- Breast cancer (Combined Oncoarray; iCOGS; GWAS meta analysis) || id:ieu-a-1128", 
                              "Breast cancer (Survival) || id:ieu-a-1165")
mr_presso_global$exposure <- c("HMGCR", "NPC1L1", "PCSK9", "LDLR", "CETP", "LPL", "APOC3")

mr_presso_global$RSSobs[1] <- HMGCR_presso1$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[1] <- HMGCR_presso1$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[2] <- HMGCR_presso2$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[2] <- HMGCR_presso2$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[3] <- HMGCR_presso3$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[3] <- HMGCR_presso3$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[4] <- HMGCR_presso4$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[4] <- HMGCR_presso4$`MR-PRESSO results`$`Global Test`$Pvalue

mr_presso_global$RSSobs[5] <- NPC1L1_presso1$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[5] <- NPC1L1_presso1$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[6] <- NPC1L1_presso2$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[6] <- NPC1L1_presso2$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[7] <- NPC1L1_presso3$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[7] <- NPC1L1_presso3$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[8] <- NPC1L1_presso4$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[8] <- NPC1L1_presso4$`MR-PRESSO results`$`Global Test`$Pvalue

mr_presso_global$RSSobs[9] <- PCSK9_presso1$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[9] <- PCSK9_presso1$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[10] <- PCSK9_presso2$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[10] <- PCSK9_presso2$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[11] <- PCSK9_presso3$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[11] <- PCSK9_presso3$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[12] <- PCSK9_presso4$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[12] <- PCSK9_presso4$`MR-PRESSO results`$`Global Test`$Pvalue

mr_presso_global$RSSobs[13] <- LDLR_presso1$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[13] <- LDLR_presso1$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[14] <- LDLR_presso2$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[14] <- LDLR_presso2$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[15] <- LDLR_presso3$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[15] <- LDLR_presso3$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[16] <- LDLR_presso4$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[16] <- LDLR_presso4$`MR-PRESSO results`$`Global Test`$Pvalue

mr_presso_global$RSSobs[17] <- CETP_presso1$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[17] <- CETP_presso1$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[18] <- CETP_presso2$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[18] <- CETP_presso2$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[19] <- CETP_presso3$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[19] <- CETP_presso3$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[20] <- CETP_presso4$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[20] <- CETP_presso4$`MR-PRESSO results`$`Global Test`$Pvalue

mr_presso_global$RSSobs[21] <- LPL_presso1$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[21] <- LPL_presso1$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[22] <- LPL_presso2$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[22] <- LPL_presso2$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[23] <- LPL_presso3$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[23] <- LPL_presso3$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[24] <- LPL_presso4$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[24] <- LPL_presso4$`MR-PRESSO results`$`Global Test`$Pvalue

mr_presso_global$RSSobs[25] <- APOC3_presso1$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[25] <- APOC3_presso1$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[26] <- APOC3_presso2$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[26] <- APOC3_presso2$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[27] <- APOC3_presso3$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[27] <- APOC3_presso3$`MR-PRESSO results`$`Global Test`$Pvalue
mr_presso_global$RSSobs[28] <- APOC3_presso4$`MR-PRESSO results`$`Global Test`$RSSobs
mr_presso_global$Pvalue[28] <- APOC3_presso4$`MR-PRESSO results`$`Global Test`$Pvalue

write.csv(mr_presso_global, "~/analysis/drug_cancer/mr_presso_global.csv")

#9. ACCOUNT FOR LD STRUCTURE
setwd("~/analysis/drug_cancer/")
dat <- dat[dat$mr_keep=="TRUE",]
HMGCR <- dat[dat$exposure=="HMGCR",]
LDmatrix(HMGCR$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "HMGCR_cor.txt")
NPC1L1 <- dat[dat$exposure=="NPC1L1",]
LDmatrix(NPC1L1$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "NPC1L1_cor.txt")
PCSK9 <- dat[dat$exposure=="PCSK9",]
LDmatrix(PCSK9$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "PCSK9_cor.txt")
LDLR <- dat[dat$exposure=="LDLR",]
LDmatrix(LDLR$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "LDLR_cor.txt")
CETP <- dat[dat$exposure=="CETP",]
LDmatrix(CETP$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "CETP_cor.txt")
APOC3 <- dat[dat$exposure=="APOC3",]
LDmatrix(APOC3$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "APOC3_cor.txt")
LPL <- dat[dat$exposure=="LPL",]
LDmatrix(LPL$SNP, pop="CEU", r2d = "r2", token = "b0a2f3f3c1f4", file = "LPL_cor.txt")

# correl <- read.table("APOC3_cor.txt", header=T)
# correl <- correl[order(correl$RS_number),]
# correl <- correl[,order(names(correl))]
# APOC3 <- APOC3[APOC3$SNP %in% names(correl),]
# correl <- within(correl, rm(RS_number))
# colnames(correl) <- NULL
# rownames(correl) <- NULL
# correl <- data.matrix(correl)
# mr_ivw_APOC3 <- mr_ivw(mr_input(bx = APOC3$beta.exposure, bxse = APOC3$se.exposure, by = APOC3$beta.outcome, byse = APOC3$se.outcome, correlation = correl))
# mr_egger_APOC3 <- mr_egger(mr_input(bx = APOC3$beta.exposure, bxse = APOC3$se.exposure, by = APOC3$beta.outcome, byse = APOC3$se.outcome, correlation = correl))
# mr_median_APOC3 <- mr_median(mr_input(bx = APOC3$beta.exposure, bxse = APOC3$se.exposure, by = APOC3$beta.outcome, byse = APOC3$se.outcome, correlation = correl))
# 
# Plot
library(corrplot)
correl <- read.table("~/analysis/drug_cancer/APOC3_cor.txt", header=T)
rownames(correl) <- correl[,1]
correl[,1] <- NULL

pdf(file = "~/analysis/drug_cancer/APOC3_cor.pdf")
corrplot(as.matrix(correl), tl.col = "brown", tl.srt = 30, bg = "White", method = "color",
         title = "\n\n Correlation plot of APOC3 SNPs",
         type = "lower")
dev.off()

# correl <- read.table("LPL_cor.txt", header=T)
# correl <- correl[order(correl$RS_number),]
# correl <- correl[,order(names(correl))]
# LPL <- LPL[LPL$SNP %in% names(correl),]
# correl <- within(correl, rm(RS_number))
# colnames(correl) <- NULL
# rownames(correl) <- NULL
# correl <- data.matrix(correl)
# mr_ivw_LPL <- mr_ivw(mr_input(bx = LPL$beta.exposure, bxse = LPL$se.exposure, by = LPL$beta.outcome, byse = LPL$se.outcome, correlation = correl))
# mr_egger_LPL <- mr_egger(mr_input(bx = LPL$beta.exposure, bxse = LPL$se.exposure, by = LPL$beta.outcome, byse = LPL$se.outcome, correlation = correl))
# mr_median_LPL <- mr_median(mr_input(bx = LPL$beta.exposure, bxse = LPL$se.exposure, by = LPL$beta.outcome, byse = LPL$se.outcome, correlation = correl))

# Plot
library(corrplot)
correl <- read.table("~/analysis/drug_cancer/LPL_cor.txt", header=T)
rownames(correl) <- correl[,1]
correl[,1] <- NULL
pdf(file = "~/analysis/drug_cancer/LPL_cor.pdf")
corrplot(as.matrix(correl), tl.col = "brown", tl.srt = 30, bg = "White", method = "color",
         title = "\n\n Correlation plot of LPL SNPs",
         type = "lower")
dev.off()

correl <- read.table("~/analysis/drug_cancer/HMGCR_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
HMGCR <- HMGCR[HMGCR$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_HMGCR <- mr_ivw(mr_input(bx = HMGCR$beta.exposure, bxse = HMGCR$se.exposure, by = HMGCR$beta.outcome, byse = HMGCR$se.outcome, correlation = correl))
mr_egger_HMGCR <- mr_egger(mr_input(bx = HMGCR$beta.exposure, bxse = HMGCR$se.exposure, by = HMGCR$beta.outcome, byse = HMGCR$se.outcome, correlation = correl))
mr_median_HMGCR <- mr_median(mr_input(bx = HMGCR$beta.exposure, bxse = HMGCR$se.exposure, by = HMGCR$beta.outcome, byse = HMGCR$se.outcome, correlation = correl))

# Plot
library(corrplot)
correl <- read.table("~/analysis/drug_cancer/HMGCR_cor.txt", header=T)
rownames(correl) <- correl[,1]
correl[,1] <- NULL
pdf(file = "~/analysis/drug_cancer/HMGCR_cor.pdf")
corrplot(as.matrix(correl), tl.col = "brown", tl.srt = 30, bg = "White", method = "color",
         title = "\n\n Correlation plot of HMGCR SNPs",
         type = "lower")
dev.off()

correl <- read.table("~/analysis/drug_cancer/NPC1L1_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
NPC1L1 <- NPC1L1[NPC1L1$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_NPC1L1 <- mr_ivw(mr_input(bx = NPC1L1$beta.exposure, bxse = NPC1L1$se.exposure, by = NPC1L1$beta.outcome, byse = NPC1L1$se.outcome, correlation = correl))
mr_egger_NPC1L1 <- mr_egger(mr_input(bx = NPC1L1$beta.exposure, bxse = NPC1L1$se.exposure, by = NPC1L1$beta.outcome, byse = NPC1L1$se.outcome, correlation = correl))
mr_median_NPC1L1 <- mr_median(mr_input(bx = NPC1L1$beta.exposure, bxse = NPC1L1$se.exposure, by = NPC1L1$beta.outcome, byse = NPC1L1$se.outcome, correlation = correl))

# Plot
correl <- read.table("~/analysis/drug_cancer/NPC1L1_cor.txt", header=T)
rownames(correl) <- correl[,1]
correl[,1] <- NULL
pdf(file = "~/analysis/drug_cancer/NPC1L1_cor.pdf")
corrplot(as.matrix(correl), tl.col = "brown", tl.srt = 30, bg = "White", method = "color",
         title = "\n\n Correlation plot of NPC1L1 SNPs",
         type = "lower")
dev.off()

correl <- read.table("~/analysis/drug_cancer/PCSK9_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
PCSK9 <- PCSK9[PCSK9$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_PCSK9 <- mr_ivw(mr_input(bx = PCSK9$beta.exposure, bxse = PCSK9$se.exposure, by = PCSK9$beta.outcome, byse = PCSK9$se.outcome, correlation = correl))
mr_egger_PCSK9 <- mr_egger(mr_input(bx = PCSK9$beta.exposure, bxse = PCSK9$se.exposure, by = PCSK9$beta.outcome, byse = PCSK9$se.outcome, correlation = correl))
mr_median_PCSK9 <- mr_median(mr_input(bx = PCSK9$beta.exposure, bxse = PCSK9$se.exposure, by = PCSK9$beta.outcome, byse = PCSK9$se.outcome, correlation = correl))

# Plot
correl <- read.table("~/analysis/drug_cancer/PCSK9_cor.txt", header=T)
rownames(correl) <- correl[,1]
correl[,1] <- NULL
pdf(file = "~/analysis/drug_cancer/PCSK9_cor.pdf")
corrplot(as.matrix(correl), tl.col = "brown", tl.srt = 30, bg = "White", method = "color",
         title = "\n\n Correlation plot of PCSK9 SNPs",
         type = "lower")
dev.off()

correl <- read.table("~/analysis/drug_cancer/LDLR_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
LDLR <- LDLR[LDLR$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_LDLR <- mr_ivw(mr_input(bx = LDLR$beta.exposure, bxse = LDLR$se.exposure, by = LDLR$beta.outcome, byse = LDLR$se.outcome, correlation = correl))
mr_egger_LDLR <- mr_egger(mr_input(bx = LDLR$beta.exposure, bxse = LDLR$se.exposure, by = LDLR$beta.outcome, byse = LDLR$se.outcome, correlation = correl))
mr_median_LDLR <- mr_median(mr_input(bx = LDLR$beta.exposure, bxse = LDLR$se.exposure, by = LDLR$beta.outcome, byse = LDLR$se.outcome, correlation = correl))

# Plot
correl <- read.table("~/analysis/drug_cancer/LDLR_cor.txt", header=T)
rownames(correl) <- correl[,1]
correl[,1] <- NULL
pdf(file = "~/analysis/drug_cancer/LDLR_cor.pdf")
corrplot(as.matrix(correl), tl.col = "brown", tl.srt = 30, bg = "White", method = "color",
         title = "\n\n Correlation plot of LDLR SNPs",
         type = "lower")
dev.off()

correl <- read.table("~/analysis/drug_cancer/CETP_cor.txt", header=T)
correl <- correl[order(correl$RS_number),]
correl <- correl[,order(names(correl))]
CETP <- CETP[CETP$SNP %in% names(correl),]
correl <- within(correl, rm(RS_number))
colnames(correl) <- NULL
rownames(correl) <- NULL
correl <- data.matrix(correl)
mr_ivw_CETP <- mr_ivw(mr_input(bx = CETP$beta.exposure, bxse = CETP$se.exposure, by = CETP$beta.outcome, byse = CETP$se.outcome, correlation = correl))
mr_egger_CETP<- mr_egger(mr_input(bx = CETP$beta.exposure, bxse = CETP$se.exposure, by = CETP$beta.outcome, byse = CETP$se.outcome, correlation = correl))
mr_median_CETP <- mr_median(mr_input(bx = CETP$beta.exposure, bxse = CETP$se.exposure, by = CETP$beta.outcome, byse = CETP$se.outcome, correlation = correl))

# Plot
correl <- read.table("~/analysis/drug_cancer/CETP_cor.txt", header=T)
rownames(correl) <- correl[,1]
correl[,1] <- NULL
pdf(file = "~/analysis/drug_cancer/CETP_cor.pdf")
corrplot(as.matrix(correl), tl.col = "brown", tl.srt = 30, bg = "White", method = "color",
         title = "\n\n Correlation plot of CETP SNPs",
         type = "lower")
dev.off()

# Combine results 
headers<-c("outcome","exposure","method","b","se","pval")
mr_results_corr <- as.data.frame(matrix(,ncol=6,nrow=15))
names(mr_results_corr)<-headers

mr_results_corr$outcome <- "cancer" 
mr_results_corr$exposure <- c("HMGCR", "HMGCR", "HMGCR", "NPC1L1", "NPC1L1", "NPC1L1", "PCSK9", "PCSK9", "PCSK9", "LDLR", "LDLR", "LDLR", "CETP", "CETP", "CETP")
mr_results_corr$method <- c("ivw", "egger", "median", "ivw", "egger", "median", "ivw", "egger", "median", "ivw", "egger", "median","ivw", "egger", "median")

mr_results_corr$b[1] <- (mr_ivw_HMGCR$Estimate)
mr_results_corr$b[2] <- (mr_egger_HMGCR$Estimate)
mr_results_corr$b[3] <- (mr_median_HMGCR$Estimate)
mr_results_corr$b[4] <- (mr_ivw_NPC1L1$Estimate)
mr_results_corr$b[5] <- (mr_egger_NPC1L1$Estimate)
mr_results_corr$b[6] <- (mr_median_NPC1L1$Estimate)
mr_results_corr$b[7] <- (mr_ivw_PCSK9$Estimate)
mr_results_corr$b[8] <- (mr_egger_PCSK9$Estimate)
mr_results_corr$b[9] <- (mr_median_PCSK9$Estimate)
mr_results_corr$b[10] <- (mr_ivw_LDLR$Estimate)
mr_results_corr$b[11] <- (mr_egger_LDLR$Estimate)
mr_results_corr$b[12] <- (mr_median_LDLR$Estimate)
mr_results_corr$b[13] <- (mr_ivw_CETP$Estimate)
mr_results_corr$b[14] <- (mr_egger_CETP$Estimate)
mr_results_corr$b[15] <- (mr_median_CETP$Estimate)

mr_results_corr$se[1] <- mr_ivw_HMGCR$StdError
mr_results_corr$se[2] <- mr_egger_HMGCR$StdError.Est
mr_results_corr$se[3] <- mr_median_HMGCR$StdError
mr_results_corr$se[4] <- mr_ivw_NPC1L1$StdError
mr_results_corr$se[5] <- mr_egger_NPC1L1$StdError.Est
mr_results_corr$se[6] <- mr_median_NPC1L1$StdError
mr_results_corr$se[7] <- mr_ivw_PCSK9$StdError
mr_results_corr$se[8] <- mr_egger_PCSK9$StdError.Est
mr_results_corr$se[9] <- mr_median_PCSK9$StdError
mr_results_corr$se[10] <- mr_ivw_LDLR$StdError
mr_results_corr$se[11] <- mr_egger_LDLR$StdError.Est
mr_results_corr$se[12] <- mr_median_LDLR$StdError
mr_results_corr$se[13] <- mr_ivw_CETP$StdError
mr_results_corr$se[14] <- mr_egger_CETP$StdError.Est
mr_results_corr$se[15] <- mr_median_CETP$StdError

mr_results_corr$pval[1] <- mr_ivw_HMGCR$Pvalue
mr_results_corr$pval[2] <- mr_egger_HMGCR$Pvalue.Est
mr_results_corr$pval[3] <- mr_median_HMGCR$Pvalue
mr_results_corr$pval[4] <- mr_ivw_NPC1L1$Pvalue
mr_results_corr$pval[5] <- mr_egger_NPC1L1$Pvalue.Est
mr_results_corr$pval[6] <- mr_median_NPC1L1$Pvalue
mr_results_corr$pval[7] <- mr_ivw_PCSK9$Pvalue
mr_results_corr$pval[8] <- mr_egger_PCSK9$Pvalue.Est
mr_results_corr$pval[9] <- mr_median_PCSK9$Pvalue
mr_results_corr$pval[10] <- mr_ivw_LDLR$Pvalue
mr_results_corr$pval[11] <- mr_egger_LDLR$Pvalue.Est
mr_results_corr$pval[12] <- mr_median_LDLR$Pvalue
mr_results_corr$pval[13] <- mr_ivw_CETP$Pvalue
mr_results_corr$pval[14] <- mr_egger_CETP$Pvalue.Est
mr_results_corr$pval[15] <- mr_median_CETP$Pvalue

# Estimate odds ratio and 95% confidence interval
mr_results_corr$or <- exp(mr_results_corr$b)
mr_results_corr$cil <- exp(mr_results_corr$b-1.96*mr_results_corr$se)
mr_results_corr$ciu <- exp(mr_results_corr$b+1.96*mr_results_corr$se)

write.csv(mr_results_corr, "~/analysis/drug_cancer/results_corr.csv")

#Hetereogenity test
headers<-c("outcome","exposure","method","Q","Q_df", "Q_pval")
mr_heterogenity_test_corr <- as.data.frame(matrix(,ncol=6,nrow=10))
names(mr_heterogenity_test_corr)<-headers
mr_heterogenity_test_corr$outcome <- "cancer" 
mr_heterogenity_test_corr$exposure <- c("HMGCR", "HMGCR", "NPC1L1", "NPC1L1", "PCSK9", "PCSK9", "LDLR", "LDLR", "CETP", "CETP")
mr_heterogenity_test_corr$method <- c("IVW", "MR Egger", "IVW", "MR Egger", "IVW", "MR Egger", "IVW", "MR Egger", "IVW", "MR Egger")

mr_heterogenity_test_corr$Q[1] <- mr_ivw_HMGCR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[2] <- mr_egger_HMGCR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[3] <- mr_ivw_NPC1L1$Heter.Stat[1]
mr_heterogenity_test_corr$Q[4] <- mr_egger_NPC1L1$Heter.Stat[1]
mr_heterogenity_test_corr$Q[5] <- mr_ivw_PCSK9$Heter.Stat[1]
mr_heterogenity_test_corr$Q[6] <- mr_egger_PCSK9$Heter.Stat[1]
mr_heterogenity_test_corr$Q[7] <- mr_ivw_LDLR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[8] <- mr_egger_LDLR$Heter.Stat[1]
mr_heterogenity_test_corr$Q[9] <- mr_ivw_CETP$Heter.Stat[1]
mr_heterogenity_test_corr$Q[10] <- mr_egger_CETP$Heter.Stat[1]

mr_heterogenity_test_corr$Q_df[1] <- mr_ivw_HMGCR$SNPs - 1
mr_heterogenity_test_corr$Q_df[2] <- mr_egger_HMGCR$SNPs - 2
mr_heterogenity_test_corr$Q_df[3] <- mr_ivw_NPC1L1$SNPs - 1
mr_heterogenity_test_corr$Q_df[4] <- mr_egger_NPC1L1$SNPs - 2
mr_heterogenity_test_corr$Q_df[5] <- mr_ivw_PCSK9$SNPs - 1
mr_heterogenity_test_corr$Q_df[6] <- mr_egger_PCSK9$SNPs - 2
mr_heterogenity_test_corr$Q_df[7] <- mr_ivw_LDLR$SNPs - 1
mr_heterogenity_test_corr$Q_df[8] <- mr_egger_LDLR$SNPs - 2 
mr_heterogenity_test_corr$Q_df[9] <- mr_ivw_CETP$SNPs - 1
mr_heterogenity_test_corr$Q_df[10] <- mr_egger_CETP$SNPs - 2 

mr_heterogenity_test_corr$Q_pval[1] <- mr_ivw_HMGCR$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[2] <- mr_egger_HMGCR$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[3] <- mr_ivw_NPC1L1$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[4] <- mr_egger_NPC1L1$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[5] <- mr_ivw_PCSK9$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[6] <- mr_egger_PCSK9$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[7] <- mr_ivw_LDLR$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[8] <- mr_egger_LDLR$Heter.Stat[2] 
mr_heterogenity_test_corr$Q_pval[9] <- mr_ivw_CETP$Heter.Stat[2]
mr_heterogenity_test_corr$Q_pval[10] <- mr_egger_CETP$Heter.Stat[2] 

write.csv(mr_heterogenity_test_corr,"~/analysis/drug_cancer/heterogeneity_results_corr.csv")

# Egger intercept 
headers<-c("outcome","exposure","egger_intercept","se","pval")
mr_pleiotropy_test_corr <- as.data.frame(matrix(,ncol=5,nrow=5))
names(mr_pleiotropy_test_corr)<-headers
mr_pleiotropy_test_corr$outcome <- "cancer" 
mr_pleiotropy_test_corr$exposure <- c("HMGCR", "NPC1L1","PCSK9","LDLR","CETP")
mr_pleiotropy_test_corr$egger_intercept[1] <- (mr_egger_HMGCR$Intercept)
mr_pleiotropy_test_corr$egger_intercept[2] <- (mr_egger_NPC1L1$Intercept)
mr_pleiotropy_test_corr$egger_intercept[3] <- (mr_egger_PCSK9$Intercept)
mr_pleiotropy_test_corr$egger_intercept[4] <- (mr_egger_LDLR$Intercept)
mr_pleiotropy_test_corr$egger_intercept[5] <- (mr_egger_CETP$Intercept)
mr_pleiotropy_test_corr$se[1] <- mr_egger_HMGCR$StdError.Int
mr_pleiotropy_test_corr$se[2] <- mr_egger_NPC1L1$StdError.Int
mr_pleiotropy_test_corr$se[3] <- mr_egger_PCSK9$StdError.Int
mr_pleiotropy_test_corr$se[4] <- mr_egger_LDLR$StdError.Int
mr_pleiotropy_test_corr$se[5] <- mr_egger_CETP$StdError.Int
mr_pleiotropy_test_corr$pval[1] <- mr_egger_HMGCR$Pvalue.Int
mr_pleiotropy_test_corr$pval[2] <- mr_egger_NPC1L1$Pvalue.Int
mr_pleiotropy_test_corr$pval[3] <- mr_egger_PCSK9$Pvalue.Int
mr_pleiotropy_test_corr$pval[4] <- mr_egger_LDLR$Pvalue.Int
mr_pleiotropy_test_corr$se[5] <- mr_egger_CETP$StdError.Int

write.csv(mr_pleiotropy_test_corr,"~/analysis/drug_cancer/pleiotropy_results_corr.csv")

########################################################
# MVMR BC # OXYSTEROL RECEPTOR
########################################################
# diagram
DiagrammeR::grViz("
      digraph mrdag {
      graph [rankdir=TB]
      node [shape=ellipse]
      U [label='Confounders']
      node [shape=plaintext, height=0.3, width=0.3]
      G [label='IVs']
      G1 [label='IVs2']
      X [label='Cholesterol-lowering drugs']
      X1 [label='Oxysterol']
      Y [label='Cancer']

      {rank = same; G; X; Y;}
      {rank = same; G1;X1[dir=back]}
    
      G -> X 
      G1 -> X1 
      X -> X1 
      U -> X
      U -> Y
      X -> Y
      X1 -> Y 
      }
      ", height = 200)%>%
  export_svg() %>%
  charToRaw %>% 
  rsvg_pdf("~/analysis/drug_cancer/mvmr.pdf")

########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
# Flip the results (to protective effect of LDL-lowering)
# i1 <- which(drug_exp_data$beta.exposure > 0, arr.ind = TRUE)
# drug_exp_data$beta.exposure[i1] <- as.numeric(drug_exp_data$beta.exposure[i1]) * -1
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_oxy<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_oxy

library(GWAS.utils)
dat_mvmr_oxy$MAF<-eaf2maf(c(dat_mvmr_oxy$eaf.exposure))

## OXY levels
oxy <- associations(dat_mvmr_oxy$SNP, "prot-a-2087", proxies = 0)

oxy <- oxy %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "oxy", sep = "."), -rsid)

mvmrdat_mvmr_oxy <- dat_mvmr_oxy %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(oxy) %>%
  mutate(harmon = case_when(ea == ea.oxy & nea == nea.oxy ~ 1,
                            ea == nea.oxy & nea == ea.oxy ~ -1,
                            TRUE ~ 0),
         beta.oxy = beta.oxy * harmon,
         z.oxy = beta.oxy / se.oxy,
         se.oxy = 1/sqrt(2 * MAF * (1 - MAF) * (n.oxy + (z.oxy^2))),
         beta.oxy = z.oxy * se.oxy) %>%
  select(-c(z.oxy, harmon, ea.oxy, nea.oxy, n.oxy))
str(mvmrdat_mvmr_oxy)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_oxy, mvmr_analysis(cbind(beta.drug, beta.oxy), beta.ca,
                                        cbind(se.drug, se.oxy), se.ca,
                                        exposures = c("drug", "oxy"))) %>%
  tibble

mvmr_res

i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_oxysterol_drug_cancer_results.csv")

########################################################
# MVMR BC # ESTROGEN RECEPTOR
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_estrogen<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_estrogen

library(GWAS.utils)
dat_mvmr_estrogen$MAF<-eaf2maf(c(dat_mvmr_estrogen$eaf.exposure))

## Estrogen receptor levels
estro <- associations(dat_mvmr_estrogen$SNP, "prot-a-991", proxies = 0)

estro <- estro %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "estro", sep = "."), -rsid)

mvmrdat_mvmr_estrogen <- dat_mvmr_estrogen %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(estro) %>%
  mutate(harmon = case_when(ea == ea.estro & nea == nea.estro ~ 1,
                            ea == nea.estro & nea == ea.estro ~ -1,
                            TRUE ~ 0),
         beta.estro = beta.estro * harmon,
         z.estro = beta.estro / se.estro,
         se.estro = 1/sqrt(2 * MAF * (1 - MAF) * (n.estro + (z.estro^2))),
         beta.estro = z.estro * se.estro) %>%
  select(-c(z.estro, harmon, ea.estro, nea.estro, n.estro))
str(mvmrdat_mvmr_estrogen)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_estrogen, mvmr_analysis(cbind(beta.drug, beta.estro), beta.ca,
                                        cbind(se.drug, se.estro), se.ca,
                                        exposures = c("drug", "estro"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_estrogen_drug_cancer_results.csv")

########################################################
# MVMR BC # Estradiol
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_estradiol<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_estradiol

library(GWAS.utils)
dat_mvmr_estradiol$MAF<-eaf2maf(c(dat_mvmr_estradiol$eaf.exposure))

## estradiol receptor levels
estradiol <- associations(dat_mvmr_estradiol$SNP, "ebi-a-GCST90020092", proxies = 0)

estradiol <- estradiol %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "estradiol", sep = "."), -rsid)

mvmrdat_mvmr_estradiol <- dat_mvmr_estradiol %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(estradiol) %>%
  mutate(harmon = case_when(ea == ea.estradiol & nea == nea.estradiol ~ 1,
                            ea == nea.estradiol & nea == ea.estradiol ~ -1,
                            TRUE ~ 0),
         beta.estradiol = beta.estradiol * harmon,
         z.estradiol = beta.estradiol / se.estradiol,
         se.estradiol = 1/sqrt(2 * MAF * (1 - MAF) * (n.estradiol + (z.estradiol^2))),
         beta.estradiol = z.estradiol * se.estradiol) %>%
  select(-c(z.estradiol, harmon, ea.estradiol, nea.estradiol, n.estradiol))
str(mvmrdat_mvmr_estradiol)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_estradiol, mvmr_analysis(cbind(beta.drug, beta.estradiol), beta.ca,
                                                      cbind(se.drug, se.estradiol), se.ca,
                                                      exposures = c("drug", "estradiol"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_estradiol_drug_cancer_results.csv")

########################################################
# MVMR BC # AGE AT MENARCHE
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_menarche<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_menarche


library(GWAS.utils)
dat_mvmr_menarche$MAF<-eaf2maf(c(dat_mvmr_menarche$eaf.exposure))


menarche <- associations(dat_mvmr_menarche$SNP, "ieu-a-1095", proxies = 0)

menarche <- menarche %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "menarche", sep = "."), -rsid)

mvmrdat_mvmr_menarche <- dat_mvmr_menarche %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(menarche) %>%
  mutate(harmon = case_when(ea == ea.menarche & nea == nea.menarche ~ 1,
                            ea == nea.menarche & nea == ea.menarche ~ -1,
                            TRUE ~ 0),
         beta.menarche = beta.menarche * harmon,
         z.menarche = beta.menarche / se.menarche,
         se.menarche = 1/sqrt(2 * MAF * (1 - MAF) * (n.menarche + (z.menarche^2))),
         beta.menarche = z.menarche * se.menarche) %>%
  select(-c(z.menarche, harmon, ea.menarche, nea.menarche, n.menarche))
str(mvmrdat_mvmr_menarche)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_menarche, mvmr_analysis(cbind(beta.drug, beta.menarche), beta.ca,
                                        cbind(se.drug, se.menarche), se.ca,
                                        exposures = c("drug", "menarche"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_menarche_drug_cancer_results.csv")

########################################################
# MVMR BC # AGE AT menopause
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_menopause<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_menopause


library(GWAS.utils)
dat_mvmr_menopause$MAF<-eaf2maf(c(dat_mvmr_menopause$eaf.exposure))

## Age at menopause levels
menopause <- associations(dat_mvmr_menopause$SNP, "ieu-b-4824", proxies = 0)

menopause <- menopause %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "menopause", sep = "."), -rsid)

mvmrdat_mvmr_menopause <- dat_mvmr_menopause %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(menopause) %>%
  mutate(harmon = case_when(ea == ea.menopause & nea == nea.menopause ~ 1,
                            ea == nea.menopause & nea == ea.menopause ~ -1,
                            TRUE ~ 0),
         beta.menopause = beta.menopause * harmon,
         z.menopause = beta.menopause / se.menopause,
         se.menopause = 1/sqrt(2 * MAF * (1 - MAF) * (n.menopause + (z.menopause^2))),
         beta.menopause = z.menopause * se.menopause) %>%
  select(-c(z.menopause, harmon, ea.menopause, nea.menopause, n.menopause))
str(mvmrdat_mvmr_menopause)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_menopause, mvmr_analysis(cbind(beta.drug, beta.menopause), beta.ca,
                                                       cbind(se.drug, se.menopause), se.ca,
                                                       exposures = c("drug", "menopause"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_menopause_drug_cancer_results.csv")

########################################################
# MVMR BC # BMI in females
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_bmi<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_bmi


library(GWAS.utils)
dat_mvmr_bmi$MAF<-eaf2maf(c(dat_mvmr_bmi$eaf.exposure))


bmi <- associations(dat_mvmr_bmi$SNP, "ieu-a-974", proxies = 0)

bmi <- bmi %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "bmi", sep = "."), -rsid)

mvmrdat_mvmr_bmi <- dat_mvmr_bmi %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(bmi) %>%
  mutate(harmon = case_when(ea == ea.bmi & nea == nea.bmi ~ 1,
                            ea == nea.bmi & nea == ea.bmi ~ -1,
                            TRUE ~ 0),
         beta.bmi = beta.bmi * harmon,
         z.bmi = beta.bmi / se.bmi,
         se.bmi = 1/sqrt(2 * MAF * (1 - MAF) * (n.bmi + (z.bmi^2))),
         beta.bmi = z.bmi * se.bmi) %>%
  select(-c(z.bmi, harmon, ea.bmi, nea.bmi, n.bmi))
str(mvmrdat_mvmr_bmi)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_bmi, mvmr_analysis(cbind(beta.drug, beta.bmi), beta.ca,
                                        cbind(se.drug, se.bmi), se.ca,
                                        exposures = c("drug", "bmi"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_bmi_drug_cancer_results.csv")

########################################################
# MVMR BC # BMI2 in All (females and males)
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_BMI2<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_BMI2


library(GWAS.utils)
dat_mvmr_BMI2$MAF<-eaf2maf(c(dat_mvmr_BMI2$eaf.exposure))

#
BMI2 <- associations(dat_mvmr_BMI2$SNP, "ieu-b-40", proxies = 0)

BMI2 <- BMI2 %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "BMI2", sep = "."), -rsid)

mvmrdat_mvmr_BMI2 <- dat_mvmr_BMI2 %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(BMI2) %>%
  mutate(harmon = case_when(ea == ea.BMI2 & nea == nea.BMI2 ~ 1,
                            ea == nea.BMI2 & nea == ea.BMI2 ~ -1,
                            TRUE ~ 0),
         beta.BMI2 = beta.BMI2 * harmon,
         z.BMI2 = beta.BMI2 / se.BMI2,
         se.BMI2 = 1/sqrt(2 * MAF * (1 - MAF) * (n.BMI2 + (z.BMI2^2))),
         beta.BMI2 = z.BMI2 * se.BMI2) %>%
  select(-c(z.BMI2, harmon, ea.BMI2, nea.BMI2, n.BMI2))
str(mvmrdat_mvmr_BMI2)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_BMI2, mvmr_analysis(cbind(beta.drug, beta.BMI2), beta.ca,
                                                 cbind(se.drug, se.BMI2), se.ca,
                                                 exposures = c("drug", "BMI2"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_BMI2_drug_cancer_results.csv")

########################################################
# MVMR BC # WC
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-61")

## 3. harmonise
dat_mvmr_WC<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_WC


library(GWAS.utils)
dat_mvmr_WC$MAF<-eaf2maf(c(dat_mvmr_WC$eaf.exposure))


WC <- associations(dat_mvmr_WC$SNP, "ieu-a-61", proxies = 0)

WC <- WC %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "WC", sep = "."), -rsid)

mvmrdat_mvmr_WC <- dat_mvmr_WC %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(WC) %>%
  mutate(harmon = case_when(ea == ea.WC & nea == nea.WC ~ 1,
                            ea == nea.WC & nea == ea.WC ~ -1,
                            TRUE ~ 0),
         beta.WC = beta.WC * harmon,
         z.WC = beta.WC / se.WC,
         se.WC = 1/sqrt(2 * MAF * (1 - MAF) * (n.WC + (z.WC^2))),
         beta.WC = z.WC * se.WC) %>%
  select(-c(z.WC, harmon, ea.WC, nea.WC, n.WC))
str(mvmrdat_mvmr_WC)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_WC, mvmr_analysis(cbind(beta.drug, beta.WC), beta.ca,
                                                 cbind(se.drug, se.WC), se.ca,
                                                 exposures = c("drug", "WC"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_WC_drug_cancer_results.csv")

########################################################
# MVMR BC # WHR
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-61")

## 3. harmonise
dat_mvmr_WHR<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_WHR


library(GWAS.utils)
dat_mvmr_WHR$MAF<-eaf2maf(c(dat_mvmr_WHR$eaf.exposure))


WHR <- associations(dat_mvmr_WHR$SNP, "ieu-a-73", proxies = 0)

WHR <- WHR %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "WHR", sep = "."), -rsid)

mvmrdat_mvmr_WHR <- dat_mvmr_WHR %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(WHR) %>%
  mutate(harmon = case_when(ea == ea.WHR & nea == nea.WHR ~ 1,
                            ea == nea.WHR & nea == ea.WHR ~ -1,
                            TRUE ~ 0),
         beta.WHR = beta.WHR * harmon,
         z.WHR = beta.WHR / se.WHR,
         se.WHR = 1/sqrt(2 * MAF * (1 - MAF) * (n.WHR + (z.WHR^2))),
         beta.WHR = z.WHR * se.WHR) %>%
  select(-c(z.WHR, harmon, ea.WHR, nea.WHR, n.WHR))
str(mvmrdat_mvmr_WHR)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_WHR, mvmr_analysis(cbind(beta.drug, beta.WHR), beta.ca,
                                                cbind(se.drug, se.WHR), se.ca,
                                                exposures = c("drug", "WHR"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_WHR_drug_cancer_results.csv")

################################################################################################################
# other outcome
########################################################
# MVMR BC +ER # OXYSTEROL RECEPTOR
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
# Flip the results (to protective effect of LDL-lowering)
# i1 <- which(drug_exp_data$beta.exposure > 0, arr.ind = TRUE)
# drug_exp_data$beta.exposure[i1] <- as.numeric(drug_exp_data$beta.exposure[i1]) * -1
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1127")

## 3. harmonise
dat_mvmr_oxy<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_oxy

library(GWAS.utils)
dat_mvmr_oxy$MAF<-eaf2maf(c(dat_mvmr_oxy$eaf.exposure))

## OXY levels
oxy <- associations(dat_mvmr_oxy$SNP, "prot-a-2087", proxies = 0)

oxy <- oxy %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "oxy", sep = "."), -rsid)

mvmrdat_mvmr_oxy <- dat_mvmr_oxy %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(oxy) %>%
  mutate(harmon = case_when(ea == ea.oxy & nea == nea.oxy ~ 1,
                            ea == nea.oxy & nea == ea.oxy ~ -1,
                            TRUE ~ 0),
         beta.oxy = beta.oxy * harmon,
         z.oxy = beta.oxy / se.oxy,
         se.oxy = 1/sqrt(2 * MAF * (1 - MAF) * (n.oxy + (z.oxy^2))),
         beta.oxy = z.oxy * se.oxy) %>%
  select(-c(z.oxy, harmon, ea.oxy, nea.oxy, n.oxy))
str(mvmrdat_mvmr_oxy)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_oxy, mvmr_analysis(cbind(beta.drug, beta.oxy), beta.ca,
                                                 cbind(se.drug, se.oxy), se.ca,
                                                 exposures = c("drug", "oxy"))) %>%
  tibble

mvmr_res

i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bcerp_oxysterol_drug_cancer_results.csv")

########################################################
# MVMR BC +ER # Estrogen 
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1127")

## 3. harmonise
dat_mvmr_estrogen<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_estrogen

library(GWAS.utils)
dat_mvmr_estrogen$MAF<-eaf2maf(c(dat_mvmr_estrogen$eaf.exposure))

## Estrogen receptor levels
estro <- associations(dat_mvmr_estrogen$SNP, "prot-a-991", proxies = 0)

estro <- estro %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "estro", sep = "."), -rsid)

mvmrdat_mvmr_estrogen <- dat_mvmr_estrogen %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(estro) %>%
  mutate(harmon = case_when(ea == ea.estro & nea == nea.estro ~ 1,
                            ea == nea.estro & nea == ea.estro ~ -1,
                            TRUE ~ 0),
         beta.estro = beta.estro * harmon,
         z.estro = beta.estro / se.estro,
         se.estro = 1/sqrt(2 * MAF * (1 - MAF) * (n.estro + (z.estro^2))),
         beta.estro = z.estro * se.estro) %>%
  select(-c(z.estro, harmon, ea.estro, nea.estro, n.estro))
str(mvmrdat_mvmr_estrogen)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_estrogen, mvmr_analysis(cbind(beta.drug, beta.estro), beta.ca,
                                                      cbind(se.drug, se.estro), se.ca,
                                                      exposures = c("drug", "estro"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bcerp_estrogen_drug_cancer_results.csv")

########################################################
# MVMR BC +ER # AGE AT MENARCHE
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1127")

## 3. harmonise
dat_mvmr_menarche<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_menarche


library(GWAS.utils)
dat_mvmr_menarche$MAF<-eaf2maf(c(dat_mvmr_menarche$eaf.exposure))


menarche <- associations(dat_mvmr_menarche$SNP, "ieu-a-1095", proxies = 0)

menarche <- menarche %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "menarche", sep = "."), -rsid)

mvmrdat_mvmr_menarche <- dat_mvmr_menarche %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(menarche) %>%
  mutate(harmon = case_when(ea == ea.menarche & nea == nea.menarche ~ 1,
                            ea == nea.menarche & nea == ea.menarche ~ -1,
                            TRUE ~ 0),
         beta.menarche = beta.menarche * harmon,
         z.menarche = beta.menarche / se.menarche,
         se.menarche = 1/sqrt(2 * MAF * (1 - MAF) * (n.menarche + (z.menarche^2))),
         beta.menarche = z.menarche * se.menarche) %>%
  select(-c(z.menarche, harmon, ea.menarche, nea.menarche, n.menarche))
str(mvmrdat_mvmr_menarche)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_menarche, mvmr_analysis(cbind(beta.drug, beta.menarche), beta.ca,
                                                      cbind(se.drug, se.menarche), se.ca,
                                                      exposures = c("drug", "menarche"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bcerp_menarche_drug_cancer_results.csv")

########################################################
# MVMR BC +ER # AGE AT menopause
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1127")

## 3. harmonise
dat_mvmr_menopause<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_menopause


library(GWAS.utils)
dat_mvmr_menopause$MAF<-eaf2maf(c(dat_mvmr_menopause$eaf.exposure))

## Age at menopause levels
menopause <- associations(dat_mvmr_menopause$SNP, "ieu-b-4824", proxies = 0)

menopause <- menopause %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "menopause", sep = "."), -rsid)

mvmrdat_mvmr_menopause <- dat_mvmr_menopause %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(menopause) %>%
  mutate(harmon = case_when(ea == ea.menopause & nea == nea.menopause ~ 1,
                            ea == nea.menopause & nea == ea.menopause ~ -1,
                            TRUE ~ 0),
         beta.menopause = beta.menopause * harmon,
         z.menopause = beta.menopause / se.menopause,
         se.menopause = 1/sqrt(2 * MAF * (1 - MAF) * (n.menopause + (z.menopause^2))),
         beta.menopause = z.menopause * se.menopause) %>%
  select(-c(z.menopause, harmon, ea.menopause, nea.menopause, n.menopause))
str(mvmrdat_mvmr_menopause)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_menopause, mvmr_analysis(cbind(beta.drug, beta.menopause), beta.ca,
                                                       cbind(se.drug, se.menopause), se.ca,
                                                       exposures = c("drug", "menopause"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bcerp_menopause_drug_cancer_results.csv")

########################################################
# MVMR BC +ER # BMI in females
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1127")

## 3. harmonise
dat_mvmr_bmi<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_bmi


library(GWAS.utils)
dat_mvmr_bmi$MAF<-eaf2maf(c(dat_mvmr_bmi$eaf.exposure))


bmi <- associations(dat_mvmr_bmi$SNP, "ieu-a-974", proxies = 0)

bmi <- bmi %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "bmi", sep = "."), -rsid)

mvmrdat_mvmr_bmi <- dat_mvmr_bmi %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(bmi) %>%
  mutate(harmon = case_when(ea == ea.bmi & nea == nea.bmi ~ 1,
                            ea == nea.bmi & nea == ea.bmi ~ -1,
                            TRUE ~ 0),
         beta.bmi = beta.bmi * harmon,
         z.bmi = beta.bmi / se.bmi,
         se.bmi = 1/sqrt(2 * MAF * (1 - MAF) * (n.bmi + (z.bmi^2))),
         beta.bmi = z.bmi * se.bmi) %>%
  select(-c(z.bmi, harmon, ea.bmi, nea.bmi, n.bmi))
str(mvmrdat_mvmr_bmi)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_bmi, mvmr_analysis(cbind(beta.drug, beta.bmi), beta.ca,
                                                 cbind(se.drug, se.bmi), se.ca,
                                                 exposures = c("drug", "bmi"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bcerp_bmi_drug_cancer_results.csv")

########################################################
# MVMR BC +ER # BMI2 in All (females and males)
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1127")

## 3. harmonise
dat_mvmr_BMI2<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_BMI2


library(GWAS.utils)
dat_mvmr_BMI2$MAF<-eaf2maf(c(dat_mvmr_BMI2$eaf.exposure))


BMI2 <- associations(dat_mvmr_BMI2$SNP, "ieu-b-40", proxies = 0)

BMI2 <- BMI2 %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "BMI2", sep = "."), -rsid)

mvmrdat_mvmr_BMI2 <- dat_mvmr_BMI2 %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(BMI2) %>%
  mutate(harmon = case_when(ea == ea.BMI2 & nea == nea.BMI2 ~ 1,
                            ea == nea.BMI2 & nea == ea.BMI2 ~ -1,
                            TRUE ~ 0),
         beta.BMI2 = beta.BMI2 * harmon,
         z.BMI2 = beta.BMI2 / se.BMI2,
         se.BMI2 = 1/sqrt(2 * MAF * (1 - MAF) * (n.BMI2 + (z.BMI2^2))),
         beta.BMI2 = z.BMI2 * se.BMI2) %>%
  select(-c(z.BMI2, harmon, ea.BMI2, nea.BMI2, n.BMI2))
str(mvmrdat_mvmr_BMI2)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_BMI2, mvmr_analysis(cbind(beta.drug, beta.BMI2), beta.ca,
                                                 cbind(se.drug, se.BMI2), se.ca,
                                                 exposures = c("drug", "BMI2"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bcerp_BMI2_drug_cancer_results.csv")

########################################################
# MVMR BC +ER # estradiol
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1127")

## 3. harmonise
dat_mvmr_estradiol<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_estradiol


library(GWAS.utils)
dat_mvmr_estradiol$MAF<-eaf2maf(c(dat_mvmr_estradiol$eaf.exposure))

## Age at estradiol levels
estradiol <- associations(dat_mvmr_estradiol$SNP, "ebi-a-GCST90020092", proxies = 0)

estradiol <- estradiol %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "estradiol", sep = "."), -rsid)

mvmrdat_mvmr_estradiol <- dat_mvmr_estradiol %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(estradiol) %>%
  mutate(harmon = case_when(ea == ea.estradiol & nea == nea.estradiol ~ 1,
                            ea == nea.estradiol & nea == ea.estradiol ~ -1,
                            TRUE ~ 0),
         beta.estradiol = beta.estradiol * harmon,
         z.estradiol = beta.estradiol / se.estradiol,
         se.estradiol = 1/sqrt(2 * MAF * (1 - MAF) * (n.estradiol + (z.estradiol^2))),
         beta.estradiol = z.estradiol * se.estradiol) %>%
  select(-c(z.estradiol, harmon, ea.estradiol, nea.estradiol, n.estradiol))
str(mvmrdat_mvmr_estradiol)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_estradiol, mvmr_analysis(cbind(beta.drug, beta.estradiol), beta.ca,
                                                       cbind(se.drug, se.estradiol), se.ca,
                                                       exposures = c("drug", "estradiol"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bcerp_estradiol_drug_cancer_results.csv")

########################################################
# MVMR BC +ER # WC
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1127")

## 3. harmonise
dat_mvmr_WC<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_WC


library(GWAS.utils)
dat_mvmr_WC$MAF<-eaf2maf(c(dat_mvmr_WC$eaf.exposure))

## Age at WC levels
WC <- associations(dat_mvmr_WC$SNP, "ieu-a-61", proxies = 0)

WC <- WC %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "WC", sep = "."), -rsid)

mvmrdat_mvmr_WC <- dat_mvmr_WC %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(WC) %>%
  mutate(harmon = case_when(ea == ea.WC & nea == nea.WC ~ 1,
                            ea == nea.WC & nea == ea.WC ~ -1,
                            TRUE ~ 0),
         beta.WC = beta.WC * harmon,
         z.WC = beta.WC / se.WC,
         se.WC = 1/sqrt(2 * MAF * (1 - MAF) * (n.WC + (z.WC^2))),
         beta.WC = z.WC * se.WC) %>%
  select(-c(z.WC, harmon, ea.WC, nea.WC, n.WC))
str(mvmrdat_mvmr_WC)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_WC, mvmr_analysis(cbind(beta.drug, beta.WC), beta.ca,
                                                       cbind(se.drug, se.WC), se.ca,
                                                       exposures = c("drug", "WC"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bcerp_WC_drug_cancer_results.csv")

########################################################
# MVMR BC +ER # WHR
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="HMGCR")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1127")

## 3. harmonise
dat_mvmr_WHR<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_WHR


library(GWAS.utils)
dat_mvmr_WHR$MAF<-eaf2maf(c(dat_mvmr_WHR$eaf.exposure))

## Age at WHR levels
WHR <- associations(dat_mvmr_WHR$SNP, "ieu-a-73", proxies = 0)

WHR <- WHR %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "WHR", sep = "."), -rsid)

mvmrdat_mvmr_WHR <- dat_mvmr_WHR %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(WHR) %>%
  mutate(harmon = case_when(ea == ea.WHR & nea == nea.WHR ~ 1,
                            ea == nea.WHR & nea == ea.WHR ~ -1,
                            TRUE ~ 0),
         beta.WHR = beta.WHR * harmon,
         z.WHR = beta.WHR / se.WHR,
         se.WHR = 1/sqrt(2 * MAF * (1 - MAF) * (n.WHR + (z.WHR^2))),
         beta.WHR = z.WHR * se.WHR) %>%
  select(-c(z.WHR, harmon, ea.WHR, nea.WHR, n.WHR))
str(mvmrdat_mvmr_WHR)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_WHR, mvmr_analysis(cbind(beta.drug, beta.WHR), beta.ca,
                                                cbind(se.drug, se.WHR), se.ca,
                                                exposures = c("drug", "WHR"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bcerp_WHR_drug_cancer_results.csv")  
 
################################################################################################################################################
################################################################################################################
# other exposure
################################################
### CETP AS DRUG THERAPY
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
# Flip the results (to protective effect of LDL-lowering)
# i1 <- which(drug_exp_data$beta.exposure > 0, arr.ind = TRUE)
# drug_exp_data$beta.exposure[i1] <- as.numeric(drug_exp_data$beta.exposure[i1]) * -1
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="CETP")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_oxy<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_oxy

library(GWAS.utils)
library(ieugwasr)
dat_mvmr_oxy$MAF<-eaf2maf(c(dat_mvmr_oxy$eaf.exposure))

## OXY levels
oxy <- associations(dat_mvmr_oxy$SNP, "prot-a-2087", proxies = 0)

oxy <- oxy %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "oxy", sep = "."), -rsid)

mvmrdat_mvmr_oxy <- dat_mvmr_oxy %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(oxy) %>%
  mutate(harmon = case_when(ea == ea.oxy & nea == nea.oxy ~ 1,
                            ea == nea.oxy & nea == ea.oxy ~ -1,
                            TRUE ~ 0),
         beta.oxy = beta.oxy * harmon,
         z.oxy = beta.oxy / se.oxy,
         se.oxy = 1/sqrt(2 * MAF * (1 - MAF) * (n.oxy + (z.oxy^2))),
         beta.oxy = z.oxy * se.oxy) %>%
  select(-c(z.oxy, harmon, ea.oxy, nea.oxy, n.oxy))
str(mvmrdat_mvmr_oxy)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_oxy, mvmr_analysis(cbind(beta.drug, beta.oxy), beta.ca,
                                                 cbind(se.drug, se.oxy), se.ca,
                                                 exposures = c("drug", "oxy"))) %>%
  tibble

mvmr_res

i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_oxysterol_drugCETP_cancer_results.csv")

########################################################
# MVMR BC # ESTROGEN RECEPTOR
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="CETP")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_estrogen<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_estrogen

library(GWAS.utils)
dat_mvmr_estrogen$MAF<-eaf2maf(c(dat_mvmr_estrogen$eaf.exposure))

## Estrogen receptor levels
estro <- associations(dat_mvmr_estrogen$SNP, "prot-a-991", proxies = 0)

estro <- estro %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "estro", sep = "."), -rsid)

mvmrdat_mvmr_estrogen <- dat_mvmr_estrogen %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(estro) %>%
  mutate(harmon = case_when(ea == ea.estro & nea == nea.estro ~ 1,
                            ea == nea.estro & nea == ea.estro ~ -1,
                            TRUE ~ 0),
         beta.estro = beta.estro * harmon,
         z.estro = beta.estro / se.estro,
         se.estro = 1/sqrt(2 * MAF * (1 - MAF) * (n.estro + (z.estro^2))),
         beta.estro = z.estro * se.estro) %>%
  select(-c(z.estro, harmon, ea.estro, nea.estro, n.estro))
str(mvmrdat_mvmr_estrogen)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_estrogen, mvmr_analysis(cbind(beta.drug, beta.estro), beta.ca,
                                                      cbind(se.drug, se.estro), se.ca,
                                                      exposures = c("drug", "estro"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_estrogen_drugCETP_cancer_results.csv")

########################################################
# MVMR BC # Estradiol
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="CETP")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_estradiol<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_estradiol

library(GWAS.utils)
dat_mvmr_estradiol$MAF<-eaf2maf(c(dat_mvmr_estradiol$eaf.exposure))

## estradiol receptor levels
estradiol <- associations(dat_mvmr_estradiol$SNP, "ebi-a-GCST90020092", proxies = 0)

estradiol <- estradiol %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "estradiol", sep = "."), -rsid)

mvmrdat_mvmr_estradiol <- dat_mvmr_estradiol %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(estradiol) %>%
  mutate(harmon = case_when(ea == ea.estradiol & nea == nea.estradiol ~ 1,
                            ea == nea.estradiol & nea == ea.estradiol ~ -1,
                            TRUE ~ 0),
         beta.estradiol = beta.estradiol * harmon,
         z.estradiol = beta.estradiol / se.estradiol,
         se.estradiol = 1/sqrt(2 * MAF * (1 - MAF) * (n.estradiol + (z.estradiol^2))),
         beta.estradiol = z.estradiol * se.estradiol) %>%
  select(-c(z.estradiol, harmon, ea.estradiol, nea.estradiol, n.estradiol))
str(mvmrdat_mvmr_estradiol)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_estradiol, mvmr_analysis(cbind(beta.drug, beta.estradiol), beta.ca,
                                                       cbind(se.drug, se.estradiol), se.ca,
                                                       exposures = c("drug", "estradiol"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_estradiol_drugCETP_cancer_results.csv")

########################################################
# MVMR BC # AGE AT MENARCHE
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="CETP")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_menarche<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_menarche


library(GWAS.utils)
dat_mvmr_menarche$MAF<-eaf2maf(c(dat_mvmr_menarche$eaf.exposure))


menarche <- associations(dat_mvmr_menarche$SNP, "ieu-a-1095", proxies = 0)

menarche <- menarche %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "menarche", sep = "."), -rsid)

mvmrdat_mvmr_menarche <- dat_mvmr_menarche %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(menarche) %>%
  mutate(harmon = case_when(ea == ea.menarche & nea == nea.menarche ~ 1,
                            ea == nea.menarche & nea == ea.menarche ~ -1,
                            TRUE ~ 0),
         beta.menarche = beta.menarche * harmon,
         z.menarche = beta.menarche / se.menarche,
         se.menarche = 1/sqrt(2 * MAF * (1 - MAF) * (n.menarche + (z.menarche^2))),
         beta.menarche = z.menarche * se.menarche) %>%
  select(-c(z.menarche, harmon, ea.menarche, nea.menarche, n.menarche))
str(mvmrdat_mvmr_menarche)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_menarche, mvmr_analysis(cbind(beta.drug, beta.menarche), beta.ca,
                                                      cbind(se.drug, se.menarche), se.ca,
                                                      exposures = c("drug", "menarche"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_menarche_drugCETP_cancer_results.csv")

########################################################
# MVMR BC # AGE AT menopause
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="CETP")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_menopause<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_menopause


library(GWAS.utils)
dat_mvmr_menopause$MAF<-eaf2maf(c(dat_mvmr_menopause$eaf.exposure))

## Age at menopause levels
menopause <- associations(dat_mvmr_menopause$SNP, "ieu-b-4824", proxies = 0)

menopause <- menopause %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "menopause", sep = "."), -rsid)

mvmrdat_mvmr_menopause <- dat_mvmr_menopause %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(menopause) %>%
  mutate(harmon = case_when(ea == ea.menopause & nea == nea.menopause ~ 1,
                            ea == nea.menopause & nea == ea.menopause ~ -1,
                            TRUE ~ 0),
         beta.menopause = beta.menopause * harmon,
         z.menopause = beta.menopause / se.menopause,
         se.menopause = 1/sqrt(2 * MAF * (1 - MAF) * (n.menopause + (z.menopause^2))),
         beta.menopause = z.menopause * se.menopause) %>%
  select(-c(z.menopause, harmon, ea.menopause, nea.menopause, n.menopause))
str(mvmrdat_mvmr_menopause)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_menopause, mvmr_analysis(cbind(beta.drug, beta.menopause), beta.ca,
                                                       cbind(se.drug, se.menopause), se.ca,
                                                       exposures = c("drug", "menopause"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_menopause_drugCETP_cancer_results.csv")

########################################################
# MVMR BC # BMI in females
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="CETP")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_bmi<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_bmi


library(GWAS.utils)
dat_mvmr_bmi$MAF<-eaf2maf(c(dat_mvmr_bmi$eaf.exposure))


bmi <- associations(dat_mvmr_bmi$SNP, "ieu-a-974", proxies = 0)

bmi <- bmi %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "bmi", sep = "."), -rsid)

mvmrdat_mvmr_bmi <- dat_mvmr_bmi %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(bmi) %>%
  mutate(harmon = case_when(ea == ea.bmi & nea == nea.bmi ~ 1,
                            ea == nea.bmi & nea == ea.bmi ~ -1,
                            TRUE ~ 0),
         beta.bmi = beta.bmi * harmon,
         z.bmi = beta.bmi / se.bmi,
         se.bmi = 1/sqrt(2 * MAF * (1 - MAF) * (n.bmi + (z.bmi^2))),
         beta.bmi = z.bmi * se.bmi) %>%
  select(-c(z.bmi, harmon, ea.bmi, nea.bmi, n.bmi))
str(mvmrdat_mvmr_bmi)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_bmi, mvmr_analysis(cbind(beta.drug, beta.bmi), beta.ca,
                                                 cbind(se.drug, se.bmi), se.ca,
                                                 exposures = c("drug", "bmi"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_bmi_drugCETP_cancer_results.csv")

########################################################
# MVMR BC # WC
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="CETP")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_WC<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_WC


library(GWAS.utils)
dat_mvmr_WC$MAF<-eaf2maf(c(dat_mvmr_WC$eaf.exposure))


WC <- associations(dat_mvmr_WC$SNP, "ieu-a-61", proxies = 0)

WC <- WC %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "WC", sep = "."), -rsid)

mvmrdat_mvmr_WC <- dat_mvmr_WC %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(WC) %>%
  mutate(harmon = case_when(ea == ea.WC & nea == nea.WC ~ 1,
                            ea == nea.WC & nea == ea.WC ~ -1,
                            TRUE ~ 0),
         beta.WC = beta.WC * harmon,
         z.WC = beta.WC / se.WC,
         se.WC = 1/sqrt(2 * MAF * (1 - MAF) * (n.WC + (z.WC^2))),
         beta.WC = z.WC * se.WC) %>%
  select(-c(z.WC, harmon, ea.WC, nea.WC, n.WC))
str(mvmrdat_mvmr_WC)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_WC, mvmr_analysis(cbind(beta.drug, beta.WC), beta.ca,
                                                 cbind(se.drug, se.WC), se.ca,
                                                 exposures = c("drug", "WC"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_WC_drugCETP_cancer_results.csv")

########################################################
# MVMR BC # WHR
########################################################
# Set working directory
setwd("~/analysis/")

#1. SNP-EXPOSURE SUMMARY DATA
drug_exp_data <- read_exposure_data("drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
head(drug_exp_data)
dim(drug_exp_data)# 36 12
## subset 
drug_exp_data<- subset(drug_exp_data, exposure=="CETP")
head(drug_exp_data)
dim(drug_exp_data)

## 2. read outcome cancer
outcome_data <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "ieu-a-1126")

## 3. harmonise
dat_mvmr_WHR<- harmonise_data(drug_exp_data, outcome_data)
dat_mvmr_WHR


library(GWAS.utils)
dat_mvmr_WHR$MAF<-eaf2maf(c(dat_mvmr_WHR$eaf.exposure))


WHR <- associations(dat_mvmr_WHR$SNP, "ieu-a-73", proxies = 0)

WHR <- WHR %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "WHR", sep = "."), -rsid)

mvmrdat_mvmr_WHR <- dat_mvmr_WHR %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF,
         beta.drug = beta.exposure, 
         se.drug = se.exposure,
         p.drug = pval.exposure,
         beta.ca = beta.outcome, 
         se.ca = se.outcome, 
         p.ca = pval.outcome) %>%
  inner_join(WHR) %>%
  mutate(harmon = case_when(ea == ea.WHR & nea == nea.WHR ~ 1,
                            ea == nea.WHR & nea == ea.WHR ~ -1,
                            TRUE ~ 0),
         beta.WHR = beta.WHR * harmon,
         z.WHR = beta.WHR / se.WHR,
         se.WHR = 1/sqrt(2 * MAF * (1 - MAF) * (n.WHR + (z.WHR^2))),
         beta.WHR = z.WHR * se.WHR) %>%
  select(-c(z.WHR, harmon, ea.WHR, nea.WHR, n.WHR))
str(mvmrdat_mvmr_WHR)

## MVMR functions
setwd("~/")
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat_mvmr_WHR, mvmr_analysis(cbind(beta.drug, beta.WHR), beta.ca,
                                                cbind(se.drug, se.WHR), se.ca,
                                                exposures = c("drug", "WHR"))) %>%
  tibble

mvmr_res
i1 <- mvmr_res$b > 0
mvmr_res$b[i1] <- -1*(mvmr_res$b[i1])
# Estimate odds ratio and 95% confidence interval
mvmr_res$or <- exp(mvmr_res$b)
mvmr_res$cil <- exp(mvmr_res$b-1.96*mvmr_res$se)
mvmr_res$ciu <- exp(mvmr_res$b+1.96*mvmr_res$se)
rio::export(mvmr_res, "~/analysis/drug_cancer/mvmr_bc_WHR_drugCETP_cancer_results.csv")

################################################
save.image("~/analysis/drug_cancer/drug_cancer_analysis.RData")
############END########################