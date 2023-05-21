library(tidyverse)
library(ieugwasr)

oxy_exp_dat <- extract_instruments(outcomes=c('prot-a-2087','prot-a-2152'))
BCAC_out_dat <- extract_outcome_data(snps = oxy_exp_dat$SNP, outcomes = c('ieu-a-1127')) #'prot-a-2087', ER+ 'ieu-a-1127', ER-'ieu-a-1128'
#
dat<- harmonise_data(exposure_dat = oxy_exp_dat,outcome_dat = BCAC_out_dat)
dat

library(GWAS.utils)
dat$MAF<-eaf2maf(c(dat$eaf.exposure))

fi <- associations(dat$SNP, c("ieu-b-115"), proxies = 0)
fi <- fi %>%
  select(-c(chr, position, id, trait, eaf)) %>%
  rename_with(~paste(.x, "fi", sep = "."), -rsid)

mvmrdat <- dat %>%
  select(rsid = SNP, 
         ea = effect_allele.exposure, 
         nea = other_allele.exposure, 
         MAF, #= eaf.exposure,
         beta.prot = beta.exposure, se.prot = se.exposure,p.prot = pval.exposure,
         beta.bcac = beta.outcome, se.bcac = se.outcome, p.bcac = pval.outcome) %>%
  inner_join(fi) %>%
  mutate(harmon = case_when(ea == ea.fi & nea == nea.fi ~ 1,
                            ea == nea.fi & nea == ea.fi ~ -1,
                            TRUE ~ 0),
         beta.fi = beta.fi * harmon,
         z.fi = beta.fi / se.fi,
         se.fi = 1/sqrt(2 * MAF * (1 - MAF) * (n.fi + (z.fi^2))),
         beta.fi = z.fi * se.fi) %>%
  select(-c(z.fi, harmon, ea.fi, nea.fi, n.fi))

mvmrdat

setwd("/Users/hu7454po/Desktop/Skeleton_projects/oxysterol")

## Clumping
# toclump <- select(mvmrdat, rsid, pval = p.prot)
# 
# ld_reflookup(mvmrdat$rsid, pop = "EUR")
# 
# clumped <- ld_clump(toclump, clump_kb = 500, clump_r2 = 0.01,
#                     plink_bin = genetics.binaRies::get_plink_binary(),
#                     bfile = "/ludc/Home/daniel_c/dva/files/1kg_ref/EUR")
# 
# mvmrdat_c <- filter(mvmrdat, rsid %in% clumped$rsid)

## MVMR functions
source(list.files(pattern = "mrfx.R$"))
mvmr_res <- with(mvmrdat, mvmr_analysis(cbind(beta.prot, beta.fi), beta.bcac,
                                          cbind(se.prot, se.fi), se.bcac,
                                          exposures = c("prot", "fi"))) %>%
  tibble

mvmr_res

rio::export(mvmr_res, "mvmr_res.tsv")
