# Using the MR sources
library(AER)
library (boot)
library(foreign)
library(rio)
library(TwoSampleMR)
library(dplyr)
library(data.table)
# https://www.alexstephenson.me/post/2022-04-02-standard-errors-and-the-delta-method/
rm(list = ls())
########
ratio <- function(x){
  x[1]*x[2]/x[3]
}
######## APOC3 ######## 
# this function by default will output this in transpose form
# estimates
apoc3_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                               c(-0.57, # first half (fixed)
                                 -0.02, # second half
                                 -0.04 # total 
                                 )))
# Variance covariance matrix 
# SEs
apoc3_bc_vcov <- diag(c(0.19, # first half (fixed)
                        0.02, # second half
                        0.05) # total 
                      ^2)

## Estimate the standard error
apoc3_bc_se_b <- sqrt(t(apoc3_bc_grad_g) %*% apoc3_bc_vcov %*%apoc3_bc_grad_g)
apoc3_bc_se_b

# 95% CIs
apoc3_bc_estimate = .33 
apoc3_bc_lwr_b = apoc3_bc_estimate - 1.96*apoc3_bc_se_b;apoc3_bc_lwr_b
apoc3_bc_upp_b = apoc3_bc_estimate + 1.96*apoc3_bc_se_b;apoc3_bc_upp_b

#p value
apoc3_SE = apoc3_bc_upp_b - apoc3_bc_lwr_b/(2*1.96);
apoc3_z = apoc3_bc_estimate/apoc3_SE;
apoc3_P = exp(-0.717*apoc3_z - 0.416*apoc3_z^2)
apoc3_P # 0.8242119

######## ER - BC ######## 
# this function by default will output this in transpose form
# estimates
apoc3_en_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                        c(-0.57, # first half (fixed)
                                          -0.02, # second half
                                          -0.01	 # total 
                                        )))
# Variance covariance matrix 
# SEs
apoc3_en_bc_vcov <- diag(c(0.19, # first half (fixed)
                        0.03, # second half
                        0.18) # total 
                      ^2)

## Estimate the standard error
apoc3_en_bc_se_b <- sqrt(t(apoc3_en_bc_grad_g) %*% apoc3_en_bc_vcov %*%apoc3_en_bc_grad_g)
apoc3_en_bc_se_b

# 95% CIs
apoc3_en_bc_estimate = .92 
apoc3_en_bc_lwr_b = apoc3_en_bc_estimate - 1.96*apoc3_en_bc_se_b;apoc3_en_bc_lwr_b
apoc3_en_bc_upp_b = apoc3_en_bc_estimate + 1.96*apoc3_en_bc_se_b;apoc3_en_bc_upp_b

#p value
apoc3_en_bc_SE = apoc3_en_bc_upp_b - apoc3_en_bc_lwr_b/(2*1.96);
apoc3_en_bc_z = apoc3_en_bc_estimate/apoc3_en_bc_SE;
apoc3_en_bc_P = exp(-0.717*apoc3_en_bc_z - 0.416*apoc3_en_bc_z^2)
apoc3_en_bc_P # 

######## ER + BC ######## 
# this function by default will output this in transpose form
# estimates
apoc3_ep_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                           c(-0.57, # first half (fixed)
                                             -0.05, # second half
                                             -0.03	# total 
                                           )))
# Variance covariance matrix 
# SEs
apoc3_ep_bc_vcov <- diag(c(0.19, # first half (fixed)
                           0.02, # second half
                           0.06) # total 
                         ^2)

## Estimate the standard error
apoc3_ep_bc_se_b <- sqrt(t(apoc3_ep_bc_grad_g) %*% apoc3_ep_bc_vcov %*%apoc3_ep_bc_grad_g)
apoc3_ep_bc_se_b

# 95% CIs
apoc3_ep_bc_estimate = .78
apoc3_ep_bc_lwr_b = apoc3_ep_bc_estimate - 1.96*apoc3_ep_bc_se_b;apoc3_ep_bc_lwr_b
apoc3_ep_bc_upp_b = apoc3_ep_bc_estimate + 1.96*apoc3_ep_bc_se_b;apoc3_ep_bc_upp_b

#p value
apoc3_ep_SE = apoc3_en_bc_upp_b - apoc3_en_bc_lwr_b/(2*1.96);
apoc3_ep_z = apoc3_ep_bc_estimate/apoc3_ep_SE;
apoc3_ep_P = exp(-0.717*apoc3_ep_z - 0.416*apoc3_ep_z^2)
apoc3_ep_P # 

######## CETP ######## 
# this function by default will output this in transpose form
# estimates
cetp_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                        c(-0.25	, # first half (fixed)
                                          -0.02, # second half
                                          -0.25 # total 
                                        )))
# Variance covariance matrix 
# SEs
cetp_bc_vcov <- diag(c(0.30, # first half (fixed)
                        0.02, # second half
                       0.08) # total 
                      ^2)

## Estimate the standard error
cetp_bc_se_b <- sqrt(t(cetp_bc_grad_g) %*% cetp_bc_vcov %*%cetp_bc_grad_g)
cetp_bc_se_b

# 95% CIs
cetp_bc_estimate = .02 
cetp_bc_lwr_b = cetp_bc_estimate - 1.96*cetp_bc_se_b;cetp_bc_lwr_b
cetp_bc_upp_b = cetp_bc_estimate + 1.96*cetp_bc_se_b;cetp_bc_upp_b

#p value
cetp_SE = cetp_bc_upp_b - cetp_bc_lwr_b/(2*1.96);
cetp_z = cetp_bc_estimate/cetp_SE;
cetp_P = exp(-0.717*cetp_z - 0.416*cetp_z^2)
cetp_P # 

######## ER - BC ######## 
# this function by default will output this in transpose form
# estimates
cetp_en_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                           c(-0.25, # first half (fixed)
                                             -0.02, # second half
                                             -0.38		 # total 
                                           )))
# Variance covariance matrix 
# SEs
cetp_en_bc_vcov <- diag(c(0.30, # first half (fixed)
                          0.03, # second half
                          0.14) # total 
                         ^2)

## Estimate the standard error
cetp_en_bc_se_b <- sqrt(t(cetp_en_bc_grad_g) %*% cetp_en_bc_vcov %*%cetp_en_bc_grad_g)
cetp_en_bc_se_b

# 95% CIs
cetp_en_bc_estimate = .01 
cetp_en_bc_lwr_b = cetp_en_bc_estimate - 1.96*cetp_en_bc_se_b;cetp_en_bc_lwr_b
cetp_en_bc_upp_b = cetp_en_bc_estimate + 1.96*cetp_en_bc_se_b;cetp_en_bc_upp_b

#p value
cetp_en_bc_SE = cetp_en_bc_upp_b - cetp_en_bc_lwr_b/(2*1.96);
cetp_en_bc_z = cetp_en_bc_estimate/cetp_en_bc_SE;
cetp_en_bc_P = exp(-0.717*cetp_en_bc_z - 0.416*cetp_en_bc_z^2)
cetp_en_bc_P # 

######## ER + BC ######## 
# this function by default will output this in transpose form
# estimates
cetp_ep_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                           c(-0.25, # first half (fixed)
                                             -0.05, # second half
                                             -0.21		# total 
                                           )))
# Variance covariance matrix 
# SEs
cetp_ep_bc_vcov <- diag(c(0.30, # first half (fixed)
                           0.02, # second half
                          0.09) # total 
                         ^2)

## Estimate the standard error
cetp_ep_bc_se_b <- sqrt(t(cetp_ep_bc_grad_g) %*% cetp_ep_bc_vcov %*%cetp_ep_bc_grad_g)
cetp_ep_bc_se_b

# 95% CIs
cetp_ep_bc_estimate = .05
cetp_ep_bc_lwr_b = cetp_ep_bc_estimate - 1.96*cetp_ep_bc_se_b;cetp_ep_bc_lwr_b
cetp_ep_bc_upp_b = cetp_ep_bc_estimate + 1.96*cetp_ep_bc_se_b;cetp_ep_bc_upp_b

#p value
cetp_ep_SE = cetp_en_bc_upp_b - cetp_en_bc_lwr_b/(2*1.96);
cetp_ep_z = cetp_ep_bc_estimate/cetp_ep_SE;
cetp_ep_P = exp(-0.717*cetp_ep_z - 0.416*cetp_ep_z^2)
cetp_ep_P # 

######## HMGCR ######## 
# this function by default will output this in transpose form
# estimates
hmgcr_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                       c(-0.09	, # first half (fixed)
                                         -0.02, # second half
                                         -0.24 # total 
                                       )))
# Variance covariance matrix 
# SEs
hmgcr_bc_vcov <- diag(c(0.25, # first half (fixed)
                       0.02, # second half
                       0.06) # total 
                     ^2)

## Estimate the standard error
hmgcr_bc_se_b <- sqrt(t(hmgcr_bc_grad_g) %*% hmgcr_bc_vcov %*%hmgcr_bc_grad_g)
hmgcr_bc_se_b

# 95% CIs
hmgcr_bc_estimate = .01
hmgcr_bc_lwr_b = hmgcr_bc_estimate - 1.96*hmgcr_bc_se_b;hmgcr_bc_lwr_b
hmgcr_bc_upp_b = hmgcr_bc_estimate + 1.96*hmgcr_bc_se_b;hmgcr_bc_upp_b

#p value
hmgcr_SE = hmgcr_bc_upp_b - hmgcr_bc_lwr_b/(2*1.96);
hmgcr_z = hmgcr_bc_estimate/hmgcr_SE;
hmgcr_P = exp(-0.717*hmgcr_z - 0.416*hmgcr_z^2)
hmgcr_P # 

######## ER - BC ######## 
# this function by default will output this in transpose form
# estimates
hmgcr_en_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                          c(-0.25, # first half (fixed)
                                            -0.02, # second half
                                            -0.17	 # total 
                                          )))
# Variance covariance matrix 
# SEs
hmgcr_en_bc_vcov <- diag(c(0.30, # first half (fixed)
                          0.03, # second half
                          0.11) # total 
                        ^2)

## Estimate the standard error
hmgcr_en_bc_se_b <- sqrt(t(hmgcr_en_bc_grad_g) %*% hmgcr_en_bc_vcov %*%hmgcr_en_bc_grad_g)
hmgcr_en_bc_se_b

# 95% CIs
hmgcr_en_bc_estimate = .01 
hmgcr_en_bc_lwr_b = hmgcr_en_bc_estimate - 1.96*hmgcr_en_bc_se_b;hmgcr_en_bc_lwr_b
hmgcr_en_bc_upp_b = hmgcr_en_bc_estimate + 1.96*hmgcr_en_bc_se_b;hmgcr_en_bc_upp_b

#p value
hmgcr_en_bc_SE = hmgcr_en_bc_upp_b - hmgcr_en_bc_lwr_b/(2*1.96);
hmgcr_en_bc_z = hmgcr_en_bc_estimate/hmgcr_en_bc_SE;
hmgcr_en_bc_P = exp(-0.717*hmgcr_en_bc_z - 0.416*hmgcr_en_bc_z^2)
hmgcr_en_bc_P # 

######## ER + BC ######## 
# this function by default will output this in transpose form
# estimates
hmgcr_ep_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                          c(-0.25, # first half (fixed)
                                            -0.05, # second half
                                            -0.28	# total 
                                          )))
# Variance covariance matrix 
# SEs
hmgcr_ep_bc_vcov <- diag(c(0.30, # first half (fixed)
                           0.02, # second half
                          0.08) # total 
                        ^2)

## Estimate the standard error
hmgcr_ep_bc_se_b <- sqrt(t(hmgcr_ep_bc_grad_g) %*% hmgcr_ep_bc_vcov %*%hmgcr_ep_bc_grad_g)
hmgcr_ep_bc_se_b

# 95% CIs
hmgcr_ep_bc_estimate = .02
hmgcr_ep_bc_lwr_b = hmgcr_ep_bc_estimate - 1.96*hmgcr_ep_bc_se_b;hmgcr_ep_bc_lwr_b
hmgcr_ep_bc_upp_b = hmgcr_ep_bc_estimate + 1.96*hmgcr_ep_bc_se_b;hmgcr_ep_bc_upp_b

#p value
hmgcr_ep_SE = hmgcr_en_bc_upp_b - hmgcr_en_bc_lwr_b/(2*1.96);
hmgcr_ep_z = hmgcr_ep_bc_estimate/hmgcr_ep_SE;
hmgcr_ep_P = exp(-0.717*hmgcr_ep_z - 0.416*hmgcr_ep_z^2)
hmgcr_ep_P # 

######## LDLR ######## 
# this function by default will output this in transpose form
# estimates
ldlr_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                        c(-0.665150720188845	, # first half (fixed)
                                          -0.0225771723876219, # second half
                                          -0.01245754458971 # total 
                                        )))
# Variance covariance matrix 
# SEs
ldlr_bc_vcov <- diag(c(0.191070875905666, # first half (fixed)
                       0.0154301328222015, # second half
                       0.0502419049219536) # total 
                      ^2)

## Estimate the standard error
ldlr_bc_se_b <- sqrt(t(ldlr_bc_grad_g) %*% ldlr_bc_vcov %*%ldlr_bc_grad_g)
ldlr_bc_se_b

# 95% CIs
ldlr_bc_estimate = 1.21
ldlr_bc_lwr_b = ldlr_bc_estimate - 1.96*ldlr_bc_se_b;ldlr_bc_lwr_b
ldlr_bc_upp_b = ldlr_bc_estimate + 1.96*ldlr_bc_se_b;ldlr_bc_upp_b

#p value
ldlr_SE = ldlr_bc_upp_b - ldlr_bc_lwr_b/(2*1.96);
ldlr_z = ldlr_bc_estimate/ldlr_SE;
ldlr_P = exp(-0.717*ldlr_z - 0.416*ldlr_z^2)
ldlr_P # 

######## ER - BC ######## 
# this function by default will output this in transpose form
# estimates
ldlr_en_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                           c(-0.665150720188845, # first half (fixed)
                                             -0.022550430586603, # second half
                                             -0.0954207813069684 # total 
                                           )))
# Variance covariance matrix 
# SEs
ldlr_en_bc_vcov <- diag(c(0.191070875905666, # first half (fixed)
                          0.0279869477723851, # second half
                          0.0968770138310482) # total 
                         ^2)

## Estimate the standard error
ldlr_en_bc_se_b <- sqrt(t(ldlr_en_bc_grad_g) %*% ldlr_en_bc_vcov %*%ldlr_en_bc_grad_g)
ldlr_en_bc_se_b

# 95% CIs
ldlr_en_bc_estimate = .016
ldlr_en_bc_lwr_b = ldlr_en_bc_estimate - 1.96*ldlr_en_bc_se_b;ldlr_en_bc_lwr_b
ldlr_en_bc_upp_b = ldlr_en_bc_estimate + 1.96*ldlr_en_bc_se_b;ldlr_en_bc_upp_b

#p value
ldlr_en_bc_SE = ldlr_en_bc_upp_b - ldlr_en_bc_lwr_b/(2*1.96);
ldlr_en_bc_z = ldlr_en_bc_estimate/ldlr_en_bc_SE;
ldlr_en_bc_P = exp(-0.717*ldlr_en_bc_z - 0.416*ldlr_en_bc_z^2)
ldlr_en_bc_P # 

######## ER + BC ######## 
# this function by default will output this in transpose form
# estimates
ldlr_ep_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                           c(-0.665150720188845, # first half (fixed)
                                             -0.0464566847330293, # second half
                                             -0.0299485588677092	# total 
                                           )))
# Variance covariance matrix 
# SEs
ldlr_ep_bc_vcov <- diag(c(0.191070875905666, # first half (fixed)
                          0.0243321227922063, # second half
                          0.0599763881157877) # total 
                         ^2)

## Estimate the standard error
ldlr_ep_bc_se_b <- sqrt(t(ldlr_ep_bc_grad_g) %*% ldlr_ep_bc_vcov %*%ldlr_ep_bc_grad_g)
ldlr_ep_bc_se_b

# 95% CIs
ldlr_ep_bc_estimate = 1.03
ldlr_ep_bc_lwr_b = ldlr_ep_bc_estimate - 1.96*ldlr_ep_bc_se_b;ldlr_ep_bc_lwr_b
ldlr_ep_bc_upp_b = ldlr_ep_bc_estimate + 1.96*ldlr_ep_bc_se_b;ldlr_ep_bc_upp_b

#p value
ldlr_ep_SE = ldlr_en_bc_upp_b - ldlr_en_bc_lwr_b/(2*1.96);
ldlr_ep_z = ldlr_ep_bc_estimate/ldlr_ep_SE;
ldlr_ep_P = exp(-0.717*ldlr_ep_z - 0.416*ldlr_ep_z^2)
ldlr_ep_P # 

######## LPL ######## 
# this function by default will output this in transpose form
# estimates
lpl_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                        c(-0.162128549293004	, # first half (fixed)
                                          -0.0225771723876219, # second half
                                          -0.0241069985574238 # total 
                                        )))
# Variance covariance matrix 
# SEs
lpl_bc_vcov <- diag(c(0.135993812257091, # first half (fixed)
                      0.0154301328222015, # second half
                      0.0352011753332954) # total 
                      ^2)

## Estimate the standard error
lpl_bc_se_b <- sqrt(t(lpl_bc_grad_g) %*% lpl_bc_vcov %*%lpl_bc_grad_g)
lpl_bc_se_b

# 95% CIs
lpl_bc_estimate = .15
lpl_bc_lwr_b = lpl_bc_estimate - 1.96*lpl_bc_se_b;lpl_bc_lwr_b
lpl_bc_upp_b = lpl_bc_estimate + 1.96*lpl_bc_se_b;lpl_bc_upp_b

#p value
lpl_SE = lpl_bc_upp_b - lpl_bc_lwr_b/(2*1.96);
lpl_z = lpl_bc_estimate/lpl_SE;
lpl_P = exp(-0.717*lpl_z - 0.416*lpl_z^2)
lpl_P # 

######## ER - BC ######## 
# this function by default will output this in transpose form
# estimates
lpl_en_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                           c(-0.162128549293004, # first half (fixed)
                                             -0.022550430586603, # second half
                                             -0.0963156328786269	 # total 
                                           )))
# Variance covariance matrix 
# SEs
lpl_en_bc_vcov <- diag(c(0.135993812257091, # first half (fixed)
                         0.0279869477723851, # second half
                         0.0702400173397047) # total 
                         ^2)

## Estimate the standard error
lpl_en_bc_se_b <- sqrt(t(lpl_en_bc_grad_g) %*% lpl_en_bc_vcov %*%lpl_en_bc_grad_g)
lpl_en_bc_se_b

# 95% CIs
lpl_en_bc_estimate = .04 
lpl_en_bc_lwr_b = lpl_en_bc_estimate - 1.96*lpl_en_bc_se_b;lpl_en_bc_lwr_b
lpl_en_bc_upp_b = lpl_en_bc_estimate + 1.96*lpl_en_bc_se_b;lpl_en_bc_upp_b

#p value
lpl_en_bc_SE = lpl_en_bc_upp_b - lpl_en_bc_lwr_b/(2*1.96);
lpl_en_bc_z = lpl_en_bc_estimate/lpl_en_bc_SE;
lpl_en_bc_P = exp(-0.717*lpl_en_bc_z - 0.416*lpl_en_bc_z^2)
lpl_en_bc_P # 

######## ER + BC ######## 
# this function by default will output this in transpose form
# estimates
lpl_ep_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                           c(-0.162128549293004, # first half (fixed)
                                             -0.0464566847330293, # second half
                                             -0.052063537316867	# total 
                                           )))
# Variance covariance matrix 
# SEs
lpl_ep_bc_vcov <- diag(c(0.135993812257091, # first half (fixed)
                         0.0243321227922063, # second half
                         0.0419743665188365) # total 
                         ^2)

## Estimate the standard error
lpl_ep_bc_se_b <- sqrt(t(lpl_ep_bc_grad_g) %*% lpl_ep_bc_vcov %*%lpl_ep_bc_grad_g)
lpl_ep_bc_se_b

# 95% CIs
lpl_ep_bc_estimate = .14
lpl_ep_bc_lwr_b = lpl_ep_bc_estimate - 1.96*lpl_ep_bc_se_b;lpl_ep_bc_lwr_b
lpl_ep_bc_upp_b = lpl_ep_bc_estimate + 1.96*lpl_ep_bc_se_b;lpl_ep_bc_upp_b


#p value
lpl_ep_SE = lpl_en_bc_upp_b - lpl_en_bc_lwr_b/(2*1.96);
lpl_ep_z = lpl_ep_bc_estimate/lpl_ep_SE;
lpl_ep_P = exp(-0.717*lpl_ep_z - 0.416*lpl_ep_z^2)
lpl_ep_P # 

######## NPC1L1 ######## 
# this function by default will output this in transpose form
# estimates
npc1l1_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                      c(-0.968294080117255	, # first half (fixed)
                                        -0.0225771723876219, # second half
                                        -0.319671261087984# total 
                                      )))
# Variance covariance matrix 
# SEs
npc1l1_bc_vcov <- diag(c(0.506117595887756, # first half (fixed)
                         0.0154301328222015, # second half
                         0.132268107024562) # total 
                    ^2)

## Estimate the standard error
npc1l1_bc_se_b <- sqrt(t(npc1l1_bc_grad_g) %*% npc1l1_bc_vcov %*%npc1l1_bc_grad_g)
npc1l1_bc_se_b

# 95% CIs
npc1l1_bc_estimate = 0.07
npc1l1_bc_lwr_b = npc1l1_bc_estimate - 1.96*npc1l1_bc_se_b;npc1l1_bc_lwr_b
npc1l1_bc_upp_b = npc1l1_bc_estimate + 1.96*npc1l1_bc_se_b;npc1l1_bc_upp_b

#p value
npc1l1_SE = npc1l1_bc_upp_b - npc1l1_bc_lwr_b/(2*1.96);
npc1l1_z = npc1l1_bc_estimate/npc1l1_SE;
npc1l1_P = exp(-0.717*npc1l1_z - 0.416*npc1l1_z^2)
npc1l1_P # 

######## ER - BC ######## 
# this function by default will output this in transpose form
# estimates
npc1l1_en_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                         c(-0.968294080117255, # first half (fixed)
                                           -0.022550430586603, # second half
                                           -0.0574168649266781	 # total 
                                         )))
# Variance covariance matrix 
# SEs
npc1l1_en_bc_vcov <- diag(c(0.506117595887756, # first half (fixed)
                            0.0279869477723851, # second half
                            0.236722556541332) # total 
                       ^2)

## Estimate the standard error
npc1l1_en_bc_se_b <- sqrt(t(npc1l1_en_bc_grad_g) %*% npc1l1_en_bc_vcov %*%npc1l1_en_bc_grad_g)
npc1l1_en_bc_se_b

# 95% CIs
npc1l1_en_bc_estimate = .38
npc1l1_en_bc_lwr_b = npc1l1_en_bc_estimate - 1.96*npc1l1_en_bc_se_b;npc1l1_en_bc_lwr_b
npc1l1_en_bc_upp_b = npc1l1_en_bc_estimate + 1.96*npc1l1_en_bc_se_b;npc1l1_en_bc_upp_b

#p value
npc1l1_en_bc_SE = npc1l1_en_bc_upp_b - npc1l1_en_bc_lwr_b/(2*1.96);
npc1l1_en_bc_z = npc1l1_en_bc_estimate/npc1l1_en_bc_SE;
npc1l1_en_bc_P = exp(-0.717*npc1l1_en_bc_z - 0.416*npc1l1_en_bc_z^2)
npc1l1_en_bc_P # 

######## ER + BC ######## 
# this function by default will output this in transpose form
# estimates
npc1l1_ep_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                         c(-0.968294080117255, # first half (fixed)
                                           -0.0464566847330293, # second half
                                           -0.414455302478785	# total 
                                         )))
# Variance covariance matrix 
# SEs
npc1l1_ep_bc_vcov <- diag(c(0.506117595887756, # first half (fixed)
                            0.0243321227922063, # second half
                            0.155899281271633) # total 
                       ^2)

## Estimate the standard error
npc1l1_ep_bc_se_b <- sqrt(t(npc1l1_ep_bc_grad_g) %*% npc1l1_ep_bc_vcov %*%npc1l1_ep_bc_grad_g)
npc1l1_ep_bc_se_b

# 95% CIs
npc1l1_ep_bc_estimate = .11
npc1l1_ep_bc_lwr_b = npc1l1_ep_bc_estimate - 1.96*npc1l1_ep_bc_se_b;npc1l1_ep_bc_lwr_b
npc1l1_ep_bc_upp_b = npc1l1_ep_bc_estimate + 1.96*npc1l1_ep_bc_se_b;npc1l1_ep_bc_upp_b

#p value
npc1l1_ep_SE = npc1l1_en_bc_upp_b - npc1l1_en_bc_lwr_b/(2*1.96);
npc1l1_ep_z = npc1l1_ep_bc_estimate/npc1l1_ep_SE;
npc1l1_ep_P = exp(-0.717*npc1l1_ep_z - 0.416*npc1l1_ep_z^2)
npc1l1_ep_P # 

######## PCSK9 ######## 
# this function by default will output this in transpose form
# estimates
pcsk9_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                         c(-0.304636488848	, # first half (fixed)
                                           -0.0225771723876219, # second half
                                           -0.121293055125405# total 
                                         )))
# Variance covariance matrix 
# SEs
pcsk9_bc_vcov <- diag(c(0.2784950004837, # first half (fixed)
                        0.0154301328222015, # second half
                        0.0590968038217517) # total 
                       ^2)

## Estimate the standard error
pcsk9_bc_se_b <- sqrt(t(pcsk9_bc_grad_g) %*% pcsk9_bc_vcov %*%pcsk9_bc_grad_g)
pcsk9_bc_se_b

# 95% CIs
pcsk9_bc_estimate = 0.06
pcsk9_bc_lwr_b = pcsk9_bc_estimate - 1.96*pcsk9_bc_se_b;pcsk9_bc_lwr_b
pcsk9_bc_upp_b = pcsk9_bc_estimate + 1.96*pcsk9_bc_se_b;pcsk9_bc_upp_b

#p value
pcsk9_SE = pcsk9_bc_upp_b - pcsk9_bc_lwr_b/(2*1.96);
pcsk9_z = pcsk9_bc_estimate/pcsk9_SE;
pcsk9_P = exp(-0.717*pcsk9_z - 0.416*pcsk9_z^2)
pcsk9_P # 

######## ER - BC ######## 
# this function by default will output this in transpose form
# estimates
pcsk9_en_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                            c(-0.304636488848, # first half (fixed)
                                              -0.022550430586603, # second half
                                              -0.191704363346305	 # total 
                                            )))
# Variance covariance matrix 
# SEs
pcsk9_en_bc_vcov <- diag(c(0.2784950004837, # first half (fixed)
                           0.0279869477723851, # second half
                           0.115170031607877) # total 
                          ^2)

## Estimate the standard error
pcsk9_en_bc_se_b <- sqrt(t(pcsk9_en_bc_grad_g) %*% pcsk9_en_bc_vcov %*%pcsk9_en_bc_grad_g)
pcsk9_en_bc_se_b

# 95% CIs
pcsk9_en_bc_estimate = .04
pcsk9_en_bc_lwr_b = pcsk9_en_bc_estimate - 1.96*pcsk9_en_bc_se_b;pcsk9_en_bc_lwr_b
pcsk9_en_bc_upp_b = pcsk9_en_bc_estimate + 1.96*pcsk9_en_bc_se_b;pcsk9_en_bc_upp_b

#p value
pcsk9_en_bc_SE = pcsk9_en_bc_upp_b - pcsk9_en_bc_lwr_b/(2*1.96);
pcsk9_en_bc_z = pcsk9_en_bc_estimate/pcsk9_en_bc_SE;
pcsk9_en_bc_P = exp(-0.717*pcsk9_en_bc_z - 0.416*pcsk9_en_bc_z^2)
pcsk9_en_bc_P # 

######## ER + BC ######## 
# this function by default will output this in transpose form
# estimates
pcsk9_ep_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                            c(-0.304636488848, # first half (fixed)
                                              -0.0464566847330293, # second half
                                              -0.0790278379781282	# total 
                                            )))
# Variance covariance matrix 
# SEs
pcsk9_ep_bc_vcov <- diag(c(0.2784950004837, # first half (fixed)
                           0.0243321227922063, # second half
                           0.0706078846426028) # total 
                          ^2)

## Estimate the standard error
pcsk9_ep_bc_se_b <- sqrt(t(pcsk9_ep_bc_grad_g) %*% pcsk9_ep_bc_vcov %*%pcsk9_ep_bc_grad_g)
pcsk9_ep_bc_se_b

# 95% CIs
pcsk9_ep_bc_estimate = .18
pcsk9_ep_bc_lwr_b = pcsk9_ep_bc_estimate - 1.96*pcsk9_ep_bc_se_b;pcsk9_ep_bc_lwr_b
pcsk9_ep_bc_upp_b = pcsk9_ep_bc_estimate + 1.96*pcsk9_ep_bc_se_b;pcsk9_ep_bc_upp_b

#p value
pcsk9_ep_SE = pcsk9_en_bc_upp_b - pcsk9_en_bc_lwr_b/(2*1.96);
pcsk9_ep_z = pcsk9_ep_bc_estimate/pcsk9_ep_SE;
pcsk9_ep_P = exp(-0.717*pcsk9_ep_z - 0.416*pcsk9_ep_z^2)
pcsk9_ep_P # 
################### ################### ################### ################### 
################### OTHER CONFIGURATIONS #############################
rm(list = ls()) 
########
ratio <- function(x){
  x[1]*x[2]/x[3]
}
#############################
# cetp -> BC ep (weighted median)
first_e = 0.245492356
first_se = 0.0299094367
second_e = 0.056612421 # second half using weighted median
second_se = 0.019827268
total_e = 0.21414186
total_se = 0.07304125 #*change 9 to 7
# any drug -> BC
any_ep_bc_grad_g <- t(numDeriv::jacobian(func = ratio,
                                          c(first_e, # first half (fixed) *HMGCR
                                            second_e, # second half *LDLR
                                            total_e		# total
                                          )))
# Variance covariance matrix
# SEs
any_ep_bc_vcov <- diag(c(first_se, # first half (fixed)
                         second_se, # second half
                         total_se) # total
                        ^2)

## Estimate the standard error
any_ep_bc_se_b <- sqrt(t(any_ep_bc_grad_g) %*% any_ep_bc_vcov %*%any_ep_bc_grad_g)
any_ep_bc_se_b #0.03269882

# 95% CIs
any_ep_bc_estimate = first_e*second_e/total_e #0.0008108263
any_ep_bc_lwr_b = any_ep_bc_estimate - 1.96*any_ep_bc_se_b;any_ep_bc_lwr_b #0.0008108263
any_ep_bc_upp_b = any_ep_bc_estimate + 1.96*any_ep_bc_se_b;any_ep_bc_upp_b #0.1289902

#p value
any_ep_SE = any_ep_bc_upp_b - any_ep_bc_lwr_b/(2*1.96);
any_ep_z = any_ep_bc_estimate/any_ep_SE;
any_ep_P = exp(-0.717*any_ep_z - 0.416*any_ep_z^2)
any_ep_P # 0.6268907

######## APOC3 ######## 
rm(list = ls())
#
ratio <- function(x){
  x[1]*x[2]/x[3]
}

# read data
drug_exp_data <- read_exposure_data("~/Desktop/Skeleton_projects/oxysterol/analysis/drug_exposure_data_carter_and_gormely_et_al.csv", sep=",",clump = FALSE) # it includes triglycerides
out_dat_comb <- import("~/Desktop/Skeleton_projects/oxysterol/analysis/drug_cancer/out_dat_comb.csv")
# filter outcome BC and exposure
out_dat_comb_bc<- filter(out_dat_comb, id.outcome == "ieu-a-1126")
dat_1 <- harmonise_data(drug_exp_data, out_dat_comb_bc, action = 2)
dat_te<- filter(dat_1, exposure == "APOC3")
dat_total_MR<- mr(dat_te)
dat_total_MR<- filter(dat_total_MR, method == "Inverse variance weighted")
TE<- abs(dat_total_MR$b)
TE_se<- abs(dat_total_MR$se)
#
out_dat_comb_lxr <- extract_outcome_data(snps = drug_exp_data$SNP, outcomes = "prot-a-2087", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat_2 <- harmonise_data(drug_exp_data, out_dat_comb_lxr, action = 2)
dat_first<- filter(dat_2, exposure == "APOC3")
dat_first_MR<- mr(dat_first)
First_e<- abs(dat_first_MR$b)
First_se<- abs(dat_first_MR$se)
#
exposure_data <- extract_instruments("prot-a-2087", p1 = 5e-08)
out_dat_comb_lxr_bc <- extract_outcome_data( snps = exposure_data$SNP,outcomes = "ieu-a-1126", proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3)
dat_3 <- harmonise_data(exposure_data, out_dat_comb_lxr_bc, action = 2)
dat_second_MR<- mr(dat_3)
dat_second_MR<- filter(dat_second_MR, method == "Inverse variance weighted")
Second_e<- abs(dat_second_MR$b)
Second_se<- abs(dat_second_MR$se)

# estimates
apoc3_bc_grad_g <- t(numDeriv::jacobian(func = ratio, 
                                        c(First_e, # first half (fixed)
                                          Second_e, # second half
                                          TE # total 
                                        )))
# Variance covariance matrix 
# SEs
apoc3_bc_vcov <- diag(c(First_se, # first half (fixed)
                        Second_se, # second half
                        TE_se) # total 
                      ^2)

## Estimate the standard error
apoc3_bc_se_b <- sqrt(t(apoc3_bc_grad_g) %*% apoc3_bc_vcov %*%apoc3_bc_grad_g)
apoc3_bc_se_b # 0.4867821

# 95% CIs
apoc3_bc_estimate = First_e*Second_e/TE
apoc3_bc_estimate # 0.3334857
apoc3_bc_lwr_b = apoc3_bc_estimate - 1.96*apoc3_bc_se_b;apoc3_bc_lwr_b #-0.6206073
apoc3_bc_upp_b = apoc3_bc_estimate + 1.96*apoc3_bc_se_b;apoc3_bc_upp_b #1.287579

#p value
apoc3_SE = apoc3_bc_upp_b - apoc3_bc_lwr_b/(2*1.96);
apoc3_z = apoc3_bc_estimate/apoc3_SE;
apoc3_P = exp(-0.717*apoc3_z - 0.416*apoc3_z^2)
apoc3_P #0.8290287

################################################
# # EXAMPLE FROM PAPER YOSHIJI ET AL. NAT MET 2023. see supplementary 
# cetp_ep_bc_grad_g <- t(numDeriv::jacobian(func = ratio,
#                                           c(0.145009643, # first half (fixed)
#                                             0.5376731, # second half
#                                             0.5615565		# total
#                                           )))
# # Variance covariance matrix
# # SEs
# cetp_ep_bc_vcov <- diag(c(0.031057256, # first half (fixed)
#                           0.0841015, # second half
#                           0.05872464) # total
#                         ^2)
# 
# ## Estimate the standard error
# cetp_ep_bc_se_b <- sqrt(t(cetp_ep_bc_grad_g) %*% cetp_ep_bc_vcov %*%cetp_ep_bc_grad_g)
# cetp_ep_bc_se_b
# 
# # 95% CIs
# cetp_ep_bc_estimate = 0.1388423
# cetp_ep_bc_lwr_b = cetp_ep_bc_estimate - 1.96*cetp_ep_bc_se_b;cetp_ep_bc_lwr_b
# cetp_ep_bc_upp_b = cetp_ep_bc_estimate + 1.96*cetp_ep_bc_se_b;cetp_ep_bc_upp_b
# 
# 
# cetp_ep_bc_grad_g <- t(numDeriv::jacobian(func = ratio,
#                                           c(0.145009643, # first half (fixed)
#                                             0.5376731, # second half
#                                             0.5615565		# total
#                                           )))
# # Variance covariance matrix
# # SEs
# cetp_ep_bc_vcov <- diag(c(0.031057256, # first half (fixed)
#                           0.0841015, # second half
#                           0.05872464) # total
#                         ^2)
# 
# ## Estimate the standard error
# cetp_ep_bc_se_b <- sqrt(t(cetp_ep_bc_grad_g) %*% cetp_ep_bc_vcov %*%cetp_ep_bc_grad_g)
# cetp_ep_bc_se_b
# 
# # 95% CIs
# cetp_ep_bc_estimate = 0.1388423
# cetp_ep_bc_lwr_b = cetp_ep_bc_estimate - 1.96*cetp_ep_bc_se_b;cetp_ep_bc_lwr_b
# cetp_ep_bc_upp_b = cetp_ep_bc_estimate + 1.96*cetp_ep_bc_se_b;cetp_ep_bc_upp_b
# 
# # proportion mediated
# library(dplyr)
# # https://www.alexstephenson.me/post/2022-04-02-standard-errors-and-the-delta-method/
# 




