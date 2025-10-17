library(dplyr)
library(geepack)

# Choose home directory where the weighted & replicated example data is stored
homedir = "C:/Users/username/homedir"
setwd(homedir)
df = read.csv("Sparsity_DataExample.csv")

##  Note: This file was generated using the following settings
# N = 500
# pNP_target=1/3 #  desired proportion of individuals expressing No Preference in stage 1
# pTheta_target=0.4 # desired proportion of individuals expressing preference for treatment A among those with a preference in stage 1
# pNP2_target=1/3 # desired proportion of patients expressing No Preference in stage 2 (among non-responders)
# pTheta2_target=0.9 # desired proportion of individuals expressing preference for treatment C among those with a preference in stage 2 (among non-responders)
# Pa=0.6
# Pa1=0.6
# Pb=0.4
# Pb1=0.4
# Other settings in type1_small.R

## Check pathways
table(df$Trial_Path)
# 2 out of 20 total pathways are empty
# RANRPD and PANRPD: non-responders to A that switch to D with a stage 2 preference

## Unrestricted Analysis -- Estimate all 16 DTRs
# Even though two pathways are empty, we can still estimate IPT weights and conduct a standard PRPP-SMART analysis
contrast_dtr_1 <- matrix(c(1,1,0,1,0,1, #AC00
                           1,1,0,0,0,0, #AD00
                           1,0,0,1,0,0, #BC00
                           1,0,0,0,0,0, #BD00
                           1,1,0,1,1,1, #AC01
                           1,1,0,0,1,0, #AD01
                           1,0,0,1,1,0, #BC01
                           1,0,0,0,1,0, #BD01
                           1,1,1,1,0,1, #AC10
                           1,1,1,0,0,0, #AD10
                           1,0,1,1,0,0, #BC10
                           1,0,1,0,0,0, #BD10
                           1,1,1,1,1,1, #AC11
                           1,1,1,0,1,0, #AD11
                           1,0,1,1,1,0, #BC11
                           1,0,1,0,1,0 #BD11
),nrow = 16, ncol=6, byrow = TRUE)
contrast_dtr_2 <- matrix(c(1,1,1,1,0,0,0,0,0, #AC00
                           1,1,0,0,0,0,0,0,0, #AD00
                           1,0,1,0,0,0,0,0,0, #BC00
                           1,0,0,0,0,0,0,0,0, #BD00
                           1,1,1,1,0,0,1,0,0, #AC01
                           1,1,0,0,0,0,0,1,0, #AD01
                           1,0,1,0,0,0,1,0,0, #BC01
                           1,0,0,0,0,0,0,1,0, #BD01
                           1,1,1,1,1,0,0,0,0, #AC10
                           1,1,0,0,1,0,0,0,0, #AD10
                           1,0,1,0,0,1,0,0,0, #BC10
                           1,0,0,0,0,1,0,0,0, #BD10
                           1,1,1,1,1,0,1,0,1, #AC11
                           1,1,0,0,1,0,0,1,1, #AD11
                           1,0,1,0,0,1,1,0,1, #BC11
                           1,0,0,0,0,1,0,1,1 #BD11
),nrow = 16, ncol=9, byrow = TRUE)
model1 <- geeglm(Y ~ T1 + S1_Preference + T2 + S2_Preference + T1:T2, 
                   id=id, weights=w, family = gaussian, corstr = "independence", data = df)
DTR_hat_1 <- round(contrast_dtr_1 %*% model1$coefficients,2)
model2 <- geeglm(Y ~ T1 + T2 + T1:T2 + S1_Preference:T1_A + S1_Preference:T1_B + S2_Preference:T2_C + S2_Preference:T2_D + S1_Preference:S2_Preference, 
                   id=id, weights=w, family = gaussian, corstr = "independence", data = df)
DTR_hat_2 <- round(contrast_dtr_2 %*% model2$coefficients,2)
DTR_hat = data.frame(DTR_hat_1, DTR_hat_2)
rownames(DTR_hat) = c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "AAD01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "AAD11", "BBC11", "BBD11")
colnames(DTR_hat) = c("Model 1 Estimates", "Model 2 Estimates")

## Restricted Analysis -- Estimate 14 DTRs (excluding AAD01 and AAD11)
# The unrestricted analysis models converge in this case, but it doesn't make sense to estimate AAD01 and AAD11 
# if we don't see any non-responders to D with a stage 2 preference.
# The restricted analysis uses the same weights and same model as the unrestricted analysis, except we don't
# use data for AD non-responders and we don't consider linear combinations of parameters to estimate AAD01 and AAD11.
contrast_dtr_1_mod <- matrix(c(1,1,0,1,0,1, #AC00
                               1,1,0,0,0,0, #AD00
                               1,0,0,1,0,0, #BC00
                               1,0,0,0,0,0, #BD00
                               1,1,0,1,1,1, #AC01
                               1,0,0,1,1,0, #BC01
                               1,0,0,0,1,0, #BD01
                               1,1,1,1,0,1, #AC10
                               1,1,1,0,0,0, #AD10
                               1,0,1,1,0,0, #BC10
                               1,0,1,0,0,0, #BD10
                               1,1,1,1,1,1, #AC11
                               1,0,1,1,1,0, #BC11
                               1,0,1,0,1,0 #BD11
),nrow = 14, ncol=6, byrow = TRUE)
contrast_dtr_2_mod <- matrix(c(1,1,1,1,0,0,0,0,0, #AC00
                               1,1,0,0,0,0,0,0,0, #AD00
                               1,0,1,0,0,0,0,0,0, #BC00
                               1,0,0,0,0,0,0,0,0, #BD00
                               1,1,1,1,0,0,1,0,0, #AC01
                               1,0,1,0,0,0,1,0,0, #BC01
                               1,0,0,0,0,0,0,1,0, #BD01
                               1,1,1,1,1,0,0,0,0, #AC10
                               1,1,0,0,1,0,0,0,0, #AD10
                               1,0,1,0,0,1,0,0,0, #BC10
                               1,0,0,0,0,1,0,0,0, #BD10
                               1,1,1,1,1,0,1,0,1, #AC11
                               1,0,1,0,0,1,1,0,1, #BC11
                               1,0,0,0,0,1,0,1,1 #BD11
),nrow = 14, ncol=9, byrow = TRUE)
df2 <- df[!((df$T1 == 1 & df$S1_Preference == 0 & df$T2 == 0 & df$S2_Preference == 1) |
               (df$T1 == 1 & df$S1_Preference == 1 & df$T2 == 0 & df$S2_Preference == 1)),]
model1_mod <- geeglm(Y ~ T1 + S1_Preference + T2 + S2_Preference + T1:T2, 
                     id=id, weights=w, family = gaussian, corstr = "independence", data = df2)
DTR_hat_1_mod <- round(contrast_dtr_1_mod %*% model1_mod$coefficients,2)
model2_mod <- geeglm(Y ~ T1 + T2 + T1:T2 + S1_Preference:T1_A + S1_Preference:T1_B + S2_Preference:T2_C + S2_Preference:T2_D + S1_Preference:S2_Preference, 
                     id=id, weights=w, family = gaussian, corstr = "independence", data = df2)
DTR_hat_2_mod <- round(contrast_dtr_2_mod %*% model2_mod$coefficients,2)
DTR_hat_mod = data.frame(DTR_hat_1_mod, DTR_hat_2_mod)
rownames(DTR_hat_mod) = c("AAC00", "AAD00", "BBC00", "BBD00", "AAC01", "BBC01", "BBD01", "AAC10", "AAD10", "BBC10", "BBD10", "AAC11", "BBC11", "BBD11")
colnames(DTR_hat_mod) = c("Model 1 Estimates", "Model 2 Estimates")
# DTR estimates are similar but not identical to the unrestricted analysis since we are no longer using the
# replicates from responders to A which were consistent with AAD01 and AAD11
