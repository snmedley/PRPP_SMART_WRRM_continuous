# PRPP_SMART_WRRM_continuous
Data generation and frequentist analytic methods for a Partially Randomized, Patient Preference, Sequential, Multiple Assignment, Randomized Trial (PRPP-SMART) with a continuous end-of-trial outcome. 

## The PRPP-SMART Design
Our PRPP-SMART design assumes there are two stages with two treatment options per stage. At the beginning of stage 1, participants are asked if they have a preference between the two stage 1 treament options (A, B). All participants with a preference are assigned to their preferred treatment while all others are randomized to one of the two treatment options. At the end of stage 1, response status is determined (i.e., responder or non-responder). Responders continue their stage 1 treatment in stage 2 while non-responders are re-assigned to two new treatment options that may be more effective (C, D). Treatment preference is elicited from non-responders at the beginning of stage 2, and non-responders with a preference receive their preferred treatment while all other non-responders are randomized. 

[add design figure here]

## Data Generation
xxx.R contains R code to generate data consistent with our PRPP-SMART design. The data generation procedure was adapted from Wang et al. (2022) by co-author Mari Wank for PRPP-SMARTs with a binary end-of-trial outcome (see https://github.com/mariwank/PRPP_SMART_WRRM). The outcome generation was modified for continuous outcomes, but the rest of the data generation procedure remains the same. 

## Analytic Methods
Notation:
- A_1 = {A, B}, A_2 = {C, D} denote stage 1 and 2 treatment, respectively
- T_1 and T_2 are binary indicators of stage 1 and 2 treatment, respectively with T_1 = 1 for A and 0 for B, T_2 = 1 for C and 0 for D
- P_1 and P_2 are binary indicators of stage 1 and 2 preference, respectively

Our goal is to estimate dynamic treatment regimens (DTRs) subject to treatment preference. We denote DTRs by [A_1A_1A_2]_{P_1P_2}, and there are 16 DTRs embedded in our PRPP-SMART design (i.e., the 4 traditional DTRs AAC, AAD, BBC, and BBD but within each preference combination). We extend frequentist weighted and replicated regression models (WRRMs) for traditional SMARTs to PRPP-SMARTs. We consider a reduced mean model (Model 1) and a more flexible mean model (Model 2):

Model 1: E[[A_1A_1A_2]_{P_1P_2} | T_1, P_1, T_2, P_2] = alpha_1 + beta_1*T_1 + theta_1*T_2 + gamma_1*T_1*T_2 + delta_1*T_1 + delta_2*T_2
Model 2: E[[A_1A_1A_2]_{P_1P_2} | T_1, P_1, T_2, P_2] = alpha_2 + beta_2*T_1 + theta_2*T_2 + gamma_2*T_1*T_2 + omega_B*T_1 + omega_A*T_1*P_1 + omega_D*T_2 + omega_C*T_2*P_2
