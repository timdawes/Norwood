# Multi-state modelling for patients after the Norwood Procedure
This is a code repository for the code used in the submitted paper "Predicting clinical trajectory after the Norwood procedure: focus on aortic and pulmonary reinterventions"
Authors: Haapanen H, Dawes TJW, Brown K, Giardini A, Dedieu N, Shetty P, Tsang V, Kostolny M.

Survival analysis after the Norwood procedure is difficult as great-vessel reinterventions are common and influence survival, but vary between patients. We apply contemporary statistical techniques to modelling these complex pathways, to identify prognostic markers, predict long-term outcomes and drive personalised care. We conducted a retrospective survival analysis of patients undergoing the Norwood procedure at a tertiary centre (October 2001 - January 2021) dividing the patients into three, equal-length eras. We developed a multi-state survival model using demographic, surgical, physiological and follow-up data to identify prognostic markers and predict transplant-free survival.

## Software and hardward requirements
All statistical analysis was performed using R version 4.0.2 (R Foundation for Statistical Computing, Vieena, Austria; www.r-project.org) using RStudio version 1.3.1073 (Boston, Mass).

## Files
### 1. Functions.R
#### Loading packages and hand-written functions used in the analysis
### 2. Setup.R
#### Re-configures the Excel sheet data into suitable formats for the analysis
### 3. Analysis.R
#### Analyses the data using multi-state modelling




