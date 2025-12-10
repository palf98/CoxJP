# CoxJP
Example code accompanying the publication "Spanish Paediatric Haematology and Oncology survival results and trends, 1999-2022"

The provided code is an implementation of the Cox proportional hazards model with joinpoints described in Mariotto et al., 2021 (Mariotto AB, Zhang F, Buckman DW, Miller D, Cho H, Feuer EJ. Characterizing Trends in Cancer Patients' Survival Using the JPSurv Software. Cancer Epidemiol Biomarkers Prev. 2021;30(11):2001-2009. doi: 10.1158/1055-9965.EPI-21-0423.). The original implementation of this model (JPSurv), is available in https://survivalstatstools.cancer.gov/jpsurv .

Our implementation is aimed at analysing long-term time trends of survival with up to 2 joinpoints, since if facilitates users to perform age-standardisation within the functions themselves. Functions are not generalised for all scenarios and assume that, given a time period of analysis, last years can be covered with observed period estimates of survival to represent graphically fitted trends over observed values.

In the provided code, users can find both all the necessary code to select the model with best fit for a dataset and determine the statistical signification of trends through a bootstrap resampling algorithm of the main trend metric used in Cox-joinpoint models: annual absolute change in survival (AACS). In Line 286, users can find an example dataset, which is completely simulated, to try the functions. 
