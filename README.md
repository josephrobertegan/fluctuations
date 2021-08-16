# fluctuations
Data and MATLAB code for reproducing the figures in the 'fluctuations' paper.

This repository relates to the paper entitled "Fluctuations in TCR and pMHC interactions regulate T cell activation" available on bioRxiv at the following link:

https://www.biorxiv.org/content/10.1101/2021.02.09.430441v2

The MATLAB files entitled "scrp_stoc_sim.m", "scrp_dose_resp_ther.m" and "scrp_dose_resp_expr.m" are scripts that reproduce figures 2, 3 and 4, respectively. Note that the first and the last of these scripts are very quick to run (of the order of seconds) whereas the second script takes much longer to run (of the order of a few hours) due to the repeated stochastic simulations that are undertaken.

The "data_andersen_2001.csv", "data_chmielewski_2004.csv" and "data_mcmahan_2006.csv" files contain data from three published papers that are reproduced in Figure 4a, 4b and 4c, respectively.  The "data_abushah_2021.xlsx" file contains previously unpublished data as shown in Figure 4d. 

In addition, the MATLAB file entitled "scrp_exam.m" is a script that includes examples of how to use the various other MATLAB functions that are provided within this repository. These functions provide computational calculations for the important equations in the paper.
