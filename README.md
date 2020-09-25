# bypass-length-sar
This repository holds the data and R code used in analyses to support Faulkner et al. (2019) and Faulkner et al. (2020).

## Data
There are two data sets used for the analyses on bypass probabilities and there are four data sets used for the analyses on adult return probabilities (smolt-to-adult returns or SAR).  These data were used in both the main paper and the formal response to comment on the paper.
* `chinook_bypass_data.csv` - bypass event data and covariates for Chinook
* `steelhead_bypass_data.csv` - bypass event data and covariates for steelhead
* `chinook_ULGR_SAR_data.csv` - adult return data and covariates for Chinook tagged upstream of Lower Granite Dam
* `chinook_LGR_SAR_data.csv` - adult return data and covariates for Chinook tagged at Lower Granite Dam
* `steelhead_ULGR_SAR_data.csv` - adult return data and covariates for steelhead tagged upstream of Lower Granite Dam
* `steelhead_LGR_SAR_data.csv` - adult return data and covariates for steelhead tagged at Lower Granite Dam

## R Code
The following scripts were used for analyses in the main paper:
* `chinook_bypass_analysis.r` - model fitting and selection for bypass probabilities for Chinook
* `steelhead_bypass_analysis.r` - model fitting and selection for bypass probabilities for steelhead
* `chinook_sar_analysis.r` - model fitting and selection and simulations for adult return probabilities for Chinook
* `steelhead_sar_analysis.r` - model fitting and selection and simulations for adult return probabilities for steelhead

The following scripts were used for analyses and plots in the response to comment:
* `chinook_power_simulations.r` - simulations to produce power estimates for Chinook used in Figure 2 and Table 1
* `steelhead_power_simulations.r`  - simulations to produce power estimates for steelhead used in Figure 2 and Table 1
* `adult_return_pred_for_response_fig1.r` - code to generate adult return probability predictions in Figure 1
* `chinook_lgs_bypass_response_fig3.r` - code to generate bypass probability predictions in Figure 3 

## References
Faulkner, JR, BL Bellerud, DL Widener, and RW Zabel. 2019. Associations among fish length, dam passage history, and survival to adulthood in two at-risk species of Pacific salmon. Transactions of the American Fisheries Society 148:1069-1087.

Faulkner, JR, BL Bellerud, DL Widener, SG Smith and RW Zabel. Associations among fish length, dam passage history, and survival to adulthood in two at-risk species of Pacific salmon: response to comment. Transactions of the American Fisheries Society -- In Review


