# Online Affinity Matrix Calibration (OAMC)

Supplementary materials for "Online Affinity Matrix Calibration"

(MLJ special issue'2022 submission)

The code has been tested on MATLAB R2016b and R2020a. It should be able to run on other recent versions.

## Main files:

- example_main.m: demo of missing data processing in Section 5.3
- example_scale.m: demo of scalable extention in Section 5.4
- example_predict.m: demo of prediction task in Section 6.3

## Baseline_MVI (missing value imputation):

- impute_zero.m: ZERO Imputation 
- impute_naive.m: Mean/Median/Mode Imputation
- impute_mean.m: MEAN Imputation
- impute_knn.m: kNN Imputation
- impute_lr.m: Linear Regression-based Imputation
- impute_grouse.m: GROUSE Imputation
- impute_kfmc.m: KFMC Imputation

## Baseline_MC (matrix calibration):

- calibrate_dmc.m: DMC Calibration
- calibrate_cmc.m: CMC Calibration
- calibrate_oamc.m: OAMC Calibration
- calibrate_oamc_block.m: Block OAMC Calibration
- calibrate_soamc.m: Scalable OAMC Calibration

## Other files:

- data_demo.mat: data file for example_main.m and example_scale.m
- data_predict.mat: data file for example_predict.m
- jaccard_approx.m: approximate Jaccard index matrix with incomplete samples
- predict_smoke.m: measure the prediction accuracy


March 30th, 2022
