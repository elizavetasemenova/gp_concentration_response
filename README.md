Example data and code for the manuscript "Flexible fitting of PROTAC concentration-response curves with changepoint Gaussian Processes": https://www.biorxiv.org/content/10.1101/2020.11.13.379883v1

The kernel implemented in the Stan file, is a changepoint kernel with a steep transition function (see Supplement).

To avoid installing rstan, one can run the model in the cloud using Kaggle R kernel: https://www.kaggle.com/lizasemenova/gp-concentration-response. One needs to create a Kaggle account in order to be able to use Kaggle's R kernel.
