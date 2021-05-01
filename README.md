Example data and code for the manuscript "Flexible fitting of PROTAC concentration-response curves with changepoint Gaussian Processes": https://www.biorxiv.org/content/10.1101/2020.11.13.379883v1

The kernel implemented in the Stan file, is a changepoint kernel with a steep transition function (see Supplement).

# Ways to run this code:

## 1. On a local machine

File 'analysis.R' provides code using the 'rstan' library, file 'analysis_cmdstanr.R' provides code using the 'cmdstanr' library.

## 2. Using Kaggle R kernel

To avoid installing the "rstan" library, one can run the model in the cloud using Kaggle R kernel, following the pre-made notebook. 

Here are the steps to follow:

- make sure to have a Kaggle account in order to be able to use Kaggle's R kernel
- follow the link to open notebook https://www.kaggle.com/lizasemenova/gp-concentration-response
- press 'copy and edit' (top right corner)
- download this GitHub repository

Data, the Stan model and a file with helper functions need to be added to Kaggle manually as explained below:

### Add data:

- make sure that the 'datagp' folder, containing 'compound.csv' and 'controls.scv', has been downloaded from GitHub

- on Kaggle, click 'File' (top left corner) -> 'Add or upload data' -> 'Ulpoad dataset' -> enter dataset title 'datagp' -> browse your files and upload 'compound.csv' and 'controls.scv' -> click 'Create'

- it might take some time for Kaggle to upload the files

- click 'File' (left top corner) -> add or upload data -> click 'Your datasets' -> chose 'datagp' -> click 'Add'

- make sure that 'datagp' folder has appeared in the 'Data' tab (at the top right corner) within the folder 'input'



### Add Stan code and helper functions:

- make sure that the 'modelgp' folder, containing 'functions.R' and 'model_invlogit_nu_g.stan', has been downloaded from GitHub

- on Kaggle, click 'File' (top left corner) -> 'Add or upload data' -> 'Upload dataset' -> enter dataset title 'modelgp' -> browse your files and upload 'functions.R' and 'model_invlogit_nu_g.stan' -> click 'Create'

- it might take some time for Kaggle to upload the files

- click 'File' (left top corner) -> add or upload data -> click 'Your datasets' -> chose 'modelgp' -> click 'Add'

- make sure that 'modelgp' folder has appeared in the 'Data' tab (at the top right corner) within the folder 'input'



- Check that both folders have been added to the right place by running `list.files(path = "../input")`. The output should be `'datagp''modelgp'`

