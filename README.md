Example data and code for the manuscript "Flexible fitting of PROTAC concentration-response curves with changepoint Gaussian Processes": https://www.biorxiv.org/content/10.1101/2020.11.13.379883v1

The kernel implemented in the Stan file, is a changepoint kernel with a steep transition function (see Supplement).


## Using Kaggle R kernel

To avoid installing the "rstan" library, one can run the model in the cloud using Kaggle R kernel, following this notebook: https://www.kaggle.com/lizasemenova/gp-concentration-response. One needs to create a Kaggle account in order to be able to use Kaggle's R kernel.

Data, the Stan model and a file with helper functions need to be added manually as explained below:

### Add data:

- download 'data' folder, containing 'compound.csv' and 'controls.scv', from this GitHub repository

- click 'File' (top left corner) -> 'Add or upload data' -> 'Ulpoad dataset' -> enter dataset title 'data' -> browse your files and upload 'compound.csv' and 'controls.scv' -> click 'Create'

- it might take some time for Kaggle to upload the files

- click 'File' (left top corner) -> add or upload data -> click 'Your datasets' -> chose 'data' -> click 'Add'

- make sure that 'data' folder has appeared in the 'Data' tab (at the top right corner) within the folder 'input'



### Add Stan code and helper functions:

- download the 'model' folder, containing 'functions.R' and 'model_invlogit.stan', from the GitHub repository

- click 'File' (top left corner) -> 'Add or upload data' -> 'Upload dataset' -> enter dataset title 'model' -> browse your files and upload 'functions.R' and 'model_invlogit.stan' -> click 'Create'

- it might take some time for Kaggle to upload the files

- click 'File' (left top corner) -> add or upload data -> click 'Your datasets' -> chose 'model' -> click 'Add'

- make sure that 'model' folder has appeared in the 'Data' tab (at the top right corner) within the folder 'input'


