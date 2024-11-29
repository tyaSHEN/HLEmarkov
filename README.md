# Modeling Disability-free Life Expectancy with Duration Dependence: A Research Note on the bias in the Markov Assumption 

***Tianyu Shen, James O'Donnell***

This is a repository for our paper in [*Demography*](https://doi.org/10.1215/00703370-11058373).

Details of the mathematical derivations and the demographic interpretation of the results are included in our paper.

The code is available in "**Modeling.R**".

It includes all the model by loading cleaned data in "output/Sex"

"Data processing.R" contains code to clean up the HRS data for the purpose of the analysis and the two imputation procedures. The cleaned data would be stored in “output/Sex”. There are 10 iterations of the imputation. These file include two duration columns to track the known and unknown duration. No need to run this script as cleaned data can be found in “output/Sex”, but if you want to you will need to download the HRS data first.

For other enquiries related to the paper please email tianyu.shen@anu.edu.au
