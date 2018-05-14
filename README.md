
# SAEM for logistic regression with missing values
Here are R codes of functions and implementations for paper
"Stochastic Approximation EM for Logistic regression with missing values (2018, Jiang W., Josse J., Lavielle M.)".

## Installation of package 
First you can install the package **misaem** from Github. The main function `miss.saem` contains the procedure of estimation for parameters, as well as their variance, and observed likelihood.
```{r}
library(devtools)
install_github("wjiang94/misaem")
 ```
## Implementation of simulation study

When you execute *Rmarkdown* codes in the following files, you can reproduce the results of the paper. You can also refer to the corresponding *html* file to see the implementations and comments. 

**convergence_saem**: 
* Brief introduction of simulation procedure
* Demonstration of convergence of SAEM (Figure 1 and Figure 8)

**simu_mcar_saem**: 
Simulation study to assess the performance of SAEM by comparing several other existing methods for missing data :
* The complete case (CC) method : all rows containing at least one unobserved data value were removed)
* Multiple imputation based on conditional modelling as implemented in the R package `mice` (with its default settings and Rubin's combining rules)
* MCEM algorithm that we implemented using adaptive rejection sampling (MCEM-AR). 
* We use the dataset without missing values (no NA) as a reference, with parameters estimated with the Newton-Raphson algorithm as implemented in the `glm` function in R.

We run repetitions of simulations. And evaluate their performance, intially in terms of estimation errors of the parameters, as well as the standard error of estimation and the coverage of confidence interval. Here we can reproduce
* Bias for estimation with 100 times of simulations (Figure 2).
* Estimated standard error with 100 times of simulations (Figure 2).
* Coverage of confidence interval with 1000 times of simulations (Table 1).
* Execution time (Table 2).

**model_selection**:  Model selection results (Table 3).
