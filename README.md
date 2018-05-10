
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

**simu_mcar_saem.R**: Simulation study on comparison of complete case analysis method (CC), multiple imputation (mice), SAEM and the classical estimation procedure on the original dataset without missing data (no NA). The missing mechanism is Missing completely at random (MCAR).
1. Bias for estimation with 100 times of simulations (Figure 2).
2. Estimated standard error with 100 times of simulations (Figure 3).
3. Coverage of confidence interval with 1000 times of simulations (Table 1).
4. Execution time (Table 2).

**model_selection.R**:  Model selection results (Table 3).
