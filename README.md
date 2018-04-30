
# Stochastic Approximated EM for Logistic regression with missing data
Codes and implementations for "Stochastic Approximated EM for Logistic regression with missing data".

**saem.R**: Functions for parameters' estimation with SAEM, estimation of their vairiance and observed log-likelihood.

**convergene_saem.R**: Implementation of convergence of SAEM (Figure 1)

**simu_mcar_saem.R**: Simulation study on comparison of complete case analysis method (CC), multiple imputation (mice), SAEM and the classical estimation procedure on the original dataset without missing data (no NA). The missing mechanism is Missing completely at random (MCAR).
1. Bias for estimation with 100 times of simulations (Figure 2).
2. Estimated standard error with 100 times of simulations (Figure 3).
3. Coverage of confidence interval with1000 times of simulations (Table 1).
4. Execution time (Table 2).

**model_selection.R**:  Model selection results (Table 3).
