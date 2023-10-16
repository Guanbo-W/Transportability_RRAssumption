# Transportability_RR

Simulation code for "Transportability under a weaker transportability assumption".

_DR_Sim.R_ is the code for comparing $\widehat \alpha_1$ with $\widehat \alpha_1^G$ and $\widehat \alpha_1^W$.

_Comp_Conv.R_ is the code for comparing $\widehat \alpha_1$ with and $\widehat \alpha'_1$.

_Efficiency_RR.R_ is the code for comparing the biases, standard deviations, coverages and RMSE of $\widehat \alpha_1$.

_sim-rates.R_ is the code for verifying the rate robustness property of the estimator. 

Across all the settings, the simulations are run 5000 times, with sample sizes: $n_0=5000, n_1=500/1000/5000$, and ${\mu_{1, 1}(X)\}$ or $\{\mu_{0, 1}(X), \tau(X)\}$ can be correctly or incorrectly specified.
