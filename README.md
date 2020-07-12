# Bayesian Parameter Estimation Approach using MCMC for Model Identification of Glucose-Insulin Dynamics for Type 1 Diabetes

Bayesian parameter estimation approach using MCMC technique is designed to identify the model of glucose-insulin dynamics
using both experimental and simulation data. Experimental data from a real patient with type 1 diabetes using subcutaneous pump therapy who participated
in a clinical research study is used. The key quantity in Bayesian parameter estimation is the a posteriori probability density function of the
unknown model parameters. From this function, which provides a complete description of the shape of the estimate
uncertainty, the 95% confidence intervals can be derived. Obtaining the a posteriori probability density function of
the unknown model parameters is a task that is analytically intractable because of the complex relationships between
parameters and data. Therefore, MCMC was used to obtain the a posteriori probability density function. In overall,
for the MCMC 6000 iterations were run. To allow the MCMC to stabilize, the first 1000 iterations were dropped out. Assessment of model fit is investigated by calculating the variance accounted for (VAF) in the all implemented
identification scenarios. Posterior distributions of parameter estimates obtained from MCMC in the experimental study are shown. The posterior distribution contains information about the mean and median as well as associated uncertainty.
Measures of uncertainty such as the 95% credible intervals can be extracted from the posterior distribution. Based on the provided results and
performance measures, the effectiveness of Bayesian approach on both simulation and experimental scenarios was
seen.
