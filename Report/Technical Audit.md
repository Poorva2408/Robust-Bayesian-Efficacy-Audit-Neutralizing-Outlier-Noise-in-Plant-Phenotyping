# Technical Audit: PlantGrowth Robust Bayesian Architecture

This repository documents the transition from legacy frequentist ANOVA to a high-fidelity Robust Bayesian Hierarchical Architecture. The audit addresses the failure of standard Gaussian models in the presence of dynamical "Black Swan" interference.

## Explanatory Data Analysis :
To isolate statistical architecture validation from agronomic latency, physical data collection was bypassed. The foundational PlantGrowth dataset was utilized as a calibrated baseline and synthetically expanded to N=42 (stress-testing) and N=125 (power optimization). This controlled scaling allowed for the deliberate injection of dynamical "Black Swan" outliers (1.2, 8.8, 9.5) and the direct benchmarking of the MCMC samplers without the confounding variables and time-costs of physical field testing. By engineering the input vectors rather than growing physical samples, the audit directly targets the resilience of the Robust-t likelihood function and subsequent power scaling strategies.

## Executive Summary 
This technical audit evaluates the transition from legacy frequentist ANOVA to a Robust Bayesian Hierarchical Architecture for analyzing dynamical growth data. Traditional ANOVA and Pooled Bayesian models ($M_1$) failed to isolate treatment effects due to "Black Swan" outlier interference (values of 1.2, 8.8, 9.5), which smeared noise across groups and produced a failing Posterior Predictive Check (PPC) p-value of 0.0475.

While a Separate Variance Gaussian model ($M_2$) identified heteroscedasticity, it remained fragile; the outliers inflated the noise floor to $\sigma=1.29$. By deploying a Robust-t Likelihood ($\nu=4$), we achieved a 35% increase in signal precision, recovering the true noise floor of $\sigma=0.90$. This model identifies Treatment 2 as a high-yield success with 99.8% Bayesian certainty.

To mitigate the 57% risk of false negatives inherent in the current $N=42$ sample, the audit identifies an Efficiency Peak at $N=125$. Scaling to this target provides 80% statistical power, ensuring maximum decision reliability and a defensible ROI for Phase III deployment.

<p align="center">
  <img src="/Figures/EDA_raincloud.png" width="600" alt="Resized Chart Description">
  <br>
  <b>Figure 1:</b> Raincloud Plot of Baseline Data
</p>


## 1. Likelihood Stress-Test Matrix: Evolutionary Benchmarking

This matrix tracks the performance of four distinct architectures across two operational phases: the Clean Baseline ($N=10$) and the High-Noise Stress Phase ($N=42$) containing dynamical "Black Swan" events (1.2, 8.8, 9.5).

| Phase | Model Architecture | Logic & Core Assumption | Outlier Handling | Trt 2 Noise ($\sigma_3$) | Signal Confidence | Audit Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Baseline ($N=10$)** | **$M_0$: Frequentist** | Identity $F$-test; $\sigma_1 = \dots = \sigma_j$ | Inflates $SS_{error}$ | 0.62 (Pooled) | $p < 0.05$ | Legacy Baseline |
| **Baseline ($N=10$)** | **$M_1$: Pooled Bayes** | Bayesian Mean; $\sigma_{shared}$ | Smears $\sigma$ globally | 0.71 (Global) | - | Optimal for clean data |
| **Baseline ($N=10$)** | **$M_2$: Separate Var** | Heteroscedastic; $\sigma_j$ unique | Distorts individual $\sigma_j$ | 0.73 | - | Overfit (DIC penalty) |
| --- | --- | --- | --- | --- | --- | --- |
| **Stress ($N=42$)** | **$M_0$: Frequentist** | Identity $F$-test; $\sigma_1 = \dots = \sigma_j$ | Fails $L_2$ norm | Inflated | Unreliable | **System Failure** |
| **Stress ($N=42$)** | **$M_1$: Pooled Bayes** | Bayesian Mean; $\sigma_{shared}$ | Smearing: High | 1.23 (Global) | Biased Low | **Fails PPC** |
| **Stress ($N=42$)** | **$M_2$: Separate Var** | Heteroscedastic; $\sigma_j$ unique | Triggers "Panic" mode | 1.29 (Panic) | Low Power | **Fragile / Boring** |
| **Stress ($N=42$)** | **$M_3$: Robust-t** | Signal Recovery; Outliers $\neq$ Errors | **Neutralized ($\nu=4$)** | **0.90 (Pure Signal)** | **99.8% Validated** | **High-Fidelity** |

The transition to $M_3$ was a mechanical necessity based on the following audit findings:

1. **Outlier Neutralization:** In $M_2$, the extreme values (1.2, 9.5) forced the Gaussian distribution to widen its "shoulders," resulting in a noise floor of 1.29. $M_3$ uses a Student-t likelihood with fixed degrees of freedom ($\nu=4$), which allows the model to categorize these spikes as "expected tail events" rather than structural variance. 
2. **Precision Recovery:** By filtering out dynamical noise, we recovered a true noise floor of 0.90 for Treatment 2, a 30.2% precision improvement over the fragile $M_2$ model. 
3. **Inhibitor Verification:** $M_3$ confirmed with 99.8% certainty that Treatment 1 is a growth inhibitor (Mean: 4.39), preventing a false-positive investment that could have occurred if noise had been "smeared" across the groups.

<table style="width: 100%; border-collapse: collapse; border: none;">
  <tr style="border: none;">
    <td align="center" style="border: none; width: 50%;">
      <img src="/Figures/PhaseII_Final_Audit_N42.png" width="100%" alt="Model Comparison Audit">
      <br>
      <b>Figure 2:</b> Model Comparison
    </td>
    <td align="center" style="border: none; width: 50%;">
      <img src="/Figures/Final_System_AuditBS242.png" width="100%" alt="Sensitivity Analysis Audit">
      <br>
      <b>Figure 3:</b> Sensitivity Analysis
    </td>
   </tr>
</table>

<p align="center">
  <img src="/Figures/Model_Comparison_BS2Audit.png" width="600" alt="Resized Chart Description">
  <br>
  <b>Figure 4:</b> Final System Audit
</p>
    
  

### 2. Black Swan Diagnosis (N = 42)
   The baseline challenge was identifying how different statistical architectures handled the $1.2, 8.8$, and $9.5$ extreme values.
<p align="center">
  <img src="/Figures/Power_sim2.png" width="75%" alt="Resized Chart Description">
  <br>
  <b>Figure 5:</b> Power Simulation
</p>

| Evaluation Phase | Architecture | Outcome | Conclusion |
| :--- | :--- | :--- | :--- |
| **Phase 1: $M_1$** | Pooled Variance | Smearing of noise across groups. | **FAILED.** Over-weighted Control noise. |
| **Phase 2: $M_2$** | Heteroscedastic Gaussian | Individual $\sigma$ revealed outliers. | **FRAGILE.** Outliers (1.2, 9.5) exploded $\sigma$. |
| **Phase 3: $M_3$** | Robust Student-t ($\nu=4$) | Mathematical "downgrading" of outliers. | **SUCCESS.** Recovered true noise floor. |

<p align="center">
   <img src="/Figures/Raincloud_SignalRecovery.png" width="600" alt="Raincloud Posterior Signal Recovery">
   <br>
   <b>Figure 6:</b> Raincloud Posterior Signal Recovery
</p>

### 3. Model 3 Performance & Outlier Neutralisation 
By transitioning to the Robust-t architecture, the model successfully isolated the true effect sizes from the dynamical noise.

<div align="center">
  <div style="display: inline-block; width: 45%; vertical-align: middle;">
    <img src="/Figures/RobustdensityBS2.png" width="600" alt="Model 3 Posterior Density Distribution">
    <br>
    <b>Figure 7:</b> Posterior Density Distribution
  </div><div style="display: inline-block; width: 45%; vertical-align: middle; margin-left: 2%;">
    <img src="/Figures/final_results.png" width="600" alt="Model 3 Final Engineering Audit">
    <br>
    <b>Figure 8:</b> Model 3 Engineering Audit
  </div>
</div>

By transitioning to the Robust-t architecture, the model successfully isolated the true effect sizes from the dynamical noise.

| Group | Gaussian $\sigma$ ($M_2$) | Robust $\sigma$ ($M_3$) | Detection Power |
| :--- | :--- | :--- | :--- |
| **Control** | 0.80 | 0.73 | Agreement (No Swans present) |
| **Treatment 1** | 1.47 | 1.09 | $M_3$ catches 3/3 Swans |
| **Treatment 2** | 1.29 | 0.90 | $M_3$ catches 3/3 Swans |

The transition from $M_2$ (Heteroscedastic Gaussian) to $M_3$ (Robust Student-t) reveals that Gaussian likelihoods were over-estimating noise by 30.2% to 34.8% in the presence of dynamical "Black Swans". 
* **Control Stability:** Both models agree where data is clean (0.80 vs 0.73), validating the baseline. 
* **Treatment 2 Recovery:** $M_3$ successfully "ignored" the 9.5 and 1.2 outliers to recover the true noise floor of 0.90, compared to the panicked 1.29 estimate in $M_2$.

## 4. Mathematical Backbone and System Specification 
The transition from $M_0$ to $M_3$ represents a shift from global averaging to local, robust signal processing.

### 4.1 Likelihood Architectures
* **$M_1$ (Homoscedastic):** $y_{i,j} \sim \mathcal{N}(\mu_j, \sigma^2)$
* **$M_2$ (Heteroscedastic):** $y_{i,j} \sim \mathcal{N}(\mu_j, \sigma_j^2)$
* **$M_3$ (Robust Shielded):** $y_{i,j} \sim \text{Student-t}(\mu_j, \sigma_j, \nu=4)$

### 4.2 Component Selection (Priors)
In this Bayesian framework, priors act as strictly defined System Constraints.

* **The "Slew Rate Limit" ($\nu=4$):** In a Gaussian model ($\nu \to \infty$), the squared-error loss forces the model to "chase" every data point. Setting $\nu=4$ defines a polynomial decay in the likelihood, acting as a slew rate limit that acknowledges Black Swans (1.2, 8.8, 9.5) without allowing them to pull the mean ($\mu$) or explode the variance ($\sigma$).
* **Noise Floor Ceiling ($\sigma \sim \text{Unif}(0, 5)$):** Defines the physical dynamic range. Prevents the MCMC from wandering into "high-noise hallucinations" during burn-in, ensuring convergence ($\hat{R} \approx 1.000$).
* **Hierarchical Tether ($\mu \sim \mathcal{N}(5, 0.01)$):** Ensures stability across groups with low sample density via a non-informative calibration point.

### 4.3 Demonstrating Outlier Capacity in Code
To show exactly how the model detects and handles outliers, we look at the JAGS implementation. The `dt` function replaces `dnorm`. Because $\nu=4$, the probability mass in the tails is vastly heavier. When $y_i=9.5$ is evaluated, `dt` assigns it a low but non-catastrophic probability, preventing the MCMC sampler from forcing $\sigma_j$ to expand to cover it.
```R
model {
  # 1. Hyper-priors (Hierarchical Tether)
  mu_grand ~ dnorm(5, 0.01) 
  tau_grand ~ dgamma(0.1, 0.1)

  # 2. Group Parameters (Heteroscedasticity)
  for (j in 1:n_groups) { 
    mu[j] ~ dnorm(mu_grand, tau_grand) 
    tau[j] <- pow(sig[j], -2) 
    sig[j] ~ dunif(0, 5) 
  }

  # 3. Robust Likelihood (Outlier Shielding)
  for (i in 1:N) { 
    # The Black Swan Filter: dt with nu=4 
    y[i] ~ dt(mu[grp[i]], tau[grp[i]], 4) 
  } 
}
```
<table style="width: 100%; border-collapse: collapse; border: none;">
  <tr style="border: none;">
    <td align="center" style="border: none; width: 50%;">
      <img src="/Figures/OutlierSensitivity.png" width="100%" alt="Outlier Sensitivity Chart">
      <br>
      <b>Figure 9:</b> Outlier Sensitivity
    </td>
    <td align="center" style="border: none; width: 100%;">
      <img src="/Figures/robust_noise_audit_N42.png" width="100%" alt="Noise Audit Chart">
      <br>
      <b>Figure 10:</b> Noise Audit
    </td>
  </tr>
</table>

## 5. Computational Hierarchy 
M_1(O(n)): Near-zero cost; useless for heteroscedastic systems.
M_2(O(k c.n)): Moderate cost; high risk of fragility.
M_3(O(k.c MCMC): Highest computational burden, but provides epistemic safety.

Audit Verdict: In our system, the "Robustness Tax" of $M_3$ is entirely justified by the recovery of a stable $0.90$ sigma. 
The MCMC converged efficiently with Information Criterion $\text{DIC} = 389.7$, $\hat{R} \approx 1.000$ (Upper C.I. $\leq 1.001$), and Effective Sample Sizes (ESS) $> 29,400$.a
<p align="center">
  <img src="/Figures/RobusttraceBS2.png" width="600" alt="Resized Chart Description">
  <br>
  <b>Figure 11:</b> MCMC Trace Plot 
</p>


6. Return On Investment & Strategic Roadmap ($N=42 \to N=125$)
   In engineering, ROI is the ratio of Net Value to Cost.
   The mathematical stability of the $M_3$ signal increases precision by approximately 30.2% relative to the volatile Gaussian baseline:
   
   <p align="center">
      <img src="/Figures/SNR.png" width="400" alt="Resized Chart Description">
      <br> 
  </p>

  6.1 Scaling to the Efficiency PeakTraditional ANOVA at $N=42$ is highly susceptible to Type II errors (false negatives).
   Current State ($N=42$): $\text{Power} \approx 43\%$.
   
   The probability of a False Negative is unacceptably high (57%).
   
   Target State ($N=125$): Reaching the Efficiency Peak provides $\sim 80\%$ Power.
   
   <p align="center">
      <img src="/Figures/ReliabilityGain.png" width="600" alt="Resized Chart Description">
      <br> 
  </p>

  6.2 The "Robustness Tax" vs. The "Error Cost"The transition to a Robust-t architecture (the "Robustness Tax") provides a defensible filter against dynamical noise that legacy frequentist models cannot replicate.

Avoided Waste: The audit identified Treatment 1 as a growth inhibitor ($0.4\%$ success probability). 
This prevents a massive R&D investment into a failing product that appeared viable under $M_1$ pooled assumptions.
Logistical Cost-Saving: Targeting $N=125$ instead of a standard $N=200$ reliability mandate saves 37.5% of the Phase III budget while maintaining an 80% detection threshold.

<p align="center">
  <img src="/Figures/BSStressTest.png" width="600" alt="Resized Chart Description">
  <br>
  <b>Figure 12:</b> Robust Model Test
</p>


Final Verdict:
The robust framework improved precision by 35%. Scaling to $N=125$ will improve operational reliability by 86% relative to the $N=42$ baseline, establishing a mathematically defensible foundation for deployment.
























