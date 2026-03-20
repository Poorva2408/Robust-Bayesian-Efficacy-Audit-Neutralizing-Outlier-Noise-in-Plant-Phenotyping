Technical Audit: PlantGrowth Robust Bayesian Architecture


This repository documents the transition from legacy frequentist ANOVA to a high-fidelity Robust Bayesian Hierarchical Architecture. The audit addresses the failure of standard Gaussian models in the presence of dynamical "Black Swan" interference.

Executive Summary
This technical audit evaluates the transition from legacy frequentist ANOVA to a Robust Bayesian Hierarchical Architecture for analyzing dynamical growth data. Traditional ANOVA and Pooled Bayesian models ($M_1$) failed to isolate treatment effects due to "Black Swan" outlier interference (values of $1.2, 8.8, 9.5$), which smeared noise across groups and produced a failing Posterior Predictive Check (PPC) p-value of $0.0475$.

While a Separate Variance Gaussian model ($M_2$) identified heteroscedasticity, it remained fragile; the outliers inflated the noise floor to $\sigma = 1.29$. By deploying a Robust-t Likelihood ($\nu=4$), we achieved a $35\%$ increase in signal precision, recovering the true noise floor of $\sigma = 0.90$. This model identifies Treatment 2 as a high-yield success with $99.8\%$ Bayesian certainty.

To mitigate the $57\%$ risk of false negatives inherent in the current $N=42$ sample, the audit identifies an Efficiency Peak at $N=125$. Scaling to this target provides $80\%$ statistical power, ensuring maximum decision reliability and a defensible ROI for Phase III deployment.

![Figure 1: Raincloud Plot of Baseline](/Figures/EDA_raincloud.png)


### 1. Likelihood Stress-Test Matrix: Evolutionary Benchmarking

This matrix tracks the performance of four distinct architectures across two operational phases: the Clean Baseline ($N=10$) and the High-Noise Stress Phase ($N=42$) containing dynamical "Black Swan" events ($1.2, 8.8, 9.5$).

| Phase | Model Architecture | Logic & Core Assumption | Outlier Handling | Trt 2 Noise ($\sigma_3$) | Signal Confidence | Audit Status |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Baseline ($N=10$)** | **$M_0$: Frequentist** | Identity $F$-test; $\sigma_1 = \dots = \sigma_j$ | Inflates $SS_{error}$ | $0.62$ (Pooled) | $p < 0.05$ | Legacy Baseline |
| **Baseline ($N=10$)** | **$M_1$: Pooled Bayes** | Bayesian Mean; $\sigma_{shared}$ | Smears $\sigma$ globally | $0.71$ (Global) | - | Optimal for clean data |
| **Baseline ($N=10$)** | **$M_2$: Separate Var** | Heteroscedastic; $\sigma_j$ unique | Distorts individual $\sigma_j$ | $0.73$ | - | Overfit (DIC penalty) |
| --- | --- | --- | --- | --- | --- | --- |
| **Stress ($N=42$)** | **$M_0$: Frequentist** | Identity $F$-test; $\sigma_1 = \dots = \sigma_j$ | Fails $L_2$ norm | Inflated | Unreliable | **System Failure** |
| **Stress ($N=42$)** | **$M_1$: Pooled Bayes** | Bayesian Mean; $\sigma_{shared}$ | Smearing: High | $1.23$ (Global) | Biased Low | **Fails PPC** |
| **Stress ($N=42$)** | **$M_2$: Separate Var** | Heteroscedastic; $\sigma_j$ unique | Triggers "Panic" mode | $1.29$ (Panic) | Low Power | **Fragile / Boring** |
| **Stress ($N=42$)** | **$M_3$: Robust-t** | Signal Recovery; Outliers $\neq$ Errors | **Neutralized ($\nu=4$)** | **$0.90$ (Pure Signal)** | **$99.8\%$ Validated** | **High-Fidelity** |

![Figure 2: Model_Comparison](Figures/PhaseII_Final_Audit_N42.png) 

![Figure 3: Sensitivity Analyis N = 42](/Figures/Final_System_AuditBS242.png)

The transition to $M_3$ was a mechanical necessity based on the following audit findings:
1. Outlier Neutralization: In $M_2$, the extreme values ($1.2, 9.5$) forced the Gaussian distribution to widen its "shoulders," resulting in a noise floor of $1.29$.
2. $M_3$ uses a Student-t likelihood with fixed degrees of freedom ($\nu=4$), which allows the model to categorize these spikes as "expected tail events" rather than structural variance.
3.Precision Recovery: By filtering out dynamical noise, we recovered a true noise floor of $0.90$ for Treatment 2, a $30.2\%$ precision improvement over the fragile $M_2$ model.
4.Inhibitor Verification: $M_3$ confirmed with $99.8\%$ certainty that Treatment 1 is a growth inhibitor (Mean: $4.39$), preventing a false-positive investment that could have occurred if noise had been "smeared" across the groups.
![Figure 4: Final N = 42](/Figures/Model_Comparison_BS2Audit.png)  

2. Black Swan Diagnosis (N = 42)
   The baseline challenge was identifying how different statistical architectures handled the $1.2, 8.8$, and $9.5$ extreme values.
   
![Figure 5: Raincloud Plot of Posterior N = 42](/Figures/Raincloud_SignalRecovery.png)

\begin{table}[h!]
\centering
\renewcommand{\arraystretch}{1.5}
\begin{tabular}{|l|l|l|l|}
\hline
\textbf{Evaluation Phase} & \textbf{Architecture} & \textbf{Outcome} & \textbf{Conclusion} \\ \hline
\textbf{Phase 1: $M_1$} & Pooled Variance & Smearing of noise across groups. & \textbf{FAILED.} Over-weighted Control noise. \\ \hline
\textbf{Phase 2: $M_2$} & Heteroscedastic Gaussian & Individual $\sigma$ revealed outliers. & \textbf{FRAGILE.} Outliers ($1.2, 9.5$) exploded $\sigma$. \\ \hline
\textbf{Phase 3: $M_3$} & Robust Student-t ($\nu=4$) & Mathematical "downgrading" of outliers. & \textbf{SUCCESS.} Recovered true noise floor. \\ \hline
\end{tabular}
\caption{Evolutionary Benchmarking of Likelihood Robustness under Black Swan Stress ($N=42$)}
\label{table:model_evolution}
\end{table}

![Figure 6: Power Simulation ](/Figures/Power_sim2.png)  







