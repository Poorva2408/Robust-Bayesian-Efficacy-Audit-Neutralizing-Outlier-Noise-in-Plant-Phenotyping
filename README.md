# Robust Bayesian Efficacy Audit: Neutralizing Outlier Noise in Plant Phenotyping

### Plant Growth Reliability Audit ($M_3$)

This repository contains the engineering audit for transitioning plant growth efficacy modeling from legacy Frequentist/Gaussian frameworks to a Robust Student-t Bayesian Architecture ($M_3$).

1. Executive Summary
   Legacy models ($M_1, M_2$) failed to account for extreme biological outliers in the growth catalog, leading to inflated variance ($\sigma = 1.29$) and high Type II error rates.
   
2. By implementing a Robust-t distribution, we successfully neutralized "Black Swan" growth spikes, recovering a stable signal and providing a mathematically defensible path to N=125 scaling.
   
3. ROI & Strategic Roadmap ($N=42 \to N=125$)The mathematical stability of the $M_3$ signal increases precision by approximately 30.2% relative to the volatile Gaussian baseline:
   <p align="center">
      <img src="/Figures/SNR.png" width="400" alt="Resized Chart Description">
      <br> 
  </p>
4. Scaling to the Efficiency Peak
      Traditional ANOVA at N=42 is highly susceptible to false negatives (missing a successful treatment).
      4.1 Current State ($N=42$): $\text{Power} \approx 43\%$.
      4.2 Target State ($N=125$): Reaching the Efficiency Peak provides $\sim 80\%$Power.

   <p align="center">
      <img src="/Figures/ReliabilityGain.png" width="600" alt="Resized Chart Description">
      <br> 
   </p>

### Final Verdict
The $M_3$ architecture is the only mathematically defensible framework for the $N=125$ scaling mandate. $M_2$ remains rejected due to catastrophic sensitivity to biological outliers. Treatment 2 is confirmed as a high-yield asset with 99.8% Bayesian certainty.

