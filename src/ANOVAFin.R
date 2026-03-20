# ANOVA 

library(ggplot2)
library(bayesplot)
library(ggdist)
library(tidyquant)
library(GGally)
library(knitr)
library(kableExtra)
data("PlantGrowth")

# ==========================================
# 1. EDA : Data Visualization 
# ==========================================
# 1.1 Box-plot 

EDA_Box <- ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) + # Transparent boxes
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) + # Shows every single plant
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "darkred") + # Highlights the Mean
  theme_minimal() + 
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Asset Distribution: PlantGrowth Weight by Treatment",
       subtitle = "Diamonds represent group means; jitter shows raw data density",
       x = "Experimental Group", y = "Weight (Units)") +
  theme(legend.position = "none")

ggsave(paste0("Figures/EDA_Box.png"), plot = EDA_Box, width = 8, height = 5)
# 1.2 Rain-cloud plot 

p_raincloud <- ggplot(PlantGrowth, aes(x = group, y = weight, fill = group)) +
  # 1. The "Cloud": Half-density plot
  stat_halfeye(
    adjust = .5, 
    width = .6, 
    .width = 0, 
    justification = -.3, 
    point_colour = NA
  ) +
  # 2. The "Umbrella": Standard boxplot
  geom_boxplot(
    width = .15, 
    outlier.shape = NA, 
    alpha = 0.5
  ) +
  # 3. The "Rain": Jittered raw data points
  stat_dots(
    side = "left", 
    justification = 1.1, 
    binwidth = .05
  ) +
  # 4. Aesthetics
  scale_fill_tq() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "darkred") + # Highlights the Mean
  theme_tq() +
  labs(
    title = "Asset Distribution: Raincloud Audit",
    subtitle = "Visualizing density, quartiles, and raw data clusters(Jitters), Mean(Dark Red Diamond)",
    x = "Experimental Group",
    y = "Weight (Units)"
  ) +
  coord_flip() # Flipping makes the "rain" effect more intuitive

# Save it for your GitHub 'figures' folder
ggsave("Figures/EDA_raincloud.png", plot = p_raincloud, width = 8, height = 6)


#1.3 Comparison Table ---
eda_comp <- data.frame(
  Criterion = c("Dist. Topology", "Raw Data Visibility", "Outlier Detection", "Feature Correlation", "Portfolio Tier"),
  Boxplot = c("Hidden", "Zero", "High", "No", "Entry"),
  Raincloud = c("Explicit (Density)", "Total (Jitter)", "Excellent", "No", "Professional"),
  Pairplot = c("Partial", "High", "Moderate", "Primary Purpose", "Engineering")
)

# This prints a clean table in your RMarkdown knit
kable(eda_comp, format = "markdown", caption = "EDA Visualization Comparison Matrix")
write.csv(eda_comp, "data/eda_comparison.csv")

# ==========================================
# 2. Model Fitting
# ==========================================

# 2.1 Frequentist Linear Model1 
lmod1 <- lm(weight ~ group, data = PlantGrowth)

summary_stats <- summary(lmod1)$coefficients # Captures Estimates, SE, t-values, and p-values

# 2. Save the ANOVA table as a clean data frame
library(broom)
anova_results <- tidy(anova(lmod1)) # Converts the messy ANOVA printout into a clean table

# 3. Export to CSV (for your 'data/' folder on GitHub)
write.csv(summary_stats, "data/frequentist_coefficients.csv")
write.csv(anova_results, "data/frequentist_anova.csv")


# 2.2 Frequentist Linear Model 

# Get Frequentist Estimates
lmod_cm <- lm(weight ~ -1 + group, data = PlantGrowth)
freq_stats <- summary(lmod_cm)$coefficients[, 1:2] # Mean and SE
freq_ci <- confint(lmod_cm)

# Bayesian Audit
# ==========================================
# 1. SETUP & UTILITIES
# ==========================================
library(rjags)
set.seed(82)

# Building the 'Power Tool' first
run_bayesian_audit <- function(mod_string, data, params, n_burn, inits, n_iter, model_name = "model" ) {
  
  # --- 1. The Plotting Engine ---
  generate_plots <- function(sim_obj, y_hat, residuals, name) {
    
    # A. MCMC Trace (Convergence Check)
    
    # vars(contains("")) Grab anything that has sig or mu : collective for mu and sig 
    
    p1 <- mcmc_trace(sim_obj, pars = vars(contains("mu"), contains("sig"))) +
      theme_minimal() +
      facet_text(size = 10) +
      scale_color_brewer(palette = "Blues")+ # Keep the pro-blue look
      labs(title = paste("MCMC Trace Diagnostics:", name))
    ggsave(paste0("Figures/trace_", name, ".png"), plot = p1, width = 8, height = 6)
    
    # B. Posterior Density (The Results)
    p2 <- mcmc_areas(sim_obj, pars = c("mu[1]", "mu[2]", "mu[3]"), prob = 0.95) +
      theme_minimal() +
      facet_text(size = 10) +
      scale_color_brewer(palette = "Blues")+ # Keep the pro-blue look
      labs(title = paste("Posterior Distributions:", name),
           subtitle = "95% Credible Intervals highlighted")
    ggsave(paste0("Figures/density_", name, ".png"), plot = p2, width = 8, height = 6)
    
    # C. Residuals (The Engineering Check)
    res_data <- data.frame(Predicted = y_hat, Residual = residuals)
    p3 <- ggplot(res_data, aes(x = Predicted, y = Residual)) +
      geom_point(size = 3, alpha = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(method = "loess", color = "blue", se = FALSE, linewidth = 0.5) +
      theme_minimal() +
      theme(plot.title = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 10)) + # This replaces facet_text
      
      labs(title = paste("Residual Diagnostics:", name),
           subtitle = "Flat blue line indicates no systematic bias")
    ggsave(paste0("Figures/resids_", name, ".png"), plot = p3, width = 8, height = 5)
  }
  
  # --- 2. Main logic  jags Implementation ---
  
  # 1. Initialize and update
  mod <- jags.model(textConnection(mod_string), data = data, inits = inits, n.chains = 3, quiet = TRUE)
  update(mod, n_burn)  # 
  
  # 2. Sample
  sim <- coda.samples(model = mod, variable.names = params, n.iter = n_iter)
  csim <- as.mcmc(do.call(rbind, sim))
  
  # 3. Diagnostics (proof of Quality)
  gelman_stats <- gelman.diag(sim)
  ess_stats <- effectiveSize(sim) # High ESS : Low Autocorrelation 
  
  # Calculate Deviance Information Criterion (native to rjags)
  dic_val <- dic.samples(mod, n.iter = 5e3, type = "pD")
  
  quality_report <- list(
    r_hat = gelman_stats$psrf,
    ess   = ess_stats,
    DIC   = dic_val
  )
  
  # 4. Analytics : Calculate our custom Engineering Metrics 
  
  # Treatment 1 Audit (The 'Risk' candidate)
  prob_better_1 <- mean(csim[, "mu[2]"] > csim[, "mu[1]"])
  # CI 95
  hpd_diff_1 <- HPDinterval(as.mcmc(csim[, "mu[2]"] - csim[, "mu[1]"]))
  
  # Treatment 2 Audit (The 'Asset' candidate)
  prob_better_2 <- mean(csim[, "mu[3]"] > csim[, "mu[1]"])
  
  # CI 95
  hpd_diff_2 <- HPDinterval(as.mcmc(csim[, "mu[3]"] - csim[, "mu[1]"]))
  
  # 5. Residuals for validation 
  
  pm_params <- colMeans(csim)
  # Dynamic indexing for mu
  mu_cols <- grep("mu", names(pm_params))
  yhat <- pm_params[mu_cols][data$grp]
  resids <- data$y - yhat
  
  #  Call the nested plotter
  
  generate_plots(sim, yhat, resids, model_name)
  

  return(list(
    sim  = sim, 
    csim = csim,
    yhat = yhat,
    resids = resids,
    summary(sim), 
    prob_better_1 = prob_better_1, prob_better_2 = prob_better_2,
    hpd_diff_1, hpd_diff_2,
    quality_report = quality_report
    ))
}
# inits declare 
# ==========================================
# 2. MODEL DEFINITIONS
# ==========================================
# 2.1 Pooled Model 

mod_string_1 = "model{
  for( i in 1: length(y)){
    y[i] ~ dnorm(mu[grp[i]], prec) # group level mean of each groups
  }
  
  for(j in 1:3){
    mu[j] ~ dnorm(0.0, 1.0/1.0e6)
  }
  
  prec ~ dgamma(5/2.0, 5*1.0/2.0)
  sig = sqrt(1.0 / prec)
  
}"

# 2.2 Separate Variance Model 

mod_string_2 = "model{
for( i in 1: length(y)){
  y[i] ~ dnorm(mu[grp[i]], prec[grp[i]]) # group level mean of each groups
}

for(j in 1:3){
  
  mu[j] ~ dnorm(0.0, 1.0/1.0e6)
  prec[j] ~ dgamma(5/2.0, 5*1.0/2.0)
  sig[j] = sqrt(1.0 / prec[j])
  
}
diff_31 = mu[3] - mu[1]  # used in Baysian Power Analysis
}"


# inits  
inits1 = function(){
  list("mu"=rnorm(3, 0.0, 100), "prec"=rgamma(1, 1.0, 1.0))
}

inits2 = function(){ 
  list("mu" = rnorm(3, 0.0, 100), "prec" = rgamma(3, 1.0, 1.0))
} 

# ==========================================
# 3. EXECUTION (The "Fitting" Happens Here)
# ==========================================

data_jags = list(y = PlantGrowth$weight, grp = as.numeric(PlantGrowth$group))

# 1. Run Audit 1 (Pooled Variance)
audit1 <- run_bayesian_audit(
  mod_string = mod_string_1, 
  data       = data_jags, 
  params     = c("mu", "sig"), 
  inits      = inits1, 
  n_burn     = 1000,
  n_iter     = 5000, 
  model_name = "Audit1_Pooled"
)

# 2. Run Audit 2 (Separate Variance)
audit2 <- run_bayesian_audit(
  mod_string = mod_string_2, 
  data       = data_jags, 
  params     = c("mu", "sig"), 
  inits      = inits2, 
  n_burn     = 1000,
  n_iter     = 5000, 
  model_name = "Audit2_Separate"
)

# 3. Save the 'Final Objects' for future reuse

saveRDS(audit1, "data/audit1_final.rds")
saveRDS(audit2, "data/audit2_final.rds")

# 4. Export the Master Summary Table for the 1-Pager 
audit1_sum <- summary(audit1$sim) 
final_table_1 <- cbind(audit1_sum$statistics, audit1_sum$quantiles) 
write.csv(final_table_1, "data/audit1_summary.csv") 

audit2_sum <- summary(audit2$sim) 
final_table_2 <- cbind(audit2_sum$statistics, audit2_sum$quantiles) 
write.csv(final_table_2, "data/audit2_summary.csv") 

# ==========================================
# 4. MODEL COMPARISON 
# ==========================================

# 4.1 Bayesian Model 1 vs Model 2 

# Using DIC Criterion 

audit1_dic <- audit1$quality_report$DIC 
print(audit1_dic) 
# Mean deviance:  58.97 
# penalty 4.102 
# Penalized deviance: 63.08 

audit2_dic <- audit2$quality_report$DIC 
print(audit2_dic) 
# Mean deviance:  61.2 
#penalty 5.691 
#Penalized deviance: 66.89 

# The Delta DIC 

Delta_DIC <- audit1_dic - audit2_dic 
print(Delta_DIC)
# Difference: -3.844633
# Sample standard error: 1.662332

# 1. Extract the actual scalar values correctly
m1_dic_val <- sum(audit1_dic$deviance) + sum(audit1_dic$penalty) 
m2_dic_val <- sum(audit2_dic$deviance) + sum(audit2_dic$penalty) 

# 2. Rebuild the DF with the correct math
model_comp <- data.frame(
  Model = c("Pooled Variance (Audit 1)", "Separate Variance (Audit 2)"),
  Mean_Deviance = c(sum(audit1_dic$deviance), sum(audit2_dic$deviance)),
  Penalty_pD = c(sum(audit1_dic$penalty), sum(audit2_dic$penalty)),
  DIC = c(m1_dic_val, m2_dic_val)
)

model_comp$Delta_DIC <- model_comp$DIC - min(model_comp$DIC)
knitr::kable(model_comp, digits = 2) 
write.csv(model_comp, "data/model_comparison.csv")

#|#Model                       | Mean_Deviance| Penalty_pD|   DIC| Delta_DIC|
#|#:---------------------------|-------------:|----------:|-----:|---------:|
#|#Pooled Variance (Audit 1)   |         58.97|       4.10| 63.08|      0.00|
#|#Separate Variance (Audit 2) |         61.20|       5.69| 66.89|      3.82|

# Model 1 Penalty: $~4.1$ (This makes sense: 3 means + 1 shared precision = ~4 parameters).
# Model 2 Penalty: $~5.8$ (This also makes sense: 3 means + 3 separate precisions = ~6 parameters).

# In Bayesian model selection, a difference of 3 to 7 points is considered "substantial" evidence. 
# Since Model 1 has the lower score, the data suggests that the extra complexity of indexing the precision 
# doesn't provide enough of a fit improvement to justify the extra parameters.

# ODI 
#1.Observation: Model 1 has a lower DIC by 3.82 units.

#2.Deduction: The separate variance model is "over-fitting" the small Dataset. The performance gain in likelihood 
# does not offset the penalty for the extra parameters. 

#3.Impact: For the final report, we should favor the Pooled Variance model for our 93.8% probability calculations, 
# as it is more parsimonious. 

#4. Frequentist Estimate vs Bayesian Model 

#5. Cell Means Model as the Primary Cross Validation tool for Bayesian Engine 

# Get Frequentist Estimates 
lmod_cm <- lm(weight ~ -1 + group, data = PlantGrowth)
freq_stats <- summary(lmod_cm)$coefficients[, 1:2] # Mean and SE
freq_ci <- confint(lmod_cm) 

# To ensure the JAGS engine was correctly initialized and that the choice of priors did not introduce 
# unintended bias, a Frequentist Cell Means Model was established as a baseline calibrator.

# The Cell Means formulation is preferred for Bayesian validation because it estimates the absolute 
# mean of each experimental group directly. This allows for a literal comparison between the Frequentist 
# 95% Confidence Intervals and the Bayesian 95% High Posterior Density (HPD) intervals.

# Technical Verdict: By using the Cell Means model as a calibrator, we prove that the 93.8% probability 
# of success for Treatment 2 is a robust signal derived from the data, not a byproduct of model parameterization 
# or prior "bullying."

# 2. Get Bayesian Estimates (from your Audit2_Separate results)

bayesian_stats_1 <- audit1_sum$statistics[grep("mu", rownames(audit1_sum$statistics)), 1:2]
bayesian_ci1 <- audit1_sum$quantiles[grep("mu", rownames(audit1_sum$quantiles)), c(1, 5)]
print(bayesian_ci1)

# Assuming 'audit2_sum' is the summary(audit2$sim) 
bayesian_stats_2 <- audit2_sum$statistics[grep("mu", rownames(audit2_sum$statistics)), 1:2]
bayesian_ci2 <- audit2_sum$quantiles[grep("mu", rownames(audit2_sum$quantiles)), c(1, 5)]
print(bayesian_ci2) 

# 3. Build the Validation Matrix 
validation_matrix <- data.frame(
  Parameter = c("Control", "Trt 1", "Trt 2"),
  Freq_Mean = freq_stats[,1],
  # Changed index from [,2] (SD) to [,1] (Mean)
  Bayes_Mean_M1 = bayesian_stats_1[,1], 
  Bayes_Mean_M2 = bayesian_stats_2[,1], 
  Freq_95_CI = paste0("[", round(freq_ci[,1], 2), ", ", round(freq_ci[,2], 2), "]"), 
  Bayes_95_HPD_M1 = paste0("[", round(bayesian_ci1[,1], 2), ", ", round(bayesian_ci1[,2], 2), "]"),
  Bayes_95_HPD_M2 = paste0("[", round(bayesian_ci2[,1], 2), ", ", round(bayesian_ci2[,2], 2), "]") 
  # Comma removed here to fix the "missing argument" error
)

knitr::kable(validation_matrix, 
             digits = 3, 
             caption = "System Calibration: Frequentist (Cell Means) vs. Bayesian Models")
write.csv(validation_matrix, "data/Validation_matrix.csv")

# 2. Combine into one data frame
plot_df_3way <- data.frame(
  Group  = rep(c("Control", "Trt1", "Trt2"), 3),
  Method = rep(c("Frequentist (lm)", "Bayesian (Mod 1)", "Bayesian (Mod 2)"), each = 3),
  Mean   = c(freq_stats[,1], bayesian_stats_1[,1], bayesian_stats_2[,1]),
  Lower  = c(freq_ci[,1], bayesian_ci1[,1], bayesian_ci2[,1]),
  Upper  = c(freq_ci[,2], bayesian_ci2[,2], bayesian_ci2[,2])
)
# 3. Plot it  
SensitivityAnalysis <- ggplot(plot_df_3way, aes(x = Group, y = Mean, color = Method)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                position = position_dodge(width = 0.6), width = 0.2) +
  theme_minimal() +
  labs(title = "Sensitivity Analysis: Weight Estimates Across Models",
       subtitle = "Comparing Pooled vs. Separate Variance Assumptions")
ggsave(paste0("figures/Sensitivity_Analysis.png"), plot = SensitivityAnalysis, width = 8, height = 5)

# power Analysis for model 1 

# Parameters from your JAGS output
M1_mu = audit1_sum$statistics[grep("mu", rownames(audit1_sum$statistics)), 1:1]
delta = M1_mu[3] - M1_mu[1] # 0.49048 : The treatment effect (mu3 - mu1)
print(delta) # 0.489
sigma = audit1_sum$statistics[grep("sig", rownames(audit1_sum$statistics)), 1:1] # 0.7119115 The pooled standard deviation from Model 1
print(sigma) # 0.71

# Risk Utility Trade-off / Jacob COhen 

# Power Analysis for 80% Power at alpha = 0.05 
# Signal to Noise Ratio calculation for Experimental design 
# alpha = 0.05 False Positive rate (Type I error) 
# It is the probability of claiming an effect exists when it doesn't. We set this low (5%) because claiming a 
# "breakthrough" that is actually noise is scientifically damaging.

# Beta = 0.2 False Negative rate (Type II Error)  probability of missing an effect that actually exists. 
# Power (1 - Beta)  = 0.8 The Probability of Success. It means if the effect is real, you have an 80 % chance of your 
# p - value being <= 0.05

# 0.20 / 0.05 = 4 this 4: 1 weighting assumes that "crying wolf" is a greater sin than "missing a discovery."

# Economics of Sampling: Power has diminishing returns. Increasing power from 80% to 95% often requires doubling 
# or tripling the sample size ($N$). For most researchers and engineers, 80% is the "sweet spot" where you have 
# a high enough chance of success without bankrupting the project on data collection.hing returns.


# This calculates how many plants per group you need to 'catch' the 0.5 difference 
power_res <- power.t.test(delta = delta, sd = sigma, sig.level = 0.05, power = 0.8) 

print(power_res$n) 
# Result: ~34 plants per group'

# cross check 
# By scaling to $N=34$, the Standard Error ($SE$) of the treatment effect is reduced enough to ensure the "Signal" 
# is never lost in the "Noise."

SE = 0.71 * sqrt(2/34)
print(SE) # 0.1722

new_upper = delta + 2*SE
print(new_upper) # 0.8348

new_lower = delta - 2*SE 
print(new_lower) # 0.1460  

new_intervals = c(new_lower, new_upper)
print(new_intervals) #   mu[3]     mu[3] 
                     #  0.1460    0.8348

# Pooled Variance Model ($\delta = 0.49$, $\sigma = 0.71$), the study requirements for a robust detection are:
# Target: 80% Power at $\alpha = 0.05$.
# Required Sample Size ($N$): 34 plants per group (102 total).
# Current Status: The pilot ($N=10$) was significantly underpowered, carrying a high risk of a Type II error (False Negative).

## Power Analysis for Model 2 
# Standard power analysis assumes a 1:1 ratio ($n_1 = n_2$). In a separate variance model, you can use Neyman 
# Allocation. If Treatment 2 has twice the standard deviation of the Control, you should allocate more samples 
# to Treatment 2 to minimize the total standard error of the difference.

# Gold standard for industrial data science 
# Instead of relying on 100-year-old math formulas, you are using your Model 2 (Separate Variances) as a 
# Generative Engine to stress-test your experimental design.

library (coda)

# 1. Ground Truth Parameters (Extracted from your Model 2 results)
# Replace these with your actual mu and sig values from audit2_sum
#true_mu  <- c(5.03, 4.66, 5.52)  # ctrl, trt1, trt2
#true_sig <- c(0.58, 0.79, 0.44)  # Separate sigmas from Model 2

M2_mu = audit2_sum$statistics[grep("mu", rownames(audit2_sum$statistics)), 1]

true_mu = c(M2_mu[1], M2_mu[2], M2_mu[3])
print(true_mu)

M2_sig = audit2_sum$statistics[grep("sig", rownames(audit2_sum$statistics)), 1] # 0.7119115 The pooled standard deviation from Model 1

true_sig = c(M2_sig[1], M2_sig[2], M2_sig[3])
print(true_sig)

Bayesian_power = function(N_per_group, n_sims = 500, model_string, true_mu, true_sig) {
  # Force numeric to prevent naming conflicts in rnorm
  true_mu <- as.numeric(true_mu)
  true_sig <- as.numeric(true_sig)
  
  success_count = 0 
  pb <- txtProgressBar(min = 0, max = n_sims, style = 3) # Progress bar initialization
  
  for(i in 1:n_sims) {
    # 1. Generate Synthetic Data
    y <- c(rnorm(N_per_group, true_mu[1], true_sig[1]),
           rnorm(N_per_group, true_mu[2], true_sig[2]),
           rnorm(N_per_group, true_mu[3], true_sig[3]))
    grp <- rep(1:3, each = N_per_group)
    
    # 2. Fit Model
    jags_sim <- jags.model(textConnection(model_string), 
                           data = list(y=y, grp=grp), quiet=TRUE)
    
    update(jags_sim, 500) 
    res <- coda.samples(jags_sim, variable.names="diff_31", n.iter=1000)
    
    # 3. Check HPD
    hpd <- HPDinterval(as.mcmc(res))[[1]]
    if(hpd[1] > 0) { 
      success_count <- success_count + 1 
    }
    setTxtProgressBar(pb, i) # Update bar
  }
  close(pb)
  return(success_count / n_sims)
}

# ---------------------------------------------------------
# EXECUTION (The Audit Standard)
# ---------------------------------------------------------
set.seed(42) # Essential for reproducibility


# 3. Execution (Test N = 10, 20, 30)
power_at_10 <- Bayesian_power(N_per_group = 10, 
                              n_sims = 500, 
                              model_string = mod_string_2, 
                              true_mu = true_mu, 
                              true_sig = true_sig)

print(paste("Bayesian Power at N=10:", power_at_10)) # 0.18

power_at_30 <- Bayesian_power(N_per_group = 30, 
                              n_sims = 500, 
                              model_string = mod_string_2, 
                              true_mu = true_mu, 
                              true_sig = true_sig)


print(paste("Bayesian Power at N=30:", power_at_30)) # 0.632

# Model Audit: DIC and Power simulations revealed that while Treatment 2 is effective, the pilot was underpowered 
# ($Power = 0.13$).Strategic Recommendation: Deploy a Phase II Validation with $N=34$ per group ($N=102$ total) 
# to achieve $80\%$ statistical power and confirm the observed $0.49$ unit biomass increase before full-scale 
# implementation.


# Phase II validation N = 34 
# Expecting: ~0.81 to 0.84 
power_curve_34 <- Bayesian_power(34, n_sims = 500, mod_string_2, true_mu, true_sig)
print(paste("Stabilized Power at N=34:", power_curve_34)) #0.712

# The "Bayesian Penalty" for having separate variances means you need slightly more data to reach the same level of certainty as a pooled model.

# Decision Matrix 

# Define the Decision Matrix
decision_matrix <- data.frame(
  Option = c("A: Accept Baseline", "B: Push to Industrial Target"),
  Action = c("Stay at N=34", "Increase to N=42"),
  New_Power = c("71%", "80% (Projected)"),
  Engineering_Logic = c(
    "5.5x improvement over pilot; balances cost vs. reliability.",
    "Accounts for the 'Uncertainty Tax' of Separate Variances."
  )
)

# Render with kableExtra
decision_matrix %>%
  kbl(caption = "Phase II Experimental Design: Decision Matrix",
      booktabs = TRUE, 
      col.names = c("Option", "Action", "Simulated Power", "Engineering Logic")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), 
                full_width = F, 
                position = "left") %>%
  column_spec(1, bold = TRUE, border_right = TRUE) %>%
  column_spec(4, width = "15em")

write.csv(decision_matrix, "data/decision_matrix.csv")
# N = 42 
power_curve_42 <- Bayesian_power(42, n_sims = 500, mod_string_2, true_mu, true_sig)
print(paste("Stabilized Power at N=34:", power_curve_42)) # 0.804


# In statistics, $N=42$ is your margin. It ensures that the heteroscedasticity (the fluctuating noise levels 
# between groups) doesn't result in a "false negative" during the final verification.


# Analyst Moves 
#1.The 3-Way Comparison (lm vs. Mod1 vs. Mod2): This confirms Consistency. It proves the growth "signal" 
# is robust and not a fluke of one specific model's assumptions.

#2.The DIC Check: This is your Occam’s Razor. It prevents you from over-complicating the story by proving 
# that Treatment 1's extra noise isn't worth a more complex model.

#3.The PPC (Posterior Predictive Check): This is your Reality Check. It proves to stakeholders that 
# your math can actually "re-create" the physical weights seen in the garden.

#4.The Forest Plot: This is your Strategic Communication. It visually explains to non-statisticians 
# exactly why Treatment 2 is the winner, even if the "error bars" (uncertainty) are still a bit too wide.

# "While Treatment 2 shows a consistent $0.5$ unit mean increase over the Control across all tested 
# models (Frequentist, Pooled Bayesian, and Heteroscedastic), the current sample size of $N=10$ provides 
# a 93.8% posterior probability of success. To achieve a 95% HPD interval entirely above zero, a power 
# analysis suggests a minimum of 27 plants per group for future iterations of this study."


## Posterior Predictive Check 
# 1. Posterior Predictive Check: Generating new data from the model

# Declare y and grp from your source data
y   <- PlantGrowth$weight  # The raw plant weights
grp <- as.numeric(PlantGrowth$group) # The group indices (1, 2, 3)

# 'y_pred' will now simulate the 'messiness' accurately


n_sim = 2000 # Increased for smoother histograms
y_pred = matrix(NA, nrow = n_sim, ncol = length(data_jags$y))
ppc_m1 = matrix(NA, n_sim, length(y))
ppc_m2 = matrix(NA, n_sim, length(y))

for (i in 1:n_sim) {
  # Model 1: Pooled (Single Sig)
  ppc_m1[i,] = rnorm(length(y), audit1$csim[i, 1:3][grp], audit1$csim[i, "sig"])
  
  # Model 2: Separate (Indexed Sig)
  mu2 = audit2$csim[i, grep("mu", colnames(audit2$csim))]
  sig2 = audit2$csim[i, grep("sig", colnames(audit2$csim))]
  ppc_m2[i,] = rnorm(length(y), mu2[grp], sig2[grp])
}

# The Litmus Test: Variance of Group 2 (Treatment 1)
obs_var_trt1 <- var(y[grp == 2])
print(obs_var_trt1) #0.6299211
pred_var_m1  <- apply(ppc_m1[, grp == 2], 1, var)
pred_var_m2  <- apply(ppc_m2[, grp == 2], 1, var)

p_val_m1 <- mean(pred_var_m1 > obs_var_trt1) # 0.265
p_val_m2 <- mean(pred_var_m2 > obs_var_trt1) # 0.6195
print(p_val_m1)
print(p_val_m2)

ppc_results <- data.frame(
  Metric = c("DIC Score (Parsimony)", "Bayesian p-value (Trt1 Var)", "Interpretation"),
  Model_1_Pooled = c("63.08 (Winner)", "0.265", "Underestimates localized noise."),
  Model_2_Separate = c("66.89", "0.6195 (Winner)", "Captures heteroscedastic profile.")
)

ppc_results %>%
  kbl(caption = "Model Integrity Audit: DIC vs. Posterior Predictive Check") %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>%
  column_spec(2, color = ifelse(ppc_results$Model_1_Pooled == "0.265", "red", "black")) %>%
  column_spec(3, bold = T, color = "darkgreen")
write.csv(ppc_results, "data/ppc_results.csv")

library(patchwork)
library(tidyr)
library(dplyr)

# --- 0. Precise Dimension Locking ---
n_to_plot <- min(100, nrow(ppc_m1)) 
df_obs <- data.frame(weight = as.numeric(PlantGrowth$weight))

# --- 1. Robust Data Prep ---
# Model 1 
sim_m1 <- as.data.frame(ppc_m1[1:n_to_plot, ]) %>%
  mutate(sim_id = row_number(), model = "Model 1") %>%
  pivot_longer(cols = starts_with("V"), names_to = "plant_idx", values_to = "weight")

# Model 2
sim_m2 <- as.data.frame(ppc_m2[1:n_to_plot, ]) %>%
  mutate(sim_id = row_number(), model = "Model 2") %>%
  pivot_longer(cols = starts_with("V"), names_to = "plant_idx", values_to = "weight")

combined_sims <- bind_rows(sim_m1, sim_m2)

# --- 2. Patchwork Construction ---

# Plot A: Spaghetti PPC 
ppc <- ggplot() + 
  # LAYER 1: The Simulations (Mapped to combined_sims)
  geom_density(data = combined_sims, 
               aes(x = weight, group = interaction(sim_id, model), color = model), 
               alpha = 0.2, linewidth = 0.3) +
  # LAYER 2: The Observed Data (Mapped to df_obs, independent of simulations)
  geom_density(data = df_obs, aes(x = weight), 
               color = "red", linewidth = 1.2) +
  facet_wrap(~model) +
  scale_color_manual(values = c("Model 1" = "gray60", "Model 2" = "skyblue")) +
  theme_minimal() + 
  theme(legend.position = "none") +
  labs(title = "A: Global Density PPC", 
       subtitle = "Red: Observed | Gray/Blue: Simulated", 
       x = "Plant Weight")

# Plot B: Variance Audit
var_df <- data.frame(
  Variance = c(pred_var_m1, pred_var_m2),
  Model = rep(c("Model 1", "Model 2"), each = length(pred_var_m1))
)

histogram <- ggplot(var_df, aes(x = Variance, fill = Model)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 40) +
  geom_vline(xintercept = obs_var_trt1, color = "black", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("Model 1" = "#E41A1C", "Model 2" = "#377EB8")) +
  theme_minimal() + 
  labs(title = "B: Treatment 1 Variance Audit", x = "Simulated Variance", y = "Count")

# Combine and Print
final_plot <- ppc / histogram
print(final_plot)
ggsave(paste0("Figures/ppc+histogram.png"), plot = final_plot, width = 8, height = 5)


# 1. Global Density Overlay (The "Spaghetti" alternative)
# This compares your actual y to 50 simulated datasets from each model
p1_m1 <- ppc_dens_overlay(y, ppc_m1[1:50, ]) + ggtitle("Model 1: Pooled")
p1_m2 <- ppc_dens_overlay(y, ppc_m2[1:50, ]) + ggtitle("Model 2: Separate")

# 2. The Variance Check (The "Litmus Test")
# This plots the distribution of the 'stat' (variance) against the real value
p2_m1 <- ppc_stat(y, ppc_m1, stat = "var") + ggtitle("M1 Variance Fit")
p2_m2 <- ppc_stat(y, ppc_m2, stat = "var") + ggtitle("M2 Variance Fit")

# Combine using patchwork
global_density <- (p1_m1 + p1_m2) / (p2_m1 + p2_m2) 
ggsave(paste0("Figures/global_density.png"), plot = global_density, width = 8, height = 5)

### Models Noise profile 
# Extract Sigma for Treatment 1 (Group 2)
# Model 1 uses one 'sig' for all; Model 2 uses 'sig[2]' for Trt1
sig_m1 <- audit1$csim[, "sig"] 
sig_m2_trt1 <- audit2$csim[, "sig[2]"] 

# Summary Statistics 
summary_noise <- data.frame( 
  Model = c("Model 1 (Pooled)", "Model 2 (Trt1 Specific)"),
  Median_Sigma = c(median(sig_m1), median(sig_m2_trt1)),
  Lower_95 = c(quantile(sig_m1, 0.025), quantile(sig_m2_trt1, 0.025)),
  Upper_95 = c(quantile(sig_m1, 0.975), quantile(sig_m2_trt1, 0.975))
)
write.csv(summary_noise, "data/noise.csv")
# Calculate the 'Underestimation Percentage'
gap = (median(sig_m2_trt1) - median(sig_m1)) / median(sig_m1) * 100
cat("Model 1 is hiding approximately", round(gap, 1), "% of the noise in Treatment 1.")
# 27.3 % noise in Trt 1

# Extract Sigmas
sig_m1 <- audit1$csim[, "sig"]

# Internal Divergence in Model 2 
sig_m2_trt1 <- audit2$csim[, "sig[2]"]
sig_m2_trt2 <- audit2$csim[, "sig[3]"]

# Function to calculate Gap
calc_gap <- function(m2_sig, m1_sig) {
  round((median(m2_sig) - median(m1_sig)) / median(m1_sig) * 100, 1)
}

# To see how Model 2 treats the treatments compared to each other, we look at the Ratio of Variances within the 
# model itself.
# In Model 1, this ratio is locked at $1:1$ (pooled). 
# In Model 2, the model "notices" the difference. If Treatment 1 is significantly noisier than Treatment 2, 
# Model 2 will reflect that in the posterior distribution of its separate sigma parameters.


# Ratio of noise: Trt1 vs Trt2
noise_ratio_m2 <- median(sig_m2_trt1) / median(sig_m2_trt2)
print(noise_ratio_m2)

# Percentage difference within Model 2
internal_diff <- (median(sig_m2_trt1) - median(sig_m2_trt2)) / median(sig_m2_trt2) * 100
print(internal_diff)
cat("Model 2 identifies that Treatment 1 is", round(internal_diff, 1), "% noisier than Treatment 2.")
# Model 2 identifies that Treatment 1 is 25.6 % noisier than Treatment 2.

# Probability that Trt 1 is noisier than Trt 2
prob_trt1_noisier <- mean(audit2$csim[, "sig[2]"] > audit2$csim[, "sig[3]"])
cat("Probability that Trt 1 > Trt 2 variance:", round(prob_trt1_noisier * 100, 2), "%\n")
# Probability that Trt 1 > Trt 2 variance: 79.75%

# Construct Audit Table

# Re-build with the 7th column to match your col.names call
audit_data <- data.frame(
  Group = c("Treatment 1", "Treatment 2"),
  Obs_Var = c(0.630, 0.196),
  M1_Sigma = c(0.7129, 0.7129),
  M2_Sigma = c(0.9256, 0.7363),
  M1_Bias = c("27.3%", "1.9%"),
  M2_Delta = c("Ref (Noisiest)", "-25.6% vs Trt1"),
  System_Risk = c("CRITICAL: Power Loss", "Negligible")
)

audit_data %>%
  kbl(caption = "Statistical Audit: Cross-Treatment Variance Suppression and Group Divergence",
      col.names = c("Group", "Obs. Var", "M1 Sigma (Pooled)", "M2 Sigma (Local)", 
                    "M1 Bias", "M2 Internal Delta", "System Risk")) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>%
  column_spec(5, bold = T, color = "white", background = "#D7261E") %>%
  column_spec(6, italic = T, color = "darkblue") %>%
  column_spec(7, bold = T, color = ifelse(audit_data$System_Risk == "Negligible", "black", "red")) %>%
  add_footnote("Model 1 assumes σ1=σ2; Model 2 recovers the 25.6% volatility gap.")

write.csv(audit_data, "data/audit_data.csv")

# Posterior Probability Distribution of delta 
# Calculate the posterior difference
# 1. Calculate the raw posterior difference
# Ensure these are numeric vectors
library(HDInterval)

# 1. Prepare Data and Stats
sigma_delta_vec <- as.numeric(audit2$csim[, "sig[2]"] - audit2$csim[, "sig[3]"])
plot_df <- data.frame(sigma_delta = sigma_delta_vec)
prob_greater <- mean(sigma_delta_vec > 0) 
z_score_noise <- mean(sigma_delta_vec) / sd(sigma_delta_vec) 
hdi_limits <- hdi(sigma_delta_vec, credMass = 0.95)

# 2. Build the Annotated Plot
post_prob <- ggplot(plot_df, aes(x = sigma_delta)) +
  # The Distribution
  geom_density(fill = "steelblue", alpha = 0.6) +
  # The HDI Segment (Visualizes the uncertainty)
  annotate("segment", x = hdi_limits[1], xend = hdi_limits[2], y = 0.05, yend = 0.05, 
           linewidth = 2, color = "darkblue") +
  annotate("text", x = mean(hdi_limits), y = 0.12, label = "95% HDI", color = "darkblue", size = 3) +
  # The Zero-Line (Reference)
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  # The Consolidated Data Box
  annotate("label", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, label_size = 0.5,
           label = paste0("AUDIT METRICS:\n",
                          "Prob(Trt1 > Trt2): ", round(prob_greater * 100, 1), "%\n",
                          "Bayesian Z-score: ", round(z_score_noise, 3), "\n",
                          "95% HDI: [", round(hdi_limits[1], 2), ", ", round(hdi_limits[2], 2), "]\n",
                          "Verdict: Moderate Heteroscedasticity\n",
                          "Recommendation: N=42 Risk-Adjusted"),
           fill = "white", alpha = 0.8, family = "mono") +
  theme_minimal() +
  scale_x_continuous(breaks = c(-1, -0.5, 0, 0.5, 1))+
  labs(title = "Structural Heteroscedasticity Audit",
       subtitle = "Testing the assumption of σ1 = σ2 (Homoscedasticity)",
       x = "Sigma Difference (Treatment 1 - Treatment 2)",
       y = "Posterior Density")
print(post_prob)
ggsave(paste0("Figures/posterior_prob.png"), plot = post_prob, width = 8, height = 5)
### conclusion 
# The 80.13% posterior probability and Z-score of 0.788 confirm a systematic rightward shift toward higher 
# entropy in Treatment 1.
# While the 95% HDI [-0.296, 0.671] overlaps zero, the median 27.7% noise underestimation in the pooled model 
# constitutes a critical structural risk to experimental power. Consequently, we reject the $N=34$ homoscedastic 
# projection in favor of a risk-adjusted $N=42$ to stabilize the signal-to-noise ratio against this observed 
# heteroscedasticity.


# Verdict 
# 1. Model 1 as a Smoothing Filter:
## Model 1 operates like a low-pass filter with an inappropriate cutoff. It observes the $0.63$ variance of 
## Treatment 1 and the $0.19$ variance of Treatment 2 and forces a compromise estimate of $\sigma = 0.70$. 
## It essentially "smooths" away the specific volatility of the treatment groups to satisfy an assumption of 
## homoscedasticity.

#2. The Risk (The Fuse Analogy):In an electronics context, this is equivalent to setting a circuit's fuse rating 
## based on the average current rather than the peak transient current. If you design your sample size ($N$) based 
## on average noise, the high-entropy bursts in Treatment 1 will "blow the fuse" of your statistical power. 
## You will fail to detect a real effect simply because your model refused to acknowledge the noise floor.

#3. The Solution:By accepting the 27.7% bias revealed by Model 2, we are forced to move from an optimistic $N=34$ 
## to a robust $N=42$. This ensures that even the "noisiest" sector of the experiment (Treatment 1) has a high 
## enough sample count to lower the Standard Error of the Mean ($SEM$) below the threshold of detection.


#######################################################################################################
# Synthetic study  Monte Carlo Power simulation (analyzing the past -> designing the future)

# Power Simulation 
# Instead of fearing a Black Swan, we simulate them. We calculate the exact % of times the "randomness" 
# will win over the "signal."

M2_mu = audit2_sum$statistics[grep("mu", rownames(audit2_sum$statistics)), 1:1]
delta_M2 = M2_mu[3] - M2_mu[1] # 0.49048 : The treatment effect (mu3 - mu1)
print(delta_M2) # 0.497

sigma_Trt1 = audit2_sum$statistics[grep("sig", rownames(audit2_sum$statistics)), 1:1][2] # 0.7119115 The pooled standard deviation from Model 1
print(sigma_Trt1) # 0.925

# --- 1. The High-Entropy Parameters ---
delta_signal <- delta_M2   # Expected Treatment Effect mu[3] - mu[1]
sigma_noise  <- sigma_Trt1  # M2-Trt1 Sigma (The "Black Swan" Generator)
iterations   <- 10000  # Law of Large Numbers (LLN) convergence 

# --- 2. Meta-Simulation Function ---
# This simulates 10,000 'universes' to see how often the Black Swan wins
run_power_sim <- function(n) {
  p_values <- replicate(iterations, {
    # Generate 'freak' data based on M2-Trt1 variance
    control <- rnorm(n, mean = 5.03, sd = 0.63) # Control is stable # can use dynamic val from audit2_sum
    treat_1 <- rnorm(n, mean = 5.03 + delta_signal, sd = sigma_noise) # Trt1 is chaotic
    t.test(control, treat_1)$p.value
  })
  return(mean(p_values < 0.05)) # Percentage of successful detection 
}

# --- 3. Execute Audit across N-Range ---
png("Figures/Black_Swan_Stress_Test.png", width = 1000, height = 800, res = 120)
n_range <- seq(10, 50, by = 5) 
power_results <- sapply(n_range, run_power_sim) 

# --- 4. The Roadmap Plot ---
#
plot(n_range, power_results, type = "b", pch = 19, col = "#0072B2",
     xlab = "Sample Size (N per group)", ylab = "Probability of Success (Power)",
     main = "Black Swan Stress Test: M2-Trt1 Entropy",
     ylim = c(0, 1))

# Add Thresholds
abline(h = 0.80, col = "red", lty = 2) # The Safety Line
text(15, 0.85, "80% Power Threshold", col = "red", font = 2)

# Mark the N=42 Requirement
points(42, run_power_sim(42), col = "darkgreen", cex = 2, lwd = 2)
text(35, 0.65, "N=42: Risk-Adjusted Goal", pos = 4, font = 2, col = "darkgreen")
#png("Black_Swan_Stress_Test.png", width = 800, height = 600, res = 120)
dev.off()

# --- 1. Data Preparation ---
# Re-using your simulation logic but wrapping in a dataframe
n_range <- seq(10, 60, by = 2) 
# Note: Using your run_power_sim function from the previous step
power_results <- sapply(n_range, run_power_sim) 
sim_df <- data.frame(N = n_range, Power = power_results)

# Calculate exact power at N=42 for the label
power_at_42 <- run_power_sim(42)

# --- 1. Data Shifting (Ensure sim_results reflects M2-Trt1 Sigma) ---
# Assuming sim_results contains N and Power from your run_power_sim(sigma=0.896)

power_sim2 <- ggplot(sim_df, aes(x = N, y = Power)) + 
  # Region 1: Gamble Zone (Now extended due to 27.7% extra noise)
  annotate("rect", xmin = 0, xmax = 25, ymin = 0, ymax = 1.1, fill = "red", alpha = 0.1) +
  # Region 2: Reliability Pivot (The climb toward N=42)
  annotate("rect", xmin = 25, xmax = 42, ymin = 0, ymax = 1.1, fill = "green", alpha = 0.1) +
  # Region 3: Efficiency Zone (Diminishing returns beyond the mandate)
  annotate("rect", xmin = 42, xmax = 60, ymin = 0, ymax = 1.1, fill = "grey", alpha = 0.1) +
  
  # The Stochastic Step Line
  geom_step(color = "#0072B2", linewidth = 1.2) +
  geom_point(color = "#0072B2", size = 2) +
  
  # 80% Power Threshold line 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgreen", linewidth = 1) +
  
  # Vertical Marker for the N=42 Mandate 
  geom_vline(xintercept = 42, color = "darkgreen", linetype = "dotted") +
  
  # Zone Labels (Placed at the top for better visibility)
  annotate("text", x = 12.5, y = 1.05, label = "GAMBLE ZONE", color = "red", fontface = "bold", size = 3) +
  annotate("text", x = 33.5, y = 1.05, label = "RELIABILITY PIVOT", color = "darkgreen", fontface = "bold", size = 3) +
  annotate("text", x = 51, y = 1.05, label = "EFFICIENCY ZONE", color = "grey40", fontface = "bold", size = 3) +
  
  # Key Milestones
  annotate("label", x = 10, y = 0.28, label = "Current N=10\n(72% Failure Risk)", color = "red", size = 3) +
  annotate("label", x = 42, y = 0.82, label = "Target N=42\n(Black Swan Shield)", color = "darkgreen", size = 3, fontface = "bold") +
  
  # Aesthetics
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1.1), labels = scales::percent) +
  scale_x_continuous(breaks = seq(0, 60, 10)) +
  theme_minimal() +
  labs(title = "Stochastic Power Roadmap: The 'Black Swan' Defense",
       subtitle = "Risk-Adjusted for Model 2 / Treatment 1 Entropy (Sigma = 0.896)",
       x = "Sample Size (n per group)",
       y = "Observed Power (1 - Beta)")
print(power_sim2)
ggsave(paste0("Figures/Power_sim2.png"), plot = power_sim2, width = 8, height = 5)

audit_table <- data.frame(
  Feature = c("Likelihood", "Priors", "Sample Strategy", "Goal"),
  Standard_Approach = c("Normal (Rigid)", "Flat/Independent", "N=10 (Gamble)", "Hope for Significance"),
  Robust_Approach = c("Student-t (Flexible)", "Hierarchical (Tethered)", "N=42 (Shielded)", "Quantified Reliability"),
  Black_Swan_Effect = c("Outliers wreck the Mean", "Overfits to noise", "72% Failure Risk", "High Discovery Risk"),
  Fix_Action = c("Masks 3-4 Sigma events", "Shrinkage toward truth", "LLN Stabilizes SNR", "Phase II Calibration")
)

audit_table %>%
  kbl(caption = "Engineering Audit: Robust vs. Standard Modeling") %>%
  kable_styling(full_width = F, position = "left") %>%
  column_spec(4, bold = T, color = "white", background = "#0072B2")
  
write.csv(audit_table, "data/audit_table.csv")
## Phase II Calibration End Game 
#  N = 42 Phase II data is processed through a Robust Hierarchical Model.

# The Decision Matrix: * If the 95% HDI of the difference clears zero: Production/Publication.
# If the 95% HDI overlaps zero: Kill the project. At $N=42$, if you can't see the signal, it's not there.

prior_audit <- data.frame(
  Parameter = c("mu_grand", "tau_grand", "sig[j]", "mu[j]", "nu (df)"),
  Value_JAGS = c("dnorm(5, 0.01)", "dgamma(0.1, 0.1)", "dunif(0, 5)", "dnorm(mu_grand, tau_grand)", "Constant = 4"),
  Strategy = c("Weakly Informative", "Non-committal Tether", "Flat Boundary", "Hierarchical / Shrinkage", "Robust / Fat-Tailed"),
  Engineering_Logic = c("Centered on expected plant height; high variance.", "Allows group means to diverge if signal is strong.", "Allows heteroscedasticity; caps noise at 100% CV.", "Tethers Treatment 1 to the 'Family Mean' to stop hype.", "Forces model to expect 3-4 sigma Black Swans.")
)

prior_audit %>%
  kbl(caption = "Phase II: Robust Hierarchical Prior Selection") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>%
  row_spec(2, bold = T, color = "white", background = "#0072B2") # Highlights the Tether

write.csv(prior_audit, "data/prior_audit.csv")


mod_string_robust = "model{

# 1. Hyper-priors (The Global Level)
mu_grand ~ dnorm(5, 0.01)       # SD of 10; very wide
tau_grand ~ dgamma(0.1, 0.1)    # Flat Gamma; weak shrinkage

# 2. Group Parameters (The Hierarchical Level)
for (j in 1:n_groups) {
  mu[j] ~ dnorm(mu_grand, tau_grand) # Hierarchical Tether
  
  # Heteroscedasticity handling
  tau[j] <- pow(sig[j], -2)
  sig[j] ~ dunif(0, 5)               # Noise floor ceiling
}

# 3. Likelihood (The Data Level)
for (i in 1:N) {
  # Robust Student-t Likelihood
  # Nu=4 ensures the model isn't 'surprised' by side-effect outliers
  y[i] ~ dt(mu[grp[i]], tau[grp[i]], 4) 
}
  
}"

robust_bayesian_audit <- function(mod_string, data, params, inits, n_iter, model_name = "model") {
  
  # --- 1. The Plotting Engine (Corrected) ---
  generate_plots <- function(sim_obj, y_hat, residuals, name) {
    
    # varnames() is the standard for mcmc.list objects
    all_params <- varnames(sim_obj)
    target_pars <- grep("^(mu|sig)\\[", all_params, value = TRUE)
    
    # Use regex to grab anything starting with mu or sig followed by brackets
    # This ensures mu[1], sig[2], etc. are all captured
    # A. MCMC Trace
    p1 <- mcmc_trace(sim_obj, pars = target_pars) +
      scale_color_brewer(palette = "Blues") +
      theme_minimal() +
      labs(title = paste("MCMC Trace Diagnostics:", name))
    ggsave(paste0("Figures/RobusttraceBS2.png"), plot = p1, width = 8, height = 6)
    
    # B. Posterior Density (Only for Means)
    mu_pars <- grep("^mu\\[", target_pars, value = TRUE)
    p2 <- mcmc_areas(sim_obj, pars = mu_pars, prob = 0.95) +
      scale_color_brewer(palette = "Blues") +
      theme_minimal() +
      labs(title = paste("Posterior Distributions:", name),
           subtitle = "95% Credible Intervals highlighted")
    ggsave(paste0("Figures/RobustdensityBS2.png"), plot = p2, width = 8, height = 6)
    
    # ... rest of your residual plotting code ...
    
    # C. Residuals (The Engineering Check)
    res_data <- data.frame(Predicted = y_hat, Residual = residuals)
    p3 <- ggplot(res_data, aes(x = Predicted, y = Residual)) +
      theme_minimal() + # Call theme first
      geom_point(size = 3, alpha = 0.6, color = "steelblue") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      geom_smooth(method = "loess", color = "darkblue", se = FALSE, linewidth = 1) +
      theme(plot.title = element_text(size = 14, face = "bold")) + 
      labs(title = paste("Residual Diagnostics:", name),
           subtitle = "Target: Random scatter with no Loess curvature")
    ggsave(paste0("Figures/RobustresidsBS2.png"), plot = p3, width = 8, height = 5)
  
  }
  
  # --- 2. JAGS Implementation ---
  mod <- jags.model(textConnection(mod_string), data = data, inits = inits, n.chains = 3, quiet = TRUE)
  update(mod, 5000) # Increased burn-in for the Student-t fat tails
  
  sim <- coda.samples(model = mod, variable.names = params, n.iter = n_iter)
  csim <- as.mcmc(do.call(rbind, sim))
  
  # 3. Quality Audit
  quality_report <- list(
    r_hat = gelman.diag(sim)$psrf,
    ess   = effectiveSize(sim),
    DIC   = dic.samples(mod, n.iter = 2000, type = "pD")
  )
  
  # 4. Robust Analytics (Using Medians for 'yhat')
  # This ensures 4-sigma side effects don't bias your 'Predicted' values
  pm_params <- apply(csim, 2, median) 
  mu_indices <- grep("mu\\[", names(pm_params))
  yhat <- pm_params[mu_indices][data$grp]
  resids <- data$y - yhat
  
  # 5. Engineering Metrics
  prob_better_1 <- mean(csim[, "mu[2]"] > csim[, "mu[1]"])
  hpd_diff_1    <- HPDinterval(as.mcmc(csim[, "mu[2]"] - csim[, "mu[1]"]))
  
  prob_better_2 <- mean(csim[, "mu[3]"] > csim[, "mu[1]"])
  hpd_diff_2    <- HPDinterval(as.mcmc(csim[, "mu[3]"] - csim[, "mu[1]"]))
  
  # Plot
  generate_plots(sim, yhat, resids, model_name)
  
  return(list(
    sim = sim, 
    csim = csim,
    quality_report = quality_report,
    stats = list(
      prob_b1 = prob_better_1, hpd_b1 = hpd_diff_1,
      prob_b2 = prob_better_2, hpd_b2 = hpd_diff_2
    )
  ))
}

# Define Initial Values for the Robust Model
# Calibrated Inits based on Audit 2 Realities
inits_robust <- list(
  # Chain 1: Centered on Audit 2 Results (The "Realistic" Start)
  list(mu = c(5.03, 4.66, 5.53), sig = c(0.80, 0.93, 0.73), mu_grand = 5.0), 
  
  # Chain 2: Low-end "Pessimistic" (Tests the lower bounds)
  list(mu = c(4.50, 4.00, 5.00), sig = c(1.20, 1.30, 1.10), mu_grand = 4.5), 
  
  # Chain 3: High-end "Optimistic" (Tests if outliers pull the chains)
  list(mu = c(5.50, 5.50, 6.00), sig = c(0.50, 0.60, 0.50), mu_grand = 5.5)
)

# Data set Generation 
#audit2_sum
# --- 1. Calibrated Parameters from Audit 2 ---
n_per_group <- 42

# Means (Posterior Means from your Audit 2)
mu_ctrl <- 5.03   # mu[1]
mu_trt1 <- 4.66   # mu[2] (The Underperformer/High Risk)
mu_trt2 <- 5.53   # mu[3] (The Clear Asset)

# Noise (Posterior Means from your Audit 2)
sig_ctrl <- 0.80  # sig[1]
sig_trt1 <- 0.93  # sig[2] (The High Entropy Generator)
sig_trt2 <- 0.73  # sig[3] (Lower Noise/More Stable)

# --- 2. Dynamic Data Generation ---
set.seed(Sys.time()) 

y_ctrl <- rnorm(n_per_group, mean = mu_ctrl, sd = sig_ctrl)
y_trt1 <- rnorm(n_per_group, mean = mu_trt1, sd = sig_trt1)
y_trt2 <- rnorm(n_per_group, mean = mu_trt2, sd = sig_trt2)

# --- 3. The Black Swan Injector (Specifically for the High-Entropy Trt1) ---
# We inject outliers here to see if the Student-t model can 
# recover the 'True' underperformance without being skewed by freak highs.
# Induce outliers in the 'Asset' group to test model resilience
outlier_indices_trt2 <- sample(1:n_per_group, 3)
outlier_indices_tr1 <- sample(1:n_per_group, 3)
y_trt1[outlier_indices_tr1] <- c(9.5, 1.2, 8.8) 
# We use the same 'Black Swan' values (9.5, 1.2, 8.8) to maintain parity
y_trt2[outlier_indices_trt2] <- c(9.5, 1.2, 8.8)

# --- 4. Package for JAGS ---
df_synth <- data.frame(
  growth = c(y_ctrl, y_trt1, y_trt2),
  grp    = rep(1:3, each = n_per_group)
)
print(df_synth)
data_list_robust <- list(
  y        = df_synth$growth,
  grp      = df_synth$grp,
  N        = nrow(df_synth),
  n_groups = 3
)
print(data_list_robust)
# 1. Run Robust Audit 
robust_audit <- robust_bayesian_audit( 
  mod_string = mod_string_robust, 
  data       = data_list_robust, 
  params     = c("mu", "sig"), 
  inits      = inits_robust, 
  n_iter     = 20000, 
  model_name = "PhaseII_Robust_Audit" 
)

saveRDS(robust_audit, "data/robust_auditBS2.rds")

robust_audit_sumBS2 <- summary(robust_audit$sim)
Robust_tableBS2 <- cbind(robust_audit_sumBS2$statistics, robust_audit_sumBS2$quantiles) 
print(Robust_tableBS2) 
write.csv(Robust_tableBS2, "data/Robust_table_summaryBS2.csv")

# DIC 
robust_audit_dic2 <- robust_audit$quality_report$DIC 
print(robust_audit_dic2) # 352.8 #390

# --- 1. Extract Metrics ---

m1_dic_val <- sum(audit1$quality_report$DIC$deviance) + sum(audit1$quality_report$DIC$penalty)
m2_dic_val <- sum(audit2$quality_report$DIC$deviance) + sum(audit2$quality_report$DIC$penalty)
m3_dic_val <- sum(robust_audit$quality_report$DIC$deviance) + sum(robust_audit$quality_report$DIC$penalty)

# --- 2. Build the Comprehensive Comparison Table ---
model_comparison_final <- data.frame(
  Model = c("Mod 1: Pooled", "Mod 2: Separate", "Mod 3: Robust-t"),
  Likelihood = c("Normal", "Normal", "Student-t (v=4)"),
  Variance_Structure = c("Homoscedastic", "Heteroscedastic", "Heteroscedastic"),
  N_per_Group = c(10, 10, 42), # Explicitly noting the Phase transition
  Penalty_pD = c(sum(audit1$quality_report$DIC$penalty), 
                 sum(audit2$quality_report$DIC$penalty), 
                 sum(robust_audit$quality_report$DIC$penalty)),
  DIC = c(m1_dic_val, m2_dic_val, m3_dic_val),
  Verdict = c("Discarded (Bias)", "Discarded (Noisy)", "SELECTED (Winner)")
)

ModelCompOpBS2 <- model_comparison_final %>%
  kable(
    format = "html", 
    digits = 2, 
    caption = "Phase I (Pilot) vs Phase II (Stress Test) Model Evolution",
    col.names = c("Model Architecture", "Likelihood Function", "Variance Logic", "Total N", "Complexity (pD)", "DIC Score", "Audit Verdict")
  ) %>% 
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "bordered"),
    full_width = F,
    position = "left"
  ) %>% 
  row_spec(0, bold = TRUE, color = "white", background = "#2c3e50") %>% # Header styling
  row_spec(3, bold = TRUE, background = "#d5f5e3") %>%                # Highlight the Winner
  column_spec(7, italic = TRUE, color = spec_color(1:3, end = 0.5, option = "viridis", direction = -1)) # Fade out discarded models

save_kable(ModelCompOpBS2, "Figures/Model_Comparison_BS2Audit.png")
print(ModelCompOpBS2)
write.csv(ModelCompOpBS2, "data/FinalModelComparisonBS2.csv")

# sensitivity Analysis 

# Extract Robust Stats
robust_statsBS2 <- summary(robust_audit$sim)$statistics[grep("mu\\[", rownames(summary(robust_audit$sim)$statistics)), 1:2]
robust_ciBS2    <- summary(robust_audit$sim)$quantiles[grep("mu\\[", rownames(summary(robust_audit$sim)$quantiles)), c(1, 5)]

# Build the 4-way Plotting DF
plot_df_4way <- data.frame(
  Group  = rep(c("Control", "Trt1", "Trt2"), 4),
  Method = rep(c("Freq (lm)", "Bayes (Mod 1)", "Bayes (Mod 2)", "Bayes (Robust)"), each = 3),
  Mean   = c(freq_stats[,1], bayesian_stats_1[,1], bayesian_stats_2[,1], robust_statsBS2[,1]),
  Lower  = c(freq_ci[,1], bayesian_ci1[,1], bayesian_ci2[,1], robust_ciBS2[,1]),
  Upper  = c(freq_ci[,2], bayesian_ci1[,2], bayesian_ci2[,2], robust_ciBS2[,2])
)

# Revised Plot
Sensitivity_Final <- ggplot(plot_df_4way, aes(x = Group, y = Mean, color = Method)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                position = position_dodge(width = 0.7), width = 0.3) +
  theme_minimal() +
  scale_color_viridis_d(option = "mako", end = 0.8) +
  labs(title = "Final System Audit: Parameter Stability Across 4 Methods",
       subtitle = "Robust Model N=42 proves maximum precision for Trt2,| Black Swan in Trt1, Tr2")
print(Sensitivity_Final)

ggsave("Figures/Final_System_AuditBS2.png", plot = Sensitivity_Final, width = 10, height = 6)


n_burn <- 5000
n_iter <- 20000

# 1. Run Audit 1 (Pooled Variance)
audit1_42 <- run_bayesian_audit(
  mod_string = mod_string_1, 
  data       = data_list_robust, 
  params     = c("mu", "sig"), 
  inits      = inits1, 
  n_burn     = 5000,
  n_iter     = 20000, 
  model_name = "Audit1_Pooled_42"
)

# 2. Run Audit 2 (Separate Variance)
audit2_42 <- run_bayesian_audit(
  mod_string = mod_string_2, 
  data       = data_list_robust, 
  params     = c("mu", "sig"), 
  inits      = inits2, 
  n_burn     = 5000,
  n_iter     = 20000, 
  model_name = "Audit2_Separate_42"
)

# 3. Save the 'Final Objects' for future reuse

saveRDS(audit1_42, "data/audit1_42_final.rds")
saveRDS(audit2_42, "data/audit2_42_final.rds")

# 4. Export the Master Summary Table for the 1-Pager 
audit1_42sum  <- summary(audit1_42$sim) 
final_table_142 <- cbind(audit1_42$statistics, audit1_42sum$quantiles) 
write.csv(final_table_142, "data/audit1_42summary.csv") 

audit2_42sum <- summary(audit2_42$sim) 
final_table_242 <- cbind(audit2_42$statistics, audit2_42sum$quantiles) 
write.csv(final_table_242, "data/audit2_42summary.csv") 


# 4.1 Bayesian Model 1 vs Model 2  with N = 42

# Using DIC Criterion 

audit1_42dic <- audit1_42$quality_report$DIC 
print(audit1_42dic) # 41

audit2_42dic <- audit2_42$quality_report$DIC 
print(audit2_42dic) 

# The Delta DIC 

Delta_42DIC <- audit1_42dic - audit2_42dic 
print(Delta_42DIC)

# 1. Extract the actual scalar values correctly
m1_42dic_val <- sum(audit1_42dic$deviance) + sum(audit1_42dic$penalty) 
m2_42dic_val <- sum(audit2_42dic$deviance) + sum(audit2_42dic$penalty) 

# 2. Rebuild the DF with the correct math N = 42 
model_comp42 <- data.frame(
  Model = c("Pooled Variance (Audit 1, N = 42)", "Separate Variance (Audit 2, N = 42)"),
  Mean_Deviance = c(sum(audit1_42dic$deviance), sum(audit2_42dic$deviance)),
  Penalty_pD = c(sum(audit1_42dic$penalty), sum(audit2_42dic$penalty)),
  DIC = c(m1_42dic_val, m2_42dic_val)
)

model_comp42$Delta_42DIC <- model_comp$DIC - min(model_comp$DIC)
knitr::kable(model_comp42, digits = 2) 
write.csv(model_comp42, "data/model_comparison42.csv")

#|#Model                       | Mean_Deviance| Penalty_pD|   DIC| Delta_DIC|
#|#:---------------------------|-------------:|----------:|-----:|---------:|
#|#Pooled Variance (Audit 1)   |         58.97|       4.10| 63.08|      0.00|
#|#Separate Variance (Audit 2) |         61.20|       5.69| 66.89|      3.82|


# Get Frequentist Estimates 
# 1. Ensure 'grp' is treated as a categorical factor
df_synth$grp <- as.factor(df_synth$grp)

# 2. Apply Cell Means Method (-1 removes the intercept)
lmod_audit_42 <- lm(growth ~ -1 + grp, data = df_synth)

# 3. Extract Estimates for the Sensitivity Matrix
freq_stats_42 <- summary(lmod_audit_42)$coefficients[, 1:2] # Mean and SE
freq_ci_42    <- confint(lmod_audit_42)                    # 95% CI

# To ensure the JAGS engine was correctly initialized and that the choice of priors did not introduce 
# unintended bias, a Frequentist Cell Means Model was established as a baseline calibrator.

# The Cell Means formulation is preferred for Bayesian validation because it estimates the absolute 
# mean of each experimental group directly. This allows for a literal comparison between the Frequentist 
# 95% Confidence Intervals and the Bayesian 95% High Posterior Density (HPD) intervals.

# Technical Verdict: By using the Cell Means model as a calibrator, we prove that the 93.8% probability 
# of success for Treatment 2 is a robust signal derived from the data, not a byproduct of model parameterization 
# or prior "bullying."

# 2. Get Bayesian Estimates (from your Audit2_Separate results)

bayesian_stats_142 <- audit1_42sum$statistics[grep("mu", rownames(audit1_42sum$statistics)), 1:2]
bayesian_ci142 <- audit1_42sum$quantiles[grep("mu", rownames(audit1_42sum$quantiles)), c(1, 5)]
print(bayesian_ci142)

# Assuming 'audit2_sum' is the summary(audit2$sim) 
bayesian_stats_242 <- audit2_42sum$statistics[grep("mu", rownames(audit2_42sum$statistics)), 1:2]
bayesian_ci242 <- audit2_42sum$quantiles[grep("mu", rownames(audit2_42sum$quantiles)), c(1, 5)]
print(bayesian_ci242) 

# 3. Build the Validation Matrix  
validation_matrix42 <- data.frame(
  Parameter = c("Control", "Trt 1", "Trt 2"),
  Freq_Mean = freq_stats_42[,1], 
  # Changed index from [,2] (SD) to [,1] (Mean)
  Bayes_Mean_M1 = bayesian_stats_142[,1], 
  Bayes_Mean_M2 = bayesian_stats_242[,1], 
  Freq_95_CI = paste0("[", round(freq_ci42[,1], 2), ", ", round(freq_ci42[,2], 2), "]"), 
  Bayes_95_HPD_M1 = paste0("[", round(bayesian_ci142[,1], 2), ", ", round(bayesian_ci142[,2], 2), "]"),
  Bayes_95_HPD_M2 = paste0("[", round(bayesian_ci242[,1], 2), ", ", round(bayesian_ci242[,2], 2), "]") 
  # Comma removed here to fix the "missing argument" error
)

knitr::kable(validation_matrix42, 
             digits = 3, 
             caption = "System Calibration: Frequentist (Cell Means) vs. Bayesian Models")
write.csv(validation_matrix42, "data/Validation_matrix.csv")

# 2. Combine into one data frame
plot_df_3way42 <- data.frame(
  Group  = rep(c("Control", "Trt1", "Trt2"), 3),
  Method = rep(c("Frequentist (lm)", "Bayesian (Mod 1)", "Bayesian (Mod 2)"), each = 3),
  Mean   = c(freq_stats_42[,1], bayesian_stats_142[,1], bayesian_stats_242[,1]),
  Lower  = c(freq_ci_42[,1], bayesian_ci142[,1], bayesian_ci242[,1]),
  Upper  = c(freq_ci_42[,2], bayesian_ci242[,2], bayesian_ci242[,2])
) 
# 3. Plot it  
SensitivityAnalysis42 <- ggplot(plot_df_3way42, aes(x = Group, y = Mean, color = Method)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                position = position_dodge(width = 0.6), width = 0.2) +
  theme_minimal() +
  labs(title = "Sensitivity Analysis: Weight Estimates Across Models",
       subtitle = "Comparing Pooled vs. Separate Variance Assumptions")

print(SensitivityAnalysis42)
ggsave(paste0("Figures/3WSensitivity_Analysis42.png"), plot = SensitivityAnalysis42, width = 8, height = 5)


# Extract Robust Stats
robust_statsBS2 <- summary(robust_audit$sim)$statistics[grep("mu\\[", rownames(summary(robust_audit$sim)$statistics)), 1:2]
robust_ciBS2    <- summary(robust_audit$sim)$quantiles[grep("mu\\[", rownames(summary(robust_audit$sim)$quantiles)), c(1, 5)]

# Build the 4-way Plotting DF
plot_df_4way42 <- data.frame(
  Group  = rep(c("Control", "Trt1", "Trt2"), 4),
  Method = rep(c("Freq (lm)", "Bayes (Mod 1)", "Bayes (Mod 2)", "Bayes (Robust)"), each = 3),
  Mean   = c(freq_stats_42[,1], bayesian_stats_142[,1], bayesian_stats_242[,1], robust_statsBS2[,1]),
  Lower  = c(freq_ci_42[,1], bayesian_ci142[,1], bayesian_ci242[,1], robust_ciBS2[,1]),
  Upper  = c(freq_ci_42[,2], bayesian_ci142[,2], bayesian_ci242[,2], robust_ciBS2[,2])
)

# Revised Plot 

Sensitivity_Final42 <- ggplot(plot_df_4way42, aes(x = Group, y = Mean, color = Method)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                position = position_dodge(width = 0.7), width = 0.3) +
  theme_minimal() +
  scale_color_viridis_d(option = "mako", end = 0.8) +
  labs(title = "Final System Audit: Parameter Stability Across 4 Methods",
       subtitle = "Robust Model N=42 proves maximum precision for Trt2,| Black Swan in Trt1, Tr2 ")
print(Sensitivity_Final42)
ggsave("Figures/Final_System_AuditBS242.png", plot = Sensitivity_Final42, width = 10, height = 6)

# conclusion 
# Outlier Resistance: While Frequentist and standard Bayesian models allow the "Black Swans" to pull the means 
# upward (see Trt2), the Robust (blue) model anchors the mean to the true signal, effectively ignoring the noise.

# Precision Optimization: The Robust-t model maintains the tightest credible intervals in the presence of chaos, 
# proving it yields the highest certainty for $N=42$.

# Asset Validation: By successfully shielding Trt2 from induced outliers, the audit confirms Trt2 as a stable 
# statistical asset that outperforms the Control regardless of environmental "freak" events.


# ##########################################

## Posterior Predictive Check 
# 1. Posterior Predictive Check: Generating new data from the model 
df_synth$grp <- as.factor(df_synth$grp) 
# Declare y and grp from your source data 
y_42   <- df_synth$growth  # The raw plant weights 
grp_42 <- as.numeric(df_synth$grp) # The group indices (1, 2, 3) 

# 'y_pred' will now simulate the 'messiness' accurately 

n_sim = 2000 # Increased for smoother histograms  

y_pred_42 = matrix(NA, nrow = n_sim, ncol = length(data_list_robust$y)) 
# PPC Model 1
ppc_m142 = matrix(NA, n_sim, length(y_42)) 
# PPC Model 2 
ppc_m242 = matrix(NA, n_sim, length(y_42)) 

for (i in 1:n_sim) {
  # Model 1: Pooled (Single Sig) 
  ppc_m142[i,] = rnorm(length(y_42), audit1_42$csim[i, 1:3][grp_42], audit1_42$csim[i, "sig"]) 
  
  # Model 2: Separate (Indexed Sig) 
  mu2 = audit2_42$csim[i, grep("mu", colnames(audit2_42$csim))]
  sig2 = audit2_42$csim[i, grep("sig", colnames(audit2_42$csim))]
  ppc_m242[i,] = rnorm(length(y_42), mu2[grp_42], sig2[grp_42])
}

# The Litmus Test: Variance of Group 2 (Treatment 1)
obs_var_trt142 <- var(y_42[grp_42 == 2])
print(obs_var_trt142) #2.27

pred_var_m1_trt142  <- apply(ppc_m142[, grp_42 == 2], 1, var)

pred_var_m2_trt142  <- apply(ppc_m242[, grp_42 == 2], 1, var)

p_val_m1_trt142 <- mean(pred_var_m1_trt142 > obs_var_trt142) # 0.0435
p_val_m2_trt142 <- mean(pred_var_m2_trt142 > obs_var_trt142) # 0.4065 

# Variance of Treatment 2 
obs_var_trt242 <- var(y_42[grp_42 == 3])
print(obs_var_trt242) # 1.728 

pred_var_m1_trt242  <- apply(ppc_m142[, grp_42 == 3], 1, var)  

pred_var_m2_trt242  <- apply(ppc_m242[, grp_42 == 3], 1, var) 

p_val_m1_trt242 <- mean(pred_var_m1_trt242 > obs_var_trt242) # 0.0435  
p_val_m2_trt242 <- mean(pred_var_m2_trt242 > obs_var_trt242) # 0.4065 
print(p_val_m1_trt242)  
print(p_val_m2_trt242) 

# Variance of Treatment 1 in Model 3 robust_audit_sumBS2

ppc_m342 = matrix(NA, n_sim, length(y_42))

for (i in 1:n_sim) {
  mu3  = robust_audit$csim[i, grep("mu", colnames(robust_audit$csim))]
  sig3 = robust_audit$csim[i, grep("sig", colnames(robust_audit$csim))]
  # Student-t simulation: mu + sigma * rt
  ppc_m342[i,] = mu3[grp_42] + sig3[grp_42] * rt(length(y_42), df = 4)
}
# model 3 trt 1 
pred_var_m3_trt1 <- apply(ppc_m342[, grp_42 == 2], 1, var)
p_val_m3_trt1    <- mean(pred_var_m3_trt1 > obs_var_trt242)
print(p_val_m3_trt1) #0.6625

# model 3 trt 2 
pred_var_m3_trt2 <- apply(ppc_m342[, grp_42 == 3], 1, var)
p_val_m3_trt2    <- mean(pred_var_m3_trt2 > obs_var_trt242)
print(p_val_m3_trt2) #0.366 

ppc_final_audit <- data.frame(
  Metric = c("Bayesian p-value (Trt1 Var)", "Bayesian p-value (Trt2 Var)", "Model Integrity"),
  Model_1_Pooled = c(round(p_val_m1_trt142, 3), round(p_val_m1_trt242, 3), "FAIL: Over-smoothed"),
  Model_2_Separate = c(round(p_val_m2_trt142, 3), round(p_val_m2_trt242, 3), "PASS: High Variance"),
  Model_3_Robust = c(round(mean(apply(ppc_m342[, grp_42 == 2], 1, var) > obs_var_trt142), 3), 
                     round(p_val_m3_trt2, 3), "WINNER: Resilient")
)

ppc_final_audit %>%
  kbl(caption = "Posterior Predictive Audit: Group-Specific Variance Recovery") %>%
  kable_styling(bootstrap_options = c("striped", "bordered")) %>%
  column_spec(4, bold = T, background = "#d5f5e3") # Highlight the Robust Model

# --- 0. Precise Dimension Locking ---##################start here tmrw
library(tidyverse)
library(patchwork)

# --- 0. Precise Dimension Locking ---
n_to_plot <- 100  # Sub-sample for spaghetti plot to keep it readable
df_obs_42 <- data.frame(weight = as.numeric(df_synth$growth))

# Helper to tidy PPC matrices
tidy_ppc <- function(ppc_matrix, model_name, n) {
  as.data.frame(ppc_matrix[1:n, ]) %>%
    mutate(sim_id = row_number(), model = model_name) %>%
    pivot_longer(cols = starts_with("V"), names_to = "plant_idx", values_to = "weight")
}

# --- 1. Robust Data Prep ---
sim_m1 <- tidy_ppc(ppc_m142, "Model 1: Pooled", n_to_plot)
sim_m2 <- tidy_ppc(ppc_m242, "Model 2: Separate", n_to_plot)
sim_m3 <- tidy_ppc(ppc_m342, "Model 3: Robust-t", n_to_plot)

combined_sims <- bind_rows(sim_m1, sim_m2, sim_m3)

# --- 2. Patchwork Construction ---

# Plot A: Spaghetti PPC (The "Shape" Audit)
plot_a <- ggplot() + 
  geom_density(data = combined_sims, 
               aes(x = weight, group = interaction(sim_id, model), color = model), 
               alpha = 0.15, linewidth = 0.3) +
  geom_density(data = df_obs_42, aes(x = weight), 
               color = "red", linewidth = 1.1) +
  facet_wrap(~model) +
  scale_color_manual(values = c("Model 1: Pooled" = "gray60", 
                                "Model 2: Separate" = "skyblue", 
                                "Model 3: Robust-t" = "#2D708EFF")) +
  theme_minimal() + 
  theme(legend.position = "none", strip.text = element_text(face="bold")) +
  labs(title = "A: Global Density PPC (N=42)", 
       subtitle = "Red: Observed Data | Faint Lines: Posterior Predictive Simulations", 
       x = "Plant Weight")

# Plot B: Variance Audit (The "Noise" Audit for Trt 1)
var_df <- data.frame(
  Variance = c(pred_var_m1_trt142, pred_var_m2_trt142, pred_var_m3_trt1),
  Model = rep(c("Model 1", "Model 2", "Model 3"), each = length(pred_var_m1_trt142))
)

plot_b <- ggplot(var_df, aes(x = Variance, fill = Model)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  geom_vline(xintercept = obs_var_trt142, color = "black", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("Model 1" = "#E41A1C", "Model 2" = "#377EB8", "Model 3" = "#4DBBD5FF")) +
  theme_minimal() + 
  labs(title = "B: Treatment 1 Variance Audit", 
       subtitle = "Dashed Line: Observed Variance (Target)",
       x = "Simulated Variance", y = "Frequency")

# Plot c: Variance Audit (The "Noise" Audit for Trt 2)
var2_df <- data.frame(
  Variance = c(pred_var_m1_trt242, pred_var_m2_trt242, pred_var_m3_trt2),
  Model = rep(c("Model 1", "Model 2", "Model 3"), each = length(pred_var_m1_trt242))
)

plot_c <- ggplot(var2_df, aes(x = Variance, fill = Model)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 50) +
  geom_vline(xintercept = obs_var_trt242, color = "black", linetype = "dashed", linewidth = 1) +
  scale_fill_manual(values = c("Model 1" = "#E41A1C", "Model 2" = "#377EB8", "Model 3" = "#4DBBD5FF")) +
  theme_minimal() + 
  labs(title = "C: Treatment 2 Variance Audit", 
       subtitle = "Dashed Line: Observed Variance (Target)",
       x = "Simulated Variance", y = "Frequency")


# Combine and Save
# Force both histograms to use the same X-axis for an honest comparison
final_audit_plot <- (plot_a / plot_b / plot_c) & xlim(0, 15)
print(final_audit_plot)

# Final Patchwork Assembly
final_production_plot <- (plot_a / plot_b / plot_c) +
  plot_layout(heights = c(1.2, 1, 1)) & 
  coord_cartesian(xlim = c(0, 15)) & # The "Goldilocks" Scale
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", colour = "darkblue"),
    panel.grid.minor = element_blank()
  )

# Adding the Global Labels
final_production_plot <- final_production_plot +
  plot_annotation(
    title = "Bayesian System Audit: Stress-Testing Outlier Resilience (N=42)",
    subtitle = "Evaluation of Pooled, Separate, and Robust-t Likelihoods under Black Swan Events",
    caption = "Data: Synthetic Plant Growth | Outliers induced at x=9.5, 1.2, 8.8 | Model 3 uses Nu=4"
  )
print(final_production_plot)
# Export for Report
ggsave("Figures/PhaseII_Final_Audit_N42.png", 
       plot = final_production_plot, 
       width = 11, height = 8.5, dpi = 300)

##########

# --- 1. Density Overlays (The "Shape" Audit) ---
# Check if the generated 'universes' overlap the real red data line
p1_m1 <- ppc_dens_overlay(y_42, ppc_m142[1:50, ]) + ggtitle("M1: Pooled (Normal)")
p1_m2 <- ppc_dens_overlay(y_42, ppc_m242[1:50, ]) + ggtitle("M2: Separate (Normal)")
p1_m3 <- ppc_dens_overlay(y_42, ppc_m342[1:50, ]) + ggtitle("M3: Robust (Student-t)")

# --- 2. Variance Stat Checks (The "Entropy" Audit) ---
# Check if the model 'understands' the total variation including Black Swans
p2_m1 <- ppc_stat(y_42, ppc_m142, stat = "var") + ggtitle("M1 Global Var")
p2_m2 <- ppc_stat(y_42, ppc_m242, stat = "var") + ggtitle("M2 Global Var")
p2_m3 <- ppc_stat(y_42, ppc_m342, stat = "var") + ggtitle("M3 Global Var")

# --- 3. Composite Assembly ---
global_audit <- (p1_m1 | p1_m2 | p1_m3) / (p2_m1 | p2_m2 | p2_m3) & 
  theme_minimal() & 
  coord_cartesian(xlim = c(0, 15)) # The Scale Lock

print(global_audit)
ggsave("Figures/global_audit_6panel.png", plot = global_audit, width = 12, height = 7)


library(ggdist)
library(tidybayes)
#### Signal Recovery under Outlier Stress (posterior distribution)

# Extract samples and label them properly
# Assume mu[1] = Control, mu[2] = Trt 1, mu[3] = Trt 2
extract_mu <- function(audit_obj, model_name) {
  as.data.frame(audit_obj$csim[, c("mu[1]", "mu[2]", "mu[3]")]) %>%
    setNames(c("1", "2", "3")) %>% # Match these to your df_synth$grp levels
    pivot_longer(cols = everything(), names_to = "Group", values_to = "mu_val") %>%
    mutate(Model = model_name)
}

raincloud_data <- bind_rows(
  extract_mu(audit1_42, "Pooled"),
  extract_mu(audit2_42, "Separate"),
  extract_mu(robust_audit, "Robust")
)

ggplot(raincloud_data, aes(x = Group, y = mu_val, fill = Model)) +
  # 1. THE CLOUDS (Posterior Expectations)
  # These show the uncertainty of the MEAN (mu). 
  # If a cloud is "tight," the model is very sure about the average.
  stat_halfeye(point_interval = median_qi, .width = c(.95), 
               position = position_dodge(width = 0.6), alpha = 0.7) +
  
  # 2. THE RAIN (Raw Observations)
  # These show the actual data points from df_synth.
  # The distance between the 'Rain' and the 'Cloud' shows the Outlier Leverage.
  geom_point(data = df_synth, aes(x = factor(grp), y = growth), 
             inherit.aes = FALSE, alpha = 0.4, size = 1.5,
             position = position_jitter(width = 0.1, height = 0)) +
  
  # Formatting
  scale_fill_manual(values = c("Pooled" = "#999999", "Separate" = "#56B4E9", "Robust" = "#0072B2")) +
  labs(
    title = "Raincloud Audit: Signal Recovery Under Outlier Stress",
    subtitle = "Clouds: Posterior distributions of Group Means (mu) | Rain: Raw observations (y)",
    x = "Experimental Group (1: Control, 2: Trt1, 3: Trt2)",
    y = "Growth Metric (Units)",
    fill = "Model Architecture"
  ) +
  theme_minimal()
ggsave("Figures/Raincloud_SignalRecovery.png", plot = global_audit, width = 12, height = 7)

### Noise profile 
# Extract Sigma for Treatment 1 (Group 2)
# Model 1 uses one 'sig' for all; Model 2 uses 'sig[2]' for Trt1
sig_m1_42 <- audit1_42$csim[, "sig"] 
sig_m2_trt1 <- audit2_42$csim[, "sig[2]"] 
sig_m3_trt1 <- robust_audit$csim[, "sig[2]"]
sig_m3_ctrl <- median(robust_audit$csim[, "sig[1]"]) # Control
print(sig_m3_ctrl) 
# Treatment 2 
sig_m2_trt2 <- audit2_42$csim[, "sig[3]"] 
sig_m3_trt2 <- robust_audit$csim[, "sig[3]"]

# Summary Statistics 
summary_noise_42 <- data.frame( 
  Model = c("Model 1 (Pooled)", "Model 2 (Trt1)", "Robust (Trt1)","Model 2 (Trt2)", "Robust (Trt2)"),
  Median_Sigma = c(median(sig_m1_42), median(sig_m2_trt1), median(sig_m3_trt1), median(sig_m2_trt2), median(sig_m3_trt2)),
  Lower_95 = c(quantile(sig_m1_42, 0.025), quantile(sig_m2_trt1, 0.025), quantile(sig_m3_trt1, 0.025), quantile(sig_m2_trt2, 0.025), quantile(sig_m3_trt2, 0.025)),
  Upper_95 = c(quantile(sig_m1_42, 0.975), quantile(sig_m2_trt1, 0.975), quantile(sig_m3_trt1, 0.975), quantile(sig_m2_trt2, 0.975), quantile(sig_m3_trt2, 0.975))
)
print(summary_noise_42)
write.csv(summary_noise_42, "data/noise_42.csv")

# Calculate the 'Underestimation Percentage'
gap_M1 = (median(sig_m2_trt1) - median(sig_m1_42)) / median(sig_m1_42) * 100
cat("Model 1 is hiding approximately", round(gap_M1, 1), "% of the noise in Treatment 1.")
# Model 1 is hiding approximately 19 % of the noise in Treatment 1.

# Percentage difference within Model 2
internal_diff_M2 <- (median(sig_m2_trt1) - median(sig_m2_trt2)) / median(sig_m2_trt2) * 100
print(internal_diff_M2)
cat("Model 2 identifies that Treatment 1 is", round(internal_diff_M2, 1), "% noisier than Treatment 2.")
# Model 2 identifies that Treatment 1 is 13.9 % noisier than Treatment 2.

# Internal difference of robust model 
internal_diff_M3 <- (median(sig_m3_trt1) - median(sig_m3_trt2)) / median(sig_m3_trt2) * 100
print(internal_diff_M3)
cat("Robust Model identifies that Treatment 1 is", round(internal_diff_M3, 1), "% noisier than Treatment 2.")

# Robust Model identifies that Treatment 1 is 20.9 % noisier than Treatment 2.
print(audit1_42$prob_better_1) # 0.02906667
print(audit1_42$prob_better_2) # 0.9919

print(audit2_42$prob_better_1) # 0.0265
print(audit1_42$prob_better_2) # 0.9919

robust_audit$stats
# prob_b1 = 0.004466667
# prob_b2 = 0.9976833

###  Probability that Trt 1 is noisier than Trt 2 
M1prob_trt1_noisier <- mean(audit1_42$csim[, "sig"] > audit1_42$csim[, "sig"])
cat("Model1 : Probability that Trt 1 > Trt 2 variance:", round(M1prob_trt1_noisier * 100, 2), "%\n")
# Probability that Trt1 > Trt2 Model assumes all groups are equally noisy 

M2prob_trt1_noisier <- mean(audit2_42$csim[, "sig[2]"] > audit2_42$csim[, "sig[3]"])
cat("Model2 : Probability that Trt 1 > Trt 2 variance:", round(M2prob_trt1_noisier * 100, 2), "%\n")
# Model2 : Probability that Trt 1 > Trt 2 variance: 81.05 %

M3prob_trt1_noisier <- mean(robust_audit$csim[, "sig[3]"] > robust_audit$csim[, "sig[2]"])
cat("Model3 : Probability that Trt 1 > Trt 2 variance:", round(M3prob_trt1_noisier * 100, 2), "%\n")
# Model3 : Probability that Trt 1 > Trt 2 variance: 81.26 %


# Posterior Probability Distribution of delta 
# Calculate the posterior difference
# 1. Calculate the raw posterior difference
# Ensure these are numeric vectors

library(bayestestR) # For hdi

# 1. Prepare Data from the WINNER (Robust Model)
# We use sig[2] (Trt1) and sig[3] (Trt2) from the robust_audit
sigma_delta_vec <- as.numeric(robust_audit$csim[, "sig[2]"] - robust_audit$csim[, "sig[3]"])
plot_df <- data.frame(sigma_delta = sigma_delta_vec)

# Stats for the Decision Box
prob_greater <- mean(sigma_delta_vec > 0) 
z_score_noise <- mean(sigma_delta_vec) / sd(sigma_delta_vec) 
hdi_limits <- hdi(sigma_delta_vec, credMass = 0.95)
median_diff <- median(sigma_delta_vec)

# 2. Build the Robust Annotated Plot
post_prob <- ggplot(plot_df, aes(x = sigma_delta)) +
  # The Distribution (The "Truth" about the Noise Gap)
  geom_density(fill = "#0072B2", alpha = 0.6) + # Dark blue for Robust
  
  # The HDI Segment (Visualizes the 95% certainty range)
  annotate("segment", x = hdi_limits$CI_low, xend = hdi_limits$CI_high, y = 0.05, yend = 0.05, 
           linewidth = 2, color = "black") +
  annotate("text", x = mean(c(hdi_limits$CI_low, hdi_limits$CI_high)), y = 0.2, 
           label = "95% Credible Interval", color = "black", fontface = "bold", size = 3.5) +
  
  # The Zero-Line (The "Plain Wrong" Assumption Line)
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
  
  # The Engineering Decision Box
  annotate("label", x = Inf, y = Inf, hjust = 1.1, vjust = 1.1, 
           label = paste0("ROBUST NOISE AUDIT (N=42):\n",
                          "Prob(Trt1 > Trt2 Noise): ", round(prob_greater * 100, 1), "%\n",
                          "Median Sigma Delta: +", round(median_diff, 3), "\n",
                          "Bayesian Z-score: ", round(z_score_noise, 3), "\n",
                          "Verdict: Structural Variance Confirmed\n",
                          "Action: Deploy Group-Specific Sigma"),
           fill = "white", alpha = 0.9, family = "mono", size = 3.5) +
  
  theme_minimal() +
  scale_x_continuous(limits = c(-0.5, 1.0), breaks = seq(-0.5, 1.0, 0.25)) +
  labs(title = "Structural Heteroscedasticity Audit: Robust-t Calibration",
       subtitle = "Posterior Difference in Sigma (Trt 1 - Trt 2) | Outliers Filtered",
       x = "Sigma Difference (Units)",
       y = "Posterior Density")

print(post_prob) 
ggsave("Figures/robust_noise_audit_N42.png", plot = post_prob, width = 8, height = 5)

# Final audit table 

library(dplyr)

# 1. Final Data Construction (Robust M3 Results)
# These values are pulled from your verified robust_audit$stats and csim
final_results_df <- data.frame(
  Group = c("Control", "Treatment 1", "Treatment 2"),
  Posterior_Mean = c(5.01, 4.39, 5.58),
  Posterior_SD = c(0.12, 0.20, 0.16),
  Lower_95_HDI = c(4.79, 4.00, 5.26),
  Upper_95_HDI = c(5.25, 4.80, 5.90),
  Noise_Floor = c(0.64, 1.09, 0.90), # Sigma (Robust)
  Prob_Success = c("-", "0.4%", "99.8%"), # P(Trt > Control)
  noisier_than  = c("-", "81.26%","18.74%"),
  Verdict = c("Baseline", "FAILURE: Inhibitor", "SUCCESS: High Yield")
)

# 2. Build the Kable
final_results_df %>%
  kbl(caption = "Final Engineering Audit: Robust Posterior Estimates and Risk Assessment (N=42)",
      col.names = c("Group", "Mean Growth (μ)","std Deviation (sd)", "95% Low", "95% High", "Noise (σ)", "P(Success)","Relative Noise", "Verdict"),
      align = "lcccccc") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = F) %>%
  # Highlight the Success/Failure for immediate scanning
  column_spec(1, bold = T) %>%
  column_spec(2, bold = T, background = "#F7F9F9") %>%
  column_spec(7, color = ifelse(final_results_df$Prob_Success == "99.8%", "green", "red"), bold = T) %>%
  column_spec(9, bold = T, color = "white", 
              background = case_when(
                final_results_df$Verdict == "SUCCESS: High Yield" ~ "#28B463",
                final_results_df$Verdict == "FAILURE: Inhibitor" ~ "#CB4335",
                TRUE ~ "#808B96"
              )) %>%
  add_footnote("Calculated using Student-t distribution (ν=4) to filter outliers. High Noise (σ) in Trt 1 correlates with poor signal stability.")
write.csv(final_results_df, "data/final_Results_Robust.csv")


# Black Swan detection (M1 + M2 + M3 )

# 1. Extract Residuals for all models
# Residual = Observed Data - Model Mean (mu)
audit1_42$csim[, c("mu[1]", "mu[2]", "mu[3]")]
audit1_42sum
calc_black_swans <- function(audit_obj, model_name) {
  # Get the posterior means for each group
  mu_samples <- audit_obj$csim[, c("mu[1]", "mu[2]", "mu[3]")]
  mu_means <- colMeans(mu_samples)
  
  
  # Map means to the raw data rows
  df_synth$expected_mu <- mu_means[df_synth$grp]
  # Calcualte Residual 
  df_synth$residual <- abs(df_synth$growth - df_synth$expected_mu)
  
  # Standardize residuals by the model's Sigma
  sigma_model = ifelse(model_name=="Model 1 (Pooled)",sigma_model <- median(audit_obj$csim[, "sig"]))
  #sigma_model <- median(audit_obj$csim[, "sig"]) # Simplified for M1; use group-sig for M2/M3
  df_synth$z_score <- df_synth$residual / sigma_model
  
  # A Black Swan is any point > 3 Standard Deviations away averaging model 1 
  df_synth$is_black_swan <- df_synth$z_score > 3
  
  return(df_synth %>% mutate(Model = model_name))
}
swan_audit <- bind_rows(
  calc_black_swans(audit1_42, "Model 1 (Pooled)"),
  #calc_black_swans(audit2_42, "Model 2 (Separate)"),
  #calc_black_swans(robust_audit, "Model 3 (Robust)")
)
# 3. Summarize the Detection Power
swan_summary <- swan_audit %>%
  group_by(Model) %>%
  summarise(
    Swans_Detected = sum(is_black_swan),
    Max_Z_Score = max(z_score),
    Avg_Residual = mean(residual)
  )
print(swan_summary)

  #Hardcoded
  # 1. Map the correct means and sigmas to each group
  # Using your verified N=42 values
  m1_map <- data.frame(grp = factor(1:3), mu_m1 = c(5.015, 4.499, 5.671), sig_m1 = c(1.242, 1.242, 1.242))
  m2_map <- data.frame(grp = factor(1:3), mu_m2 = c(5.01, 5.24, 5.68), sig_m2 = c(0.82, 1.47, 1.29))
  m3_map <- data.frame(grp = factor(1:3), mu_m3 = c(5.02, 4.39, 5.59), sig_m3 = c(0.65, 1.09, 0.90))
  
  
  # Final Audit Mapping
  swan_audit <- df_synth %>%
    left_join(m1_map, by = "grp") %>%
    left_join(m2_map, by = "grp") %>%
    left_join(m3_map, by = "grp") %>%
    mutate(
      Z_M1 = abs(growth - mu_m1) / sig_m1,
      Z_M2 = abs(growth - mu_m2) / sig_m2,
      Z_M3 = abs(growth - mu_m3) / sig_m3,
      Detection = case_when(
        # Condition 1: Only the Robust model is sensitive enough to catch it
        Z_M3 > 3 & Z_M2 <= 3 & Z_M1 <= 3 ~ "M3 Only (Success)",
        
        # Condition 2: High-intensity outliers caught by everyone
        Z_M3 > 3 & Z_M2 > 3 & Z_M1 >= 3  ~ "All Detected",
        
        # Condition 3: Caught by M2 and M3, but M1 was too blunt
        Z_M3 > 3 & Z_M2 >= 3 & Z_M1 <= 3 ~ "M2, M3 Detected",
        
        # Condition 4: If M3 says it's a swan but the others don't fit the above
        Z_M3 > 3                         ~ "M3 Detected (Misc)",
        
        # Default: It's just a normal data point
        TRUE                             ~ "Standard Point"
      )
    )
  
  # Identifying the "Invisible" Swan
  # A point that M1 and M2 miss, but M3 catches
  invisible_swans <- swan_audit %>% 
    filter(Z_M3 > 3 & Z_M2 < 3 & Z_M1 < 3)
  
  print(invisible_swans)
  
  # Standard Black Swans 
  Std_Black_Swans <- swan_audit %>%
    filter(Z_M3 > 3 & Z_M2 >3 & Z_M1>3)
  print(Std_Black_Swans)
  

  # Compare the 'Belief' of noise across models
  # Properly named data for the legend
  sigma_samples <- data.frame(
    `Model 1 (Pooled)` = audit1_42$csim[, "sig"],
    `Model 2 (Gaussian Trt1)` = audit2_42$csim[, "sig[2]"],
    `Model 3 (Robust Trt1)`   = robust_audit$csim[, "sig[2]"]
  ) %>% pivot_longer(cols = everything(), names_to = "Audit_Model", values_to = "Sigma")
  
 noise_gap <- ggplot(sigma_samples, aes(x = Sigma, fill = Audit_Model)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = 1.0, linetype = "dashed", color = "black") + 
    theme_minimal() +
    scale_fill_brewer(palette = "Set1") +
    labs(title = "Sigma Audit: The 'Noise' Detection Gap",
         subtitle = "Dashed line = True Noise Floor (1.0)",
         x = "Estimated Noise (Sigma)", 
         y = "Density",
         fill = "Model") # Legend Title Fix
  ggsave("Figures/NoiseDetGap.png", plot = noise_gap, width = 12, height = 7)

  # Final sigma plot 
  # The most powerful way to end this is a 2x3 Comparison Matrix:
  # Columns: Control, Trt 1, Trt 2.
  # Rows: Gaussian (M2) vs. Robust (M3).
  # This proves that in every single group, the Robust model found a "Lower Noise Floor" by successfully identifying and weight-adjusting the outliers.
  
  ### Outlier detection Sensitivity Plot 
  
  # 1. Create the Detection Sensitivity Data
  z_comparison <- swan_audit %>%
    select(growth, grp, Z_M2, Z_M3) %>%
    pivot_longer(cols = c(Z_M2, Z_M3), names_to = "Model", values_to = "Z") %>%
    mutate(Model = ifelse(Model == "Z_M2", "Gaussian (M2)", "Robust (M3)"))
  
  # 2. The Sigma-Sensitivity Plot
 Outlier_sens<- ggplot(z_comparison, aes(x = growth, y = Z, color = Model)) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_hline(yintercept = 3, linetype = "dashed", color = "red") + 
    facet_wrap(~grp, labeller = as_labeller(c(`1`="Control", `2`="Trt 1", `3`="Trt 2"))) +
    scale_color_manual(values = c("#E69F00", "#56B4E9")) +
    theme_minimal() +
    labs(title = "The Detection Gap: Why Robust Sigma Wins",
         subtitle = "Red Line = Black Swan Detection Threshold (3 Sigma)",
         x = "Observed Growth", y = "Detection Strength (Z-Score)")
  
 ggsave("Figures/OutlierSensitivity.png", plot = Outlier_sens, width = 12, height = 7)
  ########################################################################### 
  
  # Simulating outlier Universes (Winner Edition )
  
  library(metRology)

  # 1. Extract the Winner's Parameters (Trt 2)
  mu_w  <- mean(robust_audit$csim[, "mu[3]"])
  sig_w <- mean(robust_audit$csim[, "sig[3]"])
  nu_fixed <- 4 # Your structural assumption
  
  # 2. Simulate the 'Robust Universe' vs 'Gaussian Universe'
  set.seed(42)
  n_sim <- 10000
  sim_data <- data.frame(
    Growth = c(
      rnorm(n_sim, mu_w, sig_w),                       # What Model 2 thinks happens
      rt.scaled(n_sim, df = nu_fixed, mean = mu_w, sd = sig_w) # What Model 3 knows happens
    ),
    Universe = rep(c("Gaussian (Fragile)", "Robust (Swan-Ready)"), each = n_sim)
  )
  
  # 3. Plot the 'Tail Risk'
  Universe <- ggplot(sim_data, aes(x = Growth, fill = Universe)) +
    geom_density(alpha = 0.6) +
    geom_vline(xintercept = c(1.2, 8.8, 9.5), linetype = "dotted", color = "red", linewidth = 0.8) + 
    coord_cartesian(xlim = c(0, 12)) +
    theme_minimal() +
    labs(title = "Simulation: The Winner's Universe (Trt 2)",
         subtitle = "Red lines show your injected Black Swans; only the Blue Universe expects them.",
         x = "Growth Units", y = "Probability Density")
  print(Universe)
  ggsave("Figures/OutlierUniverse.png", plot = Universe, width = 12, height = 7)
  
###### Stress Test 
  
  #Final  Black-Swan check on Robust model
  
  # --- 1. Robust Simulation Function ---
  run_robust_sim <- function(n) {
    success_count <- replicate(500, {
      # Generate data from the 'Robust Universe' (Student-t, nu=4)
      # This includes the 6% Black Swan risk automatically
      control <- rt.scaled(n, df = 4, mean = 5.0, sd = 1.0)
      treat_2 <- rt.scaled(n, df = 4, mean = 5.5, sd = 1.0) # Signal = 0.5
      
      # Use a Robust Test (or simple check: does Trt2 > Control?)
      # A t.test is 'okay' here but rt.scaled makes it a real stress test
      t.test(control, treat_2, var.equal = FALSE)$p.value < 0.05
    })
    return(mean(success_count))
  }
  
  # --- 2. Execute across N-range ---
  n_range <- seq(10, 200, by = 5)
  robust_power <- sapply(n_range, run_robust_sim)
  sim_df_robust <- data.frame(N = n_range, Power = robust_power)
  
  # --- 3. The New 'Swan-Shield' Plot ---
  
  # Efficiency Index ( Elbow of the curve)  = del Power / del N 
  
  BS_StressTest <- ggplot(sim_df_robust, aes(x = N, y = Power)) +
    # Core Data
    geom_line(color = "#56B4E9", linewidth = 1.5) +
    geom_point(size = 3) +
    
    # Reliability Thresholds
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "darkblue", alpha = 0.6) + 
    
    # Efficiency Peak Intercepts (N=125, Power ~80%)
    geom_vline(xintercept = 125, linetype = "dotdash", color = "darkred", linewidth = 1) +
    geom_hline(yintercept = 0.80, linetype = "dotted", color = "darkred", linewidth = 0.8) +
    
    # Annotations
    annotate("text", x = 35, y = 0.98, label = "95% Reliability Goal", color = "darkblue", fontface = "bold") +
    annotate("label", x = 125, y = 0.80, label = "Efficiency Peak: N=125 (80% Power)", 
             color = "white", fill = "darkred", fontface = "bold", size = 3.5) +
    
    # Aesthetics & Scales
    theme_minimal() +
    scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Robust Efficiency Peak: Diminishing Returns Audit",
      subtitle = "Architecture: Student-t (v=4) | Signal-to-Noise Ratio: 0.5 / 1.0",
      x = "Sample Size (N per group)", 
      y = "Detection Probability (Statistical Power)",
      caption = "Note: Beyond N=125, marginal power gains require exponential sample increases."
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.title = element_text(face = "bold")
    )
  
  print(BS_StressTest)
  ggsave(paste0("Figures/BSStressTest.png"), plot = BS_StressTest, width = 8, height = 5)

  
  # 1. Construct the Milestone Data
  milestone_df <- data.frame(
    Milestone = c("Current Progress", "Efficiency Peak", "Reliability Mandate"),
    N = c(42, 125, 200),
    Power = c("~43%", "~80%", "~95%"),
    Status = c("Gamble. 1-in-2 chance of failure.", 
               "Target. Optimal signal/cost ratio.", 
               "Over-Engineered. High cost for low gain.")
  )
  
  # 2. Render Kable Table
  miles <- kable(milestone_df, 
        format = "pipe", 
        col.names = c("Milestone", "Sample Size (N)", "Reliability (Power)", "Status"),
        align = "llcl")
  write.csv(milestone_df, "data/StressTestMiles.csv")
  # |Milestone           |Sample Size (N) | Reliability (Power) |Status                                   |
  # |:-------------------|:---------------|:-------------------:|:----------------------------------------|
  # |Current Progress    |42              |        ~43%         |Gamble. 1-in-2 chance of failure.        |
  # |Efficiency Peak     |125             |        ~80%         |Target. Optimal signal/cost ratio.       |
  # |Reliability Mandate |200             |        ~95%         |Over-Engineered. High cost for low gain. |

 #### Computational Heirarchy
  # 1. Define Hierarchy Data
  hierarchy_df <- data.frame(
    Model = c("M1: Pooled", "M2: Heteroscedastic", "M3: Robust (Student-t)"),
    Parameters = c("1 mu, 1 sigma", "3 mu, 3 sigma", "3 mu, 3 sigma, fixed nu=4"),
    Compute_Cost = c("Low (O(n))", "Moderate (O(k*n))", "High (O(k*mcmc))"),
    Failure_Mode = c("Bias (Smears noise)", "Fragility (Outlier panic)", "High-N Demand (Skeptical)")
  )
  
  # 2. Render Hierarchy Table
  kable(hierarchy_df, 
        format = "pipe", 
        col.names = c("Model Architecture", "State Variables", "Compute Burden", "Primary Risk"))
  
  write.csv(hierarchy_df, "data/hierarchy.csv")
  
    #|Model Architecture     |State Variables           |Compute Burden    |Primary Risk              |
    #|:----------------------|:-------------------------|:-----------------|:-------------------------|
    #|M1: Pooled             |1 mu, 1 sigma             |Low (O(n))        |Bias (Smears noise)       |
    #|M2: Heteroscedastic    |3 mu, 3 sigma             |Moderate (O(k*n)) |Fragility (Outlier panic) |
    #|M3: Robust (Student-t) |3 mu, 3 sigma, fixed nu=4 |High (O(k*mcmc))  |High-N Demand (Skeptical) |

# Conclusion 
#1.The M3 "Robustness Tax": While M3 is more computationally expensive due to MCMC integration over heavy tails, 
# the cost is purely digital. The physical cost of a false negative at $N=42$ is significantly higher than the 
# CPU time for $N=200$.
  
#2.Signal Recovery: The Robust model successfully identified the $0.90$ noise floor for Treatment 2 by mathematically 
# "downgrading" the specific Black Swan rows (e.g., 93, 108, 113) that corrupted the Gaussian Model 2.
  
#3.Decision: Deploy Phase III with $N=125$. This provides an 80% reliability shield while avoiding the 200-sample 
# efficiency drop-off.
  
  

  