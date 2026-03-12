# ANOVA 

data("PlantGrowth")
?PlantGrowth

head(PlantGrowth)
boxplot(weight ~ group, data = PlantGrowth)

lmod1 = lm(weight ~ group, data = PlantGrowth)
summary(lmod1)
anova(lmod1)

# model fitting in jags 

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

##### model 2 with indexed Precision 
# Problem statement : Re-fit the JAGS model on plant growth from the lesson with a separate variance for 
# each of the three groups. To do so, modify the model code to index the precision in the normal 
# likelihood by group, just as we did with the mean. Use the same priors as the original model 
# (except in this case it will be three independent priors for the variances).

#Compare the estimates between the original lesson model and this model with the summary function. 
#Notice that the posterior means for the three mu parameters are essentially unchanged. However, the 
# posterior variability for these parameters has changed. The posterior for which group's mean was most 
# affected by fitting separate group variances


set.seed(82)
str(PlantGrowth)
library("rjags")

data_jags = list(y = PlantGrowth$weight, 
                  grp = as.numeric(PlantGrowth$group))

params = c("mu","sig")

# this inits for mod_string_1

inits1 = function(){
  list("mu"=rnorm(3, 0.0, 100), "prec"=rgamma(1, 1.0, 1.0))
}


# change model here 
mod1 = jags.model(textConnection(mod_string_1), data = data_jags,
                 inits = inits1, n.chains = 3)

update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                       variable.names = params,
                       n.iter = 5e3)

mod1_csim =as.mcmc(do.call(rbind, mod1_sim))

plot(mod1_sim)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
effectiveSize(mod1_sim)


pm_params1 = colMeans(mod1_csim)
print(pm_params1) 
coefficients(lmod1)

yhat1 = pm_params1[1:3][data_jags$grp]
yhat1 #predicted values 

resids1 = data_jags$y - yhat1 #actual - predicted
print(resids1)

plot(resids1)
plot(yhat1, resids1)

summary(mod1_sim)
HPDinterval(mod1_csim) #95% interval 
HPDinterval(mod1_csim, 0.9) # 90% posterior interval 

# HPD intervals mu3 - mu 1 

# 1. Create a vector of the differences from your MCMC samples
# This subtracts mu[1] from mu[3] for every single iteration
diff_mu3_mu1 <- mod1_csim[, "mu[3]"] - mod1_csim[, "mu[1]"]
print(diff_mu3_mu1)
# 2. Calculate the HPD interval for this new vector


 # (-0.2, 1.17)
# posterior prob that mu3 > mu 1 0.9378
diff_samples <- mod1_csim[, "mu[3]"] - mod1_csim[, "mu[1]"]
HPDinterval(as.mcmc(diff_samples), prob = 0.93)

# t’s more accurate to say it might reduce yield by 0.20 units.
# The Upper Bound (1.17): There is a chance Treatment 2 results in 1.17 units more growth than the 
# Control.

# The 95% HPD interval for the difference between Treatment 2 and Control is $[-0.20, 1.17]$. 
# Since this interval includes zero, the effect of Treatment 2 is not considered significant at the 95% level, despite a 93.8% posterior probability that the treatment mean is higher than the control."


# 3. Bonus: Calculate the probability that the difference is greater than zero
prob_positive <- mean(diff_mu3_mu1 > 0)
cat("Probability that mu[3] > mu[1]:", prob_positive)

head(mod1_csim)

# check draws from mu[3] >  dmu[1] then take mean to get posterior probabilities 

mean(mod1_csim[,3] > mod1_csim[,1]) # 0.9378
mean(mod1_csim[,3] > 1.1*mod1_csim[,1]) # 0.4874

# 50 : 50 odds that adopting treatment 2 will increase mean yield of a plant by  at least 10 % (Hypostheses)


# statestical modelling technique : best not to yield results 

# Iterations = 1001:6000
#Thinning interval = 1 
#Number of chains = 3 
#Sample size per chain = 5000 

#1. Empirical mean and standard deviation for each variable,
#plus standard error of the mean:
  
#  Mean      SD  Naive SE Time-series SE
#mu[1] 5.0337 0.22663 0.0018505      0.0018367
#mu[2] 4.6620 0.22733 0.0018561      0.0018277
#mu[3] 5.5259 0.22784 0.0018603      0.0018604
#sig   0.7125 0.09227 0.0007533      0.0008198

#2. Quantiles for each variable:
  
#  2.5%    25%    50%    75%  97.5%
#mu[1] 4.5855 4.8848 5.0358 5.1828 5.4810
#mu[2] 4.2056 4.5139 4.6635 4.8131 5.1029
#mu[3] 5.0749 5.3772 5.5277 5.6751 5.9783
#sig   0.5575 0.6477 0.7033 0.7668 0.9219


############################ Model 2

lmod2 = lm(weight ~ group, data = PlantGrowth)
summary(lmod2)
anova(lmod2)

mod_string_2 = "model{
for( i in 1: length(y)){
  y[i] ~ dnorm(mu[grp[i]], prec[grp[i]]) # group level mean of each groups
}

for(j in 1:3){
  
  mu[j] ~ dnorm(0.0, 1.0/1.0e6)
  prec[j] ~ dgamma(5/2.0, 5*1.0/2.0)
  sig[j] = sqrt(1.0 / prec[j])
  
}
}"

inits2 = function(){ 
  list("mu" = rnorm(3, 0.0, 100), "prec" = rgamma(3, 1.0, 1.0))
} 
# change model here 
mod2 = jags.model(textConnection(mod_string_2), data = data_jags,
                  inits = inits2, n.chains = 3)

update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                       variable.names = params,
                       n.iter = 5e3)
mod2_csim =as.mcmc(do.call(rbind, mod2_sim))

plot(mod2_sim)

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
effectiveSize(mod2_sim)

pm_params2 = colMeans(mod2_csim)
print(pm_params2)
coefficients(lmod2)

yhat2 = pm_params2[1:3][data_jags$grp]
yhat2 #predicted values 

resids2 = data_jags$y - yhat2 #actual - predicted
print(resids2)

plot(resids2)
plot(yhat2, resids2)

summary(mod2_sim)
HPDinterval(mod2_csim) #95% interval 
HPDinterval(mod2_csim, 0.9) # 90% posterior interval 

head(mod2_csim)

# check draws from mu[3] >  dmu[1] then take mean to get posterior probabilities 

mean(mod2_csim[,3] > mod2_csim[,1]) # 0.923

# Prob that Treatment 2 is at least 10% better than Control
mean(mod2_csim[,"mu[3]"] > 1.1 * mod2_csim[,"mu[1]"]) #0.488

# model 2 
#Iterations = 1001:6000
#Thinning interval = 1 
#Number of chains = 3 
#Sample size per chain = 5000 

#1. Empirical mean and standard deviation for each variable,
#plus standard error of the mean:
  
#  Mean     SD Naive SE Time-series SE
#mu[1]  5.0302 0.2582 0.002108       0.002140
#mu[2]  4.6625 0.3017 0.002463       0.002534
#mu[3]  5.5238 0.2372 0.001937       0.001953
#sig[1] 0.8030 0.1650 0.001347       0.001450
#sig[2] 0.9253 0.1914 0.001563       0.001697
#sig[3] 0.7357 0.1495 0.001221       0.001276

#2. Quantiles for each variable:
  
#  2.5%    25%    50%   75% 97.5%
#mu[1]  4.5191 4.8663 5.0297 5.194 5.548
#mu[2]  4.0580 4.4705 4.6604 4.855 5.261
#mu[3]  5.0467 5.3683 5.5251 5.677 5.991
#sig[1] 0.5561 0.6875 0.7751 0.891 1.199
#sig[2] 0.6392 0.7911 0.8944 1.027 1.382
#sig[3] 0.5108 0.6308 0.7142 0.817 1.088


# result Model 2 : Treatment 1 

### DIC Criterion 

DIC_1 = dic.samples(mod1, n.iter = 5e3)
print(DIC_1) 

DIC_2 = dic.samples(mod2, n.iter = 5e3)
print(DIC_2)


DIC_1 - DIC_2  #-3.91 model 1 wins 

# Model 1 Penalty: $~4.1$ (This makes sense: 3 means + 1 shared precision = ~4 parameters).
# Model 2 Penalty: $~5.8$ (This also makes sense: 3 means + 3 separate precisions = ~6 parameters).

# In Bayesian model selection, a difference of 3 to 7 points is considered "substantial" evidence. 
# Since Model 1 has the lower score, the data suggests that the extra complexity of indexing the precision 
# doesn't provide enough of a fit improvement to justify the extra parameters.



########## linear model with baseline mean and group effect Cell means model 

# 1. Frequentist Confidence Intervals (from lm)
lmod_cm <- lm(weight ~ -1 + group, data = PlantGrowth) # -1 removes the intercept 
freq_ci <- confint(lmod_cm)

# 2. Bayesian HPD Intervals (from JAGS Model 2)
# Selecting the mu columns from our csim object
bayesian_hpd <- HPDinterval(as.mcmc(mod2_csim[, 1:3]))

# 3. Merging for Comparison
comp_table <- data.frame(
  Group = c("Control", "Trt1", "Trt2"),
  Freq_Lower = freq_ci[, 1],
  Bayes_HPD_Lower = bayesian_hpd[, 1],
  Freq_Upper = freq_ci[, 2],
  Bayes_HPD_Upper = bayesian_hpd[, 2]
)

knitr::kable(comp_table, digits = 3, caption = "Comparison of 95% Intervals")

# Feature	  Standard Model (~ group)	      Cell Means Model (~ -1 + group)
#Intercept	Mean of the first group (Ctrl)	None
#Coefficients	Differences between groups	The actual group means
#P-value meaning	Is Trt different from Ctrl?	Is the group mean different from 0?
#JAGS Equivalent	Baseline/Offset model	Your mu[j] indexed model

# plots 
library(ggplot2)

# 1. Gather all means and intervals
# Frequentist
freq_means <- coef(lmod_cm)
freq_ci    <- confint(lmod_cm)

# Bayesian Model 1 (Pooled)
b1_means <- colMeans(mod1_csim[, 1:3])
b1_hpd   <- HPDinterval(as.mcmc(mod1_csim[, 1:3]))

# Bayesian Model 2 (Separate)
b2_means <- colMeans(mod2_csim[, 1:3])
b2_hpd   <- HPDinterval(as.mcmc(mod2_csim[, 1:3]))

# 2. Combine into one data frame
plot_df_3way <- data.frame(
  Group  = rep(c("Control", "Trt1", "Trt2"), 3),
  Method = rep(c("Frequentist (lm)", "Bayesian (Mod 1)", "Bayesian (Mod 2)"), each = 3),
  Mean   = c(freq_means, b1_means, b2_means),
  Lower  = c(freq_ci[,1], b1_hpd[,1], b2_hpd[,1]),
  Upper  = c(freq_ci[,2], b1_hpd[,2], b2_hpd[,2])
)

# 3. Plot it
ggplot(plot_df_3way, aes(x = Group, y = Mean, color = Method)) +
  geom_point(position = position_dodge(width = 0.6), size = 3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                position = position_dodge(width = 0.6), width = 0.2) +
  theme_minimal() +
  labs(title = "Sensitivity Analysis: Weight Estimates Across Models",
       subtitle = "Comparing Pooled vs. Separate Variance Assumptions")

# power Analysis for model 1 (Winner)
# Parameters from your JAGS output
delta <- 0.49 # The treatment effect (mu3 - mu1)
sigma <- 0.63 # The pooled standard deviation from Model 1

# Power Analysis for 80% Power at alpha = 0.05
# This calculates how many plants per group you need to 'catch' the 0.5 difference
power_res <- power.t.test(delta = delta, sd = sigma, sig.level = 0.05, power = 0.8)

print(power_res$n) 
# Result: ~27 plants per group'

# cross check 
SE = 0.63 /sqrt(27)
print(SE) # 0.1212

new_upper = delta +2*SE
print(new_upper) # 0.7324
new_lower = delta - 2*SE
print(new_lower) # 0.2475

new_intervals = c(new_lower, new_upper)
print(new_intervals)

# The current $N=10$ pilot study indicates a 93.8% probability of treatment success. However, due to 
# the observed variance, the model is underpowered to achieve 95% HPD separation. A subsequent study 
#with $N=27$ per group is recommended to formally validate the treatment effect size of $0.5$ units."


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
# Add 'y_pred' to your monitor list in coda.samples if you add it to the model string
# Or do it manually in R:
n_sim = nrow(mod1_csim)
y_pred = matrix(NA, nrow = n_sim, ncol = length(data_jags$y))

for (i in 1:n_sim) {
  # Extract current iteration's parameters
  current_mu = mod1_csim[i, 1:3]
  current_sig = mod1_csim[i, "sig"]
  
  # Simulate new plant weights
  y_pred[i,] = rnorm(length(data_jags$y), 
                     mean = current_mu[data_jags$grp], 
                     sd = current_sig)
}

# 2. Compare Observed vs Predicted Mean
obs_mean = mean(data_jags$y)
print(obs_mean)5.073
pred_means = rowMeans(y_pred)
print(mean(pred_means))#5.0704
hist(pred_means, main="Posterior Predictive Distribution of Mean", xlab="Weight")
abline(v=obs_mean, col="red", lwd=2) # Red line should fall inside the histogram

# Calculate the Bayesian p-value for the mean
# This is the proportion of simulated means that are greater than the observed mean
p_val_mean <- mean(pred_means > obs_mean)
cat("Bayesian P-value for the Mean:", p_val_mean) #0.4978

# Interpretation: A value near 0.5 indicates a perfect fit. 
# Values < 0.05 or > 0.95 would indicate the model is failing to capture the data.

# Results 
#1.Treatment 2 is the best performer ($93.8\%$ probability).
#2. Treatment 1 is the messiest performer (highest variance).
#3. Model 1 is the most efficient performer (best DIC score)


# conclusions 
#1.Treatment Efficacy: Treatment 2 is likely beneficial, but the small sample size ($n=10$) prevents 
# definitive 95% confidence.
#2. Noise Profile: Treatment 1 increases group variance, making mean estimates for that group less 
# reliable than ctrl or trt2.
#3. Methodology: Bayesian MCMC successfully isolated the "treatment effect" from "biological noise" 
# through the use of indexed precisions.

# Synthetic study  Monte Carlo Power simulation (analyzing the past -> designing the future)
# 1. Simulate a 'Future' Dataset with N=27
set.seed(42)
N_new <- 27

# use posterior mean and sd  from model 1 to simulate 27 plants 
y_new <- c(rnorm(N_new, mean=5.03, sd=0.63), # Simulated Control
           rnorm(N_new, mean=4.66, sd=0.63), # Simulated Trt1
           rnorm(N_new, mean=5.53, sd=0.63)) # Simulated Trt2
grp_new <- rep(1:3, each=N_new)

# 2. Run JAGS on this synthetic data
data_new <- list(y=y_new, grp=grp_new, N=length(y_new))
mod_new <- jags.model(textConnection(mod_string_1), data=data_new, n.chains=3, quiet=TRUE)
update(mod_new, 1000)
sim_new <- coda.samples(mod_new, variable.names="mu", n.iter=5000)
csim_new <- as.mcmc(do.call(rbind, sim_new))

# 3. Check the new HPD for mu[3] - mu[1]
diff_new <- csim_new[,"mu[3]"] - csim_new[,"mu[1]"]
HPD_new <- HPDinterval(as.mcmc(diff_new))
print(HPD_new)

# f you only do 100 loops, you are still vulnerable to "sampling luck." 
# This is the core of the Black Swan problem in statistics: a rare, extreme event (like a batch of 27 "freak" 
# plants) can wreck your averages and make you think a treatment works when it doesn't, or vice-versa.
# intervals can be random (Sampling luck and can subject to blackswans )
# law of large number (antidote to Blackswan) can turn this into calculated risk through meta simulation of 
# In a single run, a Black Swan event (a outlier) ruins the result.In 10,000 runs, the "truth" 
#'the $+0.5$ mean forms a stable bell curve.The Significance: We don't look at one result; we look at 
# the distribution of results. If the 95% HPD clears zero in 850 out of 1,000 simulations, we say the 
# "Power" is 85%. We have quantified the risk of the Black Swan.
#'

# Power Simulation 
# Instead of fearing a Black Swan, we simulate them. We calculate the exact % of times the "randomness" 
# will win over the "signal."

# Function to estimate power for a given N
estimate_power <- function(n_per_group, delta = 0.49, sigma = 0.63, iterations = 100) {
  successes <- 0
  for (i in 1:iterations) {
    # Generate synthetic data
    y_sim <- c(rnorm(n_per_group, 5.03, sigma), rnorm(n_per_group, 5.52, sigma))
    grp_sim <- rep(1:2, each = n_per_group)
    
    # Fast frequentist t-test as a proxy for MCMC power
    test <- t.test(y_sim ~ grp_sim, alternative = "less")
    if (test$p.value < 0.05) successes <- successes + 1
  }
  return(successes / iterations)
}

# Run for different sample sizes
n_range <- seq(5, 50, by = 1)
power_vals <- sapply(n_range, function(n) {
  power.t.test(n = n, delta = 0.49, sd = 0.63, sig.level = 0.05)$power
})

# 2. Base Plot
plot(n_range, power_vals, type = "l", lwd = 2, col = "#0072B2",
     xlab = "Sample Size (n per group)", ylab = "Statistical Power",
     main = "Risk Mitigation Roadmap", ylim = c(0, 1.1))

# 3. Label Findings
# Current State
points(10, power_vals[n_range==10], col = "red", pch = 19, cex = 1.5)
text(10, power_vals[n_range==10], labels = "Current (N=10)\n~50% Risk of Failure", pos = 1, cex = 0.8, col = "red")

# Target State
abline(h = 0.8, col = "darkgreen", lty = 2)
abline(v = 27, col = "darkgreen", lty = 2)
points(27, power_vals[n_range==27], col = "darkgreen", pch = 19, cex = 1.5)
text(27, power_vals[n_range==27], labels = "Target (N=27)\n80% Power Threshold", pos = 4, cex = 0.8, col = "darkgreen")

# Diminishing Returns
arrows(35, 0.98, 45, 0.98, length = 0.1, angle = 20)
text(40, 1.02, "Diminishing Returns Zone", cex = 0.8, font = 3)


### Stochastic (stress test  where you can anticipate sampler luck / blackswans) plot 

# 1. Generate Stochastic Power Data
# Using low iterations (n=40) to preserve 'Black Swan' sampling noise
set.seed(42) 
n_seq <- seq(5, 50, by = 2)
sim_results <- data.frame(N = n_seq)

sim_results$Power <- sapply(sim_results$N, function(n) {
  results <- replicate(40, {
    # Simulating based on Model 1 parameters
    ctrl <- rnorm(n, mean = 5.03, sd = 0.63)
    trt2 <- rnorm(n, mean = 5.53, sd = 0.63)
    t.test(trt2, ctrl, alternative = "greater")$p.value < 0.05
  })
  mean(results)
})

# 2. Plot the Stochastic Roadmap
ggplot(sim_results, aes(x = N, y = Power)) +
  # Region 1: Gamble Zone
  annotate("rect", xmin = 5, xmax = 18, ymin = 0, ymax = 1.1, fill = "red", alpha = 0.1) +
  # Region 2: Reliability Pivot
  annotate("rect", xmin = 18, xmax = 32, ymin = 0, ymax = 1.1, fill = "green", alpha = 0.1) +
  # Region 3: Efficiency Zone
  annotate("rect", xmin = 32, xmax = 50, ymin = 0, ymax = 1.1, fill = "grey", alpha = 0.1) +
  
  # The Stochastic Step Line
  geom_step(color = "#0072B2", linewidth = 1.2) +
  geom_point(color = "#0072B2", size = 2) +
  
  # Threshold line
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgreen") +
  
  # Labels and Annotations
  annotate("text", x = 11.5, y = 0.05, label = "GAMBLE ZONE", color = "red", fontface = "bold") +
  annotate("text", x = 25, y = 0.05, label = "RELIABILITY PIVOT", color = "darkgreen", fontface = "bold") +
  annotate("text", x = 41, y = 0.05, label = "EFFICIENCY ZONE", color = "grey40", fontface = "bold") +
  
  # Highlighting the 'Current' vs 'Target'
  annotate("label", x = 10, y = 0.45, label = "Current N=10\n(High Risk)", color = "red") +
  annotate("label", x = 27, y = 0.88, label = "Target N=27\n(Statistical Hedge)", color = "darkgreen") +
  
  # Aesthetics
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1.1)) +
  theme_minimal() +
  labs(title = "Stochastic Power Roadmap: The 'Black Swan' Defense",
       subtitle = "Quantifying risk through three strategic operational zones",
       x = "Sample Size (n per group)",
       y = "Observed Power (1 - Beta)")


### final verdict 
# Treatment 2 should be viewed as a High-Potential Asset. The model confirms a stable mean increase, 
# but current data carries a 6.2% risk of the observed effect being driven by sampling stochasticity. 
# To mitigate this risk and achieve 95% HPD separation (no overlap with zero), the study must be scaled 
# to 27 plants per group.


# Project Title: Stability Analysis of Treatment 2 Growth EfficacyThe 
# Result: Treatment 2 demonstrates a consistent mean growth increase of 0.5 units over the control. 
# Bayesian modeling confirms a 93.8% probability that this effect is positive.
# The Risk: At the current pilot scale ($N=10$), the study is "Underpowered." There is a 6.2% risk that 
# the observed success is a sampling fluke. Our stochastic roadmap identifies this as the "Gamble Zone," 
# where random "Black Swan" outliers can still flip the conclusion.
# The Action: To achieve 95% certainty and reach the "Reliability Pivot," the next phase must scale to 
# 27 units per group. This move hedges against sampling volatility and provides the statistical 
# integrity required for a final production rollout.



