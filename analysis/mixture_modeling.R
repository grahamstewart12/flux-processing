

# Resources
# https://cchecastaldo.github.io/BayesianShortCourse/Syllabus.html
# 

# Mixture Model

# fch4 = p_wtr * fch4_wtr + p_veg * fch4_veg + p_frs * fch4_frs
# F = phi_a * f_a + phi_b * f_b + phi_c * f_c
# - have phi, need f
# - assume f is normally distributed, but cannot assume phi is
# - cannot solve algabraically, but can solve probabilistically 


# Parameter estimation

library(rjags)
#library(coda)

site <- filter(
  site, grow_seas, !is.na(wtd_f), !is.na(fch4), night == 0, ustar > 0.1
)

# ALL DATA - patches

data <- list(
  n = nrow(site), 
  y = pull(site, fch4), 
  x1 = pull(site, wtr_ffp), 
  x2 = pull(site, wtr_ffp), 
  x3 = pull(site, frs_ffp)
)

# Are results dependent on the initial values?
# - assess by requesting independent chains of MCMC sampling
# - is it okay to set the seeds??? unsure
inits <- list(
  list(
    tau = 1, beta1 = 0, beta2 = 0, beta3 = 0, 
    .RNG.name = "base::Mersenne-Twister", .RNG.seed = 3893
  ),
  list(
    tau = 1, beta1 = 0, beta2 = 0, beta3 = 0, 
    .RNG.name = "base::Mersenne-Twister", .RNG.seed = 5629
  ),
  list(
    tau = 1, beta1 = 0, beta2 = 0, beta3 = 0, 
    .RNG.name = "base::Mersenne-Twister", .RNG.seed = 9993
  )
)

Tj_Ar_Data %>%
  rename(Sp = `Sp.`) %>%
  filter(sex != "na") %>%
  group_by(sex) %>%
  do(tidy(t.test(total_days_alive ~ Sp, data = .)))

# Early sampled params may still be distant from posterior distribution maxima
# - delete these, only take second half into account for evaluation of results
n.adapt <- 1000
# How many samples should be taken from the posterior distribution overall
n.update <- 10000
n.iter <- 10000


# Should there be an error term in the model?
sink("JAGS.R")
cat("
model {

  # priors
  
  tau ~ dgamma(0.01, 0.01)
  beta1 ~ dnorm(0, 0.001)
  beta2 ~ dnorm(0, 0.001)
  beta3 ~ dnorm(0, 0.001)
  
  s2 <- 1 / tau
  s <- sqrt(s2)
  
  
  # likelihood
  
  for (i in 1:n) {
    # mu = the 'true' flux
    mu[i] <- beta1 * x1[i] + beta2 * x2[i] + beta3 * x3[i]
    y[i] ~ dnorm(mu[i], tau)
    
  }

  # sample variance
  sv2 <- pow(sd(y[]), 2)

  # bayes R^2
  r2b <- 1 - s2 / sv2

  # expected y
  expected_y <- beta1 * mean(x1[]) + beta2 * mean(x2[]) + beta3 * mean(x3[])
  
}
", fill = TRUE)
sink()

model1_filename <- file.path(
  wd, "analysis", "mixture_modeling", "output",
  paste0("patches_alldata_", tag_out, ".txt")
)

# Write the model to a .txt file for JAGS to access
writeLines(model1, model1_filename)

fit <- jags.model(
  model1_filename, data = data, inits = inits, n.chains = length(inits), 
  n.adapt = n.adapt
)

update(fit, n.iter = n.update)

samples <- coda.samples(
  fit, c("beta1", "beta2", "beta3", "r2b", "s2"), n.iter = n.iter
)

MCMCsummary(samples, round = 4, n.eff = TRUE)

summary(window(samples, start = burnin))

# Plots

# Parameters
MCMCplot(zc.pooled, params = c("beta1", "beta2", "beta3"))

# Posterior density
MCMCplot(fit, params = "mu")

# Assess convergence
MCMCtrace(fit, pdf = FALSE)
gelman.diag(fit)
heidel.diag(fit)
raftery.diag(fit)

BCI <- MCMCpstr(
  fit, params = "mu", func = function(x) quantile(x, c(0.025, 0.5, 0.975))
)
HPDI <-  MCMCpstr(fit, params = "mu", func = function(x) hdi(x, 0.95))


# Mixture modeling
