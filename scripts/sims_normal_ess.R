library(dplyr)
library(tidyr)
library(RBesT)
library(SAMprior)
library(patchwork)
library(ggplot2)
source("source/normal_helpers.R")

# Functions ---------------------------------------------------------------
g_fcn <- function(t, a, b){
  1 / (1 + exp(a + b * log(t)))
}
t_fcn <- function(mean_h, sigma_h, n_h, mean_0, sigma_0, n_0){
  
  s <- sqrt((n_h-1) * sigma_h^2 + (n_0-1) * sigma_0^2) / 
    sqrt(n_0 + n_h - 2) * sqrt(1/n_h + 1/n_0)
  abs(mean_0 - mean_h) / s * max(n_h,n_0)^(-0.25)
}


get_a_b <- function(y_0h, n_h, csd, R, n_0, c1 = 0.99, c2 = 0.01, 
                    q0 = 0.5, q1 = 0.5, sigma){
  
  # Congruent cases
  means <- rnorm(R, mean = y_0h, sd = sigma/sqrt(n_0))
  sds <- sqrt(rchisq(R, df = n_0-1) * (sigma^2) / (n_0 - 1))
  t_list <- numeric(R)
  for (i in 1:R){
    t_list[i] <- 
      t_fcn(
        mean_h = y_0h, 
        sigma_h = sigma, 
        n_h = n_h, 
        mean_0 = means[i], 
        sigma_0 = sds[i], 
        n_0 = n_0
      )
  }
  
  # Incongruent cases
  t_list_pos <- numeric(R)
  t_list_neg <- numeric(R)
  for (i in 1:R){
    t_list_pos[i] <- 
      t_fcn(
        mean_h = y_0h, 
        sigma_h = sigma, 
        n_h = n_h, 
        mean_0 = means[i] + csd, 
        sigma_0 = sds[i], 
        n_0 = n_0
      )
    t_list_neg[i] <- 
      t_fcn(
        mean_h = y_0h, 
        sigma_h = sigma, 
        n_h = n_h, 
        mean_0 = means[i] - csd, 
        sigma_0 = sds[i], 
        n_0 = n_0
      )
  }
  
  T_q1 <- min(quantile(t_list_pos, q1), quantile(t_list_neg, q1)) # Test stat for no borrowing
  T_q0 <- quantile(t_list, q0) # Test stat for full borrowing. 
  
  a <- log(1-c1) - log(c1) -(log((1-c1)*c2/c1/(1-c2))*log(T_q0)) /(log(T_q0) - log(T_q1))
  b <- log((1-c1)*c2/c1/(1-c2)) / (log(T_q0)-log(T_q1))
  list(
    a = a, 
    b = b,
    T_q0 = T_q0,
    T_q1 = T_q1
  )
  
}


analyze <- function(y0_dat, y1_dat,
                    prior0, prior1, ep_prior, sigma, method = c("MAP","SAM", "EP"), 
                    wt = 0.20, csd = NA, a_ep, b_ep){
  
  method <- match.arg(method)
  
  if (method == "MAP"){
    
    prior0 <-
      suppressMessages(
        robustify(
          prior0, 
          weight = wt, 
          mean = summary(prior0)["mean"],
          sigma = sigma
        ))
    
  } else if (method == "SAM"){
    
    mean_use <- summary(prior0)["mean"]
    h0_prob <- pmix(prior0, mean_use + csd) - pmix(prior0, mean_use - csd)
    
    sam_weight <-
      suppressMessages(SAM_weight(
        if.prior = prior0,
        # method.w = "PPR",
        # prior.odds = h0_prob / (1-h0_prob),
        delta = csd, # as long as drift is in expected range, ok
        data = y0_dat
      ))
    # sam_weight <- max(min(sam_weight, 0.99),0.01)
    
    prior0 <-
      SAM_prior(
        if.prior = prior0,
        nf.prior = prior1, # weakly informative can go here too, that's fine
        weight = sam_weight,
        sigma = sigma
      )
  }  else if (method == "EP"){
    prior0 <- ep_prior
    
    # Determine borrowing fraction
    t_ep <- t_fcn(
      mean_h = ep_prior["m",],
      sigma_h = sigma, 
      n_h = ep_prior["s",] ^ (-2), 
      mean_0 = mean(y0_dat), 
      sigma_0 = sd(y0_dat), 
      n_0 = length(y0_dat)
    )
    g_ep <- g_fcn(t = t_ep, a = a_ep, b = b_ep)
    if(g_ep < 10^-5) g_ep <- 10^(-5)
    
    prior0["s",1] <- prior0["s",1] / g_ep
    
  }
  
  if (method == "EP"){
    post_y0 <- suppressMessages(postmix(prior0, m = mean(y0_dat), n = length(y0_dat)))
    dimnames(post_y0)[[2]] <- paste0("comp",seq_along(dimnames(post_y0)[[2]]))
    
    post_y1 <- suppressMessages(postmix(prior1, m = mean(y1_dat), n = length(y1_dat)))
    threshold <- thresholds[method]    
    
  } else{
    post_y0 <- suppressMessages(postmix(prior0, data = y0_dat))
    dimnames(post_y0)[[2]] <- paste0("comp",seq_along(dimnames(post_y0)[[2]]))
    
    post_y1 <- suppressMessages(postmix(prior1, data = y1_dat))
    threshold <- thresholds[method]    
    
  }

  
  decision <-  pmixdiff(post_y0, post_y1, 0) > threshold
  
  # Have to check weights before ess to prevent errors
  if (any(post_y0[1,] < 1e-4)){
    post_y0[1,] <- (post_y0[1,] + 1e-4) / sum(post_y0[1,]+ 1e-4)
  }
  
  if (method == "EP"){
    ess <- ess(prior0) + length(y0_dat)
  } else {
    ess <- ess(post_y0, sigma = sigma)
    
  }
  
  # ess <- ess(prior0, sigma = sigma)
  # ess <- 1
  mean <- summary(post_y0)['mean']
  out <- c(decision, ess, mean)
  names(out) <- c("decision", "ess", "mean")
  out
}



# Settings ----------------------------------------------------------------

k <- 2
lambda <- 1
sigma <- 1/sqrt(lambda)
n_1 <- 40
n_0 <- n_1 / k
prob0 <- 0.95
power_want <- 0.65
Delta <- 0
mu_0 <- 0
delta <- (qnorm(prob0) + qnorm(power_want)) * sqrt(k+1) / sqrt(n_1 * lambda)
alpha_0 <- 0.05



# Prior derivation --------------------------------------------------------

sigma <- 1
study_data <-
  data.frame(
    study = c("A","B", "C"),
    n = c(120, 150, 170),
    sd = sigma, 
    estimate  = c(0.18, 0.22,0.20)
  )
study_data$se <- with(study_data, sd / sqrt(n))

.user_mc_options <- 
  options(
    RBesT.MC.warmup = 50 * 200, 
    RBesT.MC.iter = 100 * 200,
    RBesT.MC.chains = 2, 
    RBesT.MC.thin = 5,
    RBesT.MC.control=list(adapt_delta=0.999))

set.seed(1)
map_mcmc <- gMAP(cbind(estimate, se) ~ 1 | study,
                 data = study_data,
                 weights = n,
                 beta.prior = cbind(0, sigma),
                 tau.dist = "HalfNormal",
                 tau.prior = sigma/5,
                 family = gaussian)

map <- automixfit(map_mcmc) 
map_sum <- summary(map)
max(map_sum[5] - map_sum[1], map_sum[1] - map_sum[3]) / delta

summary(map)
y_0h <- summary(map)['mean']; y_0h

weak_prior <- 
  mixnorm(c(1,y_0h,1),
          sigma = sigma, 
          param = "mn"
  )

######## ELASTIC PRIOR
# Informative prior
n_h <- sum(study_data$n)
ep_initial_prior <-
  RBesT::mixnorm(
    comp1 = c(1, y_0h, n_h),
    sigma = sigma, 
    param = "mn"
  )

ep_csd <- abs(delta) # clinically significant difference. 
q0 <- 0.025
q1 <- 0.275
R <- 5000
set.seed(1)
hparams <-
  get_a_b(
    y_0h = y_0h,
    n_h = n_h,
    csd = ep_csd,
    R = R,
    n_0 = n_0,
    c1 = 0.99,
    c2 = 0.01,
    q0 = q0,
    q1 = q1,
    sigma = sigma)
a_ep <- hparams$a
b_ep <- hparams$b

# Simulations -------------------------------------------------------------
nreps <- 1000 
methods <- c("SAM","MAP","EP")
rs <- c(-0.5,-0.3, -0.1, 0, 0.1, 0.3, 0.5)
sam_csd <- abs(delta)
thresholds <- c(0.972, 0.975,.95)
names(thresholds) <- c("SAM","MAP","EP")
thresholds <- thresholds[methods]
map_wt <- 0.20 

out_list <- list()
j <- 1
for (r in rs){
  
  mu_0 <- y_0h - r * delta
  
  for (method in methods){
    
    power <- 0
    t1e <- 0
    biases <- numeric(nreps)
    esses <- numeric(nreps)
    for (rep in seq_len(nreps)){
      
      set.seed(rep)
      y0_dat <- rnorm(n_0, mean = mu_0, sd = sigma)
      y1_dat <- rnorm(n_1, mean = mu_0, sd = sigma)
      
      # Performance under null
      out <-
        analyze(
          y0_dat = y0_dat, 
          y1_dat = y1_dat, 
          prior0 = map, 
          prior1 = weak_prior,
          ep_prior = ep_initial_prior,
          method = method,
          sigma = sigma,
          wt = map_wt, 
          csd = sam_csd,
          a_ep = a_ep,
          b_ep = b_ep
        )
      
      d0 <- out["decision"]
      ess <- out["ess"]
      mu_0_hat <- out["mean"]
      out_list[[j]] <-
        data.frame(
          method = method,
          r = r,
          rep = rep,
          scenario = 0,
          ess = ess - n_0,
          est = mu_0_hat,
          d = d0,
          mn = mean(y0_dat)
        )
      j <- j + 1
      
      # Performance under alternative
      out <-
        analyze(
          y0_dat = y0_dat, 
          y1_dat = y1_dat + delta, 
          prior0 = map, 
          prior1 = weak_prior,
          ep_prior = ep_initial_prior,
          method = method,
          sigma = sigma,
          wt = map_wt, 
          csd = sam_csd,
          a_ep = a_ep,
          b_ep = b_ep
        )
      
      d1 <- out["decision"]
      out_list[[j]] <-
        data.frame(
          method = method,
          r = r,
          rep = rep,
          scenario = 1,
          ess = ess - n_0,
          est = mu_0_hat,
          d = d1,
          mn = mean(y1_dat)
        )
      j <- j + 1
    }
  }
}

results <-
  bind_rows(out_list) %>% 
  tibble()
saveRDS(results,"results/sims_normal_esscomp.rds")

results1 <-
  results %>% 
  mutate(
    mu_0 = y_0h - r * delta,
    bias = est - mu_0
  )
results1 %>% 
  filter(scenario==1) %>% 
  group_by(r) %>% 
  summarize(
    bias = mean(bias),
    bias_sd = sd(est)
  )
results1 %>% 
  ggplot(aes(x = bias)) + 
  facet_wrap(~factor(r)) + 
  geom_histogram(binwidth = 0.05)

results1 %>% 
  ggplot(aes(x = ess)) + 
  facet_wrap(~factor(r)) + 
  geom_histogram()

results1 %>% 
  filter(scenario == 0, method == "EP") %>% 
  ggplot(aes(y = ess, x = mn)) + 
  facet_wrap(~factor(r)) + 
  geom_point(color = "dodgerblue4", alpha=0.2) + 
  labs(x = "Raw sample mean of controls", y = "ESS") + 
  theme_bw()

results1 %>% 
  filter(scenario == 0,r %in% c(-0.5, 0, 0.5)) %>%
  mutate(r = paste0("r = ",r)) %>% 
  ggplot(aes(y = ess, x = mn, color = method)) + 
  facet_wrap(~factor(r)) + 
  geom_vline(aes(xintercept = y_0h), color = "black") + 
  geom_point(alpha=0.3) + 
  scale_color_manual(values = c("dodgerblue4", "firebrick4", "darkorange4")) + 
  labs(x = "Raw sample mean of controls", y = "ESS") +
  theme_bw()
ggsave("figures/sims_normal_ess.tiff", compression = "lzw",
       width = 6.5, height = 5, units = "in", dpi = 800)

# Why does this happen?
# --Tuning the error rate requires a conservative rule.
# --By design, the EP allows the maximum borrowing to be much higher than SAM or MAP.
# --Obviously, maximum borrowing would result in extreme bias
# and inflated error rates when r < 0. 
# --To prevent this from happening, the OCs must be tuned very conservatively. 
# --However, you are still going to get maximum borrowing in rare cases when
# the sample mean closely matches the historical mean. 


