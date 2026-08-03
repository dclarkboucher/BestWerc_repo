library(dplyr)
library(tidyr)
library(RBesT)
library(SAMprior)
library(ggplot2)
library(patchwork)
source("source/binary_helpers.R")

# Functions ---------------------------------------------------------------
g_fcn <- function(t, a, b){
  1 / (1 + exp(a + b * log(t)))
}
t_fcn <- function(n_h1, n_h0, n_c1, n_c0){
  X <- matrix(c(n_h1, n_h0, n_c1, n_c0), ncol = 2)
  chisq.test(X, correct = F)$statistic * max(n_h1 + n_h0, n_c1 + n_c0)^(-0.25)
}
get_a_b <- function(n_h1, n_h0, csd, R, n_0, c1 = 0.99, c2 = 0.01, 
                    q0 = 0.5, q1 = 0.5){
  # Congruent cases
  p_h <- n_h1 / (n_h0 + n_h1)
  n_c1s <- rbinom(R, n_0, prob = p_h)
  n_c0s <- n_0 - n_c1s
  t_list <- numeric(R)
  for (i in 1:R){
    t_list[i] <- t_fcn(n_h1, n_h0, n_c1s[i], n_c0s[i])
  }
  
  # Incongruent cases
  if (p_h+csd > 1| p_h-csd < 0 ) stop("CSD too high")
  n_c1s_pos <- rbinom(R, n_0, prob = p_h + csd)
  n_c0s_pos <- n_0 - n_c1s_pos
  t_list_pos <- numeric(R)
  n_c1s_neg <- rbinom(R, n_0, prob = p_h - csd)
  n_c0s_neg <- n_0 - n_c1s_neg
  t_list_neg <- numeric(R)
  for (i in 1:R){
    t_list_pos[i] <- t_fcn(n_h1, n_h0, n_c1s_pos[i], n_c0s_pos[i])
    t_list_neg[i] <- t_fcn(n_h1, n_h0, n_c1s_neg[i], n_c0s_neg[i])
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

analyze <- function(n_0, n_1, x_0, x_1, prior0, prior1, ep_prior, method = c("EP", "MAP","SAM"), 
                    wt = 0.20, csd = NA, ep_csd = NA, a_ep = NA, b_ep = NA,
                    n_h1 = NA, n_h0 = NA){
  
  
  method <- match.arg(method)
  
  if (method == "MAP"){
    
    prior0 <-
      suppressMessages(
        robustify(
          prior0, 
          weight = wt, # CHANGED FROM 0.20
          mean = summary(prior0)["mean"], 
          n = 2
        ))
    
    
  } else if (method == "SAM"){
    
    # Try Bayesian method
    mean_use <- summary(prior0)["mean"]
    sam_weight <-
      suppressMessages(SAM_weight(
        if.prior = prior0,
        delta = csd, # as long as drift is in expected range, ok
        n = n_0,
        r = x_0
      ))
    prior0 <-
      SAM_prior(
        if.prior = prior0,
        nf.prior = prior1, # weakly informative can go here too, that's fine
        weight = sam_weight,
        n = n_0, 
        r = x_0
      )
  } else if (method == "EP"){
    prior0 <- ep_prior
    
    # Determine borrowing fraction
    t_ep <- t_fcn(n_h1 = n_h1, n_h0 = n_h0, n_c1 = x_0, n_c0 = n_0 - x_0)
    g_ep <- g_fcn(t = t_ep, a = a_ep, b = b_ep)
    if(g_ep < 10^-5) g_ep <- 10^(-5)
    
    prior0["a",1] <- prior0["a",1] * g_ep
    prior0["b",1] <- prior0["b",1] * g_ep
    
  }
  
  threshold <- thresholds[method]
  post_y0 <- suppressMessages(postmix(prior0, n = n_0, r = x_0))
  dimnames(post_y0)[[2]] <- paste0("comp",seq_along(dimnames(post_y0)[[2]]))
  
  post_y1 <- suppressMessages(postmix(prior1, n = n_1, r = x_1))
  
  decision <-  pmixdiff(post_y0, post_y1, 0) > threshold
  
  # Have to check weights before ess to prevent errors
  if (any(post_y0[1,] < 1e-4)){
    post_y0[1,] <- (post_y0[1,] + 1e-4) / sum(post_y0[1,]+ 1e-4)
  }
  
  # ess <- ess(post_y0, sigma = sigma)
  ess <- try(ess(post_y0), silent = T)
  if (inherits(ess,"try-error")) ess <- sum(post_y0[c("a","b"),])
  # ess <- 1
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
R_0 <- 0
R_0 <- 0.40
delta <- find_R1(
  power = 0.65,
  R0 = R_0,
  n_1 = n_1, 
  n_0 = n_1/k
) - R_0
alpha_0 <- 0.05

# Prior derivation --------------------------------------------------------

study_data <-
  data.frame(
    study = c("A","B", "C"),
    n = c(52, 58, 71),
    x1 = c(16, 23, 33)
  )

.user_mc_options <- 
  options(
    RBesT.MC.warmup = 50 * 200, 
    RBesT.MC.iter = 100 * 200,
    RBesT.MC.chains = 2, 
    RBesT.MC.thin = 5,
    RBesT.MC.control=list(adapt_delta=0.999))

set.seed(1)
map_mcmc <- gMAP(cbind(x1, n-x1) ~ 1 | study,
                 data = study_data,
                 beta.prior = cbind(0, 1),
                 tau.dist = "HalfNormal",
                 tau.prior = 0.2,
                 family = binomial)

map <- automixfit(map_mcmc) 
map_sum <- summary(map)
summary(map)['mean']
summary(map)
max(map_sum[5] - map_sum[1], map_sum[1] - map_sum[3]) /delta

summary(map)
R_0h <- summary(map)['mean']; R_0h

weak_prior <- 
  mixbeta(
    nf.prior = c(1,1,1)
  )

######## ELASTIC PRIOR
# Informative prior
alpha1 <- sum(study_data$x1)
beta1 <- sum(study_data$n) - alpha1

# Weakly informative prior
alpha0 <- 0.1
beta0 <- 0.1

# Initial unweighted prior
alpha01 <- alpha0 + alpha1
beta01 <- beta0 + beta1
ep_initial_prior <-
  RBesT::mixbeta(
    comp1 = c(1, alpha01, beta01),
    param = c("ab")
  )
ep_csd <- abs(delta) # clinically significant difference. 
R <- 5000
q0 <- 0.025
q1 <- 0.025
set.seed(1)
hparams <- 
  get_a_b(
    n_h1 = alpha1, 
    n_h0 = beta1, 
    csd = ep_csd,
    R = R, 
    n_0 = n_0, 
    c1 = 0.99, 
    c2 = 0.01,
    q0 = q0, 
    q1 = q1)
a_ep <- hparams$a
b_ep <- hparams$b


# Simulations -------------------------------------------------------------
nreps <- 5000
methods <- c("MAP","SAM") 
rs <- c(-0.5,-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)

sam_csd <- 0.05 
map_wt <- 0.75
thresholds <- c(0.95, 0.95)
names(thresholds) <- methods

out_list <- list()
j <- 1
for (r in rs){
  
  R_0 <- R_0h - r * delta
  
  for (method in methods){
    
    power <- 0
    t1e <- 0
    biases <- numeric(nreps)
    esses <- numeric(nreps)
    for (rep in seq_len(nreps)){
      
      # Performance under null
      set.seed(rep)
      x_0 <- rbinom(1, n_0, R_0)
      x_1 <- rbinom(1, n_1, R_0)
      out <-
        analyze(
          n_0 = n_0,
          n_1 = n_1, 
          x_0 = x_0,
          x_1 = x_1,
          prior0 = map, 
          prior1 = weak_prior,
          ep_prior = ep_initial_prior,
          method = method,
          wt = map_wt, 
          csd = sam_csd,
          ep_csd = ep_csd, 
          a_ep = a_ep, 
          b_ep = b_ep,
          n_h1 = alpha1, 
          n_h0 = beta1
        )
      d0 <- out["decision"]
      ess <- out["ess"]
      R_0_hat <- out["mean"]
      
      # Performance under alternative
      set.seed(rep)
      x_0 <- rbinom(1, n_0, R_0)
      x_1 <- rbinom(1, n_1, R_0 + delta)
      out <-
        analyze(
          n_0 = n_0,
          n_1 = n_1, 
          x_0 = x_0,
          x_1 = x_1,
          prior0 = map, 
          prior1 = weak_prior,
          ep_prior = ep_initial_prior,
          method = method,
          wt = map_wt, 
          csd = sam_csd,
          ep_csd = ep_csd, 
          a_ep = a_ep, 
          b_ep = b_ep,
          n_h1 = alpha1, 
          n_h0 = beta1
        )
      
      d1 <- out["decision"]
      
      # Update
      t1e <- t1e + d0
      power <- power + d1
      biases[rep] <- R_0_hat - R_0
      esses[rep] <- ess
      
    }
    
    
    bias <- mean(biases)
    ess <- mean(esses)
    ess_sd <- sd(esses)
    mse <- mean(biases^2)
    sd <- sd(biases)
    
    out_list[[j]] <- 
      data.frame(
        n_0 = n_0,
        n_1 = n_1,
        method = method,
        R_0 = R_0,
        r = r, 
        bias = bias, 
        bias_sd = sd,
        ess = ess, 
        ess_sd = ess_sd,
        mse = mse,
        power = power / nreps,
        type_1_error = t1e / nreps
      )
    j <- j + 1
    
    
  }
  
}
results <- 
  bind_rows(out_list) %>% 
  as_tibble()

# Load results for other methods 
results1 <- readRDS("results/sims_binary.rds") %>% 
  filter(!(method %in% c("SAM","MAP")))
results1 <- bind_rows(results, results1)
saveRDS(results1, file = "results/sims_binary_recalibrated.rds")

# Plot output -------------------------------------------------------------
# colors <- c("#fc8d62",  "#8da0cb", "#e78ac3")
colors <- c("#fda988", "#657eb9",  "#ea9acb", "firebrick1")
smallfont <- 9
largefont <- 10

results1 %>% 
  select(method, R_0, bias, mse, power, type_1_error) %>% 
  mutate(r = (R_0h - R_0) /delta
  ) %>% 
  pivot_longer(cols = c(power, type_1_error,bias, mse),
               names_to = "metric",
               values_to = "performance") %>% 
  mutate(
    metric1 = 
      factor(metric,
             levels = c("power", "type_1_error", "bias", "mse"),
             labels = c("Power", "T1E Rate", "Bias", "MSE"),
             ordered = TRUE
      ),
    method1 = 
      factor(
        method, 
        levels = c( "BW", "MAP", "SAM","EP"),
        labels = c("BestWerc", "MAP", "SAM", "EP")
      )
  ) %>% 
  ggplot(aes(x = r, y = performance, color = method1, linetype = method1)) + 
  geom_line(linewidth = 0.7) + 
  facet_wrap(~metric1, scales = "free_y") +
  scale_color_manual(values = colors) + 
  # geom_hline(yintercept = 0.1) +
  scale_linetype_manual(values = c(2, 4, 5, 3)) + 
  scale_x_continuous(breaks = round(seq(-0.6, 0.6, by = 0.2),1)) + 
  # scale_y_continuous(breaks = seq(0, 1, 0.01)) +
  labs(x = "Relative Drift (r)") +
  theme_minimal() + 
  theme(
    panel.grid = element_blank(),
    legend.key.width = unit(0.5,"in"),
    panel.background = element_rect(color = "black", linewidth = 0.5),
    legend.title = element_blank(),
    legend.position = "top",
    # legend.background = element_rect(color = "black", linewidth = 0.3),
    legend.background = element_blank(),
    legend.key.spacing.y = unit(0.01,"in"),
    legend.spacing.y = unit(0.01, "in"),
    # legend.box.margin = unit(rep(0,4),"in"),
    legend.box.spacing = unit(0,"in"),
    # legend.justification = "left",
    axis.title.y = element_blank(),
    legend.margin = margin(0.04, 0.04, 0.04, 0.04, unit="in"),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.title.x = element_text(size = largefont),
    legend.text = element_text(size = largefont),
    axis.text = element_text(size = smallfont),
    strip.text = element_text(size = largefont)
  )
ggsave("figures/sims_binary_recalibrated.tiff", compression = "lzw",
       width = 6.5, height = 3.8, units = "in", dpi = 800)


# Tables ------------------------------------------------------------------
tab <-
  results1 %>% 
  filter(round(r,2) %in% c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5)) %>% 
  mutate(
    across(c(R_0, bias, bias_sd, ess, ess_sd, mse, power, type_1_error),
           ~ as.character(format(round(.x, 3),3))),
    
    bias = ifelse(method %in% c("SAM","MAP", "EM"), paste0(bias, " (",bias_sd,")"), bias),
    ess = ifelse(method %in% c("SAM","MAP", "EM"), paste0(ess, " (",ess_sd,")"), 13.3),
    method = factor(method, levels = c("BW", "MAP", "SAM", "EP"),
                    ordered = TRUE
    )
  ) %>% 
  select(method, R_0, power, type_1_error, bias, mse, ess) %>% 
  arrange(desc(R_0)) %>% 
  pivot_longer(cols = c(power, type_1_error, bias, mse, ess),
               names_to = "metric",
               values_to = "performance") %>% 
  pivot_wider(names_from = R_0, values_from = performance) %>% 
  arrange(method)
View(tab)
write.csv(tab, "tables/sims_binary_recalibrated.csv")





