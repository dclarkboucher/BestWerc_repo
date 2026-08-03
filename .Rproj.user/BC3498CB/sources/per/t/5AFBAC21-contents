library(dplyr)
library(tidyr)
library(RBesT)
library(SAMprior)
library(patchwork)
library(ggplot2)
source("source/normal_helpers.R")

# Functions ---------------------------------------------------------------
analyze <- function(y0_dat, y1_dat,
                    prior0, prior1, sigma, method = c("MAP","SAM"), 
                    wt = 0.20, csd = NA){
  
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
    
    
    # Try Bayesian method
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
  }
  
  post_y0 <- suppressMessages(postmix(prior0, data = y0_dat))
  dimnames(post_y0)[[2]] <- paste0("comp",seq_along(dimnames(post_y0)[[2]]))
  
  post_y1 <- suppressMessages(postmix(prior1, data = y1_dat))
  
  threshold <- thresholds[method]
  
  decision <-  pmixdiff(post_y0, post_y1, 0) > threshold
  
  # Have to check weights before ess to prevent errors
  if (any(post_y0[1,] < 1e-4)){
    post_y0[1,] <- (post_y0[1,] + 1e-4) / sum(post_y0[1,]+ 1e-4)
  }
  
  ess <- ess(post_y0, sigma = sigma)
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

# Simulations -------------------------------------------------------------
nreps <- 5000 
methods <- c("SAM","MAP")
rs <- c(-0.5,-0.4,-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)
sam_csd <- 0.1
thresholds <- c(0.95, 0.95)
names(thresholds) <- methods
map_wt <- 0.7 

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
          method = method,
          sigma = sigma,
          wt = map_wt, 
          csd = sam_csd
        )
      
      d0 <- out["decision"]
      ess <- out["ess"]
      mu_0_hat <- out["mean"]
      
      # Performance under alternative
      out <-
        analyze(
          y0_dat = y0_dat, 
          y1_dat = y1_dat + delta, 
          prior0 = map, 
          prior1 = weak_prior,
          method = method,
          sigma = sigma,
          wt = map_wt, 
          csd = sam_csd
        )
      
      d1 <- out["decision"]
      
      # Update
      t1e <- t1e + d0
      power <- power + d1
      biases[rep] <- mu_0_hat - mu_0
      esses[rep] <- ess - n_0
      
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
        mu_0 = mu_0,
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
results1 <- readRDS("results/sims_normal.rds") %>% 
  filter(!(method %in% c("SAM","MAP")))
results <- bind_rows(results, results1)

saveRDS(results, file = "results/sims_normal_recalibrated.rds")
# Plot output -------------------------------------------------------------
colors <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")
colors <- c("#2d715c", "#fda988", "#657eb9",  "#ea9acb", "firebrick1")
smallfont <- 9
largefont <- 10

results %>% 
  select(method, mu_0, bias, mse, power, type_1_error) %>% 
  mutate(
    r = (y_0h - mu_0) /delta
  ) %>% 
  pivot_longer(cols = c(power, type_1_error,bias, mse),
               names_to = "metric",
               values_to = "performance") %>% 
  mutate(
    metric1 = 
      factor(metric,
             levels = c("power", "type_1_error", "bias", "mse"),
             labels = c("Power", "T1E Rate",  "Bias", "MSE"),
             ordered = TRUE
      ),
    method1 = 
      factor(
        method, 
        levels = c("BW_fixed", "BW_random", "MAP", "SAM", "EP"),
        labels = c("BestWerc with Known Variance", "BestWerc with Estimated Variance", "MAP", "SAM", "EP")
      )
  ) %>% 
  ggplot(aes(x = r, y = performance, color = method1, linetype = method1)) + 
  geom_line(linewidth = 0.7) + 
  facet_wrap(~metric1, scales = "free_y"
  ) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = c(1, 2, 4, 5,3)) +
  scale_x_continuous(breaks = round(seq(-0.6, 0.6, by = 0.2),1)) + 
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
  ) + 
  guides(color = guide_legend(nrow = 2), linetype =guide_legend(nrow = 2))
ggsave("figures/sims_normal_recalibrated.tiff", compression = "lzw",
       width = 6.5, height = 4, units = "in", dpi = 800)


# Tables ------------------------------------------------------------------
tab <-
  results %>% 
  filter(round(r,2) %in% c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5)) %>% 
  mutate(
    across(c(mu_0, bias, bias_sd, ess, ess_sd, mse, power, type_1_error),
           ~ as.character(format(round(.x, 3),3))),
    
    bias = ifelse(method %in% c("SAM","MAP","EP"), paste0(bias, " (",bias_sd,")"), bias),
    ess = ifelse(method %in% c("SAM","MAP","EP"), paste0(ess, " (",ess_sd,")"), 13.3),
    method = factor(method, levels = c("BW_random", "BW_fixed", "MAP", "SAM", "EP"),
                    ordered = TRUE
    )
  ) %>% 
  select(method, mu_0, power, type_1_error, bias, mse, ess) %>% 
  arrange(desc(mu_0)) %>% 
  pivot_longer(cols = c(power, type_1_error, bias, mse, ess),
               names_to = "metric",
               values_to = "performance") %>% 
  pivot_wider(names_from = mu_0, values_from = performance) %>% 
  arrange(method)
View(tab)
write.csv(tab, "tables/sims_normal_recalibrated.csv")
