library(dplyr)
library(tidyr)
library(RBesT)
library(SAMprior)
library(patchwork)
source("source/normal_helpers.R")



# Set up and priors -------------------------------------------------------

k <- 2
lambda <- 1
sigma <- 1/sqrt(lambda)
prob0 <- 0.95
power_want <- 0.65
Delta <- 0
mu_0 <- 0
alpha_0 <- 0.05
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
y_0h <- summary(map)['mean']; y_0h
rs <- c(-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5)


# Large sample (n = 40) ---------------------------------------------------
n_1 <- 40
n_0 <- n_1 / k
delta <- (qnorm(prob0) + qnorm(power_want)) * sqrt(k+1) / sqrt(n_1 * lambda)

# Random variance case
mu_0_use <-  y_0h - rs * delta

bw_random <-
  fast_sim_numerical_wrapper_flex(
    mu_0 = mu_0_use,
    r = rs, 
    alpha = 0.1, 
    beta = 0.1, 
    n_1 = n_1, 
    k = k, 
    a_0 = 0.4, 
    delta = delta,
    alpha_0 = 0.05,
    gridsize = 200,
    tail_error = 1e-5
  )$results %>% 
  mutate(
    method = "BW_random",
    n_0 = n_1 / k
  ) %>% 
  select(-delta, -alpha, -beta, -a_0,-k)

# Fixed variance case
bw_fixed <-
  data.frame(
    n_1 = n_1, 
    n_0 = n_1 /k,
    mu_0 = mu_0_use,
    r = rs,
    a_0 = 0.4
  ) %>% 
  mutate(
    method = "BW_fixed",
    n_0 = floor(n_1 / k),
    power = 
      pnorm(
        (delta - a_0 * r * delta - qnorm(prob0) * sigma * sqrt(1/n_1 + (1-a_0)/(n_0))) /
          (sigma * sqrt(1/n_1 + (1-a_0)^2 /n_0 ))
      ),
    type_1_error = 
      pnorm(
        (- a_0 * r * delta - qnorm(prob0) * sigma * sqrt(1/n_1 + (1-a_0)/(n_0))) /
          (sigma * sqrt(1/n_1 + (1-a_0)^2 /n_0 ))
      ),
    n_0 = n_1 / k
  ) %>% 
  select(-a_0)

results <- bind_rows(bw_fixed, bw_random)


# Small sample (n1 = 10) --------------------------------------------------
n_1 <- 10
n_0 <- n_1 / k
delta <- (qnorm(prob0) + qnorm(power_want)) * sqrt(k+1) / sqrt(n_1 * lambda)

# Random variance case
mu_0_use <-  y_0h - rs * delta

bw_random <-
  fast_sim_numerical_wrapper_flex(
    mu_0 = mu_0_use,
    r = rs, 
    alpha = 0.1, 
    beta = 0.1, 
    n_1 = n_1, 
    k = k, 
    a_0 = 0.4, 
    delta = delta,
    alpha_0 = 0.05,
    gridsize = 200,
    tail_error = 1e-5
  )$results %>% 
  mutate(
    method = "BW_random",
    n_0 = n_1 / k
  ) %>% 
  select(-delta, -alpha, -beta, -a_0,-k)

# Fixed variance case
bw_fixed <-
  data.frame(
    n_1 = n_1, 
    n_0 = n_1 /k,
    mu_0 = mu_0_use,
    r = rs,
    a_0 = 0.4
  ) %>% 
  mutate(
    method = "BW_fixed",
    n_0 = floor(n_1 / k),
    power = 
      pnorm(
        (delta - a_0 * r * delta - qnorm(prob0) * sigma * sqrt(1/n_1 + (1-a_0)/(n_0))) /
          (sigma * sqrt(1/n_1 + (1-a_0)^2 /n_0 ))
      ),
    type_1_error = 
      pnorm(
        (- a_0 * r * delta - qnorm(prob0) * sigma * sqrt(1/n_1 + (1-a_0)/(n_0))) /
          (sigma * sqrt(1/n_1 + (1-a_0)^2 /n_0 ))
      ),
    n_0 = n_1 / k
  ) %>% 
  select(-a_0)

results <- bind_rows(results, bw_fixed, bw_random)

# Plot output -------------------------------------------------------------
colors <- c("#2c7bb6","#abd9e9")
smallfont <- 9
largefont <- 10

results %>% 
  select(method, mu_0,power, type_1_error, n_1, r) %>% 
  pivot_longer(
    cols = c(power, type_1_error),
               names_to = "metric",
               values_to = "performance"
  ) %>% 
  mutate(
    metric1 = ifelse(metric == "power", "Power", "T1E~Rate"),
    n_1_use = paste0("n[1]==",n_1),
    method1 = 
      factor(
        method, 
        levels = c("BW_fixed", "BW_random", "MAP", "SAM"),
        labels = c("BestWerc with Known Variance", "BestWerc with Estimated Variance", "MAP", "SAM")
      )
  ) %>% 
  ggplot(aes(x = r, y = performance, color = method1, linetype = method1)) + 
  geom_line(linewidth = 1) + 
  facet_grid(rows = vars(metric1), cols = vars(n_1), scales = "free_y",
             labeller = "label_parsed") + 
  scale_color_manual(values = colors) + 
  scale_linetype_manual(values = c(1, 2, 4, 5)) + 
  scale_x_continuous(breaks = round(seq(-0.6, 0.6, by = 0.2),1)) + 
  labs(x = "Relative Drift (r)") +
  theme_minimal() + 
  theme(
    plot.background = element_rect(color = NA, fill = "white"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", linewidth = 0.5),
    legend.title = element_blank(),
    legend.position = "top",
    axis.title.y = element_blank(),
    legend.spacing = unit(0, units = "pt"),
    legend.margin = margin(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.title.x = element_text(size = largefont),
    legend.text = element_text(size = largefont),
    axis.text = element_text(size = smallfont),
    strip.text = element_text(size = largefont)
  )

ggsave("figures/normal_compare_n.tiff", compression = "lzw",
       width = 6.5, height = 3.8, units = "in", dpi = 800)
# 
# results %>% 
#   filter(n_1 == 10) %>% 
#   select(r, method, power) %>% 
#   pivot_wider(names_from = "method", 
#               values_from = "power") %>% 
#   mutate(dif = BW_fixed - BW_random) %>% 
#   pull(dif) %>% 
#   summary()
