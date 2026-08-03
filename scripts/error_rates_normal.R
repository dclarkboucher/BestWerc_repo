library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(forcats)
library(readr)
source("source/normal_helpers.R")

a_0s <- seq(0,0.96, by = 0.02)
rs <- c(-0.5, -0.3, 0, 0.3, 0.5)
ks <- 2
lambda <- 1
sigma <- 1/sqrt(lambda)
n_1 <- 40
prob0 <- 0.95
power1 <- 0.8
power2 <- 0.65
Delta <- 0
mu_0 <- 0
delta <- (qnorm(prob0) + qnorm(power1)) * sqrt(2) / sqrt(n_1 * lambda)
updates <- c("prior_mean")
alpha_betas <- c(0.1, 2, 20)

# Data frame of settings
settings <-
  expand_grid(
    mu_0 = mu_0,
    alpha_beta = alpha_betas,
    r = rs, 
    n_1 = n_1,
    a_0 = a_0s,
    k = ks,
    delta = delta
  )

# First do random variance case
out_random <-
  fast_sim_numerical_wrapper_flex(
    mu_0 = settings$mu_0,
    alpha = settings$alpha_beta,
    beta = settings$alpha_beta,
    n_1 = settings$n_1,
    r = settings$r,
    a_0 = settings$a_0,
    delta = settings$delta,
    k = settings$k,
    alpha_0 = 0.05,
    gridsize = 250,
    tail_error = 1e-5
  )$results %>% 
  as_tibble() %>% 
  mutate(type = "random")

# Now do fixed variance case
out_fixed <-
  settings %>% 
  filter(alpha_beta == alpha_beta[1]) %>% 
  mutate(
    type = "fixed",
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
      )
  ) %>% 
  select(-n_0)

results <-
  bind_rows(out_fixed, out_random) %>% 
  mutate(
    type_label = case_when(
      type == "fixed" ~ "plain(Fixed~variance)",
      T ~ paste0("{alpha == beta} == ",alpha)) %>%
        fct_inorder(ordered = TRUE),
    bias = a_0 * r * delta,
    mse = bias^2 + (1 - a_0)^2 * (sigma^2) * k / n_1
    
  ) %>% 
  pivot_longer(cols = c(power, type_1_error, bias, mse), 
               names_to = "metric",
               values_to = "performance") %>% 
  mutate(
    scenario = ifelse(metric == "power", "H1", "H0"),
    
    rf = factor(as.character(r),
                levels = as.character(c(-0.5, -0.3, 0, 0.3, 0.5)),
                labels = as.character(c(-0.5, -0.3, 0, 0.3, 0.5)),
                ordered = TRUE
    ),
    metric1 = (ifelse(metric == "power", "Power", "T1E Rate"))  )

colors1 <- c("#d7191c","#fdae61","gray60","#abd9e9","#2c7bb6")
colors2 <- c("gold1","dodgerblue1","firebrick1")
smallfont <- 9
largefont <- 10

# Main plot
results %>%
  filter(metric %in% c("power", "type_1_error")
  ) %>% 
  ggplot() +
  geom_hline(aes(yintercept = 0.05), color = "black",
             linewidth = 0.65) +
  geom_hline(aes(yintercept = 0.65), color = "black",
             linewidth = 0.65) +
  scale_y_continuous(limits = c(0,1),
                     expand = c(0.02,0.02),
                     breaks = seq(0, 1, by = 0.2)
  ) +
  scale_x_continuous(limits =c(0,1),
                     expand = c(0.02,0.02),
                     breaks = seq(0, 1, by = 0.2)
                     
  ) +
  geom_line(
    aes(x = a_0,
        y = performance,
        color = metric1,
        linetype = rf
        
    ),
    linewidth = 0.8
  ) +
  theme_minimal() +
  scale_linetype_manual(
    values = c(1,5,2,4,3)
  ) +
  scale_color_manual(values = c( "#2c7bb6","#abd9e9")) +
  facet_wrap(~type_label, dir = "v", labeller = "label_parsed", 
             ncol = 2) +
  theme(
    panel.grid = element_line(color = NA),
    panel.border = element_rect(color = "black", linewidth = 0.5,
                                fill = NA),
    plot.background = element_rect(color = NA, fill = "white"),
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(),
    legend.key.width = unit(0.5,"in"),
    strip.text = element_text(size = smallfont),
    plot.title = element_blank(),
    axis.text = element_text(size = smallfont),
    axis.title = element_text(size = smallfont),
    axis.title.y = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.title = element_text(size = smallfont),
    legend.text = element_text(size = smallfont)
  ) +
  annotate(
    geom = "text",
    x = 0.95,
    y = 0.075,
    label = "0.05",
    size = 3
  ) +
  annotate(
    geom = "text",
    x = 0.95,
    y = 0.68,
    label = "0.65",
    size = 3
  ) +
  labs(
    # x = expression(a[0]),
    x = expression(paste("Borrowing fraction (", a[0],")")),
    y = "Rejection Probability",
    linetype = "Relative\nDrift (r)",
    color = "Metric"
    
  )
ggsave(
  file = "figures/error_rates_normal.tiff",
  units = "in",
  height = 6,
  width = 6.5,
  dpi = 800,
  compression = "lzw"
)

results %>% 
  select(type,alpha,r, a_0, performance, metric) %>% 
  mutate(performance = as.character(format(round(performance,3),3))
  ) %>%
  filter(
    round(a_0,2) %in% c(0, 0.2, 0.4, 0.6, 0.8)
  ) %>% 
  mutate(
    metric = factor(metric, levels = c("power", "type_1_error","bias","mse"), 
                    labels = c("Power", "T1E Rate","Bias","MSE"))
  ) %>% 
  group_by(r) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  pivot_wider(names_from = r, values_from = performance) %>% 
  select(-row) %>% 
  write_csv(paste0("tables/error_rates_normal.csv"))
