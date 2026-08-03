rm(list=ls())
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(readr)
source("source/normal_helpers.R")

mu_0h <- -0.71
prior_mean <- -0.71
prior_sd <- 0.107
delta <- -0.88

r <- c(-0.24, -0.2, -0.1, 0, 0.1, 0.2, 0.24)
mu_0s <- mu_0h - r * delta 
y_0h <- prior_mean

# Generate table to choose borrowing fraction -----------------------------
r_test <- c(-0.24, 0,  0.24)
mu_0_test <-  mu_0h - r_test * delta 
alpha_0 <- 0.05
sigma <- 1.5
beta <- 2 * sigma^2
alpha <- 2
a_0 <- 0.05
n_1 <- 36
n_0 <- 18

fixed <-
  expand.grid(
    mu_0 = mu_0_test,
    n_1 = n_1,
    n_0 = n_0, 
    a_0 = seq(0, .8, 0.1)
  ) %>% 
  mutate(
    n_0h = n_0 * a_0 / (1-a_0),
    r = (mu_0h - mu_0) / delta,
    power = pnorm(
      (-delta + a_0 * (y_0h - mu_0) - qnorm(1-alpha_0) * sigma * sqrt(1/n_1 + (1-a_0)/(n_0))) /
        (sigma * sqrt(1/n_1 + (1-a_0)^2 /n_0 )) # note: signs are swapped to change direction
    ),
    t1e = pnorm(
      (a_0 * (y_0h - mu_0) - qnorm(1-alpha_0) * sigma * sqrt(1/n_1 + (1-a_0)/(n_0))) /
        (sigma * sqrt(1/n_1 + (1-a_0)^2 /n_0 ))
    ),
    method = "fixed"
  ) %>% 
  select(method, a_0, r, power, t1e) %>% 
  pivot_longer(cols = c(power, t1e), names_to = "metric", values_to = "performance") %>%
  mutate(
    performance = as.character(format(round(performance, 3), 3))
  ) %>% 
  pivot_wider(names_from = c(r, metric), values_from = performance)


random <- 
  fast_sim_numerical_wrapper(
    mu_0 = mu_0h,
    alpha = alpha, beta = beta, 
    n_1 = n_1, 
    k = 2, 
    r = r_test, 
    a_0 = seq(0, .8, 0.1), 
    delta = -delta, 
    alpha_0 = alpha_0,
    gridsize = 250,
    tail_error = 1e-5
  )$results %>% 
  mutate(
    method = "random"
  ) %>% 
  rename(t1e = type_1_error) %>% 
  select(method, a_0, r, power, t1e) %>% 
  pivot_longer(cols = c(power, t1e), names_to = "metric", values_to = "performance") %>% 
  mutate(
    performance = as.character(format(round(performance, 3), 3))
  ) %>% 
  pivot_wider(names_from = c(r, metric), values_from = performance)

design_options <- bind_rows(fixed, random)


# Make a similar table for the paper. 
design_options %>% 
  pivot_longer(cols = c(contains("t1e"), contains("power")),
               names_to = "r_setting",
               values_to = "performance") %>% 
  mutate(
    r = str_remove(r_setting,"_.*"),
    setting = str_remove(r_setting, paste0(r,"_"))
  ) %>% 
  select(-r_setting) %>% 
  pivot_wider(
    names_from = c(method, r),
    values_from = performance
  ) %>% 
  arrange(
    a_0, setting
  ) %>% 
  write_csv("tables/application.csv")
  

# Now, make a detailed plot of the designs meeting our requirements
a_0s_want <-
  design_options %>% 
  pivot_longer(cols = c(contains("t1e"), contains("power")),
               names_to = "r_setting",
               values_to = "performance") %>% 
  mutate(
    r = str_remove(r_setting,"_.*"),
    setting = str_remove(r_setting, paste0(r,"_"))
  ) %>% 
  select(-r_setting) %>% 
  pivot_wider(names_from = setting, values_from = performance) %>% 
  group_by(a_0) %>% 
  summarize(
    enough_power = all(power >= 0.70),
    low_error = all(t1e <= 0.10)
  ) %>% 
  ungroup() %>% 
  filter(enough_power & low_error) %>% 
  pull(a_0)

best_results <-
  fast_sim_numerical_wrapper(
    mu_0 = 0, 
    alpha = alpha, 
    beta = beta, 
    n_1 = n_1, 
    k = 2, 
    r = seq(-0.24, 0.24, by = 0.01),
    a_0 = a_0s_want,
    delta = -delta, 
    alpha_0 = alpha_0,
    gridsize = 250,
    tail_error = 1e-5
  )

smallfont <- 9
largefont <- 10

best_results$results %>% 
  pivot_longer(
    cols = c(power, type_1_error),
    names_to = "setting", 
    values_to = "performance"
  ) %>% 
  mutate(
    setting = ifelse(setting == "power", "Power", "T1E Rate")
  ) %>% 
  ggplot(aes(x = r, y = performance, 
             color = setting,
             linetype = factor(a_0))) + 
  geom_hline(aes(yintercept = 0.10), color = "gray30",
             linewidth = 0.65) +
  geom_hline(aes(yintercept = 0.70), color = "gray30",
             linewidth = 0.65) +
  geom_line(linewidth = 0.8) + 
  scale_color_manual(values = c( "#2c7bb6","#abd9e9")) +
  scale_linetype_manual(values = c(2,4,5,1)) + 
  annotate(
    geom = "text",
    x = -0.1,
    y = 0.66,
    label = "Minimum power: 0.70",
    size = 3
  ) +
  annotate(
    geom = "text",
    x = -0.1,
    y = 0.14,
    label = "Maximum T1E rate: 0.10",
    size = 3
  ) + 
  labs(
    linetype =expression(paste("Borrowing\nfraction", a[0]),":"),
    y = "Rejection Probability",
    x = "Relative Drift (r)",
    color = "Scenario:"
  ) + 
  facet_wrap(~setting) + 
  scale_y_continuous(
    breaks = seq(0, 1, 0.2),
    limits = c(0,1)
  ) + 
  theme_minimal() + 
  theme(
    axis.title.y = element_blank(),
    panel.grid = element_line(color = NA),
    panel.border = element_rect(color = "black", linewidth = 0.5,
                                fill = NA),
    plot.background = element_rect(color = NA, fill = "white"),
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(),
    strip.text = element_text(size = smallfont),
    # strip.background = element_rect(fill = "gray90", color = "gray90"),
    plot.title = element_blank(),
    axis.text = element_text(size = smallfont),
    axis.title = element_text(size = smallfont),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.title = element_text(size = smallfont),
    legend.text = element_text(size = smallfont)
  )



best_results$results %>% 
  pivot_longer(
    cols = c(power, type_1_error),
    names_to = "setting", 
    values_to = "performance"
  ) %>% 
  mutate(
    setting = ifelse(setting == "power", "Power", "T1E Rate")
  ) %>% 
  ggplot(aes(x = r, y = performance, 
             color = factor(a_0),
             linetype = factor(a_0))) + 
  geom_line(linewidth = 0.6) + 
  # scale_color_manual(values = c( "#2c7bb6","#abd9e9")) +
  scale_color_ordinal(
    option = "F", direction = -1,
    begin = 0.25, end = 0.75
  ) + 
  scale_linetype_manual(
    values = c(4,2,5,1)
    
  ) + 
  scale_y_continuous(n.breaks = 6) + 
  labs(
    linetype = expression(paste("Borrowing\nfraction", a[0],":")),
    y = "Rejection Probability",
    x = "Relative Drift (r)",
    color = expression(paste("Borrowing\nfraction", a[0],":"))
  ) + 
  facet_wrap(~setting, scales = "free_y") + 
  theme_minimal() + 
  theme(
    axis.title.y = element_blank(),
    panel.grid = element_line(color = NA),
    panel.border = element_rect(color = "black", linewidth = 0.5,
                                fill = NA),
    plot.background = element_rect(color = NA, fill = "white"),
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(),
    strip.text = element_text(size = smallfont),
    # strip.background = element_rect(fill = "gray90", color = "gray90"),
    plot.title = element_blank(),
    axis.text = element_text(size = smallfont),
    axis.title = element_text(size = smallfont),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    legend.title = element_text(size = smallfont),
    legend.text = element_text(size = smallfont)
  )


ggsave(
  file = "figures/application.tiff",
  units = "in",
  height = 3,
  width = 5,
  dpi = 800,
  compression = "lzw"
)

