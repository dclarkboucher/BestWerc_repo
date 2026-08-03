library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)
Rcpp::sourceCpp("source/survival_power.cpp")
gamma_0 <- c(1,2) # Average survival time
alpha_0 <- 0.05
censoring_prob <- c(0.40, 0.60)
n_1 <- 200
n_0 <- 100
set.seed(1)
settings <- 
  expand_grid(gamma_0 = gamma_0, censoring_prob = censoring_prob) %>% 
  mutate(
    gamma_c = gamma_0 * (1-censoring_prob) / censoring_prob,
    gamma_1 = 
      map2_dbl(
        gamma_0, gamma_c,
        ~find_g1(power = 0.65, gamma_0 = .x, gamma_c = .y, 
                 n_1 = n_1, n_0 = n_0, Delta = 0, alpha_0 = alpha_0, reps = 10000,
                 gridsize = 100)
      ),
    delta = gamma_1 - gamma_0,
    n_1 = n_1, 
    n_0 = n_0
  ) %>% 
  expand_grid(
    r = c(-0.5, -0.3, 0, 0.3, 0.5),
    a_0 = seq(0, 0.96, by = 0.02)
  )

set.seed(1)
system.time({
  out <-
    survival_power_vectorized(
      gamma_0_list = settings$gamma_0, 
      delta_list = settings$delta, 
      gamma_c_list = settings$gamma_c,
      n_1_list = settings$n_1, 
      n_0_list = settings$n_0,
      a_0_list = settings$a_0,
      r_list = settings$r,
      Delta = 0,
      alpha_0 = 0.05,
      reps = 10000,
      gridsize = 100,
      numerical = FALSE,
      fix_cens = FALSE
    )
})
results <- cbind(settings, as.data.frame(out))

colors1 <- c("#d7191c","#fdae61","gray60","#abd9e9","#2c7bb6")
smallfont <- 9
largefont <- 10
results %>% 
  pivot_longer(
    cols = c(power, type_1_error),
    names_to = "scenario",
    values_to = "performance"
  ) %>% 
  mutate(
    scenario = ifelse(scenario == "power", "H1", "H0"),
    metric = ifelse(scenario == "H1","Power", "T1E Rate"),
    censoring = paste0(censoring_prob * 100,"*'%'~Censoring"), # the star separates tokens without space
    # censoring = paste0("Control censoring: ", censoring_prob),
    g0_label = paste0("Control~mean~survival~time*':'~gamma[0]==",gamma_0),
    rf = factor(r)
  ) %>% 
  ggplot(aes(x = a_0, y = performance, linetype = rf, 
             # shape = scenario
             color = metric, shape = scenario
  ))  + 
  geom_hline(aes(yintercept = 0.05), color = "gray30",
             linewidth = 0.65) +
  geom_hline(aes(yintercept = 0.65), color = "gray30",
             linewidth = 0.65) +
  facet_nested_wrap( 
    vars(g0_label, censoring), labeller = "label_parsed"
  ) + 
  geom_line(linewidth = 0.8) +
  scale_x_continuous(breaks = round(seq(0, 1, 0.2),2),limits =c(0,1)) + 
  scale_y_continuous(breaks = round(seq(0, 1, 0.2),2),limits =c(0,1)) + 
  scale_linetype_manual(
    values = c(1,5,2,4,3)
  ) +
  scale_color_manual(values = c( "#2c7bb6","#abd9e9")) +
  labs(
    x = expression(paste("Borrowing fraction (", a[0],")")),
    y = "Rejection Probability",
    linetype = "Relative\nDrift (r)",
    color = "Metric"
  ) + 
  theme_minimal() + 
  theme(
    plot.background = element_rect(color = NA, fill = "white"),
    panel.grid = element_line(color = NA),
    panel.border = element_rect(color = "black", linewidth = 0.5,
                                fill = NA),
    legend.position = "right",
    legend.box = "vertical",
    legend.margin = margin(),
    legend.key.width = unit(0.5,"in"),
    strip.text = element_text(size = smallfont),
    # strip.background = element_rect(fill = "gray90", color = "gray90"),
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
    y = 0.675,
    label = "0.65",
    size = 3
  )
ggsave("figures/error_rates_survival.tiff", compression = "lzw",
       width = 6.5, height = 6, units = "in", dpi = 800
)

# add bias and mse
bias_mse <-
  survival_bias_mse_vectorized(
    gamma_0_list = results$gamma_0, 
    delta_list = results$delta, 
    gamma_c_list = results$gamma_c, 
    n_0_list = results$n_0, 
    a_0_list = results$a_0, 
    r_list = results$r
  ) %>% 
  as.data.frame()

results1 <- 
  cbind(results, bias_mse) %>% 
  mutate(
    bias = mean - gamma_0, 
    mse = var + bias^2
  ) %>% 
  select(-mean, -var)


results1 %>%
  pivot_longer(
    cols = c(power, type_1_error, bias, mse),
    names_to = "metric",
    values_to = "performance"
  ) %>% 
  mutate(
    performance = as.character(format(round(performance,3),3))
  ) %>% 
  select(gamma_0, censoring_prob, r, a_0, performance, metric) %>% 
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
  readr::write_csv(paste0("tables/error_rate_survival.csv"))
