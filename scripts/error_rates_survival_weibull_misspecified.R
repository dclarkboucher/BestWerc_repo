library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(ggh4x)
Rcpp::sourceCpp("source/survival_power.cpp")
# Settings
alpha_0 <- 0.05
n_1 <- 200
n_0 <- 100

# These simulations imagine that the Weibull shape parameter is wrongly assumed to 
# equal 1. Thus, contrary to the simulations error_rates_survival_weibull.R, p_0
# is not calibrated to control the baseline T1E rate; instead, it is fixed at 0.95.

setting_dat <-
  expand_grid(
    avg_surv_time = c(1,2),
    shape = c(0.8,1.2),
    cens_prob = c(0.4)
  ) %>%
  mutate(
    lambda = avg_surv_time / gamma(1+1/shape),
    lambda_cens = ((lambda^shape) *(1-cens_prob) / (cens_prob))^(1/shape),
    delta = c(0.39373,0.39373,0.78743,0.78743),
    alpha_0 = c(0.05, 0.05, 0.05, 0.05),
    setting = 1:n()

  )
# 
# setting_dat <-
#   expand_grid( 
#     avg_surv_time = c(1,2),
#     shape = c(1),
#     cens_prob = c(0.4)
#   ) %>%
#   mutate(
#     lambda = avg_surv_time / gamma(1+1/shape),
#     lambda_cens = ((lambda^shape) *(1-cens_prob) / (cens_prob))^(1/shape),
#     delta = c(0.39373,0.78743),
#     alpha_0 = c(0.05, 0.05),
#     setting = 1:n()
#     
#   )


reps <- 5000 # UPDATE THIS TO 5000 BEFORE SUBMITTING. 
settings <- seq_len(nrow(setting_dat))
rs <- c(-0.5, -0.3, 0, 0.3, 0.5)
# rs <- 0
scenarios <- c(0,1)
a_0s <- seq(0, 0.9, by = 0.1)
# a_0s <- 0


out_list <- list()
j <- 1
for (set in settings){
  print(set)
  lambda <- setting_dat$lambda[set]
  lambda_cens <- setting_dat$lambda_cens[set]
  delta <- setting_dat$delta[set]
  shape <- setting_dat$shape[set]
  alpha_0 <- setting_dat$alpha_0[set]
  mu_0 <- setting_dat$avg_surv_time[set]
  
  for (r in rs){
    for (scen in scenarios){
      
      mu_0h <- mu_0 + delta * r
      mu_1 <- mu_0 + delta * scen
      
      lambda0 <- mu_0 / gamma(1+1/shape)
      lambda1 <- mu_1 / gamma(1+1/shape)
      
      for( rep in seq_len(reps)){
        
        # Generate data
        set.seed(rep)
        lambda_x0 <- (lambda0^(-shape) + lambda_cens^(-shape)) ^ (-1/shape)
        lambda_x1 <- (lambda1^(-shape) + lambda_cens^(-shape)) ^ (-1/shape)
        x0_dat <- rweibull(n_0, shape = shape, scale = lambda_x0)
        x0 <- sum(x0_dat)
        x1_dat <- rweibull(n_1, shape = shape, scale = lambda_x1)
        x1 <- sum(x1_dat)
        u0 <- rbinom(1,n_0, prob = (lambda_cens^shape) / (lambda_cens^shape + lambda0^shape))
        u1 <- rbinom(1,n_1, prob = (lambda_cens^shape) / (lambda_cens^shape + lambda1^shape))
        
        for (a_0 in a_0s){
          
          n_0h <- n_0 * (a_0) / (1-a_0);
          # Determine posterior distribution for this borrowing fraction (exponential assumption)
          alpha0 <- n_0h + u0
          beta0 <- n_0h * mu_0h + x0
          alpha1 <- u1
          beta1 <- x1
          prob <- gamma_diff_prob(
            0, alpha0, beta0, alpha1, beta1, 1000, 1e-4
          )
          
          estimate <- beta0 / alpha0
          decision <- as.numeric(prob > 1 - alpha_0)
          out_list[[j]] <-
            data.frame(
              setting = set,
              r = r,
              scenario = scen,
              rep = rep,
              a_0 = a_0,
              decision = decision,
              estimate = estimate
            )
          
          j <- j + 1
          
        }
      }
    }
  }
}


results <- 
  bind_rows(out_list) %>% 
  group_by(setting, scenario, r, a_0) %>% 
  summarize(
    performance = mean(decision),
    mean = mean(estimate),
    sd = sd(estimate)
  ) %>% 
  mutate(
    scenario = ifelse(scenario == 1, "H1", "H0")
  ) %>% 
  left_join(setting_dat, by = "setting")
saveRDS(results,"results/error_rate_survival_weibull_misspecified.rds")


head(results)
smallfont <- 9
largefont <- 10
results %>% 
  mutate(
    metric = ifelse(scenario == "H1","Power", "T1E Rate"),
    censoring = paste0(cens_prob * 100,"*'%'~Censoring"),
    shape = paste0("Shape*':'~",as.character(shape)),
    g0_label = paste0("Control~mean~survival~time*':'~",avg_surv_time),
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
    vars(g0_label, shape), labeller = "label_parsed"
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
ggsave("figures/error_rates_survival_weibull_misspecified.tiff", compression = "lzw",
       width = 6.5, height = 6, units = "in", dpi = 800
)

results1 <-
  results %>% 
  ungroup() %>% 
  mutate(
    bias = mean - avg_surv_time,
    mse = bias^2 + sd^2
  ) %>% 
  select(scenario, r, a_0, performance, bias, mse, avg_surv_time, shape) %>% 
  pivot_wider(names_from = "scenario", values_from = "performance") %>% 
  rename(power = H1, type_1_error = H0)


results1 %>%
  pivot_longer(
    cols = c(power, type_1_error, bias, mse),
    names_to = "metric",
    values_to = "performance"
  ) %>% 
  mutate(
    performance = as.character(format(round(performance,3),3))
  ) %>% 
  select(avg_surv_time, shape, r, a_0, performance, metric) %>% 
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
  readr::write_csv(paste0("tables/error_rate_survival_weibull_misspecified.csv"))

