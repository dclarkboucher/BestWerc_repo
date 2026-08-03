library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
source("source/binary_helpers.R")
R_0 <- c(0.2, 0.4, 0.6, 0.8)
n_1 <- 40
k <- 2
baseline_power <- 0.65
alpha_0 <- 0.05
R_1 <- c()
for (i in seq_along(R_0)){
  R_1[i] <- find_R1(
    power = baseline_power,
    R0 = R_0[i],
    n_1 = n_1, 
    n_0 = n_1/k
  )
  
}
delta <- R_1 - R_0
a_0 <- seq(0, 0.96, by = 0.02)
r <- c(-0.5, -0.3, 0, 0.3, 0.5)
settings <- 
  expand_grid(
    n_0 = n_1 / k,
    n_1 = n_1, 
    scenario = seq_along(R_0),
    a_0 = a_0, 
    r = r
  ) %>% 
  mutate(
    R_0 = R_0[scenario],
    delta = delta[scenario]
  ) %>% 
  select(-scenario)

out <-
  binary_power_vectorized(
    R_0_list = settings$R_0,
    delta_list = settings$delta,
    n_1_list = settings$n_1,
    n_0_list = settings$n_0,
    a_0_list = settings$a_0,
    r_list = settings$r,
    alpha_0 = 0.05, 
    gridsize = 100, 
    Delta = 0
  )
results <- cbind(settings, power = out[["power"]], 
                 type_1_error = out[["type_1_error"]])

colors1 <- c("#d7191c","#fdae61","gray60","#abd9e9","#2c7bb6")

colors2 <- c("#fdae61","gray60","#abd9e9")
smallfont <- 9
largefont <- 10

# results %>% 
#   pivot_longer(
#     cols = c(power, type_1_error),
#     names_to = "scenario",
#     values_to = "performance"
#   ) %>% 
#   mutate(
#     scenario = ifelse(scenario == "power", "H1", "H0"),
#     rf = factor(r),
#     metric1 = fct_rev(ifelse(scenario == "H1", "Power", "T1E Rate")),
#     R_0_label = paste0("pi[0] == ",R_0)
#   ) %>% 
#   ggplot(aes(x = a_0, y = performance, color = rf, linetype = metric1
#   ))  + 
#   geom_hline(aes(yintercept = 0.05), color = "gray30",
#              linewidth = 0.65) +
#   geom_hline(aes(yintercept = 0.65), color = "gray30",
#              linewidth = 0.65) +
#   facet_wrap(~R_0_label, labeller = "label_parsed") + 
#   geom_line(linewidth = 0.8) +
#   # geom_point(size = 2) + 
#   scale_x_continuous(breaks = round(seq(0, 1, 0.2),2), limits = c(0,1)) + 
#   scale_y_continuous(breaks = round(seq(0, 1, 0.2),2), limits = c(0,1)) + 
#   scale_color_manual(values = rev(colors1))+
#   scale_linetype_manual(values = c(5,4)) +
#   labs(linetype = "Metric",color = "r",
#        
#        y = "Rejection Probability",
#        x = expression(a[0])) + 
#   theme_minimal() + 
#   theme(
#     plot.background = element_rect(color = NA, fill = "white"),
#     panel.grid = element_line(color = NA),
#     panel.border = element_rect(color = "black", linewidth = 0.5,
#                                 fill = NA),
#     legend.position = "right",
#     legend.box = "vertical",
#     legend.margin = margin(),
#     strip.text = element_text(size = smallfont),
#     # strip.background = element_rect(fill = "gray90", color = "gray90"),
#     plot.title = element_blank(),
#     axis.text = element_text(size = smallfont),
#     axis.title = element_text(size = smallfont),
#     axis.title.y = element_blank(),
#     axis.ticks = element_line(color = "black", linewidth = 0.5),
#     legend.title = element_text(size = smallfont),
#     legend.text = element_text(size = smallfont)
#   )
# ggsave("figures/error_rates_binary.tiff",
#        width = 6.5, height = 3.5, units = "in", dpi = 800
# )


results %>% 
  pivot_longer(
    cols = c(power, type_1_error),
    names_to = "scenario",
    values_to = "performance"
  ) %>% 
  mutate(
    scenario = ifelse(scenario == "power", "H1", "H0"),
    rf = factor(r),
    metric1 = fct_rev(ifelse(scenario == "H1", "Power", "T1E Rate")),
    R_0_label = paste0("Control~event~rate*':'~pi[0] == ",R_0)
  ) %>% 
  ggplot(aes(x = a_0, y = performance, linetype = rf, color = metric1
  ))  + 
  geom_hline(aes(yintercept = 0.05), color = "black",
             linewidth = 0.65) +
  geom_hline(aes(yintercept = 0.65), color = "black",
             linewidth = 0.65) +
  facet_wrap(~R_0_label, labeller = "label_parsed") + 
  geom_line(linewidth = 0.8) +
  # geom_point(size = 2) + 
  scale_x_continuous(breaks = round(seq(0, 1, 0.2),2), limits = c(0,1)) + 
  scale_y_continuous(breaks = round(seq(0, 1, 0.2),2), limits = c(0,1)) + 
  scale_linetype_manual(
    values = c(1,5,2,4,3)
  ) +
  scale_color_manual(values = c("#abd9e9","#2c7bb6")) +
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
    strip.text = element_text(size = smallfont),
    # strip.background = element_rect(fill = "gray90", color = "gray90"),
    plot.title = element_blank(),
    axis.text = element_text(size = smallfont),
    axis.title = element_text(size = smallfont),
    axis.title.y = element_blank(),
    legend.key.width = unit(0.5,"in"),
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
  )
ggsave("figures/error_rates_binary.tiff",
       width = 6.5, height = 6, units = "in", dpi = 800,
       compression = "lzw"
)

results %>%
  mutate(
    n_0h = n_0 / (1-a_0) * a_0,
    r_0h = delta * r + R_0,
    bias = (1 + n_0h*r_0h  - R_0 *(2 + n_0h)) / (2 + n_0h + n_0),
    mse = bias^2 + n_0 * R_0 * (1 - R_0) / (n_0 + n_0h + 2)^2
  ) %>% 
  pivot_longer(
    cols = c(power, type_1_error, bias, mse),
    names_to = "metric",
    values_to = "performance"
  ) %>% 
  mutate(
    performance = as.character(format(round(performance,3),3))
  ) %>% 
  select(R_0, r, a_0, performance, metric) %>% 
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
  readr::write_csv(paste0("tables/error_rates_binary.csv"))


