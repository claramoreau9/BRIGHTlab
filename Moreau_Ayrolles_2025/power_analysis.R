library(ggplot2)
library(dplyr)
library(pwr)
library(scales)  

# Define your effect sizes

n_range <- seq(2, 500, by = 1)

# Power analysis
df_power <- bind_rows(lapply(names(effect_sizes), function(measure) {
  data.frame(
    n1 = n_range, 
    n2 = n_range, 
    power = sapply(n_range, function(n) pwr.t2n.test(n1=n, n2=n, 
                                                     sig.level = 0.05/3, power = NULL, 
                                                     d = effect_sizes[[measure]], 
                                                     alternative = "two.sided")$power),
    measure = measure
  )
}))

# Plot power curves
ggplot(df_power, aes(x = n1, y = power, color = measure, group = measure)) +
  geom_line(linewidth = 1.2) +  
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "black") +  
  annotate("text", x = 460, y = 0.78, label = "80% Power", color = "black", size = 4, fontface = "bold") +  
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + 
  scale_color_manual(values = c("CT" = "skyblue", "SA" = "indianred", "Subcortical volume" = "mediumseagreen")) +  
  theme_minimal(base_size = 20) +
  labs(x = "Sample Size (per group)",
       y = "test power 1-ß",
       color = "Metric") +
  theme(plot.title = element_text(size = 20))
