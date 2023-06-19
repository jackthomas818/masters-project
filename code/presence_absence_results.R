# Load packages
library(tidyverse)
library(latex2exp)

# Load simulation summaries
sim12 <- read_csv("summary_stats_SimStudy12.csv") %>%
  mutate(s_d = sigma * delta)

# Parameter summaries
sim12 %>% 
  group_by(N_actual, tau, nsites, ntraps, sigma, delta) %>%
  summarize(n = n())

# Plot estimates by presence-absence detects
sim12 %>%
  filter(N_actual == 500) %>%
  mutate(p = factor(p)) %>%
  ggplot(aes(x = s_d, y = N_mean, colour = p)) + 
  geom_point() + 
  facet_grid(tau ~ .) + 
  geom_hline(aes(yintercept = N_actual)) +
  xlab(TeX("$\\sigma\\delta$")) +
  ylab("Mean Estimated Population Size") + 
  scale_y_continuous(breaks = c(250,500))

dev.copy2pdf(file="Figures/presence_absence_1.pdf", width = 6, height = 5)
dev.off()

