# ==========================================
# TITLE: Parametric vs Non-Parametric Simulation Study
# AUTHOR: Valentine Joseph Owan
# R VERSION: 4.5.2
# DESCRIPTION:
# This script reproduces all simulation results, tables, and figures.
# ==========================================

# =====================
# 1. LOAD LIBRARIES
# =====================
library(dplyr)
library(ggplot2)
library(officer)
library(flextable)

set.seed(123)

# =====================
# 2. PARAMETERS
# =====================
sample_sizes <- c(5, 10, 15, 20, 25, 30, 50, 75, 100, 200, 500, 1000, 5000, 10000)
effect_sizes <- c(0, 0.2, 0.5)
distributions <- c("normal", "skewed", "heavy_tail", "flat")

get_simulations <- function(n) {
  if (n <= 30) return(5000)
  if (n <= 100) return(3000)
  if (n <= 500) return(1500)
  if (n <= 1000) return(800)
  return(300)
}

# =====================
# 3. DATA GENERATION
# =====================
generate_2group <- function(n, dist, effect) {
  if (dist == "normal") {
    x <- rnorm(n, 0, 1)
    y <- rnorm(n, effect, 1)
  } else if (dist == "skewed") {
    x <- rlnorm(n, 0, 0.5)
    y <- rlnorm(n, effect, 0.5)
  } else if (dist == "heavy_tail") {
    x <- rt(n, df = 3)
    y <- rt(n, df = 3) + effect
  } else {
    x <- runif(n, -1, 1)
    y <- runif(n, -1 + effect, 1 + effect)
  }
  list(x = x, y = y)
}

generate_3group <- function(n, dist, effect) {
  g1 <- generate_2group(n, dist, 0)$x
  g2 <- generate_2group(n, dist, effect)$y
  g3 <- generate_2group(n, dist, 2 * effect)$y
  list(g1 = g1, g2 = g2, g3 = g3)
}

safe_cor <- function(x, y) {
  if (sd(x) == 0 || sd(y) == 0) return(NA)
  cor(x, y)
}

# =====================
# 4. SIMULATION
# =====================
results <- list()
counter <- 1

for (n in sample_sizes) {
  n_sim <- get_simulations(n)

  for (dist in distributions) {
    for (effect in effect_sizes) {

      t_sig <- mw_sig <- agree2 <- 0
      t_only <- mw_only <- 0
      aov_sig <- kw_sig <- agree3 <- 0
      aov_only <- kw_only <- 0

      t_p <- mw_p <- aov_p <- kw_p <- numeric(n_sim)

      for (i in 1:n_sim) {

        # --- 2-group ---
        d2 <- generate_2group(n, dist, effect)
        tp <- t.test(d2$x, d2$y)$p.value
        mw <- wilcox.test(d2$x, d2$y)$p.value

        t_p[i] <- tp
        mw_p[i] <- mw

        td <- tp < 0.05
        md <- mw < 0.05

        t_sig <- t_sig + td
        mw_sig <- mw_sig + md

        if (td == md) {
          agree2 <- agree2 + 1
        } else {
          if (td) t_only <- t_only + 1
          if (md) mw_only <- mw_only + 1
        }

        # --- 3-group ---
        d3 <- generate_3group(n, dist, effect)
        values <- c(d3$g1, d3$g2, d3$g3)
        group <- factor(rep(1:3, each = n))

        ap <- summary(aov(values ~ group))[[1]][["Pr(>F)"]][1]
        kw <- kruskal.test(values ~ group)$p.value

        aov_p[i] <- ap
        kw_p[i] <- kw

        ad <- ap < 0.05
        kd <- kw < 0.05

        aov_sig <- aov_sig + ad
        kw_sig <- kw_sig + kd

        if (ad == kd) {
          agree3 <- agree3 + 1
        } else {
          if (ad) aov_only <- aov_only + 1
          if (kd) kw_only <- kw_only + 1
        }
      }

      results[[counter]] <- data.frame(
        SampleSize = n,
        Distribution = dist,
        EffectSize = effect,
        T_SigRate = t_sig / n_sim,
        MW_SigRate = mw_sig / n_sim,
        Agreement_2G = agree2 / n_sim,
        T_only = t_only / n_sim,
        MW_only = mw_only / n_sim,
        ANOVA_SigRate = aov_sig / n_sim,
        KW_SigRate = kw_sig / n_sim,
        Agreement_3G = agree3 / n_sim,
        ANOVA_only = aov_only / n_sim,
        KW_only = kw_only / n_sim,
        PvalCorr_2G = safe_cor(t_p, mw_p),
        PvalCorr_3G = safe_cor(aov_p, kw_p)
      )

      counter <- counter + 1
      cat("Done:", n, dist, effect, "\n")
    }
  }
}

results_df <- do.call(rbind, results)

write.csv(results_df, "FINAL_simulation_results.csv", row.names = FALSE)

# =====================
# 5. FIGURES
# =====================
fig_data <- results_df %>% filter(EffectSize == 0.5)

fig1 <- ggplot(fig_data, aes(SampleSize, Agreement_2G, color = Distribution)) +
  geom_line(linewidth = 1.2) + geom_point() + scale_x_log10() + theme_classic()

fig2 <- ggplot(fig_data, aes(SampleSize, Agreement_3G, color = Distribution)) +
  geom_line(linewidth = 1.2) + geom_point() + scale_x_log10() + theme_classic()

fig3 <- fig_data %>%
  mutate(diff = MW_SigRate - T_SigRate) %>%
  ggplot(aes(SampleSize, diff, color = Distribution)) +
  geom_line(linewidth = 1.2) + geom_point() + scale_x_log10() + theme_classic()

ggsave("Figure1_Agreement_2G.png", fig1, width = 8, height = 5, dpi = 300)
ggsave("Figure2_Agreement_3G.png", fig2, width = 8, height = 5, dpi = 300)
ggsave("Figure3_PowerDifference.png", fig3, width = 8, height = 5, dpi = 300)

# =====================
# 6. SESSION INFO
# =====================
sessionInfo()