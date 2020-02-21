# loading libraries
library(lattice)
library(latticeExtra)
library(tidyverse)
library(Rtools)


# +++++++++++++ LOAD MODEL DATA ++++++++++++++++++++++++++++++++++++++++++++++++

# list of raw model result files to load
files <- list.files("simulations", pattern = "*.csv$", full.names = TRUE)

# Generalized function to load and combine data from multiple result tables
df <- lapply(files, function(filename) {
  read_csv(filename) %>%
  mutate(simulation = gsub("simulations.|EX_|_e|.csv", "", filename))
  }) %>% 
  dplyr::bind_rows() %>%
  rename(metabolite = X1) %>%
  separate(simulation, 
    into = c("carbon_source", "qS_carbon", "nitrogen_source", "qS_nitrogen"), 
    sep = "_") %>%
  mutate(qS_carbon = as.numeric(qS_carbon), qS_nitrogen = as.numeric(qS_nitrogen))


# +++++++++++++ YIELD AND UPTAKE RATES +++++++++++++++++++++++++++++++++++++++++

# add uptake rates to determine yield (Y = µ / qS = 
# g_bm * gDCW * h / g_subs * gDCW *h)
df_yield <- df %>% filter(metabolite == "EX_BIOMASS_c") %>%
  
  # add molecular weights to re-calculate qS to g/g*L
  mutate(
    MW_carbon = recode(carbon_source, formate = 0.04603, fru = 0.18016, succ = 0.11809),
    MW_nitrogen = 0.05349,
    qS_carbon_g = qS_carbon*MW_carbon,
    qS_nitrogen_g = qS_nitrogen*MW_nitrogen
  ) %>%
  
  mutate(
    yield_carbon = loopless/qS_carbon_g,
    yield_nitrogen = loopless/qS_nitrogen_g
  )
  
  
# +++++++++++++ PLOTTING +++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Herbet-Pirt-plot function to create plots using different parameters
HP_plot <- function(data, xvar, yvar, condvar, 
  xlimits = c(0, 0.3), ylimits = c(0, 1)
) {
  xyplot(get(yvar) ~ get(xvar) | get(condvar), data,
    par.settings = custom.lattice(), lwd = 1.5,
    xlim = xlimits, ylim = ylimits,
    ylab = expression("q"[S]*" [g h"^-1*" gDCW"^-1*"]"),
    xlab = expression('µ [h'^'-1'*']'),
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.xyplot(x, y, cex = 0.9, ...)
      # regression line through linear part of the data
      panel.lmlineq(x, y, fontfamily = "FreeSans", 
        pos = 3, offset = 7, r.squared = TRUE, cex = 0.7, col.text = grey(0.3), ...)
      coef <- lm(y ~ x, data.frame(x, y))$coeff
      # displaying maintenance and yield coefficients
      panel.abline(h = coef[[1]], lty = 2, ...)
      panel.text(0.15, coef[[1]], 
        paste("ms =", round(coef[[1]], 3), "g h-1 g_DCW-1"), 
        col = grey(0.3), pos = 3, cex = 0.7)
      panel.text(0.15, coef[[1]], paste(expression("Yx/S ="), 
          round(1/coef[[2]], 3), "g_DCW g_S-1"), 
        col = grey(0.3), pos = 1, cex = 0.7)
    }
  )
}


svg("simulations/qS_vs_growth_rate.svg", width = 6.7, height = 6.7)
print(HP_plot(data = subset(df_yield, carbon_source == "formate"), 
  xvar = "loopless", yvar = "qS_carbon_g", condvar = "carbon_source",
  ylimits = c(-5/4, 5)), split = c(1,1,2,2), more = TRUE)

print(HP_plot(data = subset(df_yield, carbon_source == "fru" & qS_carbon != 10),
  xvar = "loopless", yvar = "qS_carbon_g", condvar = "carbon_source",
  ylimits = c(-1/4, 1)), split = c(2,1,2,2), more=TRUE)

print(HP_plot(data = subset(df_yield, carbon_source == "fru" & qS_carbon == 10), 
  xvar = "loopless", yvar = "qS_nitrogen_g", condvar = "nitrogen_source",
  ylimits = c(-0.23/4, 0.15)), split = c(1,2,2,2), more = TRUE)

print(HP_plot(data = subset(df_yield, carbon_source == "succ"),
  xvar = "loopless", yvar = "qS_carbon_g", condvar = "carbon_source",
  ylimits = c(-1/4, 1)), split = c(2,2,2,2))
dev.off()