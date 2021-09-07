# loading libraries
library(lattice)
library(latticeExtra)
library(latticetools)
library(tidyverse)
library(stringi)

# set working directory
setwd("/home/michael/Documents/SciLifeLab/Resources/Models/genome-scale-models/Ralstonia_eutropha/")


# +++++++++++++ LOAD MODEL SIMULATION DATA +++++++++++++++++++++++++++++++++++++

# generalized function to load and combine data from multiple result tables
read_flux_results <- function(files) {
  lapply(files, function(filename) {
    read_csv(filename, col_types = cols()) %>%
    mutate(simulation = str_replace(filename, ".*F[BS]A_", "") %>%
      str_replace(".csv", ""))
  }) %>% 
  # combine list of files
  dplyr::bind_rows() %>%
  # rearrange result columns
  rename(reaction = X1) %>%
  separate(simulation, 
    into = c("frc", "qS_frc", "suc", "qS_suc", "form", "qS_for", "nh4", "qS_nh4"), 
    sep = "_") %>%
  select(-frc, -suc, -form, -nh4) %>%
  mutate(across(.cols = starts_with("qS"), as.numeric))
}

# list of FBA simulation result files to load
df_fba <- read_flux_results(list.files("simulations/", pattern = "^FBA_fru", full.names = TRUE))


# +++++++++++++ YIELD AND UPTAKE RATES +++++++++++++++++++++++++++++++++++++++++

# add uptake rates to determine yield (Y = µ / qS = 
# g_bm * gDCW * h / g_subs * gDCW *h)
df_summary <- df_fba %>% filter(reaction == "EX_BIOMASS_c") %>%
  # rerrange order
  arrange(qS_frc, qS_suc, qS_for) %>%
  mutate(condition = rep(c("formate", "succinate", "fructose", "ammonium"), each = 5)) %>%
  # add molecular weight in g/mmol to re-calculate qS to g/g*L
  mutate(
    MW = rep(c(0.04603, 0.11809, 0.18016, 0.05349), each = 5),
    qS_mmol_gDCW_h = c(qS_for[1:5], qS_suc[6:10], qS_frc[11:15], qS_nh4[16:20]),
    qS_g_gDCW_h = qS_mmol_gDCW_h * MW,
    yield = loopless/qS_g_gDCW_h
  )
  
  HP_plot <- function(data, xlimits = c(0, 0.3), ylimits = c(0, 1)) {
  xyplot(qS_g_gDCW_h ~ loopless | condition, data,
    par.settings = custom.colorblind(), lwd = 1.5,
    xlim = xlimits, ylim = ylimits,
    ylab = expression("q"[S]*" [g h"^-1*" gDCW"^-1*"]"),
    xlab = expression('µ [h'^'-1'*']'),
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.xyplot(x, y, cex = 0.9, ...)
      # regression line through linear part of the data (omit µ = 0.25/h)
      sel <- x != 0.25; x = x[sel]; y = y[sel]
      panel.lmline(x, y, ...)
      coef <- lm(y ~ x, data.frame(x, y))$coeff
      # displaying maintenance and yield coefficients
      panel.abline(h = coef[[1]], lty = 2, ...)
      panel.text(0.15, coef[[1]],
        paste("ms =", round(coef[[1]], 3), "g h-1 gDCW-1"),
        col = grey(0.3), pos = 3, cex = 0.7)
      panel.text(0.15, coef[[1]], paste(expression("Yx/S ="),
        round(1/coef[[2]], 3), "gDCW g-1"),
        col = grey(0.3), pos = 1, cex = 0.7)
    }
  )
}

print(HP_plot(data = filter(df_summary, condition == "formate"),
  ylimits = c(-7/4, 7)), position = c(0.017,0.5,0.5,1), more = TRUE)
print(HP_plot(data = filter(df_summary, condition == "fructose"),
  ylimits = c(-1.5/4, 1.5)), position = c(0.5,0.5,1,1), more = TRUE)
print(HP_plot(data = filter(df_summary, condition == "ammonium"),
  ylimits = c(-0.25/4, 0.25)), position = c(-0.017,0,0.5,0.5), more = TRUE)
print(HP_plot(data = filter(df_summary, condition == "succinate"),
  ylimits = c(-1.5/4, 1.5)), position = c(0.5,0,1,0.5))

