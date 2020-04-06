# loading libraries
library(lattice)
library(latticeExtra)
library(tidyverse)
library(Rtools)
library(stringi)


# +++++++++++++ LOAD MODEL SIMULATION DATA +++++++++++++++++++++++++++++++++++++

# generalized function to load and combine data from multiple result tables
read_flux_results <- function(files) {
  lapply(files, function(filename) {
    read_csv(filename) %>%
    mutate(simulation = 
      stri_extract_first(filename, 
        regex = '(succ|fru|for).*nh4_e_[0-9]*.[0-9]*') %>% 
      gsub("_e|_EX", "", .))
  }) %>% 
  
  # combine list of files
  dplyr::bind_rows() %>%
  
  # rearrange result columns
  rename(reaction = X1) %>%
  separate(simulation, 
    into = c("carbon_source", "qS_carbon", "nitrogen_source", "qS_nitrogen"), 
    sep = "_") %>%
  mutate(qS_carbon = as.numeric(qS_carbon), qS_nitrogen = as.numeric(qS_nitrogen))

}

# list of FBA simulation result files to load
fba_dir = "simulations/FBA_growth_constrained/"
fva_dir = "simulations/FVA_growth_constrained/"
df_fba <- read_flux_results(list.files(fba_dir, pattern = "_FBA.csv$", full.names = TRUE))
df_fva <- read_flux_results(list.files(fva_dir, pattern = "_FVA.csv$", full.names = TRUE))


# +++++++++++++ YIELD AND UPTAKE RATES +++++++++++++++++++++++++++++++++++++++++

# add uptake rates to determine yield (Y = µ / qS = 
# g_bm * gDCW * h / g_subs * gDCW *h)
df_yield <- df_fba %>% filter(reaction == "EX_BIOMASS_c") %>%
  
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
  
  
# +++++++++++++ PLOT YIELD +++++++++++++++++++++++++++++++++++++++++++++++++++++

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

# +++++++++++++ PLOT MAJOR FLUXES (FBA) ++++++++++++++++++++++++++++++++++++++++

# for a set of FBA simulations with different growth rate and different substrate
# limitations, plot the highest (exchnage) fluxes per condition
# if, optionally, growth rate is constrained to experimentally determined rates
# we can see where additional carbon or energy are dissipated apart from biomass
#
# First get an overview about conditions
df_fba %>% group_by(carbon_source, qS_carbon) %>% 
  summarize(n_reactions = length(reaction))

# prepare data and plot FBA fluxes for exchange reactions
plot_exchange <- df_fba %>%
  
  # select only exchnage reactions
  filter(
    grepl("^EX_", reaction),
    qS_carbon != 10,
    loopless != 0) %>%
  
  # add a nominal instead of numeric column for qS to allow
  # direct comparison
  group_by(carbon_source, reaction) %>% 
  mutate(qS_carbon_level = factor(qS_carbon) %>% as.numeric) %>%
  
  # sort by average flux per reaction, descending
  group_by(reaction) %>%
  mutate(average_flux = mean(loopless)) %>% 
  arrange(desc(average_flux)) %>%
  
  xyplot(loopless ~ factor(reaction, unique(reaction)) | carbon_source, .,
    groups = qS_carbon_level,
    par.settings = custom.lattice(), pch = 19,
    ylim = c(-10, 10), 
    ylab = expression("flux "*"[mmol h"^-1*" gDCW"^-1*"]"),
    xlab = "reaction",
    layout = c(1, 3), 
    as.table = TRUE, between = list(x = 0.5, y = 0.5),
    scales = list(alternating = FALSE, x =list(rot = 35, cex = 0.7)),
    panel = function(x, y, ...) {
      panel.grid(h = -1, v = -1, col = grey(0.9))
      panel.abline(v = 1:60, col = grey(0.9))
      panel.abline(h = 0, col = grey(0.2))
      panel.barplot(x, y, origin = 0, beside = TRUE, ewidth = 0.05, ...)
    }
  )


svg("simulations/flux_exchange_reactions.svg", width = 8, height = 6)
print(plot_exchange)
dev.off()


# +++++++++++++ FLUX VARIABILITY ANALYSIS (FVA) ++++++++++++++++++++++++++++++++

# FVA tells us how much variation is allowed to happen in the flux of certain
# reactions, without compromising the objective to e.g. 99 % of the optimal solution

svg("simulations/flux_variability_0_95.svg", width = 8, height = 6)
xyplot(minimum + maximum ~ factor(reaction, unique(reaction)) | carbon_source, 
  filter(df_fva, qS_carbon != 10),
  par.settings = custom.lattice(), pch = 19, 
  ylim = c(-6, 6), 
  ylab = expression("flux "*"[mmol h"^-1*" gDCW"^-1*"]"),
  xlab = "reaction",
  layout = c(1, 3), as.table = TRUE, between = list(x = 0.5, y = 0.5),
  scales = list(alternating = FALSE, x =list(rot = 35, cex = 0.7)),
  panel = function(x, y, ...) {
    panel.grid(h = -1, v = -1, col = grey(0.9))
    panel.rect(0, -6, 5.5, -4, col = grey(0.85), border = NA)
    panel.rect(5.5, -6, 11.5, -4, col = "white", border = grey(0.85))
    panel.rect(11.5, -6, 18.5, -4, col = grey(0.85), border = NA)
    panel.rect(18.5, -6, 27.5, -4, col = "white", border = grey(0.85))
    panel.rect(27.5, -6, 36.5, -4, col = grey(0.85), border = NA)
    panel.rect(36.5, -6, 38.5, -4, col = "white", border = grey(0.85))
    panel.text(c(3, 9, 15, 23, 32, 37), rep(-5, 6), labels = c("ED", "EMP", "CBB", "PYR", "TCA", "GLX"))
    panel.abline(h = 0, lwd = 1.5, lty = 2, col = grey(0.5))
    panel.xyplot(x, y, ...)
    panel.key(corner = c(0.98,0.95), ...)
  }
)
dev.off()
