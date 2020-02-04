# 
# Simple test script to import SBML (systems biology markup
# language) models that describe the metabolism of a cell.
# This script performs the following basic steps
# 
#  - import cobra and other required libraries
#  - import an SBML model
#  - check key data of the model
#  - run a basic flux balance analysis (FBA)
#
# All of these functions are directly taken from the COBRApy
# documentation.
#
# To install cobrapy for python3, navigate to your favored directory and 
# execute the following lines in terminal:
#
# git clone https://github.com/opencobra/cobrapy.git
# cd cobrapy-master/
# pip3 install -e
# 
# 
# LIBRARIES ------------------------------------------------------------

import numpy as np
import pandas as pd
import tabulate
import cobra
import cobra.test
import os
from os.path import join
import csv
import re


# IMPORT SBML MODEL ----------------------------------------------------

# set the path to the model's directory
data_dir = '/home/michael/Documents/SciLifeLab/Resources/Models/genome-scale-models/Ralstonia_eutropha/'
model = cobra.io.read_sbml_model(join(data_dir, "RehMBEL1391_sbml_L2V1.xml"))
# alternatively load test model for comparison
#ecoli = cobra.test.create_test_model("textbook")
salmonella = cobra.test.create_test_model()

# summary of the imported model
print('%i reactions' % len(model.reactions))
print('%i metabolites' % len(model.metabolites))
print('%i genes' % len(model.genes))


# DEALING WITH THE ENVIRONMENT -----------------------------------------
#
# the model contains two different sets of environmental variables
# that are important for simulation of growth:
# A: the model has exchange reactions that simply
#    model the extracellular availability (flux) towards the cell;
#    these have bounds = (0, 1000) when they 'produce' the metabolite
#    or e.g. (-1000, 0) when they 'consume' it
# B: model.medium is the mirror of the source exchange reaction
#    the medium has positive flux, i.e. produces a compound

# optionally export the default medium (all exchange reactions) to modify it
#with open(join(data_dir, "growth_medium.csv"), "w") as f:
#    for key in model.medium.keys():
#        f.write("%s, %s\n"%(key, model.medium[key]))

# bounds are not correctly set for all reactions during import, only
# set to standard (-1000, 1000). Many reactions are then reversible
# which is wrong. Set bounds based on info in 'notes' field, and
# set all exchange reactins to c(0, 1000))
for reaction in model.reactions:
    if reaction in model.exchanges:
        reaction.lower_bound = 0.0
    elif reaction.notes['REVERSIBILITY'] == 'true':
        reaction.lower_bound = -1000.0
    elif reaction.notes['REVERSIBILITY'] == 'false':
        reaction.lower_bound = 0.0

# check that it worked
model.exchanges.list_attr("bounds")
model.reactions.query("[A-Z]+t").list_attr("bounds")


# RUN FBA --------------------------------------------------------------

# set a growth medium: read in modified medium 
#Re_minimal_medium = pd.read_csv(join(data_dir, "growth_medium.csv"), index_col = "reaction")
#Re_minimal_medium = Re_minimal_medium.to_dict()["flux"]
#model.medium = Re_minimal_medium
model.medium = {
    'EX_mg2_e': 10.0,
    'EX_pi_e': 100.0,
    'EX_cobalt2_e': 10.0,
    'EX_cl_e': 10.0,
    'EX_k_e': 10.0,
    'EX_fe3_e': 10.0,
    'EX_so4_e': 10.0,
    'EX_fru_e': 5.0,
    'EX_nh4_e': 10.0,
    'EX_na_e': 10.0,
    'EX_o2_e': 18.5,
    'EX_mobd_e': 10.0,
    'EX_h2o_e': 1000.0,
    'EX_h_e': 100.0
    }

# set objective function
model.objective = {model.reactions.Biomass: 1}
# run FBA analysis
solution = model.optimize()

# save modified model as SBML L3 V1
cobra.io.write_sbml_model(model, join(data_dir, "RehMBEL1391_sbml_L3V1.xml"))

# REPORTING ------------------------------------------------------------
#
# print solution summary the status from the linear programming solver
print([solution, "status: ", solution.status])

# report ideal minimal medium and, indirectly, also the growth-limiting
# component of the medium (visible by component's upper_bound)
min_medium = cobra.medium.minimal_medium(
    model, 
    solution.objective_value,
    minimize_components = 4)
print(min_medium)


# export fluxes as csv
fluxes = solution.fluxes.sort_values()
# print top 20 forward and backward flux
print(fluxes[0:20])
print(fluxes[len(fluxes)-20:len(fluxes)])
#fluxes.to_csv(join(data_dir, "fluxes.csv"))

# quick summary of FBA analysis
print(model.summary())

# summary of energy balance
print(model.metabolites.atp_c.summary())

# summary of redox balance
print(model.metabolites.nadh_c.summary())

# table with shadow prices (how much does objective change by unit change of constraint?)
print(solution.shadow_prices.sort_values())


# RUN FVA --------------------------------------------------------------
#
# Flux Variability Analysis (FVA) displaying variability
# of reactions for solutions of e.g. at least 95% of optimal solution,
# and using loopless solving

cobra.flux_analysis.flux_variability_analysis(
    model, model.reactions.query("[A-Z]+t"),
    loopless = True,
    fraction_of_optimum=0.95
    )


# RUN pFBA -------------------------------------------------------------
#
# pFBA is parsimonious FBA and tries to eliminate through loops of 
# reactions. A fast and pragmatic approach is to post-process flux 
# distributions to simply set fluxes to zero wherever they can be 
# zero without changing the fluxes of any exchange reactions in the 
# model ('CycleFreeFlux')

solution = cobra.flux_analysis.loopless_solution(model)
