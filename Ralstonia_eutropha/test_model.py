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

# check some components contained in the model manually
model.reactions[1:3]
model.metabolites[1:3]
model.genes[1:3]

# DEALING WITH THE ENVIRONMENT -----------------------------------------
#
# the model contains two different sets of environmental variables
# that are important for simulation of growth:
# A: the model has exchange reactions (prefixed "EX_") that simply
#    model the extracellular availability (flux) towards the cell;
#    these have bounds = (0, 1000) when they 'produce' the metabolite
#    or e.g. (-1000, 0) when they are 'consume' it
# B: the model.medium is the mirror of the source exchange reaction
#    the medium has positive flux, i.e. produces a compound

# export the default medium (all exchange reactions) to modify it
with open(join(data_dir, "growth_medium.csv"), "w") as f:
    for key in model.medium.keys():
        f.write("%s, %s\n"%(key, model.medium[key]))

# all the model's reaction bounds were populated with (-1000, 1000)
# this means that flux can go in both directions
model.reactions.list_attr("bounds")

# we need to change bounds for most exchange reactions to (0, 1000)
# in order to make them sinks for the cell. Can do this in a loop
for reaction in model.reactions.query("EX_"):
    #if reaction.id != "EX_BIOMASS_c":
    reaction.lower_bound = 0

# check that it worked, except for biomass
model.reactions.query("EX_").list_attr("lower_bound")


# RUN FBA --------------------------------------------------------------

# set a growth medium: read in modified medium 

# set objective function
model.objective = {model.reactions.Biomass: 1}   
# run FBA analysis
solution = model.optimize()
print(solution)

# print solution summary
print(solution)
print(solution.f)                   # the value of the objective function
print(solution.status)              # the status from the linear programming solver
print(solution.y_dict)              # a dict with shadow prices indexed by reaction identifier (how much does objective change by unit change of constraint?)
print(solution.x_dict)              # a dict with fluxes indexed by reaction identifier
                                    # the flux for a reaction variable is the difference of the primal values for the forward and reverse reaction variables


