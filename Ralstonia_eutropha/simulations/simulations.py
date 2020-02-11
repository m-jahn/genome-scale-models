#
# GENOME SCALE MODEL SIMULATIONS
# for Ralstonia eutropha (aka Cupriavidus necator), a soil bacterium
# and model organism for PHB synthesis and lithoautotrophic growth
# (CO/CO2 as carbon source, hydrogen H2 as energy source) 
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


# IMPORT SBML MODEL ----------------------------------------------------

# set the path to the model's directory
data_dir = '/home/michael/Documents/SciLifeLab/Resources/Models/genome-scale-models/Ralstonia_eutropha/'
model = cobra.io.read_sbml_model(join(data_dir, "RehMBEL1391_sbml_L3V1.xml"))


# DEFINE ENVIRONMENT ---------------------------------------------------

# we define the medium as Ralstonia eutropha minimal medium with the
# following core components.
# Among those are trace elements as well as H2O, protons, and Pi
# in excess that can be used by the model.
# Some components are particularly important:
#  - the terminal e- acceptor, can be O2 or NO3, or other NO species
#  - the carbon source, we can test 
#  - the nitrogen source, mostly NH4

model.medium = {
    'EX_mg2_e': 10.0,
    'EX_pi_e': 10.0,
    'EX_cobalt2_e': 10.0,
    'EX_cl_e': 10.0,
    'EX_k_e': 10.0,
    'EX_fe3_e': 10.0,
    'EX_so4_e': 10.0,
    'EX_na_e': 10.0,
    'EX_o2_e': 18.5,
    'EX_mobd_e': 10.0,
    'EX_h2o_e': 1000.0,
    'EX_h_e': 100.0
}

# C-SOURCE
# --------
# medium does not take concentration, but requires
# an estimate of the substrate uptake rate in mmol / (gDCW x h).
# In-house measurements show the following substrate uptake rates 
# depending on dilution rate/growth rate
# 
# FRC:   qS [g / (gDCW x h)] = -0.016 + 2.246 * mu
#        qS [mmol / (gDCW x h)] = (-0.016 + 2.246 * mu) / 0.18016 g/mmol
# SUC:   qS [g / (gDCW x h)] =  0.003 + 2.246 * mu
#        qS [mmol / (gDCW x h)] = ( 0.003 + 2.246 * mu) / 0.11809 g/mmol
# FOR:   qS [g / (gDCW x h)] =  0.211 + 10.046 * mu
#        qS [mmol / (gDCW x h)] = (0.211 + 10.046 * mu) / 0.04603 g/mmol
# NH4CL: qS [g / (gDCW x h)] =  -0.014 + 0.346 * mu
#        qS [mmol / (gDCW x h)] = (-0.014 + 0.346 * mu) / 0.05349 g/mmol
# 
# These functions return uptake rate qS of the limiting substrate in 
# [mmol / (gDCW x h)] or [g / (gDCW x h)] for a given growth rate mu
def qS_FRC(mu, mmol = True):
    qS = (-0.016 + 2.246 * mu)
    if mmol: qS = qS / 0.18016
    return qS

def qS_SUC(mu, mmol = True):
    qS = ( 0.003 + 2.246 * mu)
    if mmol: qS = qS / 0.11809
    return qS

def qS_FOR(mu, mmol = True):
    qS = (0.211 + 10.046 * mu)
    if mmol: qS = qS / 0.04603
    return qS

def qS_NH4(mu, mmol = True):
    qS = (-0.014 + 0.346 * mu)
    if mmol: qS = qS / 0.05349
    return qS


# RUN SINGLE FBA SIMULATION --------------------------------------------
# 
# define medium
mm = model.medium
mm['EX_fru_e'] = qS_FRC(mu = 0.01)
mm['EX_nh4_e'] = 10
mm['EX_succ_e'] = 0#qS_SUC(mu = 0.25)
mm['EX_pi_e'] = 1000
model.medium = mm

# run FBA
#model.objective = {model.reactions.Biomass: 1}
solution = model.optimize()
print(model.summary())

# quick summary of latest FBA analysis run
# summary of energy balance
print(model.metabolites.atp_c.summary())
# summary of redox balance
print(model.metabolites.nadh_c.summary())
# 
print(model.metabolites.o2_c.summary())


# RUN SET OF SIMULATIONS -----------------------------------------------


# function to run an FBA analysis and report result in a 
# pandas data frame
# Optimizing the Biomass function is the default. The biomass 
# concentration is (per definition) balanced to 1 g/mmol, so that
# Biomass reaction flux in mmol / (gDCW x h) equals g / (gDCW x h), equals
# growth rate. This means we can also determine the yield in 
# g biomass / g substrate
def run_FBA(
    model, medium, 
    objective = "Biomass", 
    obj_maximize = True,
    print_summary = True):
    
    # assign medium to model
    model.medium = medium
    
    # assign objective function
    if obj_maximize: obj_maximize = 1
    else: obj_maximize = -1
    model.objective = {
    model.reactions.get_by_id(objective): obj_maximize
    }
    
    # run FBA analysis
    solution = model.optimize()
    if print_summary:
        print([solution, "status: ", solution.status])
        print(model.summary())
    
    # collect fluxes and return
    fluxes = solution.to_frame()
    return fluxes



# define a set of uptake rates for the limiting substrates
mu = [0.05, 0.1, 0.15, 0.2, 0.25]

qS_substrate = pd.DataFrame(
    [qS_FRC(mu = i) for i in mu] +
    [qS_SUC(mu = i) for i in mu] +
    [qS_FOR(mu = i) for i in mu] +
    [10] * 5,
    columns = ['qS_carbon'])

qS_substrate['carbon_source'] = (
    ['EX_fru_e'] * len(mu) +
    ['EX_succ_e'] * len(mu) +
    ['EX_formate_e'] * len(mu) +
    ['EX_fru_e'] * len(mu))

qS_substrate['qS_nitrogen'] = (
    [10] * 5 +
    [10] * 5 +
    [10] * 5 +
    [qS_NH4(mu = i) for i in mu])

qS_substrate['nitrogen_source'] = (
    ['EX_nh4_e'] * len(mu) +
    ['EX_nh4_e'] * len(mu) +
    ['EX_nh4_e'] * len(mu) +
    ['EX_nh4_e'] * len(mu))


# run several simulations in a loop
for index, row in qS_substrate.iterrows():
    
    # we execute FBA in a with statement in order to keep model
    # object unchanged
    # add limiting substrate uptake rate to minimal medium
    with model:
        medium = model.medium
        medium[row['carbon_source']] = row['qS_carbon']
        medium[row['nitrogen_source']] = row['qS_nitrogen']
        
        # run FBA analysis
        result = run_FBA(model, medium)
        
        # save result from pandas data frame to hdd
        result.to_csv(data_dir + 'simulations/sim_' + str(index) + '.csv')


