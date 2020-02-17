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
import csv


# IMPORT SBML MODEL ----------------------------------------------------

# set the path to the model's directory
wd = '/home/michael/Documents/SciLifeLab/Resources/Models/genome-scale-models/Ralstonia_eutropha/'
model = cobra.io.read_sbml_model(wd + "sbml/RehMBEL1391_sbml_L3V1.xml")


# DEFINE ENVIRONMENT ---------------------------------------------------

# we define the medium as Ralstonia eutropha minimal medium with the
# following core components.
# Among those are trace elements as well as H2O, protons, and Pi
# in excess that can be used by the model.
# Some components are particularly important:
#  - the terminal e- acceptor, can be O2 or NO3, or other NO species
#  - the carbon source, we can test 
#  - the nitrogen source, mostly NH4

minimal_medium = {
    'EX_mg2_e': 10.0,
    'EX_pi_e': 1000.0,
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
    loopless = cobra.flux_analysis.loopless_solution(model)
    fluxes = pd.DataFrame(
        dict(loopless = loopless.fluxes)
    )
    
    if print_summary:
        print("SOLUTION: OBJECTIVE =" + str(solution.objective_value))
        print("TOP 6 CONSUMING FLUXES:")
        print(fluxes.filter(like = 'EX_', axis = 0).sort_values("loopless").head(6))
        print("TOP 3 PRODUCING FLUXES:") 
        print(fluxes.filter(like = 'EX_', axis = 0).sort_values("loopless").tail(3))
    
    # return data frame with flux results
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
    
    # add limiting substrate uptake rate to minimal medium
    mm = minimal_medium.copy()
    mm[row['carbon_source']] = row['qS_carbon']
    mm[row['nitrogen_source']] = row['qS_nitrogen']
    
    # run FBA analysis on a copy of the model
    result = run_FBA(model = model.copy(), medium = mm)
    
    # construct filename
    filename = (
        row['carbon_source'] + '_' + 
        str(round(row['qS_carbon'], 2)) + '_' +
        row['nitrogen_source'] + '_' +
        str(round(row['qS_nitrogen']))
    )
    
    # save result from pandas data frame to hdd
    result.to_csv(wd + 'simulations/' + filename + '.csv')


# ANALYZE ONE SPECIFIC SOLUTION ----------------------------------------
#
# run FBA analysis on a copy of the model
mm = minimal_medium.copy()
mm['EX_fru_e'] = qS_FRC(0.25)#10
mm['EX_nh4_e'] = 10#qS_NH4(0.1)

model_test = model.copy()
model_test.medium = mm
model_test.optimize()

print(model_test.metabolites.fru_e.summary())
print(model_test.metabolites.co2_c.summary())
print(model_test.metabolites.atp_c.summary())
print(model_test.metabolites.nadh_c.summary())
print(model_test.metabolites.get_by_id("d6pgc_c").summary())


# Results / Interpretation:
#  - on formate, flux through PPC (phosphoenolpyruvate carboxylase)
#    is much higher than through CBB cycle (32.3/4.65 = 6.95)
#  - we know from Alegasan et al., 2018, that ratio should be
#    PPC/CBB = 0.97 / 2.79 = 0.35 (under mixotrophic growth with gly)
#  - Next question is, what is the yield on formate in the model?
#    Our own measurements reach 0.1 gDCW g-1, or 0.18 Cmol Cmol-1.
#    This is in good agreement with Grunwald et al., who hit 0.16 
#    Cmol Cmol-1
#
#    Yield on formate in model is Y = µ / qS = 
#    g_bm * gDCW * h / g_subs * gDCW *h.
#    Yield is much higher than expected. Also, no CO2 is released but
#    all CO2 is recycled by CBB or PPC reaction which is unrealistic:
#    there's only 1 NADH per CO2 produced but at least 3 required per CO2
#    for Calvin cycle to work.
print("Yield [g_bm / g_S] = " + 
    str(model_test.objective.value / qS_FRC(0.1, mmol = False)))

# Other notes to the current model implementation:
# - Rubisco and R15BP are missing: Calvin cycle is implemented as a 
#   lumped reaction CBBCYCLE
# - TCA cycle running incomplete although was reported to be functional (Alegasan paper)
# - main reasons for that is fumarate and succ-coa being replenished from amino acids 
#   (anaplerotic reactions running reverse). How realistic is that?
# - formate dehydrogenase (FDH) is making formic acid with NADH under heterotrophic conditions
#   (fructose or succinate as substrate). 
# - CO2 is never emitted but always recaptured by the CBB cycle in the model simulations.
#   some CO2 should theoretically be emitted though.  
