# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                FUNCTIONS FOR FLUX ANALYSIS IN GEMS                 +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
# Simulations with FBA, pFBA, FVA and so on, for Ralstonia eutropha 
# (aka Cupriavidus necator), a soil bacterium and model organism 
# for PHB synthesis and lithoautotrophic growth
# (CO/CO2 as carbon source, hydrogen H2 as energy source)
#
#
# LIBRARIES ------------------------------------------------------------

import numpy as np
import pandas as pd
import tabulate
import cobra
import csv


# FUNCTIONS TO DETERMINE UPTAKE RATES ----------------------------------
#
# Medium does not include concentration, but requires
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


# FLUX BALANCE ANALYSIS (FBA) ------------------------------------------

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
    print_summary = True,
    FVA = None,
    FVA_optimality = 0.99):
    
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
    result = {'FBA': fluxes}
    
    # optionally run FVA analysis and add results to df
    if FVA != None:
        fva = cobra.flux_analysis.flux_variability_analysis(
            model, 
            [model.reactions.get_by_id(id) for id in FVA],
            loopless = True,
            fraction_of_optimum = FVA_optimality
        )
        # add FVA result to dict
        result['FVA'] = fva
    
    # print summary of FBA results to terminal
    if print_summary:
        print("SOLUTION: OBJECTIVE = " + str(round(solution.objective_value, 4)))
        print("TOP 5 CONSUMING FLUXES:")
        print(fluxes.filter(like = 'EX_', axis = 0).sort_values("loopless").head(5))
        print("TOP 5 PRODUCING FLUXES:") 
        print(fluxes.filter(like = 'EX_', axis = 0).sort_values("loopless").tail(5))
    
    # return data frame with flux results
    return result


# FLUX SAMPLING ANALYSIS -----------------------------------------------

# function to run an flux sampling analysis and report result in a 
# pandas data frame. Flux sampling is performed using the OptGP
# sampler and includes only flux distributions that are within a certain
# range of feasible growth rates (such as at least 95% of mu max)

def run_flux_sampling(
    model, medium,
    n_samples = 10,
    objective = "Biomass", 
    obj_maximize = True,
    print_summary = True,
    FVA_optimality = 0.9,
    ):
    
    # assign medium to model
    model.medium = medium
    
    # assign objective function
    if obj_maximize: obj_maximize = 1
    else: obj_maximize = -1
    model.objective = {
    model.reactions.get_by_id(objective): obj_maximize
    }
    
    # run FBA analysis to predetermine optimal solution
    solution = model.optimize()
    bm_lb = solution.objective_value*FVA_optimality
    model.reactions.Biomass.bounds = (bm_lb, 1000)
    
    # run flux sampling
    result = cobra.sampling.sample(model, n = n_samples)
    
    # print summary of FBA results to terminal
    if print_summary:
        print("SOLUTION: OBJECTIVE = " + str(round(solution.objective_value, 4)))
        print("TOP 5 CONSUMING /  PRODUCING FLUXES:\n")
        print(result.mean().sort_values())
    
    # return data frame with flux results
    return result
