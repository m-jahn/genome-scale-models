# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                  GENOME SCALE MODEL SIMULATIONS                    +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
# This pipeline imports functions to run FBA or FVA and return results
# in pd dataframes
#
#
# LIBRARIES ------------------------------------------------------------

import numpy as np
import pandas as pd
import tabulate
import cobra
import csv
import re
import flux_analysis as mfa

# wrap all code into main function
def main():
    
    # IMPORT SBML MODEL ------------------------------------------------
    
    model = cobra.io.read_sbml_model("sbml/RehMBEL1391_sbml_L3V1.xml")
    
    # ADDITIONAL CONSTRAINTS -------------------------------------------
    #
    # these are set to avoid long, thermodynamically impossible cycles
    # silence glyoxylate shunt, isocitrate lyase
    model.reactions.ICL.bounds = (0.0, 0.0)
    # restrict direction of glycine hydroxymethyltransferase; should make
    # glycine from serine, not other way around although it is thermodyn possible
    # otherwise an artificial loop can occur 
    model.reactions.GHMT3.bounds = (0.0, 1000.0)
    # threonine aldolase and dehydrogenase are part of threonine degradation, 
    # but carry flux on formate for glycine+acetyl-CoA synthesis. Unlikely 
    # that this is physiologically relevant
    model.reactions.THRA.bounds = (0.0, 0.0)
    model.reactions.THRD.bounds = (0.0, 0.0)
    # optionally silence acetate secretion during N limitation
    # (on excess fructose, FBA suggests excess frc uptake and fermentation
    # but it has been shown that cells rather make PHB)
    model.reactions.EX_ac_e.bounds = (0.0, 0.0)
    
    
    # DEFINE ENVIRONMENT -----------------------------------------------
    
    # we define the medium as Ralstonia eutropha minimal medium with the
    # following core components.
    # Among those are trace elements as well as H2O, protons, and Pi
    # in excess that can be used by the model.
    # Some components are particularly important:
    #  - the terminal e- acceptor, can be O2 or NO3, or other NO species
    #  - the carbon source 
    #  - the nitrogen source, in most cases NH4
    
    minimal_medium = {
        'EX_fe2_e': 10.0,
        'EX_mg2_e': 10.0,
        'EX_pi_e': 100.0,
        'EX_cobalt2_e': 10.0,
        'EX_cl_e': 10.0,
        'EX_k_e': 10.0,
        'EX_fe3_e': 10.0,
        'EX_so4_e': 10.0,
        'EX_fru_e': 0.0,
        'EX_nh4_e': 0.0,
        'EX_na_e': 10.0,
        'EX_o2_e': 18.5,
        'EX_mobd_e': 10.0,
        'EX_h2o_e': 1000.0,
        'EX_h_e': 100.0
    }
    
    
    # define a set of growth rates
    mu = [0.05, 0.1, 0.15, 0.2, 0.25]
    
    # and calculate the experimentally determined uptake rates 
    # for the limiting substrates
    qS_substrate = pd.DataFrame({
        'EX_fru_e': [mfa.qS_FRC(mu = i) for i in mu] + [0.0] * 10 + [10.0] * 5,
        'EX_succ_e': [0.0] * 5 + [mfa.qS_SUC(mu = i) for i in mu] + [0.0] * 10,
        'EX_formate_e': [0.0] * 10 + [mfa.qS_FOR(mu = i) for i in mu] + [0.0] * 5,
        'EX_nh4_e': [10.0] * 15 + [mfa.qS_NH4(mu = i) for i in mu]
        })
    
    
    # RUN SIMULATIONS IN LOOP
    # --------------------------------------------------------------
    for index, row in qS_substrate.iterrows():
        
        # add limiting substrate uptake rate to minimal medium
        mm = minimal_medium.copy()
        for s in row.index.to_list():
            mm[s] = row[s]
        
        
        # run FBA analysis
        # --------------------------------------------------------------
        result_fba = mfa.run_FBA(
            model = model.copy(),
            medium = mm
            )
        
        
        # runs FSA analysis ('flux sampling')
        # --------------------------------------------------------------
        if index in [4, 9, 14, 19]:
            result_fsa = mfa.run_flux_sampling(
                model.copy(), mm,
                n_samples = 100, 
                objective = "Biomass", 
                print_summary = True,
                FVA_optimality = 0.95)
        
        
        # Export result tables
        suffix = ''
        for i in row.keys():
            suffix = suffix + '_' + re.sub('^EX_|_e$', '', i) + '_' + str(round(row[i], 2))
        result_fba['FBA'].to_csv('simulations/' + "FBA" + suffix + ".csv")
        if index in [4, 9, 14, 19]:
            result_fsa.to_csv('simulations/' + "FSA" + suffix + ".csv")


if __name__== "__main__" :
    main()
