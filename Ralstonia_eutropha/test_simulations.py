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
import flux_analysis as mfa

# wrap all code into main function
def main():
    
    # IMPORT SBML MODEL ----------------------------------------------------
    
    model = cobra.io.read_sbml_model("sbml/RehMBEL1391_sbml_L3V1.xml")
    
    # --- set additional constraints ---
    # silence glyoxylate shunt, isocitrate lyase
    #model.reactions.ICL.bounds = (0.0, 0.0)
    
    
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
    
    # define a set of growth rates
    mu = [0.05, 0.1, 0.15, 0.2, 0.25]
    
    # and calculate the experimentally determined uptake rates 
    # for the limiting substrates
    qS_substrate = pd.DataFrame(
        [mfa.qS_FRC(mu = i) for i in mu] +
        [mfa.qS_SUC(mu = i) for i in mu] +
        [mfa.qS_FOR(mu = i) for i in mu] +
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
        [mfa.qS_NH4(mu = i) for i in mu])
    
    qS_substrate['nitrogen_source'] = (
        ['EX_nh4_e'] * len(mu) +
        ['EX_nh4_e'] * len(mu) +
        ['EX_nh4_e'] * len(mu) +
        ['EX_nh4_e'] * len(mu))
    
    qS_substrate['growth_rate_exp'] = mu*4
    
    
    # sets of reactions for important pathways
    reactions_ED = ['PGI', 'G6PDH', 'PGL', 'EDD', 'EDA']
    reactions_EMP = ['FBA', 'TPI', 'GA3PD', 'PGK', 'PGM', 'ENO']
    reactions_CBB = ['TRKT1', 'TRKT2', 'TRADL', 'RPE', 'RPI', 'PRUK', 'RBPC']
    reactions_pyr = ['PYK', 'PPC', 'PYC', 'PPCK', 'PDH1', 'PDH2', 'ME1', 'ME2']
    reactions_TCA = ['CITS', 'ACONT1', 'ACONT2', 'ICITD', 'ICITDp', 'AKGDH', 'SUCOAS', 'SUCCD3', 'FUMR', 'MDH1']
    reactions_glx = ['ICL', 'MALS']
    reactions_FVA = (
        reactions_ED + reactions_EMP + 
        reactions_CBB + reactions_pyr + 
        reactions_TCA + reactions_glx
    )
    
    
    # run several simulations in a loop
    for index, row in qS_substrate.iterrows():
        
        # add limiting substrate uptake rate to minimal medium
        mm = minimal_medium.copy()
        mm[row['carbon_source']] = row['qS_carbon']
        mm[row['nitrogen_source']] = row['qS_nitrogen']
        
        # construct filename
        filename = (
            row['carbon_source'] + '_' + 
            str(round(row['qS_carbon'], 2)) + '_' +
            row['nitrogen_source'] + '_' +
            str(round(row['qS_nitrogen'], 2))
        )
        
        # run FBA analysis on a copy of the model
        if row['growth_rate_exp'] == 0.25:
            result = mfa.run_FBA(
                model = model.copy(), 
                medium = mm,
                FVA = reactions_FVA
                )
            # save result from pandas data frame to hdd
            result['FBA'].to_csv('simulations/' + filename + '_FBA.csv')
            result['FVA'].to_csv('simulations/' + filename + '_FVA.csv')
        else:
            result = mfa.run_FBA(
                model = model.copy(), 
                medium = mm
                )
            # save result from pandas data frame to hdd
            result['FBA'].to_csv('simulations/' + filename + '_FBA.csv')


if __name__== "__main__" :
    main()
