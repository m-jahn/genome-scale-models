# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                       ESSENTIALITY ANALYSIS                        +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
# This script has the purpose to
#   - import a genome scale model
#   - define different growth media
#   - run an essentiality analysis for genes and reactions
#
#
# LIBRARIES ------------------------------------------------------------

import numpy as np
import pandas as pd
import tabulate
import cobra
import cobra.test
import csv
import flux_analysis as mfa


# UTILITY FUNCTIONS ----------------------------------------------------

def get_minimal_medium():
    
    # Ralstonia eutropha minimal medium with the
    # following core components, not including carbon
    # or nitrogen source
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
    return(minimal_medium)


def get_lb_medium(model_reactions):
    
    # obtain LB medium from a curated Salmonella model
    lb_medium = cobra.test.create_test_model().media_compositions['LB']
    # dict comprehension to filter LB medium to include only the 
    # items present in model
    lb_medium = {i: lb_medium[i]*-1 for i in list(lb_medium.keys()) if i in model_reactions}
    # return as a dict
    return(lb_medium)


# function to run essentiality analysis on single model instance 
# and minimal medium
def run_essentiality(model, medium):
    
    model.medium = medium
    result = cobra.flux_analysis.find_essential_genes(model)
    result = [i.id for i in result]
    return(result)


# MAIN FUNCTION --------------------------------------------------------
# wrap all code into main function
def main():
    
    # IMPORT SBML MODEL ------------------------------------------------
    #
    model = cobra.io.read_sbml_model("sbml/RehMBEL1391_sbml_L3V1.xml")
    
    
    # RUN SINGLE ESSENTIALITY ANALYSIS FOR LB --------------------------
    #
    # create a new pandas df for results
    df = pd.DataFrame(
        [i.id for i in model.genes],
        columns = ['gene']
        )
    
    # test essentiality on LB complete medium
    print('...running essentiality analysis for: LB medium')
    result_lb = run_essentiality(
        model = model.copy(), 
        medium = get_lb_medium(model.reactions.list_attr('id'))
        )
    
    # assign result to new column in df
    df['LB_medium'] = [1 if i.id in result_lb else 0 for i in model.genes]
    print('...found ' + str(len(result_lb)) + ' essential genes')
    
    
    # RUN ESSENTIALITY ANALYSIS IN A LOOP ------------------------------
    #
    # to test essentiality on different growth substrates
    # that were tested experimentally
    #
    # first get predefined minimal medium
    minimal_medium = get_minimal_medium()
    
    # define different substrate species and their uptake rates
    qS_substrate = pd.DataFrame([
        mfa.qS_FRC(mu = 0.1),
        mfa.qS_SUC(mu = 0.1),
        mfa.qS_FOR(mu = 0.1)
        ], columns = ['qS_carbon'])
        
    qS_substrate['carbon_source'] = ([
        'EX_fru_e',
        'EX_succ_e',
        'EX_formate_e'
        ])
        
    qS_substrate['qS_nitrogen'] = (
        [10] * 3)
    
    qS_substrate['nitrogen_source'] = (
        ['EX_nh4_e'] * 3)
    
    
    # run several simulations in a loop
    for index, row in qS_substrate.iterrows():
        
        # add limiting substrate uptake rate to minimal medium
        mm = minimal_medium.copy()
        mm[row['carbon_source']] = row['qS_carbon']
        mm[row['nitrogen_source']] = row['qS_nitrogen']
        
        # construct name
        result_name = (
            row['carbon_source'] + '_' + 
            row['nitrogen_source']
            )
        
        # run essentiality analysis on a copy of the model
        print('...running essentiality analysis for: ' + result_name)
        result = run_essentiality(
            model = model.copy(), 
            medium = mm
            )
        
        # assign result to new column in df
        df[result_name] = [1 if i.id in result else 0 for i in model.genes]
        print('...found ' + str(len(result)) + ' essential genes')
    
    # save result from pandas data frame to hdd
    result_dir = 'simulations/essentiality/model_gene_essentiality.csv'
    df.to_csv(result_dir)
    print('...result dataframe saved to ' + result_dir)


if __name__== "__main__" :
    main()
