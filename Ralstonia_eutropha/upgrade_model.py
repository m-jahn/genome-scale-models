# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +    IMPORT, MODIFICATION, AND TESTING OF RALSTONIA EUTROPHA GEM     +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Script to import and modify an SBML (systems biology markup
# language) model which describea the metabolism of a bacterial cell.
# This script performs the following basic steps
# 
#  - import cobra and other required libraries
#  - import an SBML model
#  - check key data of the model
#  - change constraints/reactants of erroneous reactions
#  - add reactions that were missing
#  - run a basic flux balance analysis (FBA)
#
# All of these functions are based on the COBRApy documentation.
# To install cobrapy for python3, navigate to your favorite directory 
# and execute the following lines in terminal:
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
import csv
import re



# IMPORT SBML MODEL ----------------------------------------------------

def import_model(path):
    
    # read model from given path
    model = cobra.io.read_sbml_model(path)
    # alternatively load test model for comparison
    #ecoli = cobra.test.create_test_model("textbook")
    #salmonella = cobra.test.create_test_model()
    
    # summary of the imported model
    print(' ----- KEY NUMBERS ----- ')
    print('%i reactions' % len(model.reactions))
    print('%i metabolites' % len(model.metabolites))
    print('%i genes' % len(model.genes))
    
    # Fix exchnage reactions:
    # bounds are not correctly set for reactions during import, only
    # set to standard (-1000, 1000). Many reactions are then reversible
    # which is wrong. Set bounds based on info in 'notes' field, and
    # set all exchange reactions to c(0, 1000))
    for reaction in model.reactions:
        if reaction in model.exchanges:
            reaction.lower_bound = 0.0
        elif reaction.notes['REVERSIBILITY'] == 'true':
            reaction.lower_bound = -1000.0
        elif reaction.notes['REVERSIBILITY'] == 'false':
            reaction.lower_bound = 0.0
    
    # check that it worked
    print(' ----- BOUNDS FOR EXCHANGE AND TRANSPORT REACTIONS ----- ')
    print(model.exchanges.list_attr('bounds'))
    print(model.reactions.query('[A-Z]+t').list_attr('bounds'))
    
    # return imported model
    return model


# TESTING FOR ENERGY GENERATING CYCLES ---------------------------------
#
# Energy generating cycles can generate ATP, NADPH or other energy/reducant 
# 'currency' of the cell out of thin air (see Fritzemeier et a., PLOS Comp Bio, 2017). 
# Such cycles represent factual errors of the model, and they 
# can occur when reactions are not correctly mass or charge balanced. 
# For example, one reaction converts a metabolite into another and a
# second reaction converts it back but this time generates ATP. Such cycles
# can even appear when reactions are correctly balanced, but one of the two 
# reactions would not normally occur because of thermodynamic restrictions 
# (far from equilibrium,) and therefore not reversible. FBA does not take
# thermodynamics into account but we can constrain reversibility as a fix.

def test_EGC(model):
    
    # According to the workflow of Fritzemeier et al., we add dissipation reactions
    # e.g. for ATP, NADH and so on.
    diss_reactions = [s + '_diss' for s in ['ATP', 'GTP', 'CTP', 'UTP', 'ITP', 'NADH', 'NADPH', 'FADH2', 'UQ', 'ACCOA', 'GLU', 'PRO']]
    model.add_reactions([cobra.Reaction(r) for r in diss_reactions])
    
    model.reactions.ATP_diss.build_reaction_from_string('atp_c + h2o_c --> adp_c + h_c + pi_c')
    model.reactions.CTP_diss.build_reaction_from_string('ctp_c + h2o_c --> cdp_c + h_c + pi_c')
    model.reactions.GTP_diss.build_reaction_from_string('gtp_c + h2o_c --> gdp_c + h_c + pi_c')
    model.reactions.UTP_diss.build_reaction_from_string('utp_c + h2o_c --> udp_c + h_c + pi_c')
    model.reactions.ITP_diss.build_reaction_from_string('itp_c + h2o_c --> idp_c + h_c + pi_c')
    model.reactions.NADH_diss.build_reaction_from_string('nadh_c --> nad_c + h_c')
    model.reactions.NADPH_diss.build_reaction_from_string('nadph_c --> nadp_c + h_c')
    model.reactions.FADH2_diss.build_reaction_from_string('fadh2_c --> fad_c + 2.0 h_c')
    model.reactions.UQ_diss.build_reaction_from_string('uqh2_c --> uq_c + 2.0 h_c')
    model.reactions.ACCOA_diss.build_reaction_from_string('accoa_c  + h2o_c --> ac_c + coa_c + h_c')
    model.reactions.GLU_diss.build_reaction_from_string('glu_c  + h2o_c --> akg_c + nh4_c + 2.0 h_c')
    model.reactions.PRO_diss.build_reaction_from_string('h_e --> h_c')
    
    
    # Next we loop through the dissipation reactions
    print(' ----- SUMMARY OF ENERGY GENERATING CYCLES ----- ')
    for r in diss_reactions:
        
        # first set the ojective to max dissipation reaction
        model.objective = {model.reactions.get_by_id(r): 1}
        
        # run FBA analysis, all solutions should evaluate to zero
        solution = model.optimize()
        print(str(solution) + ', dissipation reaction: ' + str(r) + ', status: ' + str(solution.status))
    
    # clean up by removing all dissipation reactions
    model.remove_reactions(diss_reactions)


# MODIFY/ADD REACTIONS -------------------------------------------------
#
# several reactions were identified as erroneous and this function 
# corrects those. It also adds reactions that were missing.

def modify_reactions(model):
    
    # There's an artificial NADH generating cycle around the metabolite
    # 1-pyrroline-5-carboxylate dehydrogenase involving 3 reactions,
    # P5CD4 --> PROD4/P5CD5 --> PTO4H --> P5CD4
    # The cycle generates 2 NADH and 2 H+ per turn and carries high flux
    # According to this BiGG reaction: http://bigg.ucsd.edu/models/iIT341/reactions/4HGLSD,
    # the NADH/NAD are cofactor stoichiometry is reversed
    # (should be: h2o_c + nad_c + 4hglusa_c ⇌ 2.0 h_c + nadh_c + e4hglu_c)
    model.reactions.P5CD5.build_reaction_from_string('e4hglu_c + nadh_c + 2.0 h_c <=> 4hglusa_c + h2o_c + nad_c')
    
    # Another problem with this reaction is, it is duplicated as PROD4 
    # (including the same gene annotation --> remove it)
    model.reactions.PROD4.delete()
    
    # A reaction that seems to contain an error is NADHDH
    # a NADH dehydrogenase that converts UQ spontaneously to UQH2,
    # a reaction that requires NADH canonically (query BiGG reaction: NADHDH)
    model.reactions.NADHDH.build_reaction_from_string('3.0 h_c + nadh_c + uq_c --> nad_c + 2.0 h_e + uqh2_c')
    
    # The TCA cycle reaction succinyl-CoA-synthetase SUCOAS is importantly not set to 'reversible'
    # as it's supposed to be, and therefore constrained to wrong direction
    # regarding canonical flow of TCA
    model.reactions.SUCOAS.bounds = (-1000.0, 1000.0)
    
    # The reaction MICITL is the last step of the methyl-citrate cycle, an
    # alternative route through TCA from oaa + prop-coa --> succ + pyr.
    # It carries artificially high flux, so flux of the final reaction was constrained.
    model.reactions.MICITL.bounds = (0.0, 0.0)
    
    # PYK is allowed in the model to go in reverse direction (pyr + atp --> pep + adp)
    # but this is highly unlikely under physiological conditions (see e.g. wikipedia, or BiGG database). 
    # Standard E. coli models also exclude the reverse reaction.
    model.reactions.PYK.bounds = (0.0, 1000.0)
    # Several alternative reactions to PYK (PYK1, PYK2, PYK3) that caryy most likely very 
    # little or no flux in R. eutropha, were silenced.
    model.reactions.PYK1.bounds = (0.0, 0.0)
    model.reactions.PYK2.bounds = (0.0, 0.0)
    model.reactions.PYK3.bounds = (0.0, 0.0)
    
    # Pyruvate carboxylase (PYC) should only run in direction 
    # from pyr --> oaa, but not reverse (see E. coli reference models in BiGG).
    model.reactions.PYC.bounds = (0.0, 1000.0)
    # PEP carboxylase PPC has correct bounds but one H+ reactant too much.
    # The reaction was corrected.
    model.reactions.PPC.build_reaction_from_string('co2_c + h2o_c + pep_c --> h_c + oaa_c + pi_c')
    
    # Two different metabolites, asp_c and aspsa_c, are labelled with the 
    # name of aspartate in the model and take part in different reactions.
    # However aspsa_c is in reality L-Aspartate 4-semialdehyde (source: BiGG)
    # this was renamed in the model. The reactions are correct.
    model.metabolites.aspsa_c.name = 'L-Aspartate 4-semialdehyde'
    
    # Some reactions (e.g. CABPS, carbamoylphosphate synthase) require hco3
    # instead of co2 as substrate for phophorylation. However, a co2 <=> hco3
    # equilibration reaction is missing (see BiGG reaction HCO3E).
    # The reaction was added.
    HCO3E = cobra.Reaction('HCO3E')
    HCO3E.name = 'HCO3 equilibration reaction'
    HCO3E.lower_bound = -1000.0
    HCO3E.upper_bound = 1000.0
    model.add_reactions([HCO3E])
    model.reactions.HCO3E.build_reaction_from_string('co2_c + h2o_c <=> h_c + hco3_c')
    
    # The primary means of fructose uptake seems to be the ABC transporter
    # (ATP dependent import). A second ATP-dependent reaction, fructokinase, 
    # then phosphorylates fru --> f6p. 
    # It is not clear if PEP-PTS dependent fructose uptake and 
    # phosphorylation exists in R. eutropha. Therefore the PEP-PTS reaction
    # was silenced (more details, see Kaddor & Steinbuechel, 2011).
    model.reactions.FRUpts2.bounds = (0.0, 0.0)
    # add gene-reaction-rules for fruABC
    model.reactions.FRUabc.gene_reaction_rule = '( H16_B1498 and H16_B1499 and H16_B1500 )'
    
    # The original model only contains a lumped reaction for the CBB cycle.
    # In order to include a working CBB cycle, two reactions need to be added,
    # 1) Phosphoribulokinase (cbbP2, H16_B1389; cbbPp, PHG421) catalyzing 
    # phosphorylation of Ribulose-5-phosphate: atp_c + rl5p_c --> adp_c + h_c + rb15bp_c
    # 2) Ribulose-1,5-bisphosphate carboxylase (cbbS2, H16_B1394; 
    # cbbL2, H16_B1395; cbbSp, PHG426, cbbLp, PHG427)
    # catalyzing the addition of CO2: co2_c + h2o_c + rb15bp_c --> 2.0 h_c + 2.0 3pg_c
    # The metabolite Ribulose-1,5-bisphosphate needs to be added
    model.add_metabolites(
        cobra.Metabolite(
            id = 'rb15bp_c',
            name = 'Ribulose-1,5-bisphosphate',
            compartment = 'c'))
    
    # Reactions PRUK and RBPC (Rubisco enzyme) are added
    PRUK = cobra.Reaction('PRUK')
    PRUK.name = 'Phosphoribulokinase'
    RBPC = cobra.Reaction('RBPC')
    RBPC.name = 'Ribulose-bisphosphate carboxylase'
    model.add_reactions([PRUK, RBPC])
    model.reactions.PRUK.gene_reaction_rule = 'H16_B1389 or PHG421'
    model.reactions.RBPC.gene_reaction_rule = '( H16_B1394 and H16_B1395 ) or ( PHG426 and PHG427 )'
    
    # define reactions from string
    model.reactions.PRUK.build_reaction_from_string('atp_c + rl5p_c --> adp_c + h_c + rb15bp_c')
    model.reactions.RBPC.build_reaction_from_string('co2_c + h2o_c + rb15bp_c --> 2.0 3pg_c + 2.0 h_c')
    
    # silence original lumped reaction for CBB cycle
    model.reactions.CBBCYC.bounds = (0.0, 0.0)


# TESTING GROWTH WITH FBA ------------------------------------------
#
# run one simple FBA analysis with the model, simulating 
# growth in fructose

def test_FBA(model):
    
    # set a growth medium
    model.medium = {
        'EX_mg2_e': 10.0,
        'EX_pi_e': 100.0,
        'EX_cobalt2_e': 10.0,
        'EX_cl_e': 10.0,
        'EX_k_e': 10.0,
        'EX_fe3_e': 10.0,
        'EX_so4_e': 10.0,
        'EX_na_e': 10.0,
        'EX_o2_e': 18.5,
        'EX_mobd_e': 10.0,
        'EX_h2o_e': 1000.0,
        'EX_h_e': 100.0,
        'EX_fru_e': 3.0,
        'EX_nh4_e': 10.0
        }
    
    # set objective function
    model.objective = {model.reactions.Biomass: 1}
    # run FBA analysis
    solution = model.optimize()
    
    # print solution and status from linear programming solver
    print(' ----- SUMMARY OF FBA TEST RUN ----- ')
    print(str(solution) + ', status: ' + str(solution.status))
    
    # quick summary of FBA analysis
    print(model.summary())


# EXECUTE FUNCTIONS ----------------------------------------------------

def main():
    
    # import model
    model = import_model(path = 'sbml/RehMBEL1391_sbml_L2V1.xml')
    
    # test energy generating cycles
    test_EGC(model.copy())
    
    # modify reactions
    modify_reactions(model)
    
    # test energy generating cycles again
    test_EGC(model.copy())
    
    # test run with FBA
    test_FBA(model)
    
    # save modified model as SBML Level 3, version 1
    cobra.io.write_sbml_model(model, 'sbml/RehMBEL1391_sbml_L3V1.xml')
    cobra.io.save_json_model(model, 'escher/RehMBEL1391_sbml_L3V1.json')

if __name__== '__main__' :
    main()
