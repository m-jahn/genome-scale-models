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

# set the path to the model's directory
wd = '/home/michael/Documents/SciLifeLab/Resources/Models/genome-scale-models/Ralstonia_eutropha/'
model = cobra.io.read_sbml_model(wd + "sbml/RehMBEL1391_sbml_L2V1.xml")
# alternatively load test model for comparison
#ecoli = cobra.test.create_test_model("textbook")
#salmonella = cobra.test.create_test_model()

# summary of the imported model
print('%i reactions' % len(model.reactions))
print('%i metabolites' % len(model.metabolites))
print('%i genes' % len(model.genes))


# The model contains two different sets of environmental variables
# that are important for simulation of growth:
# A: the model has exchange reactions that simply
#    model the extracellular availability (flux) towards the cell;
#    these have bounds = (0, 1000) when they 'produce' the metabolite
#    or e.g. (-1000, 0) when they 'consume' it
# B: model.medium is the mirror of the source exchange reaction
#    the medium has positive flux, i.e. produces a compound
#
# Bounds are not correctly set for reactions during import, only
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
model.exchanges.list_attr("bounds")
model.reactions.query("[A-Z]+t").list_attr("bounds")


# set a growth
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


# MODIFY/ADD REACTIONS -------------------------------------------------
#
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
model.reactions.PYK.bounds = (0, 1000.0)
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
HCO3E.lower_bound = -1000
HCO3E.upper_bound = 1000
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


# TESTING WITH FBA -----------------------------------------------------
#
# run FBA
# run FBA analysis
solution = model.optimize()

# print solution and status from linear programming solver
print([solution, "status: ", solution.status])

# quick summary of FBA analysis
print(model.summary())

# summary of energy balance
print(model.metabolites.atp_c.summary())

# summary of redox balance
print(model.metabolites.nadh_c.summary())


# EXPORT SBML MODEL ----------------------------------------------------
#
# save modified model as SBML Level 3, version 1
cobra.io.write_sbml_model(model, wd + "sbml/RehMBEL1391_sbml_L3V1.xml")
cobra.io.save_json_model(model, wd + "escher/RehMBEL1391_sbml_L3V1.json")