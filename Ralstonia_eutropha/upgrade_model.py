# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +    IMPORT, MODIFICATION, AND TESTING OF RALSTONIA EUTROPHA GEM     +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Functions to import and modify an SBML (systems biology markup
# language) model which describes the metabolism of a bacterial cell.
# Functions perform the following steps: 
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
import cobra
import csv
import re
from cobrapy_bigg_client import client
from bioservices.uniprot import UniProt
from bioservices.kegg import KEGG
from os import path


# IMPORT SBML MODEL ----------------------------------------------------

def import_model(path):
    
    # read model from given path
    model = cobra.io.read_sbml_model(path)
    
    # summary of the imported model
    print(' ----- MODEL KEY NUMBERS ----- ')
    print('%i reactions' % len(model.reactions))
    print('%i metabolites' % len(model.metabolites))
    print('%i genes' % len(model.genes))
    
    # Fix exchange reactions:
    # bounds are not correctly set for reactions during import, only
    # set to standard (-1000, 1000). Many reactions are then reversible
    # which is wrong. Set bounds based on info in 'notes' field, and
    # set all exchange reactions to c(0, 1000))
    print('...setting bounds on exchange reactions')
    for reaction in model.reactions:
        if reaction in model.exchanges:
            reaction.lower_bound = 0.0
        elif reaction.notes['REVERSIBILITY'] == 'true':
            reaction.lower_bound = -1000.0
        elif reaction.notes['REVERSIBILITY'] == 'false':
            reaction.lower_bound = 0.0
    
    # check that it worked
    print(' ----- BOUNDS FOR EXCHANGE AND TRANSPORT REACTIONS ----- ')
    print('first 10 exchange reactions: \n' + str(model.exchanges.list_attr('bounds')[0:9]))
    print('first 10 metabolic reactions: \n' + str(model.reactions.query('[A-Z]+t').list_attr('bounds')[0:9]))
    
    # return imported model
    return model


# TESTING FOR ENERGY GENERATING CYCLES ---------------------------------
#
# Energy generating cycles can generate ATP, NADPH or other energy/reducant 
# 'currency' of the cell out of thin air (see Fritzemeier et al., PLOS Comp Bio, 2017). 
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
    model.reactions.GLU_diss.build_reaction_from_string('glu__L_c  + h2o_c --> akg_c + nh4_c + 2.0 h_c')
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
    
    # 1st STEP: remove duplicated reactions
    # -------------------------------------
    # update model IDs, names, and gene_reaction_rules the reaction that
    # is kept from a pair of duplicates
    model.reactions.TYRTA1.id = 'TYRTA'
    model.reactions.TYRTA.gene_reaction_rule = model.reactions.TYRTA2.gene_reaction_rule
    model.reactions.AMDS1.name = '2-Phenylacetamide Amidohydrolase'
    model.reactions.BZACCOAT1.id = 'BZACCOAT'
    model.reactions.BZACCOAT.gene_reaction_rule = (
        model.reactions.BZACCOAT.gene_reaction_rule + ' or ' + 
        model.reactions.BZACCOAT2.gene_reaction_rule)
    model.reactions.get_by_id('4NPHP1').id = '4NPHPP'
    model.reactions.get_by_id('4NPHPP').name = '4-Nitrophenylphosphate phosphatase'
    model.reactions.get_by_id('4NPHPP').gene_reaction_rule = (
        model.reactions.get_by_id('4NPHPP').gene_reaction_rule + ' or ' + 
        model.reactions.get_by_id('4NPHP2').gene_reaction_rule)
    model.reactions.BZACCOAT.gene_reaction_rule = (
        model.reactions.HIBD.gene_reaction_rule + ' or ' + 
        model.reactions.HACOAD3.gene_reaction_rule)
    model.reactions.G1PTT1.id = 'G1PTT'
    model.reactions.G1PTT.gene_reaction_rule = (
        model.reactions.G1PTT.gene_reaction_rule + ' or ' + 
        model.reactions.G1PTT2.gene_reaction_rule)
    model.reactions.NITRR.id = 'NTRIR2x'
    model.reactions.NTRIR2x.gene_reaction_rule = model.reactions.NITRT.gene_reaction_rule
    model.reactions.SPMS1.id = 'SPMS'
    model.reactions.P5CD2.name = 'L-Glutamate 5-semialdehyde dehydrogenase'
    model.reactions.P5CD3.name = 'Trans-4-Hydroxy-L-proline dehydrogenase'
    model.reactions.P5CD4.name = 'L-1-Pyrroline-3-hydroxy-5-carboxylate dehydrogenase'
    model.reactions.NO3RUQ1.id = 'NO3R1'
    model.reactions.UDPG4E.name = 'UDP-glucose 4-epimerase'
    model.reactions.DADNK.id = 'DADK'
    model.reactions.DADK.subtract_metabolites({'h_c': -1})
    model.reactions.UMPK.gene_reaction_rule = (
        model.reactions.UMPK.gene_reaction_rule + ' or ' + 
        model.reactions.URIDK2.gene_reaction_rule)
    model.reactions.SERD.id = 'SERD_L'
    model.reactions.SERD_L.gene_reaction_rule = (
        model.reactions.SERD_L.gene_reaction_rule + ' or ' + 
        model.reactions.SERDHT1.gene_reaction_rule)
    model.reactions.OAHSL2.id = 'AHSERL4'
    model.reactions.AHSERL4.name = 'Cysteine synthase Thiosulfate'
    model.reactions.AHSERL4.gene_reaction_rule = (
        model.reactions.AHSERL4.gene_reaction_rule + ' or ' + 
        model.reactions.CYSST3.gene_reaction_rule)
    model.reactions.GLUDH1.id = 'GLUDy'
    model.reactions.GLUDy.name = 'Glutamate dehydrogenase (NADP)'
    model.reactions.GLUDy.gene_reaction_rule = (
        model.reactions.GLUDy.gene_reaction_rule + ' or ' + 
        model.reactions.GLUDH2.gene_reaction_rule)
    model.reactions.GLUDH3.id = 'GLUDxi'
    model.reactions.GLUDxi.name = 'Glutamate dehydrogenase (NAD)'
    model.reactions.GLUDxi.gene_reaction_rule = (
        model.reactions.GLUDxi.gene_reaction_rule + ' or ' + 
        model.reactions.GLUDH4.gene_reaction_rule)
    model.reactions.CYTTS6.id = 'SLCYSS'
    model.reactions.SLCYSS.name = 'O-acetyl-L-serine sulfhydrylase'
    model.reactions.SLCYSS.gene_reaction_rule = (
        model.reactions.SLCYSS.gene_reaction_rule + ' or ' + 
        model.reactions.CYSST2.gene_reaction_rule)
    model.reactions.PPCSYN1.gene_reaction_rule = (
        model.reactions.PPCSYN1.gene_reaction_rule + ' or ' + 
        model.reactions.ACCSYN1.gene_reaction_rule)
    model.reactions.PPCSYN2.gene_reaction_rule = (
        model.reactions.PPCSYN2.gene_reaction_rule + ' or ' + 
        model.reactions.ACCSYN2.gene_reaction_rule)
    model.reactions.ASPOX1.id = 'ASPO1'
    model.reactions.ASPO1.gene_reaction_rule = (
        model.reactions.ASPO1.gene_reaction_rule + ' or ' + 
        model.reactions.LAAO1.gene_reaction_rule)
    model.reactions.ALHD3.id = 'ALDD19xr'
    model.reactions.ACCOACB.id = 'ACCOAC'
    model.reactions.ACCOAC.gene_reaction_rule = ('( ' +
        model.reactions.ACCOAC.gene_reaction_rule + ' ) or ' + 
        model.reactions.BITCB.gene_reaction_rule)
    model.reactions.ASNN.gene_reaction_rule = (
        model.reactions.ASNN.gene_reaction_rule + ' or ' + 
        model.reactions.GLNAS.gene_reaction_rule)
    
    # remove all duplicated reactions
    duplicated_reactions = (
        ['TYRTA2', 'ASPAM5', 'AMDS4', 'BZACCOAT2', '4NPHP2', 'HACOAD3', 
        'G1PTT2', 'NA1t2', 'NITRT', 'SPMS2', 'PROD3', 'PROD4', 'NO3RUQ2',
        'UDPG4E2', 'ADNK2', 'URIDK2', 'SERDHT1', 'CYSST3', 'ASPAM6',
        'GLUDH2', 'GLUDH4', 'CYSST2', 'ACFM4', 'ACCSYN1', 'ACCSYN2',
        'LAAO1', 'PHEALDD', 'BITCB', 'GLNAS'])
    
    model.remove_reactions(duplicated_reactions)
    for dupl in duplicated_reactions:
        print('...removed duplicated reaction ' + dupl)
    
    
    # 2nd STEP: add names to unnamed reactions
    # ----------------------------------------
    # overiew about reactions
    model_reactions = pd.DataFrame({
        'id': model.reactions.list_attr('id'),
        'name': model.reactions.list_attr('name'),
        'reaction': model.reactions.list_attr('reaction'),
        'gene': model.reactions.list_attr('gene_reaction_rule')
        })
    
    count = 0
    reactions_unnamed = model_reactions[model_reactions['name'] == '']['id'].to_list()
    for rea_un in reactions_unnamed:
        current_rea = model.reactions.get_by_id(rea_un)
        if current_rea.compartments == {'c', 'e'}:
            reactants = [len(i.name) for i in current_rea.reactants]
            rea_name = current_rea.reactants[reactants.index(max(reactants))].name + ' transport'
            current_rea.name = rea_name
            count = count + 1
    
    
    # 3rd STEP: correct known errors in reactions
    # -------------------------------------------
    # There's an artificial NADH generating cycle around the metabolite
    # 1-pyrroline-5-carboxylate dehydrogenase involving 3 reactions,
    # P5CD4 --> PROD4/P5CD5 --> PTO4H --> P5CD4
    # The cycle generates 2 NADH and 2 H+ per turn and carries high flux
    # According to this BiGG reaction: http://bigg.ucsd.edu/models/iIT341/reactions/4HGLSD,
    # the NADH/NAD are cofactor stoichiometry is reversed
    # (should be: h2o_c + nad_c + 4hglusa_c <-> 2.0 h_c + nadh_c + e4hglu_c)
    model.reactions.P5CD5.build_reaction_from_string('e4hglu_c + nadh_c + 2.0 h_c <=> 4hglusa_c + h2o_c + nad_c')
    model.reactions.P5CD5.name = 'L-4-Hydroxyglutamate semialdehyde dehydrogenase'
    print('...corrected cofactor usage for reaction P5CD5')
        
    # A reaction that seems to contain an error is NADHDH
    # a NADH dehydrogenase that converts UQ spontaneously to UQH2,
    # a reaction that requires NADH canonically (query BiGG reaction: NADHDH)
    model.reactions.NADHDH.build_reaction_from_string('3.0 h_c + nadh_c + uq_c --> nad_c + 2.0 h_e + uqh2_c')
    print('...corrected cofactor usage for reaction NADHDH')
    
    # The TCA cycle reaction succinyl-CoA-synthetase SUCOAS is importantly not set to 'reversible'
    # as it's supposed to be, and therefore constrained to wrong direction
    # regarding canonical flow of TCA
    model.reactions.SUCOAS.bounds = (-1000.0, 1000.0)
    print('...set reversibility for reaction SUCOAS')
    
    # The reaction MICITL is the last step of the methyl-citrate cycle, an
    # alternative route through TCA from oaa + prop-coa --> succ + pyr.
    # It carries artificially high flux, so flux of the final reaction was constrained.
    model.reactions.MICITL.bounds = (0.0, 0.0)
    print('...silenced reaction MICITL, the last step of the methyl citrate cycle')
    
    # PYK is allowed in the model to go in reverse direction (pyr + atp --> pep + adp)
    # but this is highly unlikely under physiological conditions (see e.g. wikipedia, or BiGG database). 
    # Standard E. coli models also exclude the reverse reaction.
    model.reactions.PYK.bounds = (0.0, 1000.0)
    # Several alternative reactions to PYK (PYK1, PYK2, PYK3) that caryy most likely very 
    # little or no flux in R. eutropha, were silenced.
    model.reactions.PYK1.bounds = (0.0, 0.0)
    model.reactions.PYK2.bounds = (0.0, 0.0)
    model.reactions.PYK3.bounds = (0.0, 0.0)
    print('...set irreversibility for reaction PYK')
    
    # Pyruvate carboxylase (PYC) should only run in direction 
    # from pyr --> oaa, but not reverse (see E. coli reference models in BiGG).
    model.reactions.PYC.bounds = (0.0, 1000.0)
    print('...set irreversibility for reaction PYC')
    
    # PEP carboxylase PPC has correct bounds but one H+ reactant too much.
    # The reaction was corrected.
    model.reactions.PPC.build_reaction_from_string('co2_c + h2o_c + pep_c --> h_c + oaa_c + pi_c')
    print('...corrected products for reaction PPC')
    
    # The biomass equation incorrectly contains pyridoxine as the required
    # cofactor, but the canonical metaboliyte is pyridoxal-5-phosphate (pydx5p)
    # Correcting this error prevents an infeasible cycle of reactions
    # where the same enzyme shows different reaction directionalities for
    # the same metabolites (only difference being phosphorylated or not): 
    # pdx5p --> pydx5p --> pyam5p ---|  (de-phosph. to pydam)
    # pydxn <-- pydx   <-- pydam  <--|
    # Removing pydxn from cofactors and adding pydx5p 
    # (MW = 247.142 g/mol, 0.111 g/247.142 g/mol = 0.449 mmol/g)
    model.reactions.Cofactors_and_vitamins.add_metabolites({
        model.metabolites.pydx5p_c: -0.449})
    model.reactions.Cofactors_and_vitamins.subtract_metabolites({
        model.metabolites.pydxn_c: -0.656})
    print('...corrected Cofactor reaction by replacing pyridoxine with pyridoxal-5-phophate')
    
    # Therefore correct two reactions in pyridoxal phosphate metabolism
    # catalyzed by the same enzyme (H16_A2802). These reactions are
    # the same for the phosphorylated variant of pyridoxal metabolites,
    # pdx5p --> pydx5p <-- pyam5p.
    # pyridoxamine oxidase should only run in forward direction 
    # according to http://modelseed.org/biochem/reactions/rxn01252
    # --> change name, id, and bounds
    model.reactions.PYR5OXX.id = 'PYDXNO'
    model.reactions.PYDXNO.name = 'Pyridoxine oxidase'
    model.reactions.PYDXNO.bounds = (0.0, 1000.0)
    print('...set irreversibility for reaction PYDXNO')
    
    # pyridoxal oxidase should only run in forward direction
    # according to http://modelseed.org/biochem/reactions/rxn01251
    # --> change name, id, and bounds
    model.reactions.PYR5OXM.id = 'PYDXO_1'
    model.reactions.PYDXO_1.name = 'Pyridoxamine:oxygen oxidoreductase (deaminating)'
    model.reactions.PYDXO_1.bounds = (0.0, 1000.0)
    print('...set irreversibility for reaction PYDXO_1')
    
    # Two different metabolites, asp_c and aspsa_c, are labelled with the 
    # name of aspartate in the model and take part in different reactions.
    # However aspsa_c is in reality L-Aspartate 4-semialdehyde (source: BiGG)
    # this was renamed in the model. The reactions are correct.
    model.metabolites.aspsa_c.name = 'L-Aspartate 4-semialdehyde'
    print('...corrected name for aspsa_c to L-Aspartate 4-semialdehyde')
    
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
    print('...added a HCO3 equilibration reaction')
    
    # The primary means of fructose uptake seems to be the ABC transporter
    # (ATP dependent import). A second ATP-dependent reaction, fructokinase, 
    # then phosphorylates fru --> f6p. 
    # It is not clear if PEP-PTS dependent fructose uptake and 
    # phosphorylation exists in R. eutropha. Therefore the PEP-PTS reaction
    # was silenced (more details, see Kaddor & Steinbuechel, 2011).
    model.reactions.FRUpts2.bounds = (0.0, 0.0)
    print('...inactivated FRUpts2 (fructose PEP-PTS) reaction')
    
    # add gene-reaction-rules for fruABC
    model.reactions.FRUabc.gene_reaction_rule = '( H16_B1498 and H16_B1499 and H16_B1500 )'
    print('...added gene association for fructose ABC transporter')
    
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
    model.reactions.PRUK.gene_reaction_rule = '( H16_B1389 ) or ( PHG421 )'
    model.reactions.RBPC.gene_reaction_rule = '( H16_B1394 and H16_B1395 ) or ( PHG426 and PHG427 )'
    
    # define reactions from string
    model.reactions.PRUK.build_reaction_from_string('atp_c + rl5p_c --> adp_c + h_c + rb15bp_c')
    model.reactions.RBPC.build_reaction_from_string('co2_c + h2o_c + rb15bp_c --> 2.0 3pg_c + 2.0 h_c')
    
    # silence original lumped reaction for CBB cycle
    model.reactions.CBBCYC.bounds = (0.0, 0.0)
    print('...added PRUK and RBPC (Rubisco) reactions instead of lumped reaction')
    
    # final reporting
    print(' ----- SUMMARY OF MODIFIED REACTIONS ----- ')
    print('removed duplicated reactions: ' + str(len(duplicated_reactions)))
    print('added names to reactions: ' + str(count))
    print('corrected erroneous reactions: 15')


# ADDING ANNOTATIONS FROM BIGG -----------------------------------------
# 
# helper function to convert Bigg keys to memote-compatible ones
def change_keys(input_str):
    if input_str == 'InChI Key': return 'inchikey'
    elif input_str == 'RHEA': return 'rhea'
    elif input_str == 'Reactome Compound': return 'reactome'
    elif input_str == 'Reactome Reaction': return 'reactome'
    elif input_str == 'CHEBI': return 'chebi'
    elif input_str == 'SEED Compound': return 'seed.compound'
    elif input_str == 'SEED Reaction': return 'seed.reaction'
    elif input_str == 'Human Metabolome Database': return 'hmdb'
    elif input_str == 'MetaNetX (MNX) Chemical': return 'metanetx.chemical'
    elif input_str == 'MetaNetX (MNX) Equation': return 'metanetx.reaction'
    elif input_str == 'KEGG Compound': return 'kegg.compound'
    elif input_str == 'KEGG Reaction': return 'kegg.reaction'
    elif input_str == 'BioCyc': return 'biocyc'
    elif input_str == 'EC Number': return 'ec-code'
    else: return input_str


# helper function to rename reactome IDs
def change_reactome_ID(id):
    if re.match('R-[A-Z]{3}-', id) == None:
        id = 'R-ALL-' + id
    return(id)


# function to 1) import existing Bigg annotation, and 2) add missing 
# annotation from bigg for a list of metabolites or reactions, 
# and 3) return items in a cobrapy model
def get_bigg_annotation(
    ref_path = 'data/reaction_reference.json',
    item_list = [],
    item_type = 'reaction'):
    
    
    # construct generic reference names based on input params
    ref_name = item_type + '_reference'
    
    # check if already an annotation file exists
    # if yes load, if no, create new one
    if path.exists(ref_path):
        reference = cobra.io.load_json_model(ref_path)
    else:
        reference = cobra.Model(ref_name)
    
    # determine which items are missing in reference
    ref_list = getattr(reference, item_type + 's').list_attr('id')
    item_list = [i for i in item_list if i not in ref_list]
    
    
    # download annotation for unannotated metabolites/reactions
    # and add to reference 'model'
    if len(item_list):
        
        for item in item_list:
            
            # download annotation from Bigg
            if item_type == 'reaction':
                ref_item = client.get_reaction(item)
            else: 
                ref_item = client.get_metabolite(item)[0]
            
            # process annotation
            new_annot = ref_item.annotation['database_links']
            new_annot = {change_keys(key): [i['id'] for i in new_annot[key]] for key in new_annot.keys()}
            # add the bigg identifier itself also to annotation slot
            new_annot['bigg.' + item_type] = ref_item.id
            # add a prefix to reactome to match standard regex pattern
            if 'reactome' in new_annot.keys():
                new_annot['reactome'] = [change_reactome_ID(i) for i in new_annot['reactome']]
            ref_item.annotation = new_annot
            getattr(reference, item_type + 's').add(ref_item)
            if item_type == 'reaction':
                reference.add_metabolites(ref_item.metabolites)
            reference.name = ref_name
            print('...downloaded annotation for ' + item_type + ': ' + ref_item.id)
        
        # export reference with added items as json file
        cobra.io.save_json_model(reference, ref_path)
        print('...exported file "' + ref_path + '" with ' + str(len(item_list)) + ' new ' + item_type + 's')
    
    return(reference)


# function to update metabolites or reactions with annotation
# obtained from Bigg
def update_bigg_annotation(
    model,
    item_type = 'reaction'):
    
    # fetch global list of items from BiGG
    if item_type == 'reaction':
        bigg = client.list_reactions()
    else:
        bigg = client.list_metabolites()
    
    # create pairs of metabolite trivial name and id from ref
    bigg = pd.DataFrame({
        'id': bigg.list_attr("id"),
        'name': bigg.list_attr("name")
        })
    
    # update IDs from old to new standard
    count_id = 0
    model_item_type = getattr(model, item_type + 's')
    for item in model_item_type:
        
        # capitalize first letter for reaction names
        if item_type == 'reaction':
            item.name = item.name[0:1].capitalize() + item.name[1:]
        hit = bigg['name'] == item.name
        if hit.any():
            if (sum(hit) != 1) and (item_type == 'metabolite'):
                print('...' + item_type + ' matches more than one reference, taking 1st hit')
            new_id = bigg.loc[hit, 'id'].to_list()[0]
            if item_type == 'metabolite':
                new_id = new_id + '_' + item.compartment
            # if A) new ID is not the old ID and 
            #    B) the new ID doesnt exist yet and 
            #    C) the item name is unique, otherwise more than 1 item will get same ID
            cond = new_id not in model_item_type.list_attr("id")
            #cond2 = sum([item.name == i  for i in model_item_type.list_attr("name")]) == 1
            #if (item.id != new_id) and cond1 and not cond2:
            #    print(item.id + ' | ' + new_id + ' | ' + item.name)
            #    print(sum([item.name == i  for i in model_item_type.list_attr("name")]))
            if (item.id != new_id) and cond:
                print('...converted old id: ' + item.id + ' to new id: ' + new_id)
                item.id = new_id
                count_id = count_id + 1
                # also update name of exchange reaction for respective 
                # metabolite if one exists
                if item_type == 'metabolite':
                    ex_rea = 'EX_' + item.id
                    if ex_rea in model.reactions.list_attr('id'):
                        model.reactions.get_by_id(ex_rea).id = 'EX_' + new_id
    
        
    # determine which metabolites are available in bigg
    if item_type == 'reaction':
        model_items = model.reactions.list_attr('id')
    else:
        model_items = [re.sub("_[cep]$", "", i) for i in model.metabolites.list_attr('id')]
    
    # determine intersection between model and bigg and remove duplicates
    bigg_available = [i for i in model_items if i in bigg['id'].to_list()]
    bigg_available = list(set(bigg_available))
    
        
    # import existing Bigg annotation from /data,
    # optionally downloading missing items from Bigg db
    reference = get_bigg_annotation(
        ref_path = 'data/' + item_type + '_reference.json',
        item_list = bigg_available, 
        item_type = item_type)
    
    
    # update annotation for metabolites/reactions
    count_an = 0
    for item in model_item_type:
        # remove compartment tag c, e, or p for metabolites
        if item_type == 'metabolite':
            search_id = re.sub("_[cep]$", "", item.id)
        else: 
            search_id = item.id
        if search_id in bigg_available:
            ref_item = getattr(reference, item_type + 's').get_by_id(search_id)
            item.annotation = ref_item.annotation
            if item_type == 'metabolite':
                item.charge = ref_item.charge
                item.formula = ref_item.formula
                item.notes = ref_item.notes
            count_an = count_an + 1
            print('...updated ' + item_type + 'annotation for ' + item.id)
    
    
    # final reporting
    print(' ----- SUMMARY OF ' + item_type.upper() + ' CHANGES ----- ')
    print('updated ' + item_type + ' IDs: ' + str(count_id))
    print('updated ' + item_type + ' annotation: ' + str(count_an))


# ADDING GENE ANNOTATIONS FROM UNIPROT ---------------------------------
# 
# for this purpose we use the very comprehensive bioservices
# package for python that has connectivity to all possible databases
def update_gene_annotation(model):
    
    # open client connection
    uni = UniProt(verbose = False)
    kegg = KEGG(verbose = False)
    
    # fetch uniprot annotation table for all genes as pandas df
    gene_list = model.genes.list_attr("id")
    print('...downloading annotation for ' + str(len(gene_list)) + ' genes from uniprot.org')
    df = uni.get_df(gene_list, limit = len(gene_list))
    print('...downloading annotation for ' + str(len(gene_list)) + ' genes from kegg.jp')
    dict_ncbi = {e: kegg.parse(kegg.get("reh:" + str(e)))['DBLINKS']['NCBI-ProteinID'] for e in gene_list}
    # remove entries that are not according to Reh standard
    df = df[[len(i) == 6 for i in df['Entry']]]
    
    # some rows are merged entries for 2 KEGG IDs
    # in this case we simply duplicate the rows and split entries
    kegg_id = 'Gene names  (ordered locus )'
    df_duplicated = list(filter(lambda x: len(str(x)) > 9, df[kegg_id]))
    df_duplicated = list(set(df_duplicated))
    print('...splitting duplicated entries: ' + str(df_duplicated))
    # duplicate rows for entries if necessary and split IDs
    # this currently works only if 2 IDs are merged, not more
    for i in df_duplicated:
        if sum(df[kegg_id] == i) == 1:
            df = df.append(df[(df[kegg_id] == i)])
        df.loc[df[kegg_id] == i, kegg_id] = re.split(" ", i)
    
    # replace NaN values with empty strings
    df.fillna('', inplace = True)
    
    # loop through all genes and add annotation
    for index, row in df.iterrows():
        
        # construct new dict with SBML conformity
        new_annot = {
            'uniprot': row['Entry'],
            'kegg.genes': 'reh:' + row['Gene names  (ordered locus )'],
            'ncbiprotein': dict_ncbi[row['Gene names  (ordered locus )']],
            'protein_name': row['Protein names'],
            'length': row['Length'],
            'mol_mass': row['Mass']
            }
        # assign new annotation to each gene
        gene = model.genes.get_by_id(row['Gene names  (ordered locus )'])
        gene.annotation = new_annot
    
    # final reporting
    print(' ----- SUMMARY OF GENE ANNOTATION ----- ')
    print('updated annotation for genes/proteins: ' + str(len(df)))


# ADDING SBO TERMS FOR METAB AND REACTIONS -----------------------------
# 
# SBO terms are a controlled vocabulary describing models or systems 
# biology tools. They help to characterize all entities of a model
# using well defined and fixed terms 
def update_sbo_terms(model):
    
    # ---------- METABOLITES ----------
    # add the same SBO term for 'simple chemical' for all metabolites
    # as recommended by memote
    for met in model.metabolites:
        met.annotation['sbo'] = 'SBO:0000247'
    print('...added SBO terms for ' + str(len(model.metabolites)) + ' metabolites')
    
    # ----------  REACTIONS  ----------
    # add a more diverse set of SBO terms for all reactions,
    # according to memote.io:
    #   "SBO:0000176 represents the term 'biochemical reaction'. Every 
    #   metabolic reaction that is not a transport or boundary reaction 
    #   should be annotated with this."
    #   "SBO:0000627 represents the term 'exchange reaction'. The Systems 
    #   Biology Ontology defines an exchange reaction as follows: 'A 
    #   modeling process to provide matter influx or efflux to a model'"
    #   "'SBO:0000185' and 'SBO:0000655' represent the terms 
    #   'translocation reaction' and 'transport reaction', respectively, 
    #   in addition to their children"
    for rea in model.reactions:
        # SBO term for exchange/boundary reactions
        if rea.boundary:
            rea.annotation['sbo'] = ['SBO:0000627']
        # SBO term for transport reactions
        elif rea.compartments == {'e', 'c'}:
            rea.annotation['sbo'] = ['SBO:0000655']
        # SBO term for all other biochemical reactions
        else: rea.annotation['sbo'] = ['SBO:0000176']
    print('...added SBO terms for ' + str(len(model.reactions)) + ' reactions')
    
    # ----------    GENES    ----------
    # add the SBO term "SBO:0000243" for 'gene' as recommended by memote
    for gene in model.genes:
        if 'sbo' in gene.annotation.keys():
            if not 'SBO:0000243' in gene.annotation['sbo']:
                gene.annotation['sbo'] = gene.annotation['sbo'] + ['SBO:0000243']
        else: gene.annotation['sbo'] = ['SBO:0000243']
    print('...added SBO terms for ' + str(len(model.genes)) + ' genes')
    
    # -------- BIOMASS REACTION -------
    # add the SBO term "SBO:0000243" for the biomass reaction as recommended by memote
    for biom in model.reactions.query('[Bb]iomass|BIOMASS'):
        biom.annotation['sbo'].append('SBO:0000629')


# TESTING GROWTH WITH FBA ----------------------------------------------
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

