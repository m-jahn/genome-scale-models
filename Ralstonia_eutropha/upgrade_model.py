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
    # update model IDs, names, and gene_reaction_rules for the reaction that
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
    model.reactions.SUCSD1.id = 'SSALy'
    model.reactions.SSALy.name = 'Succinate-semialdehyde dehydrogenase (NADP)'
    model.reactions.SUCCt.id = 'SUCCt2r'
    model.reactions.SUCCt2r.name = 'Succinate transport via proton symport'
    model.reactions.SUCCt2r.bounds = [-1000.0, 1000.0]
    
    # remove all duplicated reactions
    duplicated_reactions = (
        ['TYRTA2', 'ASPAM5', 'AMDS4', 'BZACCOAT2', '4NPHP2', 'HACOAD3', 
        'G1PTT2', 'NA1t2', 'NITRT', 'SPMS2', 'PROD3', 'PROD4', 'NO3RUQ2',
        'UDPG4E2', 'ADNK2', 'URIDK2', 'SERDHT1', 'CYSST3', 'ASPAM6',
        'GLUDH2', 'GLUDH4', 'CYSST2', 'ACFM4', 'ACCSYN1', 'ACCSYN2',
        'LAAO1', 'PHEALDD', 'BITCB', 'GLNAS', 'SUCSD2', 'GABAT1',
        'SUCCet'])
    
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
    
    # adjust names for important reactions in core metabolism manually
    model.reactions.TRKT1.id = 'TKT1'
    model.reactions.TRKT2.id = 'TKT2'
    model.reactions.PTA.id = 'PTAr'
    model.reactions.PTAr.name = 'Phosphotransacetylase'
    model.reactions.SUCCD1.id = 'SUCDi'
    model.reactions.SUCDi.name = 'Succinate dehydrogenase (irreversible)'
    model.reactions.SUCCD2.id = 'SUCD4'
    model.reactions.SUCCD3.id = 'SUCD1'
    model.reactions.SUCD1.name = 'Succinate dehydrogenase'
    model.reactions.BKAR1.id = 'PCNO'
    model.reactions.PCNO.name = 'Propanoyl-CoA:NADP+ 2-oxidoreductase'
    model.reactions.MNAO1.id = 'AACTOOR'
    model.reactions.AACTOOR.name = 'Aminoacetone:oxygen oxidoreductase(deaminating)(flavin-containing)'
    model.reactions.ALDRD1.id = 'LALDO2x'
    model.reactions.LALDO2x.name = 'D-Lactaldehyde:NAD+ 1-oxidoreductase'
    model.reactions.ALDRD2.id = 'LCARR'
    model.reactions.LCARR.name = 'Lacaldehyde reductase (R-propane-1,2-diol forming)'
    model.reactions.CABPS.id = 'CBPS'
    model.reactions.CBPS.name = 'Carbamoyl-phosphate synthase (glutamine-hydrolysing)'
    model.reactions.BAPYRT.id = 'APATr'
    model.reactions.APATr.name = 'B alanine pyruvate aminotransferase'
    model.reactions.DATA4.id = 'ALATA_D'
    model.reactions.DALATA.id = 'ALATA_D2'
    model.metabolites.uq_c.name = 'Ubiquinone-8'
    model.metabolites.uqh2_c.name = 'Ubiquinol-8'
    model.metabolites.bala_c.id = 'ala_B_c'
    model.metabolites.ala_B_c.name = 'Beta-Alanine'
    model.metabolites.get_by_id('3opp_c').id = 'msa'
    model.metabolites.msa.name = 'Malonate semialdehyde'
    
    
    # 3rd STEP: correct known errors in reactions
    # -------------------------------------------
    
    # Add one reaction that is missing for full heme biosynthesis, the
    # Uroporphyrinogen-III synthase. The reaction was taken from iJO1366.
    # The associated gene was annotated in uniprot and is already present in the model
    # as Uroporphyrinogen-III methyltransferase
    UPP3S = cobra.Reaction('UPP3S')
    UPP3S.name = 'Uroporphyrinogen-III synthase'
    model.add_reactions([UPP3S])
    model.reactions.UPP3S.build_reaction_from_string('hmb_c --> h2o_c + uppg3_c')
    model.reactions.UPP3S.gene_reaction_rule = 'H16_A2919'
    print('...added reaction Uroporphyrinogen-III synthase for heme biosynthesis')
    
    # add heme to the biomass sub-reaction cofactors and vitamins
    # the final concentration in mmol was adopted from E. coli reference model iJO1366
    model.reactions.Cofactors_and_vitamins.add_metabolites({'hemeA_c': -0.0074})
    print('...added heme to the Cofactor branch of biomass reaction')
    
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
    
    # reaction P5CD2 seems to have wrong reactants. There is no info on a
    # direct conversion from L-glu to L-glu-5-semialdehyde. All Bigg pathways
    # to L-glu-5-semialdehyde go via GLU5K and G5SD, the canonical pathway.
    # reaction LEUD3 seems to be wrong as well. It is annotated as L-leucine 
    # dehydrogenase, but converts (S)-3-Methyl-2-oxopentanoate to L-isoleucin,
    # although this is not the canonical pathway. The reaction does not seem to exist.
    # reaction 4AMBUAT is annotated as 4-aminobutyrate aminotransferase, but
    # catalyzes transfer of an aminogroup from malonyl semialdehyde to beta-alanine
    # it forms an artificial cycle with APAtr and does not seem to exist in nature
    model.remove_reactions(['P5CD2', 'LEUD3', '4AMBUAT'])
    
    # some reaction show infinitely high flux when doing FVA. This is often
    # related to artificial cycles that can form when reaction directionalities
    # are not correct. Here we correct some known errors, either based on
    # thermodynamic constraints (see modelSEED dG and estimated directionality),
    # or participation in cycles (compare Bigg reactions).
    model.reactions.PPCSYN2.bounds = (-1000, 0)
    model.reactions.GLUR.bounds = (-1000, 0)
    model.reactions.ALAR.bounds = (0, 1000)
    model.reactions.LAAO3.bounds = (0, 1000)
    model.reactions.ASPO1.bounds = (0, 1000)
    model.reactions.ALATA_D.bounds = (0, 1000)
    model.reactions.ALATA_D2.bounds = (0, 1000)
    
    # One reason that Succ-dehydrogenase is not essential in the model is
    # a shortcut called the methylglyoxal pathway, where DHAP (or here, via 
    # another intermediate, ac-Coa) gets metabolized to methylglyoxal, 
    # lactaldehyde, lactate, and pyruvate under glucose overflow conditions
    # However this pathway contains several very toxic intermediates and is 
    # unlikely to carry more than minimal flux (see wikipedia: methyglyoxal pathway)
    # The reaction that the model contains is completely incorrect and
    # can be removed (converts mthglx directly to pyr without neccessary 
    # intermediates). The key reaction making mthglx was silenced.
    model.remove_reactions(['LCTAD2'])
    model.reactions.AACTOOR.bounds = [0.0, 0.0]
    print('...removed erroneous reaction in methyglyoxal pathway')
    
    # Similar reactions that allow succ-DH or Akg-DH to be bypassed in the TCA
    # are the GABA shunt, that really exists but was never shown to carry major 
    # flux (a natural putrescine degradation pathway). Or the Prolyl-4 
    # hydroxylase that converts 2OG to Succ.
    model.reactions.GABAT2.bounds = [0.0, 0.0]
    model.reactions.P4HX.bounds = [0.0, 0.0]
    print('...silence reactions for GABA shunt and Prolyl hydroxylase bypass')
    
    # Succinate dehydrogenase has incorrect reaction stoichiometry regarding H+
    model.reactions.SUCD1.subtract_metabolites({model.metabolites.h_c: -1})
    print('...Succinate dehydrogenase contains one H+ too much, was removed')
    
    # The reaction MICITL is the last step of the methyl-citrate cycle, an
    # alternative route through TCA from oaa + prop-coa --> succ + pyr.
    # It carries artificially high flux, so flux of the final reaction was constrained.
    model.reactions.MICITL.bounds = (0.0, 0.0)
    print('...silenced reaction MICITL, the last step of the methyl citrate cycle')
    
    # The reaction MCTOP is annotated as 'unclear reation' and seems to be
    # erroneous compared with a Bigg reference.
    model.reactions.MCTOP.id = 'OCDOR'
    model.reactions.OCDOR.name = '3-oxopropionyl-CoA:NAD+ oxidoreductase'
    model.reactions.OCDOR.build_reaction_from_string('h_c + malcoa_c + nadph_c --> h2o_c + nadp_c + 3oppcoa_c')
    
    # The following reaction OPTHP is also annotaed as 'unclear reaction',
    # although the reaction itself is correct.
    model.reactions.OPTHP.id = 'HPCOR'
    model.reactions.HPCOR.name = '3-hydroxypropionyl-CoA:NADP+ oxidoreductase'
    
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
    
    # Ornithine cyclodeaminase is a link for a shortcut from succinate to fumarate
    # via L-proline, a very artifical reaction that bypasses succinate dehydrogenase
    # The reaction should not run in both direction according to other Bigg models
    model.reactions.ORNC.bounds = (0.0, 1000.0)
    print('...set irreversibility for reaction ORNCD')
    
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
    
    # Add a (dummy) tRNA loading reaction for asparagin that was missing;
    # in the original model, tRNA-Asn is made from tRNA-Asp directly.
    # However this does not work properly with RBA models
    model.add_metabolites(
    cobra.Metabolite(
        id = 'trnaasn_c',
        name = 'tRNA(Asn)',
        compartment = 'c',
        charge = 0,
        formula = 'C10H17O10PR2'))
    
    ASNTRS = cobra.Reaction('ASNTRS')
    ASNTRS.name = 'Asparaginyl-tRNA synthetase'
    ASNTRS_string = 'asn__L_c + atp_c + trnaasn_c --> amp_c + asntrna_c + h_c + ppi_c'
    model.add_reactions([ASNTRS])
    model.reactions.ASNTRS.build_reaction_from_string(ASNTRS_string)
    model.reactions.ASNTRS.gene_reaction_rule = 'H16_A0453'
    print('...added missing Asn-tRNA loading reaction (Asparaginyl-tRNA synthetase)')
    
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
    model.reactions.HCO3E.gene_reaction_rule = 'H16_A0169 or H16_B2270 or H16_B2403'
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
    model.reactions.PRUK.gene_reaction_rule = 'H16_B1389 or PHG421'
    model.reactions.RBPC.gene_reaction_rule = '( PHG426 and PHG427 ) or ( H16_B1394 and H16_B1395 )'
    
    # define reactions from string
    model.reactions.PRUK.build_reaction_from_string('atp_c + rl5p_c --> adp_c + h_c + rb15bp_c')
    model.reactions.RBPC.build_reaction_from_string('co2_c + h2o_c + rb15bp_c --> 2.0 3pg_c + 2.0 h_c')
    
    # silence original lumped reaction for CBB cycle
    model.reactions.CBBCYC.bounds = (0.0, 0.0)
    print('...added PRUK and RBPC (Rubisco) reactions instead of lumped reaction')
    
    # The model contains a reaction to synthesize a precursor of the essential
    # cofactor thiamin monophosphate (thmmp) using a non-exisitng pseudo reaction
    # called 'MISRXN', named in the model as 'unclear reaction'. This reaction
    # has to be removed and replaced by a correct Thiazole phosphate synthesis
    # reaction (see Bigg reference reaction THZPSN)
    model.remove_reactions(['MISRXN'])
    print('...removed erroneous thiamin precursor synthesis reaction')
    
    # For this reaction a new intermediate metabolite, 4-hydroxy benzyl 
    # alcohol, is required for which metabolizing reactions exist in the model
    M_4hba = cobra.Metabolite('4hba_c')
    M_4hba.compartment = 'c'
    M_4hba.name = '4-Hydroxy-benzyl alcohol'
    M_4hba.formula = 'C7H8O2'
    M_4hba.charge = 0
    model.add_metabolites([M_4hba])
    
    # However, one reaction that connects the metabolite with the next
    # intermediates in the model must be added, a dehydrogenase oxidizing 
    # the benzyl alcohol to benzyl aldehyde
    R_4HBADH = cobra.Reaction('4HBADH')
    R_4HBADH.name = '4 hydroxy benzyl alcohol dehydrogenase'
    R_4HBADH.lower_bound = -1000.0
    R_4HBADH.upper_bound = 1000.0
    model.add_reactions([R_4HBADH])
    model.reactions.get_by_id('4HBADH').build_reaction_from_string('nad_c + 4hba_c <=> h_c + nadh_c + 4hbzald_c')
    
    # Add reaction for Thiazole phosphate synthesis.
    # The involved genes are thiG (H16_A0238), thiF (H16_A0330), thiS (H16_A0237),
    # and thiH (?). According to Uniprot, ThiG forms heterodimers with either 
    # ThiH or ThiS.
    THZPSN = cobra.Reaction('THZPSN')
    THZPSN.name = 'Thiazole phosphate synthesis'
    THZPSN.lower_bound = 0.0
    THZPSN.upper_bound = 1000.0
    model.add_reactions([THZPSN])
    model.reactions.THZPSN.build_reaction_from_string('atp_c + cys__L_c + dx5p_c + tyr__L_c --> thzp_c + ala__L_c + amp_c + co2_c + h_c + h2o_c + ppi_c + 4hba_c')
    model.reactions.THZPSN.gene_reaction_rule = '( H16_A0238 and H16_A0237 and H16_A0330 )'
    print('...added correct thiamin precursor synthesis reaction')
    
    # final reporting
    print(' ----- SUMMARY OF MODIFIED REACTIONS ----- ')
    print('removed duplicated reactions: ' + str(len(duplicated_reactions)))
    print('added names to reactions: ' + str(count))
    print('corrected erroneous reactions: 31')


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
            # list of possible new ID
            new_id_list = bigg.loc[hit, 'id'].to_list()
            # if the old ID is present, don't change anything
            # if not, choose a new ID, and here prefer the shortest/
            # most concise one
            if item.id in new_id_list:
                new_id = item.id
            else:
                new_id_length = [len(i) for i in new_id_list]
                new_id = new_id_list[new_id_length.index(min(new_id_length))]
            if item_type == 'metabolite':
                new_id = new_id + '_' + item.compartment
            # proceed if A) new ID is not the old ID and 
            #            B) the new ID doesnt exist yet
            cond = new_id not in model_item_type.list_attr("id")
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
            print('...updated ' + item_type + ' annotation for ' + item.id)
    
    
    # final reporting
    print(' ----- SUMMARY OF ' + item_type.upper() + ' CHANGES ----- ')
    print('updated ' + item_type + ' IDs: ' + str(count_id))
    print('updated ' + item_type + ' annotation: ' + str(count_an))


# UPDATE ENZYME STOICHIOMETRIES BASED ON BARSEQ ------------------------
#
# BarSeq data available for R. eutropha in our lab contains important
# information about possible essential genes. Essential genes in the
# model were predicted using the script in 'essentiality_analysis.py'.
# The result was compared with the list of essential genes derived from
# BarSeq (cultivation of a transposon knockout library). The updates
# derived from that are implemented in the following function.

def update_stoichiometry(model):
    
    # This is the preferred way to annotate complex stoichiometry 
    # (quoted from RBApy manual, https://sysbioinra.github.io/RBApy/usage.html#sec2.3.2.3):
    # "RBApy assumes that the boolean relation is always “or”s of “and”s."
    # Example: ( STM3795  and  STM3796 )  or  ( STM3901  and  STM3902 )
    
    # STEP 1: correct false positive essential genes
    # ----------------------------------------------
    # the sulfate adenylyltransferase, reaction SADT, requires 4 
    # enzyme subunits to be present. This is a mistake, only 2 subunits are required,
    # subunit 1 is CysN (H16_A2995, H16_B0626), and subunit 2 is cysD 
    # (H16_A2996, H16_B0627), making none of the genes essential anymore
    model.reactions.SADT.gene_reaction_rule = '( H16_A2995 and H16_A2996 ) or ( H16_B0626 or H16_B0627 )'
    
    # the sulfite reductase, reaction SLFR, requires 3 enzyme subunits 
    # to be present. This is a mistake, subunits CysI1 (H16_A2999) and 
    # CysI2 (H16_A1639) seem to be isoenzymes of the beta subunit.
    # Only the alpha subunit (H16_B2500) seems to be essential (compare E.coli, uniprot)
    model.reactions.SLFR.gene_reaction_rule = '( H16_A2999 and H16_B2500) or ( H16_A1639 and H16_B2500 )'
    
    # STEP 2: correct false negative essential genes
    # ----------------------------------------------
    # On top of the list are reactions involving 7 NADH-DH subunits 
    # all turning out to be essential by BarSeq. 6 more NADH-DH subunits
    # have low BarSeq counts but *could* be essential, i.e. are at least
    # not part of the non-essential group and ahve not enough verfied
    # Tn integrations. Therefore the best assumption is to set all 
    # NADH subunits as essential
    model.reactions.NADHDH.gene_reaction_rule = (
        '( H16_A1051 and H16_A1052 and H16_A1050 and H16_A1055 and ' + 
        'H16_A1056 and H16_A1053 and H16_A1054 and H16_A1061 and ' +
        'H16_A1060 and H16_A1063 and H16_A1062 and H16_A1059 and ' +
        'H16_A1058 and H16_A1057 and H16_A0251 )'
    )
    model.reactions.NADH16.gene_reaction_rule = model.reactions.NADHDH.gene_reaction_rule
    model.reactions.NADH5.gene_reaction_rule = model.reactions.NADHDH.gene_reaction_rule
    
    # The next top candidate is ATPsynthase. Here the error is obviously a single 
    # gene that was added erroneously using an OR rule. In addition to
    # the 5 BarSeq essential genes, further 2 subunits are probably essential
    # but have too low reads/insertions. With those it'd be 7 out of 8 
    # annotated genes that are likely essential.
    model.reactions.ATPS4m.gene_reaction_rule = (
        '( H16_A3643 and H16_A3642 and H16_A3639 and H16_A3636 and ' + 
        'H16_A3637 and H16_A3638 and H16_A3640 and H16_A3641 )'
    )
    
    # The cytochrome oxidase b0/c (reactions CYTBD, CYTCOBO3) has no genes annotated. However,
    # there are several complexes known in Ralstonia that can act as terminal
    # cytochrome c oxidases (analogous to complex 3 and 4 of respiratory chain), 
    # with a high level of redundancy between complexes.
    # For details, see excellent review by R. Cramm, J Mol Microbio & Biotech, 2008.
    # The complex 3 analog (bc1 type) is the cytochrome c reductase and annotated with 3 genes,
    # qrcA, B, C. These are by far the most expressed ones judging from proteomics data.
    # (although results must be taken with care because membrane proteins are isolated/quantified
    # with low efficiency). This complex provides reduced cytochrome c for a range of
    # other terminal cyt c oxidases, of which the cta (aa3 type), cox (bb3 type) and cco (cbb3 type) operons 
    # seem to play the biggest role.
    model.reactions.CYTBD.gene_reaction_rule = (
        '( H16_A3396 and H16_A3397 and H16_A3398 and ' +        # qrcA, B, C
        'H16_A2319 and H16_A2318 and H16_A2316 and ' +          # ccoN, O, P
        'H16_B2062 and H16_B2059 and ' +                        # coxM, P
        'H16_A0342 and H16_A0343 and H16_A0347 and H16_A0345 )' # ctaC, D, E, G
    )

    # The complex 4 analog (b03 type) is another cytochrome c terminal oxidase and has three 
    # redundant complexes annotated in R. eutropha, with very low expression,
    # cyoA1, B1, C1, D1, cyoA2, B2, C2, D2, and cyoA3, B3, C3, D3. They don't
    # seem to play a major roll in any of the tested conditions
    model.reactions.CYTCOBO3.gene_reaction_rule = (
        '( H16_A1071 and H16_A1072 and H16_A1073 and H16_A1074 ) or ' +  # cyoA1,B1,C1,D1
        '( H16_A1640 and H16_A1641 and H16_A1642 and H16_A1643 ) or ' +  # cyoA2,B2,C2,D3
        '( H16_B1025 and H16_B1026 and H16_B1027 and H16_B1028 )'        # cyoA3,B3,C3,D3
    )
    
    # Succinate dehydrogenase, two reactions (SUCD1 and SUCDi) with 5 genes
    # of which 4 are essential according to BarSeq (H16_A2632, H16_A2631, 
    # H16_A2630 and H16_A2629). According to model simulations, however, 
    # they are only essential under growth on succinate, not LB, fructose,
    # or probably other C sources that don't involve later part of TCA
    # Several reactions related to unlikely TCA short cuts were removed
    # or corrected, particularly related to the methylglyoxal shunt.
    # (see modify_reactions()). This doesn't make SUCD1 and SUCDi 
    # essential but will yield considerably lower growth rate under growth 
    # on fructose ('beneficial' but not essential reaction).
    
    
    # Acetyl-Coa carboxylase (reaction id ACCOAC) is an important initial
    # step in fatty acid synthesis (Ac-Coa --> mal-Coa). In E. coli the enzyme 
    # has 4 essential subunits (AND rule), in the model the gene reaction 
    # rule is less clear. We use the four canonical subunits from Uniprot
    # (AccA, B, C2, D) that are all essential in BarSeq
    model.reactions.ACCOAC.gene_reaction_rule = (
        '( H16_A1223 and H16_A3171 and H16_A3172 and H16_A2611 )'
    )
    
    # Cytosolic lactate dehydrogenase is essential with both subunits 
    # (H16_A1681 and H16_A1682 annotated as one enzyme in Uniprot), according to
    # BarSeq. We can update the gene rules for cytosolic and cytochrome
    # associated LDH.
    model.reactions.LDHm.gene_reaction_rule = 'H16_A3091'
    model.reactions.DLDHD.gene_reaction_rule = (
        '( H16_A1681 and H16_A1682 )'
    )
    
    # reaction 2-dehydropantoate 2-reductase (DPR) has two out of three
    # associated genes marked as essential by BarSeq (H16_B1769, H16_A1715
    # but not H16_B1719). These three seem to be isoenzymes (first two identical
    # except for 1 AA). The OR rule seems to be justified.
    
    # reaction Carbamoylphosphate synthase (glutamine hydrolysing) (CBPS)
    # is essential in BarSeq but was not in the model because of an
    # erroneous OR rule. In facto both subunits seem to be required 
    # (see comparable reaction in E.coli, Bigg)
    model.reactions.CBPS.gene_reaction_rule = (
        '( H16_A2454 and H16_A2452 )'
    )
    
    # reaction CATL (2.0 3han_c + 2.0 o2_c --> cvn_c + h2_c + 2.0 h2o2_c) 
    # is associated with 3 catalase genes that should rather be 
    # associated with the canonical catalase (CAT, 2.0 h2o2_c --> 2.0 h2o_c + o2_c)
    # instead. One of them is essential (H16_B1428, katE2), the others not.
    # The rule is updated to reflect that the essential one is required.
    model.reactions.CAT.gene_reaction_rule = (
        '( H16_B1428 and H16_A2777 ) or ( H16_B1428 and H16_A3109 )'
    )
    
    # Two genes associated with Citrate synthase (CS) are essential
    # according to BarSeq (H16_A2627, H16_B0357) while 3 others or not.
    # The rule is updated to reflect that the essential ones are required.
    model.reactions.CS.gene_reaction_rule = (
        '( H16_A2627 and H16_B0357 and H16_B2211 ) or ' +
        '( H16_A2627 and H16_B0357 and H16_A1229 ) or ' +
        '( H16_A2627 and H16_B0357 and H16_B0414 )'
    )
    
    # Reaction dehydrofolate reductase (DHFR) has both associated genes
    # marked as essential in BarSeq. Uniprot query reveals that 
    # H16_A2704 (folA2) and H16_A1840 (folA1) do not align i.e. are no
    # isoenzymes but subunits. Change OR to AND rule.
    model.reactions.DHFR.gene_reaction_rule = (
        '( H16_A2704 and H16_A1840 )'
    )
    
    # 4 reactions for Ribonucleoside-diphosphate reductase for all 4 
    # nucleotides are carried out by the same 3 subunits, of which two
    # are the canonical subunits nrdA (H16_A3234) and nrdB (H16_A3235).
    # These two are essential in BarSeq data, while the third enzyme 
    # (H16_A2390, nrdJ) seems to be dispensable. Updated reaction rule
    # accordingly.
    RNDrule = '( H16_A3234 and H16_A3235 ) or ( H16_A3234 and H16_A3235 and H16_A2390 )'
    model.reactions.RNDR1.gene_reaction_rule = RNDrule
    model.reactions.RNDR2.gene_reaction_rule = RNDrule
    model.reactions.RNDR3.gene_reaction_rule = RNDrule
    model.reactions.RNDR4.gene_reaction_rule = RNDrule
    
    # reaction Ribulose 5-phosphate 3-epimerase (RPE) has 2 out of 3 
    # associated genes marked as essential by BarSeq (rpe, H16_A3317, rpe1 H16_B1391)
    # The megaplasmid borne gene PHG423 (rpe2) is an isoenzyme of rpe1
    # and not essential. Rpe1 and rpe2 are perfect isoenzymes (95.% idendity)
    # while . 
    # It is surprising that the two less-related subunits are clearly essential, pointing
    # to diverging functions. Nevertheless the model should 
    # reflect the sequence identity of the genome and plasmid isoenzymes.
    model.reactions.RPE.gene_reaction_rule = (
        '( H16_B1391 and H16_A3317 ) or ( PHG423 and H16_A3317 )'
    )
    
    # reaction transketolase has (TKT1, TKT2) has three associated genes
    # of which 2 are essential (H16_B1388, H16_A3147) and one (PHG420) 
    # is a perfect iso-enzyme of H16_B1388 located on the megaplasmid 
    # (96.6% identical). All three are isoenzymes but as with RPE, the two
    # less-similar enzymes are essential pointing to different roles. The
    # plasmid isoenzyme is not essential, maybe expression from plasmid is
    # generally lower. We treat PHG420 and H16_B1388 as perfect isoenzymes
    TKTrule = '( H16_B1388 and H16_A3147 ) or ( PHG420 and H16_A3147 )'
    model.reactions.TKT1.gene_reaction_rule = TKTrule
    model.reactions.TKT2.gene_reaction_rule = TKTrule
    
    # not directly related to stoichiometry, but anyway: 
    # correct name of dna_c to DNA_c
    model.metabolites.dna_c.id = 'DNA_c'
    model.metabolites.DNA_c.name = 'DNA'
    
    # NOTES:
    # ------
    #
    # Some reactions/genes marked as essential using BarSeq are tRNA-
    # synthetases responsible to load empty tRNAs with their corresponding
    # amino acid. The model contains a complete set of tRNA synthetase
    # reactions but currently uses amino acids directly in protein bio
    # synthesis. This could be changed later to reflect that protein
    # biosynthesis actually uses loaded tRNAs (making those reactions 
    # also essential in silico)
    #
    # Some reactions in the model contain an excessive number of genes 
    # associated with them, such as
    # - Enoyl-CoA hydratase
    # - Acetyltransferase
    # - Aldehyde dehydrogenase (butanal, NAD)
    
    # It is not clear yet what is the best way to handle those. The 
    # excessive list of genes with 'OR' rule is a consequence of very general 
    # and unspecific annotation of genes (example: putative acyl-transferase).
    # One option is to ignore this problem, the other is the minimalistic 
    # approach and weed out all non-essential genes associated with an essential 
    # reaction, thus leaving only Barseq-curated genes in.
    
    # final reporting
    print(' ----- SUMMARY OF STOICHIOMETRY CHANGES ----- ')
    print('updated gene reaction rules for 20 reactions')


# ADDING GENE ANNOTATIONS FROM UNIPROT ---------------------------------
# 
# for this purpose we use the very comprehensive bioservices
# package for python that has connectivity to all possible databases
def get_gene_annotation(
    ref_path = 'data/gene_reference.json',
    item_list = [],
    item_type = 'gene'):
    
    # open client connection
    uni = UniProt(verbose = False)
    kegg = KEGG(verbose = False)
    
    # construct generic reference names based on input params
    ref_name = item_type + '_reference'
    kegg_id = 'Gene names  (ordered locus )'
    
    # check if already an annotation file exists
    # if yes load, if no, create new one
    if path.exists(ref_path):
        reference = pd.read_json(ref_path)
        
        # determine which items are missing in reference
        ref_list = reference[kegg_id].to_list()
        item_list = [i for i in item_list if i not in ref_list]
    
    # download annotation for unannotated genes
    # and add to reference dataframe
    if len(item_list):
        
        # fetch uniprot annotation table for all genes as pandas df
        print('...downloading annotation for ' + str(len(item_list)) + ' genes from uniprot.org')
        df = uni.get_df(item_list, limit = len(item_list))
        print('...downloading annotation for ' + str(len(item_list)) + ' genes from kegg.jp')
        dict_ncbi = {e: kegg.parse(kegg.get("reh:" + str(e)))['DBLINKS']['NCBI-ProteinID'] for e in item_list}
        # remove entries that are not according to Reh standard
        df = df[[len(i) == 6 for i in df['Entry']]]
        
        # some rows are merged entries for 2 KEGG IDs
        # in this case we simply duplicate the rows and split entries
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
        
        # merge processed uniprot and kegg data
        df['ncbiprotein'] = [dict_ncbi[i] for i in df[kegg_id].to_list()]
        
        # merge newly downloaded refs with existing refs
        if path.exists(ref_path):
            reference = pd.concat([reference, df], ignore_index = True)
        else:
            reference = df
        
        # refactor index (row numbers)
        reference.index = range(0, len(reference))
        
        # export reference with added items as json file
        reference.to_json(ref_path)
        print('...exported file "' + ref_path + '" with ' + str(len(item_list)) + ' new ' + item_type + 's')
    
    return(reference)


def update_gene_annotation(model):
    
    # add names to  gene IDs where they are missing
    for i in model.genes:
        if (i.name == '') or (re.match('^G_', i.name)):
            i.name = re.sub('^G_', '', i.id)
            i.id = i.name
            print(i.name)
    
    df = get_gene_annotation(
        ref_path = 'data/gene_reference.json',
        item_list = model.genes.list_attr("id"),
        item_type = 'gene')
    
    # loop through all genes and add annotation
    for index, row in df.iterrows():
        
        # construct new dict with SBML conformity
        new_annot = {
            'uniprot': row['Entry'],
            'kegg.genes': 'reh:' + row['Gene names  (ordered locus )'],
            'ncbiprotein': row['ncbiprotein'],
            'protein_name': row['Protein names'],
            'length': row['Length'],
            'mol_mass': row['Mass']
            }
        # assign new annotation to each gene
        gene = model.genes.get_by_id(row['Gene names  (ordered locus )'])
        gene.annotation = new_annot
    
    
    # get previously downloaded genome annotation for R. eutropha
    # (source: uniprot.org)
    genome_db = pd.read_csv('data/genome_annotation.csv')
    
    # loop through all reactions, get EC number where available and
    # add genes corresponding to EC number to reaction.genes
    count_gene = 0
    # exclude some reactions that were updated manually before
    rea_excluded = ['PRUK', 'RBPC', 'HCO3E', 'TKT1', 'TKT2',
        'SADT', 'SLFR', 'NADHDH', 'NADH16', 'NADH5', 'ATPS4m',
        'ACCOAC', 'LDHm', 'DLDHD', 'CBPS', 'CAT', 'CS', 'DHFR', 
        'RNDR1', 'RNDR2', 'RNDR3', 'RNDR4', 'RPE']
    
    for rea in model.reactions:
        
        if ('ec-code' in rea.annotation) and (rea.id not in rea_excluded):
            ec_code = rea.annotation['ec-code']
            # only add gene association to reactions with unique EC number
            # to avoid gene inflation. Only match complete EC codes
            if (len(ec_code) == 1) and (re.match('[0-9]+\.[0-9]+\.[0-9]+\.([0-9]+)$', ec_code[0])):
                ec_code = ec_code[0]
                genes_db = genome_db[genome_db['EC_number'] == ec_code]['locus_tag'].to_list()
                genes_db = list(set(genes_db))
                genes_model = [re.sub('^G_', '', i.name) for i in rea.genes]
                # determine which items are missing in model
                genes_new = [i for i in genes_db if i not in genes_model]
                if len(genes_new):
                    if len(rea.gene_reaction_rule):
                        sep = ' or '
                    else:
                        sep = ''
                    rea.gene_reaction_rule = rea.gene_reaction_rule + sep + ' or '.join(genes_new)
                    print('...added genes: ' + str(genes_new) + ' to reaction: ' + rea.id + ' (' + ec_code + ')')
                    count_gene = count_gene + 1
    
    # manual changes not present in uniprot (source: kegg)
    model.reactions.GAPD.gene_reaction_rule = (
        model.reactions.GAPD.gene_reaction_rule + ' or PHG418')
    count_gene = count_gene + 1
    
    # final reporting
    print(' ----- SUMMARY OF GENE ANNOTATION ----- ')
    print('updated annotation for genes/proteins: ' + str(len(df)))
    print('added genes IDs to reactions: ' + str(count_gene))


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
        'EX_fe2_e': 10.0,
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

