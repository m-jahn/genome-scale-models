# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +    IMPORT, MODIFICATION, AND TESTING OF RALSTONIA EUTROPHA GEM     +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Script to import and modify an SBML (systems biology markup
# language) model which describes the metabolism of a bacterial cell.
# The script perform the following steps: 
# 
#  - import cobra and other required libraries
#  - import an SBML model
#  - check key data of the model
#  - change constraints/reactants of erroneous reactions
#  - add reactions that were missing
#  - run a basic flux balance analysis (FBA)
# 
# 
# LIBRARIES ------------------------------------------------------------

import upgrade_model as um
import cobra

# IMPORT, MODIFY, AND EXPORT MODEL ------------------------------------- 
def main():
    
    # import model
    model = um.import_model(path = 'sbml/RehMBEL1391_sbml_L2V1.xml')
    
    # test energy generating cycles
    um.test_EGC(model.copy())
    
    # modify reactions
    um.modify_reactions(model)
    
    # test energy generating cycles again
    um.test_EGC(model.copy())
    
    # update metabolite annotation from Bigg db
    um.update_bigg_annotation(model, item_type = 'metabolite')
    
    # update reaction annotation from Bigg db
    um.update_bigg_annotation(model, item_type = 'reaction')
    
    # update stoichiometry of complexes
    um.update_stoichiometry(model)
    
    # update gene annotation from uniprot
    um.update_gene_annotation(model)
    
    # update SBO annotation for met, react, genes
    um.update_sbo_terms(model)
    
    # test run with FBA
    um.test_FBA(model)
    
    # save modified model as SBML Level 3, version 1
    cobra.io.write_sbml_model(model, 'sbml/RehMBEL1391_sbml_L3V1.xml')
    print('...model exported as SBML file')
    cobra.io.save_json_model(model, 'escher/RehMBEL1391_sbml_L3V1.json')
    print('...model exported as JSON file')

if __name__== '__main__' :
    main()
