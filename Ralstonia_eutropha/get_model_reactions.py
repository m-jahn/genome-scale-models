#
# This script fetches the model's reactions and annotation
# and exports them as a csv file
# 
# LIBRARIES ------------------------------------------------------------

import pandas as pd
import csv
import cobra


# function to retrieve annotation objects and return empty
# string if not available
def get_annot(reaction, key):
    if key in reaction.annotation:
        return(reaction.annotation[key])
    else: return('')


# function to retrieve group membership for reaction and return empty
# string if not available
def get_groups(reaction, groups):
    # collect group memberships
    reaction_groups = []
    for g in groups:
        members = [m.id for m in list(g.members) if m.id == reaction.id]
        if members:
            reaction_groups.append(g.id)
    # return groups collapsed to one string
    if len(reaction_groups):
        return("; ".join(reaction_groups))
    else: return('')


# MAIN FUNCTION --------------------------------------------------------
def main():
    
    # adjust file paths if necessary
    input = "sbml/RehMBEL1391_sbml_L3V1.xml"
    output = 'simulations/essentiality/model_reactions.csv'
    
    # import SBML model
    model = cobra.io.read_sbml_model(input)
    
    # create a new pandas df for results
    df = pd.DataFrame(columns = ([
        'reaction_id', 
        'reaction_name',
        'reaction_reaction',
        'genes'])
        )
    
    # query information for all reactions in the model
    for r in model.reactions:
        df = df.append({
            'reaction_id': r.id,
            'reaction_name': r.name,
            'reaction_reaction': r.reaction,
            'genes': ", ".join([i.id for i in r.genes]),
            'EC_number': get_annot(r, 'ec-code'),
            'kegg_reaction': get_annot(r, 'kegg.reaction'),
            'metanetx_reaction': get_annot(r, 'metanetx.reaction'),
            'seed_reaction': get_annot(r, 'seed.reaction'),
            'groups': get_groups(r, model.groups)},
            ignore_index = True
            )
    
    # save results from pandas data frame to hdd
    df.to_csv(output)


if __name__== "__main__" :
    main()
