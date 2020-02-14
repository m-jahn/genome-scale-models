# genome-scale-models

Genome scale metabolic models in SBML format

### Getting started

The genome scale metabolic models in this repository are encoded according to [SBML standard](sbml.org) and saved as human-readable `*.xml` or `*.json` files.

The models can have well over 1000 reactions, it is therefore recommended to work with these models using a frame work such as [COBRApy](https://cobrapy.readthedocs.io/en/latest/). To install COBRApy follow the instructions on its [github page](https://github.com/opencobra/cobrapy/blob/master/INSTALL.rst). Installation of additional python 3 dependencies might be necessary for full functionality, such as `libsbml`, `numpy`, 
`scipy`, or `pandas`.

```{b}
# for linux, run the following line in terminal
sudo apt install python3-pip
pip install cobra
```

To work with the model using `COBRApy`, we can import it in a python session. We can look at the number of reaction, metabolites and genes associated with the 
model reactions.

```{python}
# load libraries
import numpy as np
import pandas as pd
import tabulate
import cobra
import os
from os.path import join

# set the path to the model's directory
data_dir = '~/genome-scale-models/Ralstonia_eutropha/'
model = cobra.io.read_sbml_model(join(data_dir, "sbml/RehMBEL1391_sbml_L3V1.xml"))

# summary of the imported model
print('%i reactions' % len(model.reactions))
print('%i metabolites' % len(model.metabolites))
print('%i genes' % len(model.genes))
```


### Structure of SBML models

The `*.xml` files containing the model definition have three (four) major slots:

- metabolite definitions, such as `nadh_c` for NADH or `fru_c` for fructose. 
- reaction definitions that link metabolites to each other. Those define the topology of
  the metabolic network (for example `PYK`: `adp_c + pep_c --> atp_c + pyr_c`)
- groups of reactions that constitute a subsystem, such as all reactions for Citrate cycle
- meta information such as gene-name reaction associations


### Parameters and constraints

Genome scale models are constrained primarily by two things:
- network topology defined by the set of reactions. Models usually include only the reaction network for low molecular weight metabolites.
- reaction bounds, e.g. ranging from `-Inf` to `+Inf` for reversible reactions. Irreversible reactions would be bounded by `0 < f < +Inf`, so that the reaction is only allowed to proceed in one direction.

In python, we can look at the bounds of all exchange reactions (the ones supplying metabolites from the environment).

```{python}
# inspect exchange reactions
model.exchanges.list_attr("bounds")

# we can also change bounds for reactions
for reaction in model.exchanges:
    reaction.lower_bound = 0.0
    reaction.upper_bound = 1000.0
```

In contrast to other types of models such as simple [resource allocation models](https://github.com/m-jahn/cell-economy-models), genome scale models usually don't include cellular processes for production of macromolecules. In other words, transcription, translation, and DNA replication are not explicitly included in the model but only appear as abstract, lumped reactions.


### Objective function

The objective function of the model is a variable similar to other variables that are optimized by the solver. However, when solving the model the prime target of the algorithm is to maximize this variable, often the `Biomass` reaction. Units of all reactions are by default in mmol per gDCW per h, but since the biomass reaction is _per definition_ formulated such that 1 mmol biomass equals 1 g biomass, it also represents the specific growth rate μ (g biomass per gDCW per hour, biomass term can be eliminated).

In COBRApy we can tell the solver to optimize (usually: maximize) the flux through any reaction of choice. To maximize growth rate that would be the Biomass equation.

```{python}
# set objective function
model.objective = {model.reactions.Biomass: 1}
```


### Solving the model

Before we can pass the model to the solver and find the optimal flux distribution towards our goal, we have to define a growth medium (a set of exchange fluxes that represent the nutrients available in the outer environment of the cell).

```{python}
model.medium = {
    'EX_mg2_e': 10.0,
    'EX_pi_e': 100.0,
    'EX_cobalt2_e': 10.0,
    'EX_cl_e': 10.0,
    'EX_k_e': 10.0,
    'EX_fe3_e': 10.0,
    'EX_so4_e': 10.0,
    'EX_fru_e': 5.0,
    'EX_nh4_e': 10.0,
    'EX_na_e': 10.0,
    'EX_o2_e': 18.5,
    'EX_mobd_e': 10.0,
    'EX_h2o_e': 1000.0,
    'EX_h_e': 100.0
    }
```

The solver will then analyze the network and find the optimal steady state flux from our input metabolites to biomass.

```{python}
# run FBA analysis
solution = model.optimize()

# print solution summary, the status from the linear programming solver
print([solution, "status: ", solution.status])

# print top 10 forward and backward flux
fluxes = solution.fluxes.sort_values()
print(fluxes[0:10])
print(fluxes[len(fluxes)-10:len(fluxes)])

# quick summary of FBA analysis
print(model.summary())

# summary of energy balance
print(model.metabolites.atp_c.summary())

# summary of redox balance
print(model.metabolites.nadh_c.summary())
```


### Visualization of results

Coming soon.
