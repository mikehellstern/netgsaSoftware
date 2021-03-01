# netgsaSoftware
netgsa software submission to PLOS Computational Biology

This repository contains the code and data to recreate the results from the NetGSA software paper submitted to PLOS Computational Biology. The code and inputs are broken up into two separate folders. The `netgsa_new` folder contains code and datasets used to compare the old NetGSA to the new NetGSA. The `netgsa_comp` folder contains code and datasets used to compare the new NetGSA to other existing methods.

# netgsa_new

## Code

There are 4 different files corresponding to different output in the paper:

1. `NetGSA_Power.R` - Estimates power for REHE for all 3 datasets
2. `NetGSA_Power_reml_pathway.R` - Estimates power for REML pathway-by-pathway for BreastCancer and ProstateCancer
3. `NetGSA_Power_reml_pathway_metab.R` - Estimates power for REML pathway-by-pathway for Metabolites
4. `NetGSA_Power_reml.R` - Estimates timing and power for REML on the entire network for BreastCancer and ProstateCancer

## Inputs

Contains the input datasets used in the power simulations

# netgsa_comp

## Code

1. `preprocess_lib.R` - Contains helper functions and data preprocessing functions. Imported in `TCGA_BCa_simulation.R`
2. `TCGA_BCa_simulation.R` - Runs comparison of NetGSA to variety of other methods

## Inputs

Contains the input BreastCancer dataset. This is in a different format and contains more R objects than the corresponding BreastCancer dataset from `netgsa_new`.

