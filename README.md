# Multiple diversification and fossil sampling models can make total-evidence phylogenies more imbalanced but only a few make them look like empirical trees

This repository includes R code for replicating the analyses and figures presented in this paper. It also includes several datasets summarising the outputs of these analyses. 

## Description of the data and file structures
R Code
1) functions.R
    A series of helper functions including tree building functions for each of the diversification models and sampling functions for each of the preservation regimes
2) simultating_trees_e100.R
   Script for simulating trees under each combination of diversification model and fossil sampling model
3) analysis_simulation.R
   Script for analysing the trees simulated using simultating_trees_e100.R
4) ClaDS_simulation.R
   Script for fitting ClaDS models to empirical total-evidence phylogenies
5) CLaDS_empirical_simulation.R
   Script for simulating trees from the parameter values estimated using ClaDS_simulation.R
6) lolipop.R
   Script for reproducing Figure 5.

Supplementary Data
CI_simulation_output_ne100.csv
    Summary of the output of the trees simulated under each combination of fiversification models and fossil sampling models produced using simultating_trees_e100.R and analysis_simulation.R. The dataset includes the number of fossil lineages on the simulated tree (Fossil_no), the number of fossil preserved by the fossil sampling models (Preservation), the length of the simulated tree (Tree_length), the starting speciation rate (lambda_init), the starting extinction rate (mu_init), the lower and upper bounds of the distributions parameterising the shifts in speciation and extinction rates (min_*_shift and max_*_shift), the number of extant species on the tree (n_extant), the diversification model under which the tree was simulated (model, see manuscript for a description of each model), the point in time in which the mass extinction event occurs (extinction_date, only relevant for models which include a mass extinction), the proportion of lineages which survive the extinction event (survivorship), the length of the root edge branch (root.edge, used to correct calculations of Tree_length for some of the diversification models), the point in time at which a universal shift in speciation rate occurs (shift_date, only relevant to the universal diversification rate shift model), the scaler for the universal speciation rate shift (shift_scaler), the limit of the size of the clade that arises from the key innovation (n_sub, only applicable to the density dependent with key innovation model), the point in time at which the key innovation arises (tinn), the fossil sampling model used to sample the extinct lineages (preservation), the Colless Index of imbalance for the simulated phylogeny (Index)

    
ClaDS_summary_all_trees.csv

4. Script for analysing imbalance of trees simulated under each combination of diversification and preservation regime
5. Script for estimating diversification parameters of empirical trees with CladDS and then simulating trees with those parameters under every combination of diversification model and preservation regime
6. Script for analysing the imbalance of trees simulated from parameter values of empirical trees estimate with CladDS
7. Script for making lolipop graph of empirical trees

CI_simulation_output_ne100.csv contains the input parameters and output values of the trees simulated under every combination of diversification model and preservation regime. 

ClaDS_summary_all_trees.csv includes the input parameters and output values of the trees simulated from CladDS estimates for empirical trees under every combination of diversification model and preservation regime.

The best_model_sim.txt includes the model summary for the glmm fit to the dataset generated from the general simulation procedure and the [study]_best_model.txt files are model sumaries for each of the glmm models fit to the dataset generated from the informed simulation procedure. 
