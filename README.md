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

ClaDS_summary_all_trees.csv

4. Script for analysing imbalance of trees simulated under each combination of diversification and preservation regime
5. Script for estimating diversification parameters of empirical trees with CladDS and then simulating trees with those parameters under every combination of diversification model and preservation regime
6. Script for analysing the imbalance of trees simulated from parameter values of empirical trees estimate with CladDS
7. Script for making lolipop graph of empirical trees

CI_simulation_output_ne100.csv contains the input parameters and output values of the trees simulated under every combination of diversification model and preservation regime. 

ClaDS_summary_all_trees.csv includes the input parameters and output values of the trees simulated from CladDS estimates for empirical trees under every combination of diversification model and preservation regime.

The best_model_sim.txt includes the model summary for the glmm fit to the dataset generated from the general simulation procedure and the [study]_best_model.txt files are model sumaries for each of the glmm models fit to the dataset generated from the informed simulation procedure. 
