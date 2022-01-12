# Opticomch
The following .py files allow reproducing the simulation results of the 
associated article "Combination of Genome-Scale Models and Bioreactor Dynamics 
to Optimize the Production of Commodity Chemicals".

## Some useful instructions ##

The genome-scale metabolic model used can be downloaded in https://www.ebi.ac.uk/biomodels/MODEL1108160000#Files

The fnyzer User Documentation is located in https://fnyzer.readthedocs.io/en/latest/
The fnyzer tool can be installed from the command line in a Linux systam with: $ pip install fnyzer

It is necessary to install a solver (CPLEX, glpk, Gurobi) in order to compute the numerical values for the 
simulation.

Installing glpk solver: $ apt-get install glpk-utils

The repository the simulation code that creates the flexible net and allows the optimization of theoretical productivity given different a glucose concentration and dilution rate which iterates over different biomass values. 
