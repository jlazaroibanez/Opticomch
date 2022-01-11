# Opticomch
The following .py files allow reproducing the simulation results of the 
associated article "Combination of Genome-Scale Models and Bioreactor Dynamics 
to Optimize the Production of Commodity Chemicals".

## Some useful instructions ##

The genome-scale metabolic model used can be downloaded in https://www.ebi.ac.uk/biomodels/MODEL1108160000#Files

The fnyzer User Documentation is located in https://fnyzer.readthedocs.io/en/latest/

It is necessary to install a solver (CPLEX, glpk, Gurobi) in order to compute the numerical values for the 
simulation.

The repository contains two folders:

1. Simul_VitrovsSilico. Includes the code necessary to reproduce the results obtained in the in vitro experiments.

2.  maxProd_optX. Contains the simulation code for the optimization of theoretical productivity given different glucose concentrations, dilution rates and biomass values. It also includes a file (comProdDG_wo.py) that shows the event of wash out when the dilution rate is too high.

