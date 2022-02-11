# opticomch
The provided .py file can be used to reproduce the simulation results of the 
paper: "Combination of Genome-Scale Models and Bioreactor Dynamics 
to Optimize the Production of Commodity Chemicals".

## Installation and use ##

- Download the file EcoliBiorFN.py from this repository.

- Download the E. coli genome-scale metabolic model MODEL1108160000 from:
https://www.ebi.ac.uk/biomodels/MODEL1108160000

- Install fnyzer (detailed documentation at https://fnyzer.readthedocs.io):

        $ pip install fnyzer

- Install a solver (GLPK, CPLEX, Gurobi) in order to compute the numerical values for the
simulation. The name of the solver (glpk, cplex, gurobi) must be set in the parameter solver of the
function loadEcolimodel(). As an example, GLPK can be easily installed as a Debian package,
e.g. in Ubuntu:

        $ apt-get install glpk-utils

- Open a Python interpreter and execute the following lines to compute the theoretical maximum
productivity for a given dilution rate, e.g. D = 0.41 h-1, and a given glucose concentration in the
medium, e.g. glc = 3.0 g L-1 (the lines below explore cell densities in the interval [0, 1.2] gdcw L-
1):

        from EcoliBiorFN import comProductivity, loadEcolimodel
    
        fnet = loadEcolimodel(filename = "MODEL1108160000_url") # Use the
        appropriate string name of the model file you downloaded. This line is for the file
        MODEL1108160000_url.xml
    
        maxProd, optX, res_glucose = comProductivity(fnet, D = 0.41, glc = 3.0, X0 = 0.0, Xf = 1.2,
        XNsamps = 13, Xint = 0.1)
    
The variables maxProd and optX contain the maximum productivity and associated cell density.
