# -*- coding: utf-8 -*-

# Functions to compute the maximum productivity of a
# flexible Net integrating the E. coli metabolic network,
# citramalate synthesis and bioreactor variables.

from __future__ import division, print_function
import numpy as np
import cobra
from cobra import Reaction, Metabolite
from fnyzer import FNFactory, cobra2fn

# Constants
glucose_mmass = 180.1577 # Glucose molar mass: 180.1577 g/mol
citra_mmass = 148.11 # Citramalic acid molar mass: 148.11 g/mol

def comProductivity(fnet, D = 0.1, glc = 5.0, X0 = 0.5, Xf = 5.0, XNsamps = 4, Xint = 0.2):
    # Computes the maximum productivity (h-1) and the optimal cell density (gdcw L-1) 
    # for which it is obtained for a given dilution rate D and glucose concentration glc
    # It returns -1 if the problem is infeasible for the given D and glc
    
    #### Parameters
    # D = 0.1/24.0 # (h-1) Dilution rate
    # glc = 5.0 # (g L-1) Glucose concentration in medium in g L-1
    # Several FNs will be optimized, each one with a different constraint for X
    # X0 = 0.5 # (gdcw L-1) Minimum X to be tested
    # Xf = 20.0 # (gdcw L-1) Maximum X to be tested
    # XNsamps = 4 # Numer of FNs to be tested
    # Xint = 0.2 # (gdcw L-1) Length of the interval in which X will be constrained
    # E.g. with the above parameters, FNs with the following intervals constraining X will be
    # optimized: [0.5, 0.7], [7.0, 7.2], [13.5, 13.7], [20.0, 20.2], 
    
    glcmM = glc*1000/glucose_mmass # Glucose concentration in medium in mM
    maxProd, optX, res_glucose = -1, -1, 1000
    for X in np.linspace(X0, Xf, XNsamps):
        genEcoliBiorFN(fnet, D = D, glcmM = glcmM, Xmin = X, Xmax = X+Xint)
        netobj = FNFactory(fnet) # Build net object
        try:
            netobj.optimize()
            prod = netobj.trans['tcout'].avl*(citra_mmass/1000)/glc
            glc_tank = netobj.places['G'].avm*(glucose_mmass/1000)
            # The above expression for productivity is equivalent to (see Section "Optimizing the productivity" of the paper):
            # mu * Yield(P/S) = D * (netobj.trans['tcout'].avl*(citra_mmass/1000))/(netobj.trans['tgin'].avl*(glucose_mmass/1000))
            print("  Productivity with X in [", f"{X:.2f}",",", f"{X+Xint:.2f}", "] gdcw L-1:", f"{prod:.4f}", "h-1")            
            if (prod > maxProd):
                maxProd = prod
                optX = netobj.places['X'].avm # Cell density for which the productivity is maximum
            if (glc_tank < res_glucose):
                res_glucose = glc_tank
            # print(netobj.trans['tcout'].avl*(citra_mmass/1000))
            print(netobj.trans['tcout'].avl*(citra_mmass/1000)/glc)
            print(glc_tank)
            #print(D * (netobj.trans['tcout'].avl*(citra_mmass/1000))/(netobj.trans['tgin'].avl*(glucose_mmass/1000)))
                        
        except:
            print("  Problem infeasible with X in [", f"{X:.2f}",",", f"{X+Xint:.2f}", "] gdcw L-1")
    
    return maxProd, optX, res_glucose


def loadEcolimodel(filename = "MODEL1108160000", name = "EcolicitFN", solver = "cplex"):
    #### Parameters
    # filename: file with the SBML model
    # name: name of the FN
    # solver: solver to be used
    biomass_reaction = 'Ec_biomass_iJO1366_core_53p95M'
    Ecolicobramodel = cobra.io.read_sbml_model(filename+'.xml')
    addCimA(Ecolicobramodel) # Add synthesis, transport and exchange of citramalate
    fnet = cobra2fn(Ecolicobramodel) # Build Flexible Net
    fnet['name'] =  name
    fnet['solver'] = solver
    return fnet


def genEcoliBiorFN(fnet, D = 0.1, glcmM = 27, Xmin = 1.8, Xmax = 1.9):
    # Adds to fnet the bioreactor dynamics
    #### Parameters
    # D = 0.1/24.0 # (h-1) Dilution rate
    # glcmM = 5 # (mM) Glucose concentration in medium in mM
    # Xmin = 1.8 # (gdcw L-1) Minimum density of cells in tank
    # Xmax = 1.9 # (gdcw L-1) Maximum density of cells in tank

    #### Variables and units
    # X: Density of cells in the tank (gdcw L-1)
    # G: Concentration of glucose in the tank (mM)
    # C: Concentration of citramalate in the tank (mM)
    # Fluxes are expressed in concentration of the reactant or product per hour
    
    ### Bioreactor(tank) concentrations
    fnet['places']['X'] = 1e-3
    fnet['mbounds'] = [str(Xmin)+"<=m['X']", "m['X']<="+str(Xmax)]
    fnet['places']['G'] = 0.0
    fnet['places']['C'] = 0.0

    ### Bioreactor reactions (transitions and handlers)
    fnet['shandlers'] = {}

    # Glucose feed
    fnet['trans']['tgin'] = {'l0': D*glcmM, 'a0': 0}
    fnet['vhandlers']['vgin'] = [{'a':('vgin','G'), 'v':('tgin','vgin')},
                                  'a == v']

    # Glucose uptake (from the tank into the cell)
    fnet['trans']['tut'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vut'] = [{'a':('G','vut'), 'v':('tut','vut')},
                                 'a == v']
    glucose_ex = 'EX_glc_LPAREN_e_RPAREN_' # Reaction ID of glucose exchange in metabolic network
    tGlucose_in = 't_'+glucose_ex+'_b'  # b is for backward reaction and f is for forward reaction
    fnet['trans'][tGlucose_in]['l0'] = 0 # Glucose uptake determined by its intensity handler
    tGlucose_out = 't_'+glucose_ex+'_f'
    fnet['trans'][tGlucose_out]['l0'] = 0 # No glucose out of the cell allowed
    fnet['shandlers']['hu'] = [{'ug':('hu',tGlucose_in), 'ut':('hu','tut')},
                                'ug*'+str(Xmin)+'<= ut', 'ut <= ug*'+str(Xmax),
                                'ug <= 10'] # This 10 is the upper bound for glucose exchange of the E. coli. Units = [mmol_per_gDW_per_hr]

    # Glucose out of the tank (to effluent)
    fnet['trans']['tgout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vgout'] = [{'a':('G','vgout'), 'v':('tgout','vgout')},
                                  'a == v']
    fnet['shandlers']['sgout'] = [{'g':('G','sgout'), 'r':('sgout','tgout')},
                                   'r == g*'+str(D)]

    # Cell growth
    fnet['trans']['txt'] = {'l0': 0, 'a0': 0}  
    fnet['vhandlers']['vxt'] = [{'a':('vxt','X'), 'v':('txt','vxt')},
                                'a == v']
    biomass_reaction = 'Ec_biomass_iJO1366_core_53p95M' # Reaction ID of biomass in metabolic network
    tBiomass = 't_'+biomass_reaction+ '_f' # transition for biomass
    fnet['trans'][tBiomass]['l0'] = 0 # Biomass production (i.e. growth rate) determined by its intensity handler. It should be equal to D in steady state.
    fnet['shandlers']['hr'] = [{'r':('hr',tBiomass), 'rt':('hr','txt')},
                                'r*'+str(Xmin)+'<= rt', 'rt <= r*'+str(Xmax)]

    # Cells out of the tank (to effluent)
    fnet['trans']['txout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vxout'] = [{'a':('X','vxout'), 'v':('txout','vxout')},
                                  'a == v']
    fnet['shandlers']['sxout'] = [{'x':('X','sxout'), 'r':('sxout','txout')},
                                   'r == x*'+str(D)]

    # Citramalate production (from cell to tank)
    
    fnet['trans']['tct'] = {'l0': 0, 'a0': 0} 
    fnet['vhandlers']['vct'] = [{'a':('vct','C'), 'v':('tct','vct')},
                                'a == v']
    tExCit = 't_EX_Citramalate_f'
    fnet['trans'][tExCit]['l0'] = 0 # Citramalate production determined by its intensity handler.
    fnet['shandlers']['hc'] = [{'z':('hc',tExCit), 'ct':('hc','tct')},
                                'z*'+str(Xmin)+'<= ct', 'ct <= z*'+str(Xmax)]

    # Citramalate out of the tank (to effluent)
    fnet['trans']['tcout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vcout'] = [{'a':('C','vcout'), 'v':('tcout','vcout')},
                                  'a == v']
    fnet['shandlers']['scout'] = [{'c':('C','scout'), 'r':('scout','tcout')},
                                   'r == c*'+str(D)]    

    # Other net options
    fnet['obj'] = {'f': "avl['tct']", 'sense': 'max'}
    fnet['exavtrans'] = 'all'
    fnet['extrans'] = 'all'
    fnet['actavplaces'] = ['X', 'G', 'C']
    fnet['actplaces'] = ['X', 'G', 'C']
    fnet['options'] = {
            'antype': 'cst',
            'savenet': False,            
            'printres': False,
            'printmodel': False,
            'writevars': {
                'avm': ['X', 'G', 'C'],
                'avl':'all'},
            'plotres': False,
            'writexls': True,
            }



#######################################################################
def addCimA(model):
    """Add CimA reaction and sink for citramalate to cobra model"""
    reaccima = Reaction('CIMA')
    reaccima.name = '(R)-Citramalate production'
    reaccima.lower_bound = 0.0
    reaccima.upper_bound = 1000.0
    # reaccima.objective_coefficient = 0.0
    
    """CIMA reaction"""
    pyr_c = model.metabolites.get_by_id("pyr_c") # Pyruvate
    accoa_c = model.metabolites.get_by_id("accoa_c") # Acetyl-CoA
    h2o_c = model.metabolites.get_by_id("h2o_c") # H2O
    citramalate_c = Metabolite(
        'citramalate_c',
        formula='C5H6O5',
        name='(R)-citramalate',
        charge=-2,
        compartment='c')
    coa_c = model.metabolites.get_by_id("coa_c") # CoA
    h_c = model.metabolites.get_by_id("h_c") # H+
    
    reaccima.add_metabolites({pyr_c: -1.0,
                              accoa_c: -1.0,
                              h2o_c: -1.0,
                              citramalate_c: 1.0,
                              coa_c: 1.0,
                              h_c: 1.0})
    reaccima.gene_reaction_rule = 'CimA37'     
#    print(reaccima.reaction)                          
#    print(reaccima.genes)                          
    model.add_reaction(reaccima)
    
    
    """Transport1 for Citramalate"""
    reactransp_1 = Reaction('CitraTransp1')
    reactransp_1.name = 'Citramalate Transport from cytoplasm to periplasm'
    reactransp_1.lower_bound = 0.0
    reactransp_1.upper_bound = 1000.0    
    # reaccisink.objective_coefficient = 0.0    
    
    citramalate_c = model.metabolites.get_by_id("citramalate_c") # Citramalate
    citramalate_p = Metabolite(
        'citramalate_p',
        formula='C5H6O5',
        name='(R)-citramalate_p',
        charge=-2,
        compartment='p')
    
    reactransp_1.add_metabolites({citramalate_c: -1.0,
                                citramalate_p: 1.0})
#    print(reacTransp.reaction)                          
#    print(reacTransp.genes)                          
    model.add_reaction(reactransp_1)
    
    """Transport2 for Citramalate"""
    reactransp_2 = Reaction('CitraTransp2')
    reactransp_2.name = 'Citramalate Transport from periplasm to the extracellular space'
    reactransp_2.lower_bound = -1000.0
    reactransp_2.upper_bound = 1000.0    
    # reaccisink.objective_coefficient = 0.0    
    
    citramalate_p = model.metabolites.get_by_id("citramalate_p") # Citramalate
    citramalate_e = Metabolite(
        'citramalate_e',
        formula='C5H6O5',
        name='(R)-citramalate_e',
        charge=-2,
        compartment='e')
    
    reactransp_2.add_metabolites({citramalate_e: -1.0,
                                citramalate_p: 1.0})
#    print(reacTransp.reaction)                          
#    print(reacTransp.genes)                          
    model.add_reaction(reactransp_2)
    
    """Exchange reaction for Citramalate"""    
    reaccEX = Reaction('EX_Citramalate')
    reaccEX.name = 'Exchange reaction to allow (R)-Citramalate to leave the system'
    reaccEX.lower_bound = 0.0
    reaccEX.upper_bound = 1000.0
    # reaccisink.objective_coefficient = 0.0
    
    reaccEX.add_metabolites({citramalate_e: -1.0})
#    print(reaccisink.reaction)                          
#    print(reaccisink.genes)                          
    model.add_reaction(reaccEX)


