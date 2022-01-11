# -*- coding: utf-8 -*-

# Flexible Net integrating the E. coli metabolic network, citramalate
# synthesis and bioreactor variables. The net is built for a given
# dilution rate, glucose concentration in the medium, and range of
# densities of cells in the tank.

from __future__ import division, print_function
import numpy as np
import cobra.test
from cobra import Reaction, Metabolite
from fnyzer import FNFactory, cobra2fn

def genEcoliCitFN(D = 0.1, l = 5, p = 1000, Xmin = 1.8, Xmax = 1.9):
    #### Parameters
    # D = 0.1/24.0 # (h-1) Dilution rate
    # l = 10 # (g L-1) Glucose concentration in medium
    # Xmin = 3.2 # (gdcw L-1) Minimum density of cells in tank
    # Xmax = 3.4 # (gdcw L-1) Maximum density of cells in tank

    #### Variables and units
    # X: Density of cells in the tank (gdcw L-1)
    # G: Concentration of glucose in the tank (mM)
    # C: Concentration of citramalate in the tank (mM)
    # Fluxes are expressed in concentration of the reactant or product per hour
    
    ### Metabolic network
    Ecolimodelfile = "MODEL1108160000" # E. coli metabolic network
    biomass_reaction = 'Ec_biomass_iJO1366_core_53p95M'
    Ecolicobramodel = cobra.io.read_sbml_model(Ecolimodelfile+'.xml')
    # Ecolicobramodel = pickle.load( open(Ecolimodelfile+'.p', "rb" ) )
    addCimA(Ecolicobramodel) # Add synthesis, transport and exchange of citramalate

    ### Build Flexible Net
    fnet = cobra2fn(Ecolicobramodel)
    fnet['name'] = 'EcoliCitFN' 
    fnet['solver'] = 'cplex'

    ### Bioreactor(tank) concentrations
    fnet['places']['X'] = 1e-3
    fnet['mbounds'] = [str(Xmin)+"<=m['X']", "m['X']<="+str(Xmax)]
    fnet['places']['G'] = 0.0
    fnet['places']['C'] = 0.0
    fnet['places']['P'] = 0.0

    ### Bioreactor reactions (transitions and handlers)
    fnet['shandlers'] = {}

    # Glucose feed
    fnet['trans']['tgin'] = {'l0': D*l, 'a0': 0}
    fnet['vhandlers']['vgin'] = [{'a':('vgin','G'), 'v':('tgin','vgin')},
                                  'a == v']
    fnet['trans']['tpin'] = {'l0': D*p, 'a0': 0}
    fnet['vhandlers']['vpin'] = [{'a':('vpin','P'), 'v':('tpin','vpin')},
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
    fnet['shandlers']['su'] = [{'u':('su',tGlucose_in), 'ut':('su','tut')},
                                'u*'+str(Xmin)+'<= ut', 'ut <= u*'+str(Xmax),
                                'u <= 10'] # This 10 is the upper bound for glucose exchange of the E. coli. Units = [mmol_per_gDW_per_hr]

    # Glucose out of the tank (to effluent)
    fnet['trans']['tgout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vgout'] = [{'a':('G','vgout'), 'v':('tgout','vgout')},
                                  'a == v']
    fnet['shandlers']['sgout'] = [{'g':('G','sgout'), 'r':('sgout','tgout')},
                                   'r == g*'+str(D)]
    
    # Phosphate uptake (from the tank into the cell)
    fnet['trans']['tput'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vput'] = [{'a':('P','vput'), 'v':('tput','vput')},
                                 'a == v']
    phosphate_ex = 'EX_pi_LPAREN_e_RPAREN_' # Reaction ID of glucose exchange in metabolic network
    tPhosphate_in = 't_'+phosphate_ex+'_b'  # b is for backward reaction and f is for forward reaction
    fnet['trans'][tPhosphate_in]['l0'] = 0 # Glucose uptake determined by its intensity handler
    tPhosphate_out = 't_'+phosphate_ex+'_f'
    fnet['trans'][tPhosphate_out]['l0'] = 0 # No glucose out of the cell allowed
    fnet['shandlers']['sup'] = [{'u':('sup',tPhosphate_in), 'ut':('sup','tput')},
                                'u*'+str(Xmin)+'<= ut', 'ut <= u*'+str(Xmax),
                                'u <= 10'] # This 10 is the upper bound for glucose exchange of the E. coli. Units = [mmol_per_gDW_per_hr]

    # Phosphate out of the tank (to effluent)
    fnet['trans']['tpout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vpout'] = [{'a':('P','vpout'), 'v':('tpout','vpout')},
                                  'a == v']
    fnet['shandlers']['spout'] = [{'g':('P','spout'), 'r':('spout','tpout')},
                                   'r == g*'+str(D)]

    # Cell growth
    fnet['trans']['txt'] = {'l0': 0, 'a0': 0}  
    fnet['vhandlers']['vxt'] = [{'a':('vxt','X'), 'v':('txt','vxt')},
                                'a == v']
    biomass_reaction = 'Ec_biomass_iJO1366_core_53p95M' # Reaction ID of biomass in metabolic network
    tBiomass = 't_'+biomass_reaction+ '_f' # transition for biomass
    fnet['trans'][tBiomass]['l0'] = 0 # Biomass production (i.e. growth rate) determined by its intensity handler. It should be equal to D in steady state.
    fnet['shandlers']['sr'] = [{'r':('sr',tBiomass), 'rt':('sr','txt')},
                                'r*'+str(Xmin)+'<= rt', 'rt <= r*'+str(Xmax)]

    # Cells out of the tank (to effluent)
    fnet['trans']['txout'] = {'l0': 0, 'a0': 0}
    fnet['vhandlers']['vxout'] = [{'a':('X','vxout'), 'v':('txout','vxout')},
                                  'a == v']
    fnet['shandlers']['sxout'] = [{'x':('X','sxout'), 'r':('sxout','txout')},
                                   'r == x*'+str(D)]

    # Citramalate production (from cell to tank) ??!! check the arrow in the sketch
    fnet['trans']['tct'] = {'l0': 0, 'a0': 0} 
    fnet['vhandlers']['vct'] = [{'a':('vct','C'), 'v':('tct','vct')},
                                'a == v']
    tExCit = 't_EX_Citramalate_f'
    fnet['trans'][tExCit]['l0'] = 0 # Citramalate production determined by its intensity handler.
    fnet['shandlers']['sz'] = [{'z':('sz',tExCit), 'zt':('sz','tct')},
                                'z*'+str(Xmin)+'<= zt', 'zt <= z*'+str(Xmax)]

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
    fnet['actavplaces'] = ['X', 'G', 'C', 'P']
    fnet['actplaces'] = ['X', 'G', 'C', 'P']
    fnet['options'] = {
            'antype': 'cst',
            'printres': False,
            'printmodel': False,
            'savenet': False,
            'writevars': {
                'avm': ['X', 'G', 'C', 'P'],
                'avl':'all'},
            'plotres': False,
            'writexls': False,
            }

    return fnet


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


EcoliCitFN = genEcoliCitFN(D = 0.1, l = 5, Xmin = 1.8, Xmax = 1.9)
