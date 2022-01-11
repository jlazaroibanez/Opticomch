from __future__ import division, print_function
import numpy as np
import cobra.test
from fnyzer import FNFactory
# from EcoliCitFNgen import genEcoliCitFN
from EcoliBiorFN import genEcoliBiorFN

def main():
    ## Parameters
    l = 27.754 # (mM) Glucose concentration in medium
    Ds = [0.1, 0.1] # (h-1) Dilution rates to be tested
    Xins = [(1.88, 1.88)] # (gdcw L-1) (minimum, maximum) densities of cells in tank to be tested for each D in Ds
    for D in Ds:
        for Xin in Xins:
            Xmin, Xmax = Xin[0], Xin[1]
            print("Parameters. glcmed(mM)=", l, "; D(h-1)=", D, "; avm[X](gdcw L-1) in [", Xmin, ",", Xmax,"]")
            print("Building net...")
            EColiCitFN = genEcoliBiorFN(l = l, D = D, Xmin = Xmin, Xmax
            = Xmax)
            netobj = FNFactory(EColiCitFN) # Build net object
            print("Optimizing...")
            try:
                netobj.optimize()
                print("Solution:", netobj.objval)
                print("avm[G](mM):", netobj.places['G'].avm)
                print("avm[X](gdcw L-1):", netobj.places['X'].avm)


                print("avm[C](mM):", netobj.places['C'].avm)
                print("avl[tgout](mM h-1):", netobj.trans['tgout'].avl)
                print("avl[txout](gdcw L-1 h-1):", netobj.trans['txout'].avl)
                print("avl[tcout](mM h-1):", netobj.trans['tcout'].avl)
                growth_reaction = 'Ec_biomass_iJO1366_core_53p95M' # Reaction ID of biomass in metabolic network
                tgrowth = 't_'+growth_reaction+ '_f' # transition for biomass
                print("avl[growth]:", netobj.trans[tgrowth].avl)
            except ValueError as inst:
                print("Problem not feasible")
            print("------------------------------")

if __name__ == '__main__':
    main()
