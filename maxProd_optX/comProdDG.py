#!/usr/bin/python
# Transform cobra model to Flexible Net (FN)
from __future__ import division, print_function
import numpy as np
import cobra.test
from EcoliBiorFN import comProductivity, loadEcolimodel
import xlwt

# Creating xls file

row = 0
xlsfile = "Results_var_D_X_glc_ok.xls"
wb = xlwt.Workbook()
bstyle = xlwt.easyxf('font: bold on')
ws = wb.add_sheet("Optimal Biomass")
ws2 = wb.add_sheet("Maximal productivity")
ws.write(row, 0, 'D (h-1)', bstyle)
ws.write(row, 1, 'glc_conc (gL-1)', bstyle)
ws.write(row, 2, 'optX (gdcw L-1)', bstyle)
ws2.write(row, 0, 'D (h-1)', bstyle)
ws2.write(row, 1, 'glc_conc (gL-1)', bstyle)
ws2.write(row, 2, 'maxProd (h-1)', bstyle)
wb.save(xlsfile)

#fnet = loadEcolimodel(filename = "MODEL1108160000")
#maxProd, optX = comProductivity(fnet, D=0.23, glc=5.0, X0 = 0.70, Xf = 1.0, XNsamps = 1, Xint = 0.02)
#
#print("maxProd:", maxProd, "optX:", optX)

# Parameters
Dmin, Dmax, DNsamps = 0.01, 0.81, 9  # (h-1) Minimal, maximal and number of dilution rates 
glcmin, glcmax, glcNsamps = 1.0, 11.0, 6  # (g L-1) Minimal, maximal and number of glucose concentrations in medium

X0, Xf, XNsamps = 0.0, 5.5, 56 # Parameters to the cell densities to which the cell density will be constrained

Xint = (Xf - X0)/(XNsamps - 1)
                                         # See function comProductivity() for details

# Load model
fnet = loadEcolimodel(filename = "MODEL1108160000")

# Optimize the net for each dilution rate and each glucose concetration in the medium
for D in np.linspace(Dmin, Dmax, DNsamps):
    for glc in np.linspace(glcmin, glcmax, glcNsamps):
       row = row + 1
       print("D:", f"{D:.4f}", "h-1; glc:", f"{glc:.4f}","g L-1")
       maxProd, optX = comProductivity(fnet, D, glc, X0 = X0, Xf = Xf, XNsamps = XNsamps, Xint = Xint)
       ws.write(row, 0, D)
       ws.write(row, 1, glc)
       ws.write(row, 2, optX)
       ws2.write(row, 0, D)
       ws2.write(row, 1, glc)
       ws2.write(row, 2, maxProd)
       
       if (maxProd<0):
           print("Infeasible conditions")
       else:
           print("  maxProd:", f"{maxProd:.4f}", "h-1; optX:", f"{optX:.4f}", "gdcw L-1")
wb.save(xlsfile)            
