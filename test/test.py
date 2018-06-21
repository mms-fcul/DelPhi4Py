from delphi4py import delphi
#from ctypes import *

#import os

f_crg = 'P.crg'
f_siz = 'DataBaseT.siz'
fpdb = 'P.pdb'

delphimol = delphi(f_crg, f_siz, fpdb, 17250, 121, 4.0895400000,
                   'single', conc=0.1, ibctyp=4, res2=0.01, nlit=500, debug=True)

natoms = delphimol.natoms
p_atpos = delphimol.p_atpos

print delphimol

delphimol.runDelPhi(scale_focus=1.0223857030, nlit_focus=50,
                    nonit_focus=5, relfac_focus=0.2,
                    relpar_focus=0.75, pbx_focus=True, pby_focus=True,
                    acent=[44.014700,44.014700,114.781000], )

#print delphimol.getSolvation() # float
#print delphimol.getSitePotential() # array

print 'exiting'