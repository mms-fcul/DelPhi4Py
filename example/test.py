import sys
sys.path.insert(0, '../')

from delphi4py import delphi

#import os

f_crg = 'P.crg'
f_siz = 'DataBaseT.siz'
fpdb = 'P.pdb'

delphimol = delphi(f_crg, f_siz, fpdb, 17250, 121, 4.0895400000,
                   'single', conc=0.1, ibctyp=4, res2=0.01, nlit=500,
                   pbx=True, pby=True, outputfile='LOG_readFiles')

natoms = delphimol.natoms
p_atpos = delphimol.p_atpos

delphimol.runDelPhi(scale_prefocus=1, scale=4, nlit_prefocus=500,
                    nonit=50, nlit=500,
                    acent=[44.014700, 44.014700, 114.781000],
                    nonit_focus=0, relfac_focus=0.0, relpar_focus=0.0,
                    relpar=0.2, relfac=0.2, pbx_focus=False,
                    pby_focus=False, debug=True,
                    outputfile='LOG_runDelPhi')

print delphimol.getSolvation() # float
#print delphimol.getSitePotential() # array

print 'exiting'
