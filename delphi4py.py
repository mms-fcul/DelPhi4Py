import readFiles
import rundelphi
from ctypes import *


# if using parallel version don't forget to set system-wide variables
# export OMP_NUM_THREADS=8
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/pedror/delphit/dependencies/NanoShaper0.7/build_lib:/home/pedror/delphit/dependencies/cppDelphi77/lib/


class delphi:
    """
    """
    def __init__(self, in_crg, in_siz, in_pdb, natoms, igrid, scale,
                 precision, perfil=0, epsin=2, epsout=80,
                 conc=0, ibctyp=2, res2=0, nlit=0, radprb=1.4,
                 relfac=0.0, relpar=0.0, nonit=0, fcrg=False,
                 pbx=False, pby=False, pbz=False, isurftype=-1,
                 parallel=False, debug=False):
        # TODO #
        # extend number of input parameters to include all/most
        # DelPhi input parameters
        #######################################################
        self.igrid        = int(igrid)
        self.scale        = float(scale)
        self.scale_focus  = None

        self.perfil = int(perfil)

        self.repsin  = float(epsin)
        self.repsout = float(epsout)
        self.conc    = float(conc)
        self.ibctyp  = int(ibctyp)
        self.res2    = float(res2)
        self.nlit    = int(nlit)

        self.acent  = []
        self.natoms = int(natoms)
        self.in_crg = str(in_crg)
        self.in_siz = str(in_siz)
        self.in_pdb = str(in_pdb)

        self.radprb = float(radprb)
        self.energy = ['s', 'c']
        self.site   = ['a', 'q', 'p']
        self.in_crg_len = len(self.in_crg)
        self.in_siz_len = len(self.in_siz)
        self.in_pdb_len = len(self.in_pdb)

        self.rmaxdim = 999.999

        self.relfac = float(relfac)
        self.relpar = float(relpar)
        self.nonit  = int(nonit)
        self.fcrg   = bool(fcrg)
        self.pbx    = bool(pbx)
        self.pby    = bool(pby)
        self.pbz    = bool(pbz)

        self.precision = str(precision)
        self.isurftype = int(isurftype)
        self.parallel  = bool(parallel)

        self.debug = bool(debug)

        if self.precision == 'double':
            self.float_type = c_double
        elif self.precision == 'single':
            self.float_type = c_float
        else:
            raise IOError('Unknown precision definition {0}. '
                          'It should be either "double" or "single"'
                          .format(self.precision))

        # -1 NanoShaper off
        # 0 connolly surface
        # 1 skin
        # 2 blobby
        # 3 mesh
        # 4 msms
        # only tested -1 and 0
        if self.isurftype not in (0, -1):
            raise IOError('Unknown precision definition {0}. '
                          'It should be either 0 to activate or -1 to deactivate'
                          .format(self.isurftype))

        # TODO #
        # check other input parameters
        #######################################################

        self.resetDelPhiData()
        self.readFiles()

        if self.perfil != 0:
            self.igrid = int(self.scale * 100 / self.perfil * self.rmaxdim)
            self.igrid_focus = self.igrid

        if self.igrid % 2 == 0:
            self.igrid += 1

        self.esolvation = 999.999

    def resetDelPhiData(self):
        """Resets all DelPhi Input data structures"""
        # internal DelPhi DataStructures
        # defined as c_type arrays and
        # passed to DelPhi as pointers
        self.atpos   = self.float_type * 3 * self.natoms
        self.p_atpos    = self.atpos()
        self.i_atpos = addressof(self.p_atpos)

        self.rad3   = self.float_type * self.natoms
        self.p_rad3    = self.rad3()
        self.i_rad3 = addressof(self.p_rad3)

        self.chrgv4   = self.float_type * self.natoms
        self.p_chrgv4    = self.chrgv4()
        self.i_chrgv4 = addressof(self.p_chrgv4)

        self.atinf   = (c_char * 15  * self.natoms)()
        self.i_atinf = addressof(self.atinf)

        self.nmedia = 1
        self.nobject = 1
        self.len_medeps = self.nmedia + self.nobject
        self.medeps  = self.float_type * self.len_medeps
        self.p_medeps    = self.medeps()
        self.i_medeps = addressof(self.p_medeps)

        self.len_iatmed  = self.natoms + 1
        self.iatmed   = c_int * self.len_iatmed
        self.p_iatmed    = self.iatmed()
        self.i_iatmmed = addressof(self.p_iatmed)

        self.dataobject   = (c_char * 96 * self.nobject * 2)()
        self.i_dataobject = addressof(self.dataobject)

    def readFiles(self):
        """ """
        self.rmaxdim = readFiles.delphi(self.igrid, self.scale,
                                        self.repsin, self.repsout,
                                        self.acent, self.in_pdb,
                                        self.in_crg, self.in_siz,
                                        self.natoms, self.nobject,
                                        self.i_atpos, self.i_rad3,
                                        self.i_chrgv4, self.i_atinf,
                                        self.i_medeps, self.i_iatmmed,
                                        self.i_dataobject,
                                        self.rmaxdim)

        if self.debug:
            print '    x        y        z     radius  charge       atinf'
            for i in range(self.natoms):
                print ('{0:8.3f} {1:8.3f} {2:8.3f} {3:7.3f} {4:7.3f} {5}'
                       .format(self.p_atpos[i][0], self.p_atpos[i][1],
                               self.p_atpos[i][2], self.p_rad3[i],
                               self.p_chrgv4[i],   self.atinf[i].value))                

    def runDelPhi(self, scale_focus=None, nlit_focus=None, acent=None,
                  nonit_focus=None, relfac_focus=None, relpar_focus=None,
                  pbx_focus=None, pby_focus=None, focusing=False):
        """
        """
        if acent:
            self.acent = acent
        focusing_trigger = False
        if scale_focus:
            scale  = float(scale_focus)
            scale_focus = scale
            self.scale_focus = scale
            ibctyp = self.ibctyp

            if nlit_focus:
                nlit = nlit_focus
            else:
                nlit = self.nlit

            nonit  = nonit_focus
            relfac = relfac_focus   
            relpar = relpar_focus
            pbx    = pbx_focus
            pby    = pby_focus

            in_frc   = 'self'            
            out_phi  = True

            self.sitpot   = self.float_type * self.natoms
            self.p_sitpot = self.sitpot()
            self.i_sitpot = addressof(self.p_sitpot)        

            self.len_phimap = self.igrid * self.igrid * self.igrid
            self.phimap4    = c_float * self.len_phimap
            self.p_phimap4  = self.phimap4()
            self.i_phimap4  = addressof(self.p_phimap4)

            focusing = True

        elif focusing:
            scale = self.scale
            scale_focus = self.scale_focus                
            ibctyp = 3                
            nlit = self.nlit               

            in_frc   = 'self'            
            out_phi   = False

            focusing = False

            nonit  = self.nonit
            relfac = self.relfac   
            relpar = self.relpar
            pbx    = self.pbx
            pby    = self.pby

        else:
            scale = self.scale
            ibctyp = self.ibctyp                
            nlit = self.nlit               
            
            in_frc    = ''            
            out_phi   = False

            nonit  = self.nonit
            relfac = self.relfac   
            relpar = self.relpar
            pbx    = self.pbx
            pby    = self.pby

            self.sitpot    = self.float_type * self.natoms
            self.p_sitpot  = self.sitpot()
            self.i_sitpot  = addressof(self.p_sitpot)

            self.len_phimap = 0
            self.phimap4    = c_float * self.len_phimap
            self.p_phimap4  = self.phimap4()
            self.i_phimap4  = addressof(self.p_phimap4)

        self.esolvation = rundelphi.delphi(self.igrid, scale,
                                           self.repsin, self.repsout,
                                           self.radprb, self.conc,
                                           ibctyp, self.res2, nlit,
                                           self.acent, self.energy,
                                           self.site, nonit, relfac,
                                           relpar, pbx, pby, in_frc,
                                           self.natoms, self.nmedia,
                                           self.nobject, self.i_atpos,
                                           self.i_rad3, self.i_chrgv4,
                                           self.i_atinf,
                                           self.i_medeps,
                                           self.i_iatmmed,
                                           self.i_dataobject,
                                           self.i_phimap4,
                                           scale_focus, out_phi,
                                           self.i_sitpot,
                                           self.esolvation,
                                           self.isurftype,
                                           self.parallel)

        if focusing:
            self.runDelPhi(focusing=True)

        if self.debug:
            #for i in range(self.natoms):
            #    if self.p_sitpot[i] != 0.0:
            #        print self.p_atpos[i][0], self.p_atpos[i][1], self.p_atpos[i][2], self.p_sitpot[i]
            print self.esolvation
            print '    x        y        z     radius  charge       atinf'
            for i in range(self.natoms):
                print ('{0:8.3f} {1:8.3f} {2:8.3f} {3:7.3f} {4:7.3f} ' \
                       '"{5}" {6:10.3f}'.format(self.p_atpos[i][0],
                                                self.p_atpos[i][1],
                                                self.p_atpos[i][2],
                                                self.p_rad3[i],
                                                self.p_chrgv4[i],
                                                self.atinf[i].value,
                                                self.p_sitpot[i]))

    def getSolvation(self):
        return self.esolvation

    def getSitePotential(self):
        """Returns site potential as a python list

        Warning: inefficient
        """
        p_sitpot_list = []
        for i in self.p_sitpot:
            p_sitpot_list.append(i)
            
        return p_sitpot_list
            

    def __str__(self):
        """Outputs the parameters used"""
        out = """
        DelPhi Parameters

        igrid       = {}
        scale       = {}
        perfil      = {}

        repsin  = {}
        repsout = {}
        conc    = {}
        ibctyp  = {}
        res2    = {}
        nlit    = {}
        
        acent  = {}
        natoms = {}
        in_crg = {}
        in_siz = {}
        in_pdb = {}
        
        radprb  = {}
        energy  = {}
        site    = {}
        in_crg_len = {}
        in_siz_len = {}
        in_pdb_len = {}

        relfac = {}
        relpar = {}
        nonit  = {}
        fcrg   = {}
        pbx    = {}
        pby    = {}
        pbz    = {}

        precision  = {}
        float_type = {}
        isurftype  = {}
        parallel   = {}
        debug      = {}

        """.format(self.igrid, self.scale, self.perfil, self.repsin,
        self.repsout, self.conc, self.ibctyp, self.res2, self.nlit,
        self.acent, self.natoms, self.in_crg, self.in_siz,
        self.in_pdb, self.radprb, self.energy, self.site,
        self.in_crg_len, self.in_siz_len, self.in_pdb_len,
        self.relfac, self.relpar, self.nonit, self.fcrg, self.pbx,
        self.pby, self.pbz, self.precision, self.float_type,
        self.isurftype, self.parallel, self.debug)
        return out

