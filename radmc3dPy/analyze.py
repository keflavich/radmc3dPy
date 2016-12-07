"""
PYTHON module for RADMC3D 
(c) Attila Juhasz 2011,2012,2013,2014

This sub-module contains classes and functions to read and write input/output data to/from RADMC3D

CLASSES:
--------
   radmc3dData
   radmc3dDustOpac
   radmc3dGrid
   radmc3dStars

FUNCTIONS:
----------
    read_data()
    readGrid()
    readMasterOpac()
    writeMasterOpac()
    readOpac()
    readParams()

"""
try:
    from numpy import *
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'


try:
    from matplotlib.pylab import *
except:
    print ' WARNING'
    print ' matploblib.pylab cannot be imported ' 
    print ' To used the visualization functionality of the python module of RADMC-3D you need to install matplotlib'
    print ' Without matplotlib you can use the python module to set up a model but you will not be able to plot things or'
    print ' display images'

from subprocess import Popen
import sys, os, copy
from radmc3dPy.natconst import *
from radmc3dPy.crd_trans import vrot, ctrans_sph2cart





class radmc3dGrid():
    """
    Class for the spatial and frequency grid used by RADMC3D

    ATTRIBUTES:
    -----------
        crd_sys    - 'car'/'cyl'/'sph' coordinate system of the spatial grid
        act_dim    - A three element vector the i-th element is 1 if the i-th dimension is active, otherwize the i-th element is zero
        nx         - Number of grid points in the x (cartesian) / r (cylindrical) / r (spherical) dimension
        ny         - Number of grid points in the y (cartesian) / theta (cylindrical) / theta (spherical) dimension
        nz         - Number of grid points in the z (cartesian) / z (cylindrical) / phi (spherical) dimension
        nxi        - Number of cell interfaces in the x (cartesian) / r (cylindrical) / r (spherical) dimension
        nyi        - Number of cell interfaces in the y (cartesian) / theta (cylindrical) / theta (spherical) dimension
        nzi        - Number of cell interfaces in the z (cartesian) / z (cylindrical) / phi (spherical) dimension
        nwav       - Number of wavelengths in the wavelength grid
        freq       - Number of frequencies in the grid (equal to nwav)
        x          - Cell centered x (cartesian) / r (cylindrical) / r (spherical)  grid points
        y          - Cell centered y (cartesian) / theta (cylindrical) / theta (spherical)  grid points
        z          - Cell centered z (cartesian) / z (cylindrical) / phi (spherical)  grid points
        xi         - Cell interfaces in the x (cartesian) / r (cylindrical) / r (spherical)  dimension
        yi         - Cell interfaces in the y (cartesian) / theta (cylindrical) / theta (spherical)  dimension
        zi         - Cell interfaces in the z (cartesian) / z (cylindrical) / phi (spherical)  dimension
        wav        - Wavelengh  grid
        freq       - Frequency  grid

    METHODS:
    --------
    """
# --------------------------------------------------------------------------------------------------

    def __init__(self):

        self.crd_sys = 'sph'
        self.act_dim = [1,1,1]
        self.nx    = -1
        self.ny    = -1
        self.nz    = -1
        self.nxi   = -1
        self.nyi   = -1
        self.nzi   = -1
        self.nwav  = -1
        self.nfreq = -1
        self.x     = -1
        self.y     = -1
        self.z     = -1
        self.xi    = -1
        self.yi    = -1
        self.zi    = -1
        self.wav   = -1
        self.freq  = -1

# --------------------------------------------------------------------------------------------------
    def makeWavelengthGrid(self, wbound=None, nw=None, ppar=None):
        """
        Function to create a wavelength/frequency grid 

        INPUT:
        ------
            wbound : list of at least two elements containing the wavelength boundaries of the wavelength grid
            nw     : list of len(wbound)-1 elements containing the number of wavelengths between the bounds
                     set by wbound
        OPTIONS:
        --------
            ppar   : parameter dictionary 
        """
        
        if ppar:
            if not wbound: wbound = ppar['wbound']
            if not nw: nw = ppar['nw']

        if (wbound==None)|(nw==None):
            if (ppar==None): 
                print 'ERROR!'
                print 'Either the boundaries or the number of gridpoints has not be specified in the wavelength grid'
                return
            
        self.nwav = nw[0]
        self.wav  = wbound[0] * (wbound[1]/wbound[0])**(arange(nw[0], dtype=float64) / nw[0])

        for ipart in range(1,len(nw)-1): 
            dum      = wbound[ipart] * (wbound[ipart+1]/wbound[ipart])**(arange(nw[ipart], dtype=float64) / nw[ipart])
            self.wav = append(self.wav, dum)

        ipart      = len(nw)-1
        dum        = wbound[ipart] * (wbound[ipart+1]/wbound[ipart])**(arange(nw[ipart], dtype=float64) / (nw[ipart]-1.))
        self.wav   = append(self.wav, dum)
        self.nwav  = self.wav.shape[0]
        self.freq  = cc / self.wav
        self.nfreq = self.nwav

# --------------------------------------------------------------------------------------------------
    def writeWavelengthGrid(self, fname=''):
        """
        Function to write the wavelength grid to a file (e.g. wavelength_micron.inp)

        OPTIONS:
        --------
            fname  - File name into which the wavelength grid should be written. If omitted 'wavelength_micron.inp' will be used
        """
        
        if fname=='':
            fname = 'wavelength_micron.inp'

        print 'Writing '+fname
        wfile = open(fname, 'w')
        wfile.write('%d\n'%self.nwav)
        for ilam in range(self.nwav):
            wfile.write('%.9e\n'%self.wav[ilam])
        wfile.close()
       
# --------------------------------------------------------------------------------------------------
    def makeSpatialGrid(self,crd_sys=None,xbound=None,ybound=None,zbound=None,nxi=None,nyi=None,nzi=None,ppar=None):
        """
        Function to create the spatial grid

        INPUT:
        ------
            crd_sys     - 'car'/'sph'  Coordinate system of the spatial grid
            xbound      - List (with at least two elements) of boundaries for the grid along the first dimension
            ybound      - List (with at least two elements) of boundaries for the grid along the second dimension
            zbound      - List (with at least two elements) of boundaries for the grid along the third dimension
            nxi         - Number of grid points along the first dimension. List with len(xbound)-1 elements with 
                            nxi[i] being the number of grid points between xbound[i] and xbound[i+1]
            nyi         - Same as nxi but for the second dimension
            nzi         - Same as nxi but for the third dimension
        
        OPTIONS:
        --------
            ppar        - Dictionary containing all input parameters of the model (from the problem_params.inp file)
                          if ppar is set all keyword arguments that are not set will be taken from this dictionary
        """

        self.act_dim = [1,1,1]
        if ppar:
            if not crd_sys : crd_sys = ppar['crd_sys']
            self.crd_sys =  crd_sys
           
            if not xbound : 
                if ppar.has_key('xbound'):
                    xbound = ppar['xbound']
                else:
                    print ' No boundary for the first dimension is given, first dimension is deactivated.'
                    self.act_dim[0] = 0
            if not nxi:
                if ppar.has_key('nx'):
                    if (type(ppar['nx']).__name__!='list'): 
                        ppar['nx'] = [ppar['nx']]
                    nxi = [i+1 for i in ppar['nx']] #nxi = ppar['nx']+1
                    if ppar['nx'][0]==0:
                        self.act_dim[0] = 0
                else:
                    self.act_dim[0] = 0


            if not ybound : 
                if ppar.has_key('ybound'):
                    ybound = ppar['ybound']
                else:
                    print ' No boundary for the second dimension is given, second dimension is deactivated.'
                    self.act_dim[1] = 0
            if not nyi:
                if ppar.has_key('ny'):
                    if (type(ppar['ny']).__name__!='list'): 
                        nyi = [ppar['ny']+1]
                        ppar['ny'] = [ppar['ny']]

                    else:
                        ppar['ny'] = ppar['ny']
                        nyi = [i for i in ppar['ny']] #ppar['ny']+1
                    
                    if ppar['ny'][0]==0:
                        self.act_dim[1] = 0
                else:
                    self.act_dim[1] = 0

            if not zbound : 
                if ppar.has_key('zbound'):
                    zbound = ppar['zbound']
                else:
                    print ' No boundary for the third dimension is given, third dimension is deactivated.'
                    self.act_dim[2] = 0
            if not nzi:
                if (ppar.has_key('nz'))&(ppar['nz']>0.):
                    if (type(ppar['nz']).__name__!='list'): 
                        ppar['nz'] = [ppar['nz']]
                    nzi = [i+1 for i in ppar['nz']] #nzi = ppar['nz']+1
                    if ppar['nz'][0]==0:
                        self.act_dim[2] = 0
                else:
                    self.act_dim[2] = 0
                    nzi = [0]

        if (crd_sys=='car'):
#
# First check whether the grid boundaries are specified
#
            if (xbound==None): 
                print 'ERROR'
                print 'Boundaries on the cartesian x-axis is not specified'
                print 'Without the boundaries no grid can be created'
                return
            
            if (ybound==None): 
                print 'ERROR'
                print 'Boundaries on the cartesian y-axis is not specified'
                print 'Without the boundaries no grid can be created'
                return
            if (zbound==None): 
                print 'ERROR'
                print 'Boundaries on the cartesian z-axis is not specified'
                print 'Without the boundaries no grid can be created'
                return
            
            if ((nxi==None)|(nyi==None)|(nzi==None)):
                print 'ERROR'
                print 'Number of grid points is not specified'
                return

#
# Type checking 
#

            if (type(nxi).__name__=='int'):  nxi = [nxi]
            if (type(nyi).__name__=='int'):  nyi = [nyi]
            if (type(nzi).__name__=='int'):  nzi = [nzi]

#
# Create the x-axis
#
            if (len(nxi)>1): 
                self.nxi = sum(nxi)
                self.nx  = self.nxi-1
                self.xi  = xbound[0] + (xbound[1] - xbound[0])*(arange(nxi[0], dtype=float64)/float(nxi[0]))
                for ipart in range(1,len(nxi)-1):
                    dum = xbound[ipart] + (xbound[ipart+1] - xbound[ipart])*(arange(nxi[ipart], dtype=float64)/float(nxi[ipart]))
                    self.xi = append(self.xi, dum)

                ipart = len(nxi)-1 
                dum = xbound[ipart] + (xbound[ipart+1] - xbound[ipart])*(arange(nxi[ipart], dtype=float64)/float(nxi[ipart]-1))
                self.xi = append(self.xi, dum)
                self.x  = 0.5*(self.xi[0:self.nx] + self.xi[1:self.nx+1])
            else:
                if self.act_dim[0]==1:
                    self.nxi = nxi[0]
                    self.xi = xbound[0] + (xbound[1] - xbound[0])*(arange(self.nxi, dtype=float64)/float(self.nxi-1.))
                    self.nx = self.nxi-1
                    self.x  = 0.5*(self.xi[0:self.nx] + self.xi[1:self.nx+1])
                else:
                    self.x = [0.]
                    self.xi = [0., 0.,]
                    self.nx = 1
                    self.nxi = 2

#
# Create the y-ayis
#
            if (len(nyi)>1): 
                self.nyi = sum(nyi)
                self.ny  = self.nyi-1
                self.yi  = ybound[0] + (ybound[1] - ybound[0])*(arange(nyi[0], dtype=float64)/float(nyi[0]))
                for ipart in range(1,len(nyi)-1):
                    dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*(arange(nyi[ipart], dtype=float64)/float(nyi[ipart]))
                    self.yi = append(self.yi, dum)

                ipart = len(nyi)-1 
                dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*(arange(nyi[ipart], dtype=float64)/float(nyi[ipart]-1))
                self.yi = append(self.yi, dum)
                self.y  = 0.5*(self.yi[0:self.ny] + self.yi[1:self.ny+1])
            else:
                if self.act_dim[0]==1:
                    self.nyi = nyi[0]
                    self.yi = ybound[0] + (ybound[1] - ybound[0])*(arange(self.nyi, dtype=float64)/float(self.nyi-1.))
                    self.ny = self.nyi-1
                    self.y  = 0.5*(self.yi[0:self.ny] + self.yi[1:self.ny+1])
                else:
                    self.y = [0.]
                    self.yi = [0., 0.,]
                    self.ny = 1
                    self.nyi = 2


#
# Create the z-azis
#
            if (len(nzi)>1): 
                self.nzi = sum(nzi)
                self.nz  = self.nzi-1
                self.zi  = zbound[0] + (zbound[1] - zbound[0])*(arange(nzi[0], dtype=float64)/float(nzi[0]))
                for ipart in range(1,len(nzi)-1):
                    dum = zbound[ipart] + (zbound[ipart+1] - zbound[ipart])*(arange(nzi[ipart], dtype=float64)/float(nzi[ipart]))
                    self.zi = append(self.zi, dum)

                ipart = len(nzi)-1 
                dum = zbound[ipart] + (zbound[ipart+1] - zbound[ipart])*(arange(nzi[ipart], dtype=float64)/float(nzi[ipart]-1))
                self.zi = append(self.zi, dum)
                self.z  = 0.5*(self.zi[0:self.nz] + self.zi[1:self.nz+1])
            else:
                if self.act_dim[0]==1:
                    self.nzi = nzi[0]
                    self.zi = zbound[0] + (zbound[1] - zbound[0])*(arange(self.nzi, dtype=float64)/float(self.nzi-1.))
                    self.nz = self.nzi-1
                    self.z  = 0.5*(self.zi[0:self.nz] + self.zi[1:self.nz+1])
                else:
                    self.z = [0.]
                    self.zi = [0., 0.0]
                    self.nz = 1
                    self.nzi = 2



        if (crd_sys=='sph'): 
#
# r->x, theta->y, phi-z            
#
            if (xbound==None): 
                print 'ERROR'
                print 'Boundaries on the radius is not specified'
                print 'Without the boundaries no grid can be created'
                return

            if (ybound==None): ybound = [0.0, pi]
            if (zbound==None): zbound = [0.0, 2.0*pi]

            if ((nxi==None)|(nyi==None)|(nzi==None)):
                print 'ERROR'
                print 'Number of grid points is not specified'
                return

#
# Type checking (what is in the dimension numbers)
#

            if (type(nxi).__name__=='int'):  nxi = [nxi]
            if (type(nyi).__name__=='int'):  nyi = [nyi]
            if (type(nzi).__name__=='int'):  nzi = [nzi]
#
# Create the x axis
#
            if (len(nxi)>1): 
                self.nxi = sum(nxi)
                self.nx  = self.nxi-1
                self.xi  = xbound[0] * (xbound[1] / xbound[0])**(arange(nxi[0], dtype=float64)/float(nxi[0]))
                for ipart in range(1,len(nxi)-1):
                    dum = xbound[ipart] * (xbound[ipart+1] / xbound[ipart])**(arange(nxi[ipart], dtype=float64)/float(nxi[ipart]))
                    self.xi = append(self.xi, dum)

                ipart = len(nxi)-1 
                dum = xbound[ipart] * (xbound[ipart+1] / xbound[ipart])**(arange(nxi[ipart], dtype=float64)/float(nxi[ipart]-1))
                self.xi = append(self.xi, dum)
                self.x  = sqrt(self.xi[0:self.nx] * self.xi[1:self.nx+1])
            else:
                if self.act_dim[0]==1:
                    self.nxi = nxi[0]
                    self.xi = xbound[0] * (xbound[1] / xbound[0])**(arange(self.nxi, dtype=float64)/float(self.nxi-1.))
                    self.nx = self.nxi-1
                    self.x  = sqrt(self.xi[0:self.nx] * self.xi[1:self.nx+1])
                else:
                    self.x = [0.]
                    self.xi = [0., 0.,]
                    self.nx = 1
                    self.nxi = 2
                
            ## This has to be done properly
            #if ppar.has_key('xres_nlev'):
                #ri_ext = array([self.xi[0], self.xi[ppar['xres_nspan']]])
                #for i in range(ppar['xres_nlev']):
                    #dum_ri = ri_ext[0] + (ri_ext[1]-ri_ext[0]) * arange(ppar['xres_nstep']+1, dtype=float64) / float(ppar['xres_nstep'])
                    #print ri_ext[0:2]/au
                    #print dum_ri/au
                    #ri_ext_old = array(ri_ext)
                    #ri_ext = array(dum_ri)
                    #ri_ext = append(ri_ext,ri_ext_old[2:])
                    #print ri_ext/au
                    #print '----------'
                    
                #r_ext = (ri_ext[1:] + ri_ext[:-1]) * 0.5

                #self.xi = append(ri_ext, self.xi[ppar['xres_nspan']+1:])
                #self.x = append(r_ext, self.x[ppar['xres_nspan']:])
                #self.nx = self.x.shape[0]
                #self.nxi = self.xi.shape[0]

            # Refinement of the inner edge of the grid
            # This has to be done properly
            if ppar.has_key('xres_nlev'):
                ri_ext = array([self.xi[0], self.xi[ppar['xres_nspan']]])
                for i in range(ppar['xres_nlev']):
                    dum_ri = ri_ext[0] + (ri_ext[1]-ri_ext[0]) * arange(ppar['xres_nstep']+1, dtype=float64) / float(ppar['xres_nstep'])
                    #print ri_ext[0:2]/au
                    #print dum_ri/au
                    ri_ext_old = array(ri_ext)
                    ri_ext = array(dum_ri)
                    ri_ext = append(ri_ext,ri_ext_old[2:])
                    
                r_ext = (ri_ext[1:] + ri_ext[:-1]) * 0.5

                self.xi = append(ri_ext, self.xi[ppar['xres_nspan']+1:])
                self.x = append(r_ext, self.x[ppar['xres_nspan']:])
                self.nx = self.x.shape[0]
                self.nxi = self.xi.shape[0]
                
#
# Create the y axis
#
            if (len(nyi)>1):
                
                # Check if we go to the full [0,pi] interval or only use the upper half-plane [0, pi/2]
                
                if ybound[len(ybound)-1]!=pi/2.:
                    self.nyi = sum(nyi)+1
                    self.ny  = self.nyi-1
                    self.yi  = ybound[0] + (ybound[1] - ybound[0])*(arange(nyi[0], dtype=float64)/float(nyi[0]))                
                    for ipart in range(1,len(nyi)-1):
                        # Now make sure that pi/2 will be a cell interface
                        # 
                        # BUGFIX! 16-05-2012
                        # The grid was not symmetric to pi/2 when the grid contained multiple sections (i.e. len(nyi)>1)
                        # This is now fixed
                        if (ybound[ipart]<pi/2.):
                            dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*(arange(nyi[ipart], dtype=float64)/float(nyi[ipart]))
                        else:
                            if (ybound[ipart]==pi/2.):
                                dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((arange(nyi[ipart]+1, dtype=float64))/(float(nyi[ipart])))
                            else:
                                dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((arange(nyi[ipart], dtype=float64)+1.)/float(nyi[ipart]))

                        self.yi = append(self.yi, dum)

                    ipart   = len(nyi)-1 
                    if len(nyi)==2:
                        dum     = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((arange(nyi[ipart], dtype=float64))/(float(nyi[ipart])-1.))
                    else:
                        dum     = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((arange(nyi[ipart], dtype=float64)+1.)/float(nyi[ipart]))

                else:
                    self.nyi = sum(nyi)+1
                    self.ny  = self.nyi-1
                    self.yi  = ybound[0] + (ybound[1] - ybound[0])*(arange(nyi[0], dtype=float64)/float(nyi[0]))                
                    for ipart in range(1,len(nyi)-1):
                        # Now make sure that pi/2 will be a cell interface
                        # 
                        # BUGFIX! 16-05-2012
                        # The grid was not symmetric to pi/2 when the grid contained multiple sections (i.e. len(nyi)>1)
                        # This is now fixed
                        if (ybound[ipart]<pi/2.):
                            dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*(arange(nyi[ipart], dtype=float64)/float(nyi[ipart]))
                        else:
                            dum = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((arange(nyi[ipart]+1, dtype=float64))/(float(nyi[ipart])))
                        
                        self.yi = append(self.yi, dum)

                    ipart   = len(nyi)-1 

                    if len(nyi)==2:
                        dum     = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((arange(nyi[ipart]+1, dtype=float64))/(float(nyi[ipart])))
                    else:
                        dum     = ybound[ipart] + (ybound[ipart+1] - ybound[ipart])*((arange(nyi[ipart], dtype=float64)+1.)/float(nyi[ipart]))



            
                self.yi = append(self.yi, dum)
                self.y  = 0.5*(self.yi[0:self.ny] + self.yi[1:self.ny+1])

            else:
                if self.act_dim[1]==1:
                    self.nyi = nyi[0]
                    self.yi = ybound[0] + (ybound[1] - ybound[0])*(arange(self.nyi, dtype=float64)/float(self.nyi-1.))
                    self.ny = self.nyi-1
                    self.y  = 0.5*(self.yi[0:self.ny] + self.yi[1:self.ny+1])
                else:
                    self.y = [0.]
                    self.yi = [0., 0.,]
                    self.ny = 1
                    self.nyi = 2
#
# Create the z axis

            if (len(nzi)>1):
                self.nzi = sum(nzi)
                self.nz  = self.nzi-1

                self.zi  = zbound[0] + (zbound[1] - zbound[0])*(arange(nzi[0], dtzpe=float64)/float(nzi[0]))                
                for ipart in range(1,len(nzi)-1):
                    dum = zbound[ipart] + (zbound[ipart+1] - zbound[ipart])*(arange(nzi[ipart], dtzpe=float64)/float(nzi[ipart]))
                    self.zi = append(self.zi, dum)
                ipart   = len(nzi)-1 
                dum     = zbound[ipart] + (zbound[ipart+1] - zbound[ipart])*(arange(nzi[ipart], dtzpe=float64)/float(nzi[ipart]-1))
                self.zi = append(self.zi, dum)
                self.z  = 0.5*(self.zi[0:self.nz] + self.zi[1:self.nz+1])
            else:
                if self.act_dim[2]==1:
                    self.nzi = nzi[0]
                    self.zi = zbound[0] + (zbound[1] - zbound[0])*(arange(self.nzi, dtype=float64)/float(self.nzi-1))
                    self.nz = self.nzi-1
                    self.z  = 0.5*(self.zi[0:self.nz] + self.zi[1:self.nz+1])
                else:
                    self.z = [0.]
                    self.zi = [0., np.pi*2.]
                    self.nz = 1
                    self.nzi = 2
                

        #if (crd_sys!='sph'):
            #print 'WARNING:'
            #print 'Currently only spherical coordinate system is supported!'
            #return

# --------------------------------------------------------------------------------------------------
    def writeSpatialGrid(self, fname=''):
        """
        Function to write the wavelength grid to a file (e.g. amr_grid.inp)

        OPTIONS:
        --------
            fname - File name into which the spatial grid should be written. If omitted 'amr_grid.inp' will be used. 
        """
        


        if fname=='':
            fname = 'amr_grid.inp'

        print 'Writing '+fname
        wfile = open(fname, 'w')
        wfile.write('%d\n'%1)                    # Format number
        wfile.write('%d\n'%0)                    # AMR self.style (0=regular self. NO AMR)
        if self.crd_sys=='car':
            wfile.write('%d\n'%0)                  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
        if self.crd_sys=='sph':
            wfile.write('%d\n'%100)                  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
        if self.crd_sys=='cyl':
            wfile.write('%d\n'%200)                  # Coordinate system (0-99 cartesian, 100-199 spherical, 200-299 cylindrical)
        wfile.write('%d\n'%0)                    # Gridinfo
        
        wfile.write('%d %d %d \n'%(self.act_dim[0], self.act_dim[1], self.act_dim[2]))       # Which dimension is active
        wfile.write('%d %d %d \n'%(self.nx,self.ny,self.nz))    # Grid size (x,y,z or r,phi,theta, or r,phi,z)
        for i in range(self.nxi): wfile.write('%.9e\n'%self.xi[i])
        for i in range(self.nyi): wfile.write('%.9e\n'%self.yi[i])
        for i in range(self.nzi): wfile.write('%.9e\n'%self.zi[i])
        wfile.close()


# --------------------------------------------------------------------------------------------------
    def readGrid(self, fname=''):
        """
        Function to read the spatial (amr_grid.inp) and frequency grid (wavelength_micron.inp).
        
        OPTIONS:
        --------
            fname - File name from which the spatial grid should be read. If omitted 'amr_grid.inp' will be used. 
        """
#
# Natural constants
#

        cc = 29979245800.

        if fname=='':
            fname = 'amr_grid.inp'

# 
# Read the spatial grid 
#
        try :
            rfile = open(fname, 'r')
        except:
            print 'Error!' 
            print 'amr_grid.inp was not found!'
            return 
    
        form        = float(rfile.readline())
        grid_style  = float(rfile.readline())
        crd_system  = int(rfile.readline())
        if crd_system<100:
            self.crd_sys = 'car'
        elif ((crd_system>=100)&(crd_system<200)):
            self.crd_sys = 'sph'
        elif ((crd_system>=200)&(crd_system<300)):
            self.crd_sys = 'cyl'
        else:
            rfile.close()
            print 'ERROR'
            print ' unsupported coordinate system in the amr_grid.inp file'
            print crd_system
            return

        grid_info   = float(rfile.readline())
        dum         = rfile.readline().split()
        self.act_dim = [int(dum[i]) for i in range(len(dum))]
        dum         = rfile.readline().split()
        self.nx,self.ny,self.nz    = int(dum[0]), int(dum[1]), int(dum[2])
        self.nxi,self.nyi,self.nzi = self.nx+1, self.ny+1, self.nz+1

        self.xi           = zeros(self.nx+1, dtype=float64)
        self.yi           = zeros(self.ny+1, dtype=float64)
        self.zi           = zeros(self.nz+1, dtype=float64)
       
        for i in range(self.nxi): self.xi[i] = float(rfile.readline())
        for i in range(self.nyi): self.yi[i] = float(rfile.readline())
        for i in range(self.nzi): self.zi[i] = float(rfile.readline())

        if self.crd_sys=='car':
            self.x = (self.xi[0:self.nx] +  self.xi[1:self.nx+1]) * 0.5
            self.y = (self.yi[0:self.ny] +  self.yi[1:self.ny+1]) * 0.5
            self.z = (self.zi[0:self.nz] +  self.zi[1:self.nz+1]) * 0.5
        else: 
            self.x = sqrt(self.xi[0:self.nx] * self.xi[1:self.nx+1])
            self.y = (self.yi[0:self.ny] +  self.yi[1:self.ny+1]) * 0.5
            self.z = (self.zi[0:self.nz] +  self.zi[1:self.nz+1]) * 0.5

        rfile.close()

# 
# Read the frequency grid 
#

        try :
            rfile = open('wavelength_micron.inp', 'r')
        except:
            print 'Error!' 
            print 'wavelength_micron.inp was not found!'
            return 

        self.nwav = int(rfile.readline())
        self.nfreq = self.nwav
        self.wav  = zeros(self.nwav, dtype=float64)

        for i in range(self.nwav): self.wav[i] = float(rfile.readline())

        self.freq = cc / self.wav

        rfile.close()

# --------------------------------------------------------------------------------------------------
    def getCellVolume(self):
        """
        Function to calculate the volume of grid cells
        """

        if self.crd_sys=='sph':

            if self.act_dim[0]==0:
                print '----------------------------------------------------------'
                print 'ERROR'
                print 'The r-dimension of a spherical grid is switched off'
                print 'This model (ppdisk) is not perpared for such grid style'
                print '----------------------------------------------------------'
            elif self.act_dim[1]==0:
                print '----------------------------------------------------------'
                print 'ERROR'
                print 'The theta-dimension of a spherical grid is switched off'
                print 'This model (ppdisk) is not perpared for such grid style'
                print '----------------------------------------------------------'
            elif self.act_dim[2]==0:
                vol = zeros([self.nx, self.ny, self.nz], dtype=float64)
                diff_r3   = self.xi[1:]**3 - self.xi[:-1]**3
                diff_cost = cos(self.yi[:-1]) - cos(self.yi[1:])
                diff_phi  = 2.*pi
                for ix in range(self.nx):
                    for iy in range(self.ny):
                        vol[ix,iy,:] = 1./3. * diff_r3[ix] * diff_cost[iy] * diff_phi

            else:
                vol = zeros([self.nx, self.ny, self.nz], dtype=float64)
                diff_r3   = self.xi[1:]**3 - self.xi[:-1]**3
                diff_cost = cos(self.yi[:-1]) - cos(self.yi[1:])
                diff_phi  = self.zi[1:] - self.zi[:-1] 
                for ix in range(self.nx):
                    for iy in range(self.ny):
                        vol[ix,iy,:] = 1./3. * diff_r3[ix] * diff_cost[iy] * diff_phi
        else:
            print 'ERROR!'
            print "coordinate system '" + self.crd_sys+ "' is not yet supported"
            return 0

        return vol
# --------------------------------------------------------------------------------------------------
class radmc3dData():
    """
    RADMC3D data class
        reading and writing dust density/temperature, gas density/temperature/velocity,
        generating a legacy vtk file for visualization

    ATTRIBUTES:
    -----------
        rhodust   -  Dust density in g/cm^3 
        dusttemp  -  Dust temperature in K 
        rhogas    -  Gas density in g/cm^3
        ndens_mol -  Number density of the molecule [molecule/cm^3]
        ndens_cp  -  Number density of the collisional partner [molecule/cm^3]
        gasvel    -  Gas velocity in cm/s 
        gastemp   -  Gas temperature in K
        vturb     -  Mictroturbulence in cm/s
        taux      -  Optical depth along the x (cartesian) / r (cylindrical) / r (spherical) dimension
        tauy      -  Optical depth along the y (cartesian) / theta (cylindrical) / theta (spherical) dimension
        tauz      -  Optical depth along the z (cartesian) / z (cylindrical) / phi (spherical) dimension
        sigmadust -  Dust surface density in g/cm^2
        sigmagas  -  Gas surface density in molecule/cm^2 (or g/cm^2 depending on the dimension of rhogas)
    """
    
    def __init__(self, grid=None):

        if grid:
            self.grid = copy.deepcopy(grid)
        else:
            self.grid = radmc3dGrid()

        self.rhodust   = -1
        self.dusttemp  = -1
        self.rhogas    = -1
        self.ndens_mol = -1
        self.ndens_cp  = -1
        self.gasvel    = -1
        self.gastemp   = -1
        self.vturb     = -1
        self.taux      = -1
        self.tauy      = -1
        self.tauz      = -1
        self.sigmadust = -1
        self.sigmagas  = -1
# --------------------------------------------------------------------------------------------------
    def _scalarfieldWriter(self, data=None, fname='', binary=True):
        """
        Function to write a scalar field to a file

        INPUT:
        ------
            data   - Scalar variable to be written
            fname  - Name of the file containing a scalar variable
            binary - If True the file will be in binary format, if False the file format is formatted ASCII text
        
        """

        wfile = open(fname, 'w')
        if binary:
            if len(data.shape)==3:
                hdr = array([1, 8, self.grid.nx*self.grid.ny*self.grid.nz], dtype=int)
            elif len(data.shape)==4:
                hdr = array([1, 8, self.grid.nx*self.grid.ny*self.grid.nz,  self.rhodust.shape[3]], dtype=int)
            hdr.tofile(wfile)
            # Now we need to flatten the dust density array since the Ndarray.tofile function writes the 
            # array always in C-order while we need Fortran-order to be written
            if len(data.shape)==4:
                data = swapaxes(data,0,3)
                data = swapaxes(data,1,2)
                data.tofile(wfile)
            elif len(data.shape)==3:
                data = swapaxes(data,0,2)
                data.tofile(wfile)
            else:
                print 'ERROR'
                print 'Unknown array shape  : '
                print data.shape
                return
        else:
            
            if len(data.shape)==3:
                hdr = array([1, self.grid.nx*self.grid.ny*self.grid.nz], dtype=int)
                hdr.tofile(wfile, sep=" ", format="%d\n")
                # Now we need to flatten the dust density array since the Ndarray.tofile function writes the 
                # array always in C-order while we need Fortran-order to be written
                data = swapaxes(data,0,2)
                data.tofile(wfile, sep=" ", format="%.9e\n")


            elif len(data.shape)==4:
                hdr = array([1, self.grid.nx*self.grid.ny*self.grid.nz,  self.rhodust.shape[3]], dtype=int)
                hdr.tofile(wfile, sep=" ", format="%d\n")
                # Now we need to flatten the dust density array since the Ndarray.tofile function writes the 
                # array always in C-order while we need Fortran-order to be written
                data = swapaxes(data,0,3)
                data = swapaxes(data,1,2)
                data.tofile(wfile, sep=" ", format="%.9e\n")
            else:
                print 'ERROR'
                print 'Unknown array shape  : '
                print data.shape
                return
            
        wfile.close()

# --------------------------------------------------------------------------------------------------
    def _scalarfieldReader(self, fname='', binary=True):
        """
        Function to read a scalar field from file

        INPUT:
        ------
            fname  - Name of the file containing a scalar variable
            binary - If True the file is in binary format, if False the file format is formatted ASCII text
        
        OUTPUT:
        -------
            Returns a numpy Ndarray with the scalar field
        """

        if binary:
            # hdr[0] = format number
            # hdr[1] = data precision (4=single, 8=double)
            # hdr[2] = nr of cells
            # hdr[3] = nr of dust species
            hdr = fromfile(fname, count=4, dtype=int)
            if hdr[2]!=(self.grid.nx*self.grid.ny*self.grid.nz):
                print ' ERROR'
                print ' Number of grid points in '+fname+' is different from that in amr_grid.inp'
                print npoints
                print hdr[2]
                return

            if hdr[1]==8:
                data = fromfile(fname, count=-1, dtype=float64)
            elif hdr[1]==4:
                data = fromfile(fname, count=-1, dtype=float)
            else:
                print 'ERROR'
                print 'Unknown datatype in '+fname
                return
            

            if data.shape[0]==(hdr[2]+3):
                data = reshape(data[3:], [1, self.grid.nz,self.grid.ny,self.grid.nx])
            elif data.shape[0]==(hdr[2]*hdr[3]+4):
                data = reshape(data[4:], [hdr[3],self.grid.nz,self.grid.ny,self.grid.nx])

            
            #data = reshape(data, [hdr[3],self.grid.nz,self.grid.ny,self.grid.nx])
            # We need to change the axis orders as Numpy always writes binaries in C-order while RADMC3D
            # uses Fortran-order
            data = swapaxes(data,0,3)
            data = swapaxes(data,1,2)

        else:
            rfile = -1
            try :
                rfile = open(fname, 'r')
            except:
                print 'Error!' 
                print fname+' was not found!'
                
             
            if (rfile!=(-1)):

                hdr = fromfile(fname, count=3, sep="\n", dtype=int)
                
                if ((self.grid.nx * self.grid.ny * self.grid.nz)!=hdr[1]):
                    print 'Error!'
                    print 'Number of grid points in amr_grid.inp is not equal to that in '+fname
                else:

                    data = fromfile(fname, count=-1, sep="\n", dtype=float64)
                    if data.shape[0]==hdr[1]+2:
                        data = reshape(data[2:], [1, self.grid.nz,self.grid.ny,self.grid.nx])
                    elif data.shape[0]==hdr[1]*hdr[2]+3:
                        data = reshape(data[3:], [hdr[2],self.grid.nz,self.grid.ny,self.grid.nx])
                    # We need to change the axis orders as Numpy always reads  in C-order while RADMC3D
                    # uses Fortran-order
                    data = swapaxes(data,0,3)
                    data = swapaxes(data,1,2)
            
            else:
                data = -1

            if rfile!=(-1):
                rfile.close()
        return data

# --------------------------------------------------------------------------------------------------
    def  getTauOneDust(self, idust=0, axis='', kappa=0.):
        """
        Function to calculate the optical depth of a single dust species along any given combination of the axes 

        INPUT:
        ------
            idust - Index of the dust species whose optical depth should be calculated
            axis  - Name of the axis/axes along which the optical depth should be calculated 
                    (e.g. 'x' for the first dimension or 'xyz' for all three dimensions)
            kappa - Mass extinction coefficients of the dust species at the desired wavelength
        
        OUTPUT:
        -------
            Returns a dictionary with the following keys;
            
            taux  - optical depth along the first dimension
            tauy  - optical depth along the second dimension
            
            (tauz is not yet implemented)
        """
    
        # Check along which axis should the optical depth be calculated
        do_taux = False
        do_tauy = False

        if axis.find('x')>=0 : do_taux = True
        if axis.find('y')>=0 : do_tauy = True
       
        # Calculate the optical depth along the x-axis (r in spherical coordinates)
        if do_taux:
            taux = zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=float64)
            diff_x    = self.grid.xi[1:] - self.grid.xi[:-1]
            taux[0,:,:] = self.rhodust[0,:,:,idust] * kappa * diff_x[0] 
            for ix in range(1,self.grid.nx):
                taux[ix,:,:] = taux[ix-1,:,:] + self.rhodust[ix,:,:,idust] * kappa * diff_x[ix] 
        else:
            taux = [-1.]
        
        # Calculate the optical depth along the theta in spherical coordinates
        # Warning the formulation below is valid only in spherical coordinate sytem

        dum_x = zeros([self.grid.nx, self.grid.nz], dtype=float64)
        for iz in range(self.grid.nz):
            dum_x[:,iz] = self.grid.x

        if do_tauy:
            tauy = zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=float64)
            diff_y    = self.grid.yi[1:] - self.grid.yi[:-1]
            tauy[:,0,:] = self.rhodust[:,0,:,idust] * kappa * diff_y[0] * dum_x
            for iy in range(1,self.grid.ny):
                tauy[:,iy,:] = tauy[:,iy-1,:] + self.rhodust[:,iy,:,idust] * kappa * diff_y[iy] * dum_x
        else:
            tauy = [-1.]
        return {'taux':taux, 'tauy':tauy}

# --------------------------------------------------------------------------------------------------
    def  getTau(self, idust=[], axis='xy', wav=0., kappa=None):
        """
        Function to calculate the optical depth along any given combination of the axes 

        INPUT:
        ------
            idust : List of dust component indices whose optical depth should be calculated
                    If multiple indices are set the total optical depth is calculated summing 
                    over all dust species in idust
            axis  - Name of the axis/axes along which the optical depth should be calculated 
                    (e.g. 'x' for the first dimension or 'xyz' for all three dimensions)
            wav   : Wavelength at which the optical depth should be calculated
            kappa : If set it should be a list of mass extinction coefficients at the desired wavelength
                    The number of elements in the list should be equal to that in the idust keyword

        """
        # Check if the input idust indices can be found in rhoudust 
        if len(self.rhodust.shape)==3: 
            ndust = 1
        else:
            ndust = self.rhodust.shape[3]

        if len(idust)==0:
            idust = arange(ndust)

        if max(idust)>ndust:
            print 'ERROR'
            print ' There are less number of dust species than some of the indices in idust'
            return -1


        # If the kappa keyword is set it should be used during the optical depth calculation
        if kappa:
            # Safety check
            if len(kappa)!=len(idust):
                print 'ERROR'
                print ' The number of kappa values should be identical to the number of specified dust species '
                return -1
        else:  
            # Read the master opacity file to get the dustkappa file name extensions
            dum = radmc3dDustOpac()
            dummy_ext = dum.readMasterOpac()['ext']
            if len(dummy_ext)<=max(idust):
                print 'ERROR'
                print 'There are less dust species specified in dustopac.inp than some of the specified idust indices'
                return -1
            else:
                ext = [dummy_ext[i] for i in idust]

        if axis.find('x')>=0:
            self.taux = zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=float64) 
        if axis.find('y')>=0:
            self.tauy = zeros([self.grid.nx, self.grid.ny, self.grid.nz], dtype=float64) 

        for i in idust:

            if kappa==None:
                opac = readOpac(ext=ext[i])
                if opac.ext==[]:
                    return -1
                else:
                    kabs = 10.**interp(log10(array(wav)), log10(opac.wav[0]), log10(opac.kabs[0]))
                if opac.ksca[0][0]>0:
                    ksca = 10.**interp(log10(array(wav)), log10(opac.wav[0]), log10(opac.ksca[0]))
                else:
                    ksca = array(kabs)*0.
            
                print ' Opacity at '+("%.2f"%wav)+'um : ', kabs+ksca
                dum  = self.getTauOneDust(i, axis=axis, kappa=kabs + ksca)
            else:
                dum  = self.getTauOneDust(i, axis=axis, kappa=kappa[i])

            if axis.find('x')>=0:
                self.taux = self.taux + dum['taux']
            if axis.find('y')>=0:
                self.tauy = self.tauy + dum['tauy']

# --------------------------------------------------------------------------------------------------
    def readDustDens(self, fname='', binary=True):
        """
        Function to read the dust density

        OPTIONS:
        --------
            fname - Name of the file that contains the dust density. If omitted 'dust_density.inp' is used
                    (or if binary=True the 'dust_density.binp' is used).
            binary - If true the data will be read in binary format, otherwise the file format is ascii
        """
    
        if (self.grid.nx==-1):
            self.grid.readGrid()
            
        print 'Reading dust density'

        if binary:
            if fname=='':
                fname = 'dust_density.binp'
        else:
            if fname=='':
                fname = 'dust_density.inp'
            
        self.rhodust = self._scalarfieldReader(fname=fname, binary=binary)
# --------------------------------------------------------------------------------------------------
    def readDustTemp(self, fname='', binary=True):
        """
        Function to read the dust temperature

        OPTIONS:
        --------
            fname - Name of the file that contains the dust temperature. 
                    If omitted 'dust_temperature.dat' (if binary=True 'dust_temperature.bdat')is used.
            binary - If true the data will be read in binary format, otherwise the file format is ascii
        """
       

        if (self.grid.nx==-1):
            self.grid.readGrid()
            
        print 'Reading dust temperature'

        if binary:
            if fname=='':
                fname = 'dust_temperature.bdat'
        else:
            if fname=='':
                fname = 'dust_temperature.dat'
            
        self.dusttemp = self._scalarfieldReader(fname=fname, binary=binary)

# --------------------------------------------------------------------------------------------------
    def readGasVel(self, fname='', binary=True):
        """
        Function to read the gas velocity.  
        
        OPTIONS:
        --------
            fname - Name of the file that contains the gas velocity
                    If omitted 'gas_velocity.inp' (if binary=True 'gas_velocity.binp')is used.
            binary - If true the data will be read in binary format, otherwise the file format is ascii

        """

        if binary:
            if fname=='':
                fname = 'gas_velocity.binp'
            if (self.grid.nx==-1):
                self.grid.readGrid()

            print 'Reading gas velocity'
            
            try :
                rfile = open(fname, 'r')
            except:
                print 'Error!' 
                print fname+' was not found!'

            if (rfile!=(-1)):            
                hdr = fromfile(fname, count=3, dtype=int)
                if (hdr[2]!=self.grid.nx*self.grid.ny*self.grid.nz):
                    print 'ERROR'
                    print 'Number of grid points in '+fname+' is different from that in amr_grid.inp'
                    print self.grid.nx, self.grid.ny, self.grid.nz
                    print hdr[1]
                    return

                if hdr[1]==8:
                    self.gasvel = fromfile(fname, count=-1, dtype=float64)
                elif hdr[1]==4:
                    self.gasvel = fromfile(fname, count=-1, dtype=float)
                else:
                    print 'ERROR'
                    print 'Unknown datatype in '+fname
                    return
                self.gasvel = reshape(self.gasvel[3:], [self.grid.nz,self.grid.ny,self.grid.nx,3])
                self.gasvel = swapaxes(self.gasvel, 0, 2)

            else:
                self.gasvel=-1
            

        else:
            if fname=='':
                fname = 'gas_velocity.inp'

            if (self.grid.nx==-1):
                self.grid.readGrid()

            print 'Reading gas velocity'

            rfile = -1

            try :
                rfile = open(fname, 'r')
            except:
                print 'Error!' 
                print fname+' was not found!'
                
            if (rfile!=(-1)):            
                dum = rfile.readline()
                dum = int(rfile.readline())
                
                if ((self.grid.nx * self.grid.ny * self.grid.nz)!=dum):
                    print 'Error!'
                    print 'Number of self.grid.points in amr_grid.inp is not equal to that in gas_velocity.inp'
                else:
                    
                    self.gasvel = zeros([self.grid.nx, self.grid.ny, self.grid.nz, 3], dtype=float64)
                    
                    for k in range(self.grid.nz):
                        for j in range(self.grid.ny):
                            for i in range(self.grid.nx):
                                dum = rfile.readline().split()
                                self.gasvel[i,j,k,0] = float(dum[0])
                                self.gasvel[i,j,k,1] = float(dum[1])
                                self.gasvel[i,j,k,2] = float(dum[2])
    #                            self.gasvel[i,j,k,:] = [float(dum[i]) for i in range(3)]

            else:
                self.gasvel = -1                            

            rfile.close()
# --------------------------------------------------------------------------------------------------
    def readVTurb(self, fname='', binary=True):
        """
        Function to read the turbulent velocity field. 
        
        OPTIONS:
        --------
            fname - Name of the file that contains the turbulent velocity field
                    If omitted 'microturbulence.inp' (if binary=True 'microturbulence.binp') is used.
            binary - If true the data will be read in binary format, otherwise the file format is ascii
        """
        
        if (self.grid.nx==-1):
            self.grid.readGrid()
            
        print 'Reading microturbulence'

        if binary:
            if fname=='':
                fname = 'microturbulence.binp'
        else:
            if fname=='':
                fname = 'microturbulence.inp'
            
        self.vturb = self._scalarfieldReader(fname=fname, binary=binary)
       
# --------------------------------------------------------------------------------------------------
    def readGasDens(self,ispec='',binary=True):
        """
        Function to read the gas density

        INPUT:
        ------
            ispec - File name extension of the 'numberdens_ispec.inp' (or if binary=True 'numberdens_ispec.binp') file.


        OPTIONS:
        --------
            binary - If true the data will be read in binary format, otherwise the file format is ascii

        """
        
        if (self.grid.nx==-1):
            self.grid.readGrid()
            
        print 'Reading gas density (numberdens_'+ispec+'.inp)'

        if binary:
            fname = 'numberdens_'+ispec+'.binp'
        else:
            fname = 'numberdens_'+ispec+'.inp'
            
        self.ndens_mol = self._scalarfieldReader(fname=fname, binary=binary)
       
# --------------------------------------------------------------------------------------------------

    def readGasTemp(self, fname='', binary=True):
        """
        Function to read the gas temperature

        OPTIONS:
        --------
            fname - Name of the file that contains the gas temperature. If omitted 'gas_temperature.inp' 
                    (or if binary=True 'gas_tempearture.binp') is used.
            binary - If true the data will be read in binary format, otherwise the file format is ascii
        """
      
        if (self.grid.nx==-1):
            self.grid.readGrid()
            
        print 'Reading gas temperature'

        if binary:
            if fname=='':
                fname = 'gas_temperature.binp'
        else:
            if fname=='':
                fname = 'gas_temperature.inp'
            
        self.gastemp = self._scalarfieldReader(fname=fname, binary=binary)

# --------------------------------------------------------------------------------------------------
    def writeDustDens(self, fname='', binary=True):
        """
        Function to write the dust density

        OPTIONS:
        --------
            fname - Name of the file into which the dust density should be written. If omitted 'dust_density.inp' is used.
            binary - If true the data will be written in binary format, otherwise the file format is ascii
        """
        
        if fname=='':
            if binary:
                fname = 'dust_density.binp'
            else:
                fname = 'dust_density.inp'

        print 'Writing '+fname

        self._scalarfieldWriter(data=self.rhodust, fname=fname, binary=binary)


# --------------------------------------------------------------------------------------------------
    def writeDustTemp(self, fname='', binary=True):
        """
        Function to write the dust density

        OPTIONS:
        --------
            fname - Name of the file into which the dust density should be written. If omitted 'dust_density.inp' is used.
            binary - If true the data will be written in binary format, otherwise the file format is ascii
        """
        if fname=='':
            if binary:
                fname = 'dust_temperature.bdat'
            else:
                fname = 'dust_temperature.dat'

        print 'Writing '+fname
        self._scalarfieldWriter(data=self.dusttemp, fname=fname, binary=binary)
    
# --------------------------------------------------------------------------------------------------
    def writeGasDens(self, fname='', ispec='',binary=True):
        """
        Function to write the gas density

        INPUT:
        ------
        fname  - Name of the file into which the data will be written. If omitted "numberdens_xxx.inp" and
                 "numberdens_xxx.binp" will be used for ascii and binary format, respectively (xxx is the name of the molecule).
        ispec  - File name extension of the 'numberdens_ispec.inp' (if binary=True 'numberdens_ispec.binp') 
                 file into which the gas density should be written
        binary - If true the data will be written in binary format, otherwise the file format is ascii
        """
        if ispec=='':
            print 'ERROR'
            print 'ispec keyword was not specified. This keyword is required to generate the '
            print "output file name 'numberdens_ispec.dat'" 
            return -1
        else:
            if fname=='':
                if binary:
                    fname = 'numberdens_'+ispec+'.binp'
                else:
                    fname = 'numberdens_'+ispec+'.inp'

            print 'Writing '+fname
            self._scalarfieldWriter(data=self.ndens_mol, fname=fname, binary=binary)
        
       
# --------------------------------------------------------------------------------------------------
    def writeGasTemp(self, fname='', binary=True):
        """
        Function to write the gas temperature

        OPTIONS:
        --------
            fname - Name of the file into which the gas temperature should be written. If omitted 
                    'gas_temperature.inp' (if binary=True 'gas_tempearture.binp') is used.
            binary - If true the data will be written in binary format, otherwise the file format is ascii
        """
        if fname=='':
            if binary:
                fname = 'gas_temperature.binp'
            else:
                fname = 'gas_temperature.inp'

        print 'Writing '+fname
        self._scalarfieldWriter(data=self.gastemp, fname=fname, binary=binary)
   
# --------------------------------------------------------------------------------------------------
    def writeGasVel(self, fname='', binary=True):
        """
        Function to write the gas velocity

        OPTIONS:
        --------
            fname  - Name of the file into which the gas temperature should be written. 
                    If omitted 'gas_velocity.inp' (if binary=True 'gas_velocity.binp') is used.
            binary - If true the data will be written in binary format, otherwise the file format is ascii
        """
   
        if binary:
            if fname=='':
                fname = 'gas_velocity.binp'

            wfile = open(fname, 'w')
            hdr = array([1, 8, self.grid.nx*self.grid.ny*self.grid.nz], dtype=int)
            hdr.tofile(wfile)
            # Now we need to change the axis orders since the Ndarray.tofile function writes the 
            # array always in C-order while we need Fortran-order to be written
            self.gasvel = swapaxes(self.gasvel,0,2)
            self.gasvel.tofile(wfile)

            # Switch back to the original axis order
            self.gasvel = swapaxes(self.gasvel,0,2)
            wfile.close()
        else:
            if fname=='':
                fname = 'gas_velocity.inp'

            wfile = open(fname, 'w')
            
            wfile.write('%d\n'%1)
            wfile.write('%d\n'%(self.grid.nx*self.grid.ny*self.grid.nz))

            for iz in range(self.grid.nz):
                for iy in range(self.grid.ny):
                    for ix in range(self.grid.nx):
                        wfile.write("%9e %9e %9e\n"%(self.gasvel[ix,iy,iz,0], self.gasvel[ix,iy,iz,1], self.gasvel[ix,iy,iz,2]))
                    
            wfile.close()
        print 'Writing '+fname
# --------------------------------------------------------------------------------------------------
    def writeVTurb(self, fname='', binary=True):
        """
        Function to write the microturbulence file

        OPTIONS:
        --------
            fname - Name of the file into which the turubulent velocity field should be written. 
                    If omitted 'microturbulence.inp' (if binary=True 'microturbuulence.binp') is used.
            binary - If true the data will be written in binary format, otherwise the file format is ascii
        """
   
        if fname=='':
            if binary:
                fname = 'microturbulence.binp'
            else:
                fname = 'microturbulence.inp'

        print 'Writing '+fname
        self._scalarfieldWriter(data=self.vturb, fname=fname, binary=binary)


# --------------------------------------------------------------------------------------------------
    def writeVTK(self, vtk_fname='', ddens=False, dtemp=False, idust=[0], \
                          gdens=False, gvel=False, gtemp=False):
        """
        Function to dump all physical variables to a legacy vtk file 

        INPUT:
        ------
            vtk_fname : name of the file to be written, if not specified 'radmc3d_data.vtk' will be used
            ddens     : if set to True the dust density will be written to the vtk file
            dtemp     : if set to True the dust temperature will be written to the vtk file
            idust     : a list of indices that specifies which dust component should be written 
                        if not set then the first dust species (zero index) will be used
            gdens     : if set to True the gas density will be written to the vtk file
            gtemp     : if set to True the gas temperature will be written to the vtk file
            gvel      : if set to True the gas velocity will be written to the vtk file
        """

        if (vtk_fname==''):
            vtk_fname = 'radmc3d_data.vtk'
        else:
            vtk_fname = str(vtk_fname)

#
# Get the grid 
#
        
        x  = self.grid.xi
        # For the theta axis I leave out the poles
        #  The current cell type is hexahedron and the cells near the pole are
        #    rather tetrahedra than hexahedra and this is not yet implemented
        y  = array(self.grid.yi[1:self.grid.nyi-1])
        z  = self.grid.zi
        nxi = x.shape[0]
        nyi = y.shape[0]
        nzi = z.shape[0]


#
# Gas velocity field (Should be corner centered)
# TODO
#  The lines below should be double checked and re-implemented as
#  the re-mapping of the cell centered velocity field to the cell corners
#  is physically not correct and very messy... 
#
        if gvel:
            vgas = zeros([nxi,nyi,nzi,3], dtype=float64)
            vgas[0:nxi-1,0:nyi,0:nzi-1,:]  = self.gasvel[:,1:nyi+1,:,:]
            vgas[nxi-1,:,:,:] = vgas[nxi-2,:,:,:]
            vgas[:,nyi-1,:,:] = vgas[:,nyi-2,:,:]
            vgas[:,:,nzi-1,:] = vgas[:,:,nzi-2,:]

# 
# Header 
# 

        
        wfile = open(vtk_fname, 'w')
        wfile.write('%s\n'%'# vtk DataFile Version 3.0')
        wfile.write('%s\n'%'RADMC3D Data')
        wfile.write('%s\n'%'ASCII')
        wfile.write('%s\n'%'DATASET UNSTRUCTURED_GRID')

#
# Write out the coordinates of the cell corners 
# 
        wfile.write('%s\n'%('POINTS '+str(nxi*nyi*nzi).strip()+' double'))
        print 'Writing POINTS: '
        for ix in range(nxi):
            print ix, nxi
            for iy in range(nyi):
                for iz in range(nzi):
                    crd = ctrans_sph2cart([x[ix],z[iz],y[iy]])
                    wfile.write('%.9e %9e %9e\n'%(crd[0], crd[1], crd[2]))
            
# ---------------------------------------------------------------------------------------------
# Write out the indices of the cell interface mesh that define a
# hexahedron (VTK cell type #12)
# 
# The indexing of a hexahedron is as follows
#
#                  7________6
#                 /|      / |               
#                / |     /  |               
#               4_------5   |             z ^   ^ y
#               |  3____|___2               |  /
#               | /     |  /                | /
#               |/      | /                 |/
#               0-------1                   0-----> x
#
# ---------------------------------------------------------------------------------------------
 
        wfile.write('%s %d %d\n'%('CELLS ', ((nxi-1)*(nyi-1)*(nzi-1)), ((nxi-1)*(nyi-1)*(nzi-1))*9))


        for ix in range(nxi-1):
            print 'Writing CELL COORDINATES: ', ix, self.grid.nxi-2
            for iy in range(nyi-1):
                for iz in range(nzi-1):                
                
                    id1 = nzi*nyi*ix     + nzi*iy     + iz
                    id2 = nzi*nyi*ix     + nzi*(iy+1) + iz
                    id4 = nzi*nyi*ix     + nzi*iy     + ((iz+1) % (nzi-1))
                    id3 = nzi*nyi*ix     + nzi*(iy+1) + ((iz+1) % (nzi-1))
                    id5 = nzi*nyi*(ix+1) + nzi*iy     + iz
                    id6 = nzi*nyi*(ix+1) + nzi*(iy+1) + iz
                    id7 = nzi*nyi*(ix+1) + nzi*(iy+1) + ((iz+1) % (nzi-1))
                    id8 = nzi*nyi*(ix+1) + nzi*iy     + ((iz+1) % (nzi-1))
                
                
                    line = array([8,id1,id2,id3,id4,id5,id6,id7,id8])
                    line.tofile(wfile, sep=' ', format='%d')
                    wfile.write('\n')
#
# Now write out the type of each cell (#12)
#
        wfile.write('%s %d\n'%('CELL_TYPES', ((nxi-1)*(nyi-1)*(nzi-1))))

        for ix in range(nxi-1):
            for iy in range(nyi-1):
                for iz in range(nzi-1):      
                    wfile.write('%d\n'%12)
# 
# Now write out the corner centered velocities
#
                
        if gvel:
            wfile.write('%s %d\n'%('POINT_DATA', (nxi*nyi*nzi)))
            wfile.write('%s\n'%'VECTORS gas_velocity double')
            for ix in range(nxi):
                print 'Writing velocity : ', ix, nxi-1
                for iy in range(nyi):
                    for iz in range(nzi):      
                        vsph = array([vgas[ix,iy,iz,0],vgas[ix,iy,iz,2],vgas[ix,iy,iz,1]])
                        vxyz = vtrans_sph2cart([x[ix],z[iz],y[iy]], vsph)
                
                        wfile.write('%.9e %.9e %.9e\n'%(vxyz[0], vxyz[1], vxyz[2]))
                      
# 
# Write out the cell centered scalars
# 
        wfile.write('%s %d\n'%('CELL_DATA', ((nxi-1)*(nyi-1)*(nzi-1))))

    
        if ddens:
            for ids in idust:
                wfile.write('%s\n'%('SCALARS dust_density_'+str(int(ids))+' double'))
                wfile.write('%s\n'%'LOOKUP_TABLE default')

                for ix in range(nxi-1):
                    print 'Writing dust density : ', ix, nxi-2
                    for iy in range(nyi-1):
                        for iz in range(nzi-1):
                            wfile.write('%.9e\n'%self.rhodust[ix,iy,iz,ids])

        if dtemp:
            for ids in idust:
                wfile.write('%s\n'%('SCALARS dust_temperature_'+str(int(ids))+' double'))
                wfile.write('%s\n'%'LOOKUP_TABLE default')

                for ix in range(nxi-1):
                    print 'writing dust temperature : ', ix, nxi-2
                    for iy in range(nyi-1):
                        for iz in range(nzi-1):
                            wfile.write('%.9e\n'%self.dusttemp[ix,iy,iz,ids])


        if gdens:
            wfile.write('%s\n'%'SCALARS gas_numberdensity double')
            wfile.write('%s\n'%'LOOKUP_TABLE default')

            for ix in range(nxi-1):
                print 'writing gas density : ', ix, nxi-2
                for iy in range(nyi-1):
                    for iz in range(nzi-1):
                        wfile.write('%.9e\n'%self.rhogas[ix,iy,iz])

        if gtemp:
            for ids in idust:
                wfile.write('%s\n'%('SCALARS gas_temperature double'))
                wfile.write('%s\n'%'LOOKUP_TABLE default')

                for ix in range(nxi-1):
                    print 'writing dust temperature : ', ix, nxi-2
                    for iy in range(nyi-1):
                        for iz in range(nzi-1):
                            wfile.write('%.9e\n'%self.gastemp[ix,iy,iz])
                            
# --------------------------------------------------------------------------------------------------
# Close the file
# --------------------------------------------------------------------------------------------------
        wfile.close()

# --------------------------------------------------------------------------------------------------
    def getSigmaDust(self, idust=0):
        """
        Function to calculate dust surface density 
        
        OPTIONS:
        ------
            idust - index of the dust species for which the surface density should be calculated 
                    if omitted the calculated surface density will be the sum over all dust species
        """

        # Calculate the volume of each grid cell
        vol  = self.grid.getCellVolume()
        # Dustmass in each grid cell
        if len(self.rhodust)>3:
            if idust>=0:
                mass = vol * self.rhodust[:,:,:,idust]
            else:
                mass = vol * self.rhodust.sum(3)
        else:
            mass = vol * self.rhodust

        # Calculate the surface of each grid facet in the midplane
        surf     = zeros([self.grid.nx, self.grid.nz], dtype=float64)
        diff_r2  = (self.grid.xi[1:]**2 - self.grid.xi[:-1]**2) * 0.5
        diff_phi = self.grid.zi[1:] - self.grid.zi[:-1]
        for ix in range(self.grid.nx):
            surf[ix,:] = diff_r2[ix] * diff_phi

        
        # Now get the surface density 
        dum = squeeze(mass.sum(1))
        self.sigmadust = dum / squeeze(surf)

# --------------------------------------------------------------------------------------------------
    def getSigmaGas(self):
        """
        Function to calculate gas surface density 
        This function uses radmc3dData.rhogas to calculate the surface density, thus the 
        unit of surface density depends on the unit of radmc3dData.rhogas (g/cm^2 or molecule/cm^2)
        """

        # Calculate the volume of each grid cell
        vol  = self.grid.getCellVolume()
        # Total number of molecules / gas mass in each grid cell
        mass = vol * self.rhogas
        # Calculate the surface are of each grid facet in the midplane
        surf     = zeros([self.grid.nx, self.grid.nz], dtype=float64)
        diff_r2  = (self.grid.xi[1:]**2 - self.grid.xi[:-1]**2) * 0.5
        diff_phi = self.grid.zi[1:] - self.grid.zi[:-1]  
        for ix in range(self.grid.nx):
            surf[ix,:] = diff_r2[ix] * diff_phi


        # Now get the surface density 
        dum = squeeze(mass.sum(1))
        self.sigmagas = dum / squeeze(surf)

# --------------------------------------------------------------------------------------------------
class radmc3dStars():
    """
    Class of the radiation sources (currently only stars)

    ATTRIBUTES:
    -----------
        mstar - List of stellar masses
        tstar - List of stellar effective temperatures
        rstar - List of stellar radii
        lstar - List of stellar luminosities 
        nstar - Number of stars
        pstar - Locations (coordinates) of the stars
        wav   - Wavelength for the stellar spectrum
        freq  - Frequency for the stellar spectrum
        fnu   - Stellar spectrum (flux@1pc)
        nwav  - Number of wavelenghts in the stellar spectrum
        nfreq - Number of frequencies in the stellar spectrum
    """
    def __init__(self, ppar=None):

        self.mstar    = []  
        self.tstar    = []  
        self.rstar    = []  
        self.lstar    = []  
        self.nstar    = 0
        self.pstar    = []
        self.wav      = []
        self.freq     = []
        self.fnu      = []
        self.nwav     = 0
        self.nfreq    = 0

        if ppar:
            if type(ppar['mstar']).__name__=='list':
                self.mstar = ppar['mstar']
            else:
                self.mstar = [ppar['mstar']]
                
            if type(ppar['tstar']).__name__=='list':
                self.tstar = ppar['tstar']
            else:
                self.tstar = [ppar['tstar']]

            if type(ppar['rstar']).__name__=='list':
                self.rstar = ppar['rstar']
            else:
                self.rstar = [ppar['rstar']]

            self.nstar = len(self.rstar)
            for istar in range(self.nstar):
                self.lstar.append(4.*pi*self.rstar[istar]**2. * ss* self.tstar[istar]**4.)
            self.pstar = ppar['pstar']

# --------------------------------------------------------------------------------------------------

    def findPeakStarspec(self):

        """
        Function to  calculate the peak wavelength of the stellar spectrum
       
        OUTPUT:
        -------
            Returns the peak wavelength of the stellar spectrum in nu*Fnu for all 
                stars as a list
        """
   
        pwav = []
   
        for istar in range(self.nstar):
            ii = (self.fnu[:,istar]*self.freq).argmax()
            pwav.append(self.wav[ii])

            #nufnu = self.fnu[:,istar] * self.freq
            #dpwav  = 0.0
            #dflux = 0.0

            #for iw in range(self.nwav):
                #if nufnu[iw]>dflux:
                    #dflux = nufnu[iw]
                    #dpwav  = self.wav[iw]


            #pwav.append(dpwav)

        return pwav

# --------------------------------------------------------------------------------------------------
    def readStarsinp(self, fname=''):
        """
        Function to read the stellar data from the stars.inp file

        OPTIONS:
        --------
            fname - File name of the file that should be read (if omitted stars.inp will be used)
        """
        
        if fname=='':
            fname = 'stars.inp'

        try:
            rfile = open(fname, 'r')
        except:
            print ' ERROR '
            print fname+' cannot be opened '
            return

        dum = rfile.readline()
        iformat = int(dum)
        if iformat!=2:
            print ' ERROR '
            print ' Unknown file format '
            print ' Format number : ', iformat
            rfile.close()
            return

        dum = rfile.readline().split()
        self.nstar = int(dum[0])
        self.nwav  = int(dum[1])
        self.nfreq = self.nwav
        for istar in range(self.nstar):
            dum = rfile.readline().split()
            self.rstar.append(float(dum[0]))
            self.mstar.append(float(dum[1]))
            self.pstar.append([float(dum[2]), float(dum[3]), float(dum[4])])
        
        dum = rfile.readline()
        wav = []
        for ilam in range(self.nwav):
            dum = rfile.readline()
            wav.append(float(dum))

        self.wav = array(wav, dtype=float)
        self.freq = cc/self.wav*1e4
        dum = rfile.readline()
        for istar in range(self.nstar):
            dum = rfile.readline()
            self.tstar.append(-float(dum))

        rfile.close()

        # 
        # Now calculates the stellar spectrum
        #

        self.getStellarSpectrum()

# --------------------------------------------------------------------------------------------------
    def writeStarsinp(self, wav=None, freq=None, pstar=None, tstar=None):
        """
        Writes the stars.inp file

        INPUT:
        ------
            wav   - Wavelength grid for the stellar spectrum
            freq  - Frequency grid for the stellar spectrum (either freq or wav should be set)
            pstar - List of the cartesian coordinates of the stars (each element of pstar should be a list of three elements
                    with the [x,y,z] coordinate of the individual stars)
            tstar - List containing the effective temperature of the stars
        """

        if freq!=None:
            self.wav  = cc/array(freq)
            self.freq = array(freq)
            self.nwav = self.wav.shape[0]
            self.nfreq = self.nwav

        if wav!=None:
            self.wav = array(wav)
            self.freq = cc/self.wav
            self.nwav = self.wav.shape[0]
            self.nfreq = self.nwav

        self.nstar = len(self.rstar)
        self.pstar = pstar

        # if we don't have
        print 'Writing stars.inp'
        wfile = open('stars.inp', 'w')
        wfile.write('%d\n'%2)
        wfile.write('%d %d\n'%(self.nstar,self.nwav))
        
        if (self.nstar>1):
            for istar in range(self.nstar):
                wfile.write('%.9e %.9e %.9e %.9e %.9e\n'%(self.rstar[istar], self.mstar[istar],
                    self.pstar[istar][0],self.pstar[istar][1],self.pstar[istar][2]))
        else:
            wfile.write('%.9e %.9e %.9e %.9e %.9e\n'%(self.rstar[0], self.mstar[0],
                self.pstar[0],self.pstar[1],self.pstar[2]))

        wfile.write('%s\n'%' ')
        for ilam in range(self.nwav):
            wfile.write('%.9e\n'%self.wav[ilam])
        wfile.write('%s\n'%' ')
        for istar in range(self.nstar):
            wfile.write('%.9e\n'%(-self.tstar[istar]))
        wfile.close()

# --------------------------------------------------------------------------------------------------
    def getStellarSpectrum(self, tstar=None, rstar=None, lstar=None, nu=None, wav=None):
        """
        Function to calculate a blackbody stellar spectrum

        INPUT:
        ------
            tstar : Effective temperature of the star in [K]
            rstar : Radius of the star in [cm]
            lstar : Bolometric luminosity of the star [erg/s] (either rstar or lstar should be given)
            nu    : frequency grid on which the spectrum should be calculated [Hz] 
            wav   : wavelength grid on which the spectrum should be calculated [micron] 
        """
#
# Check the input which parameters are set and which should be calculated
#

        if nu and wav:
            print 'ERROR'
            print ' Either nu or wav keyword should be set but not both!'
            return 0 
        elif wav!=None:
            self.wav   = array(wav)
            self.nwav  = self.wav.shape[0]
            self.freq  = cc/self.wav * 1e4
            self.nfreq = self.freq.shape[0]
        elif nu!=None:
            self.freq  = array(nu)
            self.nfreq = self.freq.shape[0]
            self.wav   = cc/self.freq*1e4
            self.nwav  = self.wav.shape[0]


        if tstar:
            if type(tstar).__name__!='list':
                tstar = [tstar]
            dum1 = len(tstar)
            if lstar and rstar: 
                print 'ERROR'
                print ' Only two of the input variables tstar, rstar, lstar should be set not all three'
                return 0
            elif lstar:
                if len(lstar)!=dum1:
                    print 'ERROR'
                    print 'lstar and tstar have different number of elements'
                    return 0
                else:
                    self.tstar = array(tstar)
                    self.lstar = array(lstar)
                    self.nstar = self.lstar.shape[0]
                    self.rstar = sqrt(self.lstar / (4.*pi*ss*self.tstar**4.))
        else:
            if lstar and rstar:
                if len(lstar)!=len(rstar):
                    print 'ERROR'
                    print 'lstar and rstar have different number of elements'
                    return 0
                else:
                    self.lstar = array(lstar)
                    self.rstar = array(rstar)
                    self.nstar = self.rstar.shape[0]
                    self.tstar = (self.lstar / (4.*pi*ss*self.rstar**2.))**0.25


        self.fnu   = zeros([self.nwav, self.nstar], dtype=float64)
        for istar in range(len(self.tstar)):
            self.fnu[:,istar]   = 2.*hh*self.freq**3./(exp(hh*self.freq/kk/self.tstar[istar])-1.0) * pi * self.rstar[istar]**2. / pc**2.

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class radmc3dDustOpac():
    """
    Dust opacity class

    ATTRIBUTES:
    -----------
        wav     - wavelength grid
        freq    - frequency grid
        nwav    - number of wavelengths
        kabs    - absorption coefficient per unit mass
        ksca    - scattering coefficient per unit mass
        phase_g - phase function
        ext     - if set it contains the file name extension of the duskappa_ext.Kappa file
        therm   - if False the dust grains are quantum-heated (default: True)
        idust   - index of the dust species in the dust density distribution array

        NOTE: Each attribute is a list with each element containing the corresponding data for
              a given dust species

    METHODS:
    --------
        readOpac()          - Read the dust opacity files
        readMasterOpac()   - Read the master opacity file
        writeMasterOpac()  - Write the master opacity file
        makeOpac()          - Calculates opacities with the Mie-code that comes with RADMC-3D (using the runMakedust() function)
        runMakedust()      - Runs the Mie-code to calculate dust opacities

    """
# --------------------------------------------------------------------------------------------------
    def __init__(self):

        self.wav     = []
        self.freq    = []
        self.nwav    = []
        self.nfreq   = []
        self.kabs    = []
        self.ksca    = []
        self.phase_g = []
        self.ext     = []
        self.idust   = []
        self.therm   = []
         
# --------------------------------------------------------------------------------------------------
    def  readOpac(self, ext=[''], idust=None):
        """
        Function to read the dust opacity files

        INPUT:
        ------
            ext  : file name extension (file names should look like 'dustkappa_ext.inp')
            idust: index of the dust species in the master opacity file (dustopac.inp') - starts at 0 
        """
        
        if (type(ext).__name__=='str'):  ext = [ext]
        if idust!=None:
            if (type(idust).__name__=='int'):  idust = [idust]

        if (len(ext)==1)&(ext[0]!=''):
            if idust!=None:
                print 'ERROR'
                print 'Either idust or ext should be specified, but not both'
                print idust
                print ext
                return [-1]
        
        # Read the master dust opacity file to get the dust indices and dustkappa file name extensions
        mopac = self.readMasterOpac()

        # Find the file name extensions in the master opacity file if idust is specified instead of ext
        if idust:
            ext = []
            for ispec in idust:
                if (ispec+1)>len(mopac['ext']):    
                    print 'ERROR'
                    print 'No dust species found at index ', ispec
                    return [-1]
                else:
                    ext.append(mopac['ext'][ispec])

        # If only the extension is specified look for the master opacity file and find the index of this dust species
        #  or set the index to -1 if no such dust species is present in the master opacity file
        else:
            idust = []
            for iext in ext:
                try:
                    dum2 = mopac['ext'].index(iext)
                except:
                    dum2 = -1
                idust.append(dum2)
        
        # Now read all dust opacities
        for i in range(len(ext)):
            try:
                rfile = open('dustkappa_'+ext[i]+'.inp', 'r')
            except:
                print 'ERROR'
                print ' No dustkappa_'+ext[i]+'.inp file was found'
                return -1

            self.ext.append(ext[i])

            # Read the file format
            iformat = int(rfile.readline())
            if (iformat<1)|(iformat>3):
                print 'ERROR'
                print 'Unknown file format in the dust opacity file'
                rfile.close()
                return -1


            # Read the number of wavelengths in the file
            dum = rfile.readline()
            self.nwav.append(int(dum))
            self.nfreq.append(int(dum))
            self.idust.append(idust[i])
            idu = len(self.nwav)-1

            # If only the absorption coefficients are specified
            if iformat==1:
                wav = zeros(self.nwav[idu], dtype=float64)
                kabs = zeros(self.nwav[idu], dtype=float64)
                for ilam in range(self.nwav[idu]):
                    dum = rfile.readline().split()
                    wav[ilam] = dum[0] 
                    kabs[ilam] = dum[1] 
                self.wav.append(wav)
                self.freq.append(cc/wav*1e4)
                self.kabs.append(kabs)
                self.ksca.append([-1])
                self.phase_g.append([-1])
            # If the absorption and scattering coefficients are specified
            elif iformat==2:
                wav = zeros(self.nwav[idu], dtype=float64)
                kabs = zeros(self.nwav[idu], dtype=float64)
                ksca = zeros(self.nwav[idu], dtype=float64)
                for ilam in range(self.nwav[idu]):
                    dum = rfile.readline().split()
                    wav[ilam] = dum[0] 
                    kabs[ilam] = dum[1] 
                    ksca[ilam] = dum[2] 
                self.wav.append(wav)
                self.freq.append(cc/wav*1e4)
                self.kabs.append(kabs)
                self.ksca.append(ksca)
                self.phase_g.append([-1])
            
            # If the absorption and scattering coefficients and also the scattering phase function are specified
            elif iformat==3:
                wav = zeros(self.nwav[idu], dtype=float64)
                kabs = zeros(self.nwav[idu], dtype=float64)
                ksca = zeros(self.nwav[idu], dtype=float64)
                phase_g = zeros(self.nwav[idu], dtype=float64)
                for ilam in range(self.nwav[idu]):
                    dum = rfile.readline().split()
                    wav[ilam] = dum[0] 
                    kabs[ilam] = dum[1] 
                    ksca[ilam] = dum[2] 
                    phase_g[ilam] = dum[3] 
                
                self.wav.append(wav)
                self.freq.append(cc/wav*1e4)
                self.kabs.append(kabs)
                self.ksca.append(ksca)
                self.phase_g.append(phase_g)
       
            rfile.close()
        return 0 
#--------------------------------------------------------------------------------------------------------------------
    def makeOpac(self, ppar=None, wav=None):
        """
        Function to create dust opacities for RADMC3D using MIE calculation 
        
        INPUT:
        ------
            ppar  - dictionary containing all parameter of the simulation
        
        OPTIONS:
        --------
            wav  - numpy.ndarray containing the wavelength grid on which the mass absorption coefficients should be calculated
        """

    #
    # Create the wavelength grid if it is not specified
    #
        if wav==None:
            grid = radmc3dGrid()
            grid.makeWavelengthGrid(ppar=ppar)
            wav = grid.wav

    #
    # Do we need to mix the opacities?
    #
        if type(ppar['lnk_fname']).__name__=='str':
            ppar['lnk_fname'] = [ppar['lnk_fname']]

        if len(ppar['lnk_fname'])>1:
            ext = []
            for idust in range(len(ppar['lnk_fname'])):
                
                # makedust needs the lnk file to be sorted in wavelength so create a dummy file 
                # which contains the sorted optical constants 
                try:
                    rfile = open(ppar['lnk_fname'][idust], 'r')
                except:
                    print 'ERROR'
                    print ppar['lnk_fname'][idust] + ' could not be opened'
                    return

                try:
                    w = []
                    n = []
                    k = []
                    dum = rfile.readline()
                    while len(dum)>0:
                        dum = dum.split()
                        w.append(dum[0])
                        n.append(dum[1])
                        k.append(dum[2])
                        dum = rfile.readline()

                    rfile.close()
                except:
                    print 'ERROR'
                    print ppar['lnk_fname'][idust] + ' could not be read'
                    return

                w = array(w, dtype=float)
                n = array(n, dtype=float)
                k = array(k, dtype=float)

                if float(w[0])>float(w[w.shape[0]-1]):
                    w = w[::-1]
                    n = n[::-1]
                    k = k[::-1]

                #Write out the dummy file containing the sorted optical constants
                wfile = open('opt_const.dat', 'w')
                for iwav in range(w.shape[0]):
                    wfile.write("%s %s %s \n"%(w[iwav], n[iwav], k[iwav]))
                wfile.close()

                # Run makedust
                self.runMakedust(freq=cc/wav*1e4, gmin=ppar['gsmin'], gmax=ppar['gsmax'], ngs=ppar['ngs'], \
                        lnk_fname='opt_const.dat', gdens=ppar['gdens'][idust])

                # Change the name of makedust's output
                for igs in range(ppar['ngs']):
                    dum = Popen('mv dustkappa_'+str(igs+1)+'.inp dustkappa_idust_'+str(idust+1)+'_igsize_'+str(igs+1)+'.inp', shell=True).wait()
                    ext.append('idust_'+str(idust+1)+'_igsize_'+str(igs+1))

                os.remove('opt_const.dat')

            # Mix the opacity of different dust species for a given grain size if mixing is requested
            if ppar.has_key('mixabun'):
                if len(ppar['mixabun'])==len(ppar['lnk_fname']):
                    ext = []
                    for igs in range(ppar['ngs']):
                        mixnames = ['dustkappa_igsize_'+str(igs+1)+'.inp']
                        mixspecs = [['dustkappa_idust_'+str(idust+1)+'_igsize_'+str(igs+1)+'.inp' for idust in range(len(ppar['lnk_fname']))]]
                        self.mixOpac(mixnames=mixnames, mixspecs=mixspecs, mixabun=[ppar['mixabun']])
                    
                        ext.append('igsize_'+str(igs+1))
                else:
                    print 'ERROR'
                    print ' mixabun and lnk_fname should have the same number of elements.'
                    print ' To disable mixing either set mixabun to an empty list ([]) or comment it out in the problem_params.inp file'
                    return
            
            therm = [True for i in range(len(ext))]
            self.writeMasterOpac(ext=ext, therm=therm, scattering_mode_max=ppar['scattering_mode_max'])

        else:
            # makedust needs the lnk file to be sorted in wavelength so create a dummy file 
            # which contains the sorted optical constants 
            try:
                rfile = open(ppar['lnk_fname'][0], 'r')
            except:
                print 'ERROR'
                print ppar['lnk_fname'][0] + ' could not be opened'
                return

            try:
                w = []
                n = []
                k = []
                dum = rfile.readline()
                while len(dum)>0:
                    dum = dum.split()
                    w.append(dum[0])
                    n.append(dum[1])
                    k.append(dum[2])
                    dum = rfile.readline()

                rfile.close()
            except:
                print 'ERROR'
                print ppar['lnk_fname'][0] + ' could not be read'
                return

            w = array(w, dtype=float)
            n = array(n, dtype=float)
            k = array(k, dtype=float)

            if float(w[0])>float(w[w.shape[0]-1]):
                w = w[::-1]
                n = n[::-1]
                k = k[::-1]
            
            # Write out the dummy file containing the sorted optical constants
            wfile = open('opt_const.dat', 'w')
            for iwav in range(w.shape[0]):
                wfile.write("%s %s %s \n"%(w[iwav], n[iwav], k[iwav]))
            wfile.close()

            # Run makedust
            self.runMakedust(freq=cc/wav*1e4, gmin=ppar['gsmin'], gmax=ppar['gsmax'], ngs=ppar['ngs'], \
                    lnk_fname='opt_const.dat', gdens=ppar['gdens'][0])

            # Change the name of makedust's output
            ext = []
            therm = []
            for igs in range(ppar['ngs']):
                dum = Popen('mv dustkappa_'+str(igs+1)+'.inp dustkappa_idust_1_igsize_'+str(igs+1)+'.inp', shell=True).wait()
                ext.append('idust_1_igsize_'+str(igs+1))
                therm.append(True)
#            # Change the name of makedust's output 
#            dum = Popen('mv dustkappa_1.inp dustkappa_idust_1_igsize_1.inp', shell=True).wait()
#            os.remove('opt_const.dat')

            self.writeMasterOpac(ext=ext, therm=therm, scattering_mode_max=ppar['scattering_mode_max'])
        
        # Clean up and remove dust.inp and frequency.inp
        os.remove('dust.inp')
        os.remove('frequency.inp')
# --------------------------------------------------------------------------------------------------
    def mixOpac(self, ppar=None, mixnames=[], mixspecs=[], mixabun=[], writefile=True):
        """
        Function to mix opacities

        INPUT:
        ------
            ppar     - A dictionary containing all parameters of the actual model setup
                        If any keyword is set besides ppar, the value of the separate keyword
                        will be taken instead of that in ppar. If mixname, mixspecs, and mixabun are all set
                        ppar is completely omitted and not necessary to set when mixOpac is called.
            mixnames  - Names of the files into which the mixed dust opacities will be written (not needed if writefile=False)
            mixspecs  - Names of the files from which the dust opacities are read (not needed if readfile=False)
            mixabun   - Abundances of different dust species
            writefile - If False the mixed opacities will not be written out to files given in mixnames.  
           
        """

        if writefile:
            if len(mixnames)==0:
                if ppar!=None:
                    mixnames = ppar['mixnames']
                else:
                    print 'ERROR'
                    print ' Neither ppar nor mixnames are set in mixOpac '
                    return

        if len(mixspecs)==0:
            if ppar!=None:
                mixspecs = ppar['mixspecs']
            else:
                print 'ERROR'
                print ' Neither ppar nor mixspecs are set in mixOpac '
                return
            
        if len(mixabun)==0:
            if ppar!=None:
                mixabun = ppar['mixabun']
            else:
                print 'ERROR'
                print ' Neither ppar nor mixabun are set in mixOpac '
                return

        mwav  = []
        mcabs = []
        mcsca = []
        for i in range(len(mixnames)):
            #
            # Read the dust opacities to be mixed for composite dust species #1
            #
            ocabs = []
            ocsca = []
            ogsym = []
            oform = 0
            for j in range(len(mixspecs[i])):
                try:
                    rfile=open(mixspecs[i][j], 'r')
                    form   = int(rfile.readline())
                    nwav   = int(rfile.readline())
                    dw     = zeros(nwav, dtype     = float)
                    dcabs  = zeros(nwav, dtype     = float)
                    dcsca  = zeros(nwav, dtype     = float)
                    gsym   = zeros(nwav, dtype     = float)
                    if form==1:
                        if ((oform==0)|(oform==1)):
                            oform = 1
                        else:
                            print ' '
                            print 'WARNING'
                            print ' You are trying to mix opacity tables with different formats'
                            print ' Some of the tables contain scattering coefficients while (format>=2) while other do not (format=1)'
                            print ' If you wish to continue mixing will only be done for the absorption and the output opacity table'
                            print ' will have a format number of 1.'
                            dum = raw_input('Do you wish to continue (1-yes, 0-no) ?')
                            if dum.strip()!='1':
                                return

                        for iwav in range(nwav):
                            dum = rfile.readline().split()
                            dw[iwav], dcabs[iwav] = float(dum[0]), float(dum[1])
                    if form==2:
                        if ((oform==0)|(oform==2)):
                            oform=2
                        else:
                            print ' '
                            print 'WARNING'
                            print ' You are trying to mix opacity tables with different formats'
                            print ' Some of the tables contain scattering coefficients while (format>=2) while other do not (format=1)'
                            print ' If you wish to continue mixing will only be done for the absorption and the output opacity table'
                            print ' will have a format number of 1.'
                            dum = raw_input('Do you wish to continue (1-yes, 0-no) ?')
                            if dum.strip()!='1':
                                return
                        for iwav in range(nwav):
                            dum = rfile.readline().split()
                            dw[iwav], dcabs[iwav], dcsca[iwav] = float(dum[0]), float(dum[1]), float(dum[2])
                    if form==3:
                        if ((oform==0)|(oform==3)):
                            oform=3
                        else:
                            print ' '
                            print 'WARNING'
                            print ' You are trying to mix opacity tables with different formats'
                            print ' Some of the tables contain scattering coefficients while (format>=2) while other do not (format=1)'
                            print ' If you wish to continue mixing will only be done for the absorption and the output opacity table'
                            print ' will have a format number of 1.'
                            dum = raw_input('Do you wish to continue (1-yes, 0-no) ?')
                            if dum.strip()!='1':
                                return
                        for iwav in range(nwav):
                            dum = rfile.readline().split()
                            dw[iwav], dcabs[iwav], dcsca[iwav], gsym[iwav] = float(dum[0]), float(dum[1]), float(dum[2]), float(dum[3])
                    if form>3:
                        print ' '
                        print ' ERROR'
                        print ' Unsupported dust opacity table format (format number: '+form+')'
                        print ' Currently only format number 1 and 2 are supported'
                        return
                    rfile.close()

                    if dw[1]<dw[0]:
                        print ' Dust opacity table seems to be sorted in frequency instead of wavelength'
                        print ' Reversing the arrays'
                        dw = dw[::-1]
                        dcabs = dcabs[::-1]
                        dcsca = dcsca[::-1]
                except:
                    print 'ERROR'
                    print mixspecs[i][j]+ ' could not be read'
                    return

                if j==0:
                    ocabs = array(dcabs) * mixabun[i][j]
                    ocsca = array(dcsca) * mixabun[i][j]
                    ogsym = array(gsym) * mixabun[i][j]
                    nwav0 = dw.shape[0]
                    owav  = array(dw)
                else:
                    #
                    # Interpolate dust opacities to the wavelength grid of the first dust species
                    #
                    ii = ( (owav>=dw[0])&(owav<=dw[nwav-1]) )
                    il = (owav<dw[0]) 
                    ih = (owav>dw[nwav-1])
                    dum = zeros(nwav0, dtype=float)
                    dum[ii] = 10.**interp(log10(owav[ii]), log10(dw), log10(dcabs))
                   
                    # Edwtrapolate the absorption coefficients using linear fit in log-log space (i.e. fitting a polinomial) for short wavelengths
                    der = log10(dcabs[1]/dcabs[0]) / log10(dw[1]/dw[0])
                    dum[il] = 10.**(log10(dcabs[0]) + log10(dw[0]/owav[il]))
                    
                    # Edwtrapolate the absorption coefficients using linear fit in log-log space (i.e. fitting a polinomial) for long wavelengths
                    der = log10(dcabs[nwav-1]/dcabs[nwav-2]) / log10(dw[nwav-1]/dw[nwav-2])
                    dum[ih] = 10.**(log10(dcabs[nwav-1]) + log10(owav[il]/dw[nwav-1]))
                 
                    ocabs = ocabs + array(dum) * mixabun[i][j]
                    
                    if oform==2:
                        # Do the inter-/extrapolation of for the scattering coefficients
                        dum = zeros(nwav0, dtype=float)
                        dum[ii] = 10.**interp(log10(owav[ii]), log10(dw), log10(dcsca))
                       
                        der = log10(dcsca[1]/dcsca[0]) / log10(dw[1]/dw[0])
                        dum[il] = 10.**(log10(dcsca[0]) + log10(dw[0]/owav[il]))
                        
                        der = log10(dcsca[nwav-1]/dcsca[nwav-2]) / log10(dw[nwav-1]/dw[nwav-2])
                        dum[ih] = 10.**(log10(dcsca[nwav-1]) + log10(owav[il]/dw[nwav-1]))
                       
                        ocsca = ocsca + array(dum) * mixabun[i][j]

                    if oform==3:
                        # Do the inter-/extrapolation of for the scattering phase function
                        dum = zeros(nwav0, dtype=float)
                        dum[ii] = 10.**interp(log10(owav[ii]), log10(dw), log10(gsym))
                       
                        der = log10(gsym[1]/gsym[0]) / log10(dw[1]/dw[0])
                        dum[il] = 10.**(log10(gsym[0]) + log10(dw[0]/owav[il]))
                        
                        der = log10(gsym[nwav-1]/gsym[nwav-2]) / log10(dw[nwav-1]/dw[nwav-2])
                        dum[ih] = 10.**(log10(gsym[nwav-1]) + log10(owav[il]/dw[nwav-1]))
                       
                        ogsym = ogsym + array(dum) * mixabun[i][j]


       
            #
            # Write out the mixed dust opacities
            #
            wfile = open(mixnames[i], 'w')
            wfile.write("%d\n"%oform) 
            wfile.write("%d\n"%owav.shape[0])
            if oform==1:
                for iwav in range(owav.shape[0]):
                    wfile.write("%.9e %.9e\n"%(owav[iwav], ocabs[iwav]))
            if oform==2:
                for iwav in range(owav.shape[0]):
                    wfile.write("%.9e %.9e %.9e\n"%(owav[iwav], ocabs[iwav], ocsca[iwav]))
            if oform==3:
                for iwav in range(owav.shape[0]):
                    wfile.write("%.9e %.9e %.9e %.9e\n"%(owav[iwav], ocabs[iwav], ocsca[iwav], ogsym[iwav]))

        return 
# --------------------------------------------------------------------------------------------------
    def  readMasterOpac(self):
        """
        Function to read the master opacity file 'dustopac.inp' 
        it reads the dustkappa filename extensions (dustkappa_ext.inp) corresponding to dust species indices

        OUTPUT:
        -------
        Returns a dictionary with the following keys:
            'ext'   - list of dustkappa file name extensions
            'therm' - a list of integers specifying whether the dust grain is thermal or quantum heated 
                      (0 - thermal, 1 - quantum heated)
        """
        
        try: 
            rfile = open('dustopac.inp', 'r')
        except:
            print 'Error'
            print ' No dustopac.inp file was found'
            return -1

       
        # file format
        dum = rfile.readline()
        # nr of dust species
        ndust = int(rfile.readline().split()[0])
        # Comment line
        dum = rfile.readline()

        ext = []
        therm= []
        for idust in range(ndust):
            dum = rfile.readline()
            # Check if the dust grain is thermal or quantum heated
            dum = int(rfile.readline().split()[0])
            if dum==0:
                therm.append(True)
            else:
                therm.append(False)
            # Dustkappa filename extension
            dum = rfile.readline().split()[0]
            ext.append(dum)
            #Comment line
            dum = rfile.readline()
        rfile.close()

        return {'ext':ext, 'therm':therm}
# --------------------------------------------------------------------------------------------------
    def  writeMasterOpac(self, ext=None, therm=None, scattering_mode_max=1):
        """
        Function to write the master opacity file 'dustopac.inp' 

        INPUT:
        ------
            ext : list of dustkappa file name extensions
            therm : list of integers specifying whether the dust grain is thermal or quantum heated
                    (0-thermal, 1-quantum)
        """

        print 'Writing dustopac.inp'
       
        if not ext:
            print 'ERROR'
            print 'No file name extension is specified. Without it dustopac.inp cannot be written'
            return -1
        else:
            if (type(ext).__name__=='str'):  ext = [ext]

        if therm:
            if (type(therm).__name__=='int'): therm = [therm]
            if (len(ext)!=len(therm)):
                print 'ERROR'
                print ' The number of dust species in ext and in therm are different'
                return -1
        else:
            # If therm is not specified it is assumed that all grains are thermal, no quantum heating

            therm = [True for i in range(len(ext))]

        wfile = open('dustopac.inp', 'w')

        # File format
        wfile.write('%-15s %s\n'%('2', 'Format number of this file'))
        # Number of dust species
        wfile.write('%-15s %s\n'%(str(len(ext)), 'Nr of dust species'))
        # Separator
        wfile.write('%s\n'%'============================================================================')

        for idust in range(len(ext)):
            # Dust opacity will be read from a file
            if scattering_mode_max<5:
                wfile.write('%-15s %s\n'%('1', 'Way in which this dust species is read'))
            else:
                wfile.write('%-15s %s\n'%('10', 'Way in which this dust species is read'))

            # Check if the dust grain is thermal or quantum heated
            if therm:
                if therm[idust]:
                    wfile.write('%-15s %s\n'%('0', '0=Thermal grain, 1=Quantum heated'))
            else:
                wfile.write('%-15s %s\n'%('1', '0=Thermal grain, 1=Quantum heated'))

            # Dustkappa filename extension
            wfile.write('%s %s %s\n'%(ext[idust], '    ', 'Extension of name of dustkappa_***.inp file'))
            # Separator
            wfile.write('%s\n'%'----------------------------------------------------------------------------')
            
        wfile.close()
#
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# --------------------------------------------------------------------------------------------------
    def runMakedust(self, freq=None, gmin=None, gmax=None, ngs=None, lnk_fname=None, gdens=None):
        """
        Interface function to the F77 code makedust to calculate mass absorption
        coefficients from the optical constants using Mie-theory

        INPUT:
        ------
            freq       - numpy.ndarray containing the frequency grid on which the opacities should be calculated
            gmin       - minimum grain size
            gmax       - maximum grain size
            ngs        - number of grain sizes
            gdens      - density of the dust grain in g/cm^3
            lnk_faname - name of the file in which the optical constants are stored

        OUTPUT:
        -------
            result         - numpy.ndarray[nfreq,ngs] containing the resulting opacities

        FILE OUTPUT:
        ------------
            dustopac_i.inp - Contains the dust opacities in radmc3d format
            dustopac.inp   - Master dust opacity file

        """

#
# Calculate the grain sizes
#
        if ngs>1:
            gsize = gmin * (gmax/gmin)**(arange(ngs, dtype=float64)/(float(ngs)-1.))
        else:
            gsize = [gmin]

#
# Write the frequency.inp file
#
        wfile = open('frequency.inp', 'w')
        wfile.write("%d\n"%freq.shape[0])
        wfile.write("  \n")
        for i in range(freq.shape[0]):
            wfile.write("%.10e\n"%freq[i])
        wfile.close()

#
# Write the dust.inp file (makedust main control file)
#
        wfile = open('dust.inp', 'w')
        for igs in range(ngs):
            wfile.write("%s\n"%lnk_fname)
            wfile.write("%s\n"%"MIE")
            wfile.write("%d %f %f %f %d %f %f %f\n"%(1,0.0,log10(gsize[igs]), log10(gsize[igs]),1.,-3.5,gdens,-2.0))
        wfile.close()

#
# Run the Mie-code
#
        dum = Popen('makedust', shell=True).wait()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

class radmc3dPar():
    """
    Class for parameters in a RADMC-3D model

    ATTRIBUTES:
    -----------
        ppar   : Dictionary containing parameter values with parameter names as keys 
        pdesc  : Disctionary containing parameter description (comments in the parameter file) with parameter names as keys
        pblock : Dictionary containing the block names in the parameter file and parameter names as values 
        pvalstr: Dictionary containing parameter values as strings with parameter names as keys
    
    """

    def __init__(self):

        self.ppar = {}
        self.pdesc = {}
        self.pblock = {}
        self.pvalstr = {}
# --------------------------------------------------------------------------------------------------
    def readPar(self, fname=''):
        """
        Function to read a parameter file 
        The parameters in the files should follow the python syntax


        INPUT:
        ------
            fname  : file name to be read (if omitted problem_params.inp is used)

        OUTPUT:
        -------
            Returns a dictionary with the parameter names as keys
            
        """

        if fname=='':
            fname = 'problem_params.inp'

        try:
            rfile = open(fname, 'r')
        except:
            return 

        cchar  = '#'
        lbchar = ""

    # ------------------------------------------------------------------------------------------------------------------------
    # First read every line that is not commented (i.e. does not begin with a comment character)
    # ------------------------------------------------------------------------------------------------------------------------
        dumlist = []
        dumline = '-'

        dumline = rfile.readline()
        while dumline!='':
            # First check if the line is commented out, in which case ignore it
            comment = False
            if dumline[0]==cchar:
                if dumline.find('Block')<0:
                    comment = True

            # OK, the line is not commented out, now check if it contains a '=' sign (if not ignore the line)
            if not comment:
                # Check if we have an empty line in which case also ignore it
                if dumline.strip()!='':
                    dumlist.append(dumline)

            # Read the next line
            dumline = rfile.readline()
        
        rfile.close()

    # ------------------------------------------------------------------------------------------------------------------------
    # After every line in the file was read try to decode the lines to 
    #  [variable name] = [variable value] # [comment]
    # also try to catch if an expression has been broken into multiple lines
    # ------------------------------------------------------------------------------------------------------------------------

        varlist = []
        iline = 0
        while iline<len(dumlist):
            # First check if the line contains an '=' sign if not we have a problem 
            #  expression broken into multiple lines are should already be dealt with
            ind = dumlist[iline].find('=')
            if ind<=0:
                if dumlist[iline].find('Block')<=0:
                    print 'ERROR'
                    print ' Invalid expression in line ', iline
                    print dumlist[iline]
                    print dumlist[iline+1]
                    return
                else:
                    if dumlist[iline].find(':')<=0:
                        print 'ERROR'
                        print 'Invalid block identified'
                        print 'The syntax of the block name field is :'
                        print ' # Block : Blockname '
                        return
                    else:
                        blockname = dumlist[iline].split(':')[1].strip()

            else:
                # The line contains a '=' sign and a variable name, so let's check if the
                #  value expression is broken into multiple lines
                vlist = dumlist[iline].split('=')
                lbind = vlist[1].find('\\')
                cind  = vlist[1].find('#')

                # The line is full not broken
                if lbind==-1:
                    # Check if there is a comment field
                    if cind>=0:
                        vlist = [vlist[0], vlist[1][:cind], vlist[1][cind+1:], blockname]
                    else:
                        vlist = [vlist[0], vlist[1][:cind], ' ', blockname]
                    
                    varlist.append(vlist)
                # The value expression is broken into multiple lines; take all lines and join the pieces
                else:
                    # Check if there is any comment in the line 
                    inBrokenLine = False
                    if cind>=0:
                        # Part of the line is commented, now check if the line break is before or after the comment character
                        if lbind>cind:
                            # The line break is in the comment field so there is no real line break
                            vlist = [vlist[0], vlist[1][:cind], vlist[1][cind+1:], blockname]
                        else:
                            # The line break is before the comment character 
                            inBrokenLine = True
                            expr = vlist[1][:lbind]
                            com  = vlist[1][cind+1:]
                    else: 
                        inBrokenLine = True
                        expr  = vlist[1][:lbind]
                        com   = ' '

                    if inBrokenLine:
                        # Now gather all other pieces of this line
                        iline2 = 0
                        while inBrokenLine:
                            iline2 = iline2 + 1
                            dummy = dumlist[iline + iline2]
                            # Search for comments
                            cind2 = dummy.find('#')
                            # Search for another line break
                            lbind2 = dummy.find('\\')

    # TODO:
    # At the moment I neglect the possiblity that the second line in a broken long line begins
    # with a linebreak or commented out

                            # There is comment
                            if cind2>0:

                                # There is line break
                                if lbind2>0:
                                    # The line break is commented out
                                    if lbind2>cind:
                                        expr = expr + dummy[:cind2].strip()
                                        com  = com  + dummy[cind2+1:]
                                        inBrokenLine = False
                                    else:
                                        # The line break is not commented out
                                        expr = expr + dummy[:lbind2].strip()
                                        com  = com + dummy[cind2+1:]
                                else:
                                    #There is no line break
                                    expr = expr + dummy[:cind2].strip()
                                    com  = com  + dummy[cind2+1:]
                                    inBrokenLine = False

                            # There is no comment
                            else:
                                # There is a line break
                                if lbind2>0:
                                    expr = expr + dummy[:lbind2].strip()
                                    com  = com + dummy[cind2+1:]
                                    
                                #There is no line break
                                else:
                                    expr = expr + dummy[:cind2].strip()
                                    com  = com  + ' '
                                    inBrokenLine = False
                        iline = iline + iline2 
                        vlist = [vlist[0], expr, com, blockname]
                        varlist.append(vlist)

            iline = iline + 1
    # ------------------------------------------------------------------------------------------------------------------------
    # Now evaluate the expressions in the value field and make the final dictionary
    # ------------------------------------------------------------------------------------------------------------------------
        self.ppar = {}
        glob = globals()
        loc  = locals()
        for i in range(len(varlist)):
            try:
                val= eval(varlist[i][1], glob)
                glob[varlist[i][0].strip()] = val
            except:
                try:
                    val= eval(varlist[i][1], loc)
                    loc[varlist[i][0].strip()] = val
                except:
                    print 'Unknown expression "'+varlist[i][1]+'"'
            self.ppar[varlist[i][0].strip()] = val
            self.pvalstr[varlist[i][0].strip()] = varlist[i][1].strip()
            self.pdesc[varlist[i][0].strip()] = varlist[i][2].strip()
            self.pblock[varlist[i][0].strip()] = varlist[i][3].strip()
        return

# --------------------------------------------------------------------------------------------------
    def setPar(self,parlist=[]):
        """
        Function to add parameter to the radmc3DPar parameter class
        If the paramter is already defined its value will be modified

        INPUT:
        ------
            parlist - If the parameter is already defined parlist should be a two element
                      list 1) parameter name, 2) parameter expression/value as a string

                      If the parameter is not yet defined parlist should be a four element
                      list 1) parameter name, 2) parameter expression/value as a string
                      3) Parameter description (= comment field in the parameter file)
        """

        parname = parlist[0].strip()

        # 
        # Check whether or not the parameter is already defined
        #
        new_par = False
        if len(parlist)==2:
            if not self.ppar.keys().__contains__(parname):
                print ' ERROR'
                print ' The argument of radmc3dPar.setPar() should be a four element list if a new'
                print ' parameter is defined 1) parameter name, 2) parameter expression/value as a string'
                print ' 3) Parameter description (= comment field in the parameter file)'
                print ' 4) The name of the block in which the parameter must be placed in the problem_params.inp file'
                return
        else:
            new_par = True

        # 
        # Add the parameter to the dictionaries /change its value
        #
        glob = globals()
        loc = locals()

        try:
            self.ppar[parname] = eval(parlist[1].strip(), glob)
            glob[parname] = self.ppar[parname]
        except Exception, e:
            print e
            try:
                self.ppar[parname] = eval(parlist[1].strip(), loc)
                loc[parname] = self.ppar[parname]
            except Exception, e:
                print 'Unknown expression '+parlist[1].strip()
                print e
                return

        self.pvalstr[parname] = parlist[1].strip()
        
        if new_par:
            if not self.pdesc.has_key(parname):
                self.pdesc[parname] = parlist[2].strip()
            if len(parlist)==4:
                if not self.pblock.has_key(parname):
                    self.pblock[parname] = parlist[3].strip()


# --------------------------------------------------------------------------------------------------
    def loadDefaults(self, model='', ppar={}, reset=True):
        """
        Function to fill up the classs attributes with default values

        OPTIONS:
        ------
            model - Model name whose paraemters should also be loaded
            ppar - Dictionary containing parameter values as string and parameter names as keys
                   Default values will be re-set to the values in this dictionary

            reset - If True the all class attributes will be re-initialized before
                    the default values would be loaded. I.e. it will remove all entries
                    from the dictionary that does not conain default values either in this
                    function or in the optional ppar keyword argument
        """

        if reset:
            self.ppar = {}
            self.pvarstr = {}
            self.pdesc = {}
            self.pblock = {}

        #
        # Radiation sources
        #
        self.setPar(['mstar', '[1.0*ms]', '# Mass of the star(s)', 'Radiation sources'])
        self.setPar(['rstar','[2.0*rs]', '# Radius of the star(s)', 'Radiation sources'])
        self.setPar(['tstar','[4000.0]', '# Effective temperature of the star(s) [K]', 'Radiation sources'])
        self.setPar(['pstar','[0.0, 0.0, 0.0]', '# Position of the star(s) (cartesian coordinates)', 'Radiation sources'])
        #
        # Grid parameters
        #
        self.setPar(['crd_sys', "'sph'", '  Coordinate system used (car/cyl)', 'Grid parameters']) 
        self.setPar(['nx', '50', '  Number of grid points in the first dimension', 'Grid parameters']) 
        self.setPar(['ny', '30', '  Number of grid points in the second dimension', 'Grid parameters'])
        self.setPar(['nz', '36', '  Number of grid points in the third dimension', 'Grid parameters'])
        self.setPar(['xbound', '[1.0*au, 100.*au]', '  Boundaries for the x grid', 'Grid parameters'])
        self.setPar(['ybound', '[0.0, pi]', '  Boundaries for the y grid', 'Grid parameters'])
        self.setPar(['zbound', '[0.0, 2.0*pi]', '  Boundraries for the z grid', 'Grid parameters'])
        self.setPar(['xres_nlev', '3', 'Number of refinement levels (spherical coordinates only', 'Grid parameters'])
        self.setPar(['xres_nspan', '3', 'Number of the original grid cells to refine (spherical coordinates only)', 'Grid parameters'])
        self.setPar(['xres_nstep', '3', 'Number of grid cells to create in a refinement level (spherical coordinates only)', 'Grid parameters'])
        self.setPar(['wbound', '[0.1, 7.0, 25., 1e4]', '  Boundraries for the wavelength grid', 'Grid parameters'])
        self.setPar(['nw', '[19, 50, 30]', '  Number of points in the wavelength grid', 'Grid parameters'])

        #
        # Dust opacity
        #
        self.setPar(['lnk_fname', "['/disk2/juhasz/Data/JPDOC/astrosil/astrosil_WD2001_new.lnk', '/disk2/juhasz/Data/JPDOC/carbon/A/cel600.lnk']", ' ', 'Dust opacity'])
        self.setPar(['gdens', '[3.6, 1.8]', ' Bulk density of the materials in g/cm^3', 'Dust opacity'])
        self.setPar(['gsmin', '0.1', ' Minimum grain size', 'Dust opacity'])
        self.setPar(['gsmax', '10.0', ' Maximum grain size', 'Dust opacity'])
        self.setPar(['ngs', '1', ' Number of grain sizes', 'Dust opacity'])
        self.setPar(['gsdist_powex', '-3.5', ' Grain size distribution power exponent', 'Dust opacity'])
        self.setPar(['mixabun',       '[0.75, 0.25]', ' Mass fractions of the dust componetns to be mixed', 'Dust opacity'])
        self.setPar(['dustkappa_ext',"['silicate']", ' ', 'Dust opacity'])
        
        #
        # Gas line RT 
        #
        self.setPar(['gasspec_mol_name', "['co']", '  Name of the gas species - the extension of the molecule_EXT.inp file', 'Gas line RT'])
        self.setPar(['gasspec_mol_abun', '[1e-4]', '  Abundance of the molecule', 'Gas line RT']) 
        self.setPar(['gasspec_mol_dbase_type',"['leiden']", '  leiden or linelist', 'Gas line RT'])
        self.setPar(['gasspec_colpart_name', "['h2']", '  Name of the gas species - the extension of the molecule_EXT.inp file', 'Gas line RT'])
        self.setPar(['gasspec_colpart_abun', '[1e0]', '  Abundance of the molecule', 'Gas line RT']) 
        self.setPar(['gasspec_vturb', '0.1e5', '  Microturbulence', 'Gas line RT'])
        #self.setPar(['writeGasTemp', 'False', '  Whether or not to write a separate gas temperature file (gas_temperature.inp) if such function exists in the model', 'Gas line RT'])
        #
        # Code parameters
        #
        self.setPar(['nphot', 'long(1e5)', '  Nr of photons for the thermal Monte Carlo', 'Code parameters'])
        self.setPar(['nphot_scat','long(3e4)', '  Nr of photons for the scattering Monte Carlo (for images)', 'Code parameters'])
        self.setPar(['nphot_spec','long(1e5)', '  Nr of photons for the scattering Monte Carlo (for spectra)', 'Code parameters'])
        self.setPar(['scattering_mode_max','1', '  0 - no scattering, 1 - isotropic scattering, 2 - anizotropic scattering', 'Code parameters'])
        self.setPar(['lines_mode', '-1', '  Line raytracing mode', 'Code parameters'])
        self.setPar(['istar_sphere', '0', '  1 - take into account the finite size of the star, 0 - take the star to be point-like', 'Code parameters'])
        self.setPar(['itempdecoup', '1', '  Enable for different dust components to have different temperatures', 'Code parameters'])
        self.setPar(['tgas_eq_tdust', '1', '  Take the dust temperature to identical to the gas temperature', 'Code parameters'])
        self.setPar(['rto_style', '1', '  Format of outpuf files (1-ascii, 2-unformatted f77, 3-binary', 'Code parameters'])
        #
        # Model parameters
        #
        if model!='':
            try:
                mdl = __import__('model_'+model)
            except:
                try:
                    mdl  = __import__('radmc3dPy.model_'+model, fromlist=['']) 
                except:
                    print 'ERROR'
                    print ' model_'+model+'.py could not be imported'
                    print ' The model files should either be in the current working directory or'
                    print ' in the radmc3d python module directory'
                    return

            modpar = mdl.getDefaultParams()
            for i in range(len(modpar)):
                dum = modpar[i]
                dum.append('Model '+model)
                self.setPar(dum)
        
# --------------------------------------------------------------------------------------------------
    def printPar(self):
        """
        Print the parameters of the current model
        
        """
        
        #
        # First get the unique block names 
        #

        blocknames = ['Radiation sources', 'Grid parameters', 'Dust opacity', 'Gas line RT', 'Code parameters']
        for key in self.pblock.keys():
            dum = self.pblock[key]
            if not blocknames.__contains__(dum):
                blocknames.append(dum)

       
        #
        # Get the parameter block names and distionary keys
        #
        par_keys = self.pblock.keys()
        par_block = self.pblock.values()

        #
        # Print the parameters by blocks 
        #
        for iblock in blocknames:
            print ('%s'%'# -------------------------------------------------------------------------------------------------------------------------')
            txt = '# Block: '+iblock
            print ('%s'%txt)
            print ('%s'%'# -------------------------------------------------------------------------------------------------------------------------')
           

            keys = []
            for i in range(len(par_block)):
                if par_block[i]==iblock:
                    keys.append(par_keys[i])

            keys.sort()
            for key in keys:
                print (key.ljust(25) + ' = ' + self.pvalstr[key].strip() + '  # ' + self.pdesc[key].strip())
# --------------------------------------------------------------------------------------------------
    def writeParfile(self, fname=''):
        """
        Function to write a parameter file 


        INPUT:
        ------
            fname  : File name to be read (if omitted problem_params.inp is used)

        """
        
        if fname=='':
            fname = 'problem_params.inp'

        print 'Writing '+fname
    
        #
        # First get the uniq block names 
        #

        blocknames = ['Radiation sources', 'Grid parameters', 'Dust opacity', 'Gas line RT', 'Code parameters']
        for key in self.pblock.keys():
            dum = self.pblock[key]
            if not blocknames.__contains__(dum):
                blocknames.append(dum)

        
        try :
            wfile = open(fname, 'w')
        except:
            print ' ERROR '
            print ' Cannot create '+fname 
            return
        #
        # Write header
        #

        wfile.write('%s\n'%'###########################################################################################################################')
        wfile.write('%s\n'%'# RADMC-3D PARAMETER SETUP')
        wfile.write('%s\n'%'# Created by the python module of RADMC-3D')
        wfile.write('%s\n'%'###########################################################################################################################')
       
        #
        # Get the parameter block names and distionary keys
        #
        par_keys = self.pblock.keys()
        par_block = self.pblock.values()

        #
        # Write the parameterfile
        #
        for iblock in blocknames:
            wfile.write('%s\n'%'# -------------------------------------------------------------------------------------------------------------------------')
            txt = '# Block: '+iblock
            wfile.write('%s\n'%txt)
            wfile.write('%s\n'%'# -------------------------------------------------------------------------------------------------------------------------')
           

            keys = []
            for i in range(len(par_block)):
                if par_block[i]==iblock:
                    keys.append(par_keys[i])

            keys.sort()
            for key in keys:
                wfile.write(key.ljust(25) + ' = ' + self.pvalstr[key].strip() + '  # ' + self.pdesc[key].strip() + '\n')

# --------------------------------------------------------------------------------------------------
# Functions for an easy compatibility with the IDL routines
# --------------------------------------------------------------------------------------------------
def readOpac(ext=[''], idust=None):
    """
    Function to read the dust opacity files 
    This function is an interface to radmc3dDustOpac.readOpac()

    INPUT:
    ------
        ext  : file name extension (file names should look like 'dustkappa_ext.inp')
        idust: index of the dust species in the master opacity file (dustopac.inp')

    OUTPUT:
    -------
        Returns an instance of the radmc3dDustOpac class with the following attributes:
        
        wav     - wavelength grid
        freq    - frequency grid
        nwav    - number of wavelengths
        kabs    - absorption coefficient per unit mass
        ksca    - scattering coefficient per unit mass
        phase_g - phase function
        ext     - if set it contains the file name extension of the duskappa_ext.Kappa file
        therm   - if False the dust grains are quantum-heated (default: True)
        idust   - index of the dust species in the dust density distribution array
    
    """


    res = radmc3dDustOpac()
    res.readOpac(ext=ext, idust=idust)
    
    return res
# --------------------------------------------------------------------------------------------------
# Functions for an easy compatibility with the IDL routines
# --------------------------------------------------------------------------------------------------
def readData(ddens=False, dtemp=False, gdens=False, gtemp=False, gvel=False, ispec=None, vturb=False, binary=True):
    """
    Function to read the model data (e.g. density, velocity, temperature)

    INPUT:
    ------
        ddens - If True dust density will be read (all dust species and grain sizes)
        dtemp - If True dust temperature will be read (all dust species and grain sizes)
        gdens - If True gas density will be read (NOTE: the gas density will be number density in 1/cm^3)
        gtemp - If True gas temperature will be read (all dust species and grain sizes)
        gvel  - If True the velocity field will be read
        ispec - Name of the molecule in the 'molecule_ispec.inp' filename

    OUTPUT:
    ------
        Returns an instance of the radmc3dData class with the following attributes:
            rhodust   -  Dust density in g/cm^3 
            dusttemp  -  Dust temperature in K 
            rhogas    -  Gas density in molecule/cm^3
            gasvel    -  Gas velocity in cm/s 
            gastemp   -  Gas temperature in K
            vturb     -  Mictroturbulence in cm/s
            taux      -  Optical depth along the x (cartesian) / r (cylindrical) / r (spherical) dimension
            tauy      -  Optical depth along the y (cartesian) / theta (cylindrical) / theta (spherical) dimension
            tauz      -  Optical depth along the z (cartesian) / z (cylindrical) / phi (spherical) dimension
            sigmadust -  Dust surface density in g/cm^2
            sigmagas  -  Gas surface density in molecule/cm^2 (or g/cm^2 depending on the dimension of rhogas)
    """

    res = radmc3dData()
    if ddens: res.readDustDens(binary=binary)
    if dtemp: res.readDustTemp(binary=binary)
    if gvel: res.readGasVel(binary=binary)
    if gtemp: res.readGasTemp(binary=binary)
    if vturb: res.readVTurb(binary=binary)
    if gdens:
        if not ispec:
            print 'ERROR'
            print 'No gas species is specified!'
            print 'The ispec input keyword should be set to the name of the gas species as it appears in '
            print ' numberdens_gasspecname.inp'
            return 0
        else:
            res.readGasDens(ispec=ispec,binary=binary)

    return res

# --------------------------------------------------------------------------------------------------
def readGrid():
    """
    Function to read the spatial and frequency grid

    OUTPUT
    ------

        Returns an instance of the radmc3dGrid class with the following attributes:

        crd_sys    - 'car'/'cyl'/'sph' coordinate system of the spatial grid
        act_dim    - A three element vector the i-th element is 1 if the i-th dimension is active, otherwize the i-th element is zero
        nx         - Number of grid points in the x (cartesian) / r (cylindrical) / r (spherical) dimension
        ny         - Number of grid points in the y (cartesian) / theta (cylindrical) / theta (spherical) dimension
        nz         - Number of grid points in the z (cartesian) / z (cylindrical) / phi (spherical) dimension
        nxi        - Number of cell interfaces in the x (cartesian) / r (cylindrical) / r (spherical) dimension
        nyi        - Number of cell interfaces in the y (cartesian) / theta (cylindrical) / theta (spherical) dimension
        nzi        - Number of cell interfaces in the z (cartesian) / z (cylindrical) / phi (spherical) dimension
        nwav       - Number of wavelengths in the wavelength grid
        freq       - Number of frequencies in the grid (equal to nwav)
        x          - Cell centered x (cartesian) / r (cylindrical) / r (spherical)  grid points
        y          - Cell centered y (cartesian) / theta (cylindrical) / theta (spherical)  grid points
        z          - Cell centered z (cartesian) / z (cylindrical) / phi (spherical)  grid points
        xi         - Cell interfaces in the x (cartesian) / r (cylindrical) / r (spherical)  dimension
        yi         - Cell interfaces in the y (cartesian) / theta (cylindrical) / theta (spherical)  dimension
        zi         - Cell interfaces in the z (cartesian) / z (cylindrical) / phi (spherical)  dimension
        wav        - Wavelengh  grid
        freq       - Frequency  grid
    """

    grid = radmc3dGrid()
    grid.readGrid()

    return grid

# --------------------------------------------------------------------------------------------------
def readParams():
    """
    Function to read the problem_params.inp file (interface function to radmc3dPar.readPar())

    OUTPUT:
    -------
        Returns an instance of the radmc3dPar class with the following attributes:

        ppar   : Dictionary containing parameter values with parameter names as keys 
        pdesc  : Disctionary containing parameter description (comments in the parameter file) with parameter names as keys
        pblock : Dictionary containing the block names in the parameter file and parameter names as values 
        pvalstr: Dictionary containing parameter values as strings with parameter names as keys
    """

    dum = radmc3dPar()
    dum.readPar()
    return dum
# --------------------------------------------------------------------------------------------------
def writeDefaultParfile(model='', fname=''):
    """
    Function to write a parameter file (problem_params.inp) with default parameters for a given model

    INPUT:
    ------
        model - Name of the model whose parameter should be written to the file

    OPTIONS:
    --------
        fname - Name of the parameter file to be written (if omitted problem_params.inp will be used)


    """
    
    if model=='':
        print ' ERROR '
        print ' No model name is given '
        return

    dum  = radmc3dPar()
    dum.loadDefaults(model=model)
    dum.writeParfile()
# --------------------------------------------------------------------------------------------------
def readSpectrum(fname=''):
    """
    Function to read the spectrum / SED


    OPTIONS:
    --------
        fname - Name of the file to be read


    OUTPUT:
    -------
        Returns a two dimensional Numpy array with [Nwavelength, 2] dimensions 
        [Nwavelength,0] is the wavelength / velocity and
        [Nwavelength,1] is the flux density
        
    """
   
    if fname.strip()=='':
        fname = 'spectrum.out'

    
    rfile = open(fname, 'r')
    # Read the format number
    dum = rfile.readline()
    # Read the number of wavelengths 
    nwav = int(rfile.readline())
    # Read a blank line
    dum = rfile.readline()
    
    res = zeros([nwav, 2], dtype=float64)
    for iwav in range(nwav):
        dum = rfile.readline().split()
        res[iwav,0] = float(dum[0])
        res[iwav,1] = float(dum[1])
    rfile.close()
    return res

# --------------------------------------------------------------------------------------------------
def getDensVstruct(data=None, vmean_temp=False, ispec_tgas=0, gsize=[], idust=None, mstar=0.):
    """
    Calculates the vertical hydrostatic equilibrium

    INPUT:
    ------
        data        - An instance of the radmc3DData class
        vmean_temp  - If True (T(z) = T(-z) = 0.5*(T(z) + T(-z))) if False (T(z)!=T(-z)) 
        idust       - List of dust indices whose structure must be calculated
        mstar       - Stellar mass
    
    OPTIONS:
    --------
        ispec_tgas  - Index of dust species whose temperature is taken to be the gas temperature
        gsize       - Dust grain sizes - If specified, the gas temperature is calculated as the average temperature
                      of all dust grains in the grid cell weighted by the total surface area of dust grains with given
                      size - NOTE: this approach assumes that all dust grains of a given size have the same bulk density

    OUTPUT:
    -------
        Returns a Numpy array with the dust density
    """
        
    # Fix the mean molecular weight to 2.3
    mu = 2.3

    # Pre-calculate some constants
    A  = mu*mp*gg*mstar / kk
    cost  = cos(data.grid.y)
    costi = cos(data.grid.yi)

    if not mstar:
        print 'ERROR'
        print ' You should specify the stellar mass (mstar = ??)'
        return

    if idust==None:
        print ' No dust index was given for which the vertical structure should be calculated'
        print ' So we do for all dust species'
        idust = range(data.rhodust.shape[3])
    else:
        if (type(idust).__name__=='int')| (type(idust).__name__=='float'):
            idust = [idust]
    # To improve the smoothness of the temperature structure, if the density structure is
    #  symmetric to the disk midplane we use T_new(theta) = T_new(pi-theta) = 0.5 * (T(theta) + T(pi-theta))
    if vmean_temp:

        if abs(data.grid.yi[data.grid.nyi-1]-pi/2.)<1e-8:
            print 'ERROR'
            print "Cannot average temperature in the vertical direction if theta mirroring is active"
            return None
        else:
            print ' Smoothing the vertical temperature structure by averaging the temperature of the '
            print " two half planes above and below the disk midplane"
            dusttemp = zeros(data.dusttemp.shape, dtype=float64)
            for iy in range(data.grid.ny/2):
                print iy
                dusttemp[:,iy,:,:] = 0.5 * (data.dusttemp[:,iy,:,:] + data.dusttemp[:,data.grid.ny-1-iy,:,:])
                dusttemp[:,data.grid.ny-1-iy,:,:] = dusttemp[:,iy,:,:]
    # Calculate the vertical hydrostatic equilibrium for the two half space (z<0, z>0) separately
    else:
        dusttemp = data.dusttemp
      
    #rho_new = zeros(data.rhodust.shape, dtype=float64)
    rho_new = array(data.rhodust)
    if len(gsize)!=0:
        mean_dusttemp = zeros([data.grid.nx, data.grid.ny, data.grid.nz], dtype=float64)
        w             = zeros(data.rhodust.shape, dtype=float64)
        for ispec in idust:
            w[:,:,:,ispec] = gsize[ispec]**2 * (data.rhodust[:,:,:,ispec] / gsize[ispec]**3) 
        
        wnorm = w.sum(3)
        for ispec in idust:
            w[:,:,:,ispec] = w[:,:,:,ispec]/wnorm

        for ispec in idust:
            mean_dusttemp = mean_dusttemp + data.dusttemp[:,:,:,ispec] * w[:,:,:,ispec]

    # Loop over all dust species where we should calculate the vertical structure
    for ispec in idust:
        rho_new[:,:,:,ispec] = 0.
        for ir in range(data.grid.nx):
            print ir, data.grid.nx-1
            r     = data.grid.x[ir]
            z     = r * cost
            zi    = r * costi
            dz    = z[:-1] - z[1:]
            const = A / r**3

            # Do we have theta mirroring active?
            if abs(data.grid.yi[data.grid.nyi-1]-pi/2.)<1e-8:
                for ip in range(data.grid.nz):
                    dlgrho  = log(data.rhodust[ir,1:,ip,ispec]) - log(data.rhodust[ir,:-1,ip,ispec])
                    if len(gsize)!=0:
                        temp    = mean_dusttemp[ir,:,ip] 
                    else:
                        temp    = dusttemp[ir,:,ip,ispec]
                    
                    it = data.grid.ny-1
                    temp[it] = 0.5 * (temp[it] + temp[it-1])    

                    dlgtemp = log(temp[1:]) - log(temp[:-1])
                    zpt     = z/temp
                    zpt     = 0.5*(zpt[1:] + zpt[:-1])


                    # Calculate the normalized (rho[z=0] = 1.0) density
                    rho_new[ir,data.grid.ny-1,ip,ispec] = 1.0

                    for it in range(data.grid.ny-1)[::-1]:
                        rho_new[ir,it,ip,ispec] = rho_new[ir,it+1,ip,ispec] * exp(-(const*zpt[it] + dlgtemp[it]/dz[it])*dz[it])
                
                    rho_new = rho_new.clip(1e-90, 1e90)
                  
                    # Now re-normalize the surface density to the input value
                    sigma = (data.rhodust[ir,:,ip,ispec] * (zi[1:] - zi[:-1])).sum()
                    sigma_new = (rho_new[ir,:,ip,ispec] * (zi[1:] - zi[:-1])).sum()

                    rho_new[ir,:,ip,ispec] = rho_new[ir,:,ip,ispec] * sigma / sigma_new

            else:
                for ip in range(data.grid.nz):
                    dlgrho  = log(data.rhodust[ir,1:,ip,ispec]) - log(data.rhodust[ir,:-1,ip,ispec])
                    if len(ispec_weights)!=0:
                        temp    = (dusttemp[ir,:,ip,ispec]*ispec_weights).sum()
                    else:
                        temp    = dusttemp[ir,:,ip,ispec]
                    dlgtemp = log(temp[1:]) - log(temp[:-1])
                    zpt     = z/temp
                    zpt     = 0.5*(zpt[1:] + zpt[:-1])

                    # Calculate the normalized (rho[z=0] = 1.0) density
                    rho_new[ir,data.grid.ny/2-1,ip,ispec] = 1.0
                    rho_new[ir,data.grid.ny/2,ip,ispec] = 1.0

                    for it in range(data.grid.ny/2)[::-1]:
                        rho_new[ir,it-1,ip,ispec] = rho_new[ir,it,ip,ispec] * exp(-(const*zpt[it] + dlgtemp[it]/dz[it])*dz[it])
                    for it in range(data.grid.ny/2, data.grid.ny-1)[::1]:
                        rho_new[ir,it,ip,ispec] = rho_new[ir,it-1,ip,ispec] * exp((const*zpt[it-1] + dlgtemp[it-1]/dz[it-1])*dz[it-1])
                   
                    rho_new = rho_new.clip(1e-90, 1e90)

                    # Now re-normalize the surface density to the input value
                    sigma = (data.rhodust[ir,:,ip,ispec] * (zi[1:] - zi[:-1])).sum()
                    sigma_new = (rho_new[ir,:,ip,ispec] * (zi[1:] - zi[:-1])).sum()

                    rho_new[ir,:,ip,ispec] = rho_new[ir,:,ip,ispec] * sigma / sigma_new
            
            print rho_new[ir,data.grid.ny/2-1,ip,ispec]

    return rho_new

