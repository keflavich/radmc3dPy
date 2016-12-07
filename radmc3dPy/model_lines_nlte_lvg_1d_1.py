"""
PYTHON module for RADMC3D 
(c) Attila Juhasz, Kees Dullemond 2011,2012,2013,2014

Original IDL model by Kees Dullemond, Python translation by Attila Juhasz
"""
try:
    import numpy as np
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'

from radmc3dPy.natconst import *



# ============================================================================================================================
#
# ============================================================================================================================
def getModelDesc():
    """
    Function to provide a brief description of the model
    """

    return "Example model: A 1D simple velocity gradient model to calculate lines with the LVG method"
           

# ============================================================================================================================
#
# ============================================================================================================================
def getDefaultParams():
    """
    Function to provide default parameter values 

    OUTPUT:
    -------

    Returns a list whose elements are also lists with three elements:
    1) parameter name, 2) parameter value, 3) parameter description
    All three elements should be strings. The string of the parameter
    value will be directly written out to the parameter file if requested,
    and the value of the string expression will be evaluated and be put
    to radmc3dData.ppar. The third element contains the description of the
    parameter which will be written in the comment field of the line when
    a parameter file is written. 
    """

    defpar = {}

    defpar = [['mstar', '1.0*ms', 'Mass of the star(s)'],
    ['pstar', '[0., 0., 0.]', 'Position of the star(s) (cartesian coordinates)'],
    ['rstar', '1.0*rs', 'Radius of the star(s)'],
    ['tstar', '1.0*ts', 'Effective temperature of the star(s)'],
    ['crd_sys', "'car'", 'Coordinate system used (car/sph)'],
    ['nx', '10', 'Number of grid points in the first dimension'],
    ['ny', '1', 'Number of grid points in the second dimension'],
    ['nz', '1', 'Number of grid points in the third dimension'],
    ['xbound', '[-1000.0*au, 1000.0*au]', 'Boundaries for the x-grid'], 
    ['ybound', '[-1000.0*au/nx, 1000.0*au/nx]', 'Boundaries for the y-grid'], 
    ['zbound', '[-1000.0*au/nx, 1000.0*au/nx]', 'Boundaries for the z-grid'], 
    ['nw', '[20,100,30]', 'Number of points in the wavelength grid'],
    ['wbound', '[0.1, 7., 25., 1e4]', 'Boundaries for the wavelength grid'],
    ['dustkappa_ext', "['silicate']", 'Dust opacity file name extension'],  
    ['nphot', '1000000', 'Number of photons in the thermal Monte Carlo simulation'],
    ['lines_mode', '3', ''],
    ['scattering_mode_max', '1', '0 - no scattering, 1 - isotropic scattering, 2 - anizotropic scattering'],
    ['gasspec_mol_name', "['co']", ''],
    ['gasspec_mol_abun', '[1e-4]', ''],
    ['gasspec_mol_dbase_type', "['leiden']", ''],
    ['gasspec_colpart_name', "['h2']", ''],
    ['gasspec_colpart_abun', '[1e0]', ''],
    ['gasspec_vturb', '1e5', ''],
    ['abun_h2', '0.5', ''],
    ['abun_he', '0.1', ''],
    ['nh2', '1e5', ''],
    ['temp0', '30.', ''],
    ['tdust0', '30.', ''],
    ['dusttogas', '1e-2', ''],
    ['dvdau', '1e-2*1e5', '']]

    return defpar

# ============================================================================================================================
#
# ============================================================================================================================
def getGasTemperature(grid=None, ppar=None):
    """
    Function to calcualte/set the gas temperature
    
    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model 
    
    OUTPUT:
    -------
        returns the gas temperature in K
    """


    tgas = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + ppar['temp0']
    return tgas
# ============================================================================================================================
#
# ============================================================================================================================
def getDustTemperature(grid=None, ppar=None):
    """
    Function to calcualte/set the dust temperature
    
    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model 
    
    OUTPUT:
    -------
        returns the dust temperature in K
    """


    tdust = np.zeros([grid.nx, grid.ny, grid.nz, 1], dtype=np.float64) + ppar['tdust0']
    return tdust
# ============================================================================================================================
#
# ============================================================================================================================
def getGasAbundance(grid=None, ppar=None, ispec=''):
    """
    Function to create the conversion factor from volume density to number density of molecule ispec.
    The number density of a molecule is rhogas * abun 
   
    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model 
        ispec - The name of the gas species whose abundance should be calculated

    OUTPUT:
    -------
        returns the abundance as a Numpy array
    """
    
    # Mass of gas per H2-molecule
    #mgas    = mp*(2.0*ppar['abun_h2']+4*ppar['abun_he'])/ppar['abun_h2']

    gasabun = -1
    if ppar['gasspec_mol_name'].__contains__(ispec):
        gasabun = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) 
        ind = ppar['gasspec_mol_name'].index(ispec)
        gasabun[:,:,:] = ppar['gasspec_mol_abun'][ind]#/mgas
 
    elif ppar['gasspec_colpart_name'].__contains__(ispec):
        gasabun = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) 
        ind = ppar['gasspec_colpart_name'].index(ispec)
        gasabun[:,:,:] = ppar['gasspec_colpart_abun'][ind]#/mgas
    else:
        print 'ERROR'
        print ' The abundance of "'+ispec+'" is not specified in the parameter file'
   
    return gasabun

# ============================================================================================================================
#
# ============================================================================================================================
def getGasDensity(grid=None, ppar=None):
    """
    Function to create the total gas density distribution 
    
    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model 
    
    OUTPUT:
    -------
        returns the volume density in g/cm^3
    """
    # Mass of gas per H2-molecule
    mgas    = mp*(2.0*ppar['abun_h2']+4*ppar['abun_he'])/ppar['abun_h2']
    rhogas = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + ppar['nh2'] * mgas
    return rhogas

# ============================================================================================================================
#
# ============================================================================================================================
def getDustDensity(grid=None, ppar=None):
    """
    Function to create the dust density distribution 
    
    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model 
    
    OUTPUT:
    -------
        returns the volume density in g/cm^3
    """

    rhogas  = getGasDensity(grid=grid, ppar=ppar)
    rhodust = np.zeros([grid.nx, grid.ny, grid.nz, 1], dtype=np.float64) 
    rhodust[:,:,:,0] = rhogas * ppar['dusttogas']
    return rhodust

# ============================================================================================================================
#
# ============================================================================================================================
def getVTurb(grid=None, ppar=None):
    """
    Function to create the turbulent velocity field
    
    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model 
    
    OUTPUT:
    -------
        returns the turbulent velocity in cm/s
    """

    vturb = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64) + ppar['gasspec_vturb']
    return vturb

# ============================================================================================================================
#
# ============================================================================================================================
def getVelocity(grid=None, ppar=None):
    """
    Function to create the turbulent velocity field
    
    INPUT:
    ------
        grid - An instance of the radmc3dGrid class containing the spatial and wavelength grid
        ppar - Dictionary containing all parameters of the model 
    
    OUTPUT:
    -------
        returns the turbulent velocity in cm/s
    """

    vel = np.zeros([grid.nx, grid.ny, grid.nz, 3], dtype=np.float64)
    for ix in range(grid.nx):
        vel[ix,:,:,0] = ppar['dvdau']*grid.x[ix]/au
    
    return vel



