"""
PYTHON module for RADMC3D 
(c) Attila Juhasz 2011,2012,2013

This sub-module contains functions for coordinate transformations (e.g. rotation)
"""
try:
    import numpy as np
except:
    print 'ERROR'
    print ' Numpy cannot be imported '
    print ' To use the python module of RADMC-3D you need to install Numpy'


#import scipy.linalg as la

def ctrans_sph2cyl(crd=None, theta=None, reverse=False):
    """
    Function to transform coordinates between spherical to cylindrical systems

    INPUT : 
    -------
            r,phi,theta : numpy arrays containing the spherical coordinates

    OPTIONS : 
    ----------
            reverse=False : Calculates the inverse trasnformation
                       (cartesian -> spherical). In this case crd should be [r,phi,theta]

    OUTPUT : 
    --------
            result   : a numpy array of [Nr,Nphi,Ntheta,3] dimensions containing the cylindrical 
                       coordinates [rcyl, z, phi]
    """
    
    nr     = r.shape[0]
    nt     = theta.shape[0]
    result = np.zeros([nr, nt], dtype=np.float64)
    for ix in range(grid.nr):
        result[ix,:,0] = r[ix] * np.sin(theta)
        result[ix,:,1] = r[ix] * np.cos(theta)

    return result


def ctrans_sph2cart(crd=[0,0,0], reverse=False):
    """
    Function to transform coordinates between spherical to cartesian systems

    INPUT : 
    -------
            crd      : Three element array containing the input
                       coordinates [x,y,z] or [r,phi,theta] by default
                       the coordinates assumed to be in the cartesian system
    OPTIONS :
    ---------
            reverse=False : Calculates the inverse trasnformation
                       (cartesian -> spherical). In this case crd should be [r,phi,theta]

    OUTPUT : 
    --------
            result   : A three element array containg the output
                       coordinates [r,phi,theta] or [x,y,z]
    """
    if (reverse==False):
        r     = crd[0]
        phi   = crd[1]
        theta = crd[2] + 1e-50

        x = np.sin(theta) * np.cos(phi) * r
        y = np.sin(theta) * np.sin(phi) * r
        z = np.cos(theta) * r

        crdout = [x, y, z]

    else:

        x = crd[0]
        y = crd[1]
        z = crd[2]

        r     = np.sqrt(x**2 + y**2 + z**2)
        phi   = np.arccos(x / np.sqrt(x**2 + y**2) + 1e-90)
        theta = np.arccos(z / r)

        if (y<0.): phi = 2.0 * np.pi - phi
        
        crdout = [r, phi, theta]
        

    return crdout

def vtrans_sph2cart(crd=[0,0,0], v=[0,0,0], reverse=False):
    """
    Function to transform velocities between spherical to cartesian systems

    INPUT : 
    -------
            crd      : Three element array containing the input
                       coordinates [x,y,z] or [r,phi,theta] by default
                       the coordinates assumed to be in the cartesian system

            v        : Three element array containing the input
                       velocities in the same coordinate system as crd


    OPTIONS :
    ---------
            reverse=False : Calculates the inverse trasnformation (cartesian -> spherical)

    OUTPUT : 
    --------
            result   : A three element array containg the output
                       velocities [vr,vphi,vtheta] or [vx,vy,vz]


    NOTE!!!!! The velocities in the spherical system are not angular velocities!!!!
    v[1] = dphi/dt * r
    v[2] = dtheta/dt * r
    """
    if (reverse==False):
        r      = crd[0]
        phi    = crd[1]
        theta  = crd[2]
        
        vr     = v[0]
        vphi   = v[1]
        vtheta = v[2]
        
        vx     = vr*np.sin(theta)*np.cos(phi) - vphi*np.sin(phi) + vtheta*np.cos(theta)*np.cos(phi)
        vy     = vr*np.sin(theta)*np.sin(phi) + vphi*np.cos(phi) + vtheta*np.cos(theta)*np.sin(phi)
        vz     = vr*np.cos(theta) - vtheta*np.sin(theta)

        vout   = [vx,vy,vz]

    else:
        
        crd_sph = ctrans_sph2cart(crd, reverse=True)
        r       = crd_sph[0]
        phi     = crd_sph[1]
        theta   = crd_sph[2]
        
        a       = [[np.sin(theta)*np.cos(phi), -np.sin(phi), np.cos(theta)*np.cos(phi)],\
                   [np.sin(theta)*np.sin(phi), np.cos(phi), np.cos(theta)*np.sin(phi)],\
                   [np.cos(theta), 0., -np.sin(theta)]]

        a       = np.array(a, dtype=np.float64)

        vout = la.solve(a,v)

    return vout

def csrot(crd=None, ang=None, xang=0.0, yang=0.0, zang=0.0, deg=False):
    """
    Function to make coordinate system rotation
 
    INPUT : 
    -------
           crd  : three element vector containing the coordinates of a
                  given point in a cartesian system

           ang  : three element array, angles of rotation around the x,y,z axes

    OPTIONS : 
    ---------
           deg=True : if this keyword is set angles should be given in
                  angles instead of radians (as by default)
 
    Rotation matrices :
    -------------------
    X-axis

     |      1               0            0        |
     |      0          cos(alpha)    -sin(alpha)  | 
     |      0          sin(alpha)     cos(alpha)  |

    Y-axis

     |   cos(beta)          0         sin(beta)   |
     |      0               1            0        |
     |   -sin(beta)         0         cos(beta)   |

    Z-axis
 
     |   cos(gamma)     -sin(gamma)       0        |
     |  sin(gamma)       cos(gamma)       0        |
     |      0                0            1        |

    """

    crd_new = np.zeros(len(crd), dtype=np.float64)

    if (ang!=None):
        xang = ang[0]
        yang = ang[1]
        zang = ang[2]

#
# Convert degree into radian if the angles are given in degree
#

    if (deg==True):
        xang = xang / 180.0 * np.pi
        yang = yang / 180.0 * np.pi
        zang = zang / 180.0 * np.pi

#
# Rotation around the x axis
#

    if (xang!=0.0):
        dumx = crd[0]
        dumy = np.cos(xang)*crd[1]  - np.sin(xang)*crd[2]
        dumz = np.sin(xang)*crd[1] + np.cos(xang)*crd[2] 
    
        crd_new = [dumx, dumy, dumz]

#
# Rotation around the y axis
#

    if (yang!=0.0):
        dumx = np.cos(yang)*crd[0] + np.sin(yang)*crd[2]
        dumy = crd[1]
        dumz = -np.sin(yang)*crd[0] + np.cos(yang)*crd[2]
        
        crd_new = [dumx, dumy, dumz] 

#
# Rotation around the z axis
#
  
    if (zang!=0.0):
        dumx = np.cos(zang)*crd[0] - np.sin(zang)*crd[1] + 0.0
        dumy = np.sin(zang)*crd[0] + np.cos(zang)*crd[1] + 0.0
        dumz = crd[2]
        
        crd_new = [dumx, dumy, dumz]

    return crd_new

def vrot(crd=None, v=None, ang=None):
    """
    Function to rotate a vector in spherical coordinate system

    First transform the vector to cartesian coordinate system do the rotation then make the
     inverse transformation

    INPUT : 
    -------
           crd  : three element vector containing the coordinates of a
                  given point in a cartesian system

           v    : three element array, angles of rotation around the x,y,z axes

           ang  : angle around the x, y, z, axes with which the vector should be rotated

    """
# Convert the position vector to cartesian coordinate system
    crd_xyz   = ctrans_sph2cart(crd=crd)
# Convert the velocity vector to cartesian coordinate system
    v_xyz     = vtrans_sph2cart(crd=crd, v=v)
# Rotate the vecto
    v_xyz_rot = csrot(crd=v_xyz, ang=ang)
# Transform the rotated vector back to the spherical coordinate system
    v_rot     = vtrans_sph2cart(crd=crd_xyz, v=v_xyz_rot, reverse=True) 
    
    return v_rot
