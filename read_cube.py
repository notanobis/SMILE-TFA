import numpy as np
import netCDF4 as nc
from scipy.interpolate import RegularGridInterpolator
from scipy.io import savemat
from matplotlib import pyplot as plt
from matplotlib import ticker

#Constants
Re = 6371.0 # km
emission = 568.4900 # eV
abundance = 1.1236e-04 # dimensionless
ds = 0.05 * Re * 10**5 # step in cm

def read(name):
    load_path = f"./Data/{name}.nc"
    data = nc.Dataset(load_path,'r')
    print("The available variables are:", data.variables.keys())
    return data

def cart2pol(x, y, center):
    x = x - center[0]
    y = y - center[1]
    rho = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x)
    return(rho, theta)

def pol2cart(rho, theta):
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return(x, y)

def readProdGSO(filename):
    # Read data
    ds = read(filename)
    Prod_GSO = ds.variables['Prod'][:]  # shape (nz, ny, nx) 1/cm^3/s
    ds.close()
    # Trasnpose Prog_GSO to get [X,Y,Z] format
    Prod_GSO = np.transpose(Prod_GSO, (1, 2, 0)) 
    return Prod_GSO

def GSOtoGSE(Prod_GSO):
    # Build GSE grid
    x_GSOatGSE = np.arange(20,  - 0.05, -0.05)   # 20 down to 0
    y_GSOatGSE = np.arange(10, -10.05,  -0.05)   # 10 down to -10
    z_GSOatGSE = np.arange(-10, 10.05,  0.05)   # -10 up to 10
    x_GSE = np.sort(x_GSOatGSE)
    y_GSE = np.sort(y_GSOatGSE)
    z_GSE = np.sort(z_GSOatGSE)
    [X_GSE, Y_GSE, Z_GSE] = np.meshgrid(x_GSE,y_GSE,z_GSE)

    # Build interpolation function to extract the values of Prod_GSO in the GSE grid
    f_interp = RegularGridInterpolator(
            (x_GSOatGSE, y_GSOatGSE, z_GSOatGSE),
            Prod_GSO,
            method='linear',
            bounds_error=False,
            fill_value=None
        )
    
    # Flatten points to interpolate. 
    pts = np.array([X_GSE.ravel(), Y_GSE.ravel(), Z_GSE.ravel()]).T
    P_GSE_flat = f_interp(pts)
    P_GSE = P_GSE_flat.reshape(X_GSE.shape) 
    return x_GSE,y_GSE,z_GSE,P_GSE


def GSOtoGSE_Vector(Vector_GSO):
    # Build GSO grid as before
    x_GSOatGSE = np.arange(20, -0.05, -0.05)
    y_GSOatGSE = np.arange(10, -10.05, -0.05)
    z_GSOatGSE = np.arange(-10, 10.05, 0.05)
    
    x_GSE = np.sort(x_GSOatGSE)
    y_GSE = np.sort(y_GSOatGSE)
    z_GSE = np.sort(z_GSOatGSE)

    [X_GSE, Y_GSE, Z_GSE] = np.meshgrid(x_GSE, y_GSE, z_GSE)

    # Prepare flattened points for interpolation
    pts = np.array([X_GSE.ravel(), Y_GSE.ravel(), Z_GSE.ravel()]).T

    # Interpolate each component of the vector field
    Vx_interp = RegularGridInterpolator((x_GSOatGSE, y_GSOatGSE, z_GSOatGSE), Vector_GSO[0], bounds_error=False, fill_value=None)
    Vy_interp = RegularGridInterpolator((x_GSOatGSE, y_GSOatGSE, z_GSOatGSE), Vector_GSO[1], bounds_error=False, fill_value=None)
    Vz_interp = RegularGridInterpolator((x_GSOatGSE, y_GSOatGSE, z_GSOatGSE), Vector_GSO[2], bounds_error=False, fill_value=None)

    Vx_GSE = Vx_interp(pts).reshape(X_GSE.shape)
    Vy_GSE = Vy_interp(pts).reshape(X_GSE.shape)
    Vz_GSE = Vz_interp(pts).reshape(X_GSE.shape)

    return X_GSE, Y_GSE, Z_GSE, Vx_GSE, Vy_GSE, Vz_GSE

def Q_GSE(Prog_GSE):
    return Prog_GSE * emission * abundance

def shue_parameters(bz, vx, n):
    """
    Computes the magnetopause standoff distance (r0) and flaring parameter (alpha)
    using the Shue et al. (1997) model.

    Parameters:
    bz : float
        Interplanetary magnetic field Bz component (nT).
    vx : float
        Solar wind velocity (km/s).
    n  : float
        Solar wind density (cm^-3).

    Returns:
    r0 : float
        Magnetopause standoff distance (Re).
    alpha : float
        Flaring parameter.
    """
    dp = 1.15 * n * vx**2 * 1.67e-6  # Dynamic pressure in nPa

    if bz >= 0:
        r0 = (11.4 + 0.013 * bz) * dp**(-1 / 6.6)
    else:
        r0 = (11.4 + 0.14 * bz) * dp**(-1 / 6.6)

    alpha = (0.58 - 0.01 * bz) * (1 + 0.01 * dp)

    return r0, alpha

def shue_polar(r0, alpha, theta):
    return r0 * (2 / (1 + np.cos(theta)))**alpha

