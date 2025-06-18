import numpy as np
import netCDF4 as nc
from scipy.interpolate import RegularGridInterpolator

#Constants
Re = 6371.0 # km
emission = 568.4900 # eV
abundance = 1.1236e-04 # dimensionless

class Cube:
    def __init__(self, cube_name, step):
        self.step = step
        self.ds = step * Re * 10**5 # step in cm
        self.name = cube_name
        # Define GSE coordinate grid
        self.x_GSE = np.arange(- 0.05, 20, self.step)   # 20 down to 0
        self.y_GSE = np.arange(-10.05, 10,  self.step)   # 10 down to -10
        self.z_GSE = np.arange(-10, 10.05,  self.step)   # -10 up to 10
        [self.X_GSE, self.Y_GSE, self.Z_GSE] = np.meshgrid(self.x_GSE,self.y_GSE,self.z_GSE)

    def get_Production_GSO(self):  
        '''
        Load NetCDF cube that is saved in the ./Data folder
        '''
        load_path = f"./Data/{self.name}.nc"
        cube_data = nc.Dataset(load_path,'r')
        Prod_GSO = cube_data.variables['Prod'][:]  # order [Z,Y,X] 1/cm^3/s
        cube_data.close()
        # Trasnpose Prog_GSO to get [X,Y,Z] format
        self.Prod_GSO = np.transpose(Prod_GSO, (1, 2, 0)) 

    def get_Production_GSE(self):
        # Build interpolation function to extract the values of Prod_GSO in the GSE grid
        f_interp = RegularGridInterpolator(
                (self.x_GSE, self.y_GSE, self.z_GSE),
                self.Prod_GSO,
                method='linear',
                bounds_error=False,
                fill_value=None
            )
        # Flatten points to interpolate. 
        pts = np.array([self.X_GSE.ravel(), self.Y_GSE.ravel(), self.Z_GSE.ravel()]).T
        Prod_GSE_flat = f_interp(pts)
        self.Prod_GSE = Prod_GSE_flat.reshape(self.X_GSE.shape) 

    def Q_GSE(self):
        self.get_Production_GSO()
        self.get_Production_GSE()
        self.Q_GSE =  self.Prod_GSE * emission * abundance





'''
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

'''