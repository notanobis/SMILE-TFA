import numpy as np
import netCDF4 as nc
from scipy.interpolate import RegularGridInterpolator
from matplotlib import pyplot as plt

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
        self.x_GSOatGSE = np.arange(20,  - 0.05, -0.05)   # 20 down to 0
        self.y_GSOatGSE = np.arange(10, -10.05,  -0.05)   # 10 down to -10
        self.z_GSOatGSE = np.arange(-10, 10.05,  0.05)   # -10 up to 10
        self.x_GSE = np.sort(self.x_GSOatGSE)
        self.y_GSE = np.sort(self.y_GSOatGSE)
        self.z_GSE = np.sort(self.z_GSOatGSE)
        [self.X_GSE, self.Y_GSE, self.Z_GSE] = np.meshgrid(self.x_GSE,self.y_GSE,self.z_GSE)

    def Production_GSO(self):  
        '''
        Load NetCDF cube that is saved in the ./Data folder
        '''
        load_path = f"./Data/{self.name}.nc"
        cube_data = nc.Dataset(load_path,'r')
        Prod_GSO = cube_data.variables['Prod'][:]  # order [Z,Y,X] 1/cm^3/s
        cube_data.close()
        # Trasnpose Prog_GSO to get [X,Y,Z] format
        self.Prod_GSO = np.transpose(Prod_GSO, (1, 2, 0)) 

    def Production_GSE(self):
        # Build interpolation function to extract the values of Prod_GSO in the GSE grid
        f_interp = RegularGridInterpolator(
                (self.x_GSOatGSE, self.y_GSOatGSE, self.z_GSOatGSE),
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
        self.Production_GSO()
        self.Production_GSE()
        self.Q_GSE =  self.Prod_GSE * emission * abundance
        return self.Q_GSE
    
    def get_intengration(self, axis):
        '''
        Sum the emissivity over axis
        Choose axis: 0, 1, 2 ==> x, y, z
        '''
        return np.sum(self.Q_GSE, axis=axis) * self.ds/(4*np.pi)

