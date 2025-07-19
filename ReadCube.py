import numpy as np
import netCDF4 as nc
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interpn 


from astropy.io.votable import parse
import glob
# from scipy.interpolate import interp1d
# from scipy.optimize import fsolve
# import sympy as sp
# from sympy.vector import CoordSys3D



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
    
        
    def cartesian_to_spherical(self, rmax=15, dr=0.2, dtheta=1, dphi=1):
        '''
        Parameters
        ----------
        x_1d - 1D array of x positions
        y_1d - 1D array of y positions
        z_1d - 1D array of z positions 
        eta_3d - 3D array of emissivity values
        rmax - Max radius of spherical grid. def = 20 RE. 
        dr - Radial grid spacing. def = 0.2 RE. 
        dtheta - Theta grid spacing. def = 1 degree. 
        dphi - Phi grid spacing. def = 1 degree. 
        
        Returns
        -------
        self.R, self.THETA, self.PHI, self.eta_sph - 3D arrays. 
        
        '''
  
        #Define spherical grid. 
        self.define_spherical_shells(dr=dr, dtheta=dtheta, dphi=dphi, rmax=rmax)
        
        #Put original grid into this tuple. 
        self.points_original = (self.x_GSE, self.y_GSE, self.z_GSE)
        
        #Calculate emissivity by interpolation at these new points. 
        new_points = (self.X_sph, self.Y_sph, self.Z_sph)
        
        #Calculate the emissivity on the resampled grid. 
        self.eta_sph = interpn(self.points_original, self.Q_GSE, new_points, method='linear', bounds_error=False, fill_value=0).reshape(self.R.shape)
        
        return self.R, self.THETA, self.PHI, self.eta_sph
        
        
         
    def define_spherical_shells(self, dr=0.2, dtheta=1, dphi=1, rmax=15):
        '''This defines a default grid.
        
        Parameters
        ----------
        dr - radial separation in grid. def = 0.2 RE.
        dtheta - angular separation in theta direction. def = 1 deg. 
        dphi - angular separation in phi direction. def = 1 deg. 
        rmax - max distance to take grid out to. def = 10 RE. 
        '''
        
        
        #Save grid resolutions. 
        self.dr = dr 
        self.dtheta = dtheta
        self.dphi = dphi 
        self.rmax = rmax 
        
        #Get number in each array. 
        nr = int((rmax-1)/dr + 1)
        nth = int(180/dtheta + 1)
        nph = int(360/dphi + 1)
        
        #Define a spherical grid. 
        self.r = np.linspace(1, rmax, nr)
        
        theta = np.linspace(-np.deg2rad(65), np.deg2rad(65) , nth)
        phi = np.linspace(0, 2*np.pi, nph) 
        
        self.R, self.THETA, self.PHI = np.meshgrid(self.r, theta, phi) 
        
        #Convert this grid to Cartesian. 
        self.X_sph, self.Y_sph, self.Z_sph = self.convert_shue_to_xyz_coords(self.R, self.THETA, self.PHI) 

    @staticmethod
    def convert_xyz_to_shue_coords(x, y, z):
        '''This will convert the x,y,z coordinates to those used in the Shue model 
         of the magnetopause and bowshock. 

        Parameters
        ---------
        x, y, z - now 3D. Must be entered as numpy arrays.  

        Returns
        -------
        r, theta (rad) and phi (rad)
        '''
        
        # r 
        r = (x**2 + y**2 + z**2)**0.5
       
        # theta - only calc. where coordinate singularities won't occur. 
        theta = np.zeros(r.shape)
        i = np.where(r != 0)
        theta[i] =  np.arccos(x[i]/r[i])

        # phi - only calc. where coordinate singularities won't occur. 
        phi = np.zeros(r.shape)
        phi[i] = np.arctan2(z[i], y[i]) 
     
        return r, theta, phi

    def convert_shue_to_xyz_coords(self, r, theta, phi):
        '''This will convert the spherical coordinates to GSE cartesian. 

        Parameters
        ---------
        r, theta (rad) and phi (rad) - now 3D. Must be entered as numpy arrays.  

        Returns
        -------
        x, y, z - 3D arrays. 
        '''
        
        #Convert to Cartesian coordinates. 
        x = r * np.cos(theta)
        y = r * np.sin(theta) * np.cos(phi)
        z = r * np.sin(theta) * np.sin(phi)

        return x, y, z

    @staticmethod
    def cart2pol(x, y):
        rho = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        return(rho, theta)

    @staticmethod
    def pol2cart(rho, theta):
        x = rho * np.cos(theta)
        y = rho * np.sin(theta)
        return(x, y)
