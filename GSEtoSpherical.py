#Code for written by Samuel Wharton. 

import numpy as np 
from scipy.interpolate import interpn 


class cartesian_to_spherical():
    '''This takes in an emissivity cube in cartesian coordinates and returns it in Shue spherical coordinates. '''
    
    def __init__(self, x_1d, y_1d, z_1d, eta_3d, rmax=20, dr=0.2, dtheta=1, dphi=1):
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
        
        
        
        self.x_1d = x_1d
        self.y_1d = y_1d
        self.z_1d = z_1d
        self.eta_3d = eta_3d
        
        #Define spherical grid. 
        self.define_spherical_shells(dr=dr, dtheta=dtheta, dphi=dphi, rmax=rmax)
        
        #Put original grid into this tuple. 
        self.points_original = (self.x_1d, self.y_1d, self.z_1d)
        
        #Calculate emissivity by interpolation at these new points. 
        new_points = (self.X_sph, self.Y_sph, self.Z_sph)
        
        #Calculate the emissivity on the resampled grid. 
        self.eta_sph = interpn(self.points_original, self.eta_3d, new_points, method='linear', bounds_error=False, fill_value=0).reshape(self.R.shape)
        
        
         
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
 