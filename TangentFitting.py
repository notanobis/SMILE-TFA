from astropy.io.votable import parse
import numpy as np
import glob
from scipy.interpolate import interp1d
import sympy as sp
from sympy.vector import CoordSys3D



Re = 6371  # Earth radius in km
N = CoordSys3D('N')
theta = sp.symbols('theta')
phi = sp.symbols('phi')
a_symb = sp.symbols('a_symb')
r0_symb = sp.symbols('r0_symb')

theta_vals = np.linspace(-np.pi/2, np.pi/2, 500)
phi_vals = np.linspace(0, np.pi, 500)
Theta_grid, Phi_grid = np.meshgrid(theta_vals, phi_vals)

def shue_analytical():
    """
    Returns the analytical expression for the Shue model in polar coordinates 
    """
    return r0_symb * (2 / (1 + sp.cos(theta)))**a_symb

def shue_analytical_GSE():
    shue = shue_analytical()
    x = shue * sp.cos(theta)
    y = shue * sp.sin(theta) * sp.cos(phi)
    z = shue * sp.sin(theta) * sp.sin(phi)
    return x, y, z

shue_num = sp.lambdify((theta, phi, a_symb, r0_symb), shue_analytical(), modules='numpy')
shue_GSE_num = sp.lambdify((theta,phi,a_symb,r0_symb), shue_analytical_GSE(), modules='numpy')

def tangent_equation(model, satellite_pos):
    """
    Return the dot product - equation = 0 which is the condition for the tangent points
    """
    if model == "shue":
        x, y, z = shue_analytical_GSE()
        R_shue = x*N.i + y*N.j + z*N.k
        satellite_pos_vector = satellite_pos[0]*N.i + satellite_pos[1]*N.j + satellite_pos[2]*N.k
        R_theta = R_shue.diff(theta).simplify()
        R_phi = R_shue.diff(phi).simplify()
        n = R_theta.cross(R_phi).simplify()
    else:
        print("Model not yet implemented.")
    return n.dot(R_shue - satellite_pos_vector)

def read_images(folder_name, subfolder_name):
    """
    Read images from XML files in the specified folder and subfolder.
    Returns 
    """
    folder_path = f"./Data/FOV_data/{folder_name}/{subfolder_name}/"
    xml_files = sorted(glob.glob(folder_path + "*.xml"))
    # List to store data from all files
    all_images = []
    # Loop through each XML file and parse it
    for file in xml_files:
        VOTable = parse(file)
        table = VOTable.get_first_table()
        data = table.array
        values = data['values']
        # Normlize images
        normalized = (values - np.min(values)) / (np.max(values) - np.min(values))
        all_images.append(normalized)
    elevation = data['elevation']
    azimuth = VOTable.get_field_by_id_or_name('azimuth').value
    return all_images, azimuth, elevation, xml_files

def read_satellite_pos_att(folder_name, xml_file_name):
    """
    Read satellite position and attitude data from XML files.
    """
    folder_path = f"./Data/FOV_data/{folder_name}/{xml_file_name}.xml"
    VOTable = parse(folder_path)
    table = VOTable.get_first_table()
    data = table.array
    # Extract position and attitude data
    positions = np.vstack((data['X']/Re, data['Y']/Re, data['Z']/Re)).T # Positions in RE
    # Extract satellite base unit vectors in GSE
    vx = np.vstack((data['Att-vx1'], data['Att-vx2'], data['Att-vx3'])).T
    vy = np.vstack((data['Att-vy1'], data['Att-vy2'], data['Att-vy3'])).T
    vz = np.vstack((data['Att-vz1'], data['Att-vz2'], data['Att-vz3'])).T
    return positions, vx, vy, vz

def solve_tangent_GSE(r0, a, theta, phi, tangent_eq_num):
    """
    Solve the tangent equation in GSE coordinates
    """
    dot_product_num = tangent_eq_num(theta, phi, a, r0)
    mask = np.abs(dot_product_num) < 0.1  
    theta_sol = theta[mask]
    phi_sol = phi[mask]
    
    r_sol = shue_num(theta_sol, phi_sol, a, r0) 

    x_sol = r_sol * np.cos(theta_sol)
    y_sol = r_sol * np.sin(theta_sol) * np.cos(phi_sol)
    z_sol = r_sol * np.sin(theta_sol) * np.sin(phi_sol)

    return x_sol, y_sol, z_sol


def GSE_to_Sat(x_GSE, y_GSE, z_GSE, satellite_pos, vx, vy, vz):
    """
    Convert GSE coordinates to satellite coordinates
    """
    R = np.array([vx, vy, vz])
    GSE_vector = np.array([x_GSE, y_GSE, z_GSE]) - satellite_pos
    sat_vector = np.dot(R, GSE_vector.T )  # Perform matrix multiplication
    return sat_vector[0], sat_vector[1], sat_vector[2]

def Sat_to_SXI(x_sat, y_sat, z_sat):
    """
    Convert satellite coordinates to SXI FOV angles.
    """
    azim = np.arctan2(x_sat, z_sat)  # X-Z plane: up/down
    elev = np.arctan2(-y_sat, z_sat)  # Y-Z plane: left/right

    # FOV limits
    azim_mask = (azim >= np.radians(-7.8)) & (azim<= np.radians(7.7))
    elev_mask = (elev >= np.radians(-13.2)) & (elev <= np.radians(13.2))
    # Combine masks: keep values inside either azim or elev FOV
    fov_mask = azim_mask & elev_mask  # <- this means inside both FOV limits
    # Apply mask
    azim = azim[fov_mask]
    elev = elev[fov_mask]

    return azim, elev

def interpolate(az_tangent, el_tangent, Az, El):
    if len(az_tangent) < 5 or len(el_tangent) < 5:
        return None, None

    # interpolate = interp1d(np.degrees(el_tangent),np.degrees(az_tangent), kind='cubic', fill_value="extrapolate")
    az_deg = np.degrees(az_tangent)
    el_deg = np.degrees(el_tangent)

    unique_el_deg, unique_indices = np.unique(el_deg, return_index=True)
    unique_az_deg = az_deg[unique_indices]

    interpolate = interp1d(unique_el_deg, unique_az_deg, kind='cubic', fill_value="extrapolate")

    az_deg = interpolate(np.array(El[:,0]))
    el_deg = np.array(El[:, 0])

    az_max = Az[-1, -1]
    az_min = Az[0, 0]
    el_max = El[-1, -1]
    el_min = El[0, 0]

    az_ind = ((az_deg - az_min) / (az_max - az_min) * (Az.shape[1] - 1))
    el_ind = ((el_deg - el_min) / (el_max - el_min) * (El.shape[0] - 1))
    az_ind = np.clip(az_ind.round().astype(int), 0, Az.shape[1] - 1)
    el_ind = np.clip(el_ind.round().astype(int), 0, El.shape[0] - 1)
    return az_ind, el_ind 

def shue_tangent_SXI(r0, a, Az, El, satellite_pos, vx, vy, vz):
    theta_vals = np.linspace(-np.pi/2, np.pi/2, 500)
    phi_vals = np.linspace(0, np.pi, 500)
    Theta, Phi = np.meshgrid(theta_vals, phi_vals)
    # Solve the tangent equation in GSE
    tangent_eq_num = sp.lambdify((theta, phi, a_symb, r0_symb), tangent_equation("shue", satellite_pos), modules='numpy')
    x_tangent_GSE, y_tangent_GSE, z_tangent_GSE = solve_tangent_GSE(r0, a, Theta, Phi, tangent_eq_num)
    # Convert GSE coordinates to satellite coordinates
    x_tangent_sat , y_tangent_sat, z_tangent_sat = np.zeros((len(x_tangent_GSE))), np.zeros((len(x_tangent_GSE))), np.zeros((len(x_tangent_GSE)))
    for i in range(len(x_tangent_GSE)):
        x_tangent_sat[i], y_tangent_sat[i], z_tangent_sat[i] = GSE_to_Sat(x_tangent_GSE[i], y_tangent_GSE[i], z_tangent_GSE[i], satellite_pos, vx, vy, vz)
    # Convert satellite coordinates to SXI FOV angles
    azim_SXI, elev_SXI = Sat_to_SXI(x_tangent_sat, y_tangent_sat, z_tangent_sat)
    # interpolate the azimuth and elevation angles to the image grid
    az_ind, el_ind = interpolate(azim_SXI, elev_SXI, Az, El)
    return az_ind, el_ind

def hough(image, satellite_pos, vx, vy, vz, Az, El, r0_array, a_array):
    ### Here we set an arbitrary limit to entering the magnetopause. With the complete satellite data we can do that by using the analyser's data to determine in which region we are in and wether we should process the image.
    d_satellite = np.sqrt(satellite_pos[0]**2+satellite_pos[1]**2+satellite_pos[2]**2)
    if d_satellite<7:
        return None,None,None
    hough = np.zeros((len(r0_array), len(a_array)))
    normalized = (image - np.min(image)) / (np.max(image) - np.min(image))
    tangent_eq_num = sp.lambdify((theta, phi, a_symb, r0_symb), tangent_equation("shue", satellite_pos), modules='numpy')

    for i in range(len(r0_array)):
        for j in range(len(a_array)):
            # Remesh inside FOV
            # Theta, Phi  = remesh(r0_array[i], a_array[j], satellite_pos, vz)
            Theta, Phi = Theta_grid, Phi_grid
            # Get tagent curve
            x_tangent_GSE, y_tangent_GSE, z_tangent_GSE = solve_tangent_GSE(r0_array[i], a_array[j], Theta, Phi, tangent_eq_num)
            x_tangent_sat , y_tangent_sat, z_tangent_sat = np.zeros((len(x_tangent_GSE))), np.zeros((len(x_tangent_GSE))), np.zeros((len(x_tangent_GSE)))
            for k in range(len(x_tangent_GSE)):
                x_tangent_sat[k], y_tangent_sat[k], z_tangent_sat[k] = GSE_to_Sat(x_tangent_GSE[k], y_tangent_GSE[k], z_tangent_GSE[k], satellite_pos, vx, vy, vz)
            azim_SXI, elev_SXI = Sat_to_SXI(x_tangent_sat, y_tangent_sat, z_tangent_sat)
            az_ind, el_ind = interpolate(azim_SXI, elev_SXI, Az, El)
            hough[i,j] = np.mean(normalized[el_ind,az_ind])
    max_index = np.unravel_index(np.argmax(hough), hough.shape)
    return hough, r0_array[max_index[0]], a_array[max_index[1]]


def tangent_points(R_numerical, Theta, Phi, satellite_pos):
    x = R_numerical * np.cos(Theta)
    y = R_numerical * np.sin(Theta) * np.cos(Phi)
    z = R_numerical * np.sin(Theta) * np.sin(Phi)

    dx_dtheta = np.gradient(x, Theta[:,0], axis=0)
    dy_dtheta = np.gradient(y, Theta[:,0], axis=0)
    dz_dtheta = np.gradient(z, Theta[:,0], axis=0)

    dx_dphi = np.gradient(x, Phi[0], axis=1)
    dy_dphi = np.gradient(y, Phi[0], axis=1)
    dz_dphi = np.gradient(z, Phi[0], axis=1)

    dtheta = np.stack((dx_dtheta,dy_dtheta,dz_dtheta), axis=-1)
    dphi = np.stack((dx_dphi, dy_dphi, dz_dphi), axis=-1)

    normals = np.cross(dtheta, dphi)
    norms = np.linalg.norm(normals, axis=-1, keepdims=True)
    normals /= norms    

    x_sat = x - satellite_pos[0]
    y_sat = y - satellite_pos[1]
    z_sat = z - satellite_pos[2]
    vec_to_sat = np.stack((x_sat, y_sat, z_sat), axis=-1)  # shape: (n_theta, n_phi, 3)
    norm_vec = np.linalg.norm(vec_to_sat, axis=-1, keepdims=True)
    vec_to_sat_unit = vec_to_sat / norm_vec  # shape: (n_theta, n_phi, 3)

    dot_product_num = np.sum(normals * vec_to_sat_unit, axis=-1)  # shape: (n_theta, n_phi)
    # dot_product_num = normals[..., 0] * x_sat + normals[..., 1] * y_sat + normals[..., 2] * z_sat

    tolerance = np.deg2rad(0.05)  # e.g. 5 degrees grazing
    mask = np.abs(dot_product_num) < np.cos(np.pi/2 - tolerance)
    # mask = np.abs(dot_product_num) < 0.1
    x_sol = x[mask]
    y_sol = y[mask]
    z_sol = z[mask]

    return x_sol, y_sol, z_sol

def surface_tangent_SXI(R_numerical,Theta,Phi, Az, El, satellite_pos, vx, vy, vz):
    # Solve the tangent equation in GSE
    x_tangent_GSE, y_tangent_GSE, z_tangent_GSE = tangent_points(R_numerical, Theta, Phi, satellite_pos)
    # Convert GSE coordinates to satellite coordinates
    x_tangent_sat , y_tangent_sat, z_tangent_sat = np.zeros((len(x_tangent_GSE))), np.zeros((len(x_tangent_GSE))), np.zeros((len(x_tangent_GSE)))
    for i in range(len(x_tangent_GSE)):
        x_tangent_sat[i], y_tangent_sat[i], z_tangent_sat[i] = GSE_to_Sat(x_tangent_GSE[i], y_tangent_GSE[i], z_tangent_GSE[i], satellite_pos, vx, vy, vz)
    # Convert satellite coordinates to SXI FOV angles
    azim_SXI, elev_SXI = Sat_to_SXI(x_tangent_sat, y_tangent_sat, z_tangent_sat)
    # interpolate the azimuth and elevation angles to the image grid
    # az_ind, el_ind = interpolate(azim_SXI, elev_SXI, Az, El)
    return azim_SXI, elev_SXI, z_tangent_sat
# model = "shue"
# folder_name = "01.08.2025-03.08.2025-Res0.5"
# xml_positions = "satellite_pos_att"
# images_folder = "latep"

# positions, vx, vy, vz = read_satellite_pos_att(folder_name, xml_positions)
# images, az_image, el_image, xmls  = read_images(folder_name, images_folder)
# Az, El = np.meshgrid(az_image, el_image)

# i = 10

# az_shue, el_shue = shue_tangent_SXI(7.818, 0.386, Az, El, positions[i], vx[i], vy[i], vz[i])
# plt.pcolormesh(Az, El, images[i], cmap='viridis', shading='auto')
# plt.scatter(Az[el_shue,az_shue], El[el_shue,az_shue], s=2, color = 'black', marker = 's', label='Interpolated values')
# plt.show()
        