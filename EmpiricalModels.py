# Code from Nguyen et al. 2022

import numpy as np

'''
For all the spherical functions here, we take the following convetion:

X = r cos(theta)
Y = r sin(theta) sin(phi)
Z = r sin(theta) cos(phi)
'''

def get_cartesian(r, theta, phi):
    x = r*np.cos(theta)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.sin(theta)*np.cos(phi)

    return x, y, z


def get_spherical(x, y, z):
    r = np.sqrt(x**2+y**2+z**2)
    if z != 0:
        phi = np.arctan2(y, z)
    else:
        phi = np.sign(y)*np.pi/2
    theta = np.arccos(x/r)
    return r, theta, phi


def MP_Lin2010(phi_in ,th_in, Pd, Pm, Bz, tilt=0.):
    ''' The Lin 2010 Magnetopause model. Returns the MP distance for a given
    azimuth (phi), zenith (th), solar wind dynamic and magnetic pressures (nPa)
    and Bz (in nT).
    * th: Zenith angle from positive x axis (zenith) (between 0 and pi)
    * phi: Azimuth angle from y axis, about x axis (between 0 abd 2 pi)
    * Pd: Solar wind dynamic pressure in nPa
    * Pm: Solar wind magnetic pressure in nPa
    * tilt: Dipole tilt
    '''
    a = [12.544,
         -0.194,
         0.305,
         0.0573,
         2.178,
         0.0571,
         -0.999,
         16.473,
         0.00152,
         0.382,
         0.0431,
         -0.00763,
         -0.210,
         0.0405,
         -4.430,
         -0.636,
         -2.600,
         0.832,
         -5.328,
         1.103,
         -0.907,
         1.450]

    arr = type(np.array([]))

    if(type(th_in) == arr):
        th = th_in.copy()
    else:
        th = th_in

    if(type(phi_in) == arr):
        phi = phi_in.copy()
    else:
        phi = phi_in

    el = th_in < 0.
    if(type(el) == arr):
        if(el.any()):
            th[el] = -th[el]

            if(type(phi) == type(arr)):
                phi[el] = phi[el]+np.pi
            else:
                phi = phi*np.ones(th.shape)+np.pi*el
    else:
        if(el):
            th = -th
            phi = phi+np.pi

    P = Pd+Pm

    def exp2(i):
        return a[i]*(np.exp(a[i+1]*Bz)-1)/(np.exp(a[i+2]*Bz)+1)

    def quad(i, s):
        return a[i]+s[0]*a[i+1]*tilt+s[1]*a[i+2]*tilt**2

    r0 = a[0]*P**a[1]*(1+exp2(2))

    beta = [a[6] + exp2(7),
            a[10],
            quad(11, [1, 0]),
            a[13]]

    f = np.cos(0.5*th)+a[5]*np.sin(2*th)*(1-np.exp(-th))
    s = beta[0]+beta[1]*np.sin(phi)+beta[2]*np.cos(phi)+beta[3]*np.cos(phi)**2
    f = f**(s)

    c = {}
    d = {}
    TH = {}
    PHI = {}
    e = {}
    for i, s in zip(['n', 's'], [1, -1]):
        c[i] = a[14]*P**a[15]
        d[i] = quad(16, [s, 1])
        TH[i] = quad(19, [s, 0])
        PHI[i] = np.cos(th)*np.cos(TH[i])
        PHI[i] = PHI[i] + np.sin(th)*np.sin(TH[i])*np.cos(phi-(1-s)*0.5*np.pi)
        PHI[i] = np.arccos(PHI[i])
        e[i] = a[21]
    r = f*r0

    Q = c['n']*np.exp(d['n']*PHI['n']**e['n'])
    Q = Q + c['s']*np.exp(d['s']*PHI['s']**e['s'])

    return r+Q


def MP_Shue1998(th, DP, Bz):
    ''' Shue 1998 Magnetopause model. Returns the MP distance for given
    theta (r), dynamic pressure (in nPa) and Bz (in nT).
    * theta: Angle from the x axis (model is cylindrical symmetry)
    * PD: Dynamic Pressure in nPa
    * Bz: z component of IMF in nT'''
    r0 = (10.22+1.29*np.tanh(0.184*(Bz+8.14)))*DP**(-1./6.6)

    a = (0.58-0.007*Bz)*(1+0.024*np.log(DP))
    return r0*(2./(1+np.cos(th)))**a


def MP_Liu2015(phi_in, th_in, Pd, Pm, Bx, By, Bz, tilt=0):

    arr = type(np.array([]))

    if(type(th_in) == arr):
        th = th_in.copy()
    else:
        th = th_in

    if(type(phi_in) == arr):
        phi = phi_in.copy()
    else:
        phi = phi_in

    el = th_in < 0.
    if(type(el) == arr):
        if(el.any()):
            th[el] = -th[el]

            if(type(phi) == type(arr)):
                phi[el] = phi[el]+np.pi
            else:
                phi = phi*np.ones(th.shape)+np.pi*el
    else:
        if(el):
            th = -th
            phi = phi+np.pi

    P = Pd+Pm

    r0 = (10.56+0.956*np.tanh(0.1795*(Bz+10.78)))*P**(-0.1699)

    alpha_0 = (0.4935+0.1095*np.tanh(0.07217*(Bz+6.882)))*(1+0.01182*np.log(Pd))

    alpha_z = 0.06263*np.tanh(0.0251*tilt)

    alpha_phi = (0.06354+0.07764*np.tanh(0.07217*(abs(Bz)+4.851)))*(1+0.01182*np.log(Pd))

    delta_alpha = 0.02582*np.tanh(0.0667*Bx)*np.sign(Bx)

    if Bz !=0:
        omega = np.arctan2(0.1718*By*(By**2+Bz**2)**0.194, Bz)
    else:
        omega = np.sign(By)*np.pi/2

    alpha= alpha_0+alpha_z*np.cos(phi)+(alpha_phi+delta_alpha*np.sign(np.cos(phi)))*np.cos(2*(phi-omega))

    l_n = (0.822+0.2921*np.tanh(0.08792*(Bz+10.12)))*(1-0.01278*tilt)
    l_s = (0.822+0.2921*np.tanh(0.08792*(Bz+10.12)))*(1+0.01278*tilt)
    w = (0.2382+0.005806*np.log(Pd))*(1+0.0002335*tilt**2)

    C = np.exp(-abs(th-l_n)/w)*(1+np.sign(np.cos(phi)))+np.exp(-abs(th-l_s)/w)*(1+np.sign(-np.cos(phi)))

    r = (r0*(2/(1+np.cos(th)))**alpha)*(1-0.1*C*np.cos(phi)**2)
    return r


def Jelinek_boundary(phi, theta, Pdyn, boundary='mp'):
    if boundary =='bs':
        lamb = 1.17
        R = 15.02
        epsilon = 6.55
    else:
        lamb = 1.54
        R = 12.82
        epsilon = 5.26

    R0 = 2*R*Pdyn**(-1/epsilon)
    return R0/(np.cos(theta)+np.sqrt(np.cos(theta)**2+np.sin(theta)*np.sin(theta)*lamb**2))


def MSH_MF_Kobel1994(x, y, z, x0, x1, xc, B0x, B0y, B0z):
    '''
    given a meshgrid (x, y), confocal bow shocks and magnetopause with noses
    x0 and x1 and focii xc, and IMF condition
    '''
    r = np.sqrt((x-xc)**2+y**2+z**2)
    sigma = np.sqrt((x-xc)+r)
    sigma_0 = np.sqrt(2*(x0-xc))
    sigma_1 = np.sqrt(2*(x1-xc))
    C = (sigma_1**2)/(sigma_1**2-sigma_0**2)


    # FIrst, put B0x terms

    Bx = B0x*C*(1-(sigma_0**2)/(2*r))
    By = -B0x*C*sigma_0*sigma_0*y/(2*r*sigma**2)
    Bz = -B0x*C*sigma_0*sigma_0*z/(2*r*sigma**2)

    # Then add B0y terms

    Bx += -C*sigma_0*sigma_0*B0y*y/(r*sigma**2)
    By += C*B0y*(1+sigma_0*sigma_0/sigma**2 - y*y*sigma_0*sigma_0/(r*sigma**4))
    Bz += -C*B0y*y*z*sigma_0*sigma_0/(r*sigma**4)

    # Then add B0z terms

    Bx += -C*sigma_0*sigma_0*B0z*z/(r*sigma**2)
    By += -C*B0z*y*z*sigma_0*sigma_0/(r*sigma**4)
    Bz += C*B0z*(1+sigma_0*sigma_0/sigma**2 - z*z*sigma_0*sigma_0/(r*sigma**4))


    if mesh==True:
        Bx[sigma<sigma_0] = 0
        Bx[sigma>sigma_1] = 0

        By[sigma<sigma_0] = 0
        By[sigma>sigma_1] = 0

        Bz[sigma<sigma_0] = 0
        Bz[sigma>sigma_1] = 0

    return Bx, By, Bz


def MP_Nguyen_2021_a(phi_in, th_in, Pd, Pm, Bz, omega, gamma):
    def inv_cos(t):
        return 2/(1+np.cos(t))
    
    ## return an array of phi and theta in the right ranges of values
    arr = type(np.array([]))

    if(type(th_in) == arr):
        th = th_in.copy()
    else:
        th = th_in

    if(type(phi_in) == arr):
        phi = phi_in.copy()
    else:
        phi = phi_in

    el = th_in < 0.
    if(type(el) == arr):
        if(el.any()):
            th[el] = -th[el]

            if(type(phi) == type(arr)):
                phi[el] = phi[el]+np.pi
            else:
                phi = phi*np.ones(th.shape)+np.pi*el
    else:
        if(el):
            th = -th
            phi = phi+np.pi
            
    #coefficients
    a = [10.73,
         -0.150,
         0.0208,
         0.38 ,
         2.09  ,
         0.55 ,
         0.088,
         0.0150,
         -0.087]
    
    
    P = Pd+Pm
    r0 = a[0]*(1+a[2]*np.tanh(a[3]*Bz+a[4]))*P**(a[1])
    
    alpha0 = a[5]
    alpha1 = a[6]*gamma
    alpha2 = a[7]*np.cos(omega)
    alpha3 = a[8]*np.cos(omega)
    
    alpha = alpha0+alpha1*np.cos(phi)+alpha2*np.sin(phi)**2+alpha3*np.cos(phi)**2
    return r0*inv_cos(th)**alpha


def MP_Liu2015(phi_in, th_in, Pd, Pm, Bx, By, Bz, tilt=0):

    arr = type(np.array([]))

    if(type(th_in) == arr):
        th = th_in.copy()
    else:
        th = th_in

    if(type(phi_in) == arr):
        phi = phi_in.copy()
    else:
        phi = phi_in

    el = th_in < 0.
    if(type(el) == arr):
        if(el.any()):
            th[el] = -th[el]

            if(type(phi) == type(arr)):
                phi[el] = phi[el]+np.pi
            else:
                phi = phi*np.ones(th.shape)+np.pi*el
    else:
        if(el):
            th = -th
            phi = phi+np.pi


    P = Pd+Pm

    r0 = (10.56+0.956*np.tanh(0.1795*(Bz+10.78)))*P**(-0.1699)

    alpha_0 = (0.4935+0.1095*np.tanh(0.07217*(Bz+6.882)))*(1+0.01182*np.log(Pd))

    alpha_z = 0.06263*np.tanh(0.0251*tilt)

    alpha_phi = (0.06354+0.07764*np.tanh(0.07217*(abs(Bz)+4.851)))*(1+0.01182*np.log(Pd))

    delta_alpha = 0.02582*np.tanh(0.0667*Bx)*np.sign(Bx)

    if Bz != 0:
        omega = np.arctan2(0.1718*By*(By**2+Bz**2)**0.194, Bz)
    else:
        omega = np.sign(By)*np.pi/2

    alpha = alpha_0+alpha_z*np.cos(phi)+(alpha_phi+delta_alpha*np.sign(np.cos(phi)))*np.cos(2*(phi-omega))

    l_n = (0.822+0.2921*np.tanh(0.08792*(Bz+10.12)))*(1-0.01278*tilt)
    l_s = (0.822+0.2921*np.tanh(0.08792*(Bz+10.12)))*(1+0.01278*tilt)
    w = (0.2382+0.005806*np.log(Pd))*(1+0.0002335*tilt**2)

    C = np.exp(-abs(th-l_n)/w)*(1+np.sign(np.cos(phi)))+np.exp(-abs(th-l_s)/w)*(1+np.sign(-np.cos(phi)))

    r = (r0*(2/(1+np.cos(th)))**alpha)*(1-0.1*C*np.cos(phi)**2)
    return r

def Q_MHD(a, phi_in, th_in, P, Bz, tilt=0):
   # Returns the indentation term following the expression established in Liu et al. (2015) and adapted in Nguyen et al. (2021)
    arr = type(np.array([]))

    if(type(th_in) == arr):
        th = th_in.copy()
    else:
        th = th_in

    if(type(phi_in) == arr):
        phi = phi_in.copy()
    else:
        phi = phi_in

    el = th_in < 0.
    if(type(el) == arr):
        if(el.any()):
            th[el] = -th[el]

            if(type(phi) == type(arr)):
                phi[el] = phi[el]+np.pi
            else:
                phi = phi*np.ones(th.shape)+np.pi*el
    else:
        if(el):
            th = -th
            phi = phi+np.pi




    l_n = (a[1]+a[2]*np.tanh(a[3]*(Bz+a[4])))*(1-a[5]*tilt)
    l_s = (a[1]+a[2]*np.tanh(a[3]*(Bz+a[4])))*(1+a[5]*tilt)
    w = (a[6]+a[7]*np.log(P))*(1+a[8]*tilt**2)


    C = np.exp(-abs(th-l_n)/w)*(1+np.sign(np.cos(phi)))+np.exp(-abs(th-l_s)/w)*(1+np.sign(-np.cos(phi)))

    Q =a[0]*C*np.cos(phi)**2
    return Q


def MP_GMN_2021_b(phi_in, th_in, Pd, Pm, Bz, omega, gamma, Q='non_indented'):
    a = [ 10.85,
         -0.15,
         2.70e-02,
         2.96e-01,
         2.14,
         5.49e-01,
         7.45e-02,
         1.0e-02,
        -7.13e-02,
         1.23e-01,
         8.77e-01,
         3.29e-01,
         2.11e-01,
         1.013e+01,
         4.64e-01,
         3.26e-01,
         8.35e-02,
         7.21e-03]
    
    P = Pd+Pm

    if Q == 'non_indented':
        q_value = 0
    elif Q == 'indented':
        q_value = Q_MHD(a[9:], phi_in, th_in, P, Bz, gamma)
    else:
        raise ValueError('the only possible values of Q are "indented" and "non indented"')
        
        
    def inv_cos(t):
        return 2/(1+np.cos(t))

    ## return an array of phi and theta in the right ranges of values
    arr = type(np.array([]))

    if(type(th_in) == arr):
        th = th_in.copy()
    else:
        th = th_in

    if(type(phi_in) == arr):
        phi = phi_in.copy()
    else:
        phi = phi_in

    el = th_in < 0.
    if(type(el) == arr):
        if(el.any()):
            th[el] = -th[el]

            if(type(phi) == type(arr)):
                phi[el] = phi[el]+np.pi
            else:
                phi = phi*np.ones(th.shape)+np.pi*el
    else:
        if(el):
            th = -th
            phi = phi+np.pi
            
    #coefficients
    r0 = a[0]*(1+a[2]*np.tanh(a[3]*Bz+a[4]))*P**(a[1])
    
    alpha0 = a[5]
    alpha1 = a[6]*gamma
    alpha2 = a[7]*np.cos(omega)
    alpha3 = a[8]*np.cos(omega)
    
    alpha = alpha0+alpha1*np.cos(phi)+alpha2*np.sin(phi)**2+alpha3*np.cos(phi)**2
    return (1-q_value)*r0*inv_cos(th)**alpha
