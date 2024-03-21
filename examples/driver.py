import numpy as np
import matplotlib.pyplot as plt



def main():
    #Set the known values

    #Layer Density
    rho_1 = 1800 # [kg/m^3]
    rho_2 = 2500 # [kg/m^3]

    #Shear wave velocity
    beta_1 = 1900 # [m/s]
    beta_2 = 3200 # [m/s]

    #Layer thickness
    H = 4000 # [m]

    #Frequency values *Not sure if these values actually work*
    f = np.array([-10, -5, -1, 0, 1, 5, 10])
    

    #Equation 1
    g_c = (rho_2/rho_1)*np.sqrt((H ** 2)*((beta_1 ** -2) - (beta_1 ** -2)) - C ** 2) / C

    #Equation 2 (cL -> Love wave velocity [m/s])
    cL = beta_1 ** -2 - (C / H) ** 2

    #Equation 3 (lam_L -> Love wave wavelength [m])
    lam_L = cL / f

