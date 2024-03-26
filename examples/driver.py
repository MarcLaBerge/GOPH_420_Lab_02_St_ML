import numpy as np
import matplotlib.pyplot as plt
from goph420_lab02.functions import root_newton_raphson


def main():
    """
    Main driver for lab 2, root finding using Newton-Raphson and secant method
    """

    #Set the known values

    #Layer Density
    rho_1 = 1800 # [kg/m^3]
    rho_2 = 2500 # [kg/m^3]

    #Shear wave velocity
    beta_1 = 1900 # [m/s]
    beta_2 = 3200 # [m/s]

    #Layer thickness
    H = 4000 # [m]
    
    #Set some frequencies
    f = np.array([0.01, 0.04, 0.09, 0.16, 0.25, 0.36, 0.49, 0.64, 0.81, 1, 2, 5, 10])

    #Equation 1 root finding equation g(C) = 0
        #Lambda can take allows user to take multiple arguments (f, C) in one function
    f = lambda f, C: (np.tan(2 * np.pi * f * C) 
                      - (rho_2/rho_1) * np.sqrt((H ** 2)*((beta_1 ** -2) - (beta_1 ** -2)) - C ** 2) / C)

    #Derivative of Equation 1
    dfdx = lambda f, C: ((2 * np.pi * f) * (1 / np.cos(2 * np.pi * f * C) ** 2))

if __name__ == "__main__":
    main()

