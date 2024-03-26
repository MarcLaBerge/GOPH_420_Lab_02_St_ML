import numpy as np
import matplotlib.pyplot as plt
from goph420_lab02.functions import root_secant_modified


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
    frequency = np.array([0.01, 0.25, 1, 2, 5, 10])

    #Zeta Max
    C_max = np.sqrt(H ** 2 * ((beta_1 ** -2) - (beta_2 ** -2)))

    #Equation 1 root finding equation g(C) = 0
        #Lambda can take allows user to take multiple arguments (f, C) in one function
    f = lambda f, C: (np.tan(2 * np.pi * f * C) 
                      - (rho_2/rho_1) * np.sqrt((H ** 2)*((beta_1 ** -2) - (beta_1 ** -2)) - C ** 2) / C)
    dx = 0.01

    # Creating empty array to store roots
    roots = np.zeros(len(frequency))
    # Creating empty array to store interations
    iterations = np.zeros(len(frequency))

    # Finding the initial guess
    for i in range(len(frequency)):
        # Finding the mode
        k = 0  
        while k < 3:  
            C_asy = (0.25 / frequency[i]) * (2 * k + 1)
            if C_asy > C_max:
                break
            else:
                C_initial = C_asy - 1e-4
                roots[i, k], iterations[i,k], error = root_secant_modified(C_initial, dx, f)
                k += 1


    
if __name__ == "__main__":
    main()

