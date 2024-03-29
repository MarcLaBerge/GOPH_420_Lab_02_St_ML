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
    frequency = np.array([0.01, 0.2, 0.50, 0.75, 1, 2, 5, 10])

    #Zeta Max
    C_max = np.sqrt(H ** 2 * ((beta_1 ** -2) - (beta_2 ** -2)))

    #Equation 1 root finding equation g(C) = 0
    def fx(f, C):

        fx = ((rho_2 / rho_1) 
         * np.sqrt((C_max ** 2 - C ** 2)) / C
         - np.tan(2 * np.pi * f * C))
        return fx
    
    #Derviative of equation 1, big function
    def dfdx(f, C):
        dfdx = -((2 * np.pi * f) / np.cos(2 * f * np.pi * C) ** 2 
                + (rho_2 / rho_1) * np.sqrt(C_max ** 2 - C ** 2) / C ** 2 
                + (rho_2 / rho_1) / (np.sqrt(C_max ** 2 - C ** 2)) )
          
        return dfdx
    
    modes = [0, 1, 2, 3]

    # Creating empty array to store roots or zeta
    zeta = np.zeros([len(frequency), len(modes)])
    # Creating empty array to store velocity
    velocity = np.zeros([len(frequency), len(modes)])
    # Creating empty array to store wavelength
    wavelength = np.zeros([len(frequency), len(modes)])
    
    # Loop of frequencies to find zeta (C)
    for i, f in enumerate(frequency):

        # Create the figure
        plt.figure()
        # Create line on the x = 0 axis
        plt.plot([0, C_max], [0,0], '-k')

        # Mode count starting at 0
        k = 0
        
        # Function for the given frequency and derivative
        function = lambda zeta: fx(f, zeta)
        derivative = lambda zeta: dfdx(f, zeta)

        # Initialise the minimum boundary for plotting the real roots
        z0 = 1e-4

        # Set initial guesses by finding the asymptotes of tan(2 * pi * f * C)
        # Keep watch of certain places that N-R will diverge
        x0 = ((2 * k) + 1) / (4 * f) - 1e-4

    # While the initial guess is less than the C_max
        while x0 < C_max:
            # The max bound is around the asymptotes
            z_max = x0

            # Create the range that will be plotted
            z_range = np.linspace(z0, z_max)

            # Plotting function to find the root
            plt.plot(z_range, function(z_range), '-g')
            # Plotting asymptotes
            plt.plot([z_max + 1e-4, z_max + 1e-4], [-50, 50], '-.c')
            plt.plot(x0, 0.0, 'co')
            # Calculating root using modified secant method
            roots, iterations, error = root_newton_raphson(x0, function, derivative)
            # Plotting the values of the root
            plt.plot(roots, function(roots), 'go')

            # Storing the values of the modes into the zeta, velocity and wavelength arrays
            if k < len(modes):
                zeta[i, k] = roots
                # Use functions given in the lab to find the values for the other arrays
                    # Equation 2
                velocity[i, k] = 1 / (np.sqrt((beta_1 ** -2) - (zeta[i, k] / H ** 2)))
                    # Equation 3
                wavelength[i , k] = velocity[i, k] / f
            k += 1
            # Re-Initialise the guess and bounds for the next mode
            x0 = ((2 * k) + 1) / (4 * f) - 1e-4

            # z0 has to be big enough to not show the function going back up to inf (1e-4, isn't big enough)
            z0 = z_max + 2e-4

        if not k and function(x0:= C_max - 1e-4) < 0:
            roots, itr, error = root_newton_raphson(x0, function, derivative)
            zeta[i, k] = roots
            velocity[i, k] = 1 / np.sqrt(beta_1 ** -2 - (zeta[i, k] / H ** 2))
            wavelength[i, k] = velocity[i, k] / f
            plt.plot(roots, function(roots), 'go')

        # Plotting the frequency vs zeta for a given frequency
            # Helps us undertsand what is going on with the different frequencies
        z_max = C_max - 1e-4
        z_range = np.linspace(z0, z_max)
        plt.plot(z_range, function(z_range), '-g')
        plt.ylim((-10,10))
        plt.grid()
        plt.legend(['_nolegend_', '_nolegend_','Asymptotes'])
        plt.xlabel('Zeta (C)')
        plt.ylabel("Frequency(Zeta) [Hz]")
        plt.title(f"Frequency of {f} [Hz]")
        plt.savefig(f"figures/frequency_plots/Frequency_{f}.png")
        
        print(f"\nIt took the Newton-Raphson method {itr} iterations to converge with the relative error of each iterations being {error}")
    # Stop editing the last figure
    plt.close("all")
    


# Plot zeta vs frequency for each mode
    plt.plot(frequency, zeta[:,0], label = 'Mode 0')

    # Finding frequencies for where roots exist - mode 1
    mode1 = np.argwhere(zeta[:,1] > 0).flatten() 
    plt.plot(frequency[mode1[0]:], zeta[mode1[0]:,1], label = 'Mode 1')

    # Finding frequencies for where roots exist - mode 2
    mode2 = np.argwhere(zeta[:,2] > 0).flatten() 
    plt.plot(frequency[mode2[0]:], zeta[mode2[0]:,2], label = 'Mode 2')

    # Finding frequencies for where roots exist - mode 3
    mode3 = np.argwhere(zeta[:,3] > 0).flatten() 
    plt.plot(frequency[mode3[0]:], zeta[mode3[0]:,3], label = 'Mode 3')

    plt.title('Zeta as a Function of Frequency')
    plt.legend()
    plt.xlabel('Frenquency [Hz]')
    plt.ylabel('Zeta (C)')
    plt.savefig('figures/Zeta_vs_Frequency.png')
    
    # Stop editing this figure
    plt.close('all')
    

# Plot velocity vs frequency for each mode
    plt.plot(frequency, velocity[:,0], label = 'Mode 0')

    # Finding frequencies for where roots exist - mode 1
    mode1 = np.argwhere(velocity[:,1] > 0).flatten() 
    plt.plot(frequency[mode1[0]:], velocity[mode1[0]:,1], label = 'Mode 1')

    # Finding frequencies for where roots exist - mode 2
    mode2 = np.argwhere(velocity[:,2] > 0).flatten() 
    plt.plot(frequency[mode2[0]:], velocity[mode2[0]:,2], label = 'Mode 2')

    # Finding frequencies for where roots exist - mode 3
    mode3 = np.argwhere(velocity[:,3] > 0).flatten() 
    plt.plot(frequency[mode3[0]:], velocity[mode3[0]:,3], label = 'Mode 3')

    plt.title('Love Wave Velocity as a Function of Frequency')
    plt.legend()
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Love wave Velocity [m/s]')
    plt.savefig('figures/Velocity_vs_Frequency.png')

    # Stop editing this figure
    plt.close("all")


# Plot wavelength vs frequency for each mode
    # Limit the size of the plot since is soooo big, mode 1 is very big at small frequencies
    plt.ylim(-1, 5000)
    plt.plot(frequency, wavelength[:,0], label = 'Mode 0', linewidth = '0.75')

    # Finding frequencies for where roots exist - mode 1
    mode1 = np.argwhere(wavelength[:,1] > 0).flatten() 
    plt.plot(frequency[mode1[0]:], wavelength[mode1[0]:,1], label = 'Mode 1', linewidth = '0.75')

    # Finding frequencies for where roots exist - mode 2
    mode2 = np.argwhere(wavelength[:,2] > 0).flatten() 
    plt.plot(frequency[mode2[0]:], wavelength[mode2[0]:,2], label = 'Mode 2', linewidth = '0.75')

    # Finding frequencies for where roots exist - mode 3
    mode3 = np.argwhere(wavelength[:,3] > 0).flatten() 
    plt.plot(frequency[mode3[0]:], wavelength[mode3[0]:,3], label = 'Mode 3', linewidth = '0.75')

    plt.title('Love Wave Wavelength as a Function of Frequency')
    plt.legend()
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Love Wave Wavelength [m]')
    plt.savefig('figures/Wavelength_vs_Frequency.png')

    # Stop editing this figure
    plt.close("all")
            
if __name__ == "__main__":
    main()

