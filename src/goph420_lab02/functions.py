
import numpy as np



def root_newton_raphson(x0, f, dfdx):
    """
    
    Peramters
    ---------
    x0: Float
        Initial guess

    f:  Function
        The entered function

    dfdx:   Function
        The first derivative of f

    Returns
    -------
    Final estimate of the root:
        Float

    Number of iterations it takes to converge:
        Int

    One-dimensional vector of the approximate relative error at each iteration:
        Numpy.ndarray
    
    """
    
    # Variables needed for the loop
        # Max iteration
    max_itr = 300
        # Starting iteration
    itr = 1
        # Approximate relative error tolerance
    tol = 1e-7
        # Starting approximate relative error
    eps_a = 2 * tol
        #Error array at each iteration
    rel_error = np.array([])
    
    # While loop, will stop once error is smaller than tolerance and has completed less than the max iterations
    while eps_a > tol and itr < max_itr:
        # Implementing the Newton-Raphson formula
        x1 = x0 - (f(x0) / dfdx(x0))
        # Updating the relative approximate error for each iteration
        eps_a = np.abs((x1 - x0) / x1)
        
        # Updating the approximate relative error vector
        rel_error = np.append(rel_error, eps_a)
        # Update the guess
        x0 = x1
        # Update the iteration
        itr += 1

        # Error message if the relative error isn't small enough by the time we complete the max iterations
        if itr >= max_itr:
            raise RuntimeError(f"After {itr} iterations, the Newton-Raphson method has not converged with the initial guess of {x0}")

    return x0 , itr , rel_error

def root_secant_modified(x0, dx, f):
    '''
    
    Parameters
    ----------
    x0: Float
        Initial guess
    dx: Float
        Step size for derivative estimation
    f:  Function
        The entered function
        

    Returns
    -------
    Final estimate of the root:
        Float

    Number of iterations it takes to converge:
        Int

    One-dimensional vector of the approximate relative error at each iteration:
        Numpy.ndarray
    '''
    # Variables needed for the loop
        # Max iteration
    max_itr = 100
        # Starting iteration
    itr = 1
        # Approximate relative error tolerance
    tol = 1e-5
        # Starting approximate relative error
    eps_a = 2 * tol
        # Error array at each iteration
    rel_error = np.array([])
        # Zeta max
    x0 = 2
    
    # While loop, will stop once error is smaller than tolerance and has completed less than the max iterations
    while eps_a > tol and itr < max_itr:
        x_change = x0 + dx
        x1 = x0 - (f(x0) * dx) / (f(x_change) - f(x0))
        print("Root is", x1)
        # Updating relative approximate error
        eps_a = np.abs((x1 - x0) / x1)
        # Updating the approximate relative error vector
        rel_error = np.append(rel_error, eps_a)
        # Update the guess
        x0 = x1
        # Update the iteration
        itr += 1
        if itr >= max_itr:
            raise RuntimeError(f"Modified Secant Method couldn't converge fast enough with the initial guess{x0}")
        
    return x0 , itr , rel_error
    