import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as scp

#Constants
G = 6.67408e-11     #Gravitational constant, [m^3/kgs^2]

N = 100000          #Number of steps

#Arrays
sigv_vals = np.linspace(1e-14,1e-7,N/10+1)  #Sigma*v, [GeV^-2]
W = np.zeros((N+1,N/10+1))
x_arr = np.linspace(1,1e3,N+1)      #Matrix set up with W values along x, for a given value of sigma*v
W_eq_arr = np.zeros_like(W)

def solver(sigv,ind):
    """
    Solving the ODE by built in initial value problem solving from scipy that
    uses Radau (Runge-Kutta 4/5) method. Function takes a value of sigma v and
    an index value that saves the W for correct sigma*v. Returns the
    y-values.
    """
    W_eqf = lambda x_val: np.log(9.35e9*(sigv/1e-10)*x_val**(3/2.)*np.exp(-x_val))
    dW_dx = lambda x_val,W_val: (np.exp(2.*W_eqf(x_val) - W_val) - np.exp(W_val))*1./x_val**2

    W[0,ind] = W_eqf(x_arr[0])

    WW = scp.solve_ivp(dW_dx, t_span=(1.0,1e3), y0=[W_eqf(x_arr[0])], t_eval=x_arr, method='Radau')
    return WW.y[0]

def xf_finder(sigv,ind):    #ind is y value = dimension of sigma v
    """
    Finds xf for sigma v value for an array y. Ind tells us which dimension of y
    for sigma*v values we are on for values of y (Not at all a confusing
    sentence).
    """
    y = np.exp(solver(sigv,ind))
    for i in range(N):
        if y[i] <= y[0]*0.9:
            xf = x_arr[i]
            return xf


def abbdm2():
    """
    Calculating the abbundance of dark matter as previously. Checking whether
    the value is within the error 0.12 \pm 0.05, and appending value as well as
    the index to their respective lists. The index again tells us the sigma*v
    value. Since the value of the abbundance of DM decreases, we break the loop
    once the value is below 0.05. Found the best value
    """
    ind = 0             #indice for sigma v value
    abbdm_inds = []     #List for abbundance of dm indices
    abbdms = []         #List for the values of the abbundance of dm
    for sigv in sigv_vals:
        abbdm = 1.69*xf_finder(sigv,ind)/20.*(1e-10/sigv)
        if abbdm <= 0.17 and abbdm >= 0.07: #Checking if abb. of dm is within \pm 0.05.
            abbdms.append(abbdm)
            abbdm_inds.append(ind)
        elif abbdm < 0.05:
            break
        ind += 1

    print abbdms[5], abbdm_inds[5], sigv_vals[abbdm_inds[5]]
#0.115408394408 18 1.80009982e-10
abbdm2()
