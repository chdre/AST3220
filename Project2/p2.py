import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as scp

N = 10000           #Number of steps

#Arrays
x_arr = np.linspace(1,1e3,N+1)
sigv_vals = np.array([1e-9, 1e-10, 1e-11])   #Sigma*v, [GeV^-2]
W = np.zeros((N+1,3))       #Matrix set up with W values along x, for a given value of sigma*v
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

    W[0,ind] = W_eqf(x_arr[0])  #Initial condition
    W_eq_arr[:,ind] = W_eqf(x_arr)

    WW = scp.solve_ivp(dW_dx, t_span=(1.0,1e3), y0=[W_eqf(x_arr[0])], t_eval=x_arr, method='Radau')
    return WW.y[0]

#Finding the W values.
W[:,0] = solver(sigv_vals[0],0)
W[:,1] = solver(sigv_vals[1],1)
W[:,2] = solver(sigv_vals[2],2)

def xf_finder(ind):
    """
    Finds xf for sigma v value for an array y. Ind tells us which dimension of y
    for sigma*v values we are on for values of y (Not at all a confusing
    sentence).
    """
    y_arr = np.exp(W[:,ind])
    for i in range(N):
        if y_arr[i] <= y_arr[0]*0.9:
            xf = x_arr[i]
            return xf

def abbdm():
    """
    Finding the abbundance of dark matter for each sigma*v and printing the
    values. Calling "abbdm()" will simply print the values.
    """
    abbdm_vals = []
    ind = 0     #counter for which sigv value
    for sigv in sigv_vals:
        abbdm_vals.append(1.69*xf_finder(ind)/20.*(1e-10/sigv))
        ind += 1
    abbdm_vals = np.array(abbdm_vals)
    print abbdm_vals

#[ 0.02111232  0.21112325  2.1112325 ]

"""
Below one can find how the plots were made. A simple call of "plot()" will show
them.
"""
def plot():
    plt.figure()
    plt.grid()
    plt.loglog(x_arr, np.exp(W[:,1]),x_arr, np.exp(W_eq_arr[:,1]))
    plt.xlabel('log(x)'); plt.ylabel('log(y)')
    plt.title('$\\left< \\sigma v \\right> = 1e-9$')
    plt.legend(['$Y$', '$Y_{eq}$'], loc='best')

    plt.figure()
    plt.grid()
    plt.loglog(x_arr, np.exp(W[:,1]),x_arr, np.exp(W_eq_arr[:,1]))
    plt.xlabel('log(x)'); plt.ylabel('log(y)')
    plt.title('$\\left< \\sigma v \\right> = 1e-10$')
    plt.legend(['$Y$', '$Y_{eq}$'], loc='best')

    plt.figure()
    plt.grid()
    plt.loglog(x_arr, np.exp(W[:,2]),x_arr, np.exp(W_eq_arr[:,2]))
    plt.xlabel('log(x)'); plt.ylabel('log(y)')
    plt.title('$\\left< \\sigma v \\right> = 1e-11$')
    plt.legend(['$Y$', '$Y_{eq}$'], loc='best')

    plt.show()
