import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

if __name__ =='__main__':
    #Setup
    z0 = 2.0e7
    n = 3000000             #Number of steps

    #Constants
    G = 6.67408e-11         #Gravitational constant, [m^3/kgs^2]
    kappa = np.sqrt(8.*np.pi*G)
    alpha = 1.
    C = 3./2                #Curly C

    #Matrices
    N0 = np.log(1./(1. + z0))
    N = np.linspace(N0,0,n+1)    #Since we are making dN out of N we make the N's out of linspace and z's out of this
    z = np.exp(-N)-1.
    x1 = np.zeros(n+1)
    x2 = np.zeros(n+1)
    x3 = np.zeros(n+1)
    lda = np.zeros(n+1)   #Lambda

    #Initial conditions
    inverse_pl = False       #If running the inverse power law
    if inverse_pl == True:
        gamma = (alpha + 1.)/alpha
        x1[0] = 5e-5
        x2[0] = 1e-8
        x3[0] = 0.9999
        lda[0] = 1e9
    else:
        lda[0] = C
        gamma = 1.
        x1[0] = 0.
        x2[0] = 5e-13
        x3[0] = 0.9999

    #Integration loop
    """
    Using Eulers integration method. Since we choose to create the N's out of
    linspace, we can now use this to create our dN.
    """
    for i in range(n):
        dN = N[i+1]-N[i]
        fac = 3. + 3.*x1[i]**2 - 3.*x2[i]**2 + x3[i]**2    #To simplify dx/dN

        dlda_dN = -np.sqrt(6)*lda[i]**2*(gamma - 1)*x1[i]    #dlambda/dN
        dx1_dN = -3.*x1[i] + 0.5*np.sqrt(6)*lda[i]*x2[i]**2 + 0.5*x1[i]*fac
        dx2_dN = -0.5*np.sqrt(6)*lda[i]*x1[i]*x2[i] + 0.5*x2[i]*fac
        dx3_dN = -2.*x3[i] + 0.5*x3[i]*fac

        lda[i+1] = lda[i] + dlda_dN*dN
        x1[i+1] = x1[i] + dx1_dN*dN
        x2[i+1] = x2[i] + dx2_dN*dN
        x3[i+1] = x3[i] + dx3_dN*dN


    #Density parameters
    omega_r = x3**2
    omega_phi = x1**2 + x2**2
    omega_m = 1 - omega_r - omega_phi
    omega_r0 = 8.4e-5
    omega_m0 = 0.3
    omega_p0 = 1 - omega_r0 - omega_m0

    #Equation of state
    w_phi = (x1**2 - x2**2)/(x1**2 + x2**2)

    """
    Following are the functions for calculating the definite integrals.
    """
    def H_int(w_phi_val, N_val):
        """
        Integral in the Hubble rate.
        """
        result = integrate.simps(-3.*(1. + w_phi_val), N_val)
        return result

    def Hubble_rate(z_val, w_phi_val, N_val):
        result = np.sqrt(omega_r0*(1 + z_val)**4 + omega_m0*(1 + z_val)**3 + omega_p0*np.exp(H_int(w_phi_val,N_val)))
        return result

    def H0t0():
        LCDM = False    #Set to true to calculate for LCDM model
        if LCDM == False:
            result = integrate.simps(1./Hubble_rate(z,w_phi,N), N)
            return result
        else:
            array = np.sqrt(Hubble_LCDM(z))
            result = integrate.simps(1./array, N)
            return result

    def dL():
        """
        Calculating the luminosity distance. We set a new z array, as a limited
        version of the original z array, to consider the z = 0.83.
        """
        fac = 1. + z[2892157]
        result = fac*integrate.simps(1./Hubble_rate(z[:2892157],w_phi[:2892157],N[:2892157]), z[:2892157])
        result2 = fac*integrate.simps(1./Hubble_LCDM(z[:2892157]), z[:2892157])
        return result


    """
    Following are the Hubble rate for lambda CDM model and a function to find
    the index for z = 0.83.
    """
    def Hubble_LCDM(z_val):
        return omega_m0*(1+z_val)**3 + (1 - omega_m0)

    def index_z_finder():
        """
        Found value is 0.83 at index 2892157 for 3000000. 4820264 for 5000000.
        """
        for i in range(n):
            if z[i] <= 0.83 and z[i-1] > 0.83:
                print z[i], i, z[i-1], z[i+1]   #So not to limit the index, the i-1 proved to be best
                break

    """
    Following are functions for plotting.
    """
    def plot_omega():
        plt.figure()
        plt.grid('on')
        plt.xscale('log')
        plt.xlabel('$log(z)$'); plt.ylabel('$\\Omega$')

        plt.plot(z,omega_r, label='$\\Omega_r$')
        plt.plot(z,omega_phi, label='$\\Omega_{\\phi}$')
        plt.plot(z,omega_m, label='$\\Omega_m$')

        plt.legend(['$\\Omega_r$', '$\\Omega_{\\phi}$', '$\\Omega_m$'])
        plt.show()

    def plot_eos():
        plt.plot(z,w_phi, label='$w_{\\phi}$')
        plt.xscale('log')
        plt.grid('on')
        plt.xlabel('$log(z)$'); plt.ylabel('$w_{\\phi}$')
        plt.show()

    def plot_Hubble_rate():
        plt.plot(z,Hubble_rate(z,w_phi,N), label='$\\frac{H}{H_0}$')
        plt.xlabel('$z$'); plt.ylabel('$\\frac{H}{H_0}$')
        plt.grid('on')
        plt.show()

    def plot_Hubble_LCDM():
        plt.plot(z,np.sqrt(Hubble_LCDM(z)))
        plt.xlabel('$z$'); plt.ylabel('$\\frac{H_{\\Lambda CDM}}{H_0}$')
        plt.grid('on')
        plt.show()

"""
Below are statements to run various functions to prove same results as in
the report.
"""
plot_omega()
#plot_eos()
#plot_Hubble_rate()
#plot_Hubble_LCDM()

#print H0t0()
#print dL()

"""
Following is the infamous failed script.
"""
if __name__ == '__TheLastJedi__':
    def integral_calc(w_phi_val,z_val):
        """
        Swapping limits on integral to deal with the negative sign from
        dz/dN = -1e^(-N).
        """
        result = integrate.quad(lambda N_val: -3.*(1 + w_phi_val), 0, np.log(1./(1. + z_val)))[0]
        return result

    def H_square(z_val,integral):
        """
        Used in function Hubble_rate below to calculate H^2/H0^2.
        """
        #result = (omega_r0*np.exp(-4*N_val) + omega_m0*np.exp(-3*N_val) + omega_p0*np.exp(func))
        #integral = integrate.quad(lambda N_val: -3.*(1 + w_phi_val), 0, np.log(1./(1. + z_val)))[0]
        result = omega_r0*(1. + z_val)**4 + omega_m0*(1. + z_val)**3 + omega_p0*np.exp(integral)
        return result

    def Hubble_rate():
        """
        Calculating H/H0. Returns array of H(z)/H0 for all z.
        """
        H_rate = np.zeros(n+1)        #Calculated Hubble parameter squared
        for i in range(n+1):
            integral = integral_calc(w_phi[i],z[i])         #Calculates integral in hubble rate
            H_rate[i] = np.sqrt(H_square(z[i], integral))
        return H_rate

    def H0t02():
        """
        Calculating H0t0 by integrating 1/H_rate = H0/H.
        """
        result = np.zeros(n+1)
        for i in range(n):
            integral = integral_calc(w_phi[i],z[4820264])
            H_rate = np.sqrt(H_square(z[i],integral))
            result[i] = integrate.quad(lambda N_val: 1./H_rate, N[0], 0)[0]
        return np.sum(result)


    def dL():
        """
        Calculating the luminosity distance.
        """
        result = np.zeros(n+1)
        for i in range(n+1):
            #result[i] = integrate.quad(lambda : 1./H_rate[i], 0, N[2892159])[0]
            #result[i] = integrate.quad(lambda N_val: 1./H_rate[i], N[4820264], 0)[0]
            #result[i] = integrate.quad(lambda z_val: 1./(omega_r0*(1. + z[4820264])**4 + omega_m0*(1. + z[4820264])**3 + omega_p0*np.exp(integral_calc(w_phi[i],z[i]))), 0,z[4820264])[0]
            integral = integral_calc(w_phi[i],z[4820264])
            H_rate = np.sqrt(H_square(z[i],integral))
            result[i] = integrate.quad(lambda z_val: 1./H_rate, 0, z[4820264])[0]
        return (1 + z[4820264])*np.sum(result)


    """
    Quad og for-loop
    """
    def integral_quad():
        integral = np.zeros(n+1)
        integral[0] = integrate.quad(lambda w_val: 3.*(1. + w_phi[0]), N[0], N[0])[0]
        for i in range(1,n):
            integral[i] = integrate.simps(3.*(1. + w_phi[:i]), N[:i])
        return integral

    def Hubble_rate_quad(z_val):
        integral = integral_quad()      #Matrix of all integrals up to element i
        result = omega_r0*(1 + z_val)**4 + omega_m0*(1 + z_val)**3 + omega_p0*np.exp(np.sum(integral))
        return result

    def H0t0_quad():
        integral = np.zeros(n+1)
        result = np.zeros(n+1)
        for i in range(n):
            integral[i] = integrate.quad(lambda w_val: -3.*(1. + w_phi[i]), 0, N[i])[0]
            result[i] = integrate.quad(lambda z_val: 1./np.sqrt(omega_r0*(1 + z[i])**4 + omega_m0*(1 + z[i])**3 + omega_p0*np.exp(integral[i])), N[0], 0)[0]
        return np.sum(result)


    """
    Simps
    """
    #H_H0_integral = integrate.simps(3.*(1. + w_phi), N)
    H_H0_integral = lambda w_phi_val: 3.*(1. + w_phi_val)
    H_H0_square = lambda z_val, integral: omega_r0*(1 + z_val)**4 + omega_m0*(1 + z_val)**3 + omega_p0*np.exp(integral)


    def H0t0():
        LCDM = False
        if LCDM == False:
            array = np.sqrt(H_H0_square(z,H_H0_integral))
            result = integrate.simps(1./array, N)
            return result
        else:
            array = np.sqrt(Hubble_LCDM(z))
            result = integrate.simps(1./array, N)
            return result

    def dL():
        z_lim = z[:4820264]
        fac = 1 + z[4820264]
        array = np.sqrt(omega_r0*(1 + z_lim)**4 + omega_m0*(1 + z_lim)**3 + omega_p0*np.exp(H_H0_integral))
        result = fac*integrate.simps(1./array, z_lim)
        return fac*np.sum(result)


    def test_int():
        integral = np.zeros(n+1)
        for i in range(n):
            integral[i] = integrate.quad(lambda vals: -3.*(1. + w_phi[i]), 0, N[i])[0]
        return integral
