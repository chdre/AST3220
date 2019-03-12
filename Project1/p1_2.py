import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate


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
