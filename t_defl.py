# file which contains all the standard physical constants for the computation
import numpy as np

m_1 = 9.1093897e-28 # grams

e = 4.8032e-10 # esu (erg^1/2 cm^1/2)

k = 1.380658e-16 # erg/K

T = 1e4 # K 

Z_1 = -1 # electron charge

Z_2 = 1 # proton charge

z = 20 # redshift 


#---------
KE = 1.655 * k * T # kinetic energy of the electron

v_1 = np.sqrt(2 * KE / m_1) # cm/s (note: this is the computation from the exam problem)

print("Electron velocity %2.3e cm/s" % v_1)


#---------
n_2 = (5.17e-3 / (1 + 30)**3) * (1 + z)**3 # cm^-3

print("H number density %2.3e cm^-3" % n_2)

n_e = n_2 # cm^-3 (assuming fully ionized hydrogen)


#---------
# calculating the values for t_defl
Lambda = ((k**3 * T**3)/(np.pi * n_e ))**0.5 * (3/(2*np.abs(Z_1 * Z_2) * e**3))

t_defl = (m_1**2 * v_1**3) / (8 * np.pi * Z_1**2 * Z_2**2 * e**4 * n_2 * np.log(Lambda))

t_defl_yr = t_defl*(1/(np.pi * 1e7))

print("Deflection time is %3.3e yr" % t_defl_yr)
