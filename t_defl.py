# short file to compute the deflection time for an electron in a plasma at a given temperature and redshift
import numpy as np

m_1 = 9.1093897e-28 # grams (mass of electron)

e = 4.8032e-10 # esu (cm^3/2 g^1/2 s^-1) (electron charge)

k = 1.380658e-16 # Boltzmann constant (erg/K)

mu_0 = 1 # unitless in cgs units (4 * np.pi * 1e-7 in SI units)

T = 1e5 # K (is this the temparature that we want to use (from Pop II.1 MS temperature)? Higher T would give longer t_defl)

Z_1 = -1 # electron charge

Z_2 = 1 # proton charge

z = 20 # redshift

# CHECK THIS VALUE
Psi = 2.70118 # value from original problem: 1.655


#---------
KE = Psi * k * T # kinetic energy of the electron

v_1 = np.sqrt(2 * KE / m_1) # cm/s (note: this is the computation from the exam problem)
print("Electron velocity: %2.3e cm/s" % v_1)
#print("Electron velocity: %2.3e km/s" % (v_1 / 1e5))


#---------
n_0 = 5.17e-3 / (1 + 30)**3  # based on z=30 numbers
#print("H number density now: %2.3e cm^-3" % n_0)

n_2 = n_0 * (1 + z)**3 # cm^-3 (H number density at redshift z=20)
print("H number density at z=20: %2.3e cm^-3" % n_2)

n_e = n_2 # cm^-3 (assuming fully ionized hydrogen)


#---------
# calculating the values for t_defl
# Lambda = ((k**3 * T**3)/(np.pi * n_e ))**0.5 * (3/(2*np.abs(Z_1 * Z_2) * e**3)) # approx bmax / bmin
# print("Lambda value: %3.3e" % Lambda)

Lambda = ((k * T)/(np.pi * n_e))**0.5 * ((m_1 * v_1**2) / 2) * 1/(np.abs(Z_1 * Z_2) * e**3)
# print("Lambda test value: %3.3e" % Lambda_test)

t_defl = (m_1**2 * v_1**3) / (8 * np.pi * Z_1**2 * Z_2**2 * e**4 * n_2 * np.log(Lambda)) # s

t_defl_yr = t_defl*(1/(np.pi * 1e7)) # convert the deflection time to years
print("Deflection time is %3.3e yr" % t_defl_yr)


#---------
# distance over which the region stretchs (assuming sharp cuttoff)
d = v_1 * t_defl
#print("Distance", d, "cm")
d_pc = d / 3.0857e18
print("Distance", d_pc, "pc")

#---------
# get the mean free path of ionizing photons at z=20
mfp = 3.4e17 * (1/n_2)
mfp_pc = mfp / 3.0857e18
print("Mean free path is %3.3e pc" % mfp_pc)

#---------
# determine the total distance over which the ionization front can travel
d_total = d_pc + mfp_pc
print("Total distance over which the ionization front can travel is %3.3e pc" % d_total)

#---------
#CHECK THESE VALUES
J = n_e * e * v_1 # currrent density (cm^-1/2 g^1/2 s^-2)
#get the magnitude of the magnetic field generated at z=20
curl_B = mu_0 * J # Gauss/cm
print("Curl of B at z=20: %3.3e G/cm" % curl_B)

B = (4 * np.pi * n_2 * m_1 * v_1**2)**0.5 # Gauss
print("Magnetic field magnitude at z=20: %3.3e G" % B)
