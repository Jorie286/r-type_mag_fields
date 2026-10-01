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
print(r"Electron velocity in r: %2.3e cm/s" % v_1)
#print("Electron velocity: %2.3e km/s" % (v_1 / 1e5))

# theta hat velocity
v_t = 1e5 # cm/s (10 km/s)

# assume no phi hat velocity?

#---------
n_0 = 5.17e-3 / (1 + 30)**3  # based on z=30 numbers
#print("H number density now: %2.3e cm^-3" % n_0)

n_2 = n_0 * (1 + z)**3 # cm^-3 (H number density at redshift z=20)
print("H number density at z=20: %2.3e cm^-3" % n_2)

n_e = n_2 # cm^-3 (assuming fully ionized hydrogen)


#---------
# calculating the values for t_defl (does not depend on theta hat velocity)
# Lambda = ((k**3 * T**3)/(np.pi * n_e ))**0.5 * (3/(2*np.abs(Z_1 * Z_2) * e**3)) # approx bmax / bmin

Lambda = ((k * T)/(np.pi * n_e))**0.5 * ((m_1 * v_1**2) / 2) * 1/(np.abs(Z_1 * Z_2) * e**3)

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

print("Mean free path to deflection length", mfp_pc/d_pc)

#---------
#CHECK THESE VALUES
J_r = n_e * e * v_1 # currrent density in r (cm^-1/2 g^1/2 s^-2)
print("Current density in r", J_r)

J_t = n_e * e * v_t # current density in theta (cm^-1/2 g^1/2 s^-2)
print("Current density in theta", J_t)

# get the ratio of the mean free path to the deflection length (dilution of current density)
f = mfp_pc / d_pc
print("Mean free path to deflection length", f)

#get the magnitude of the magnetic field generated at z=20
# R = 1e11 # cm (arbitrary choice for the size of the region that the star occupies)
# n_o = n_e * (r/R)**2 # get the number density of electrons at the chosen radius

r = 0.1 * 1e6 * 3.086e18  # cm (10% of 1 comoving Mpc in cm, rough size of overdensity at z=20)
r_1 = 0.65 * 1e3* 3.086e18 # cm (60% of a comoving kpc, size of a smaller overdensity at z=20)
f_ion = 0.5 # fraction of ionized hydrogen

# B = ((R**2)/r) * n_o * f_ion * e * v_t * np.log(R/r) # magnetic field in the phi direction (Gauss)
# print("Magnetic field magnitude at z=20 when R=r: %3.3e G" % np.abs(B))
# print("Magnetic field magnitude at z=20 with dilution: %3.3e G" % np.abs(B / f))

# curl_B = mu_0 * J # Gauss/cm
# print("Curl of B at z=20: %3.3e G/cm" % curl_B)

# B = (4 * np.pi * n_2 * m_1 * v_1**2)**0.5 # Gauss
# print("Magnetic field magnitude at z=20: %3.3e G" % B)

# magnetic field based on biot-savart law using constant current density and a spherical region of radius r
B = (4 * np.pi * J_r * f_ion * r) / (3 * 3e10) # Gauss
B_diluted = B / f # Gauss
print("Magnetic field magnitude at z=20 using comoving Mpc: %3.3e G" % np.abs(B))
print("Magnetic field magnitude at z=20 with dilution using comoving Mpc: %3.3e G" % np.abs(B_diluted))

B_1 = (4 * np.pi * J_r * f_ion * r_1) / (3 * 3e10) # Gauss
B_diluted_1 = B_1 / f # Gauss
print("Magnetic field magnitude at z=20 using comoving kpc: %3.3e G" % np.abs(B_1))
print("Magnetic field magnitude at z=20 with dilution using comoving kpc: %3.3e G" % np.abs(B_diluted_1))

# find the magnetic field magnitude that was generated at z=20 and is now at z=0
z0 = 0 # redshift now
B_z0 = (B_diluted) * ((1 + z)/(1 + z0))**-2 # Gauss
print("Magnetic field magnitude at z=0 with dilution for comoving Mpc: %3.3e G" % np.abs(B_z0))

B_z0_1 = (B_diluted_1) * ((1 + z)/(1 + z0))**-2 # Gauss
print("Magnetic field magnitude at z=0 with dilution for comoving kpc: %3.3e G" % np.abs(B_z0_1))
