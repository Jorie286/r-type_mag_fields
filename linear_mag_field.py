# short file to compute the deflection time for an electron in a plasma at a given temperature and redshift
import numpy as np

m_1 = 9.1093897e-28 # grams (mass of electron)

e = 4.8032e-10 # esu (cm^3/2 g^1/2 s^-1) (electron charge)

k = 1.380658e-16 # Boltzmann constant (erg/K)

mu_0 = 1 # unitless in cgs units (4 * np.pi * 1e-7 in SI units)

T = 1e5 # K (is this the temparature that we want to use (from Pop III.1 MS temperature)? Higher T would give longer t_defl)

Z_1 = -1 # electron charge

Z_2 = 1 # proton charge

z = 20 # redshift

c = 3e10 # cm/s

epsilon_0 = 1 # in cgs units

# CHECK THIS VALUE
Psi = 2.70118 # value from original problem: 1.655


#---------
KE = Psi * k * T # kinetic energy of the electron

v_1 = np.sqrt(2 * KE / m_1) # cm/s (note: this is the computation from the exam problem)
print(r"Electron velocity in r: %2.3e cm/s" % v_1)
#print("Electron velocity: %2.3e km/s" % (v_1 / 1e5))

# theta hat velocity
v_t = 1e6 # cm/s (10 km/s)

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
d = v_1 * t_defl # cm
#print("Distance", d, "cm")
d_pc = d / 3.0857e18 # convert to pc
print("Distance", d_pc, "pc")

#---------
# get the mean free path of ionizing photons at z=20
mfp = 3.4e17 * (1/n_2) # cm
mfp_pc = mfp / 3.0857e18 # convert to pc
print("Mean free path of photons is %3.3e pc" % mfp_pc)

v_th = np.sqrt((k * T) / m_1) # cm s^-1 (thermal velocity of electron)

# get mean free path of electrons in this medium
mfp_e = v_th * t_defl # cm

#---------
# determine the total distance over which the ionization front can travel
d_total = d_pc + mfp_pc
print("Total distance over which the ionization front can travel is %3.3e pc" % d_total)

print("Mean free path to deflection length", mfp_pc/d_pc)

# get the ratio of the mean free path to the deflection length (dilution of current density)
f = mfp_pc / d_pc

#---------
#CHECK THESE VALUES
J_r = n_e * e * v_1 # currrent density in r (cm^-1/2 g^1/2 s^-2)
print("\nCurrent density in r", J_r)

J_t = n_e * e * v_t # current density in theta (cm^-1/2 g^1/2 s^-2)
print("Current density in theta", J_t)



#---------------
#get the magnitude of the magnetic field generated at z=20
# R = 1e11 # cm (arbitrary choice for the size of the region that the star occupies)
# n_o = n_e * (r/R)**2 # get the number density of electrons at the chosen radius

# R = 1e6 * 3.086e18 * (1/(1+z)) # cm (10% of 1 comoving Mpc in cm, rough size of overdensity at z=20)
# r=R+1 # location to measure the field
# R_1 = 1e3* 3.086e18 * (1/(1+z)) # cm (60% of a comoving kpc, size of a smaller overdensity at z=20)
# r_1=R_1+1 # location to measure the field
f_ion = 0.5 # fraction of ionized hydrogen

# B = ((R**2)/r) * n_o * f_ion * e * v_t * np.log(R/r) # magnetic field in the phi direction (Gauss)
# print("Magnetic field magnitude at z=20 when R=r: %3.3e G" % np.abs(B))
# print("Magnetic field magnitude at z=20 with dilution: %3.3e G" % np.abs(B / f))

# curl_B = mu_0 * J # Gauss/cm
# print("Curl of B at z=20: %3.3e G/cm" % curl_B)

# B = (4 * np.pi * n_2 * m_1 * v_1**2)**0.5 # Gauss
# print("Magnetic field magnitude at z=20: %3.3e G" % B)

print("\nFor B field over volume:")

# CHECK THESE NUMBERS

# magnetic field based on biot-savart law using constant current density and a spherical region of radius r
# B = (4 * np.pi * J_r * f_ion * 1.1 * R**3) / (3 * c * r**2) # Gauss
# B_diluted = B / f # Gauss

B_d = (4 * np.pi * J_r * f_ion * 1.1 * d_total**3) / (3 * c * f * (1)**2) # Gauss
B_d_n = (4 * np.pi * J_r * 1.1 * d_total**3) / (3 * c * (1)**2) # Gauss
# print("Magnetic field magnitude at z=20 using comoving Mpc: %3.3e G" % np.abs(B))
# print("Magnetic field magnitude at z=20 with dilution using comoving Mpc: %3.3e G" % np.abs(B_diluted))
print("Magnetic field magnitude at z=20 without dilution using d_total = mfp + d_tdefl: %3.3e G" % np.abs(B_d_n))
print("Magnetic field magnitude at z=20 with dilution using d_total = mfp + d_tdefl: %3.3e G" % np.abs(B_d))

# B_1 = (4 * np.pi * J_r * f_ion * 1.1 * R_1**3) / (3 * c * r_1**2) # Gauss
# B_diluted_1 = B_1 / f # Gauss

# print("Magnetic field magnitude at z=20 using comoving kpc: %3.3e G" % np.abs(B_1))
# print("Magnetic field magnitude at z=20 with dilution using comoving kpc: %3.3e G" % np.abs(B_diluted_1))

# find the magnetic field magnitude that was generated at z=20 and is now at z=0
z0 = 0 # redshift now
# B_z0 = (B_diluted) * ((1 + z)/(1 + z0))**-2 # Gauss
# print("Magnetic field magnitude at z=0 with dilution for comoving Mpc: %3.3e G" % np.abs(B_z0))

# B_z0_1 = (B_diluted_1) * ((1 + z)/(1 + z0))**-2 # Gauss
# print("Magnetic field magnitude at z=0 with dilution for comoving kpc: %3.3e G" % np.abs(B_z0_1))

B_z0_d = (B_d) * ((1 + z)/(1 + z0))**-2 # Gauss
print("Magnetic field magnitude at z=0 with dilution for d_total = mfp + d_tdefl: %3.3e G" % np.abs(B_z0_d))


print("\nFor B field in 'wire':")

# total_current = J_r * 1.1 * 0.5 * (np.pi * R**2) # with 0.1 overdensity in current

# magnetic field based on biot-savart law using constant current density and a spherical region of radius r
# B_magnitude = (2 * total_current) / (c * r) # Gauss
# B_diluted_magnitude = B_magnitude / f # Gauss
# print("Magnetic field magnitude at z=20 using comoving Mpc wire: %3.3e G" % np.abs(B_magnitude))
# print("Magnetic field magnitude at z=20 with dilution using comoving Mpc wire: %3.3e G" % np.abs(B_diluted_magnitude))

# total_current_1 = J_r * 1.1 * 0.5 * (np.pi * R_1**2) # with 0.1 overdensity in current

# B_1_magnitude = (2 * total_current_1) / (c * r_1) # Gauss
# B_diluted_1_magnitude = B_1_magnitude / f # Gauss
# print("Magnetic field magnitude at z=20 using comoving kpc wire: %3.3e G" % np.abs(B_1_magnitude))
# print("Magnetic field magnitude at z=20 with dilution using comoving kpc wire: %3.3e G" % np.abs(B_diluted_1_magnitude))

total_current_d = J_r * 1.1 * 0.5 * (np.pi * d_total**2) # with 0.1 overdensity in current
B_diluted_d_magnitude = (2 * total_current_d) / (c * (1e-1) * f) # Gauss
total_current_d_n = J_r * 1.1 * (np.pi * d_total**2) # with 0.1 overdensity in current
B_diluted_d_magnitude_n = (2 * total_current_d) / (c * (1e-1)) # Gauss
print("Magnetic field magnitude at z=20 without dilution using d_total = mfp + d_tdefl wire: %3.3e G" % np.abs(B_diluted_d_magnitude_n))
print("Magnetic field magnitude at z=20 with dilution using d_total = mfp + d_tdefl wire: %3.3e G" % np.abs(B_diluted_d_magnitude))

# find the magnetic field magnitude that was generated at z=20 and is now at z=0
# z0 = 0 # redshift now
# B_z0_magnitude = (B_diluted_magnitude) * ((1 + z)/(1 + z0))**-2 # Gauss
# print("Magnetic field magnitude at z=0 with dilution for comoving Mpc wire: %3.3e G" % np.abs(B_z0_magnitude))

# B_z0_1_magnitude = (B_diluted_1_magnitude) * ((1 + z)/(1 + z0))**-2 # Gauss
# print("Magnetic field magnitude at z=0 with dilution for comoving kpc wire: %3.3e G" % np.abs(B_z0_1_magnitude))

B_z0_d_magnitude = (B_diluted_d_magnitude) * ((1 + z)/(1 + z0))**-2 # Gauss
print("Magnetic field magnitude at z=0 with dilution for d_total = mfp + d_tdefl wire: %3.3e G" % np.abs(B_z0_d_magnitude))


# -------------
# approximate the ionization front velocity and find the fraction of the field that will be locked into the background

v_i = 1e53/(4 * np.pi * n_2 * (1e6 * 3.8e18 / (1 + z))**2)
print("\nIonization front velocity: %3.3e cm/s" % v_i)

v_ion = J_r / (e * n_2) #cm/s
print("Rough estimate of the ionization front velocity: %3.3e cm/s" % v_ion)

v_ionization = 2000*1e5 #km/s

t_ioniztion = (mfp/v_ion) * np.pi*10**7
print(t_ioniztion) #30,000 yrs

f_transverse = v_t / v_ionization # fraction of the field that will be transverse in front crossing time
print("\nf transverse", f_transverse)

B_lock_vol = B_d * f_transverse # fraction of the field that is locked in
print("Magnitude of field that is twisted for a volume: %3.3e G" % B_lock_vol)

B_lock_wire = B_diluted_d_magnitude * f_transverse # fraction of the field that is locked in
print("Magnitude of field that is twisted for a wire: %3.3e G" % B_lock_wire)

B_z0_wire = (B_lock_wire) * ((1 + z)/(1 + z0))**-2 # Gauss
print("Magnitude of twisted field today without amplification: %3.3e G" % B_z0_wire)

#---------------
#get the small scale eddy circulation time
t_eddy = (0.01 * d_total * 3.0857e18) / v_t # cm
print("Small scale eddy circulation time: %3.3e yr" % (t_eddy * (1/(np.pi * 1e7))))

#---------------
# get the magnitude of the electric field in the IGM at z=20
E = J_r / ((n_e * e**2 * t_defl)/m_1) # Gauss (CHECK THIS EQUATION)
print("\nElectric field magnitude at z=20: %3.3e G" % np.abs(E))

#---------------
# get the timescale for electron restoration in the reionization front using random walk

D_reg = (1/6) * mfp_e * v_th # cm^2 s^-1, diffusion coeffeicient

t_reg = (d**2)/(2 * 3 * D_reg) # s, drifting/random walk timescale for electrons

print("\nRandom walk timescale without magnetic field %3.3e s" % t_reg)

#get average velocity of random walk
v_rand_reg = d/t_reg # cm/s
print("Velocity during random walk %3.3e cm/s" % v_rand_reg)

J_reg = n_e * e * v_rand_reg # currrent density in r for drifting electrons (cm^-1/2 g^1/2 s^-2)
print("Current density in r for drifting electrons", J_reg)

# get the magnitude of the electric field in the IGM at z=20
E_reg = J_reg / ((n_e * e**2 * t_reg)/m_1) # Gauss (CHECK THIS EQUATION)
print("Electric field magnitude at z=20 for drifting electrons: %3.3e G" % np.abs(E_reg))



B_max = 1.197e-09 # Gauss, undiulted wire mag field at z=20

B_min = 4.801e-13 # Gauss, diluted mag field at z=20

r_L = (m_1*v_1)/(e*B_min)# Larmor radius

D_mag = (1/3) * (r_L**2) * (v_th/mfp_e) # cm^2 s^-1, diffusion coeffeicient

t_mag = (d**2)/(2 * 3 * D_mag) # s, drifting/random walk timescale for electrons
print("\nRandom walk timescale with magnetic field %3.3e s" % t_mag)

#get average velocity of random walk
v_rand_mag = d/t_mag # cm/s
print("Velocity during random walk %3.3e cm/s" % v_rand_mag)

J_mag = n_e * e * v_rand_mag # currrent density in r for drifting electrons (cm^-1/2 g^1/2 s^-2)
print("Current density in r for drifting electrons with magnetic field", J_mag)

# get the magnitude of the electric field in the IGM at z=20
E_mag = J_mag / ((n_e * e**2 * t_mag)/m_1) # Gauss (CHECK THIS EQUATION)
print("Electric field magnitude at z=20 for drifting electrons with magnetic field: %3.3e G" % np.abs(E_mag))
