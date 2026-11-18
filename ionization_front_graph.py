# plot the ionization front velocity evolution relative to radius
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

# get a range of radii
r = np.logspace(20.5, 21.5, 100) # (cm) get distances out to ~ 100 pc

E_phot = [2, 10, 20] # (eV) a range of post-reionization energies for electrons

cm_to_pc = 1/3.0857e18 # conversion factor to pc

ev_to_erg = 1.6e-12 # conversion factor to erg

m_e = 9.1093897e-28 # grams (mass of electron)

M_star = 1e5 # solar masses

T_star = 1e5 # K (temperature of the star)

A_sun = 6.082e22 # cm^2 (surface area of the sun)

L_sun = 3.826e33 # erg/s (luminosity of the sun)

R_sun = 6.96e10 # cm (radius of the sun)

h = 6.626e-27 # erg s (Planck's constant)

c = 2.998e10 # cm/s (speed of light)

k_B = 1.38e-16 # erg/K (Boltzmann's constant)


# -----------

# get the number of photons emitted per second given the energy of the emitted photon and the star's mass
def N_phot_v1(E, M):
    """
    Function to calculate the number of photons emitted per second given the energy of the emitted photon. This assumes a blackbody spectrum.

    Inputs:
        E, the energy of the emitted photon (eV)
        M, the mass of the star (M_sun)

    Returns:
        N_phot, number of photons emitted per second
    """

    # convert energy to erg
    E_erg = E * ev_to_erg

    # calculate the luminosity of the star
    L = L_sun * 32000 * M # erg/s (mass-luminosity relation for high mass stars)

    # calculate the number of photons emitted per second
    N_phot = L / E_erg

    return N_phot

# calculate the number of photons emitted per second for different post-ionization electron energies
for E in E_phot:
    E_photon = E + 13.6 # eV (add the ionization energy of hydrogen)
    N_phot = N_phot_v1(E_photon, M_star)
    print(f"Number of photons emitted per second for E = {E_photon} eV: {N_phot:.2e} photons/s")


# -----------

# integrate over all energies that can ionize hydrogen to get the total number of photons emitted per second
def N_phot_total(M, T, E_min):
    """
    Function to calculate the total number of photons emitted per second given the mass and temperature of the star. This assumes a blackbody spectrum.

    Inputs:
        M, the mass of the star (M_sun)
        T, the temperature of the star (K)
        E_min, the minimum energy of the emitted photon (eV)

    Returns:
        N_phot_total, the total number of photons emitted per second
    """

    # calculate the radius of the star
    R = 10 * R_sun * (M/400)**(5/11) # cm (mass-radius relation for high mass stars from Yoon et al. 2015)

    # calculate the surface area of the star
    A = 4 * np.pi * R**2

    # Dimensionless minimum energy
    x_min = (E_min * ev_to_erg) / (k_B * T)

    # Dimensionless integrand
    def integrand(x):
        return x**2 / np.expm1(x)

    # integrate from minimum frequency to infinity
    integral, error = quad(integrand, x_min, np.inf)

    # Convert back to photons per second
    N_phot_total = ((2 * np.pi * A / c**2) * (k_B * T / h)**3 * integral)

    return N_phot_total, error

N_phot_total_value, error = N_phot_total(M_star, T_star, 13.6)
print(f"Total number of photons emitted per second: {N_phot_total_value:.2e} photons/s (error: {error:.2e})")

# using table 1 of Bromm et al. 2001, find the number of photons emitted per second (claims to be applicable for above 300 M_sun)
print(f"\nNumber of ionizing photons emitted per second for 13.6 eV: {1.6e48 * 1e5:.2e}")
print(f"Number of ionizing photons emitted per second for 24.6 eV: {1.1e48 * 1e5:.2e}")
print(f"Number of ionizing photons emitted per second for 54.4 eV: {3.8e47 * 1e5:.2e}")

    
# -----------

# define the function for the ionization front velocity
def v_ion(R, N_phot, n_e):
    """
    Function to calculate the ionization front velocity based on photon conservation. This assumes all photons ionize electrons, this will need to be adapted if this is not the case.

    Inputs:
        R, the radius reached by the ionization front (cm)
        N_phot, number of ionizing photons released per second (1/s)
        n_e, number density of electrons within the ionization front (1/cm^3)

    Returns:
        v_ion, the ionization front velocity
    """

    return N_phot/(4 * np.pi * (R**2) * n_e)


v_i = np.empty_like(r)
N_phot = 1e53 # photons/s
f = 0.5 # ionization fraction
p_f = 1e-4 # fraction of the atoms that were already ionized before the front moves through
n_H = 1.607e-3 # number density of protons
n_e = f* (1 - p_f) * 1.607e-3 # number density of electrons correcting for ionization fraction

# calculate the ionization front velocity at each radius
for (i, R) in enumerate(r):
    v_i[i] = v_ion(R, N_phot, n_H)


E_list = [2, 10, 20] # eV
# calculate the electron velocity for a given energy
v_e = []
for E in E_list:
    v_e.append(np.sqrt(2 * (E * ev_to_erg) / m_e)) # cm/s

fig, ax = plt.subplots()
ax.plot(r*cm_to_pc, v_i*1e-5, label = "Ionization Front Velocity")
ax.set_xlabel("Radius (pc)")
ax.set_ylabel("Velocity (km/s)")
ax.hlines(v_e[0], min(r*cm_to_pc), max(r*cm_to_pc), "r", linestyle = "--", label = "Electron Energy = 2 eV")
ax.hlines(v_e[1], min(r*cm_to_pc), max(r*cm_to_pc), "g", linestyle = "--", label = "Electron Energy = 10 eV")
ax.hlines(v_e[2], min(r*cm_to_pc), max(r*cm_to_pc), "m", linestyle = "--", label = "Electron Energy = 20 eV")
ax.legend()
plt.show()
