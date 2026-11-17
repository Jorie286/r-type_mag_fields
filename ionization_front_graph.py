# plot the ionization front velocity evolution relative to radius
import matplotlib.pyplot as plt
import numpy as np

# get a range of radii
r = np.logspace(20.5, 21.5, 100) # (cm) get distances out to ~ 100 pc
cm_to_pc = 1/3.0857e18 # conversion factor to pc
ev_to_erg = 1.6e-12 # conversion factor to erg
m_e = 9.1093897e-28 # grams (mass of electron)

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
