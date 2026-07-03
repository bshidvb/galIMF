from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker

masses = np.logspace(-1.3, 2.5, 500)
alpha3 = 2.3
def kroupa_imf_unnormalized(mass, alpha3):  # there is no time dependence for Kroupa IMF
        if mass < 0.08:
            return 0
        elif mass < 0.5:
            return 2*mass**(-1.3)
        elif mass < 1:
            return mass**(-2.3)
        elif mass < 150:
            return mass**(-alpha3)
        else:
            return 0

def mass_function(mass):
    return kroupa_imf_unnormalized(mass, alpha3) * mass

integrated_mass = quad(mass_function, 0.1, 150, limit=50)[0]

def kroupa_imf_canon(mass):  # normalized to a population with mass = 1 Msun
    if mass < 0.5:
        return 2*mass**(-1.3)/integrated_mass
    elif mass < 1:
        return mass**(-2.3)/integrated_mass
    elif mass < 150:
        return mass**(-2.3)/integrated_mass

plt.figure(figsize=(8, 6))
plt.rc('font', family='serif')
matplotlib.rcParams['xtick.minor.size'] = 0
matplotlib.rcParams['xtick.minor.width'] = 0
matplotlib.rcParams['ytick.minor.size'] = 0
matplotlib.rcParams['ytick.minor.width'] = 0
cmap = plt.get_cmap('tab10')

# collect handles but don't add labels here
plt.loglog(masses, [kroupa_imf_canon(m) for m in masses],
                 color="red", linestyle='--', zorder=10, lw=2, label='Canonical Kroupa IMF')

def kroupa_imf(mass):  # normalized to a population with mass = 1 Msun
    if mass < 0.5:
        return mass**(0)/integrated_mass
    elif mass < 1:
        return 0.7*mass**(-0.5)/integrated_mass
    elif mass < 150:
        return 0.7*mass**(-1.4)/integrated_mass
 
plt.loglog(masses, [kroupa_imf(m) for m in masses],
                    color='green', zorder=1, lw=2, label='Different slope')

# plot the red dashed line last so it's on top / visible

# create legend with the order you want (red dashed in the middle)
# plt.legend([handles[0], hd], # handles[2]],
#            [f"Kroupa IMF (α3={alpha3[0]})",
#             "Kroupa IMF (α3=2.3)"])
            # f"Kroupa IMF (α3={alpha3[2]})"])

plt.xlabel("Stellar Mass ($\\mathrm{M_{\\odot}}$)", fontsize=12)
plt.ylabel("IMF ($\\mathrm{\\xi}$)", fontsize=12)
# plt.savefig("./IMF_Kroupa_3a.pdf", dpi=300)
plt.tight_layout()
plt.show()

def kroupa_imf_2001(mass):  # normalized to a population with mass = 1 Msun
    # if mass < 0.08:
    #     return 25*mass**(-0.3)/integrated_mass
    if mass < 0.5:
        return 2*mass**(-1.3)/integrated_mass
    elif mass < 1:
        return mass**(-2.3)/integrated_mass
    elif mass < 150:
        return mass**(-2.3)/integrated_mass

plt.loglog(masses, [kroupa_imf_2001(m) for m in masses], color="purple", linestyle='-', zorder=10, lw=2, label="Kroupa IMF")
plt.xlabel("Stellar Mass ($\\mathrm{M_{\\odot}}$)", fontsize=12)
plt.ylabel("IMF ($\\mathrm{\\xi}$)", fontsize=12)
plt.legend()
plt.savefig("./IMF_Kroupa_2001_without_bd.png", dpi=300)
plt.tight_layout()
# plt.show()