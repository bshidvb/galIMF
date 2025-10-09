from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt

masses = np.logspace(-2, 2.5, 500)
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

integrated_mass = quad(mass_function, 0.08, 150, limit=50)[0]

def kroupa_imf(mass, alpha3, time=0):  # normalized to a population with mass = 1 Msun
    if mass < 0.08:
        return 0
    elif mass < 0.5:
        return 2*mass**(-1.3)/integrated_mass
    elif mass < 1:
        return mass**(-2.3)/integrated_mass
    elif mass < 150:
        return mass**(-alpha3)/integrated_mass
    else:
        return 0
    
alpha3 = [1.3, 2.3, 3.0]

plt.figure(figsize=(8, 6))
plt.rc('font', family='serif')
cmap = plt.get_cmap('tab10')
colors = [cmap(i) for i in range(len(alpha3))]

# collect handles but don't add labels here
handles = []
for i, alpha_value in enumerate(alpha3):
    h, = plt.loglog(masses, [kroupa_imf(m, alpha_value) for m in masses],
                    color=colors[i], zorder=1)
    handles.append(h)

# plot the red dashed line last so it's on top / visible
hd, = plt.loglog(masses, [kroupa_imf(m, 2.3) for m in masses],
                 color="red", linestyle='--', zorder=10)

# create legend with the order you want (red dashed in the middle)
plt.legend([handles[0], hd, handles[2]],
           [f"Kroupa IMF (α3={alpha3[0]})",
            "Kroupa IMF (α3=2.3)",
            f"Kroupa IMF (α3={alpha3[2]})"])

plt.xlabel("Stellar Mass (M$_\odot$)", fontsize=14)
plt.ylabel("IMF ($\\xi$)", fontsize=14)
plt.title("Kroupa IMF with different high-mass slopes", fontsize=16)
plt.savefig("./figs/IMF_comparison.png", dpi=300)
plt.tight_layout()
plt.show()