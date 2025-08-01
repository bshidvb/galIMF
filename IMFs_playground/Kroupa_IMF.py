from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt

masses = np.logspace(-2, 2.5, 500)
alpha3 = 2.5
def kroupa_imf_unnormalized(mass):  # there is no time dependence for Kroupa IMF
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
    return kroupa_imf_unnormalized(mass) * mass

integrated_mass = quad(mass_function, 0.08, 150, limit=50)[0]

def kroupa_imf(mass, time=0):  # normalized to a population with mass = 1 Msun
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