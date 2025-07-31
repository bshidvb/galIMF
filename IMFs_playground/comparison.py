import given_IMF
from diet_Salpeter_IMF import diet_salpeter_imf
from Kroupa_IMF import kroupa_imf, kroupa_imf_unnormalized
from given_IMF import given_imf
from Salpeter_IMF import custom_imf, custom_imf_unnormalized
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt

alpha3 = 2.3

# Generate a range of stellar masses
masses = np.logspace(-2, 2.5, 500)  # Mass range from 0.01 to ~300 solar masses

# Evaluate the IMF for each mass
salpeter_imf = [diet_salpeter_imf(m, None) for m in masses]

# Evaluate the Kroupa IMF
kroupa_imf_values = [kroupa_imf(m) for m in masses]

# Evaluate the given IMF
given_imf_values = [given_imf(m, 0) for m in masses]

# Evaluate the custom Salpeter IMF
custom_imf_values = [custom_imf(m) for m in masses]

# Plot the IMF
plt.rc('font', family='serif')
plt.figure(figsize=(8, 6))
#plt.loglog(masses, salpeter_imf, label="Diet Salpeter IMF", color="blue")
plt.loglog(masses, kroupa_imf_values, label="Kroupa IMF", color="green")
# plt.loglog(masses, given_imf_values, label="Time dependent IMF", color="purple")
# plt.loglog(masses, custom_imf_values, label="Salpeter IMF", color="red")
plt.xlabel("Stellar Mass ($M_{\\odot}$)")
plt.ylabel("IMF ($\\xi$)")
# plt.title("Comparison of various IMFs")
plt.legend()
plt.savefig("IMF_comparison.png")
plt.show()

# Cumulative integration for the total number of stars
def cumulative_number_of_stars(imf_function, masses):
    cumulative = []
    for m in masses:
        result, _ = quad(imf_function, 0.08, 150)
        cumulative.append(result)
    return cumulative

# Cumulative integration for the total mass of stars
def cumulative_mass_of_stars(imf_function, masses):
    cumulative = []
    for m in masses:
        result, _ = quad(lambda x: x * imf_function(x), 0.08, 150)
        cumulative.append(result)
    return cumulative

# Compute cumulative integrals for each IMF
diet_salpeter_cumulative_number = cumulative_number_of_stars(lambda m: diet_salpeter_imf(m, None), masses)
diet_salpeter_cumulative_mass = cumulative_mass_of_stars(lambda m: diet_salpeter_imf(m, None), masses)

kroupa_cumulative_number = cumulative_number_of_stars(kroupa_imf, masses)
kroupa_cumulative_mass = cumulative_mass_of_stars(kroupa_imf, masses)

given_cumulative_number = cumulative_number_of_stars(lambda m: given_imf(m, 0), masses)
given_cumulative_mass = cumulative_mass_of_stars(lambda m: given_imf(m, 0), masses)

custom_cumulative_number = cumulative_number_of_stars(custom_imf, masses)
custom_cumulative_mass = cumulative_mass_of_stars(custom_imf, masses)


# Plot cumulative number of stars
# plt.figure(figsize=(10, 7))
# plt.loglog(masses, diet_salpeter_cumulative_number, label="Diet Salpeter IMF", color="blue")
# plt.loglog(masses, kroupa_cumulative_number, label="Kroupa IMF", color="green")
# plt.loglog(masses, given_cumulative_number, label="Given IMF", color="purple")
# plt.loglog(masses, custom_cumulative_number, label="Custom Salpeter IMF", color="orange")
# plt.xlabel("Stellar Mass (M☉)")
# plt.ylabel("Cumulative Number of Stars")
# plt.title("Cumulative Number of Stars for Various IMFs")
# plt.legend()
# plt.savefig("Cumulative_number_of_stars.png")

# Plot cumulative mass of stars
# plt.figure(figsize=(10, 7))
# plt.loglog(masses, diet_salpeter_cumulative_mass, label="Diet Salpeter IMF", color="blue")
# plt.loglog(masses, kroupa_cumulative_mass, label="Kroupa IMF", color="green")
# plt.loglog(masses, given_cumulative_mass, label="Given IMF", color="purple")
# plt.loglog(masses, custom_cumulative_mass, label="Custom Salpeter IMF", color="orange")
# plt.xlabel("Stellar Mass (M☉)")
# plt.ylabel("Cumulative Mass of Stars")
# plt.title("Cumulative Mass of Stars for Various IMFs")
# plt.legend()
# plt.savefig("Cumulative_mass_of_stars.png")
