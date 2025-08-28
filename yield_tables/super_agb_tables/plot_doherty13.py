import matplotlib.pyplot as plt
import re

# File path
file_path = "yield_tables/super_agb_doherty13.txt"

# Initialize lists to store data
initial_masses = []
metallicities = []
yields = []

# Specify the isotope to filter
chosen_isotope = "&H-1"

# Read the file and process lines
with open(file_path, 'r') as file:
    for line in file:
        if line.startswith('H'):  # Process header lines for mass and metallicity
            match_mass = re.search(r'(\d+\.\d+)M', line)
            match_metallicity = re.search(r'Z=(\d+\.\d+)', line)
            if match_mass and match_metallicity:
                initial_masses.append(float(match_mass.group(1)))
                metallicities.append(float(match_metallicity.group(1)))
        elif line.startswith(chosen_isotope):  # Process lines starting with the chosen isotope
            columns = line.split()
            yield_value = float(columns[1].lstrip('&'))  # Remove '&' symbol and convert to float
            yields.append(yield_value)
            print(yield_value)
print(yields)
# Organize data into a dictionary by metallicity
data_by_metallicity = {}
for mass, metallicity, yield_value in zip(initial_masses, metallicities, yields):
    if metallicity not in data_by_metallicity:
        data_by_metallicity[metallicity] = {"masses": [], "yields": []}
    data_by_metallicity[metallicity]["masses"].append(mass)
    data_by_metallicity[metallicity]["yields"].append(yield_value)

# Plot the data
plt.figure(figsize=(12, 6))
for metallicity, data in data_by_metallicity.items():
    plt.plot(data["masses"], data["yields"], label=f"Z={metallicity}", marker='o')

# Add labels, title, and legend
plt.xlabel("Initial Mass (M☉)")
plt.ylabel("Mg-24 Yields")
plt.title("Mg-24 Yields vs Initial Mass for Different Metallicities")
plt.legend()

# Adjust layout and show the plot
plt.tight_layout()
plt.show()