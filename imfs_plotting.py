from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt

def read_and_plot_multiple(files_and_labels):
    """
    Reads gas_NO and gas_OH from multiple files and plots them.

    Parameters:
        files_and_labels (list of tuples): Each tuple contains the file path and the label for the plot.
    """
    plt.figure(figsize=(8, 6))
    plt.rc('font', family='serif')

    for file_path, label in files_and_labels:
        logNO = []
        logOH = []

        # Open the file and read line by line
        with open(file_path, 'r') as file:
            lines = file.readlines()
            for i, line in enumerate(lines):
                # Check for the relevant sections
                if line.startswith("Stellar mass-weighted [N/O]:"):
                    logNO = list(map(float, lines[i + 1].strip().split()))
                elif line.startswith("Stellar mass-weighted [O/H]:"):
                    logOH = list(map(float, lines[i + 1].strip().split()))
    logOH = [value + 12 - 8.69/12 for value in logOH]
    plt.plot(logOH, logNO, label=label)
    plt.xlabel('12+log(O/H)')
    plt.ylabel('log(N/O)')
    plt.savefig('figs/NO_OH_plot.png')
    plt.legend()
    plt.show()

# Example usage
files_and_labels = [
    ("/home/adriana/python/galIMF/simulation_results_from_galaxy_evol/imfKroupaSTF0.5log_SFR0.1SFEN9.0Z_0-15.98/chemical_and_SN_evolution.txt", "$\\alpha=2.35$"),
    #("/home/adriana/python/galIMF/simulation_results_from_galaxy_evol/imfigimfSTF0.5log_SFR0.1SFEN9Z_0-15.98/chemical_and_SN_evolution.txt", "igimf"),
    #("simulation_results_from_galaxy_evol/imfigimfSTF0.5log_SFR0.3SFEN9.0Z_0-15.98/chemical_and_SN_evolution.txt", "Scenario 3"),
]

read_and_plot_multiple(files_and_labels)


