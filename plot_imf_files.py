#!/usr/bin/env python3
"""
Script to read multiple IMF files from a plots directory and plot them all together.
Each file contains mass_list, xi_each_time, and xi_Kroupa data.
"""

import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def parse_imf_file(filepath):
    """
    Parse an IMF file and extract mass_list, xi_each_time, and xi_Kroupa.
    
    Returns:
        tuple: (mass_list, xi_each_time, xi_Kroupa) or (None, None, None) if parsing fails
    """
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        mass_list = None
        xi_each_time = None
        xi_Kroupa = None
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            if line == "# mass_list":
                # Next line should be the data
                if i + 1 < len(lines):
                    mass_list = np.array([float(x) for x in lines[i+1].split()])
                    i += 2
                else:
                    i += 1
            elif line == "# xi_each_time":
                if i + 1 < len(lines):
                    xi_each_time = np.array([float(x) for x in lines[i+1].split()])
                    i += 2
                else:
                    i += 1
            elif line == "# xi_Kroupa":
                if i + 1 < len(lines):
                    xi_Kroupa = np.array([float(x) for x in lines[i+1].split()])
                    i += 2
                else:
                    i += 1
            else:
                i += 1
        
        return mass_list, xi_each_time, xi_Kroupa
    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
        return None, None, None


def plot_imf_files(plots_dir):
    """
    Plot all IMF files from a directory.
    
    Args:
        plots_dir: Path to the plots directory
        output_file: Optional output file path for the plot
        max_files: Maximum number of files to plot (for testing)
    """
    # Find all imf_at_time files
    imf_files = sorted(glob.glob(os.path.join(plots_dir, "imf_at_time_*.txt")))
    
    # Create figure
    plt.rc('font', family='serif')
    plt.figure(figsize=(8, 6))
    
    # Color map for different files
    colors = plt.cm.viridis(np.linspace(0, 1, len(imf_files)))
    
    # Plot each file
    for idx, filepath in enumerate(imf_files):
        mass_list, xi_each_time, xi_Kroupa = parse_imf_file(filepath)
        
        if mass_list is None:
            print(f"Skipping {os.path.basename(filepath)} - parsing failed")
            continue
        
        # Extract time from filename (e.g., "imf_at_time_10_Myr.txt" -> "10 Myr")
        filename = os.path.basename(filepath)
        time_str = filename.replace("imf_at_time_", "").replace(".txt", "")
        
        # Plot xi_Kroupa with solid line
        if xi_Kroupa is not None and len(xi_Kroupa) == len(mass_list):
            plt.plot(mass_list, xi_Kroupa, '-', color='blue', alpha=0.7, linewidth=1.5)
        
        # Plot xi_each_time with dashed line
        if xi_each_time is not None and len(xi_each_time) == len(mass_list):
            plt.plot(mass_list, xi_each_time, '--', color='grey', alpha=0.7, linewidth=1.5)
    
    # Configure plot
    plt.xlabel('Mass', fontsize=14)
    plt.ylabel('xi', fontsize=14)
    # plt.set_title('IMF Evolution Over Time', fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(0, 1), loc="upper left", fontsize=10, frameon=False)
    plt.tight_layout()
    plt.show()

plot_imf_files("/Users/adriana_work/Desktop/galIMF/simulation_results_from_galaxy_evol/igimf/20260701/imfigimfSTF-4.15alpha3.0log_SFR<module 'igimf_SFR_-190459_Fe_over_H_25849' from '/Users/adriana_work/Desktop/galIMF/Generated_IGIMFs/igimf_SFR_-190459_Fe_over_H_25849.py'>SFEN1.2SFE0.006Z_015infall0.009/plots/", output_file=None, max_files=15)