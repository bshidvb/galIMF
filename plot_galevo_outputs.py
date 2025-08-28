from scipy.integrate import quad
import glob
import numpy as np
import matplotlib.pyplot as plt
import os

def load_data_with_names(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    data_dict = {}
    for i in range(len(lines)):
        line = lines[i].strip()
        if line.startswith("#"):
            row_name = line[1:].strip()
            if i + 1 < len(lines):
                data = list(map(float, lines[i + 1].strip().split()))
                data_dict[row_name] = data
    return data_dict

def recalculate_NO(N_over_O): # Recalculate logOH
    N_over_O = [i - 0.86 for i in N_over_O]
    return N_over_O

def recalculate_OH(O_over_H): # Recalculate logOH
    O_over_H = [i + 8.696 for i in O_over_H]
    return O_over_H

def plot_NO_over_time(file_paths):
    data = load_data_with_names(file_paths)
    log_time, N_over_O = data['log_time_axis'], data['gas_N_over_O_list']
    N_over_O = recalculate_NO(N_over_O)
    plt.scatter(log_time[4:], N_over_O[4:], marker='.')
    plt.plot(log_time[4:], N_over_O[4:], lw=0.8)
    plt.xlabel('log time', fontsize=14)
    plt.ylabel('log(N/O)', fontsize=14)
    plt.title('N over O over time', fontsize=12)
    #plt.savefig('./figs/galevo/weird_N_over_O_over_time.png', bbox_inches='tight', dpi=300)
    plt.show()

def plot_OH_over_time(file_paths):
    data = load_data_with_names(file_paths)
    log_time, O_over_H = data['log_time_axis'], data['gas_O_over_H_list']
    O_over_H = recalculate_OH(O_over_H)
    plt.scatter(log_time[4:], O_over_H[4:], marker='.')
    plt.plot(log_time[4:], O_over_H[4:], lw=0.8)
    plt.xlabel('log time', fontsize=14)
    plt.ylabel('log(O/H)', fontsize=14)
    plt.title('O over H over time', fontsize=12)
    #plt.savefig('./figs/galevo/O_over_H_over_time.png', bbox_inches='tight', dpi=300)
    plt.show()

def plot_sfh(file_paths):
    data = load_data_with_names(file_paths)
    age_list, sfr = data['age_list'], data['SFR_list']
    plt.plot(age_list, sfr, lw=0.8)
    plt.xlabel('Time (Gyr)')
    plt.ylabel('log(SFR)')
    plt.title('Star Formation History')
    plt.savefig('./figs/galevo_output_plots/sfh.png', bbox_inches='tight', dpi=300)
    plt.show()

def plot_mass_evolution(file_paths):
    data = load_data_with_names(file_paths)
    time_axis, total_gas_mass_list = data['time_axis'], data['total_gas_mass_list']
    ejected_gas_mass_list = data['ejected_gas_mass_list']
    stellar_mass_list = data['stellar_mass_list']
    remnant_mass_list = data['remnant_mass_list']
    BH_mass_list = data['BH_mass_list']
    NS_mass_list = data['NS_mass_list']
    WD_mass_list = data['WD_mass_list']
    plt.plot(time_axis, total_gas_mass_list, lw=1.5, label='gas', ls='dotted', c='k')
    plt.plot(time_axis, ejected_gas_mass_list, lw=0.8, label='ejected gas')
    plt.plot(time_axis, stellar_mass_list, lw=0.8, label='living stars', c='k')
    print('plot stellar_mass final', stellar_mass_list[-1])
    plt.plot(time_axis, remnant_mass_list, lw=0.8, label='stellar remnants', ls='dashed', c='k')
    plt.plot(time_axis, BH_mass_list, lw=0.8, label='black holes')
    plt.plot(time_axis, NS_mass_list, lw=0.8, label='neutron stars')
    plt.plot(time_axis, WD_mass_list, lw=0.8, label='white dwarfs')
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel(r'log$_{10}$(Mass [$M_\odot$])')
    plt.title('Mass evolution', fontsize=10)
    plt.legend(fontsize=6, loc='upper left')
    plt.savefig('./figs/galevo_output_plots/mass_evolution.png', bbox_inches='tight', dpi=300)
    plt.show()

#plot_mass_evolution("./simulation_results_from_galaxy_evol/solution/correct agb/imfKroupaSTF-4.15alpha2.1log_SFR<module 'IMFs.Kroupa_IMF' from '/Users/adriana_work/Desktop/galIMF/IMFs/Kroupa_IMF.py'>SFEN1.3SFE0.0072Z_0100infall0.008/plots/mass_evolution.txt")
plot_sfh("./simulation_results_from_galaxy_evol/solution/correct agb/imfKroupaSTF-4.15alpha2.1log_SFR<module 'IMFs.Kroupa_IMF' from '/Users/adriana_work/Desktop/galIMF/IMFs/Kroupa_IMF.py'>SFEN1.3SFE0.0072Z_0100infall0.008/plots/SFH.txt")
