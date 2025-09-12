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

def plot_high_redshift_gal_data():
    cmap = plt.get_cmap('tab20')

#   points = [
#         #list from Ji+25, https://ui.adsabs.harvard.edu/abs/2025arXiv250512505J/abstract
#         #(log(age in Myr, logNO, 'name', colour, age_err, logNO_err)
#         (np.log10(3*1e6), -1.43, 'Mrk 996 (low density)', cmap(1), np.log10(0.5*1e6), 0.14), #James+09
#         (np.log10(3.5*1e6), -0.21, 'LyC', cmap(2), np.log10(0.9*1e6), 0.11), #Pascale+23
#         # (8.00, -0.88, 'ID150880', cmap(4), 0.10, 0.15), # no information in Stiavelli+25
#         (np.log10(400*1e6), 0.20, 'UNCOvER-45924', cmap(5), 0, 0.06), #Labbe+24
#         # (8.26, -0.93, 'ID1665', cmap(6), 0.15, 0.15), # no information in Stiavelli+25
#         # (7.92, -0.85, 'ID1746', cmap(7), 0.13, 0.15), # no information in Stiavelli+25
#         # (7.72, -0.62, 'ID1477', cmap(8), 0.09, 0.11), # no information in Stiavelli+25
#         # (7.75, -0.76, 'ID60001', cmap(9), 0.03, 0.03), # no information in Stiavelli+25
#         (np.log10(100*1e6), -0.86, 'EXCELS-121806 - uncertain age', cmap(10), 0, 0.10), #Arellano-Cordova+25
#         (np.log10(2.7*1e6), -1.10, 'GS_3073 (low density)', cmap(12), np.log10(0.1*1e6), 0.12), #Ji+24b -> Barchiesi+23 (mentions only GS-14)
#         (np.log10(50*1e6), -0.85, 'GS_9422 (tentative)', cmap(13), 0, 0), #Tacchella+24
#         #(7.96, -0.67, 'ID397', cmap(14), 0.10, 0.14), # no information in Stiavelli+25
#         (np.log10(1.8*1e6), -0.39, 'RXJ2248-ID', cmap(15), 0, 0.10), #Topping+24
#         # (7.65, -0.40, 'GLASS_150008', cmap(16), 0.12, 0.08), #no information in Isobe+23
#         (np.log10(1.6*1e6), -0.6, 'A1703-zdk6', cmap(17), np.log10(0.45*1e6), 0.3), #Topping+25
#         (np.log10(10*1e6), -0.44, 'GN-z8-LAE', cmap(18), 0, 0.36), #Navarro-Carrera+24
#         # (8.37, -0.01, 'CEERS_01019 (AGN)', cmap(19), 0.13, 0), #not sure, Isobe+23 says that they were using 10 Myr
#         # (7.38, -0.59, 'GN-z9p4', cmap(0), 0.15, 0.24), #no information in Schaerer+24
#         (np.log10(70*1e6), -0.25, 'GN-z11', cmap(6), np.log10(40*1e6), 0.05), #Tacchella+23
#         (np.log10(12*1e6), -0.93, 'GS-z9-0', cmap(11), np.log10(1.5*1e6), 0.37), #oldest stars, Curti+25
#         # ((6.69+7.69)/2, (-0.08-0.12)/2, 'GHZ9p', cmap(3), 0, 0), #no information in Napolitano+25
#         (np.log10(28*1e6), -0.25, 'GHZ2r', cmap(4), np.log10(12*1e6), 0.05) #Castellano+24
#     ]
    points = [
        # now (age_yr, logNO, 'name', colour, age_err_yr, logNO_err)
        (3e6, -1.43, 'Mrk 996 (low density)', cmap(1), 0.5e6, 0.14),
        (3.5e6, -0.21, 'LyC', cmap(2), 0.9e6, 0.11),
        (400e6, 0.20, 'UNCOvER-45924', cmap(5), 0.0, 0.06),
        (100e6, -0.86, 'EXCELS-121806 - uncertain age', cmap(10), 0.0, 0.10),
        (2.7e6, -1.10, 'GS_3073 (low density)', cmap(12), 0.1e6, 0.12),
        (50e6, -0.85, 'GS_9422 (tentative)', cmap(13), 0.0, 0.0),
        (1.8e6, -0.39, 'RXJ2248-ID', cmap(15), 0.0, 0.10),
        (1.6e6, -0.6, 'A1703-zdk6', cmap(17), 0.45e6, 0.3),
        (10e6, -0.44, 'GN-z8-LAE', cmap(18), 0.0, 0.36),
        (70e6, -0.25, 'GN-z11', cmap(6), 40e6, 0.05),
        (12e6, -0.93, 'GS-z9-0', cmap(11), 1.5e6, 0.37),
        (28e6, -0.25, 'GHZ2r', cmap(4), 12e6, 0.05)
    ]

    handles = []
    for age, y, label, color, age_err, yerr in points:
        x = np.log10(age)
        # compute asymmetric log errors from linear age errors (guard against age_err==0)
        if age_err and age_err > 0 and age > age_err:
            low = max(age - age_err, 1e-8)
            xerr_minus = x - np.log10(low)
            xerr_plus = np.log10(age + age_err) - x
            # matplotlib expects shape (2,) for asymmetric single-point xerr or (2, N) for many points
            xerr_plot = np.array([[xerr_minus], [xerr_plus]])
        else:
            xerr_plot = None

        sc = plt.scatter(x, y, color=color, label=label, s=20, marker='1', zorder=3)
        plt.errorbar(x, y, xerr=xerr_plot, yerr=(None if yerr == 0 else yerr),
                     fmt='none', elinewidth=0.8, ecolor=color, zorder=2)
        handles.append(sc)

    lgd = plt.legend(handles=handles, labels=[p[2] for p in points],
                     bbox_to_anchor=(1.01, 1), loc='upper left', fontsize=8,
                     title="Observational Data")
    return lgd

def plot_NO_over_time(file_paths):
    data = load_data_with_names(file_paths)
    log_time, N_over_O = data['log_time_axis'], data['gas_N_over_O_list']
    N_over_O = recalculate_NO(N_over_O)
    plot_high_redshift_gal_data()
    plt.scatter(log_time, N_over_O, marker='.')
    plt.plot(log_time, N_over_O, lw=0.8)
    plt.xlabel('log time', fontsize=14)
    plt.ylabel('log(N/O)', fontsize=14)
    plt.title('N over O over time', fontsize=12)
    plt.savefig('./figs/galevo/weird_N_over_O_over_time.png', bbox_inches='tight', dpi=300)
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

def plot_masses(file_paths):
    data = load_data_with_names(file_paths)
    log_time_axis = np.log(data['time step list:'])
    Y_list = data['Y_list:']
    X_list = data['X_list:']
    Z_list = data['Z_list:']
    plt.rc('font', family='serif')
    # plt.rc('xtick', labelsize='x-small')
    # plt.rc('ytick', labelsize='x-small')
    # fig = plt.figure(61, figsize=(3, 2.5))
    # fig.add_subplot(1, 1, 1)
    plt.stackplot(log_time_axis, Y_list, X_list, Z_list, labels=["Y", "X", "Z"])
    plt.title('gas-phase H, He, and metal mass fraction', fontsize=10)
    plt.xlim(7, log_time_axis[-1])
    plt.ylim(0, 1)
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel('stacked mass fraction')
    plt.legend(loc='lower left', prop={'size': 7})
    plt.tight_layout()
    plt.show()
    stellar_Y_list = data['stellar_Y_list:']
    stellar_X_list = data['stellar_X_list:']
    stellar_Z_list = data['stellar_Z_list:']
    # if plot_save is True:
    #     plt.savefig('XYZ_gas_phase.pdf', dpi=250)
    #     fig = plt.figure(62, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    plt.stackplot(log_time_axis, stellar_Y_list, stellar_X_list, stellar_Z_list, labels=["Y", "X", "Z"])
    #     if plot_save is not True:
    #         plt.title('stellar H, He, and metal mass fraction', fontsize=10)
    plt.xlim(7, log_time_axis[-1])
    plt.ylim(0, 1)
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel('stacked mass fraction')
    plt.legend(loc='lower left', prop={'size': 7})
    plt.tight_layout()
    plt.show()
    #     if plot_save is True:
    #         plt.savefig('XYZ_star_MW.pdf', dpi=250)
    #     fig = plt.figure(63, figsize=(3, 2.5))
    #     fig.add_subplot(1, 1, 1)
    stellar_Y_list_luminosity_weighted = data['stellar_Y_list_luminosity_weighted:']
    stellar_X_list_luminosity_weighted = data['stellar_X_list_luminosity_weighted:']
    stellar_Z_list_luminosity_weighted = data['stellar_Z_list_luminosity_weighted:']
    plt.stackplot(log_time_axis, stellar_Y_list_luminosity_weighted, stellar_X_list_luminosity_weighted, stellar_Z_list_luminosity_weighted, labels=["Y", "X", "Z"])
    #     if plot_save is not True:
    #         plt.title('stellar luminosity-weighted H, He, and metal mass fraction', fontsize=10)
    plt.xlim(7, log_time_axis[-1])
    plt.ylim(0, 1)
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel('stacked mass fraction')
    plt.legend(loc='lower left', prop={'size': 7})
    plt.tight_layout()
    plt.show()
    #     if plot_save is True:
    #         plt.savefig('XYZ_star_LW.pdf', dpi=250)

def nitrogen_oxygen_mass_fraction_evolution(file_paths):
    data = load_data_with_names(file_paths)
    time_axis = data['time step list:']
    gas_NO = data['Gas [N/O]:']
    stellar_NO = data['Stellar mass-weighted [N/O]:']
    stellar_NO_luminosity_weighted = data['Stellar luminosity-weighted [N/O]:']
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    plt.plot(time_axis, gas_NO, label='gas')
    plt.plot(time_axis, stellar_NO, label='stellar MW')
    plt.plot(time_axis, stellar_NO_luminosity_weighted, label='stellar LW')
    plt.xlabel(r'log$_{10}$(time [yr])')
    plt.ylabel('N/O')
    plt.title('Nitrogen mass fraction evolution', fontsize=10)
    plt.legend(prop={'size': 7})
    plt.tight_layout()
    plt.show()

def nitrogen_hydrogen_mass_fraction_evolution(file_paths):
    data = load_data_with_names(file_paths)
    time_axis = data['time step list:']
    gas_NH = data['Gas [N/H]:']
    stellar_NH = data['Stellar mass-weighted [N/H]:']
    stellar_NH_luminosity_weighted = data['Stellar luminosity-weighted [N/H]:']
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    plt.plot(time_axis, gas_NH, label='gas')
    plt.plot(time_axis, stellar_NH, label='stellar MW')
    plt.plot(time_axis, stellar_NH_luminosity_weighted, label='stellar LW')
    plt.xlabel('time [yr]')
    plt.ylabel('log(N/H)')
    plt.title('Nitrogen mass fraction evolution', fontsize=10)
    plt.legend(prop={'size': 7})
    plt.tight_layout()
    plt.show()

#plot_NO_over_time("./simulation_results_from_galaxy_evol/solution/correct agb/imfKroupaSTF-4.15alpha2.1log_SFR<module 'IMFs.Kroupa_IMF' from '/Users/adriana_work/Desktop/galIMF/IMFs/Kroupa_IMF.py'>SFEN1.3SFE0.0072Z_0100infall0.008/plots/N_over_O_time.txt")
#plot_mass_evolution("./simulation_results_from_galaxy_evol/solution/correct agb/imfKroupaSTF-4.15alpha2.1log_SFR<module 'IMFs.Kroupa_IMF' from '/Users/adriana_work/Desktop/galIMF/IMFs/Kroupa_IMF.py'>SFEN1.3SFE0.0072Z_015infall0.008/plots/mass_evolution.txt")
#plot_sfh("./simulation_results_from_galaxy_evol/solution/correct agb/imfKroupaSTF-4.15alpha2.1log_SFR<module 'IMFs.Kroupa_IMF' from '/Users/adriana_work/Desktop/galIMF/IMFs/Kroupa_IMF.py'>SFEN1.3SFE0.0072Z_015infall0.008/plots/SFH.txt")
plot_NO_over_time("./simulation_results_from_galaxy_evol/solution/correct agb/imfKroupaSTF-4.15alpha2.1log_SFR<module 'IMFs.Kroupa_IMF' from '/Users/adriana_work/Desktop/galIMF/IMFs/Kroupa_IMF.py'>SFEN1.3SFE0.0072Z_0100infall0.008/plots/N_over_O_time.txt")