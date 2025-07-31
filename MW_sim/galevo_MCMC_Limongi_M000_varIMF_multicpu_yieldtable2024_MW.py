# A python3 code
# This is a single-zone closed-box galaxy chemical evolution module.
# It is coupled with a variable galaxy-wide IMF that depends on the galactic property at the time of star formation.
# The stellar population forms at every 10 Myr (the shortest time step) over 10 Gyr;
# with each stellar population a different galaxy-wide IMF calculated using the IGIMF theory (the galimf.py model).
# The parameters assumed for the simulation are specified at the end of this file or imported from other files,
# i.e., element_weight_table.py, element_abundances_solar.py, element_abundances_primordial.py.

# import ctypes
# from scipy import LowLevelCallable
#
# # Load the shared library
# lib = ctypes.CDLL('IMFs/function_get_target_mass_interpolation.so')  # Update with the correct path
# # Define the function signature for igimf_xi_function
# lib.function_get_target_mass_interpolation.argtypes = [
#     ctypes.c_double, ctypes.POINTER(ctypes.c_double),
#     ctypes.POINTER(ctypes.c_double), ctypes.c_int,
#     ctypes.c_int, ctypes.c_int
# ]
# lib.function_get_target_mass_interpolation.restype = ctypes.c_double
# # def wrapper_function_get_target_mass_interpolation(initial_mass, mass_grid_table2, Mtarget_table, low, high, len_mass_grid_table2):
# #     initial_mass_c = ctypes.c_double(initial_mass)
# #     mass_grid_table2_c = mass_grid_table2.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
# #     Mtarget_table_c = Mtarget_table.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
# #     low_c = ctypes.c_int(low)
# #     high_c = ctypes.c_int(high)
# #     len_mass_grid_table2_c = ctypes.c_int(len_mass_grid_table2)
# #     return lib.function_get_target_mass_interpolation(initial_mass_c, mass_grid_table2_c, Mtarget_table_c, low_c, high_c, len_mass_grid_table2_c)
# function_get_target_mass_interpolation_C_callable = LowLevelCallable(lib.function_get_target_mass_interpolation)

from functools import partial
from multiprocessing import Pool
import glob
import numpy as np
import math
from scipy.integrate import quad
# from scipy.interpolate import UnivariateSpline
from scipy import stats
import sys
import warnings
from astropy.table import Table
import pandas as pd
# import dill
import os
warnings.filterwarnings("ignore")
sys.path.insert(0, 'yield_tables')
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import element_weight_table, element_abundances_solar, element_abundances_primordial
from yield_tables import SNIa_yield

SNIa_yield_table = 'Iwamoto1999_W7'
# SNIa_yield_table = 'Iwamoto1999_W70'

def igimf_xi_function(mass, alpha1, alpha3, integrate_mass_):  # Kroupa 01 IMF, normalized to a population with mass = 1 Msun assuming mass limits 0.08 to 150
    if mass > 50.1:
        return 0
    elif mass > 1:
        return mass ** (-alpha3) / integrate_mass_
    elif mass > 0.08:
        return mass ** (-alpha1-1) / integrate_mass_
    else:
        return 0.5 ** (-alpha1-1) * (mass / 0.5) ** (-alpha1) / integrate_mass_


# output_Fe_over_H_values = np.linspace(-4, 0.6, 100)


def process_epoches(age_of_this_epoch_list, epoch_info_list, alpha1, alpha3, integrate_mass):
    observable_star_number_of_different_epochs_at_this_time_localcpu = []
    ejected_gas_mass_till_this_time_localcpu = 0
    ejected_metal_mass_till_this_time_localcpu = 0
    ejected_H_mass_till_this_time_localcpu = 0
    ejected_O_mass_till_this_time_localcpu = 0
    ejected_N_mass_till_this_time_localcpu = 0
    ejected_C_mass_till_this_time_localcpu = 0
    ejected_Si_mass_till_this_time_localcpu = 0
    ejected_Fe_mass_till_this_time_localcpu = 0
    for i in range(len(age_of_this_epoch_list)):
        age_of_this_epoch = age_of_this_epoch_list[i]
        epoch_info_n = epoch_info_list[i]
        SFR_of_this_epoch = epoch_info_n[0]
        M_tot_of_this_epoch = epoch_info_n[1]
        mass_grid_table = epoch_info_n[4]
        lifetime_table = epoch_info_n[5]
        mass_grid_table2 = epoch_info_n[7]
        M_element_table = epoch_info_n[9]
        epoch_info_n[10] = age_of_this_epoch
        SNIa_number_prob = epoch_info_n[11]
        metal_in_gas = epoch_info_n[12]
        if SFR_of_this_epoch > 0:
            mass_boundary = fucntion_mass_boundary(age_of_this_epoch, mass_grid_table, lifetime_table)
            # mass_boundary_giant_stars = fucntion_mass_boundary(age_of_this_epoch * 1.04, mass_grid_table, lifetime_table)
            observable_star_number = quad(igimf_xi_function, 0.6, max(min(mass_boundary, 1), 0.6),
                                               args=(alpha1, alpha3, integrate_mass), limit=30)[0]  # normalized mass
            observable_star_number_of_a_epoch_at_a_time_step = M_tot_of_this_epoch * observable_star_number  # real number
            # ejected_ :
            len_mass_grid_table2 = len(mass_grid_table2)
            metal_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range(mass_boundary, alpha1, alpha3, integrate_mass,
                                                                                                     len_mass_grid_table2, M_element_table[0])
            H_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range(mass_boundary, alpha1, alpha3, integrate_mass,
                                                                                             len_mass_grid_table2, M_element_table[1])
            He_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range(mass_boundary, alpha1, alpha3, integrate_mass,
                                                                                               len_mass_grid_table2, M_element_table[2])
            ejected_gas_mass_of_this_epoch = H_mass_of_this_epoch + He_mass_of_this_epoch + metal_mass_of_this_epoch
            O_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range(mass_boundary, alpha1, alpha3, integrate_mass,
                                                                                             len_mass_grid_table2, M_element_table[3])
            N_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range(mass_boundary, alpha1, alpha3, integrate_mass,
                                                                                             len_mass_grid_table2, M_element_table[4])
            C_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range(mass_boundary, alpha1, alpha3, integrate_mass,
                                                                                             len_mass_grid_table2, M_element_table[5])
            Si_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range(mass_boundary, alpha1, alpha3, integrate_mass,
                                                                                               len_mass_grid_table2, M_element_table[6])
            Fe_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range(mass_boundary, alpha1, alpha3, integrate_mass,
                                                                                               len_mass_grid_table2, M_element_table[7])
            # if consider SNIa
            Fe_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Fe')
            Si_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Si')
            N_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'N')
            C_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'C')
            O_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'O')
            S_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'S')
            total_mass_eject_per_SNIa = Fe_mass_eject + Si_mass_eject + O_mass_eject + S_mass_eject + N_mass_eject + C_mass_eject
            SNIa_number_from_this_epoch_till_this_time = function_number_SNIa_power_law(0, age_of_this_epoch, SNIa_number_prob, M_tot_of_this_epoch)
            ejected_gas_mass_of_this_epoch += total_mass_eject_per_SNIa * SNIa_number_from_this_epoch_till_this_time
            metal_mass_of_this_epoch += 1.2 * SNIa_number_from_this_epoch_till_this_time  # Chandrasekhar_mass = 1.44
            Fe_mass_of_SNIa = Fe_mass_eject * SNIa_number_from_this_epoch_till_this_time
            O_mass_of_SNIa = O_mass_eject * SNIa_number_from_this_epoch_till_this_time
            Si_mass_of_SNIa = Si_mass_eject * SNIa_number_from_this_epoch_till_this_time
            O_mass_of_this_epoch += O_mass_of_SNIa
            Si_mass_of_this_epoch += Si_mass_of_SNIa
            Fe_mass_of_this_epoch += Fe_mass_of_SNIa
            Fe_over_H_of_an_epoch = function_element_abundunce("Asplund2009", "Fe", "H", metal_in_gas[6], metal_in_gas[1], False)
            observable_star_number_of_different_epochs_at_this_time_localcpu.append([Fe_over_H_of_an_epoch, observable_star_number_of_a_epoch_at_a_time_step])
            ejected_gas_mass_till_this_time_localcpu += ejected_gas_mass_of_this_epoch
            ejected_metal_mass_till_this_time_localcpu += metal_mass_of_this_epoch
            ejected_H_mass_till_this_time_localcpu += H_mass_of_this_epoch
            ejected_O_mass_till_this_time_localcpu += O_mass_of_this_epoch
            ejected_N_mass_till_this_time_localcpu += N_mass_of_this_epoch
            ejected_C_mass_till_this_time_localcpu += C_mass_of_this_epoch
            ejected_Si_mass_till_this_time_localcpu += Si_mass_of_this_epoch
            ejected_Fe_mass_till_this_time_localcpu += Fe_mass_of_this_epoch
        # else:
        #     observable_star_number_of_different_epochs_at_this_time_localcpu.append([0, 0])
    return observable_star_number_of_different_epochs_at_this_time_localcpu, \
           ejected_gas_mass_till_this_time_localcpu, ejected_metal_mass_till_this_time_localcpu, \
           ejected_H_mass_till_this_time_localcpu, ejected_O_mass_till_this_time_localcpu, \
           ejected_N_mass_till_this_time_localcpu, ejected_C_mass_till_this_time_localcpu, \
           ejected_Si_mass_till_this_time_localcpu, ejected_Fe_mass_till_this_time_localcpu


def extract_inner_lists(lst):
    flattened_list = []
    for item in lst:
        if isinstance(item, list):
            flattened_list.extend(item)
    return flattened_list


def galaxy_evol(logOGM=0.5, Z_0=0.000000134, solar_mass_component='Anders1989_mass',
                SF_timescale_limit=400, SFE=0.1, solar_abu_table='Asplund2009',
                outflow=10.0,
                tau_infalle9=0.1, original_gas_mass_fraction=1/100, alpha1=1.3, alpha3=2.3, integrate_mass=4.6, maximum_core_number=14):
    ######################
    # If imf='igimf', the model will use variable IMF, imf='Kroupa' will use Kroupa IMF
    # A 1 in SFH.txt stand for SFR = 1 [solar mass/year] in a 10 Myr epoch.
    # logOGM is the total stellar mass/total gas mass in 13Gyr, which determines the initial gas mass. See Yan et al. 2019
    # Z_0 is the initial metallicity
    ######################
    global mass_grid_table, mass_grid_table2, M_element_table, time_axis, \
        stellar_number_of_different_epochs_list, infall_mass_till_this_time_list, \
        observable_star_number_of_different_epochs_list, original_gas_mass, all_sfr, Z_solar

    ###################
    ### preparation ###
    ###################

    # get all avaliable metallicity from stellar evolution table
    # (Z_table_list, Z_table_list_2) = function_get_avaliable_Z()
    Z_table_list = [0.0004, 0.0008, 0.0012, 0.0016, 0.002, 0.0024, 0.0028, 0.0032, 0.0036, 0.004, 0.008, 0.012]
    Z_table_list_2 = [2e-05, 0.0002, 0.002, 0.02]
    # Z_table_list_2 = [0.008, 0.004, 0.02, 0.0001, 0.001]

    # constrain SFH
    total_gas_mass = 10**logOGM  # in solar mass unit
    original_gas_mass = total_gas_mass * original_gas_mass_fraction
    time_axis = [10**6, 10*10**9, 138*10**8]  ############################################# galaxy age
    time_axis_for_SFH_input = []

    jj = 1
    while jj < SF_timescale_limit:
        time_axis_for_SFH_input += [jj * 10 ** 7]
        (jj) = (jj + 1)

    # the final time axis is the sorted combination of the two
    time_axis = sorted(list(set(time_axis + time_axis_for_SFH_input)))
    length_list_time_step = len(time_axis)

    ###################
    ###  main loop  ###
    ###################
    # define an array save SF event informations that will be used in every latter time steps
    all_sfr = []
    epoch_info = []  # This array saves the properties of stellar populations formed at different time steps
    observable_star_number_of_different_epochs_list = []
    Z_solar = element_abundances_solar.function_solar_element_abundances(solar_mass_component, 'Metal')
    primary_H_mass_fraction = element_abundances_primordial.function_element_mass_primary_fraction(solar_abu_table, "H", Z_0, Z_solar)
    # do calculation for each time start from time 0
    time_step = 0
    gas_infall = True
    infall_mass_till_this_time_list = [1e-9]
    # do calculation for each time to the end time
    while time_step < length_list_time_step:
        # get time
        this_time = time_axis[time_step]
        # calculated the array index (line number in SFH.txt) this_time has reached
        epoch_index_limit = (this_time + 1) / 10 ** 7
        if epoch_index_limit > SF_timescale_limit:
            epoch_index_limit = SF_timescale_limit
        observable_star_number_of_different_epochs_at_this_time = []
        if time_step == 0:
            total_gas_mass_at_this_time = original_gas_mass
            ejected_gas_mass_till_last_time = 0
            ejected_metal_mass_till_last_time = 0
            ejected_H_mass_till_last_time = 0
            ejected_O_mass_till_last_time = 0
            ejected_N_mass_till_last_time = 0
            ejected_C_mass_till_last_time = 0
            ejected_Si_mass_till_last_time = 0
            ejected_Fe_mass_till_last_time = 0
            ejected_gas_mass_till_this_time = 0
            ejected_metal_mass_till_this_time = 0
            ejected_H_mass_till_this_time = 0
            ejected_Fe_mass_till_this_time = 0
            ejected_O_mass_till_this_time = 0
            ejected_N_mass_till_this_time = 0
            ejected_C_mass_till_this_time = 0
            ejected_Si_mass_till_this_time = 0
            Z_gas_this_time_step = Z_0
            total_H_mass_at_last_time = original_gas_mass * primary_H_mass_fraction
            total_metal_mass_in_gas_at_last_time = original_gas_mass * Z_0
            total_gas_mass_at_last_time = original_gas_mass
            ###### element_mass_fraction_after_popiii ######
            total_O_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(solar_abu_table, "O", Z_0, Z_solar)
            total_N_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(solar_abu_table, "N", Z_0, Z_solar)
            total_C_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(solar_abu_table, "C", Z_0, Z_solar)
            total_Si_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(solar_abu_table, "Si", Z_0, Z_solar)
            total_Fe_mass_at_last_time = original_gas_mass * element_abundances_primordial.function_element_mass_primary_fraction(solar_abu_table, "Fe", Z_0, Z_solar)
            metal_mass_in_gas = [Z_0, total_H_mass_at_last_time, total_O_mass_at_last_time, total_N_mass_at_last_time, total_C_mass_at_last_time, total_Si_mass_at_last_time, total_Fe_mass_at_last_time]
        else:
            total_gas_mass_at_last_time = total_gas_mass_at_this_time
            total_metal_mass_in_gas_at_last_time = total_metal_mass_at_this_time
            total_H_mass_at_last_time = total_H_mass_at_this_time
            total_O_mass_at_last_time = total_O_mass_at_this_time
            total_N_mass_at_last_time = total_N_mass_at_this_time
            total_C_mass_at_last_time = total_C_mass_at_this_time
            total_Si_mass_at_last_time = total_Si_mass_at_this_time
            total_Fe_mass_at_last_time = total_Fe_mass_at_this_time
            ejected_gas_mass_till_last_time = ejected_gas_mass_till_this_time
            ejected_metal_mass_till_last_time = ejected_metal_mass_till_this_time
            ejected_H_mass_till_last_time = ejected_H_mass_till_this_time
            ejected_O_mass_till_last_time = ejected_O_mass_till_this_time
            ejected_N_mass_till_last_time = ejected_N_mass_till_this_time
            ejected_C_mass_till_last_time = ejected_C_mass_till_this_time
            ejected_Si_mass_till_last_time = ejected_Si_mass_till_this_time
            ejected_Fe_mass_till_last_time = ejected_Fe_mass_till_this_time
            ejected_gas_mass_till_this_time = 0
            ejected_metal_mass_till_this_time = 0
            ejected_H_mass_till_this_time = 0
            ejected_O_mass_till_this_time = 0
            ejected_N_mass_till_this_time = 0
            ejected_C_mass_till_this_time = 0
            ejected_Si_mass_till_this_time = 0
            ejected_Fe_mass_till_this_time = 0
            Z_gas_this_time_step = total_metal_mass_in_gas_at_last_time / total_gas_mass_at_last_time
            metal_mass_in_gas = [Z_gas_this_time_step, total_H_mass_at_last_time, total_O_mass_at_last_time, total_N_mass_at_last_time, total_C_mass_at_last_time, total_Si_mass_at_last_time, total_Fe_mass_at_last_time]

        age_of_this_epoch = this_time - len(epoch_info) * 10 ** 7
        # SFR_of_this_epoch = (total_gas_mass_at_this_time / 10 ** 6) ** 0.72 * SFE  # SFE = 0.79e-3 de los Reyes 2022
        SFR_of_this_epoch = total_gas_mass_at_this_time / 10 ** 9 * SFE
        # if SFR_of_this_epoch < 10 ** (-9):
        #     SFR_of_this_epoch = 0
        M_tot_of_this_epoch = SFR_of_this_epoch * 10 ** 7
        if SFR_of_this_epoch > 0:
            # Choose the closest metallicity
            Z_select_in_table = function_select_metal(Z_gas_this_time_step, Z_table_list)
            Z_select_in_table_2 = function_select_metal(Z_gas_this_time_step, Z_table_list_2)
            log_Z_select_in_table_2 = [math.log(x) for x in Z_select_in_table_2[1:]]
            # read in interpolated stellar lifetime table
            (mass, lifetime_table) = function_read_lifetime(Z_select_in_table)
            # read in interpolated stellar ejected metal mass
            (mass2, Mmetal_table) = function_read_Mmetal(Z_select_in_table_2, log_Z_select_in_table_2)
            MH_table = function_read_M_element_H(Z_select_in_table_2, log_Z_select_in_table_2)
            MHe_table = function_read_M_element_He(Z_select_in_table_2, log_Z_select_in_table_2)
            MO_table = function_read_M_element_O(Z_select_in_table_2, log_Z_select_in_table_2)
            MN_table = function_read_M_element_N(Z_select_in_table_2, log_Z_select_in_table_2)
            MC_table = function_read_M_element_C(Z_select_in_table_2, log_Z_select_in_table_2)
            MSi_table = function_read_M_element_Si(Z_select_in_table_2, log_Z_select_in_table_2)
            MFe_table = function_read_M_element_Fe(Z_select_in_table_2, log_Z_select_in_table_2)
            M_element_table = [Mmetal_table, MH_table, MHe_table, MO_table, MN_table, MC_table, MSi_table, MFe_table]
            mass_grid_table = mass
            mass_grid_table2 = mass2
            number_in_SNIa_boundary = quad(igimf_xi_function, 2, 8, args=(alpha1, alpha3, integrate_mass), limit=50)[0]
            number_all = quad(igimf_xi_function, 0.08, 50.1, args=(alpha1, alpha3, integrate_mass), limit=50)[0]
            SNIa_number_prob = number_in_SNIa_boundary ** 2 / number_all
            last_time_age = age_of_this_epoch
            all_sfr.append(SFR_of_this_epoch)
            epoch_info.append(
                [SFR_of_this_epoch, M_tot_of_this_epoch, 0, 0,
                 mass_grid_table, lifetime_table, 0, mass_grid_table2, 0, M_element_table,
                 last_time_age, SNIa_number_prob, metal_mass_in_gas, 1])
        else:  # if SFR == 0
            all_sfr.append(1e-11)
            epoch_info.append(
                [0, 0, 0, 0, 0, 0, 0, 0, 0, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 0, 0, [0, 0, 0, 0, 0], 0])

        if epoch_index_limit > max(maximum_core_number, 60):
            iterable_age_of_this_epoch = []
            iterable_epoch_info_n = []
            # processes_num = 4
            processes_num = max(1, min(int(epoch_index_limit-1), maximum_core_number))
            # Calculate the number of tasks per CPU core
            number_of_tasks = int(epoch_index_limit)

            tasks_per_core = math.ceil(number_of_tasks / processes_num)
            start_index = 1
            for i in range(processes_num):
                end_index = min(start_index + tasks_per_core, number_of_tasks + 1)  # Ensure end_index does not exceed N
                subset_indices = list(range(start_index, end_index))
                subset_ages = [this_time - epoch_indexs * 10 ** 7 for epoch_indexs in
                               list(range(start_index, end_index))]
                subset_epoch_info = [epoch_info[index - 1] for index in
                                     subset_indices]  # Adjust indices to 0-based indexing
                iterable_age_of_this_epoch.append(subset_ages)
                iterable_epoch_info_n.append(subset_epoch_info)
                start_index = end_index

            subpool = Pool(processes=processes_num)
            partial_task_function = partial(process_epoches, alpha1=alpha1, alpha3=alpha3, integrate_mass=integrate_mass)
            results = subpool.starmap(partial_task_function, zip(iterable_age_of_this_epoch, iterable_epoch_info_n))
            subpool.close()
            subpool.join()
            # observable_star_number_of_different_epochs_at_this_time = [item for sublist[0] in results for item in sublist[0]]
            observable_star_number_of_different_epochs_at_this_time.extend(extract_inner_lists([result[0] for result in results]))
            ejected_gas_mass_till_this_time += sum(result[1] for result in results)
            ejected_metal_mass_till_this_time += sum(result[2] for result in results)
            ejected_H_mass_till_this_time += sum(result[3] for result in results)
            ejected_O_mass_till_this_time += sum(result[4] for result in results)
            ejected_N_mass_till_this_time += sum(result[5] for result in results)
            ejected_C_mass_till_this_time += sum(result[6] for result in results)
            ejected_Si_mass_till_this_time += sum(result[7] for result in results)
            ejected_Fe_mass_till_this_time += sum(result[8] for result in results)
        else:
            epoch_index = 0
            while epoch_index < epoch_index_limit:
                age_of_this_epoch = this_time - epoch_index * 10 ** 7
                SFR_of_this_epoch = epoch_info[epoch_index][0]
                M_tot_of_this_epoch = epoch_info[epoch_index][1]
                mass_grid_table = epoch_info[epoch_index][4]
                lifetime_table = epoch_info[epoch_index][5]
                mass_grid_table2 = epoch_info[epoch_index][7]
                M_element_table = epoch_info[epoch_index][9]
                epoch_info[epoch_index][10] = age_of_this_epoch
                SNIa_number_prob = epoch_info[epoch_index][11]
                metal_in_gas = epoch_info[epoch_index][12]
                if SFR_of_this_epoch > 0:
                    mass_boundary = fucntion_mass_boundary(age_of_this_epoch, mass_grid_table, lifetime_table)
                    # mass_boundary_giant_stars = fucntion_mass_boundary(age_of_this_epoch * 1.04, mass_grid_table, lifetime_table)
                    # observable_star_number = quad(igimf_xi_function, max(0.08, mass_boundary_giant_stars), mass_boundary, args=(alpha1, alpha3, integrate_mass), limit=30)[0]  # normalized mass
                    # observable_star_number_of_a_epoch_at_a_time_step = M_tot_of_this_epoch * observable_star_number  # real number
                    observable_star_number = quad(igimf_xi_function, 0.6, max(min(mass_boundary, 1), 0.6), args=(alpha1, alpha3, integrate_mass), limit=30)[0]  # normalized mass
                    observable_star_number_of_a_epoch_at_a_time_step = M_tot_of_this_epoch * observable_star_number  # real number
                    # ejected_ :
                    len_mass_grid_table2 = len(mass_grid_table2)
                    metal_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range_metal(mass_boundary, alpha1, alpha3, integrate_mass, len_mass_grid_table2)
                    H_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range_H(mass_boundary, alpha1, alpha3, integrate_mass, len_mass_grid_table2)
                    He_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range_He(mass_boundary, alpha1, alpha3, integrate_mass, len_mass_grid_table2)
                    ejected_gas_mass_of_this_epoch = H_mass_of_this_epoch + He_mass_of_this_epoch + metal_mass_of_this_epoch
                    O_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range_O(mass_boundary, alpha1, alpha3, integrate_mass, len_mass_grid_table2)
                    N_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range_N(mass_boundary, alpha1, alpha3, integrate_mass, len_mass_grid_table2)
                    C_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range_C(mass_boundary, alpha1, alpha3, integrate_mass, len_mass_grid_table2)
                    Si_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range_Si(mass_boundary, alpha1, alpha3, integrate_mass, len_mass_grid_table2)
                    Fe_mass_of_this_epoch = M_tot_of_this_epoch * function_get_target_mass_in_range_Fe(mass_boundary, alpha1, alpha3, integrate_mass, len_mass_grid_table2)
                    # if consider SNIa
                    Fe_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Fe')
                    Si_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'Si')
                    N_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'N')
                    C_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'C')
                    O_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'O')
                    S_mass_eject = SNIa_yield.function_mass_ejected(SNIa_yield_table, 'S')
                    total_mass_eject_per_SNIa = Fe_mass_eject + Si_mass_eject + O_mass_eject + S_mass_eject + N_mass_eject + C_mass_eject
                    # SNIa_energy_release_per_event = 0.8 * 10 ** 51  # in the unit of 10^51 erg
                    SNIa_number_from_this_epoch_till_this_time = function_number_SNIa_power_law(0, age_of_this_epoch, SNIa_number_prob, M_tot_of_this_epoch)
                    ejected_gas_mass_of_this_epoch += total_mass_eject_per_SNIa * SNIa_number_from_this_epoch_till_this_time
                    metal_mass_of_this_epoch += 1.2 * SNIa_number_from_this_epoch_till_this_time  # Chandrasekhar_mass = 1.44
                    Fe_mass_of_SNIa = Fe_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    O_mass_of_SNIa = O_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    Si_mass_of_SNIa = Si_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    #N_mass_of_SNIa = N_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    #C_mass_of_SNIa = C_mass_eject * SNIa_number_from_this_epoch_till_this_time
                    O_mass_of_this_epoch += O_mass_of_SNIa
                    Si_mass_of_this_epoch += Si_mass_of_SNIa
                    Fe_mass_of_this_epoch += Fe_mass_of_SNIa
                    #N_mass_of_this_epoch += N_mass_of_SNIa
                    #C_mass_of_this_epoch += C_mass_of_SNIa
                    Fe_over_H_of_an_epoch = function_element_abundunce(solar_abu_table, "Fe", "H", metal_in_gas[6], metal_in_gas[1], False)
                    observable_star_number_of_different_epochs_at_this_time.append([Fe_over_H_of_an_epoch, observable_star_number_of_a_epoch_at_a_time_step])
                    ejected_gas_mass_till_this_time += ejected_gas_mass_of_this_epoch
                    ejected_metal_mass_till_this_time += metal_mass_of_this_epoch
                    ejected_H_mass_till_this_time += H_mass_of_this_epoch
                    ejected_O_mass_till_this_time += O_mass_of_this_epoch
                    ejected_N_mass_till_this_time += N_mass_of_this_epoch
                    ejected_C_mass_till_this_time += C_mass_of_this_epoch
                    ejected_Si_mass_till_this_time += Si_mass_of_this_epoch
                    ejected_Fe_mass_till_this_time += Fe_mass_of_this_epoch
                # Goes to the next SF epoch until all SF event before this time step is accounted:
                (epoch_index) = (epoch_index + 1)
                # output of this time step
        ### yields at this time step from all SF epoch:
        ejected_gas_mass_at_this_time = ejected_gas_mass_till_this_time - ejected_gas_mass_till_last_time
        ejected_metal_mass_at_this_time = ejected_metal_mass_till_this_time - ejected_metal_mass_till_last_time
        ejected_H_mass_at_this_time = ejected_H_mass_till_this_time - ejected_H_mass_till_last_time
        ejected_O_mass_at_this_time = ejected_O_mass_till_this_time - ejected_O_mass_till_last_time
        ejected_N_mass_at_this_time = ejected_N_mass_till_this_time - ejected_N_mass_till_last_time
        ejected_C_mass_at_this_time = ejected_C_mass_till_this_time - ejected_C_mass_till_last_time
        ejected_Si_mass_at_this_time = ejected_Si_mass_till_this_time - ejected_Si_mass_till_last_time
        ejected_Fe_mass_at_this_time = ejected_Fe_mass_till_this_time - ejected_Fe_mass_till_last_time
        if total_gas_mass_at_last_time > M_tot_of_this_epoch * (outflow + 1):
            outflow_mass = M_tot_of_this_epoch * outflow
        elif total_gas_mass_at_last_time < M_tot_of_this_epoch:
            outflow_mass = 0
        else:
            outflow_mass = total_gas_mass_at_last_time - M_tot_of_this_epoch
        lockup_and_outflow_mass = M_tot_of_this_epoch + outflow_mass
        total_gas_mass_at_this_time = total_gas_mass_at_last_time - lockup_and_outflow_mass + ejected_gas_mass_at_this_time
        if total_gas_mass_at_this_time < 0:
            total_gas_mass_at_this_time = 0
            lockup_and_outflow_mass = total_gas_mass_at_last_time - ejected_gas_mass_at_this_time
            # print("SFE {} too high at time step {}".format(SFE, time_step))
        gas_mass_minimum = 0.0001
        if total_gas_mass_at_this_time < gas_mass_minimum:
            total_gas_mass_at_this_time = gas_mass_minimum
        total_metal_mass_at_this_time = total_metal_mass_in_gas_at_last_time - lockup_and_outflow_mass * \
                                        Z_gas_this_time_step + ejected_metal_mass_at_this_time
        Z_minimum = gas_mass_minimum * element_abundances_solar.function_solar_element_abundances(solar_mass_component, 'Metal')
        if total_metal_mass_at_this_time < Z_minimum:
            total_metal_mass_at_this_time = Z_minimum
        total_H_mass_at_this_time = total_H_mass_at_last_time - lockup_and_outflow_mass * (
            total_H_mass_at_last_time / total_gas_mass_at_last_time) + ejected_H_mass_at_this_time
        H_minimum = gas_mass_minimum * primary_H_mass_fraction
        if total_H_mass_at_this_time < H_minimum:
            total_H_mass_at_this_time = H_minimum

        total_O_mass_at_this_time = total_O_mass_at_last_time - lockup_and_outflow_mass * (
                total_O_mass_at_last_time / total_gas_mass_at_last_time) + ejected_O_mass_at_this_time
        O_minimum = gas_mass_minimum * element_abundances_primordial.function_element_mass_primary_fraction(solar_abu_table, "O", Z_0, Z_solar)
        if total_O_mass_at_this_time < O_minimum:
            total_O_mass_at_this_time = O_minimum
        total_N_mass_at_this_time = total_N_mass_at_last_time - lockup_and_outflow_mass * (
                total_N_mass_at_last_time / total_gas_mass_at_last_time) + ejected_N_mass_at_this_time
        N_minimum = gas_mass_minimum * element_abundances_primordial.function_element_mass_primary_fraction(solar_abu_table, "N", Z_0, Z_solar)
        if total_N_mass_at_this_time < N_minimum:
            total_N_mass_at_this_time = N_minimum
        total_C_mass_at_this_time = total_C_mass_at_last_time - lockup_and_outflow_mass * (
                total_C_mass_at_last_time / total_gas_mass_at_last_time) + ejected_C_mass_at_this_time
        C_minimum = gas_mass_minimum * element_abundances_primordial.function_element_mass_primary_fraction(solar_abu_table, "C", Z_0, Z_solar)
        if total_C_mass_at_this_time < C_minimum:
            total_C_mass_at_this_time = C_minimum
        total_Si_mass_at_this_time = total_Si_mass_at_last_time - lockup_and_outflow_mass * (
                total_Si_mass_at_last_time / total_gas_mass_at_last_time) + ejected_Si_mass_at_this_time
        Si_minimum = gas_mass_minimum * element_abundances_primordial.function_element_mass_primary_fraction(solar_abu_table, "Si", Z_0, Z_solar)
        if total_Si_mass_at_this_time < Si_minimum:
            total_Si_mass_at_this_time = Si_minimum
        total_Fe_mass_at_this_time = total_Fe_mass_at_last_time - lockup_and_outflow_mass * (
            total_Fe_mass_at_last_time / total_gas_mass_at_last_time) + ejected_Fe_mass_at_this_time
        Fe_minimum = gas_mass_minimum * element_abundances_primordial.function_element_mass_primary_fraction(solar_abu_table, "Fe", Z_0, Z_solar)
        if total_Fe_mass_at_this_time < Fe_minimum:
            total_Fe_mass_at_this_time = Fe_minimum
        if gas_infall == True:
            tau_infall = tau_infalle9 * 1e9
            # total_infall_mass = 100  # the total infall mass at t=infinity is total_infall_mass * 1e9 Msun
            # A_in = total_infall_mass / tau_infalle9 / tau_infalle9 / 1e9
            A_in = total_gas_mass / tau_infall / tau_infall
            infall_mass_till_this_time = A_in * (tau_infall ** 2 - (tau_infall * this_time + tau_infall ** 2) * math.exp(-this_time / tau_infall))
            infall_mass_at_this_time = infall_mass_till_this_time - infall_mass_till_this_time_list[-1]
            total_gas_mass_at_this_time += infall_mass_at_this_time
            total_H_mass_at_this_time += infall_mass_at_this_time * 0.75

        ##########################################
        ##### luminosity weighted abundances #####
        ##########################################

        observable_star_number_of_different_epochs_list += [observable_star_number_of_different_epochs_at_this_time]
        infall_mass_till_this_time_list += [infall_mass_till_this_time]
        (time_step) = (time_step + 1)

    ###################
    ###     end     ###
    ###################

    giant_number_Fe_over_H_sim_hist = []
    N_over_O_evolution_list = []
    C_over_O_evolution_list = []
    O_over_Fe_evolution_list = []
    Si_over_Fe_evolution_list = []
    Fe_over_H_evolution_list = []
    this_time = time_axis[-1]
    epoch_index_limit = (this_time + 1) / 10 ** 7
    if epoch_index_limit > SF_timescale_limit:
        epoch_index_limit = SF_timescale_limit
    epoch_index = 0
    while epoch_index < epoch_index_limit:
        SFR_of_this_epoch = epoch_info[epoch_index][0]
        if SFR_of_this_epoch != 0:
            mass_grid_table = epoch_info[epoch_index][4]
            metal_in_gas = epoch_info[epoch_index][12]
            Fe_over_H_of_an_epoch = function_element_abundunce(solar_abu_table, "Fe", "H", metal_in_gas[6], metal_in_gas[1], False)
            Si_over_Fe_of_an_epoch = function_element_abundunce(solar_abu_table, "Si", "Fe", metal_in_gas[5], metal_in_gas[6], False)
            O_over_Fe_of_an_epoch = function_element_abundunce(solar_abu_table, "O", "Fe", metal_in_gas[2], metal_in_gas[6], False)
            N_over_O_of_an_epoch = function_element_abundunce(solar_abu_table, "N", "O", metal_in_gas[3], metal_in_gas[2], False)
            C_over_O_of_an_epoch = function_element_abundunce(solar_abu_table, "C", "O", metal_in_gas[4], metal_in_gas[2], False)
            Fe_over_H_evolution_list.append(Fe_over_H_of_an_epoch)
            Si_over_Fe_evolution_list.append(Si_over_Fe_of_an_epoch)
            O_over_Fe_evolution_list.append(O_over_Fe_of_an_epoch)
            N_over_O_evolution_list.append(N_over_O_of_an_epoch)
            C_over_O_evolution_list.append(C_over_O_of_an_epoch)
        (epoch_index) = (epoch_index + 1)

    n__ = len(observable_star_number_of_different_epochs_list[-1])
    max_value = max(observable_star_number_of_different_epochs_list[-1][i][1] for i in range(n__))
    for i in range(n__):
        giant_number_Fe_over_H_sim_hist += round(observable_star_number_of_different_epochs_list[-1][i][1]/max_value*100) * [
            observable_star_number_of_different_epochs_list[-1][i][0]]
    xerr__ = 0.030375089141004858  # 0.096
    length__errors = len(giant_number_Fe_over_H_sim_hist)
    random_errors = np.random.normal(0, xerr__, length__errors)
    model_error_list = [sum(x) for x in zip(random_errors, giant_number_Fe_over_H_sim_hist)]
    kde_model_with_error__ = stats.gaussian_kde(model_error_list)
    output_Fe_over_H_evolution_list = []
    output_Si_over_Fe_evolution_list = []
    output_O_over_Fe_evolution_list = []
    output_N_over_O_evolution_list = []
    output_C_over_O_evolution_list = []
    output_observable_star_number_list = []
    length__Fe_over_H_evolution_list = len(Fe_over_H_evolution_list)
    k = 0
    while k < length__Fe_over_H_evolution_list:
        j = 0
        observable_star_number_at_this_output_step = 0
        while k + j < length__Fe_over_H_evolution_list and abs(Fe_over_H_evolution_list[k + j] - Fe_over_H_evolution_list[k]) < 0.01 and abs(
                O_over_Fe_evolution_list[k + j] - O_over_Fe_evolution_list[k]) < 0.01:
            observable_star_number_at_this_output_step += observable_star_number_of_different_epochs_list[-1][k + j][1]
            (j) = (j + 1)
        output_Fe_over_H_evolution_list.append(Fe_over_H_evolution_list[k+round(j/2)])
        output_Si_over_Fe_evolution_list.append(Si_over_Fe_evolution_list[k+round(j/2)])
        output_O_over_Fe_evolution_list.append(O_over_Fe_evolution_list[k+round(j/2)])
        output_N_over_O_evolution_list.append(N_over_O_evolution_list[k+round(j/2)])
        output_C_over_O_evolution_list.append(C_over_O_evolution_list[k+round(j/2)])
        output_observable_star_number_list.append(observable_star_number_at_this_output_step)
        (k) = (k + j)
    return kde_model_with_error__, all_sfr, output_Fe_over_H_evolution_list, output_Si_over_Fe_evolution_list, output_O_over_Fe_evolution_list, output_observable_star_number_list, output_N_over_O_evolution_list, output_C_over_O_evolution_list


def function_number_SNIa_power_law(last_delay_time, this_delay_time, SNIa_number_prob__, M_tot_of_this_epoch):
    SNIa_number_per_solar_mass = quad(function_SNIa_DTD, last_delay_time, this_delay_time, limit=40)[0]
    SNIa_number = M_tot_of_this_epoch * SNIa_number_per_solar_mass / 0.0013510468127287789/3*2 * SNIa_number_prob__
    # return SNIa_number * 2.2
    return SNIa_number


def function_SNIa_DTD(delay_time):
    if delay_time < 4 * 10 ** 7:  # [yr] #  2.3 * 10 ** 7 for a burst of star formation from Greggio 1983
    # if delay_time < 10 * 10 ** 7:  # [yr] #  2.3 * 10 ** 7 for a burst of star formation from Greggio 1983
        number = 0
    else:
        number = 10 ** (-4) * delay_time ** (-1) * 4  # 4.31 # Kroupa IMF  60 Myr
    return number


def function_read_lifetime(Z_select_in_table):
    #### if apply instantaneous recycling approximation ####
    if Z_select_in_table[0] == 'out':
        file_lifetime = open(
            'yield_tables/rearranged___/setllar_lifetime_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table[1]),
            'r')
        data = file_lifetime.readlines()
        mass_1 = data[3]
        lifetime_ = data[5]
        file_lifetime.close()
        mass = [float(x) for x in mass_1.split()]
        lifetime_table = [float(x) for x in lifetime_.split()]
    else:
        file_lifetime_low = open(
            'yield_tables/rearranged___/setllar_lifetime_from_portinari98/portinari98_Z={}.txt'.format(Z_select_in_table[1]),
            'r')
        data_low = file_lifetime_low.readlines()
        mass_1 = data_low[3]
        lifetime_low = data_low[5]
        file_lifetime_low.close()
        file_lifetime_high = open(
            'yield_tables/rearranged___/setllar_lifetime_from_portinari98/portinari98_Z={}.txt'.format(
                Z_select_in_table[3]),
            'r')
        data_high = file_lifetime_high.readlines()
        lifetime_high = data_high[5]
        file_lifetime_high.close()
        mass = [float(x) for x in mass_1.split()]
        lifetime_table_low = [float(x) for x in lifetime_low.split()]
        lifetime_table_high = [float(x) for x in lifetime_high.split()]
        x1 = Z_select_in_table[1]
        x2 = Z_select_in_table[2]
        x3 = Z_select_in_table[3]
        lifetime_table = [y1+(y3-y1)*(x2-x1)/(x3-x1) for y1, y3 in zip(lifetime_table_low, lifetime_table_high)]
    return mass, lifetime_table


def function_read_Mmetal(Z_select_in_table_2, log_Z_select_in_table_2):
    global mm, zz
    if Z_select_in_table_2[0] == 'out':
        Metal_eject_table = []
        file_Metal_eject = open('yield_tables__2024/rearranged___/setllar_Metal_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
            Z_select_in_table_2[1]), 'r')
        data = file_Metal_eject.readlines()
        mass_2 = data[3]
        Metal_eject_ = data[5]
        file_Metal_eject.close()
        mass = [float(x) for x in mass_2.split()]

        for i in range(len(Metal_eject_.split())):
            Metal_eject_table.append(float(Metal_eject_.split()[i]))
    else:
        file_Metal_eject = open('yield_tables__2024/rearranged___/setllar_Metal_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
            Z_select_in_table_2[1]), 'r')
        data_ = file_Metal_eject.readlines()
        mass_2 = data_[3]
        Metal_eject_low_ = data_[5]
        file_Metal_eject.close()
        mass = [float(x) for x in mass_2.split()]
        Metal_eject_table_low_ = [float(x) for x in Metal_eject_low_.split()]
        file_Metal_eject = open('yield_tables__2024/rearranged___/setllar_Metal_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
            Z_select_in_table_2[3]), 'r')
        data_ = file_Metal_eject.readlines()
        Metal_eject_high_ = data_[5]
        file_Metal_eject.close()
        Metal_eject_table_high_ = [float(x) for x in Metal_eject_high_.split()]
        Metal_eject_table_low = []
        Metal_eject_table_high = []
        for i in range(len(Metal_eject_table_high_)):
            Metal_eject_table_low.append(Metal_eject_table_low_[i])
            Metal_eject_table_high.append(Metal_eject_table_high_[i])
        x1 = log_Z_select_in_table_2[0]
        x2 = log_Z_select_in_table_2[1]
        x3 = log_Z_select_in_table_2[2]
        Metal_eject_table = [y1 + (y3 - y1) * (x2 - x1) / (x3 - x1) for y1, y3 in
                             zip(Metal_eject_table_low, Metal_eject_table_high)]
    return mass, Metal_eject_table
def function_read_M_element_H(Z_select_in_table_2, log_Z_select_in_table_2):
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_H_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[1]), 'r')
    data = file_M_eject.readlines()
    M_eject_low2 = data[5]
    file_M_eject.close()
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_H_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[3]), 'r')
    data = file_M_eject.readlines()
    M_eject_high2 = data[5]
    file_M_eject.close()
    M_eject_table_low_2 = [float(x) for x in M_eject_low2.split()]
    M_eject_table_high_2 = [float(x) for x in M_eject_high2.split()]
    M_eject_table_low = []
    M_eject_table_high = []
    for i in range(len(M_eject_table_low_2)):
        M_eject_table_low.append(M_eject_table_low_2[i])
        M_eject_table_high.append(M_eject_table_high_2[i])
    x1 = log_Z_select_in_table_2[0]
    x2 = log_Z_select_in_table_2[1]
    x3 = log_Z_select_in_table_2[2]
    if x3 == x1:
        M_eject_table = M_eject_table_high
    else:
        # interpolate yields between two metallicity
        M_eject_table = [y1 + (y3 - y1) * (x2 - x1) / (x3 - x1) for y1, y3 in
                         zip(M_eject_table_low, M_eject_table_high)]
    return M_eject_table
def function_read_M_element_He(Z_select_in_table_2, log_Z_select_in_table_2):
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_He_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[1]), 'r')
    data = file_M_eject.readlines()
    M_eject_low2 = data[5]
    file_M_eject.close()
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_He_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[3]), 'r')
    data = file_M_eject.readlines()
    M_eject_high2 = data[5]
    file_M_eject.close()
    M_eject_table_low_2 = [float(x) for x in M_eject_low2.split()]
    M_eject_table_high_2 = [float(x) for x in M_eject_high2.split()]
    M_eject_table_low = []
    M_eject_table_high = []
    for i in range(len(M_eject_table_low_2)):
        M_eject_table_low.append(M_eject_table_low_2[i])
        M_eject_table_high.append(M_eject_table_high_2[i])
    x1 = log_Z_select_in_table_2[0]
    x2 = log_Z_select_in_table_2[1]
    x3 = log_Z_select_in_table_2[2]
    if x3 == x1:
        M_eject_table = M_eject_table_high
    else:
        # interpolate yields between two metallicity
        M_eject_table = [y1 + (y3 - y1) * (x2 - x1) / (x3 - x1) for y1, y3 in
                         zip(M_eject_table_low, M_eject_table_high)]
    return M_eject_table
def function_read_M_element_O(Z_select_in_table_2, log_Z_select_in_table_2):
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_O_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[1]), 'r')
    data = file_M_eject.readlines()
    M_eject_low2 = data[5]
    file_M_eject.close()
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_O_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[3]), 'r')
    data = file_M_eject.readlines()
    M_eject_high2 = data[5]
    file_M_eject.close()
    M_eject_table_low_2 = [float(x) for x in M_eject_low2.split()]
    M_eject_table_high_2 = [float(x) for x in M_eject_high2.split()]
    M_eject_table_low = []
    M_eject_table_high = []
    for i in range(len(M_eject_table_low_2)):
        M_eject_table_low.append(M_eject_table_low_2[i])
        M_eject_table_high.append(M_eject_table_high_2[i])
    x1 = log_Z_select_in_table_2[0]
    x2 = log_Z_select_in_table_2[1]
    x3 = log_Z_select_in_table_2[2]
    if x3 == x1:
        M_eject_table = M_eject_table_high
    else:
        # interpolate yields between two metallicity
        M_eject_table = [y1 + (y3 - y1) * (x2 - x1) / (x3 - x1) for y1, y3 in
                         zip(M_eject_table_low, M_eject_table_high)]
    return M_eject_table
def function_read_M_element_N(Z_select_in_table_2, log_Z_select_in_table_2):
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_N_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[1]), 'r')
    data = file_M_eject.readlines()
    M_eject_low2 = data[5]
    file_M_eject.close()
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_N_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[3]), 'r')
    data = file_M_eject.readlines()
    M_eject_high2 = data[5]
    file_M_eject.close()
    M_eject_table_low_2 = [float(x) for x in M_eject_low2.split()]
    M_eject_table_high_2 = [float(x) for x in M_eject_high2.split()]
    M_eject_table_low = []
    M_eject_table_high = []
    for i in range(len(M_eject_table_low_2)):
        M_eject_table_low.append(M_eject_table_low_2[i])
        M_eject_table_high.append(M_eject_table_high_2[i])
    x1 = log_Z_select_in_table_2[0]
    x2 = log_Z_select_in_table_2[1]
    x3 = log_Z_select_in_table_2[2]
    if x3 == x1:
        M_eject_table = M_eject_table_high
    else:
        # interpolate yields between two metallicity
        M_eject_table = [y1 + (y3 - y1) * (x2 - x1) / (x3 - x1) for y1, y3 in
                         zip(M_eject_table_low, M_eject_table_high)]
    return M_eject_table
def function_read_M_element_C(Z_select_in_table_2, log_Z_select_in_table_2):
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_C_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[1]), 'r')
    data = file_M_eject.readlines()
    M_eject_low2 = data[5]
    file_M_eject.close()
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_C_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[3]), 'r')
    data = file_M_eject.readlines()
    M_eject_high2 = data[5]
    file_M_eject.close()
    M_eject_table_low_2 = [float(x) for x in M_eject_low2.split()]
    M_eject_table_high_2 = [float(x) for x in M_eject_high2.split()]
    M_eject_table_low = []
    M_eject_table_high = []
    for i in range(len(M_eject_table_low_2)):
        M_eject_table_low.append(M_eject_table_low_2[i])
        M_eject_table_high.append(M_eject_table_high_2[i])
    x1 = log_Z_select_in_table_2[0]
    x2 = log_Z_select_in_table_2[1]
    x3 = log_Z_select_in_table_2[2]
    if x3 == x1:
        M_eject_table = M_eject_table_high
    else:
        # interpolate yields between two metallicity
        M_eject_table = [y1 + (y3 - y1) * (x2 - x1) / (x3 - x1) for y1, y3 in
                         zip(M_eject_table_low, M_eject_table_high)]
    return M_eject_table
def function_read_M_element_Si(Z_select_in_table_2, log_Z_select_in_table_2):
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_Si_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[1]), 'r')
    data = file_M_eject.readlines()
    M_eject_low2 = data[5]
    file_M_eject.close()
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_Si_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[3]), 'r')
    data = file_M_eject.readlines()
    M_eject_high2 = data[5]
    file_M_eject.close()
    M_eject_table_low_2 = [float(x) for x in M_eject_low2.split()]
    M_eject_table_high_2 = [float(x) for x in M_eject_high2.split()]
    M_eject_table_low = []
    M_eject_table_high = []
    for i in range(len(M_eject_table_low_2)):
        M_eject_table_low.append(M_eject_table_low_2[i])
        M_eject_table_high.append(M_eject_table_high_2[i])
    x1 = log_Z_select_in_table_2[0]
    x2 = log_Z_select_in_table_2[1]
    x3 = log_Z_select_in_table_2[2]
    if x3 == x1:
        M_eject_table = M_eject_table_high
    else:
        # interpolate yields between two metallicity
        M_eject_table = [y1 + (y3 - y1) * (x2 - x1) / (x3 - x1) for y1, y3 in
                         zip(M_eject_table_low, M_eject_table_high)]
    return M_eject_table
def function_read_M_element_Fe(Z_select_in_table_2, log_Z_select_in_table_2):
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_Fe_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[1]), 'r')
    data = file_M_eject.readlines()
    M_eject_low2 = data[5]
    file_M_eject.close()
    file_M_eject = open('yield_tables__2024/rearranged___/setllar_Fe_eject_mass_from_Limongi_M000/Limongi_M000_Z={}.txt'.format(
        Z_select_in_table_2[3]), 'r')
    data = file_M_eject.readlines()
    M_eject_high2 = data[5]
    file_M_eject.close()
    M_eject_table_low_2 = [float(x) for x in M_eject_low2.split()]
    M_eject_table_high_2 = [float(x) for x in M_eject_high2.split()]
    M_eject_table_low = []
    M_eject_table_high = []
    for i in range(len(M_eject_table_low_2)):
        M_eject_table_low.append(M_eject_table_low_2[i])
        M_eject_table_high.append(M_eject_table_high_2[i])
    x1 = log_Z_select_in_table_2[0]
    x2 = log_Z_select_in_table_2[1]
    x3 = log_Z_select_in_table_2[2]
    if x3 == x1:
        M_eject_table = M_eject_table_high
    else:
        # interpolate yields between two metallicity
        M_eject_table = [y1 + (y3 - y1) * (x2 - x1) / (x3 - x1) for y1, y3 in
                         zip(M_eject_table_low, M_eject_table_high)]
    return M_eject_table


def function_get_target_mass_in_range(lower_mass_limit, alpha1, alpha3, integrate_mass, len_mass_grid_table2, M_element_table):
    return quad(integrator_for_function_get_target_mass_in_range2, lower_mass_limit, 50.1,
                              args=(M_element_table, alpha1, alpha3, integrate_mass, len_mass_grid_table2), limit=40)[0]

def integrator_for_function_get_target_mass_in_range2(initial_mass, Mtarget_table, alpha1, alpha3, integrate_mass, len_mass_grid_table2):
    mass_grid_table2 = [0.08, 1.0, 1.25, 1.5, 1.75, 1.9, 2.0, 2.25, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 120.0, 150]
    if initial_mass > 50.1:
        return 0
    elif initial_mass > 1:
        number = initial_mass ** (-alpha3) / integrate_mass
    elif initial_mass > 0.08:
        number = initial_mass ** (-alpha1-1) / integrate_mass
    else:
        number = 0.5**(-alpha1-1) * (initial_mass / 0.5)**(-alpha1) / integrate_mass
    high = len_mass_grid_table2 - 1
    low = 0
    while low < high:
        mid = (low + high) // 2
        if mass_grid_table2[mid] < initial_mass:
            low = mid + 1
        else:
            high = mid
    if low == 0:
        mass_eject_fraction_per_star = Mtarget_table[low]
    elif low == len_mass_grid_table2:
        mass_eject_fraction_per_star = Mtarget_table[low - 1]
    else:
        Mtarget_table_low = Mtarget_table[low]
        mass_grid_table_2_low = mass_grid_table2[low]
        mass_eject_fraction_per_star = Mtarget_table_low + (Mtarget_table[low - 1] - Mtarget_table_low) * (
                initial_mass - mass_grid_table_2_low) / (mass_grid_table2[low - 1] - mass_grid_table_2_low)
    integrator = number * mass_eject_fraction_per_star * initial_mass
    return integrator

def function_get_target_mass_in_range_metal(lower_mass_limit, alpha1, alpha3, integrate_mass, len_mass_grid_table2):
    return quad(integrator_for_function_get_target_mass_in_range, lower_mass_limit, 50.1,
                              args=(M_element_table[0], alpha1, alpha3, integrate_mass, len_mass_grid_table2), limit=40)[0]
def function_get_target_mass_in_range_H(lower_mass_limit, alpha1, alpha3, integrate_mass, len_mass_grid_table2):
    return quad(integrator_for_function_get_target_mass_in_range, lower_mass_limit, 50.1,
                              args=(M_element_table[1], alpha1, alpha3, integrate_mass, len_mass_grid_table2), limit=40)[0]
def function_get_target_mass_in_range_He(lower_mass_limit, alpha1, alpha3, integrate_mass, len_mass_grid_table2):
    return quad(integrator_for_function_get_target_mass_in_range, lower_mass_limit, 50.1,
                              args=(M_element_table[2], alpha1, alpha3, integrate_mass, len_mass_grid_table2), limit=40)[0]
def function_get_target_mass_in_range_O(lower_mass_limit, alpha1, alpha3, integrate_mass, len_mass_grid_table2):
    return quad(integrator_for_function_get_target_mass_in_range, lower_mass_limit, 50.1,
                              args=(M_element_table[3], alpha1, alpha3, integrate_mass, len_mass_grid_table2), limit=40)[0]
def function_get_target_mass_in_range_N(lower_mass_limit, alpha1, alpha3, integrate_mass, len_mass_grid_table2):
    return quad(integrator_for_function_get_target_mass_in_range, lower_mass_limit, 50.1,
                              args=(M_element_table[4], alpha1, alpha3, integrate_mass, len_mass_grid_table2), limit=40)[0]
def function_get_target_mass_in_range_C(lower_mass_limit, alpha1, alpha3, integrate_mass, len_mass_grid_table2):
    return quad(integrator_for_function_get_target_mass_in_range, lower_mass_limit, 50.1,
                              args=(M_element_table[5], alpha1, alpha3, integrate_mass, len_mass_grid_table2), limit=40)[0]
def function_get_target_mass_in_range_Si(lower_mass_limit, alpha1, alpha3, integrate_mass, len_mass_grid_table2):
    return quad(integrator_for_function_get_target_mass_in_range, lower_mass_limit, 50.1,
                              args=(M_element_table[6], alpha1, alpha3, integrate_mass, len_mass_grid_table2), limit=40)[0]
def function_get_target_mass_in_range_Fe(lower_mass_limit, alpha1, alpha3, integrate_mass, len_mass_grid_table2):
    return quad(integrator_for_function_get_target_mass_in_range, lower_mass_limit, 50.1,
                              args=(M_element_table[7], alpha1, alpha3, integrate_mass, len_mass_grid_table2), limit=40)[0]

def integrator_for_function_get_target_mass_in_range(initial_mass, Mtarget_table, alpha1, alpha3, integrate_mass, len_mass_grid_table2):
    if initial_mass > 50.1:
        return 0
    elif initial_mass > 1:
        number = initial_mass ** (-alpha3) / integrate_mass
    elif initial_mass > 0.08:
        number = initial_mass ** (-alpha1-1) / integrate_mass
    else:
        number = 0.5 ** (-alpha1-1) * (initial_mass / 0.5)**(-alpha1) / integrate_mass
    high = len_mass_grid_table2 - 1
    low = 0
    while low < high:
        mid = (low + high) // 2
        if mass_grid_table2[mid] < initial_mass:
            low = mid + 1
        else:
            high = mid
    if low == 0:
        mass_eject_fraction_per_star = Mtarget_table[low]
    elif low == len_mass_grid_table2:
        mass_eject_fraction_per_star = Mtarget_table[low - 1]
    else:
        Mtarget_table_low = Mtarget_table[low]
        mass_grid_table_2_low = mass_grid_table2[low]
        mass_eject_fraction_per_star = Mtarget_table_low + (Mtarget_table[low - 1] - Mtarget_table_low) * (
                initial_mass - mass_grid_table_2_low) / (mass_grid_table2[low - 1] - mass_grid_table_2_low)
    integrator = number * mass_eject_fraction_per_star * initial_mass
    return integrator


def function_element_abundunce(solar_abu_table, element_1_name, element_2_name, metal_1_mass, metal_2_mass, instant_ejection):
    if metal_2_mass == 0:
        if metal_1_mass == 0:
            metal_1_over_2 = -6
        elif metal_1_mass > 0:
            metal_1_over_2 = None
        elif metal_1_mass < 0:
            if instant_ejection is False:
                print("Warning: current {} mass < 0. See galevo.py".format(element_1_name))
            metal_1_over_2 = -6
    elif metal_2_mass < 0:
        if instant_ejection is False:
            print("Warning: current {} mass < 0. See galevo.py".format(element_2_name))
        if metal_1_mass == 0:
            metal_1_over_2 = -6
        elif metal_1_mass > 0:
            metal_1_over_2 = -6
        elif metal_1_mass < 0:
            if instant_ejection is False:
                print("Warning: current {} mass < 0. See galevo.py".format(element_1_name))
            metal_1_over_2 = None
    else:
        if metal_1_mass == 0:
            metal_1_over_2 = None
        elif metal_1_mass < 0:
            if instant_ejection is False:
                print("Warning: current {} mass < 0. See galevo.py".format(element_1_name))
            metal_1_over_2 = None
        else:
            solar_metal_1_logarithmic_abundances = element_abundances_solar.function_solar_element_abundances(
                solar_abu_table, element_1_name)
            solar_metal_2_logarithmic_abundances = element_abundances_solar.function_solar_element_abundances(
                solar_abu_table, element_2_name)
            metal_1_element_weight = element_weight_table.function_element_weight(element_1_name)
            metal_2_element_weight = element_weight_table.function_element_weight(element_2_name)
            metal_1_over_2 = math.log(metal_1_mass / metal_2_mass / metal_1_element_weight * metal_2_element_weight, 10) - (solar_metal_1_logarithmic_abundances - solar_metal_2_logarithmic_abundances)
    return metal_1_over_2


def function_get_avaliable_Z():
    # list 1
    file_names_setllar_lifetime_from_str_yield_table = os.listdir(
        'yield_tables/rearranged___/setllar_lifetime_from_portinari98')
    Z_table_list = []
    for name in file_names_setllar_lifetime_from_str_yield_table:
        length_file_name = len(name)
        i = 0
        i_start = 0
        i_end = 0
        while i < length_file_name:
            if name[i] == '=':
                i_start = i
            if name[i] == '.':
                i_end = i
            (i) = (i + 1)
        i = i_start + 1
        Z = ''
        while i < i_end:
            Z += name[i]
            (i) = (i + 1)
        Z_table_list += [float(Z)]
    sorted_Z_table_list = sorted(Z_table_list)
    # list 2
    file_names_setllar_lifetime_from_str_yield_table = os.listdir(
        'yield_tables__2024/rearranged___/setllar_Mg_eject_mass_from_Limongi_M000')
    Z_table_list_2 = []
    for name in file_names_setllar_lifetime_from_str_yield_table:
        length_file_name = len(name)
        i = 0
        i_start = 0
        i_end = 0
        while i < length_file_name:
            if name[i] == '=':
                i_start = i
            if name[i] == '.':
                i_end = i
            (i) = (i + 1)
        i = i_start + 1
        Z = ''
        while i < i_end:
            Z += name[i]
            (i) = (i + 1)
        if Z != '':
            Z_table_list_2 += [float(Z)]
    sorted_Z_table_list_2 = sorted(Z_table_list_2)
    return sorted_Z_table_list, sorted_Z_table_list_2


def function_select_metal(Z, Z_table_list):
    # the list for stellar lifetime is
    # [0.0004, 0.0008, 0.0012, 0.0016, 0.002, 0.0024, 0.0028, 0.0032, 0.0036, 0.004, 0.008, 0.012]
    # the list for stellar metallicity is
    # [0.0004, 0.004, 0.008, 0.0127] or [0, 0.004, 0.02] for Kobayashi2006 massive star table
    if Z <= Z_table_list[0]:
        Z_select__ = Z_table_list[0]
        return 'out', Z_select__, Z_select__, Z_select__
        # The 'out' flag means the current gas metallicity is outside the range of provided stellar yield table.
    elif Z >= Z_table_list[-1]:
        Z_select__ = Z_table_list[-1]
        return 'out', Z_select__, Z_select__, Z_select__
    else:
        i = 1
        while i < len(Z_table_list):
            if Z < Z_table_list[i]:
                Z_select__low = Z_table_list[i - 1]
                Z_select__high = Z_table_list[i]
                return 'in', Z_select__low, Z, Z_select__high
            (i) = (i + 1)


def fucntion_mass_boundary(time, mass, lifetime):
    # length_list_lifetime = len(lifetime)
    # length_list_lifetime = 2901
    # x = length_list_lifetime // 2
    x = 1450
    # loop_number_fucntion_mass_boundary = math.ceil(math.log2(length_list_lifetime))
    loop_number_fucntion_mass_boundary = 12
    if lifetime[x] == time:
        mass_boundary = mass[x]
    else:
        i = 0
        low = 0
        # high = length_list_lifetime
        high = 2901
        while i < loop_number_fucntion_mass_boundary:
            if lifetime[x] > time:
                low = x
                x = x + (high - x) // 2
            else:
                high = x
                x = x - (x - low) // 2
            (i) = (i + 1)
        # if x == length_list_lifetime - 1:
        if x == 2900:
            mass_boundary = mass[x]
        else:
            if lifetime[x - 1] > time > lifetime[x]:
                x = x - 1
            mass_boundary = mass[x] + (mass[x + 1] - mass[x]) * (lifetime[x] - time) / (lifetime[x] - lifetime[x + 1])
    return mass_boundary


if __name__ == '__main__':

    def unnormalized_mass_function(mass, alpha1, alpha3):
        if mass > 50.1:
            return 0
        elif mass > 1:
            return mass ** (1-alpha3)
        elif mass > 0.08:
            return mass ** (1-alpha1 - 1)
        else:
            return 0.5 ** (1-alpha1 - 1) * (mass / 0.5) ** (1-alpha1)


    def average_list(lst):
        return sum(lst) / len(lst)

    ######################################

    # import cProfile
    # cProfile.run('galaxy_evol()', sort='tottime')
    import matplotlib.pyplot as plt
    import time
    start_time = time.time()

    # log_prob -2.0044127568870516, SFE 0.09532987740311372, out 0.05926706797683223, tau_in 5.199789554313394, alpha1 0.644202050415464, alpha3 2.3007932997065152, Time: 49 m 18 s
    maximum_core_number = 25
    SFT = 800
    SFE = 0.12
    outflow = 8
    tau_infalle9 = 0.4
    alpha1__ = 1.08
    alpha3__ = 2

    integrate_mass = quad(unnormalized_mass_function, 0.08, 50.1, args=(alpha1__, alpha3__), limit=50)[0]

    kde_model_with_error, all_sfr, Fe_over_H_evolution_list, Si_over_Fe_evolution_list, O_over_Fe_evolution_list, observable_star_number_list, N_over_O_evolution_list, C_over_O_evolution_list = \
        galaxy_evol(logOGM=11, Z_0=1e-9,
                    solar_mass_component="Asplund2009_mass",
                    SF_timescale_limit=SFT, SFE=SFE,  # 0.001 0.0005 0.0002 (0.003 canonical)
                    solar_abu_table='Asplund2009',
                    outflow=outflow,  # 2697 outflow_mass = SN_number_per_century[-1]*1e5 * outflow
                    tau_infalle9=tau_infalle9, original_gas_mass_fraction=1/1000,
                    alpha1=alpha1__, alpha3=alpha3__, integrate_mass=integrate_mass, maximum_core_number = maximum_core_number)

    GCE_time = time.time()
    computation_time_seconds = round((GCE_time - start_time), 2)
    minutes, seconds = divmod(computation_time_seconds, 60)
    hours, minutes = divmod(minutes, 60)
    print("SFT: {}*10Myr. SFE={}. outflow={}. tau_infalle9={}. Use {} cores. GCE_time: {} m {} s".format(SFT, SFE, outflow, tau_infalle9, maximum_core_number, minutes, seconds))

    print(len(all_sfr), len(Fe_over_H_evolution_list), len(O_over_Fe_evolution_list))
    # print(all_sfr)
    # print(Fe_over_H_evolution_list)
    # print(Si_over_Fe_evolution_list)
    # print(O_over_Fe_evolution_list)

    ### import data from Bensby 2014
    from astropy.io import fits
    filename = "MW_sim/MW_data/J_A+A_562_A71_tablec3.dat.fits"
    with fits.open(filename) as hdul:
        data = hdul[1].data
        Fe_over_H_MW = []
        O_over_Fe_MW = []
        Si_over_Fe_MW = []
        for i in data:
            if i[10] < 1 and i[4] > 4.1 and i[34] > 0 and i[30] > 0:
            # if i[10] < 1 and i[4] > 4.1 and i[34] > 0 and i[30] > 0 and i[-10] > 2:
                # if mass[i] < 1 and logg[i]>4.1 and e_si_fe[i]>0 and e_o_fe[i]>0 and td/d<0.5 (thin disc stars):
                Fe_over_H_MW.append(i[16])
                O_over_Fe_MW.append(i[17])
                Si_over_Fe_MW.append(i[21])
    Si_over_Fe_error = 0.074047619047619
    O_over_Fe_error = 0.1694444444444444
    Fe_over_H_error = 0.0696031746031746

    #load data from Israelian 2004, https://www.aanda.org/articles/aa/pdf/2004/26/aa0132-04.pdf
    mw_data = pd.read_csv("MW_sim/MW_data/Israelian2004_metal-poor-stars.txt", delimiter="|")
    FeH_israel = mw_data['Fe/H']
    NH_israel = mw_data['N/H']
    OH_israel = mw_data['O/H']
    NO_israel = mw_data['N/O']
    kde_obs_israel = stats.gaussian_kde(FeH_israel)

    sampled_data = pd.read_csv("MW_sim/MW_data/Xiang2019_filtered_sampled_data.txt", delimiter="|")
    FeH_xiang19 = sampled_data['    [Fe/H]']
    OFe_xiang19 = sampled_data['    [O/Fe]']
    CFe_xiang19 = sampled_data['    [C/Fe]']
    NFe_xiang19 = sampled_data['    [N/Fe]']
    SiFe_xiang19 = sampled_data['   [Si/Fe]']
    # e_FeH_xiang19 = sampled_data['    e_[Fe/H]']
    # e_OFe_xiang19 = sampled_data['   e_[O/Fe]']
    # e_SiFe_xiang19 = sampled_data['  e_[Si/Fe]']
    # print(sum(e_FeH_xiang19)/len(e_FeH_xiang19))  # 0.030375089141004858
    # print(sum(e_OFe_xiang19)/len(e_OFe_xiang19))  # 0.06767662884927061
    # print(sum(e_SiFe_xiang19)/len(e_SiFe_xiang19))  # 0.0483514262560778
    kde_obs_xiang19 = stats.gaussian_kde(FeH_xiang19)  # Kernel Density Estimation (KDE)
    # Fe_over_H_error = 0.030375089141004858  # 0.06973
    # O_over_Fe_error = 0.06767662884927061  # 0.15573
    # Si_over_Fe_error = 0.0483514262560778  # 0.06253
    # length_O = len(O_over_Fe_MW)
    # length_Si = len(Si_over_Fe_MW)
    # plt.plot(kde__x, kde_obs(kde__x), c="b", lw=1.3, label='Observations\' kernel-density estimate')
    # plt.show()

    plt.figure(0, figsize=(6, 4))
    time_list_Gyr = np.linspace(0, len(all_sfr) / 100, len(all_sfr))
    plt.plot(time_list_Gyr, all_sfr, alpha=0.5)
    plt.scatter(time_list_Gyr, all_sfr, alpha=0.5)
    plt.xlabel("Time [Gyr]")
    plt.ylabel(r"Arbitrary normalized SFR [$M_\odot$/yr]")
    plt.savefig("MW_sim/plots/SFH.png", dpi=300)
    plt.tight_layout()

    # time_list_Gyr = np.linspace(0, len(all_sfr) / 100, len(N_over_C_evolution_list))
    # plt.figure(1, figsize=(6, 4))
    # plt.plot(time_list_Gyr, N_over_C_evolution_list, alpha=0.5)
    # plt.scatter(time_list_Gyr, N_over_C_evolution_list, alpha=0.5)
    # plt.xlabel("Time [Gyr]")
    # plt.ylabel(r"[$^{17}$O/$^{18}$O]")
    # plt.tight_layout()

    plt.figure(10, figsize=(6, 4))
    ### plot data
    plt.scatter(Fe_over_H_MW, O_over_Fe_MW, s=10, color='r', alpha=0.2)
    ### plot data error
    # plt.scatter([0.35], [0.75], s=20, color='tab:blue', alpha=0.1)
    # plt.errorbar([0.35], [0.75], yerr=O_over_Fe_error, xerr=Fe_over_H_error, color='tab:blue', fmt='none', lw=0.5, zorder=0)
    ### plot model
    plt.plot(Fe_over_H_evolution_list, O_over_Fe_evolution_list, alpha=0.5)
    plt.scatter(Fe_over_H_evolution_list, O_over_Fe_evolution_list, alpha=0.5)
    plt.scatter(FeH_xiang19, OFe_xiang19, alpha=0.1)
    plt.xlabel("[Fe/H]")
    plt.ylabel("[O/Fe]")
    plt.xlim(-4, 1)
    plt.ylim(-1, 2)
    plt.savefig("MW_sim/plots/OFe-FeH.png", dpi=300)
    plt.tight_layout()

    plt.figure(11, figsize=(6, 4))
    ### plot data
    plt.scatter(Fe_over_H_MW, Si_over_Fe_MW, s=10, color='r', alpha=0.2)
    ### plot data error
    # plt.scatter([0.35], [0.42], s=20, color='tab:blue', alpha=0.1, label='MW Bensby14')
    # plt.errorbar([0.35], [0.42], yerr=Si_over_Fe_error, xerr=Fe_over_H_error, color='tab:blue', fmt='none', lw=0.5, zorder=0)
    ### plot model
    plt.plot(Fe_over_H_evolution_list, Si_over_Fe_evolution_list, alpha=0.5)
    plt.scatter(Fe_over_H_evolution_list, Si_over_Fe_evolution_list, alpha=0.5)
    plt.scatter(FeH_xiang19, SiFe_xiang19, alpha=0.1)
    plt.xlabel("[Fe/H]")
    plt.ylabel("[Si/Fe]")
    plt.xlim(-4, 1)
    plt.ylim(-1, 1.5)
    plt.savefig("MW_sim/plots/SiFe-FeH.png", dpi=300)
    plt.tight_layout()

    # plt.figure(12, figsize=(6, 4))
    # plt.plot(Fe_over_H_evolution_list, N_over_C_evolution_list, alpha=0.5)
    # plt.scatter(Fe_over_H_evolution_list, N_over_C_evolution_list, alpha=0.5)
    # plt.xlabel("[Fe/H]")
    # plt.ylabel(r"[$^{17}$O/$^{18}$O]")
    # plt.tight_layout()
    kde__x = np.linspace(-3, 0.5, 100)
    plt.figure(12, figsize=(6, 4))
    plt.hist(Fe_over_H_MW, density=True, bins=25, color='k', histtype='step', linestyle='--', label='MW Bensby14')
    plt.hist(FeH_israel, density=True, bins=25, color='orange', histtype='step', linestyle='--', label='MW Israelian04')
    plt.plot(kde__x, kde_model_with_error(kde__x), alpha=0.5)
    plt.plot(kde__x, kde_obs_xiang19(kde__x), ls='dashed', alpha=0.5)
    plt.xlabel("[Fe/H]")
    plt.ylabel("PDF")
    plt.xlim(-4, 1)
    plt.ylim(0, 2)
    plt.legend()
    plt.savefig("MW_sim/plots/PDF-FeH.png", dpi=300)
    plt.tight_layout()

    O_over_H_evolution_list = [8.69 + a + b for a, b in zip(Fe_over_H_evolution_list, O_over_Fe_evolution_list)]
    NO_xiang19 = [a - b for a, b in zip(NFe_xiang19, OFe_xiang19)]
    CO_xiang19 = [a - b for a, b in zip(CFe_xiang19, OFe_xiang19)]
    OH_xiang19 = [8.69 + a + b for a, b in zip(OFe_xiang19, FeH_xiang19)]
    OH_israel = 8.69 + OH_israel

    plt.figure(13, figsize=(6, 4))
    ### plot data
    plt.scatter(OH_xiang19, NO_xiang19, color='orange', alpha=0.3)
    plt.scatter(OH_israel, NO_israel, color='red', alpha=0.3)
    ### plot model
    plt.plot(O_over_H_evolution_list, N_over_O_evolution_list, alpha=0.5, color='tab:blue')
    plt.scatter(O_over_H_evolution_list, N_over_O_evolution_list, alpha=0.5, color='tab:blue')
    plt.xlabel("12+log(O/H)")
    plt.ylabel("[N/O]")
    # plt.xlim(-4, 1)
    # plt.ylim(-1, 1.5)
    plt.savefig("MW_sim/plots/NO-12logOH.png", dpi=300)
    plt.tight_layout()
    
    plt.figure(14, figsize=(6, 4))
    ### plot data
    plt.scatter([], [], s=10, color='r', alpha=0.2)
    ### plot model
    plt.plot(O_over_H_evolution_list, C_over_O_evolution_list, alpha=0.5)
    plt.scatter(O_over_H_evolution_list, C_over_O_evolution_list, alpha=0.5)
    plt.scatter(OH_xiang19, CO_xiang19, alpha=0.1)
    plt.xlabel("12+log(O/H)")
    plt.ylabel("[C/O]")
    # plt.xlim(-4, 1)
    # plt.ylim(-1, 1.5)
    plt.tight_layout()
    plt.savefig("MW_sim/plots/CO-12logOH.png", dpi=300)
    plt.show()