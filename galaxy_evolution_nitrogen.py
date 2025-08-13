import galevo_nitrogen as galevo

    # For gas_mass_dependent SFH model, we need to provide the gas mass evolution.
    # This is an example of how to generate a gas mass evolution.
Log_SFR_list = [1.3]  # 20 M_sun/yr (Tacchella et al. 2023)
SFEN_list = [15]       # 70 ± 40 Myr (Tacchella et al. 2023)
STF_list = [0.04]
SFE_list = [0.004]
tau_infall_list = [0.008] 
# Log_SFR = 1.5 # for 20 M_sun/yr - real value
location = 0
skewness = 20
sfr_tail = 5
# SFEN = 5 # for 66 Myr

for SFEN in SFEN_list:
    for Log_SFR in Log_SFR_list:
        for STF in STF_list:
            for SFE in SFE_list:
                for tau_infall in tau_infall_list:
                    galevo.generate_SFH('flat', Log_SFR, SFEN, sfr_tail, skewness, location)
                    galevo.galaxy_evol(imf='Kroupa', STF=STF, SFEN=SFEN, Z_0=1e-6, solar_mass_component="Asplund2009_mass",
                                    str_yield_table='Limongi_R300', IMF_name='Kroupa', steller_mass_upper_bound=150, #try 100
                                    time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
                                    SFH_model='gas_mass_dependent', SFE=SFE, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
                                    solar_abu_table='Asplund2009',
                                    high_time_resolution=None, plot_show=None, plot_save=True, outflow=200, check_igimf=None, tau_infalle9=tau_infall)

# Bekki input
# galevo.galaxy_evol(imf='Kroupa', STF=0.5, SFEN=SFEN, Z_0=0.015*1e-6, solar_mass_component="Asplund2009_mass",
#                 str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=120,
#                 time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
#                 SFH_model='provided', SFE=0.04, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
#                 solar_abu_table='Asplund2009',
#                 high_time_resolution=None, plot_show=None, plot_save=None, outflow=100, check_igimf=None)