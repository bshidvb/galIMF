# Python3 code

# An example


import galevo_nitrogen as galevo

Log_SFR = 3
location = 0
skewness = 20
sfr_tail = 5
SFEN = 0.66

galevo.generate_SFH('boxy', Log_SFR, SFEN, sfr_tail, skewness, location)
# galevo.galaxy_evol(imf='Kroupa', STF=0.5, SFEN=SFEN, Z_0=0.015*1e-16, solar_mass_component="Asplund2009_mass",
#                 str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=150,
#                 time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
#                 SFH_model='provided', SFE=0.04, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
#                 solar_abu_table='Asplund2009',
#                 high_time_resolution=None, plot_show=None, plot_save=None, outflow=100, check_igimf=None)

# Bekki input
galevo.galaxy_evol(imf='Kroupa', STF=0.097, SFEN=SFEN, Z_0=0.015*1e-16, solar_mass_component="Asplund2009_mass",
                str_yield_table='Kobayashi06', IMF_name='Kroupa', steller_mass_upper_bound=120,
                time_resolution_in_Myr=1, mass_boundary_observe_low=1.5, mass_boundary_observe_up=8,
                SFH_model='provided', SFE=0.04, SNIa_ON=True, SNIa_yield_table='Iwamoto1999',
                solar_abu_table='Asplund2009',
                high_time_resolution=None, plot_show=None, plot_save=None, outflow=100, check_igimf=None)
    