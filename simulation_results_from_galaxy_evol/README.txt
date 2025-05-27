Notes on structure of files:

16.5/ - alpha = 2.35, variations in STF (0.09, 0.15), SFR (10, 50, 100) and SFEN (10, 15)
      - alpha = 1.55, variations in STF (0.09, 0.15), SFR (10, 50, 100) and SFEN (10, 15)
      - alpha = 1.15, variations in STF (0.09, 0.15), SFR (10, 50, 100) and SFEN (10, 15)
      - fixed parameters: Z_0=1e-6, no outflows, flat SFH, str_yield_table='Kobayashi06' (Kob), SFH_model='provided', steller_mass_upper_bound=25

19.5/ - variation of STF, fixed SFEN = 15, SFR = 15, alpha = 1.15

20.5/ - alpha = 2.35, variations in STF (0.09, 0.15), SFR (10, 50, 100) and SFEN (10, 15)
      - alpha = 1.55, variations in STF (0.09, 0.15), SFR (10, 50, 100) and SFEN (10, 15)
      - alpha = 1.15, variations in STF (0.09, 0.15), SFR (10, 50, 100) and SFEN (10, 15)
      - fixed parameters: Z_0=1e-6, no outflows, flat SFH, str_yield_table='Limongi_R000' (Lim), SFH_model='provided', steller_mass_upper_bound=25

21.5/ - alpha = 2.35, variations in STF (0.09, 0.12, 0.15, 0.18, 0.24), SFEN fixed at 150, SFR (15, 20, 30)
      - alpha = 1.55, variations in STF (0.09, 0.12, 0.15, 0.18, 0.24), SFEN fixed at 150, SFR (15, 20, 30)
      - alpha = 1.15, variations in STF (0.09, 0.12, 0.15, 0.18, 0.24), SFEN fixed at 150, SFR (15, 20, 30)
      - fixed parameters: Z_0=1e-6, no outflows, flat SFH, str_yield_table='Limongi_R000' (Lim), SFH_model='provided', steller_mass_upper_bound=25

23.5/ - trying lower SFEN, modifying steller_mass_upper_bound = 100
      - variation of alpha (1.15, 1.55, 2.35), LC18 yields, no outflow
      - SFR fixed at 100, SFEN (2.5, 5, 10), STF (0.09, 0.12, 0.15)

experiment/ - testing whether steller_mass_upper_bound works right
            - same setup as 23.5/ but with steller_mass_upper_bound=2 (=1 giving errors)

24.5/ - lower STF (0.02, 0.05, 0.07, 0.09), variation of SFEN (25, 50, 100) and SFR (15, 30, 50, 100)
      - LC18 yields, no outflows, low metallicity Z_0=1e-6,  steller_mass_upper_bound (100)
      - trying different alpha
      
27.5/ - lower STF (0.02, 0.05, 0.07, 0.09), variation of SFEN (25, 50, 100) and SFR (15, 30, 50, 100)
      - LC18 yields, no outflows, low metallicity Z_0=1e-6,  steller_mass_upper_bound (150)
      - trying different alpha

27.5/mod_SFH/ - lower STF (0.02, 0.05, 0.07, 0.09), 20-step SFH, till 100, sfen = 100
              - LC18 yields, no outflows, low metallicity Z_0=1e-6,  steller_mass_upper_bound (150)
              - trying different alpha

28.5/ - lower STF (0.02, 0.05, 0.07, 0.09), variation of SFEN (25, 50, 100) and SFR (100)
      - LC18 yields, no outflows, low metallicity Z_0=1e-6,  steller_mass_upper_bound (150)
      - SFH_model='gas_mass_dependent', SFE=0.04
      - trying different alpha