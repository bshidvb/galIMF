#!/bin/bash
# chmod +x push_origin-master.sh # Give execute permission to the script
# sh push_origin-master.sh # To run the script

git add plotting/plot_stellar_lifetimes_comparison.py

git add plotting/plots_for-thesis.ipynb

git commit -m "Created new plotting functions for thesis and stellar lifetimes from yields"

git add plotting/plots_nitrogen_correct.ipynb

git commit -m "Marked filepaths to solutions"

git add plotting/

git commit -m "Modified existing plotting functions"

git add yield_tables__2024/rearranged___/setllar_lifetime_from_Limongi_M300

git commit -m "Generated lifetimes for Limongi M300 yields"

git add simulation_results_from_galaxy_evol/

git commit -m "Generated new solutions for thesis"

git add sfh_comparison.ipynb

git commit -m "Created plotting functions for SFHs models"

git add read_yield_table_output_fraction_2024.py

git commit -m "Modified input to generate yields"

git add galevo_nitrogen.py

git add galaxy_evolution_nitrogen.py

git commit -m "Modified functions to get new solutions"

git push origin master
