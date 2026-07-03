#!/bin/bash
# chmod +x push_origin-master.sh # Give execute permission to the script
# sh push_origin-master.sh # To run the script

git add plotting/igimf_animation/

git add plotting/igimf_evolution.ipynb

git commit -m "Created function that plots evolution of IGIMF"

git add IMFs_playground/kroupa-imfs-different-a modified.py

git commit -m "Created function to plot wild variations of Kroupa IMF"

git add figs/yield_plots/

git add plot_yields_LC18-N:O.ipynb

git add IMFs_playground/

git commit -m "Made new plots"

git add figs/galevo/

git add simulation_results_from_galaxy_evol/

git commit -m "Found new models"

git add galaxy_evolution_nitrogen.py

git add galevo_nitrogen.py

git add IMFs/

git commit -m "Did parameter modifications"

git push origin master
