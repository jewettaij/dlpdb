#!/bin/sh

# Different statistics are generated depending upon the desired 
# resolution of the coarse-grained model.
# Specifically: do we want one base pair of the coarse-grained 
# model to represent 1, 2, or 3 base pairs of real DNA?
# For each of these choices, we invoke different scripts for 
# and copy the results into separate directories.

mkdir statistics_keeping_every_base_pair
cd statistics_keeping_every_base_pair
cp -f ../pdbs_FINAL_USE_THESE_FOR_ANALYSIS/*.pdb .
./calc_distances_angles.sh
rm -f *.pdb
cd ../

mkdir statistics_keeping_every_2nd_base_pair
cd statistics_keeping_every_2nd_base_pair
cp -f ../pdbs_FINAL_USE_THESE_FOR_ANALYSIS/*.pdb .
./calc_distances_angles.sh
rm -f *.pdb
cd ../

mkdir statistics_keeping_every_3rd_base_pair
cd statistics_keeping_every_3rd_base_pair
cp -f ../pdbs_FINAL_USE_THESE_FOR_ANALYSIS/*.pdb .
./calc_distances_angles.sh
rm -f *.pdb
cd ../

