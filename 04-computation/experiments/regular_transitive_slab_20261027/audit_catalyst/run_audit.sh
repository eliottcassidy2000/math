#!/bin/sh
# Reproduce the THM-4820 audit (about one minute on one core).  Run from this directory.
set -e
cc -O3 -march=native -o audit_dp audit_dp.c
cc -O3 -march=native -o audit_lattice audit_lattice.c
python3 audit_algebra.py > audit_algebra.out
python3 audit_small_trees.py > audit_small_trees.out
./audit_dp 46 35 35 verify > audit_dp_46_35_35.out 2> audit_dp_46_35_35.err
./audit_dp 47 36 36 verify tree 47 36 36 > audit_dp_47_36_36.out 2> audit_dp_47_36_36.err; mv audit_tree_out.txt audit_own_tree_47_36_36.txt
./audit_dp 47 36 36 lastmin tree 47 36 36 > /dev/null 2>&1; mv audit_tree_out.txt audit_own_tree_47_36_36_lastmin.txt
./audit_dp 76 62 62 > audit_dp_76_62_62.out 2> audit_dp_76_62_62.err
python3 audit_tree.py ../catalyst/tree_47_36_36.txt 47,36,36 > audit_tree.out
python3 audit_tree.py audit_own_tree_47_36_36.txt 47,36,36 >> audit_tree.out
python3 audit_tree.py audit_own_tree_47_36_36_lastmin.txt 47,36,36 >> audit_tree.out
./audit_lattice 2500 2600 6 dump > audit_lat_W2500.txt 2> audit_lat_W2500.log
./audit_lattice 2500 0 3 dump > audit_lat_W2500_m3.txt 2> audit_lat_W2500_m3.log
./audit_lattice 400 0 6 > audit_lat_W400.txt 2> audit_lat_W400.log
./audit_lattice 900 0 6 > audit_lat_W900.txt 2> audit_lat_W900.log
./audit_lattice 3500 0 6 dump > audit_lat_W3500.txt 2> audit_lat_W3500.log
./audit_lattice 7500 7500 6 dump > audit_lat_W7500.txt 2> audit_lat_W7500.log
./audit_lattice 10500 10500 6 dump tree 101 0 > audit_lat_W10500.txt 2> audit_lat_W10500.log; mv audit_value_tree.txt audit_value_tree_101.txt
python3 audit_compare.py > audit_compare.out
python3 audit_value_tree.py audit_value_tree_101.txt 101 0 > audit_value_tree_101.out
python3 audit_excess2.py 10500 101 0 > audit_excess2_W10500.out
python3 audit_slack_levels.py 3000 > audit_slack_levels.out
echo "audit outputs regenerated"
