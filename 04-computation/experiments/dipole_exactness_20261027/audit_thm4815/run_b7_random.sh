#!/bin/sh
# random battery for B7 (audit); seeds fixed for reproducibility
for mode in 0 1 2 3 4; do
  ./audit_b7_random 11 100000 $mode 11
  ./audit_b7_random 12 100000 $mode 12
  ./audit_b7_random 13 40000 $mode 13
  ./audit_b7_random 14 40000 $mode 14
  ./audit_b7_random 15 10000 $mode 15
  ./audit_b7_random 16 10000 $mode 16
  ./audit_b7_random 18 500 $mode 18
  ./audit_b7_random 20 100 $mode 20
done
