#!/bin/sh
# Reproduce every output of the microcosm packet (mac-mini-2026-11-06-microcosm).
set -e
cd "$(dirname "$0")"
python3 gate_basics.py > gate_basics.out 2>&1
python3 gate_descent_face_formula.py > gate_descent_face_formula.out 2>&1
python3 gate_resolvent_exact.py > gate_resolvent_exact.out 2>&1
python3 tower_limits.py > tower_limits.out 2>&1
python3 tower_spectrum.py > tower_spectrum.out 2>&1
cc -O2 -o tower_k1_paths tower_k1_paths.c && ./tower_k1_paths > tower_k1_paths.out 2>&1 && rm -f tower_k1_paths
cc -O2 -o gated_spectrum gated_spectrum.c && ./gated_spectrum > gated_spectrum.out 2>&1 && rm -f gated_spectrum
python3 tower_ocf.py > tower_ocf.out 2>&1
cc -O2 -o hc_gate_free hc_gate_free.c && for n in 5 6 7 8; do "${GENTOURNG:-/opt/homebrew/bin/gentourng}" -c -q $n | ./hc_gate_free $n; done > hc_gate_free.out 2>&1 && rm -f hc_gate_free
python3 c3_first_order_law.py > c3_first_order_law.out 2>&1
python3 cover_tower_values.py > cover_tower_values.out 2>&1
for n in 5 6 7 8 9; do python3 phi_strong_minimum.py $n; done > phi_strong_minimum.out 2>&1
echo "all done"
