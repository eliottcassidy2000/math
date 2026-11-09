#!/bin/sh
# S19 audit: reproduce every output in out/. BUILD must be a directory OUTSIDE the repository.
set -e
BUILD=${BUILD:-/tmp/s19audit}; mkdir -p "$BUILD"
gcc -O2 -o "$BUILD/tbrute" tbrute.c
gcc -O2 -o "$BUILD/cores" cores.c
gcc -O2 -o "$BUILD/cover" cover.c
gcc -O2 -o "$BUILD/cover2" cover2.c
gcc -O2 -o "$BUILD/q4" q4.c
gcc -O2 -o "$BUILD/verify_cert" verify_cert.c
gcc -O2 -o "$BUILD/find_core" find_core.c
python3 naive_check.py "$BUILD/tbrute" 2000 11 > out/naive_check.out
python3 check_transfer.py "$BUILD/tbrute" 2026 600 > out/check_transfer.out
python3 check_transfer.py "$BUILD/tbrute" 99 1500 > out/check_transfer_seed99.out
"$BUILD/cores" classes 7 > "$BUILD/class_stats_le7.txt"
for n in 1 2 3 4 5 6 7; do "$BUILD/cores" labelled $n > "$BUILD/labelled_$n.txt" 2>> "$BUILD/labelled.log"; done
grep '^order 7' "$BUILD/labelled.log" | tail -1 > out/labelled_order7.out
python3 dedup_cores.py "$BUILD/class_stats_le7.txt" "$BUILD" "$BUILD/labelled_5.txt" "$BUILD/labelled_6.txt" "$BUILD/labelled_7.txt" > out/core_states.out
grep '^#' "$BUILD/class_stats_le7.txt" > out/class_counts.out
for f in cores_b2_le6 cores_b2_le7 cores_b3_le6 cores_b3_le7; do cp "$BUILD/$f.txt" out/; done
"$BUILD/cover" "$BUILD/cores_b2_le5.txt" 2 20000 > out/cover_b2_le5_2e4.out
"$BUILD/cover" "$BUILD/cores_b3_le5.txt" 3 20000 > out/cover_b3_le5_2e4.out
"$BUILD/cover" "$BUILD/cores_b3_le6.txt" 3 20000 > out/cover_b3_le6_2e4.out
"$BUILD/cover" "$BUILD/cores_b2_le6.txt" 2 200000 > out/cover_b2_le6_2e5.out
"$BUILD/cover" "$BUILD/cores_b2_le6.txt" 2 1000000 > out/cover_b2_le6_1e6.out
"$BUILD/cover" "$BUILD/cores_b2_le7.txt" 2 1000000 "$BUILD/hits_b2_le7_1e6.bin" > out/cover_b2_le7_1e6.out
"$BUILD/cover" "$BUILD/cores_b3_le7.txt" 3 1000000 "$BUILD/hits_b3_le7_1e6.bin" > out/cover_b3_le7_1e6.out
"$BUILD/cover2" "$BUILD/cores_b2_le7.txt" 2 1000000 "$BUILD/hits_b2_le7_1e6.bin" > out/cover2_b2_le7_1e6.out
"$BUILD/cover2" "$BUILD/cores_b3_le7.txt" 3 1000000 "$BUILD/hits_b3_le7_1e6.bin" > out/cover2_b3_le7_1e6.out
"$BUILD/cover2" "$BUILD/cores_b2_le6.txt" 2 1000000 > out/cover2_b2_le6_1e6.out
"$BUILD/cover2" "$BUILD/cores_b3_le6.txt" 3 1000000 > out/cover2_b3_le6_1e6.out
"$BUILD/cover" "$BUILD/cores_b2_le7.txt" 2 10000000 "$BUILD/hits_b2_le7_1e7.bin" > out/cover_b2_le7_1e7.out
"$BUILD/cover" "$BUILD/cores_b3_le7.txt" 3 10000000 "$BUILD/hits_b3_le7_1e7.bin" > out/cover_b3_le7_1e7.out
"$BUILD/cover2" "$BUILD/cores_b2_le7.txt" 2 10000000 "$BUILD/hits_b2_le7_1e7.bin" > out/cover2_b2_le7_1e7.out
"$BUILD/cover2" "$BUILD/cores_b3_le7.txt" 3 10000000 "$BUILD/hits_b3_le7_1e7.bin" > out/cover2_b3_le7_1e7.out
"$BUILD/cover2" "$BUILD/cores_b2_le6.txt" 2 10000000 > out/cover2_b2_le6_1e7.out
python3 density_windows.py out/cover2_b2_le6_1e7.out > out/b2_le6_density_1e7.out
"$BUILD/q4" "$BUILD/cores_b2_le7.txt" 2 1001 200000 ../out/q4_greedy_b2.out "$BUILD/q4_sparse_b2.txt" 80 > out/q4_b2.out
"$BUILD/q4" "$BUILD/cores_b3_le7.txt" 3 1000 200000 ../out/q4_greedy_b3.out "$BUILD/q4_sparse_b3.txt" 80 > out/q4_b3.out
(python3 q4_lp.py "$BUILD/q4_sparse_b2.txt" 45 0 out/q4_lp_b2_cert_K45.txt; python3 q4_lp.py "$BUILD/q4_sparse_b2.txt" 80) > out/q4_lp_b2.out
(python3 q4_lp.py "$BUILD/q4_sparse_b3.txt" 45 0 out/q4_lp_b3_cert_K45.txt; python3 q4_lp.py "$BUILD/q4_sparse_b3.txt" 80) > out/q4_lp_b3.out
# optional, slow: python3 q4_lp.py "$BUILD/q4_sparse_b2.txt" 45 120 > out/q4_milp_b2_K45.out
(python3 cert_targets_check.py out/q4_lp_b2_cert_K45.txt 2 1001; python3 cert_targets_check.py out/q4_lp_b3_cert_K45.txt 3 1000) > out/q4_cert_targets_check.out
"$BUILD/verify_cert" "$BUILD/cores_b2_le7.txt" 2 200000 out/q4_lp_b2_cert_K45.txt > out/q4_cert_verify_b2.out
"$BUILD/verify_cert" "$BUILD/cores_b3_le7.txt" 3 200000 out/q4_lp_b3_cert_K45.txt > out/q4_cert_verify_b3.out
python3 oneapex.py 22 > out/oneapex_share_B22.out
python3 capacity.py "$BUILD/cores_b2_le5.txt" "$BUILD/cores_b2_le6.txt" "$BUILD/cores_b2_le7.txt" > out/capacity_b2.out
python3 dyadic_rays.py "$BUILD/cores_b2_le7.txt" > out/dyadic_rays.out
python3 dyadic_bridge_example.py "$BUILD/tbrute" > out/dyadic_ray_bridges_bruteforce.out
"$BUILD/find_core" 7 71 287 -5 31    # prints the core used by dyadic_bridge_example.py
python3 dictionary_check.py > out/dictionary_check.out
