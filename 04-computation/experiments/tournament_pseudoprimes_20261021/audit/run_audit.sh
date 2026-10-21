#!/bin/sh
# AUDIT reproduction (independent auditor, 2026-10-21).  Run from this directory:
#   cd 04-computation/experiments/tournament_pseudoprimes_20261021/audit && sh run_audit.sh
# Needs: clang, python3 (sympy), gentourng + dreadnaut (nauty).  Total time ~2-3 min (census n = 9 ~45 s, circulants n <= 23 ~3 min).
set -e
clang -O3 -o a2_smod a2_smod.c
clang -O3 -o a6_egf_prime a6_egf_prime.c
clang -O3 -o a11_circ a11_circ.c
clang -O3 -o a12_circ_brute a12_circ_brute.c
clang -O3 -o a14_zero_classes a14_zero_classes.c
clang -O3 -o a15_census a15_census.c
python3 a1_exact_bruteforce.py > a1_exact_bruteforce.out              # exact s, D (r<=200), brute force n<=6, factorisations
./a2_smod period 100 4 > a2_period_100.out                              # Thm 1, odd p<=100, m<=4p
./a2_smod period 400 4 > a2_period_400.out                              # Thm 1, odd p<=400, m<=4p
python3 a3_groupring_check.py > a3_groupring_check.out                  # replay of the group-ring proof, p<=31
./a2_smod direct 4 2000 > a2_direct_4_2000.out                          # direct s(n) mod n, composite n<=2000
./a2_smod direct 2001 2500 > a2_direct_2001_2500.out
python3 a4_korselt_check.py > a4_korselt_check.out                      # criterion vs direct
python3 a4b_session_counts.py > a4b_session_counts.out
./a2_smod single 262143 > a2_single_262143.out                          # 2^18-1 direct (~1 min)
python3 a5_mersenne.py > a5_mersenne.out                                # all composite 2^k-1, k<=64
./a2_smod wief 3000 > a2_wief_3000.out                                  # Camion-Wieferich primes p<=3000
python3 a7_prime_square_identity.py > a7_prime_square_identity.out      # s(p^2)-c(p^2) = eps D(p) mod p^2, p<=97
./a2_smod seq 1401 45 > a2_seq_mod45.txt
python3 a8_mod3_mod5_mod9.py > a8_mod3_mod5_mod9.out                    # statement 3; mod-9 law to 1401
python3 a9_count_test.py > a9_count_test.out                            # statement 6 + Burnside m=3,5,7
python3 a10_controls.py > a10_controls.out                              # statement 7 controls
python3 a10b_connected_extended.py > a10b_connected_extended.out        # connected graphs ARE quasi-periodic (m>=2)
for n in 3 5 7 9 11 13; do ./a11_circ $n all; done > a11_circ_3_13.out
for n in 15 17; do ./a11_circ $n all; done > a11_circ_15_17.out
for n in 19 21 23; do ./a11_circ $n all; done > a11_circ_19_23.out
cat a11_circ_3_13.out a11_circ_15_17.out a11_circ_19_23.out > a11_circ_all.out
(for n in 3 5 7 9; do ./a12_circ_brute $n; done) > a12_circ_brute.out; ./a12_circ_brute 11 > a12_circ_brute_11.out
python3 a12_join.py > a12_join.out                                      # brute force vs subset DP, all S, n <= 11
for n in 2465 2821 4623 6601 8911 10585 12099 15841 16383 19599; do ./a2_smod single $n; done > a2_single_families.out
python3 a18_thm4747_checks.py > a18_thm4747_checks.out                  # THM-4747 congruences, Paley data, HYP-9357 table
python3 a13_list20000.py > a13_list20000.out                            # squarefree list to 20000
python3 a13b_exotic_classify.py > a13b_exotic_classify.out
./a14_zero_classes 400 > a14_zero_classes_400.out
python3 a14_zero_classes_check.py > a14_zero_classes_check.out          # zero classes, Poisson table
python3 a17_mod_p2_quasiperiod.py > a17_mod_p2_quasiperiod.out          # mod p^2 quasi-periodicity (m>=3)
python3 a17b_mod_pk_quasiperiod.py > a17b_mod_pk_quasiperiod.out
python3 a19_rotation_proof.py > a19_rotation_proof.out                  # Cor 2 rotation proof
python3 a20_designs.py > a20_designs.out                                # HYP-9357 designs (dreadnaut)
python3 a21_large_primes_period.py > a21_large_primes_period.out
python3 a22_nine_criterion.py > a22_nine_criterion.out
python3 a23_square_criterion.py > a23_square_criterion.out              # criterion for cube-free odd parts
for k in 3 4 5 6 7 8 9; do gentourng -q $k | ./a15_census $k > a15_census_n$k.txt; done
python3 a16_census_analyze.py > a16_census_analyze.out                  # census spectra
echo AUDIT RUN COMPLETE
