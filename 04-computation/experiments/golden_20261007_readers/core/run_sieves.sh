#!/bin/bash
# sequential queue of exact sieve runs (each < ~10 min, ~1-30 MB)
cd "$(dirname "$0")"
./sieve 36 0 0 0            > out_desc_K36.txt
./sieve 32 2048 0 0         > out_join_K32.txt          # Roosendaal-type: descent + all translation joins
./sieve 30 2048 1 0 0 1     > out_branch_mergeafter_K30.txt   # branch first, then merges on survivors (subsumption test)
./sieve 32 0 1 0 1          > out_pred1_K32.txt         # descent + depth-1 predecessor (Angeltveit-type)
./sieve 30 2048 1 0 1       > out_pred1_join_K30.txt    # depth-1 predecessor + joins
./sieve 34 0 1 0            > out_branch_K34.txt        # maximal 2-adic (descent + full backward trees)
./sieve 28 0 1 2            > out_branch_J2_K28.txt     # with n mod 9
echo DONE > queue_done.txt
