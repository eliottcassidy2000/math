#!/bin/bash
# Parallel driver for the S54 prediction tables (opus-2026-10-22-S54).  Splits the class/sample files into chunks,
# runs predict_nonregular.py on each chunk and diag_exact.py on the diagnostic subsets, then merges the per-chunk
# prediction files.  Usage: bash run_all.sh   (from this directory; ~40 min on 12 cores)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p chunks
head -1500 ar12_samples.txt > chunks/ar12_1500.txt
head -1500 ar14_samples.txt > chunks/ar14_1500.txt
head -600 ar16_samples.txt > chunks/ar16_600.txt
head -1000 ar12_samples.txt > ar12_sub.txt
head -1000 ar14_samples.txt > ar14_sub.txt
head -300 ar16_samples.txt > ar16_sub.txt
split -n l/6 -d ar10_classes.txt chunks/c10_
split -n l/2 -d chunks/ar12_1500.txt chunks/c12_
split -n l/2 -d chunks/ar14_1500.txt chunks/c14_
cp chunks/ar16_600.txt chunks/c16_00
cp ar20_samples.txt chunks/c20_00
for f in chunks/c10_* chunks/c12_* chunks/c14_* chunks/c16_* chunks/c20_*; do
    case "$f" in *.pred|*.log) continue;; esac
    n=$(basename "$f" | sed 's/^c\([0-9]*\)_.*/\1/')
    python predict_nonregular.py "$n" "$f" "$f.pred" > "$f.log" 2>&1 &
done
for n in 12 14 16; do
    python diag_exact.py "$n" "ar${n}_sub.txt" > "diag${n}.out" 2>&1 &
done
wait
for n in 10 12 14 16 20; do
    first=1
    for p in chunks/c${n}_*.pred; do
        if [ $first = 1 ]; then cat "$p"; first=0; else tail -n +2 "$p"; fi
    done > "pred${n}.txt"
done
echo ALLDONE
