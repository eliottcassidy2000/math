#!/bin/sh
cd "$(dirname "$0")"
U=..
python zeros_big.py $U/ar10_classes.txt 10 > zb_ar10.out 2>&1
python zeros_big.py $U/ar10_classes.txt 10 --complete > zb_reg11.out 2>&1
python zeros_big.py random 12 20000 2026 > zb_rand12.out 2>&1
python zeros_big.py $U/ar12_samples.txt 12 > zb_ar12.out 2>&1
python zeros_big.py $U/rand12_samples.txt 12 > zb_rand12_file.out 2>&1
python zeros_big.py random 14 3000 2027 > zb_rand14.out 2>&1
python zeros_big.py $U/ar14_samples.txt 14 > zb_ar14.out 2>&1
python zeros_big.py $U/ar16_samples.txt 16 > zb_ar16.out 2>&1
echo done > zb_done.flag
