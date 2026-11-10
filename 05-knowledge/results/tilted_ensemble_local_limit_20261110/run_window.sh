#!/bin/sh
# Two-dimensional deletion-pair census (section 6).  Needs gentourng and cc.  About 15 minutes on 8 cores.
set -e
W=${TMPDIR:-/tmp}/tilted_window_$$; mkdir -p "$W"
cc -O2 -o "$W/pairs" pairs.c
for n in 2 3 4 5 6 7 8 9; do gentourng -q $n 2>/dev/null | "$W/pairs" $n "$W/pA9.bin" "$W/pZ9.bin"; done
for r in 0 1 2 3 4 5 6 7; do ( gentourng -q 10 $r/8 2>/dev/null | "$W/pairs" 10 "$W/pA10_$r.bin" "$W/pZ10_$r.bin" ) & done
wait
python3 merge_parts.py bits "$W/pA_le10.bin" "$W/pA9.bin" "$W"/pA10_[0-7].bin
python3 merge_parts.py bits "$W/pZ_le10.bin" "$W/pZ9.bin" "$W"/pZ10_[0-7].bin
python3 analyze_window.py "$W/pA_le10.bin" "$W/pZ_le10.bin"
