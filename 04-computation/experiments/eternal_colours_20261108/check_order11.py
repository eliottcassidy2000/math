#!/usr/bin/env python3
# mac-mini-2026-11-08-firstrealizers: consistency check of sc_odd11.out (all self-converse 11-tournaments with H < 625).
# Padding a self-converse T of order n by a new source and sink (source => T => sink) is self-converse with the same H,
# so every self-converse value of order 9 (and, padded twice, of order 7, 5, 3) below 625 must occur at order 11.
# Also the order-11 values must all be odd and below 625, and each strong flag must be compatible (H >= f(11) = 225).
import sys, os
D = sys.argv[1] if len(sys.argv) > 1 else "."
v11 = {}
for line in open(os.path.join(D, "sc_odd11.out")):
    _, h, st = line.split(); v11[int(h)] = int(st)
ok = True
for n in (3, 5, 7, 9):
    sc = set()
    for line in open(os.path.join(D, f"cls_{n}.txt")):
        hd = line.split("|")[0].split()
        if int(hd[4]) and int(hd[1]) < 625: sc.add(int(hd[1]))
    miss = sorted(sc - set(v11))
    print(f"order-{n} self-converse values < 625: {len(sc)}; missing at order 11: {miss}")
    ok &= not miss
bad = [h for h in v11 if h % 2 == 0 or h >= 625]
weak = [h for h, s in v11.items() if s and h < 225]
print(f"order-11 values: {len(v11)}; strong ones: {sum(v11.values())}; even or >= 625: {bad}; strong below f(11) = 225: {weak}")
ok &= not bad and not weak
print("CHECK", "PASSED" if ok else "FAILED")
