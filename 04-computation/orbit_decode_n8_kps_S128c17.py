#!/usr/bin/env python3
"""orbit_decode_n8_kps_S128c17.py -- kind-pasteur S128 cont.17.
The n=8 fourth entries: recompute the 404 non-gridsym quasi-fixed tilings (two-pass engine),
form the 101 free Klein orbits, type each by carrier SC-ness -> orbitsSC(8) + orbitsNS(8)."""
import sys
from math import comb
from collections import Counter
sys.stdout.reconfigure(line_buffering=True)
exec(open('04-computation/selfline_sequence_atlas_kps_S128c15.py',encoding='utf-8').read().split('# (A) Burnside-affine')[0])
n = 8
tiles, tidx, m = setup(n)
gmap = [tidx[(n - y + 1, n - x + 1)] for (x, y) in tiles]
print("pass 1+2 over 2^%d tilings..." % m)
X = []
for t in range(1 << m):
    tv = [(t >> i) & 1 for i in range(m)]
    kv = [1 - b for b in tv]
    B1 = beats(tv, n, tidx)
    B2 = beats(kv, n, tidx)
    if inv_pair(B1, n) == inv_pair(B2, n) and iso(B1, B2, n):
        if not all(tv[i] == tv[gmap[i]] for i in range(m)):
            X.append(t)
    if t % (1 << 19) == 0:
        print("  ...%d, X so far %d" % (t, len(X)), flush=True)
print("non-gridsym quasi-fixed: %d (expect 404)" % len(X))
seen = set(); oSC = oNS = 0
for t in X:
    if t in seen:
        continue
    tv = [(t >> i) & 1 for i in range(m)]
    gt = sum(tv[gmap[i]] << i for i in range(m))
    kt = t ^ ((1 << m) - 1); gkt = gt ^ ((1 << m) - 1)
    seen |= {t, gt, kt, gkt}
    B1 = beats(tv, n, tidx)
    BT = [[B1[u][v] for u in range(n)] for v in range(n)]
    sc = inv_pair(B1, n) == inv_pair(BT, n) and iso(B1, BT, n)
    if sc:
        oSC += 1
    else:
        oNS += 1
print("n=8: orbits=%d = orbitsSC %d + orbitsNS %d" % (oSC + oNS, oSC, oNS))
print("SEQUENCES: orbitsSC = 2, 0, 9, %d ; orbitsNS = 0, 3, 13, %d" % (oSC, oNS))
print("DONE")
