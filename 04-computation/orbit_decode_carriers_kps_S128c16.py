#!/usr/bin/env python3
"""orbit_decode_carriers_kps_S128c16.py -- kind-pasteur S128 cont.16.
DECODE 2,3,22,101: per-carrier tables + orbit-type split orbitsSC/orbitsNS at n=5,6,7.
An orbit {t,gt,kt,gkt} has cls(t)=cls(kt)=K, cls(gt)=K^op: SC-carrier orbit iff K iso K^op."""
import sys
from math import comb
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)
exec(open('04-computation/selfline_sequence_atlas_kps_S128c15.py',encoding='utf-8').read().split('# (A) Burnside-affine')[0])
def gm(n, tiles, tidx):
    return [tidx[(n - y + 1, n - x + 1)] for (x, y) in tiles]
for n in [5, 6, 7]:
    tiles, tidx, m = setup(n)
    gmap = gm(n, tiles, tidx)
    X = []
    for t in range(1 << m):
        tv = [(t >> i) & 1 for i in range(m)]
        kv = [1 - b for b in tv]
        B1 = beats(tv, n, tidx); B2 = beats(kv, n, tidx)
        if inv_pair(B1, n) == inv_pair(B2, n) and iso(B1, B2, n):
            gs = all(tv[i] == tv[gmap[i]] for i in range(m))
            if not gs:
                X.append(t)
    seen = set(); oSC = oNS = 0; mults = Counter()
    for t in X:
        if t in seen: continue
        tv = [(t >> i) & 1 for i in range(m)]
        gt = sum(tv[gmap[i]] << i for i in range(m))
        kt = t ^ ((1 << m) - 1); gkt = gt ^ ((1 << m) - 1)
        seen |= {t, gt, kt, gkt}
        B1 = beats(tv, n, tidx)
        BT = [[B1[u][v] for u in range(n)] for v in range(n)]
        sc = inv_pair(B1, n) == inv_pair(BT, n) and iso(B1, BT, n)
        if sc: oSC += 1
        else: oNS += 1
    print("n=%d: orbits=%d = orbitsSC %d + orbitsNS %d" % (n, oSC + oNS, oSC, oNS))
print("DONE")
