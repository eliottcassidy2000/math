#!/usr/bin/env python3
"""Turn audit F's numerical rank >= 3 evidence (HYP-9244, 04-computation/experiments/reframes_20261007/audit_F/partC_maps.out)
into proofs: adaptive block certificates (block_adaptive.adaptive_certify) for the maps with their exact constants r_i.
Usage: python3 auditF_maps_cert.py NAME K0 KMAX"""
import math, sys, time
from block_balance import int_coords, ratio_group
from block_adaptive import adaptive_certify
MAPS = {
    'Z5_1_8_3_7_12': (5, [1, 8, 3, 7, 12], [0, -3, 4, 4, 2]),
    'Z5_1_2_3_7_6': (5, [1, 2, 3, 7, 6], [0, 3, 4, 4, 1]),
    'Z7_1_2_3_5_1_1_1': (7, [1, 2, 3, 5, 1, 1, 1], [0, 5, 1, 6, 3, 2, 1]),
}
name, K0, KMAX = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
d, m, r = MAPS[name]
assert all((m[i] * i + r[i]) % d == 0 for i in range(d))
rho, v = int_coords(m); G = ratio_group(m, d)
Lam = sum(math.log(x / d) for x in m) / d
bad = [(a, b) for a in G for b in range(d) if (a, b) != (1, 0) and all(m[(a * j + b) % d] == m[j] for j in range(d))]
t0 = time.time()
print(f"{name}: d = {d}, m = {m}, r = {r}, rank {rho}, Lambda = {Lam:+.4f}, coupling group {G}, property (i) {'holds' if not bad else 'FAILS ' + str(bad)}", flush=True)
res = adaptive_certify(d, m, r, K0, KMAX, 0.02, log=lambda x: print(x, flush=True))
if res:
    Qi, nU, counts, mgl = res
    print(f"   EXACT: integer form Q = {Qi} balances all {nU} distinct leaf blocks (block lengths {K0}..{KMAX}; nodes refined per level {counts}; leaf margin {mgl:+.4f})  [{time.time() - t0:.0f}s]", flush=True)
else:
    print(f"   not certified with lengths {K0}..{KMAX}  [{time.time() - t0:.0f}s]", flush=True)
