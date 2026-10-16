#!/usr/bin/env python3
"""Integrity checks for a single (unreplicated) run of the dense resumable census (S38 audit E6).

The per-block invariants of collatz_census_dense_resumable_20261016.py do not cover the weighted moments, the energy sums,
the parity table or the grid pass.  This script checks the published outputs of one run against bounds that any
corrupted value of significant size would violate:
  * grid values: |S(N') - S(N)| <= N' - N and |H(N') - H(N)| <= 2 (N' - N) for consecutive grid points, including the
    step from the resume point (S, H at 2^K0 from the previous run's JSON);
  * prefix sums: sum_{n<=N} H(n)^2 non-decreasing, and |Delta sum H| <= Delta n * (max |H| allowed by Lipschitz);
  * moments: per sub-block |Sigma t| <= 2^17 and |Sigma (n - lo) t| <= 2^17 (2^17 - 1)/2, and the total Sigma t over the run
    plus the start value equals S at the end (minus t(2^K) as stored).
It cannot detect small (low-bit) corruptions of the moments; the replication run does that.
Session opus-2026-10-16-S38.  Usage: python collatz_census_single_run_integrity_20261016.py run.npz run.json run_moments.npz
       previous_run.json
"""
import sys, json
import numpy as np

NPZ, JS, MOM, PREV = sys.argv[1:5]
Z = np.load(NPZ)
g = np.array(Z['grid'], dtype=np.int64)
S = np.array(Z['S'], dtype=np.int64)
H = np.array(Z['H'], dtype=np.int64)
p1 = [int(x) for x in Z['preH1_str']]
p2 = [int(x) for x in Z['preH2_str']]
J = json.load(open(JS))
P = json.load(open(PREV))
k0 = int(J['k0'])
N0 = 1 << k0
S0, H0 = int(P['S']), int(P['H'])
out = {}
dn = np.diff(np.concatenate([[N0], g]))
dS = np.diff(np.concatenate([[S0], S]))
dH = np.diff(np.concatenate([[H0], H]))
out['lipschitz_S_ok'] = bool(np.all(np.abs(dS) <= dn))
out['lipschitz_H_ok'] = bool(np.all(np.abs(dH) <= 2 * dn))
out['preH2_nondecreasing'] = all(p2[i + 1] >= p2[i] for i in range(len(p2) - 1))
# |Delta preH1| <= Delta n * (|H(N)| + 2 Delta n)
ok1 = True
for i in range(len(g) - 1):
    d = int(g[i + 1]) - int(g[i])
    bound = d * (abs(int(H[i])) + 2 * d)
    if abs(p1[i + 1] - p1[i]) > bound:
        ok1 = False
        break
out['preH1_increments_bounded'] = ok1
M = np.load(MOM)
mom = np.array(M['moments'], dtype=np.int64)
SUB = int(M['B'])
out['moment0_bounded'] = bool(np.all(np.abs(mom[:, 0]) <= SUB))
out['moment1_bounded'] = bool(np.all(np.abs(mom[:, 1]) <= SUB * (SUB - 1) // 2))
# total t over the run: S(2^K) - S(2^K0) = sum mom0 + t(2^K) - t(2^K0)
K = int(J['logmax'])
tK = -1 if K % 2 else 1
tK0 = -1 if k0 % 2 else 1
tot = int(mom[:, 0].sum())
out['moment_total_matches_S'] = bool(int(J['S']) - S0 == tot + tK - tK0)
out['n_grid'] = int(len(g))
out['n_subblocks'] = int(len(mom))
out['all_ok'] = all(v for k, v in out.items() if isinstance(v, bool))
print(json.dumps(out, indent=1))
if not out['all_ok']:
    raise SystemExit('integrity check failed')
