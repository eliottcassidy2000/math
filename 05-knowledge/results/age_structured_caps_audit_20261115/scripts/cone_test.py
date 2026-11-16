#!/usr/bin/env python3
"""cone_test.py: is the history-refined feature-cap cone an invariant certificate for literal span-3 growth?
usage: cone_test.py mats3.txt frag.txt K rounds [lo:hi] caps_file [caps_file ...]
Caps files are hratio5 outputs ("k key feature num den"); the cap of (k, key, feature) is the maximum over the files.
For history length K the cone of a key (the last K letters) is
    C(key) = { x >= 0 : f_j . x <= cap[key][j] * h . x for every feature j },
with f_j the coordinates (j < N), the fragment-class sums (N..N+2), the two-fragment functionals Z_l (N+3..N+10) and
the next-step counts H(W l) = (h A_l) . x (N+11..N+18).
Soundness: if every word of length K lies in its cone (true by construction: caps are maxima over data that include
those words) and A_l C(key) is contained in C(key l) for all key, l, then every reachable state lies in its cone, so the
caps are theorems; in particular the lemma Z_l <= 38 H follows when every Z_l cap is <= 38.
Each containment is checked by linear programmes: for every (key, l, feature j) maximise
    (f_j A_l) . x - cap[key l][j] * (h A_l) . x   subject to  x in C(key), h . x = 1.
A positive maximum is a violation; the cap of the target is then raised to the maximum (value iteration) and the
round repeated, up to `rounds` rounds.  Prints per round: violations, largest excess, the largest Z cap."""
import sys, numpy as np
from fractions import Fraction
from collections import defaultdict
from scipy.optimize import linprog
mats, fragf, K, rounds = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
# optional: restrict the source keys to a slice KEYLO:KEYHI (for parallel single-round runs), given as 'lo:hi' first
keyslice = None
if ':' in sys.argv[5]:
    lo, hi = map(int, sys.argv[5].split(':')); keyslice = (lo, hi); sys.argv.pop(5)
tok = open(mats).read().split(); N = int(tok[0]); p = 1
A = []
for b in range(8):
    A.append(np.array(tok[p:p + N * N], dtype=float).reshape(N, N)); p += N * N
v0 = np.array(tok[p:p + N], dtype=float); p += N; h = np.array(tok[p:p + N], dtype=float)
nfrag = list(map(int, open(fragf).read().split()))
hA = [h @ A[b] for b in range(8)]
NF = N + 3 + 8 + 8
F = np.zeros((NF, N))
for s in range(N): F[s, s] = 1
for s in range(N):
    if 2 <= nfrag[s] <= 4: F[N + nfrag[s] - 2, s] = 1
for b in range(8):
    for s in range(N):
        if nfrag[s] == 2: F[N + 3 + b, s] = hA[b][s]
        F[N + 11 + b, s] = hA[b][s]
caps = defaultdict(lambda: np.zeros(NF))
for fn in sys.argv[5:]:
    for line in open(fn):
        if line.startswith('#'): continue
        k, key, j, num, den = map(int, line.split())
        if k == K:
            caps[key][j] = max(caps[key][j], num / den)
keys = sorted(caps)
if keyslice: keys = [k for k in keys if keyslice[0] <= k < keyslice[1]]
print(f"K = {K}: {len(caps)} history keys with data ({len(keys)} tested in this run); features {NF}")
def shift(key, b): return ((key * 8) + b) % (8 ** K) if K > 0 else 0
tol = 1e-7
for rnd in range(1, rounds + 1):
    viol = 0; worst = 0.0; raised = defaultdict(dict); checked = 0
    for key in keys:
        c = caps[key]
        # cone constraints: F x <= c (since h.x = 1); coordinates as bounds, group rows as A_ub
        bounds = [(0, c[s]) for s in range(N)]
        A_ub = F[N:]; b_ub = c[N:]
        for b in range(8):
            key2 = shift(key, b)
            if key2 not in caps: continue
            c2 = caps[key2]
            for j in range(NF):
                obj = -(F[j] @ A[b] - c2[j] * hA[b])        # minimise the negative
                if not np.any(obj): continue
                checked += 1
                res = linprog(obj, A_ub=A_ub, b_ub=b_ub, A_eq=h.reshape(1, -1), b_eq=[1.0], bounds=bounds, method='highs')
                if res.status != 0: continue
                excess = -res.fun
                if excess > tol:
                    viol += 1; worst = max(worst, excess)
                    # the new cap: max of (f_j A_b x)/(hA_b x) at the optimum
                    x = res.x; denom = hA[b] @ x
                    if denom > 1e-12:
                        newcap = (F[j] @ A[b] @ x) / denom
                        raised[key2][j] = max(raised[key2].get(j, 0.0), newcap)
    zmax = max(caps[k][N + 3:N + 11].max() for k in keys)
    print(f"round {rnd}: LPs {checked}, violations {viol}, largest excess {worst:.6g}, largest Z cap {zmax:.6g}, largest coordinate cap {max(caps[k][:N].max() for k in keys):.6g}")
    if viol == 0:
        print("INVARIANT: the cone is a certificate; Z caps per key <= 38:", all(caps[k][N + 3:N + 11].max() <= 38 + 1e-9 for k in keys)); break
    for key2, d in raised.items():
        for j, v in d.items(): caps[key2][j] = max(caps[key2][j], v)
