#!/usr/bin/env python3
"""
The [14,48] x duty-structure residue check, RUN at bounded height (the finite
check in action), + the exact H(41) impersonation floor.

D1: for all 78 two-far complements, scan every pair 14 <= x < y <= B and test
    the pure residue condition: (x mod q, y mod q) in E_q for ALL q in [14,48].
    Zero survivors ==> certified: no two-far family with fars <= B escapes the
    strip+window residue system (independent of, and beyond, the exact sweeps'
    exit-q observations).
D2: exact x-floor of the impersonation channel for {1,2} at ceiling Q=41
    (x mod q in {+-1,+-2} for all q in [14,41]) -- the H(41) of the ladder law,
    to compare with kps's constructed witness x = 1 + L41 ~ 2.19e17.
"""
import time
from math import gcd
from itertools import combinations
from functools import reduce
import numpy as np

T0 = time.time()

def tq(q): return -(-3*q // 41)
def dist(v, j, q):
    r = (v*j) % q
    return min(r, q-r)

def escape_bool(S, q):
    """q x q boolean array: esc[u,v] True iff (u,v) escapes q for S u {x,y}."""
    t = tq(q)
    A = [j for j in range(1, q) if all(dist(s, j, q) >= t for s in S)]
    if not A:
        return np.ones((q, q), dtype=bool)
    u = np.arange(q)
    clear = np.zeros((len(A), q), dtype=bool)   # clear[ji, u]: dist(u*j) >= t
    for i, j in enumerate(A):
        r = (u*j) % q
        clear[i] = np.minimum(r, q-r) >= t
    esc = np.ones((q, q), dtype=bool)
    for i in range(len(A)):
        cu = clear[i]
        esc &= ~np.outer(cu, cu)   # j witnesses (u,v) if both clear
    return esc

B = 20000
QS = list(range(14, 49))
print(f"D1: residue check over q in [14,48], all 78 complements, pairs up to B={B}")
total_surv = 0
yidx = {q: (np.arange(B+1) % q) for q in QS}
for ab in combinations(range(1, 14), 2):
    S = [v for v in range(1, 14) if v not in ab]
    tables = {q: escape_bool(S, q) for q in QS}
    # order moduli by escape density ascending (strongest filter first)
    order = sorted(QS, key=lambda q: tables[q].mean())
    surv = []
    for x in range(14, B):
        ok = None
        for q in order[:6]:
            row = tables[q][x % q]
            if not row.any():
                ok = False; break
        if ok is False:
            continue
        mask = np.ones(B+1, dtype=bool)
        mask[:x+1] = False
        for q in order:
            mask &= tables[q][x % q][yidx[q]]
            if not mask.any():
                break
        else:
            ys = np.nonzero(mask)[0]
            for y in ys:
                surv.append((x, int(y)))
    if surv:
        print(f"  removed {ab}: SURVIVORS {surv[:5]}{'...' if len(surv)>5 else ''}")
        total_surv += len(surv)
print(f"  => {'ZERO pairs' if total_surv==0 else str(total_surv)+' pairs'} with "
      f"14 <= x < y <= {B} escape the full [14,48] residue system, over all 78 "
      f"complements   [{time.time()-T0:.1f}s]")

# ---------------------------------------------------------------------------
print(f"\nD2: exact impersonation x-floor at ceiling Q=41 (the H(41) datum)")
cluster = [32, 27, 25, 7, 11, 13, 17, 19, 23]
iso41   = [29, 31, 37, 41]
Mc = reduce(lambda a, b: a*b, cluster)
Mi = reduce(lambda a, b: a*b, iso41)
comp_qs = [q for q in range(14, 42) if q not in iso41]
combos = [(0, 1)]
for m in cluster:
    R = sorted({1 % m, (m-1) % m, 2 % m, (m-2) % m})
    new = []
    for (r, M) in combos:
        inv = pow(M % m, -1, m)
        for c in R:
            new.append((r + M*((c - r) * inv % m), M*m))
    combos = new
surv = [r for (r, M) in combos if all(r % q in (1, 2, q-1, q-2) for q in comp_qs)]
iso_res = [(0, 1)]
for p in iso41:
    new = []
    for (r, M) in iso_res:
        inv = pow(M % p, -1, p)
        for c in (1, 2, p-1, p-2):
            new.append((r + M*((c - r) * inv % p), M*p))
    iso_res = new
invMc = pow(Mc % Mi, -1, Mi)
best = None
for r in surv:
    for (s, _) in iso_res:
        x = r + Mc * (((s - r) * invMc) % Mi)
        if x >= 14 and (best is None or x < best):
            best = x
L41 = Mc * Mi
assert all(best % q in (1, 2, q-1, q-2) for q in range(14, 42))
print(f"  coherent cluster residues: {len(surv)}; channel classes mod L41: "
      f"{len(surv)*len(iso_res):,}")
print(f"  EXACT H(41) x-floor = {best} ~ {float(best):.3e}")
print(f"  (kps c96 witness sits at x = 1+L41 = {1+L41} ~ {float(1+L41):.3e}; "
      f"the channel floor is {float(1+L41)/float(best):.0f}x below it)")
print(f"\nTOTAL {time.time()-T0:.1f}s")
