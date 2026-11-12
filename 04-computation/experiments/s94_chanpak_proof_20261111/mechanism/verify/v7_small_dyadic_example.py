#!/usr/bin/env python3
"""Verifier: explicit order-7 arc cubes of dimension 3 and 4 covering small dyadic intervals (X < 24), checked by
brute force over permutations.  Shows that 'no 4-cube covers a dyadic interval' needs the qualifier X >= 24."""
import itertools, numpy as np
from vlib import *
n = 7
pi, pj = pair_arrays(n)
H, S = build_table(n, pi, pj)
prs = pairs_of(n); na = len(prs)
def find(k, X, need_strong_count=False):
    want = set(v for v in range(X, 2 * X + 1) if v % 2 == 1)
    for arcs in itertools.combinations(range(na), k):
        allb = sum(1 << a for a in arcs)
        offs = [sum(1 << arcs[t] for t in range(k) if (s >> t) & 1) for s in range(1 << k)]
        free = ((1 << na) - 1) ^ allb
        # vectorised over bases
        bases = np.arange(1 << na, dtype=np.int64)
        bases = bases[(bases & allb) == 0]
        ok = np.ones(len(bases), bool)
        vals = np.stack([H[bases | o] for o in offs], axis=1)
        for w in want:
            ok &= (vals == w).any(axis=1)
        idx = np.where(ok)[0]
        if len(idx):
            m = int(bases[idx[0]])
            return arcs, m, offs
    return None
for k, X in ((3, 8), (3, 9), (3, 10), (4, 10), (5, 22)):
    r = find(k, X)
    arcs, m, offs = r
    members = [m | o for o in offs]
    bf = [brute_H(n, adjm_from_mask(n, mm)) for mm in members]
    st = [int(S[mm]) for mm in members]
    want = [v for v in range(X, 2 * X + 1) if v % 2 == 1]
    print(f"k = {k}, X = {X}: arcs {[prs[a] for a in arcs]}, base mask {m}; member H (brute force) {sorted(bf)}; "
          f"strong flags {st}; contains all odd numbers of [{X},{2*X}] = {want}: {set(want) <= set(bf)}")
