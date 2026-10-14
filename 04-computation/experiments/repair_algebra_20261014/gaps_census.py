#!/usr/bin/env python3
"""Distribution of ladder gaps (2i), length differences and representations per primitive K among canonical increasing representations with exponents <= E (no leading exponent 1).  Usage: python3 gaps_census.py E"""
import numpy as np, sys, time
from collections import Counter
E = int(sys.argv[1])
K = np.array([-1], dtype=np.int64); L = np.array([0], dtype=np.int8); M = np.array([0], dtype=np.int8); F = np.array([0], dtype=np.int8)  # F = first exponent
for e in range(1, E + 1):
    K = np.concatenate([K, 3 * K + (1 << e)]); L = np.concatenate([L, np.full(len(L), e, dtype=np.int8)])
    F = np.concatenate([F, np.where(M == 0, e, F).astype(np.int8)]); M = np.concatenate([M, M + 1])
# canonical reps: strip the root twin (first exponent 1): drop sequences with F == 1 (they duplicate the rest)
keep = F != 1
K, L, M = K[keep], L[keep], M[keep]
o = np.lexsort((L, K)); Ks, Ls, Ms = K[o], L[o], M[o]
bnd = np.nonzero(np.concatenate([[True], Ks[1:] != Ks[:-1], [True]]))[0]
gap = Counter(); lend = Counter(); nrep = Counter(); nprim = 0
for a, b in zip(bnd[:-1], bnd[1:]):
    if b - a < 2: continue
    ls = Ls[a:b].astype(int); ms = Ms[a:b].astype(int)
    if len(set(ls.tolist())) < 2: continue
    nprim += 1; nrep[b - a] += 1
    # pairs with distinct last exponents
    for i in range(b - a):
        for j in range(i + 1, b - a):
            if ls[i] != ls[j]:
                gap[abs(int(ls[j]) - int(ls[i]))] += 1; lend[abs(int(ms[j]) - int(ms[i]))] += 1
print(f"E <= {E}: {nprim} primitive K (canonical reps); reps per primitive K: {dict(sorted(nrep.items()))}")
print("  last-exponent gaps (2i) over pairs with distinct last:", dict(sorted(gap.items())))
print("  length differences over those pairs:", dict(sorted(lend.items())))
