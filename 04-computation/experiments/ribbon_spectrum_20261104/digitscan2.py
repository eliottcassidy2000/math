#!/usr/bin/env python3
"""digitscan2.py -- extended digit-block search (numpy): blocks B of length <= 5 over the 3-letter
alphabet and over all 4 letters (length <= 4); multipliers alpha in 2..9; residual H(uB) - alpha H(u)
over all prefix vectors u of length 3..9 (all 4 letters).  A digit block would have a bounded residual
set; we report the smallest residual-set sizes."""
import itertools, numpy as np
from ribbon import transfer_system, letters, word_str
w = 2
order, states, mats, final = transfer_system(w)
d = len(order); L = letters(w)
F = frozenset
M = {}
for a in L:
    A = np.zeros((d, d), dtype=np.int64)
    for (j, i), m in mats[a].items(): A[j, i] = m
    M[a] = A
phi = np.zeros(d, dtype=np.int64); phi[final] = 1
# prefix vectors
cur = [np.eye(d, dtype=np.int64)[0]]
prefixes = []
for n in range(1, 10):
    cur = [M[a] @ v for v in cur for a in L]
    # dedupe
    uniq = {tuple(v) for v in cur}
    cur = [np.array(t, dtype=np.int64) for t in uniq]
    if n >= 3: prefixes.extend(cur)
P = np.array(prefixes).T   # d x K
K = P.shape[1]
h = phi @ P
print(f"{K} prefix vectors")
results = []
alph3 = [F([1]), F([2]), F([1, 2])]
blocks = [B for k in range(1, 6) for B in itertools.product(alph3, repeat=k)] + [B for k in range(1, 5) for B in itertools.product(L, repeat=k) if any(not a for a in B)]
for B in blocks:
    MB = np.eye(d, dtype=np.int64)
    for a in B: MB = M[a] @ MB
    hB = phi @ (MB @ P)
    for alpha in range(2, 10):
        res = hB - alpha * h
        nres = len(set(res.tolist()))
        results.append((nres, alpha, word_str(list(B), w), sorted(set(res.tolist()))[:6]))
results.sort()
for r in results[:15]:
    print(f"  {r[0]:4d} distinct residuals: block {r[2]:>14} alpha={r[1]} (smallest residuals {r[3]})")
