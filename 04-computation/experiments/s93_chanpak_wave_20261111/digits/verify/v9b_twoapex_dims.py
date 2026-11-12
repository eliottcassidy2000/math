#!/usr/bin/env python3
"""Verifier: (a1, a2, dim C) of the category (iii) pairs among the words in the four mu = 0 letters (two-apex family)."""
import json, os
from collections import Counter
from fractions import Fraction
import numpy as np
from v9_twomark_classes import classify, OPS
from vlin import charpoly_int
from vroots import int_roots
here = os.path.dirname(os.path.abspath(__file__))
rep = json.load(open(here + "/../twomark_rep_5_3.json"))
d = rep["d"]
A = [[[int(Fraction(x)) for x in row] for row in M] for M in rep["mats"]]
mm = lambda X, Y: [[sum(X[i][k] * Y[k][j] for k in range(d)) for j in range(d)] for i in range(d)]
stay = [o for o in range(12) if OPS[o][1] == 0]
names = [(o,) for o in stay] + [(o1, o2) for o1 in stay for o2 in stay]
mats = [A[w[0]] if len(w) == 1 else mm(A[w[1]], A[w[0]]) for w in names]
eigs = [[a for a in int_roots(charpoly_int(M)) if a >= 1] for M in mats]
print("eigenvalues >= 1 of the four letters:", [eigs[i] for i in range(4)], "; all integer eigenvalues of the letters:", [int_roots(charpoly_int(mats[i])) for i in range(4)])
print("two-letter words: distinct eigenvalue lists:", Counter(tuple(e) for e in eigs[4:]))
res = {}
for p in (2147483647, 2147483629):
    Mp = [np.array([[x % p for x in r] for r in M], dtype=np.int64) for M in mats]
    inc = [(i, a) for i in range(len(names)) for a in eigs[i]]
    lam = {(i, a): [mats[i][0][k] - (a if k == 0 else 0) for k in range(d)] for (i, a) in inc}
    out = Counter()
    for x in range(len(inc)):
        i, a1 = inc[x]
        for y in range(x + 1, len(inc)):
            j, a2 = inc[y]
            if i == j:
                continue
            kap = np.array([((a2 - 1) * u - (a1 - 1) * v) % p for u, v in zip(lam[(i, a1)], lam[(j, a2)])], dtype=np.int64)
            nb, inV = classify(Mp[i], Mp[j], np.array([v % p for v in lam[(i, a1)]], dtype=np.int64), np.array([v % p for v in lam[(j, a2)]], dtype=np.int64), kap, p)
            if nb < d and not inV:
                out[(a1, a2, d - nb)] += 1
    res[p] = dict(sorted(out.items()))
    print(f"prime {p}: (iii) pairs {sum(out.values())}; by (a1, a2, dim C): {res[p]}")
print("primes agree:", len(set(str(v) for v in res.values())) == 1)
