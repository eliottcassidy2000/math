#!/usr/bin/env python3
"""exploration: truly bounded functionals (flags) via blocks of length <= 3, with a closure certificate."""
import itertools
import sys
from fractions import Fraction

import minrep
from lattice_tools import rank_q
from explore1 import nullspace_q, matmul, matpow

rep = minrep.build()
R = rep["R"]
A = rep["A"]
L0 = rep["L0"]
FIN = rep["FIN"]
I = [[int(i == j) for j in range(R)] for i in range(R)]

cons = []
for L in (1, 2, 3):
    for d in itertools.product(range(4), repeat=L):
        Ad = I
        for c in d:
            Ad = matmul(Ad, A[c])
        AN = matpow(Ad, R)
        AmI = [[Ad[i][j] - I[i][j] for j in range(R)] for i in range(R)]
        cons += matmul(AN, AmI)
    F = nullspace_q(cons, R)
    print(f"blocks of length <= {L}: dim of bounded-candidate functionals = {len(F)}")

# closure of the flag values under the letters
k = len(F)
# express the action on flag vectors: need A_c F in span F
def flag(v):
    return tuple(sum(v[j] * lam[j] for j in range(R)) for lam in F)

seen = {}
frontier = [((), L0)]
seen[flag(L0)] = ()
level = 0
while frontier and level < 12:
    nxt = []
    for w, v in frontier:
        for c in range(4):
            v2 = [sum(v[i] * A[c][i][j] for i in range(R)) for j in range(R)]
            fl = flag(v2)
            if fl not in seen:
                seen[fl] = w + (c,)
                nxt.append((w + (c,), v2))
    frontier = nxt
    level += 1
    print(f"  level {level}: new flag vectors {len(nxt)}, total {len(seen)}")
print("flag vectors:")
for fl, w in sorted(seen.items(), key=lambda t: (len(t[1]), t[1])):
    print("  ", "".join(map(str, w)) or "(empty)", [str(x) for x in fl])
