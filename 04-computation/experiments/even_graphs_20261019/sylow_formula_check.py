#!/usr/bin/env python3
"""Theorem S closed form: tau_T(P) = tau_E(P) = (1/|N(P)|) sum_{h in N(P), |<P,h>| odd} 2^{c_2(<P,h>)} for Sylow p-subgroups P,
compared with the census type counts (subgroups_fix*.out)."""
import itertools, sys
sys.path.insert(0, '.')
from subgroups_fix import cyc, mul, closure, pair_orbits
def inv(a):
    r = [0] * len(a)
    for i, x in enumerate(a): r[x] = i
    return tuple(r)
cases = []
for n in (6, 7, 8):
    cases.append((n, "Z3xZ3 (Sylow-3)", closure([cyc([[0,1,2]], n), cyc([[3,4,5]], n)], n)))
    cases.append((n, "Z5 (Sylow-5)", closure([cyc([[0,1,2,3,4]], n)], n)))
for n in (7, 8):
    cases.append((n, "Z7 (Sylow-7)", closure([cyc([[0,1,2,3,4,5,6]], n)], n)))
for n, name, P in cases:
    perms = list(itertools.permutations(range(n)))
    N = [s for s in perms if frozenset(mul(mul(s, g), inv(s)) for g in P) == P]
    tot = 0
    for h in N:
        G = closure(list(P) + [h], n)
        if len(G) % 2 == 0: continue
        tot += 2 ** len(pair_orbits(G, n))
    assert tot % len(N) == 0
    print(f"n={n} P={name} |N(P)|={len(N)}: formula = {tot // len(N)}")
