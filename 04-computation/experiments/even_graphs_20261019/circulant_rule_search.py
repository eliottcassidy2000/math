#!/usr/bin/env python3
"""Search for a uniform quantifier-free rule S -> D(S) between circulant tournaments and even circulant graphs of prime
order: D(S) = symmetrisation of {x in Z_p^* : f(x in S, a1 x + b1 in S, a2 x + b2 in S)} for boolean f and atoms from a
small dictionary; a rule passes p if it induces a bijection from tournament types onto even circulant types."""
import itertools, sys
from sympy import isprime
sys.path.insert(0, '.')
from circulant_prime import brute, order
PR = [5, 7, 11, 13, 17, 19, 23]
data = {}
for p in PR:
    tours, tbs, evens, ebs = brute(p)
    units = list(range(1, p))
    def canon(S, p=p, units=units): return min(tuple(sorted((a * x) % p for x in S)) for a in units)
    data[p] = (tours, evens, canon)
atoms = [(1, 0), (2, 0), (-2, 0), (4, 0), (1, 1), (1, -1), (2, 1), (-1, 1), (3, 0), (-3, 0)]
def test(rule_atoms, f):
    ok_ps = []
    for p in PR:
        tours, evens, canon = data[p]
        images = set(); ok = True
        for T in tours:
            S = frozenset(T)
            D = set()
            for x in range(1, p):
                bits = tuple(((a * x + b) % p) in S if (a * x + b) % p else False for a, b in rule_atoms)
                if f[bits]: D.add(x); D.add(p - x)
            img = canon(frozenset(D))
            if img not in evens or img in images: ok = False; break
            images.add(img)
        if ok and len(images) == len(evens): ok_ps.append(p)
    return ok_ps
best = []
for k in (1, 2, 3):
    for rule_atoms in itertools.combinations(atoms, k):
        if (1, 0) not in rule_atoms: continue
        keys = list(itertools.product([False, True], repeat=k))
        for vals in itertools.product([0, 1], repeat=len(keys)):
            f = dict(zip(keys, vals))
            if not any(vals) or all(vals): continue
            ok = test(rule_atoms, f)
            if len(ok) >= 3: best.append((len(ok), rule_atoms, vals, ok))
best.sort(reverse=True)
print("rules passing >= 3 of", PR, ":", len(best))
for b in best[:25]: print(b)
