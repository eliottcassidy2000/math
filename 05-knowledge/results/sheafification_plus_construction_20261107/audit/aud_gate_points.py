#!/usr/bin/env python3
"""Independent audit of THM-4892 (S85) (F): how many points (irreducible objects = nonempty strong subsets)
do the strong-component sites of T, Sigma(T) (padding) and Gamma(T) (gate) have?

Claim checked:  #points(Sigma T) = #points(T) + 2   and   #points(Gamma T) = #points(T) + 2^n + 1
(the extra points of Gamma T are {a}, {b} and {a, b} u S for every nonempty S subset V).  So the gate does
not 'collapse the site to one point'; it makes the top object irreducible while adding 2^n + 1 points.
Auditor's own helpers (aud_sites_colimit.py); nothing imported from the session.
"""
import random
from collections import Counter

from aud_sites_colimit import components, random_tour, tour_from_bits, submasks
from aud_ocf_q_gate import gate, pad


def n_points(out):
    n = len(out)
    return sum(1 for S in submasks((1 << n) - 1) if S and len(components(out, S)) == 1)


rep = Counter()
rng = random.Random(61)
cases = [tour_from_bits(n, b) for n in range(1, 5) for b in range(1 << (n * (n - 1) // 2))]
cases += [random_tour(n, rng) for n in (5, 6, 7) for _ in range(40)]
examples = []
for out in cases:
    n = len(out)
    p, ps, pg = n_points(out), n_points(pad(out)), n_points(gate(out))
    rep['#points(Sigma T) = #points(T) + 2: ' + str(ps == p + 2)] += 1
    rep['#points(Gamma T) = #points(T) + 2^n + 1: ' + str(pg == p + 2 ** n + 1)] += 1
    top_irreducible = len(components(gate(out), (1 << (n + 2)) - 1)) == 1
    rep['Gamma T top object irreducible: ' + str(top_irreducible)] += 1
    if len(examples) < 4 and n >= 3:
        examples.append((n, p, ps, pg))
for k in sorted(rep):
    print(f'{k}: {rep[k]}')
print('examples (n, #points T, #points Sigma T, #points Gamma T):', examples)
