#!/usr/bin/env python3
"""Exhaustive Paley7 action/path check, using only the Python standard library."""
import hashlib
import itertools as it
import json
from pathlib import Path


def need(condition, label):
    if not condition:
        raise RuntimeError(label)


V = tuple(range(7))
Q = {1, 2, 4}
arcs = {(x, y) for x in V for y in V if (y - x) % 7 in Q}
perms = tuple(it.permutations(V))


def ordinary(p):
    return all((p[x], p[y]) in arcs for x, y in arcs)


def anti(p):
    return all((p[y], p[x]) in arcs for x, y in arcs)


def word(p):
    return all((x, y) in arcs for x, y in zip(p, p[1:]))


def normalize(p):
    a = pow((p[1] - p[0]) % 7, -1, 7)
    need(a in Q, 'normalizing slope must be a square')
    return tuple(a * (x - p[0]) % 7 for x in p)


def reversal_map(p):
    r = [None] * 7
    for x, y in zip(p, p[::-1]):
        r[x] = y
    return tuple(r)


G = tuple(p for p in perms if ordinary(p))
A = tuple(p for p in perms if anti(p))
F = {tuple((a*x+b) % 7 for x in V) for a in range(1, 7) for b in V}
Fplus = {tuple((a*x+b) % 7 for x in V) for a in Q for b in V}
need(len(G) == 21 and set(G) == Fplus, 'ordinary group')
need(len(A) == 21 and set(G) | set(A) == F, 'signed group')
need(set(G).isdisjoint(A), 'orientation cosets')
paths = tuple(p for p in perms if word(p))
need(len(paths) == 189, 'Hamilton word count')

# Independent subset DP, with terminal vertex retained.
dp = {(1 << v, v): 1 for v in V}
for mask in range(1, 1 << 7):
    for v in V:
        count = dp.get((mask, v), 0)
        if count:
            for w in V:
                if not mask >> w & 1 and (v, w) in arcs:
                    key = (mask | 1 << w, w)
                    dp[key] = dp.get(key, 0) + count
need(sum(dp.get((127, v), 0) for v in V) == len(paths), 'independent DP')

presentations = tuple(sorted({normalize(p) for p in paths}))
need(len(presentations) == 9, 'ordinary presentation count')
need(presentations == tuple(p for p in paths if p[:2] == (0, 1)), 'unique affine normalization')
fibers = {p: [w for w in paths if normalize(w) == p] for p in presentations}
need(all(len(f) == 21 for f in fibers.values()), 'free word action')
blue = tuple(p for p in presentations if anti(reversal_map(p)))
need(len(blue) == 3, 'blue presentation count')
need(sum(anti(reversal_map(p)) for p in paths) == 63, 'blue raw word count')

# Converse acts on directed words only after reversing their positions.
def converse_mate(p):
    return normalize(tuple(-x % 7 for x in p[::-1]))


mates = {p: converse_mate(p) for p in presentations}
need(all(mates[mates[p]] == p for p in presentations), 'converse involution')
need({p for p in presentations if mates[p] == p} == set(blue), 'fixed equals blue')
pairs = sorted({tuple(sorted((p, mates[p]))) for p in presentations if mates[p] != p})
need(len(pairs) == 3, 'three black pairs')
signed_orbits = []
unseen = set(paths)
while unseen:
    p = min(unseen)
    orbit = {tuple(g[x] for x in p) for g in G}
    orbit |= {tuple(g[x] for x in p[::-1]) for g in A}
    need(orbit <= unseen, 'signed orbits partition')
    unseen -= orbit
    signed_orbits.append(len(orbit))
need(sorted(signed_orbits) == [21, 21, 21, 42, 42, 42], 'signed path orbit sizes')

# Direct regular-subset test of the prior converse-exception bridge.
subsets = [frozenset(v for v in V if mask >> v & 1) for mask in range(128)]
regular = [s for s in subsets if sum(frozenset(p[x] for x in s) == s for p in G) == 1]
signed_regular = [s for s in subsets if sum(frozenset(p[x] for x in s) == s for p in F) == 1]
need(len(regular) == 84 and not signed_regular, 'ordinary/signed regular subsets')

result = {
    'status': 'PASS: finite-exact Paley7 actions and presentations',
    'universe': 'all 5040 permutations of F7; all 128 vertex subsets',
    'orientation': 'x -> y iff y-x modulo7 is in {1,2,4}',
    'Aut': len(G), 'anti_automorphisms': len(A), 'signed_group': len(F),
    'H': len(paths), 'H_subset_DP': sum(dp.get((127,v), 0) for v in V),
    'h': len(presentations), 'b': len(blue), 'raw_blue_words': 63,
    'deficit_h_minus_b': 6, 'black_converse_pairs': len(pairs),
    'signed_word_orbit_sizes': sorted(signed_orbits),
    'normalized_presentations': presentations, 'blue_presentations': blue,
    'black_pairs': pairs, 'ordinary_regular_subsets': len(regular),
    'signed_regular_subsets': len(signed_regular),
    'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
print(json.dumps(result, indent=2, sort_keys=True))
