#!/usr/bin/env python3
"""Power-sum detector at p = 61 (p-1 = 60 = 4*3*5): meet-in-the-middle over the 2^30 tournament sets S (one of each pair
{i,-i}) for the condition sigma_k(S) = 0 for every k coprime to 60; count traceless sets and compare with sets having a
nontrivial multiplier stabiliser; exhibit a traceless set with trivial stabiliser (the odd 'slice-mixed' construction)."""
from math import gcd
from collections import defaultdict
p = 61; units = [k for k in range(1, p - 1) if gcd(k, p - 1) == 1]   # 16 exponents
classes = list(range(1, (p + 1) // 2))                                 # 30 classes
half = 15
def vec(signs, cls):
    return tuple(sum((x if s == 0 else -x) * pow(x, k - 1, p) for x, s in zip(cls, signs)) % p for k in units)
# sigma_k(S) = sum over classes of (+x^k if x in S else (-x)^k = -x^k)  (k odd) 
import itertools
left = defaultdict(list)
for signs in itertools.product([0, 1], repeat=half):
    left[vec(signs, classes[:half])].append(signs)
traceless = []
for signs in itertools.product([0, 1], repeat=half):
    v = vec(signs, classes[half:]); need = tuple((-a) % p for a in v)
    for sl in left.get(need, []): traceless.append(sl + signs)
print("traceless tournament sets:", len(traceless))
def stab(S):
    return [a for a in range(1, p) if frozenset(a * x % p for x in S) == S]
sets = [frozenset((x if s == 0 else p - x) for x, s in zip(classes, signs)) for signs in traceless]
nontriv = sum(1 for S in sets if len(stab(S)) > 1)
print("of which with nontrivial multiplier stabiliser:", nontriv, "; with trivial stabiliser:", len(sets) - nontriv)
# explicit construction: Z_61^* = <2>; t = log_2(x) in Z_60 = (t mod 4, t mod 3, t mod 5)
g = 2; log = {pow(g, t, p): t for t in range(p - 1)}
phi = {0: 1, 1: 1, 2: -1, 3: 1, 4: -1}; psi = {0: 1, 1: -1, 2: 1}
def f(t):
    t4, t3, t5 = t % 4, t % 3, t % 5
    if t4 == 0: return phi[t5]
    if t4 == 1: return psi[t3]
    if t4 == 2: return -phi[t5]
    return -psi[t3]
S = frozenset(x for x in range(1, p) if f(log[x]) == 1)
assert all((x in S) != ((p - x) in S) for x in range(1, p))
print("constructed set: traceless?", all(sum(pow(s, k, p) for s in S) % p == 0 for k in units), "stabiliser:", stab(S))
