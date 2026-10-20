#!/usr/bin/env python3
"""Power sums as uniform invariants of circulant tournaments: sigma_k(S) = sum_{s in S} s^k mod p transforms by a^k under the
multiplier a. For k coprime to p-1 the map a -> a^k is a bijection, so sigma_k(S) != 0 allows a canonical normalisation
of the multiplier orbit. K-invariant S (K a nontrivial odd multiplier subgroup) have sigma_k(S) = 0 whenever |K| does not
divide k. Question: are the "traceless" tournaments (sigma_k(S) = 0 for every k coprime to p-1) exactly the K_max-invariant
ones? Checked for primes p <= 31."""
import itertools
from math import gcd
from sympy import isprime
for p in [q for q in range(5, 32) if isprime(q)]:
    units = list(range(1, p)); classes = list(range(1, (p + 1) // 2))
    ks = [k for k in range(1, p) if gcd(k, p - 1) == 1]
    def order(a):
        x = a % p; k = 1
        while x != 1: x = x * a % p; k += 1
        return k
    Kmax = [a for a in units if order(a) % 2 == 1]
    traceless = []; types = set(); traceless_types = set(); kinv_types = set()
    for signs in itertools.product([0, 1], repeat=len(classes)):
        S = frozenset((x if s == 0 else p - x) for x, s in zip(classes, signs))
        canon = min(tuple(sorted((a * x) % p for x in S)) for a in units); types.add(canon)
        if all(sum(pow(s, k, p) for s in S) % p == 0 for k in ks): traceless_types.add(canon)
        if all(frozenset((a * x) % p for x in S) == S for a in Kmax): kinv_types.add(canon)
    # types with nontrivial stabiliser
    nontriv = set()
    for canon in types:
        S = frozenset(canon)
        if any(frozenset((a * x) % p for x in S) == S for a in units if a != 1): nontriv.add(canon)
    print(f"p = {p}: |K_max| = {len(Kmax)}, types {len(types)}, traceless types {len(traceless_types)}, "
          f"K_max-invariant types {len(kinv_types)}, nontrivial-stabiliser types {len(nontriv)}; "
          f"traceless == nontrivial-stabiliser? {traceless_types == nontriv}; traceless == K_max-invariant? {traceless_types == kinv_types}")
