#!/usr/bin/env python3
"""Census of Petersen-graph copies in G_q = Cay(F_(q^2), mu_(q+1)), q an odd prime (mac-mini-2026-10-18-unitpath).
Counts all injective homomorphisms of the Petersen graph with vertex 0 -> 0 and vertex 1 -> 1 (every copy of a directed
edge is equivalent to (0, 1) under translations and multiplication by mu_(q+1)). F_(q^2) = F_q[t]/(t^2 - n), n a
non-residue (n = -1 when q = 3 mod 4); norm N(a + b t) = a^2 - n b^2.
Usage: python3 petersen_census.py QMAX [count|exist]"""
import sys
from sympy import primerange
QMAX = int(sys.argv[1]); mode = sys.argv[2] if len(sys.argv) > 2 else 'count'
for q in primerange(3, QMAX + 1):
    n = next(x for x in range(2, q) if pow(x, (q - 1) // 2, q) == q - 1)
    mu = [(a, b) for a in range(q) for b in range(q) if (a * a - n * b * b) % q == 1]
    assert len(mu) == q + 1
    V = q * q
    enc = lambda a, b: a * q + b
    nb = [None] * V
    for a in range(q):
        for b in range(q):
            nb[enc(a, b)] = frozenset(enc((a + x) % q, (b + y) % q) for x, y in mu)
    z0, one = enc(0, 0), enc(1, 0)
    count = 0; first = None
    for v2 in nb[one] - {z0}:
        for v3 in nb[v2] - {z0, one}:
            for v4 in (nb[v3] & nb[z0]) - {one, v2}:
                used = {z0, one, v2, v3, v4}
                for v5 in nb[z0] - used:
                    for v7 in (nb[v2] & nb[v5]) - used:
                        for v9 in (nb[v4] & nb[v7]) - used - {v5}:
                            for v6 in (nb[one] & nb[v9]) - used - {v5, v7}:
                                for v8 in (nb[v3] & nb[v6] & nb[v5]) - used - {v5, v7, v9}:
                                    count += 1
                                    if first is None: first = (v2, v3, v4, v5, v6, v7, v8, v9)
                                    if mode == 'exist': break
                                if mode == 'exist' and count: break
                            if mode == 'exist' and count: break
                        if mode == 'exist' and count: break
                    if mode == 'exist' and count: break
                if mode == 'exist' and count: break
            if mode == 'exist' and count: break
        if mode == 'exist' and count: break
    print('q = %3d (q mod 20 = %2d, q mod 12 = %2d, q mod 5 = %d): Petersen copies with 0 -> 0, 1 -> 1: %s'
          % (q, q % 20, q % 12, q % 5, count if mode == 'count' else (count > 0)), flush=True)
