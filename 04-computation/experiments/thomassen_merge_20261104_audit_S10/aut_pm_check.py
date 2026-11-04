#!/usr/bin/env python3
"""Aut+-(T) (automorphisms and anti-automorphisms) for C_3, Z_5{1,2}, P_7: order, abelian?, element orders,
and whether it is isomorphic to the dihedral group of the same order (cyclic subgroup of index 2 + an involution
inverting it), i.e. a Maschke (planar) group; AGL(1,7) is not on Maschke's list (no element of order 21 or 42)."""
import itertools
from math import gcd
def circ(n, S):
    return [{(v + s) % n for s in S} for v in range(n)]
def autpm(out):
    n = len(out); G = []
    for p in itertools.permutations(range(n)):
        auto = all((p[w] in out[p[v]]) == (w in out[v]) for v in range(n) for w in range(n) if v != w)
        anti = all((p[v] in out[p[w]]) == (w in out[v]) for v in range(n) for w in range(n) if v != w)
        if auto or anti:
            G.append((p, auto))
    return G
def order(p):
    n = len(p); seen = [False] * n; o = 1
    for i in range(n):
        if not seen[i]:
            L = 0; j = i
            while not seen[j]:
                seen[j] = True; j = p[j]; L += 1
            o = o * L // gcd(o, L)
    return o
for name, T in (("C_3", circ(3, {1})), ("Z_5{1,2}", circ(5, {1, 2})), ("P_7", circ(7, {1, 2, 4}))):
    G = autpm(T)
    N = len(G); na = sum(1 for _, a in G if a)
    orders = sorted(order(p) for p, _ in G)
    comp = lambda p, q: tuple(p[q[i]] for i in range(len(p)))
    abelian = all(comp(p, q) == comp(q, p) for p, _ in G for q, _ in G)
    maxo = max(orders)
    dihedral_like = (maxo == N // 2) and (orders.count(2) == N // 2 + (1 if (N // 2) % 2 == 0 else 0))
    print("%-9s |Aut+-| = %2d, |Aut| = %2d, abelian = %s, max element order = %d, #involutions = %d, dihedral = %s"
          % (name, N, na, abelian, maxo, orders.count(2), dihedral_like))
