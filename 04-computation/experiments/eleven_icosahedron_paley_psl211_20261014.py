#!/usr/bin/env python3
"""The eleven-square torus, the icosahedron and the Paley tournament P_11 inside PSL(2, 11) (exact computations).

Checks:
 (1) Golden numerology made exact: N(phi^n - 1) = (-1)^n + 1 - L_n, so R_n = Z[phi]/(phi^n - 1) has order
     |L_n - 1 - (-1)^n|; R_5 has order 11 = L_5.  L_n = 2n + 1 only for n = 5 (the paper's 'exactly when n = 5'), and
     C(L_n, 2) = F_(2n) iff F_(n-3) = 1 iff n in {2, 4, 5}: L = 3, 7, 11 and C(L, 2) = 3, 21, 55 (Luo's triangular
     Fibonacci numbers); F_n = n only for n in {0, 1, 5}.
 (2) PSL(2, 11) as Moebius maps on P^1(F_11) (12 points): order 660; the Borel subgroup B = {x -> a x + b : a a nonzero
     square} = Aut(P_11) = the eleven-square group G_5 has order 55 and acts regularly on the 55 arcs (unordered pairs)
     of P_11 on F_11.
 (3) An icosahedral subgroup A ~ A_5 containing the 5-fold rotation x -> 4x: |A| = 60, A acts on the 12 points with
     stabilisers of order 5, A meets B in C_5, and PSL(2,11) = A B (exact factorisation, 60 * 55 / 5 = 660).
     The A-invariant graph of valency 5 on P^1(F_11) is the icosahedron, with infinity antipodal to 0 and the two
     pentagons around them equal to the squares {1,3,4,5,9} and the non-squares {2,6,7,8,10} (the gold and teal
     classes of the eleven-square torus figure).
Session opus-2026-10-14-S36.
"""
from itertools import product
from math import comb

# (1) golden numerology
F = [0, 1]
for _ in range(130):
    F.append(F[-1] + F[-2])
L = [2, 1]
for _ in range(60):
    L.append(L[-1] + L[-2])
print('(1) |N(phi^n - 1)| for n = 1..12:', [abs((-1) ** n + 1 - L[n]) for n in range(1, 13)])
print('    n with L_n = 2n + 1 (n <= 60):', [n for n in range(1, 61) if L[n] == 2 * n + 1])
print('    n with C(L_n, 2) = F_2n (n <= 60):', [n for n in range(1, 61) if comb(L[n], 2) == F[2 * n]],
      ' values', [(L[n], comb(L[n], 2)) for n in (2, 4, 5)])
print('    n with F_n = n (n <= 60):', [n for n in range(0, 61) if F[n] == n])

# (2)-(3) PSL(2, 11) on P^1(F_11); points 0..10 and 11 = infinity
p = 11
INF = p
sq = {x * x % p for x in range(1, p)}


def mob(a, b, c, d):
    def f(x):
        if x == INF:
            return INF if c == 0 else a * pow(c, -1, p) % p
        den = (c * x + d) % p
        if den == 0:
            return INF
        return (a * x + b) * pow(den, -1, p) % p
    return tuple(f(x) for x in range(p + 1))


G = set()
for a, b, c, d in product(range(p), repeat=4):
    if (a * d - b * c) % p == 1:
        G.add(mob(a, b, c, d))
print('(2) |PSL(2,11)| =', len(G))
B = {mob(a, b, 0, pow(a, -1, p)) for a in range(1, p) for b in range(p)}  # x -> a^2 x + ab: squares as multipliers
print('    |Borel B| =', len(B), '; B equals {x -> s x + t : s square}:',
      B == {tuple((s * x + t) % p if x < p else INF for x in range(p + 1)) for s in sq for t in range(p)})
arcs = [(i, j) for i in range(p) for j in range(p) if i != j and (j - i) % p in sq]
regular = all(len({(g[i], g[j]) for g in B}) == 55 for (i, j) in arcs[:1]) and len(arcs) == 55
print('    P_11 has 55 arcs and B acts regularly on them (orbit of one arc has size 55):', regular)


def compose(g, h):  # g after h
    return tuple(g[h[x]] for x in range(p + 1))


def generate(gens):
    e = tuple(range(p + 1))
    S = {e}
    frontier = [e]
    while frontier:
        nf = []
        for x in frontier:
            for g in gens:
                y = compose(g, x)
                if y not in S:
                    S.add(y); nf.append(y)
        frontier = nf
    return S


rot5 = mob(2, 0, 0, pow(2, -1, p))  # x -> 4x
assert rot5 == tuple((4 * x) % p if x < p else INF for x in range(p + 1))
A = None
for g in G:
    H = generate([rot5, g])
    if len(H) == 60:
        A = H
        break
print('(3) found a subgroup A of order', len(A), 'containing x -> 4x')
orb = {g[INF] for g in A}
stab = [g for g in A if g[INF] == INF]
print('    A is transitive on the 12 points:', len(orb) == 12, '; point stabiliser order:', len(stab), '; |A cap B| =', len(A & B))
AB = {compose(a, b) for a in A for b in B}
print('    PSL(2,11) = A B (exact factorisation):', AB == G)
# A-invariant graphs of valency 5: orbitals of A on ordered pairs
pairs = [(x, y) for x in range(p + 1) for y in range(p + 1) if x != y]
seen = set(); orbitals = []
for (x, y) in pairs:
    if (x, y) in seen:
        continue
    o = {(g[x], g[y]) for g in A}
    seen |= o
    orbitals.append(o)
sizes = sorted(len(o) for o in orbitals)
print('    orbital sizes of A on ordered pairs:', sizes)
for o in orbitals:
    if len(o) == 60 and all((y, x) in o for (x, y) in o):  # symmetric orbital, valency 5
        nbr_inf = sorted(y for (x, y) in o if x == INF)
        nbr_0 = sorted(y for (x, y) in o if x == 0)
        antipode = [y for y in range(p) if (INF, y) not in o]
        print('    symmetric valency-5 orbital: neighbours of infinity', nbr_inf, '; neighbours of 0', nbr_0)
