#!/usr/bin/env python3
"""Exact checks for the six-seven note (mac-mini, 2026-10-06).

  A. automorphism groups of the knight torus G_n = Cay(Z_n^2, knight moves), n = 5..14 (nauty
     dreadnaut), twin classes, vertex-stabiliser orders;
  B. G_6 = C_4 x Paley(9) (tensor product) by the explicit CRT map Z_6^2 -> Z_2^2 x Z_3^2;
  C. G_7 inside F_49 = F_7[i]: the 8 knight moves are the elements of norm 5 (a coset of mu_8,
     all non-squares); nightrider = {N in NQR_7} = complement of Paley(49); queen = {N in QR_7}
     = squares = Paley(49); both srg(49,24,11,12); omega = 5+5i (order 8) preserves the knight set;
  D. the four knight direction classes are pairwise transversal iff gcd(n,30) = 1 (n = 5..60);
  E. the point-stripping ladder PSL(2,7) on P^1(F_7) (8 points) > Borel = Aut(P_7) (21)
     > split torus = Aut(P_7 - 0) = <x -> 2x> (3); B-invariant tournaments on F_7 = {P_7, reverse};
  F. octonions from the Fano lines {x, x+1, x+3}: alternative algebra check; the 3 lines through 0
     pair q with 3q (q in QR_7), so J = L_{e_0} maps e_q -> e_{3q}; <x -> 2x> permutes them cyclically;
     P_7 - 0 = two cyclic triangles O = QR_7, I = NQR_7, the antipodal arcs run I -> O;
  G. Paley P_7 Hamiltonian-path blocking census: beta = 6 = hall, 63 minimum sets by type.
Run: python3 sixseven_20261006_structure.py   (needs nauty's dreadnaut on PATH)
"""
import re, subprocess, itertools
from math import gcd
from collections import Counter
from fractions import Fraction

MOVES = [(1, 2), (2, 1), (-1, 2), (-2, 1), (1, -2), (2, -1), (-1, -2), (-2, -1)]
OK = True


def check(cond, msg):
    global OK
    print(("  ok   " if cond else "  FAIL ") + msg)
    OK &= bool(cond)


def knight(n):
    return {n * i + j: {n * ((i + a) % n) + (j + b) % n for a, b in MOVES} for i in range(n) for j in range(n)}


def dread(adj, fix0=False):
    V = len(adj)
    s = "n=%d g\n" % V + ";\n".join(" ".join(map(str, sorted(adj[u]))) for u in range(V)) + ".\n"
    if fix0:
        s += "f=[0|1:%d]\n" % (V - 1)
    s += "x\nq\n"
    out = subprocess.run(["dreadnaut"], input=s, capture_output=True, text=True).stdout
    return int(float(re.search(r"grpsize=([0-9.e+]+)", out).group(1)))


print("A. knight-torus automorphism groups (nauty)")
table = {}
for n in range(5, 15):
    G = knight(n)
    aut = dread(G)
    stab = dread(G, fix0=True)
    twins = sum(1 for c in Counter(tuple(sorted(G[u])) for u in G).values() if c > 1)
    table[n] = (aut, stab, twins)
    print(f"  n={n:2d}  |Aut| = {aut:>9}  stabiliser of a square = {stab:>7}  twin classes = {twins}")
check(table[6] == (2 ** 18 * 144, 2 ** 18 * 144 // 36, 18), "n=6: |Aut| = 2^18*144, 18 twin pairs")
check(table[7][1] == 16, "n=7: vertex stabiliser of order 16 (generic n: 8)")
check(all(table[n][1] == 8 for n in (9, 11, 12, 13, 14)), "n = 9, 11, 12, 13, 14: stabiliser = D4 (order 8)")
check(table[5][0] == 2 * 120 ** 2, "n=5: |Aut| = 28800 = |S5 wr S2| (G_5 = K5 box K5)")

print("B. G_6 = C_4 x Paley(9)")
G6 = knight(6)
lab = {(0, 0): 0, (1, 0): 1, (1, 1): 2, (0, 1): 3}
C4 = {0: {1, 3}, 1: {0, 2}, 2: {1, 3}, 3: {0, 2}}
P9 = {3 * a + b: {3 * ((a + s) % 3) + (b + t) % 3 for s, t in [(1, 1), (-1, -1), (1, -1), (-1, 1)]} for a in range(3) for b in range(3)}
T = {9 * g + h: {9 * g2 + h2 for g2 in C4[g] for h2 in P9[h]} for g in C4 for h in P9}
phi = {6 * i + j: 9 * lab[(i % 2, j % 2)] + 3 * (i % 3) + (j % 3) for i in range(6) for j in range(6)}
check(len(set(phi.values())) == 36 and all({phi[v] for v in G6[u]} == T[phi[u]] for u in G6),
      "CRT map (i,j) -> ((i,j) mod 2, (i,j) mod 3) is an isomorphism G_6 -> C_4 x P_9")
S2 = {((a % 2, b % 2)) for a, b in MOVES}
S3 = {((a % 3, b % 3)) for a, b in MOVES}
check(S2 == {(1, 0), (0, 1)} and S3 == {(1, 1), (2, 2), (1, 2), (2, 1)} and
      {((a % 2, b % 2), (a % 3, b % 3)) for a, b in MOVES} == {(x, y) for x in S2 for y in S3},
      "knight set = S_2 x S_3 exactly (8 = 2 * 4): the mod-2 rook step times the mod-3 bishop step")
check(all(G6[6 * i + j] == G6[6 * ((i + 3) % 6) + (j + 3) % 6] for i in range(6) for j in range(6)),
      "x and x + (3,3) have the same 8 neighbours (twins)")
check(P9 == {3 * a + b: {3 * ((a + s) % 3) + (b + t) % 3 for s, t in [(1, 1), (2, 2), (1, 2), (2, 1)]} for a in range(3) for b in range(3)},
      "Paley(9) = Cay(Z_3^2, {+-(1,1), +-(1,-1)}) = K3 box K3 (lines of slope 1 and -1)")

print("C. G_7 inside F_49 = F_7[i]")
QR, NQR = {1, 2, 4}, {3, 5, 6}
N = lambda a, b: (a * a + b * b) % 7
mul = lambda z, w: ((z[0] * w[0] - z[1] * w[1]) % 7, (z[0] * w[1] + z[1] * w[0]) % 7)
Km = {(a % 7, b % 7) for a, b in MOVES}
check(Km == {(a, b) for a in range(7) for b in range(7) if N(a, b) == 5}, "knight moves = the 8 elements of norm 5")
mu8 = {(a, b) for a in range(7) for b in range(7) if N(a, b) == 1}
z0 = (1, 2)
check(len(mu8) == 8 and {mul(z0, u) for u in mu8} == Km, "norm-5 set = (1+2i) * mu_8, a coset of the norm-one torus mu_8")
sq = {mul((a, b), (a, b)) for a in range(7) for b in range(7) if (a, b) != (0, 0)}
check(all(z not in sq for z in Km), "every knight move is a non-square in F_49 (5 in NQR_7)")
night = {(a, b) for a in range(7) for b in range(7) if (a, b) != (0, 0) and any(((k * m[0]) % 7, (k * m[1]) % 7) == (a, b) for k in range(1, 7) for m in MOVES)}
queen = {(a, b) for a in range(7) for b in range(7) if (a, b) != (0, 0) and (a == 0 or b == 0 or a == b or (a + b) % 7 == 0)}
check(night == {(a, b) for a in range(7) for b in range(7) if N(a, b) in NQR}, "nightrider steps = {N in NQR_7}")
check(queen == {(a, b) for a in range(7) for b in range(7) if N(a, b) in QR} == sq, "queen steps = {N in QR_7} = squares of F_49*")
check(not (night & queen) and len(night | queen) == 48, "queen + nightrider partition the 48 nonzero steps (P^1(F_7) = 4 + 4 slopes)")


def cay(S):
    return {7 * a + b: {7 * ((a + s) % 7) + (b + t) % 7 for s, t in S} for a in range(7) for b in range(7)}


def srg(adj):
    lam, mu = set(), set()
    for u in adj:
        for v in adj:
            if u < v:
                (lam if v in adj[u] else mu).add(len(adj[u] & adj[v]))
    return {len(x) for x in adj.values()}, lam, mu


check(srg(cay(queen)) == ({24}, {11}, {12}) == srg(cay(night)), "queen and nightrider torus graphs are srg(49,24,11,12) (Paley(49) and its complement)")
om = (5, 5)
p = (1, 0)
order = 0
while True:
    p = mul(p, om); order += 1
    if p == (1, 0):
        break
check(order == 8 and {mul(om, m) for m in Km} == Km, "omega = 5+5i has order 8 and maps knight moves to knight moves (the 45-degree rotation)")
# linear automorphisms of Z_n^2 preserving the knight set (the D4 point group always does)
def lin_order(M, n):
    P, k = ((1, 0), (0, 1)), 0
    while True:
        P = tuple(tuple(sum(P[i][t] * M[t][j] for t in range(2)) % n for j in range(2)) for i in range(2)); k += 1
        if P == ((1, 0), (0, 1)):
            return k


stab_lin = {}
for n in range(5, 27):
    K = {(a % n, b % n) for a, b in MOVES}
    Ms = [((a_, b_), (c_, d_)) for a_, b_, c_, d_ in itertools.product(range(n), repeat=4)
          if gcd((a_ * d_ - b_ * c_) % n, n) == 1 and {((a_ * x + b_ * y) % n, (c_ * x + d_ * y) % n) for x, y in K} == K]
    stab_lin[n] = (len(Ms), max(lin_order(M, n) for M in Ms))
print("  linear stabiliser of the knight set (size, max element order), n = 5..26:", stab_lin)
check({n: v[0] for n, v in stab_lin.items() if v[0] != 8} == {5: 32, 6: 16, 7: 16, 8: 32, 10: 16},
      "extra linear symmetries beyond D4 exactly at n = 5, 6, 7, 8, 10 (n <= 26)")
check({n for n, v in stab_lin.items() if v[1] == 8} == {5, 7},
      "an order-8 linear symmetry (a '45-degree rotation') exists only at n = 5 (degenerate, K5 box K5) and n = 7 (omega in F_49)")

print("D. transversality of the four knight direction classes")
dirs = [(1, 2), (2, 1), (1, -2), (2, -1)]
for n in range(5, 61):
    trans = all(gcd(abs(a * d - b * c), n) == 1 for (a, b), (c, d) in itertools.combinations(dirs, 2))
    if trans != (gcd(n, 30) == 1):
        check(False, f"n={n}: transversality {trans} vs gcd(n,30)=1 {gcd(n, 30) == 1}")
check(sorted(abs(a * d - b * c) for (a, b), (c, d) in itertools.combinations(dirs, 2)) == [3, 3, 4, 4, 5, 5],
      "pairwise determinants of the knight directions are 3, 3, 4, 4, 5, 5; transversal iff gcd(n,30) = 1 (n = 5..60): 7 is the least n > 1")

print("E. the point-stripping ladder PSL(2,7) > Borel > torus")
INF = 7


def mob(a, b, c, d, x):
    if x == INF:
        return INF if c == 0 else (a * pow(c, -1, 7)) % 7
    den = (c * x + d) % 7
    return INF if den == 0 else ((a * x + b) * pow(den, -1, 7)) % 7


G = set()
for a, b, c, d in itertools.product(range(7), repeat=4):
    if (a * d - b * c) % 7 == 1:
        G.add(tuple(mob(a, b, c, d, x) for x in range(8)))
check(len(G) == 168, "PSL(2,7) acting on P^1(F_7): 168 permutations")
check(len({(g[x], g[y]) for g in G for x in [0] for y in [INF]}) == 56, "2-transitive on the 8 points (so it preserves no tournament)")
B = [g for g in G if g[INF] == INF]
arc = lambda u, v: (v - u) % 7 in QR
autP7 = [p for p in itertools.permutations(range(7)) if all(arc(p[i], p[j]) == arc(i, j) for i in range(7) for j in range(7) if i != j)]
check(len(B) == 21 and {g[:7] for g in B} == set(autP7), "stabiliser of infinity = Aut(P_7) (order 21, x -> ax+b, a in QR_7)")
orbitals = {frozenset(((g[x], g[y]) for g in B)) for x in range(7) for y in range(7) if x != y}
check(len(orbitals) == 2, "the Borel has exactly 2 orbitals on ordered pairs of F_7: the invariant tournaments are P_7 and its reverse")
Tor = [g for g in B if g[0] == 0]
V6 = [1, 2, 3, 4, 5, 6]
autP70 = [p for p in itertools.permutations(V6) if all(arc(p[i], p[j]) == arc(V6[i], V6[j]) for i in range(6) for j in range(6) if i != j)]
check(len(Tor) == 3 and {tuple(g[x] for x in V6) for g in Tor} == set(autP70) and
      {tuple(g[x] for x in V6) for g in Tor} == {tuple((a * x) % 7 for x in V6) for a in QR},
      "stabiliser of 0 and infinity = Aut(P_7 - 0) = {x -> ax : a in QR_7} = <x -> 2x> (order 3): the split torus")

print("F. octonions from the Fano plane")
lines = [((x) % 7, (x + 1) % 7, (x + 3) % 7) for x in range(7)]
prod = {}
for (i, j, k) in lines:
    for (a, b, c) in [(i, j, k), (j, k, i), (k, i, j)]:
        prod[(a, b)] = (1, c); prod[(b, a)] = (-1, c)


def omul(x, y):
    # x, y: lists of 8 Fractions, index 7 = real unit, 0..6 = e_0..e_6
    z = [Fraction(0)] * 8
    for p_ in range(8):
        for q_ in range(8):
            if x[p_] == 0 or y[q_] == 0:
                continue
            c = x[p_] * y[q_]
            if p_ == 7:
                z[q_] += c
            elif q_ == 7:
                z[p_] += c
            elif p_ == q_:
                z[7] -= c
            else:
                s, r = prod[(p_, q_)]
                z[r] += s * c
    return z


import random
rng = random.Random(1)
alt = True
for _ in range(30):
    x = [Fraction(rng.randint(-3, 3)) for _ in range(8)]
    y = [Fraction(rng.randint(-3, 3)) for _ in range(8)]
    alt &= omul(x, omul(x, y)) == omul(omul(x, x), y) and omul(omul(y, x), x) == omul(y, omul(x, x))
check(alt, "the Fano table e_x e_{x+1} = e_{x+3} defines an alternative algebra (octonions), 30 random exact tests")
e = lambda k: [Fraction(int(t == k)) for t in range(8)]
J = {q: omul(e(0), e(q)) for q in range(1, 7)}
pairs = sorted(tuple(sorted((q, [t for t in range(7) if J[q][t] != 0][0]))) for q in range(1, 7))
check(set(pairs) == {tuple(sorted((q, (3 * q) % 7))) for q in QR} and all(J[q] == e((3 * q) % 7) for q in QR),
      "J = left multiplication by e_0 on span(e_1..e_6): e_q -> e_{3q} for q in QR_7 (QR -> NQR along the 3 Fano lines through 0)")
lines0 = [frozenset(l) - {0} for l in map(set, lines) if 0 in l]
check({frozenset((2 * x) % 7 for x in L) for L in lines0} == set(lines0) and all(frozenset((2 * x) % 7 for x in L) != L for L in lines0),
      "x -> 2x permutes the 3 complex lines of T_{e_0} S^6 cyclically")
O_, I_ = [1, 2, 4], [3, 5, 6]
check(arc(1, 2) and arc(2, 4) and arc(4, 1) and arc(3, 5) and arc(5, 6) and arc(6, 3), "P_7 - 0: O = QR_7 and I = NQR_7 are cyclic triangles")
check(sorted((i, o) for i in I_ for o in O_ if arc(i, o)) == [(3, 4), (5, 2), (6, 1)], "the only I -> O arcs are the antipodal arcs u -> -u")

print("G. Paley P_7: Hamiltonian-path blocking")
arcs = [(u, v) for u in range(7) for v in range(7) if arc(u, v)]


def has_hp(A):
    out = {u: [] for u in range(7)}
    for u, v in A:
        out[u].append(v)
    dp = [0] * 128
    for v in range(7):
        dp[1 << v] |= 1 << v
    for S in range(1, 128):
        if dp[S]:
            for v in range(7):
                if dp[S] >> v & 1:
                    for w in out[v]:
                        if not S >> w & 1:
                            dp[S | 1 << w] |= 1 << w
    return dp[127] != 0


beta = next(k for k in range(1, 8) if any(not has_hp([arcs[i] for i in range(21) if i not in D]) for D in itertools.combinations(range(21), k)))
mins = [D for D in itertools.combinations(range(21), beta) if not has_hp([arcs[i] for i in range(21) if i not in D])]


def kind(D):
    rem = [a for i, a in enumerate(arcs) if i not in D]
    indeg = Counter(b for a, b in rem); outdeg = Counter(a for a, b in rem)
    if any(indeg[v] == 0 and outdeg[v] == 0 for v in range(7)):
        return "star"
    if sum(indeg[v] == 0 for v in range(7)) >= 2:
        return "two sources"
    if sum(outdeg[v] == 0 for v in range(7)) >= 2:
        return "two sinks"
    outs = {v: frozenset(b for a, b in rem if a == v) for v in range(7)}
    ins = {v: frozenset(a for a, b in rem if b == v) for v in range(7)}
    for S in itertools.combinations(range(7), 3):
        if len(frozenset().union(*(outs[v] for v in S))) == 1 or len(frozenset().union(*(ins[v] for v in S))) == 1:
            return "three vertices, one common out- or in-neighbour"
    return "other"


census = Counter(kind(D) for D in mins)
print("  beta(P_7) =", beta, "; minimum blocking sets:", len(mins), dict(census))
check(beta == 6 and len(mins) == 63 and census == Counter({"star": 7, "two sources": 21, "two sinks": 21,
                                                          "three vertices, one common out- or in-neighbour": 14}),
      "beta(P_7) = 6 = hall = N - 1; 63 minimum sets, all Hall obstructions: 7 stars, 21 + 21 two sources/sinks, 14 triples")
print("H. spectra: the 7x7 knight torus is the finite Euclidean graph E_7(5), a Kloosterman / Ramanujan graph")
import cmath, math
R = 2 * math.sqrt(7)


def eig(n, a, b):
    return sum(math.cos(2 * math.pi * (a * x + b * y) / n) for x, y in MOVES)


ram = []
for n in range(5, 61):
    lam = max(abs(eig(n, a, b)) for a in range(n) for b in range(n)
              if abs(eig(n, a, b) - 8) > 1e-9 and abs(eig(n, a, b) + 8) > 1e-9)
    if lam <= R + 1e-9:
        ram.append(n)
check(ram == [5, 6, 7, 8, 10], f"Ramanujan knight tori (all nontrivial |lambda| <= 2 sqrt 7), n = 5..60: {ram}")
check(all(4 * math.cos(2 * math.pi / n) + 4 * math.cos(4 * math.pi / n) > R for n in range(12, 2000)),
      "for n >= 12 the character (1,0) alone gives 4cos(2pi/n) + 4cos(4pi/n) >= 2 + 2 sqrt 3 > 2 sqrt 7 (checked to 2000; monotone)")
e7 = lambda t: cmath.exp(2j * math.pi * t / 7)
Kl = {c: sum(e7(x + c * pow(x, -1, 7)) for x in range(1, 7)).real for c in range(1, 7)}
okk = all(abs(eig(7, a1, a2) + Kl[(3 * (a1 * a1 + a2 * a2)) % 7]) < 1e-9 for a1 in range(7) for a2 in range(7) if (a1, a2) != (0, 0))
check(okk, "every nontrivial eigenvalue of G_7 at the character (a1,a2) is -Kl_7(3 N(a)), Kl_7(c) = sum_x e((x + c/x)/7)")
check(max(abs(v) for v in Kl.values()) <= R, f"Weil bound |Kl_7(c)| <= 2 sqrt 7 (max {max(abs(v) for v in Kl.values()):.4f}): G_7 is Ramanujan because d - 1 = 7 = q")
check(sorted(set(n for n in range(5, 61) if len({(a % n, b % n) for a, b in MOVES}) == 8 and
                 {(a % n, b % n) for a, b in MOVES} == {(x, y) for x in range(n) for y in range(n) if (x * x + y * y) % n == 5 % n}
                 and n in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59))) == [7],
      "among primes p <= 59 the knight set is the whole 'sphere' {x^2 + y^2 = 5} only for p = 7 (p + 1 = 8)")

print("I. the modular-curve reading of the ladder (PARI/GP)")
gpout = subprocess.run(["gp", "-q"], input='E=ellinit([1,-1,0,-2,-1]); print(ellglobalred(E)[1]); print(E.j); print(mfdim([49,2],1)); print(mfdim([7,2],1)); print(Set(apply(p->p%7, select(p->kronecker(-7,p)==1, primes(200)))));quit\n',
                       capture_output=True, text=True).stdout.split()
check(gpout[:4] == ["49", "-3375", "1", "0"], "49a1: conductor 49, j = -3375 (CM by Q(sqrt(-7))); genus X_0(49) = 1, genus X_0(7) = 0")
check("".join(gpout[4:]) == "[1,2,4]", "primes splitting in Q(sqrt(-7)) are 1, 2, 4 mod 7 = QR_7 = the trivial cycle's code")

print("ALL CHECKS PASSED" if OK else "SOME CHECK FAILED")
