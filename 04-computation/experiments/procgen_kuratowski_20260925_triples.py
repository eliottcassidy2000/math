#!/usr/bin/env python3
"""procgen_kuratowski_20260925, part T: exact checks behind the triple dictionary.

The Kuratowski-Tutte pattern has two exact shapes:
  KT-a "twins + container": two incomparable minimal obstructions X, Y (K5, K33 for planarity) and a
       third Z >= X, Y in the same order (Petersen, minimal obstruction to 3-edge-colouring of
       bridgeless cubic graphs / the 4-flow conjecture);
  KT-b "dual pair + self-dual": a duality D with D(X) = Y, D(Z) = Z (Tutte's regular matroids:
       F7 <-> F7* and the self-dual U_{2,4}).
Sections (every check raises on failure):
  T1  KT-a in graphs: K5 and K33 are minors of Petersen (explicit), incomparable; Delta-Y on a triangle of
      K5 gives K33 + e; the Y-Delta/Delta-Y class of K6 has 7 graphs and contains Petersen (Petersen family).
  T2  KT-b in matroids: ranks of F7, F7*, U_{2,4} and its dual, M(K5), M(K33) and their duals.
  T3  sheets: nu T_b nu = T_{-b} on Z (b = +-1, +-3, +-5), and nu fixes the central b = 0 map on Z[1/2].
  T4  means: QM^2 + GM^2 = 2 AM^2 (parallelogram law); s -> 2AM^2 - s swaps GM^2, QM^2 and fixes AM^2.
  T5  the negative cycles -1, -5, -17: parity words, clocks K/L, Stern-Brocot positions on the path to
      log2(3), mediants, balancedness (Christoffel test).
  T6  four-vertex tournaments: the one-arc-reversal graph on the 4 classes, converse action, H values.
  T7  Kohl colours and the places: the colour of a doubling edge {m,2m} is the mod-3 sign of m, flipped by
      doubling (2 = -1 mod 3) and by negation (nu b nu = c).
Run: python3 04-computation/experiments/procgen_kuratowski_20260925_triples.py   (seconds, < 100 MB)
"""
import itertools
import math
from fractions import Fraction as Fr

import networkx as nx
import sympy as sp


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


print("=" * 100)
print("T1  KT-a in graphs")
print("=" * 100)
P = nx.petersen_graph()          # outer 0..4, inner 5..9, spokes i -- i+5
M = [(i, i + 5) for i in range(5)]
Q = nx.Graph()
for u, v in P.edges():
    cu, cv = u % 5, v % 5
    if cu != cv:
        Q.add_edge(cu, cv)
check(nx.is_isomorphic(Q, nx.complete_graph(5)), "Petersen / spokes = K5")
print("  contracting the 5 spokes (a perfect matching) of Petersen gives K5")
ok, kur = nx.check_planarity(P, counterexample=True)
dk = dict(kur.degree())
br = [v for v in dk if dk[v] == 3]
check(not ok and len(br) == 6 and max(dk.values()) == 3, "Petersen contains a K33 subdivision")
print("  Petersen contains a K33 subdivision (branch vertices", sorted(br), ") -- never a K5 subdivision (cubic)")
check(nx.complete_graph(5).number_of_nodes() < 6 and nx.complete_bipartite_graph(3, 3).number_of_edges() < 10,
      "K5, K33 incomparable")
print("  K5 and K33 are minor-incomparable (5 < 6 vertices; 9 < 10 edges)")
# Delta-Y on K5
G = nx.complete_graph(5)
G.remove_edges_from([(0, 1), (0, 2), (1, 2)])
G.add_edges_from([("x", 0), ("x", 1), ("x", 2)])
found = None
for e in G.edges():
    H = G.copy()
    H.remove_edge(*e)
    if nx.is_isomorphic(H, nx.complete_bipartite_graph(3, 3)):
        found = e
        break
check(found is not None, "Delta-Y(K5) = K33 + e")
print("  Delta-Y on a triangle of K5 gives K33 + one edge (the extra edge", found, "): the twins are one")
print("  star-triangle (Kennelly 1899, electrical) move apart, up to an edge")


def dy_moves(G):
    out = []
    for tri in (c for c in nx.enumerate_all_cliques(G) if len(c) == 3):
        H = G.copy()
        H.remove_edges_from(itertools.combinations(tri, 2))
        z = max(H.nodes()) + 1
        H.add_edges_from((z, t) for t in tri)
        out.append(H)
    for v in G.nodes():
        if G.degree(v) == 3:
            nb = list(G.neighbors(v))
            if any(G.has_edge(p, q) for p, q in itertools.combinations(nb, 2)):
                continue
            H = G.copy()
            H.remove_node(v)
            H.add_edges_from(itertools.combinations(nb, 2))
            out.append(nx.convert_node_labels_to_integers(H))
    return out


fam = [nx.complete_graph(6)]
frontier = [fam[0]]
while frontier:
    new = []
    for G_ in frontier:
        for H in dy_moves(G_):
            H = nx.convert_node_labels_to_integers(H)
            if not any(nx.is_isomorphic(H, F) for F in fam):
                fam.append(H)
                new.append(H)
    frontier = new
check(len(fam) == 7 and any(nx.is_isomorphic(F, P) for F in fam), "Petersen family: 7 graphs incl. Petersen")
print("  the Delta-Y / Y-Delta class of K6 has", len(fam), "graphs (vertices:", sorted(F.number_of_nodes() for F in fam),
      ") and contains Petersen: the Petersen family")
print("  (CITED, Wikipedia pages read 2026-09-25: these 7 are the forbidden minors for linkless embedding [RST],")
print("   and for Colin de Verdiere mu <= 4 [Lovasz-Schrijver 1998]; mu <= 3 <=> planar <=> no K5, K33 minor.)")

print()
print("=" * 100)
print("T2  KT-b in matroids (Tutte 1958/1965: regular <=> no U24, F7, F7*; graphic adds M*(K5), M*(K33))")
print("=" * 100)


def gf2_rank(vectors):
    rows = [int("".join(map(str, v)), 2) for v in vectors]
    rank = 0
    for bit in reversed(range(max(len(v) for v in vectors))):
        piv = next((r for r in rows if r >> bit & 1), None)
        if piv is None:
            continue
        rows = [r ^ piv if (r >> bit & 1) else r for r in rows if r != piv]
        rank += 1
    return rank


F7 = [v for v in itertools.product((0, 1), repeat=3) if any(v)]
rF7 = gf2_rank(F7)
print(f"  F7 (7 nonzero vectors of GF(2)^3): rank {rF7}; its dual F7* has rank 7 - {rF7} = {7 - rF7}: a dual PAIR")
print("  U_{2,4}: dual is U_{4-2,4} = U_{2,4}: SELF-DUAL")
for name, Gm in (("K5", nx.complete_graph(5)), ("K33", nx.complete_bipartite_graph(3, 3))):
    r = Gm.number_of_nodes() - 1
    print(f"  M({name}): rank {r}, corank {Gm.number_of_edges() - r}; M*({name}) is regular, not graphic")
check(rF7 == 3, "F7 rank")
print("  => Tutte's regular-matroid triple {F7, F7*, U24} is exactly KT-b: a dual pair plus a self-dual third")

print()
print("=" * 100)
print("T3  sheets b in {-1, 0, +1}")
print("=" * 100)


def Tb(b):
    return lambda x: (3 * x + b) // 2 if x % 2 else x // 2


for b in (1, 3, 5):
    for x in range(-5000, 5001):
        check(-Tb(b)(-x) == Tb(-b)(x), "nu T_b nu = T_-b")
print("  nu T_b nu = T_{-b} on Z for b = 1, 3, 5 (|x| <= 5000): nu swaps the sheets +-1")
T0 = lambda x: x / 2 if (x.denominator == 1 and x.numerator % 2 == 0) else 3 * x / 2
check(all(-T0(-Fr(k, 2 ** j)) == T0(Fr(k, 2 ** j)) for k in range(-50, 51) for j in range(4)), "T0 odd")
print("  the central sheet b = 0 (x -> x/2 or 3x/2) commutes with nu: SELF-DUAL; so {-1,0,+1} is KT-b exactly.")
print("  It is not KT-a: b = 0 is the common scaling limit (x -> lambda x, lambda -> inf) of BOTH sheets, i.e. BELOW")
print("  them (a degeneration), not a container above them.")

print()
print("=" * 100)
print("T4  means GM, AM, QM")
print("=" * 100)
x, y = sp.symbols("x y", positive=True)
AM, GM2, QM2 = (x + y) / 2, x * y, (x ** 2 + y ** 2) / 2
check(sp.simplify(QM2 + GM2 - 2 * AM ** 2) == 0, "QM^2 + GM^2 = 2 AM^2")
check(sp.simplify((2 * AM ** 2 - GM2) - QM2) == 0 and sp.simplify((2 * AM ** 2 - QM2) - GM2) == 0, "involution")
a_, b_ = sp.symbols("a b")
check(sp.expand((a_ + b_) ** 2 + (a_ - b_) ** 2 - 2 * a_ ** 2 - 2 * b_ ** 2) == 0, "parallelogram")
print("  QM^2 + GM^2 = 2 AM^2; the involution s -> 2 AM^2 - s swaps GM^2 <-> QM^2 and fixes AM^2: KT-b exactly")
print("  (the parallelogram law (a+b)^2 + (a-b)^2 = 2a^2 + 2b^2: the sign of the cross term +-2ab is the swap;")
print("   THM-4471 puts the sign law in exactly this cross term)")

print()
print("=" * 100)
print("T5  the negative cycles -1, -5, -17")
print("=" * 100)


def T(x):
    return (3 * x + 1) // 2 if x % 2 else x // 2


words = {}
for s in (-1, -5, -17):
    w = []
    x = s
    while True:
        w.append(x % 2)
        x = T(x)
        if x == s:
            break
    words[s] = w
    K, L = len(w), sum(w)
    print(f"  cycle of {s:4d}: T-word {''.join(map(str, w)):12s} K = {K:2d} steps, L = {L} odd, clock K/L = {K}/{L}, "
          f"2^K - 3^L = {2 ** K - 3 ** L}")
check([len(words[s]) for s in (-1, -5, -17)] == [1, 3, 11], "lengths 1, 3, 11")
alpha = math.log2(3)
# Stern-Brocot path to log2(3)
lo, hi = (0, 1), (1, 0)
path = []
for _ in range(14):
    med = (lo[0] + hi[0], lo[1] + hi[1])
    path.append(med)
    if med[0] / med[1] < alpha:
        lo = med
    else:
        hi = med
print("  Stern-Brocot path to log2 3 (p/q = K/L):", ", ".join(f"{p}/{q}{'<' if p / q < alpha else '>'}" for p, q in path))
pos = {kl: path.index(kl) for kl in ((1, 1), (3, 2), (11, 7))}
check(all(kl in path for kl in ((1, 1), (3, 2), (11, 7))), "clocks on the SB path")
check((3 + 8, 2 + 5) == (11, 7), "11/7 = mediant(3/2, 8/5)")
print(f"  the clocks 1/1, 3/2, 11/7 are SB-path nodes at depths {pos[(1, 1)]}, {pos[(3, 2)]}, {pos[(11, 7)]}, all LOWER")
print("  approximations: a CHAIN on one path, not twins + container; 11/7 = mediant(3/2, 8/5), and 8/5 is an upper node")


def balanced(w):
    n = len(w)
    ww = w * 3
    for L_ in range(1, n + 1):
        cnt = {sum(ww[i:i + L_]) for i in range(n)}
        if max(cnt) - min(cnt) > 1:
            return False
    return True


for s in (-1, -5, -17):
    print(f"  cycle {s}: cyclic word balanced (Christoffel/Sturmian) = {balanced(words[s])}")
check(balanced(words[-5]) and not balanced(words[-17]), "balancedness")
chris = [(math.floor((i + 1) * 7 / 11) - math.floor(i * 7 / 11)) for i in range(11)]
print("  the balanced word of slope 7/11 is", "".join(map(str, chris)), "-- the -17 word 11110111000 is not a rotation of it:")
check(all(chris[i:] + chris[:i] != words[-17] for i in range(11)), "not Christoffel")
print("  so '11/7 = 3/2 (+) 8/5' holds for the clocks only; no mediant structure on the cycle words (NUMEROLOGY)")

print()
print("=" * 100)
print("T6  four-vertex tournaments")
print("=" * 100)
V4 = range(4)
pairs = list(itertools.combinations(V4, 2))


def score(arcs):
    s = [0] * 4
    for (u, v) in arcs:
        s[u] += 1
    return tuple(sorted(s))


NAMES = {(0, 1, 2, 3): "TT", (0, 2, 2, 2): "C3+sink", (1, 1, 1, 3): "source+C3", (1, 1, 2, 2): "strong"}


def ham(arcs):
    A = set(arcs)
    return sum(1 for p in itertools.permutations(V4) if all((p[i], p[i + 1]) in A for i in range(3)))


flip = set()
Hval = {}
for bits in itertools.product((0, 1), repeat=6):
    arcs = [(u, v) if b else (v, u) for (u, v), b in zip(pairs, bits)]
    c1 = NAMES[score(arcs)]
    Hval.setdefault(c1, set()).add(ham(arcs))
    for k in range(6):
        arcs2 = list(arcs)
        arcs2[k] = (arcs[k][1], arcs[k][0])
        c2 = NAMES[score(arcs2)]
        if c1 != c2:
            flip.add(tuple(sorted((c1, c2))))
conv = {c: NAMES[tuple(sorted(3 - s for s in k))] for k, c in NAMES.items()}
print("  one-arc-reversal adjacencies between classes:", sorted(flip))
print("  converse:", conv, "  H:", {k: sorted(v) for k, v in Hval.items()})
check(conv["C3+sink"] == "source+C3" and conv["TT"] == "TT" and conv["strong"] == "strong", "converse action")
check(("C3+sink", "source+C3") not in flip, "the diamonds are not one reversal apart")
print("  => KT-b with TWO self-dual members (TT, strong) around the converse pair; H(TT) + H(strong) = 1 + 5 =")
print("     2 H(diamond): the diamonds sit at the arithmetic mean, like AM between GM and QM (a coincidence of")
print("     small numbers: THM-4472 reads the pair as time reversal)")

print()
print("=" * 100)
print("T7  Kohl colours and the places {inf, 2, 3}")
print("=" * 100)
for m in range(-3000, 3001):
    if m % 3:
        col = "b" if m % 3 == 1 else "c"
        col2 = "b" if (2 * m) % 3 == 1 else "c"
        coln = "b" if (-m) % 3 == 1 else "c"
        check(col != col2 and col != coln, "doubling and negation flip the mod-3 sign")
print("  the colour of the doubling edge {m, 2m} is the mod-3 sign (m | 3) in (Z/3)^x = {+1, -1} (b = +1, c = -1);")
print("  doubling flips it because 2 = -1 mod 3 (so doubling rays are properly 2-coloured), and nu flips it")
print("  (nu b nu = c): the archimedean sign acts on the Tait colouring through the 3-adic sign.")
print("  multiples of 3 have no mod-3 sign: they are exactly Kohl's deleted/leaf vertices 0(6), 3(6).")
print()
print("ALL CHECKS PASSED (triples)")
