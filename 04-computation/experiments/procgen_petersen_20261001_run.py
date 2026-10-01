#!/usr/bin/env python3
"""procgen_petersen_20261001_run.py -- re-verifies every computational claim of
05-knowledge/results/procgen_petersen_20261001_petersen_heawood_paley.md and prints to stdout only.
Ends with 'ALL CHECKS PASSED'.

Sections
  A  Delta-Y families (Petersen, Heawood, K3311): sizes, Delta-Y descendants, Hasse arrows, invariants
  B  Paley/Fano coordinates: Heawood = split(P7), Fano-line Delta-Y orbits, Pasch construction of the Petersen family,
     K44-e from the two codes of the trivial cycle, K331, Petersen = Heawood - point, vertex-deletion shadows,
     order-7 symmetry, Aut(P7 - 0) = Stab(0) = Stab(D)
  C  Fano colourings of the Petersen graph in Paley coordinates (Frobenius-equivariant ones)
  D  arc-HP parity over all orientations of the 7 Petersen-family graphs; Redei-complete lemma; transit lemma;
     Frobenius-symmetric oriented Petersen graphs; C engine cross-checked against an independent Python count
  E  Conway-Gordon-Sachs parity on random PL embeddings of all 7 graphs
  F  Theorem L: linear K6 and the Gale tournament, exhaustive over all order/sign patterns (exact arithmetic);
     F2: the walk theorem for linear K_{n+3} in R^n, n = 3, 5, 7, 9 (random configurations)
  G  Collatz side: Cayley arc-evenness, code digraphs and punctures, the Mersenne-Legendre circulant T_L,
     Legendre signs of cycle codes, q-interpolation, small-q cycle search
  H  two-code tournaments; the all-odd 14-vertex tournament QR_127[mu14] (three engines); anti-circulant tournaments
     (antipodal theorem, all-odd census N <= 22, prime-derived classes, the m = 3, 5 criteria, N = 26 checks)
Requirements: python3 with networkx and numpy, nauty's labelg and gentourng on PATH, gcc. Runtime about 3 minutes; the
largest child process (the N = 26 anti-circulant engine) uses about 260 MB.
The full N = 26 census is the separate script procgen_petersen_20261001_census26.py.
"""
import itertools
import math
import os
import random
import subprocess
import sys
import time
from collections import Counter, defaultdict
from fractions import Fraction

import numpy as np
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import procgen_petersen_20261001_lib as L  # noqa: E402

NCHECK = 0


def check(cond, msg):
    global NCHECK
    NCHECK += 1
    if not cond:
        print("CHECK FAILED:", msg)
        sys.exit(1)


def section(title):
    print()
    print("=" * 100)
    print(title)
    print("=" * 100)
    sys.stdout.flush()


T0 = time.time()
WORK = os.path.join(os.path.dirname(os.path.dirname(HERE)), "scratch", "procgen_petersen", "build")
EXE = L.build_engine(os.path.join(HERE, "procgen_petersen_20261001_hp.c"), WORK)

# =============================================================================================== A
section("A. Delta-Y families")
FAM = {}
for name, G0, exp_size, exp_m, exp_dy in [("K6", nx.complete_graph(6), 7, 15, 6),
                                          ("K7", nx.complete_graph(7), 20, 21, 14),
                                          ("K3311", nx.complete_multipartite_graph(3, 3, 1, 1), 58, 22, 26)]:
    fam, arrows, skips = L.closure(G0)
    c0 = L.canon(G0)
    dyd = L.dy_descendants(c0, arrows)
    ms = sorted(set(G.number_of_edges() for G in fam.values()))
    nxc = L.nx_class_count(list(fam.values()))
    print("%-6s family size %d, edge counts %s, Delta-Y descendants of the root %d, Y-Delta moves creating "
          "multi-edges %d, networkx class count %d" % (name, len(fam), ms, len(dyd), skips, nxc))
    check(len(fam) == exp_size and ms == [exp_m] and len(dyd) == exp_dy and skips == 0 and nxc == exp_size,
          "family " + name)
    FAM[name] = (fam, arrows, dyd, c0)

fam6, arrows6, dyd6, c06 = FAM["K6"]
INFO6 = {c: L.describe(G) for c, G in fam6.items()}
PET = L.canon(nx.petersen_graph())
K331 = L.canon(nx.complete_multipartite_graph(3, 3, 1))
K44e = nx.complete_bipartite_graph(4, 4)
K44e.remove_edge(0, 4)
K44E = L.canon(K44e)
check(PET in fam6 and K331 in fam6 and K44E in fam6, "named members present")
NAME6 = {c06: "K6", PET: "P10", K331: "K331", K44E: "K44-e"}
for c, d in INFO6.items():
    if c not in NAME6:
        NAME6[c] = "P%d" % d['n']
check(sorted(NAME6.values()) == sorted(["K6", "P7", "K331", "P8", "K44-e", "P9", "P10"]), "names")
ORDER6 = sorted(fam6, key=lambda c: (INFO6[c]['n'], NAME6[c]))
print("\nPetersen family invariants:")
print("  %-6s %3s %-24s %4s %6s %4s %5s %6s %4s" % ("name", "n", "degrees", "tri", "girth", "bip", "|Aut|", "HP_u", "HC"))
HPU = {}
for c in ORDER6:
    G = fam6[c]
    d = INFO6[c]
    # undirected Hamiltonian paths by DFS
    nodes = sorted(G.nodes())
    cntp = 0

    def dfs(v, seen, depth):
        global cntp
        if depth == len(nodes):
            cntp += 1
            return
        for w in G.neighbors(v):
            if w not in seen:
                seen.add(w)
                dfs(w, seen, depth + 1)
                seen.remove(w)
    for s in nodes:
        dfs(s, {s}, 1)
    HPU[c] = cntp // 2
    hc = L.count_hc_undirected(G)
    print("  %-6s %3d %-24s %4d %6s %4s %5d %6d %4d" % (NAME6[c], d['n'], "".join(str(x) for x in d['degs']), d['tri'],
                                                    d['girth'], d['bip'], d['aut'], HPU[c], hc))
    INFO6[c]['hc'] = hc
exp_aut = {"K6": 720, "P7": 36, "K331": 72, "P8": 8, "K44-e": 72, "P9": 12, "P10": 120}
exp_hpu = {"K6": 360, "P7": 324, "K331": 324, "P8": 264, "K44-e": 324, "P9": 192, "P10": 120}
exp_hc = {"K6": 60, "P7": 36, "K331": 36, "P8": 20, "K44-e": 36, "P9": 8, "P10": 0}
for c in ORDER6:
    check(INFO6[c]['aut'] == exp_aut[NAME6[c]] and HPU[c] == exp_hpu[NAME6[c]] and INFO6[c]['hc'] == exp_hc[NAME6[c]],
          "invariants " + NAME6[c])
dy_arrows = sorted(set((NAME6[a], NAME6[b]) for (a, b, k) in arrows6 if k == 'DY' and a != b))
print("\nDelta-Y Hasse arrows:", dy_arrows)
check(dy_arrows == sorted([("K6", "P7"), ("P7", "P8"), ("P7", "K44-e"), ("K331", "P8"), ("P8", "P9"), ("P9", "P10")]),
      "Hasse diagram")
check(set(NAME6[c] for c in dyd6) == {"K6", "P7", "P8", "K44-e", "P9", "P10"}, "Delta-Y descendants of K6")

fam7, arrows7, dyd7, c07 = FAM["K7"]
INFO7 = {c: L.describe(G) for c, G in fam7.items()}
HEA = L.canon(nx.heawood_graph())
check(HEA in fam7 and HEA in dyd7, "Heawood graph in the Heawood family, as a Delta-Y descendant")
print("\nHeawood family (n, degrees, |Aut|, Delta-Y descendant of K7):")
for c in sorted(fam7, key=lambda c: (INFO7[c]['n'], INFO7[c]['degs'])):
    d = INFO7[c]
    print("  %2d %-44s %5d %s%s" % (d['n'], "".join(str(x) for x in d['degs']), d['aut'], c in dyd7,
                                   "  <- Heawood graph" if c == HEA else ""))

# =============================================================================================== B
section("B. Paley / Fano coordinates")
K7 = nx.complete_graph(7)
P7 = L.paley(7)
for t in range(7):
    a, b, c = t + 1, t + 2, t + 4
    check(P7[a % 7][b % 7] and P7[b % 7][c % 7] and P7[c % 7][a % 7], "line is a directed 3-cycle")
    check(set(j for j in range(7) if P7[t][j]) == set(L.LINES[t]), "line = out-neighbourhood")
check(sorted(tuple(sorted(e)) for t in range(7) for e in itertools.combinations(sorted(L.LINES[t]), 2))
      == sorted(itertools.combinations(range(7), 2)), "lines partition K7")
H_all = L.dy_on_lines(K7, [L.LINES[t] for t in range(7)])
check(L.canon(H_all) == HEA, "Delta-Y on the 7 lines is the Heawood graph")
check(L.canon(L.split_digraph(P7)) == HEA, "split(P7) is the Heawood graph")
print("Heawood graph = Delta-Y of K7 along the 7 translates of {1,2,4} = split(P7): OK")

COL = L.collineations()
check(len(COL) == 168, "|GL(3,2)| = 168")
orbits = defaultdict(list)
for r in range(8):
    for S in itertools.combinations(range(7), r):
        key = min(tuple(sorted(tuple(sorted(p[x] for x in L.LINES[t])) for t in S)) for p in COL)
        orbits[key].append(S)
check(len(orbits) == 10, "10 orbits of line sets")
fano_graphs = []
print("Delta-Y along Fano line sets (one per GL(3,2)-orbit):")
for key, Ss in sorted(orbits.items(), key=lambda kv: (len(kv[0]), len(kv[1]))):
    S = Ss[0]
    G = L.dy_on_lines(K7, [L.LINES[t] for t in S])
    c = L.canon(G)
    check(c in dyd7, "Fano Delta-Y graph is a Delta-Y descendant")
    # orbit independence: every set in the orbit gives the same graph
    check(all(L.canon(L.dy_on_lines(K7, [L.LINES[t] for t in S2])) == c for S2 in Ss), "orbit independence")
    fano_graphs.append(c)
    print("   |S| = %d, orbit size %2d, rep %s -> n = %d, degrees %s%s" % (
        len(S), len(Ss), S, INFO7[c]['n'], "".join(map(str, INFO7[c]['degs'])), "  (Heawood)" if c == HEA else ""))
check(len(set(fano_graphs)) == 10, "10 distinct Fano Delta-Y graphs")
check(sorted(len(v) for v in orbits.values()) == sorted([1, 7, 21, 7, 28, 28, 7, 21, 7, 1]), "orbit sizes")

K6 = nx.complete_graph(7)
K6.remove_node(0)
pasch = [t for t in range(7) if 0 not in L.LINES[t]]
through0 = [t for t in range(7) if 0 in L.LINES[t]]
check(pasch == [0, 1, 2, 4] and through0 == [3, 5, 6], "lines avoiding / through 0 are t in {0} u QR / NQR")
check(all(len(L.LINES[s] & L.LINES[t]) == 1 for s, t in itertools.combinations(pasch, 2)), "Pasch: lines meet in 1 point")
print("Pasch lines (avoiding 0): t =", pasch, "; lines through 0: t =", through0)
for r in range(5):
    names = set()
    for S in itertools.combinations(pasch, r):
        names.add(NAME6[L.canon(L.dy_on_lines(K6, [L.LINES[t] for t in S]))])
    print("   Delta-Y on %d Pasch lines -> %s" % (r, names))
    check(names == {["K6", "P7", "P8", "P9", "P10"][r]}, "Pasch Delta-Y %d" % r)
Gp = L.dy_on_lines(K6, [L.LINES[t] for t in pasch])
left = sorted(tuple(sorted(e)) for e in Gp.edges() if e[0] < 7 and e[1] < 7)
check(left == sorted(tuple(sorted(L.LINES[t] - {0})) for t in through0) == [(1, 3), (2, 6), (4, 5)], "matching")
check(left == sorted(tuple(sorted((u, 3 * u % 7))) for u in L.D), "matching = {u, 3u}")
check(nx.girth(Gp) == 5 and all(d == 3 for _, d in Gp.degree()) and Gp.number_of_nodes() == 10, "cubic girth 5")
print("   matching left of K6:", left, "= {u, 3u : u in D}")
G = L.dy_on_lines(K6, [L.QR7, L.NQR7])
check(NAME6[L.canon(G)] == "K44-e", "K44-e from the two codes")
yQ, yN = [v for v in G.nodes() if v >= 7]
check(not G.has_edge(yQ, yN) and nx.is_bipartite(G), "missing edge joins the two centres")
G2 = L.dy_on_lines(K6, [L.LINES[0], L.LINES[1]])
common = list(L.LINES[0] & L.LINES[1])[0]
check(common == 2 and G2.degree(2) == 3 and sorted(G2.neighbors(2)) == [6, 7, 8], "common point 2: neighbours 6, y0, y1")
G3 = L.y_delta(G2, 2)
check(NAME6[L.canon(G3)] == "K331" and G3.degree(6) == 6, "K331 with apex 6")
print("K44-e = Delta-Y on D and -D (missing edge = the two centres); K331 = Y-Delta at point 2 of Delta-Y_{L0,L1}: OK")

S = L.split_digraph(P7)
for v in [(0, 'in'), (0, 'out')]:
    S2 = S.copy()
    S2.remove_node(v)
    check(L.canon(L.smooth(S2)) == PET, "Heawood minus a vertex, smoothed = Petersen")
GP11 = L.dy_on_lines(K7, [L.LINES[t] for t in pasch])
check(GP11.degree(0) == 6, "deg of 0 in the Pasch member")
G4 = GP11.copy()
G4.remove_node(0)
check(L.canon(G4) == PET, "Pasch member minus 0 = Petersen")
print("Petersen = Heawood - point (suppressed) = Delta-Y_{Pasch}(K7) - 0: OK")

shadow = Counter()
for c, G in fam7.items():
    for v in G.nodes():
        H = G.copy()
        H.remove_node(v)
        Hs = L.smooth(H)
        ch = L.canon(Hs) if Hs.number_of_nodes() else None
        if ch in NAME6:
            shadow[NAME6[ch]] += 1
print("vertex-deletion shadows of the Heawood family in the Petersen family:", dict(shadow))
check(set(shadow) == {"K6", "P7", "P8", "K44-e", "P9", "P10"}, "every member but K331 is a Heawood shadow")
fam3311 = FAM["K3311"][0]
sh2 = Counter()
for c, G in fam3311.items():
    for v in G.nodes():
        H = G.copy()
        H.remove_node(v)
        ch = L.canon(L.smooth(H))
        if ch in NAME6:
            sh2[NAME6[ch]] += 1
print("vertex-deletion shadows of the K3311 family in the Petersen family:", dict(sh2))
check("K331" in sh2, "K331 is a K3311 shadow")

aut6 = {NAME6[c]: INFO6[c]['aut'] for c in fam6}
aut7 = sorted(INFO7[c]['aut'] for c in fam7)
print("Aut orders, Petersen family:", aut6)
print("Aut orders, Heawood family:", aut7)
check(all(a % 7 for a in aut6.values()), "no order-7 symmetry in the Petersen family")
check(sorted(c for c in fam7 if INFO7[c]['aut'] % 7 == 0) == sorted([c07, HEA]), "order 7 only for K7 and Heawood")

AUTP7 = [(a, b) for a in range(1, 7) for b in range(7)
         if all(P7[(a * i + b) % 7][(a * j + b) % 7] == P7[i][j] for i in range(7) for j in range(7))]
check(len(AUTP7) == 21 and set(a for a, b in AUTP7) == {1, 2, 4}, "Aut(P7) affine of order 21")
stab0 = [(a, b) for a, b in AUTP7 if b == 0]
stabD = [(a, b) for a, b in AUTP7 if set((a * x + b) % 7 for x in L.D) == set(L.D)]
check(stab0 == stabD == [(1, 0), (2, 0), (4, 0)], "Stab(0) = Stab(D) = <x2>")
P7m = [row[1:] for row in P7[1:]]
autm = [p for p in itertools.permutations(range(6))
        if all(P7m[p[i]][p[j]] == P7m[i][j] for i in range(6) for j in range(6))]
check(len(autm) == 3 and all(any(all((a * (i + 1)) % 7 == p[i] + 1 for i in range(6)) for a in (1, 2, 4)) for p in autm),
      "Aut(P7 - 0) = multipliers by QR")
print("Aut(P7) = 21 affine maps; Stab(0) = Stab(D) = Aut(P7 - 0) = {x, 2x, 4x}: OK")
# arc classes of P7 - 0 and the HP counts per class
arcs = L.tour_arcs(P7m)
Hm, cm = L.py_hp_counts(6, arcs)
cls = defaultdict(list)
for (i, j) in arcs:
    u, v = i + 1, j + 1
    if u in L.QR7 and v in L.QR7:
        k = "D-internal"
    elif u in L.NQR7 and v in L.NQR7:
        k = "(-D)-internal"
    elif u in L.QR7 and v == 3 * u % 7:
        k = "u->3u (matching)"
    elif u in L.QR7 and v == 5 * u % 7:
        k = "u->5u"
    elif u in L.NQR7 and v == (7 - u):
        k = "antipodal -u->u"
    else:
        k = "?"
    cls[k].append(cm[(i, j)])
print("P7 - 0: H = %d, arc classes:" % Hm, {k: v for k, v in cls.items()})
check(Hm == 45 and "?" not in cls and all(len(v) == 3 for v in cls.values()), "five arc classes of size 3")
check(cls["antipodal -u->u"] == [23, 23, 23] and all(x == 13 for k, v in cls.items() if k != "antipodal -u->u" for x in v),
      "c = 23 on antipodal arcs, 13 elsewhere")
inpasch = set()
for t in pasch:
    for (a, b) in itertools.permutations(L.LINES[t], 2):
        if P7[a][b]:
            inpasch.add((a, b))
anti = [(7 - u, u) for u in L.D]
check(all(a in inpasch for a in anti) and all((u, 3 * u % 7) not in inpasch for u in L.D), "antipodal in Pasch, matching outside")

# =============================================================================================== C
section("C. Fano colourings of the Paley-coordinate Petersen graph")
V = ['p%d' % p for p in range(1, 7)] + ['y%d' % t for t in pasch]
E = [('p%d' % u, 'p%d' % (3 * u % 7)) for u in (1, 2, 4)]
for t in pasch:
    for p in sorted(L.LINES[t]):
        E.append(('y%d' % t, 'p%d' % p))
check(L.canon(nx.Graph(E)) == PET, "Paley-coordinate Petersen")
inc = {v: [i for i, e in enumerate(E) if v in e] for v in V}
LINESET = set(L.LINES.values())
sols = []
col = [None] * len(E)


def okv(v):
    cs = [col[i] for i in inc[v] if col[i] is not None]
    if len(set(cs)) < len(cs):
        return False
    if len(cs) == 2:
        return sum(1 for ln in LINESET if cs[0] in ln and cs[1] in ln) == 1
    if len(cs) == 3:
        return frozenset(cs) in LINESET
    return True


def bt(i):
    if i == len(E):
        sols.append(tuple(col))
        return
    for c in range(7):
        col[i] = c
        if okv(E[i][0]) and okv(E[i][1]):
            bt(i + 1)
    col[i] = None


bt(0)


def lines_used(s):
    return frozenset(frozenset(s[i] for i in inc[v]) for v in V)


hist = Counter(len(lines_used(s)) for s in sols)
print("Fano colourings:", len(sols), "lines-used histogram", sorted(hist.items()))
check(len(sols) == 28560 and dict(hist) == {4: 3360, 5: 10080, 6: 10080, 7: 5040}, "28560 Fano colourings (S15)")


def mul2(v):
    return v[0] + str((2 * int(v[1:])) % 7)


perm = [next(k for k, (x, y) in enumerate(E) if {x, y} == {mul2(a), mul2(b)}) for (a, b) in E]
eq = [s for s in sols if all(s[perm[i]] == (2 * s[i]) % 7 for i in range(len(E)))]
eq4 = [s for s in eq if len(lines_used(s)) == 4]
target = frozenset([L.LINES[0]] + [L.LINES[t] for t in through0])
print("Frobenius-equivariant colourings:", len(eq), "; using 4 lines:", len(eq4),
      "; their line sets:", set(tuple(sorted(tuple(sorted(l)) for l in lines_used(s))) for s in eq4))
check(len(eq) == 24 and len(eq4) == 6 and all(lines_used(s) == target for s in eq4), "pencil(0) + D forced")
check(all(len(lines_used(s)) in (4, 7) for s in eq), "equivariant colourings use 4 or 7 lines")

# =============================================================================================== D
section("D. Arc-HP parity over all orientations of the Petersen family")
exp_D = {"K6": (32768, 240, 0, 0, 737280), "P7": (12960, 0, 3752, 0, 331776), "K331": (11808, 0, 4592, 144, 331776),
         "P8": (12064, 0, 10628, 0, 135168), "K44-e": (9504, 0, 11156, 0, 165888), "P9": (7552, 0, 19872, 4, 49152),
         "P10": (3968, 0, 26480, 0, 15360)}
print("  %-6s %9s %8s %9s %12s %8s" % ("member", "H odd", "all-odd", "all-even", "even,H odd", "sum H"))
rng = random.Random(20261001)
for c in ORDER6:
    G = nx.convert_node_labels_to_integers(fam6[c])
    n = G.number_of_nodes()
    edges = sorted(tuple(sorted(e)) for e in G.edges())
    st, allodd = L.engine_all(EXE, n, edges)
    got = (st["HODD"], st["ALLODD"], st["ALLEVEN"], st["ALLEVEN_HODD"], st["HSUM"])
    print("  %-6s %9d %8d %9d %12d %8d" % ((NAME6[c],) + got))
    check(got == exp_D[NAME6[c]], "census " + NAME6[c])
    check(st["HSUM"] == 2 * HPU[c] * 2 ** (15 - (n - 1)), "sum of H vs undirected HPs " + NAME6[c])
    if NAME6[c] == "K6":
        cans = set(L.labelg([L.digraph6(6, L.orient(edges, mk)) for mk, _ in allodd]))
        p7m_can = L.labelg([L.digraph6(6, L.tour_arcs(P7m))])[0]
        check(cans == {p7m_can} and all(h == 45 for _, h in allodd), "all-odd K6 orientations = P7 - v")
    # Redei-complete lemma: non-complete -> an orientation with H = 0 (two non-adjacent sources)
    if NAME6[c] != "K6":
        u, v = next((a, b) for a, b in itertools.combinations(range(n), 2) if not G.has_edge(a, b))
        mask = 0
        for k, (a, b) in enumerate(edges):
            if b in (u, v):
                mask |= 1 << k
        (H0, _), = L.engine_list(EXE, n, edges, [mask])
        check(H0 == 0, "two sources give H = 0 in " + NAME6[c])
    # independent cross-check of the C engine on random orientations
    masks = [rng.randrange(1 << len(edges)) for _ in range(40)]
    for mk, (H, cc) in zip(masks, L.engine_list(EXE, n, edges, masks)):
        arcsO = L.orient(edges, mk)
        Hp, cp = L.py_hp_counts(n, arcsO)
        check(H == Hp and cc == [cp[a] for a in arcsO], "engine vs python " + NAME6[c])
print("C engine agrees with the independent Python enumeration on 280 random orientations")

# transit lemma (undirected)


def ham_paths(G):
    nodes = sorted(G.nodes())
    out = []

    def rec(path, seen):
        if len(path) == len(nodes):
            if path[0] < path[-1]:
                out.append(tuple(path))
            return
        for w in G.neighbors(path[-1]):
            if w not in seen:
                seen.add(w)
                path.append(w)
                rec(path, seen)
                path.pop()
                seen.remove(w)
    for s in nodes:
        rec([s], {s})
    return out


ntri = 0
for c in ORDER6:
    G = fam6[c]
    hps = ham_paths(G)
    for tri in L.triangles(G):
        ntri += 1
        TE = {frozenset(e) for e in itertools.combinations(tri, 2)}
        used = [sum(1 for i in range(len(P) - 1) if frozenset((P[i], P[i + 1])) in TE) for P in hps]
        Gy = L.delta_y(G, tri, max(G.nodes()) + 1)
        hpy = ham_paths(Gy)
        pred = sum(1 for u in used if u == 1) + sum((P[0] in tri) + (P[-1] in tri) for P, u in zip(hps, used) if u == 0)
        check(len(hpy) == pred, "transit lemma (paths)")
        nodes = sorted(G.nodes())
        hcs = [P for P in hps if G.has_edge(P[0], P[-1])]
        # Hamiltonian cycles of G using exactly one triangle edge, counted as cycles
        cyc_set = set()
        for P in hps:
            if G.has_edge(P[0], P[-1]):
                cyc = P + (P[0],)
                k = sum(1 for i in range(len(cyc) - 1) if frozenset((cyc[i], cyc[i + 1])) in TE)
                if k == 1:
                    cyc_set.add(frozenset(frozenset((cyc[i], cyc[i + 1])) for i in range(len(cyc) - 1)))
        check(L.count_hc_undirected(Gy) == len(cyc_set), "transit lemma (cycles)")
print("transit lemma verified on all %d triangles of the family" % ntri)

# Frobenius-symmetric orientations of the Paley-coordinate Petersen graph
pid = {p: p - 1 for p in range(1, 7)}
yid = {t: 6 + i for i, t in enumerate(pasch)}
results = []
for o0 in (0, 1):
    for o1 in range(8):
        arcsO = [(pid[u], pid[3 * u % 7]) for u in (1, 2, 4)]
        for p in sorted(L.LINES[0]):
            arcsO.append((yid[0], pid[p]) if o0 else (pid[p], yid[0]))
        base = sorted(L.LINES[1])
        for kk in range(3):
            t = pow(2, kk, 7)
            for j, p in enumerate(base):
                q = (p * pow(2, kk, 7)) % 7
                arcsO.append((yid[t], pid[q]) if (o1 >> j) & 1 else (pid[q], yid[t]))
        check(L.canon(nx.Graph(arcsO)) == PET, "oriented Petersen underlying graph")
        H, cc = L.py_hp_counts(10, arcsO)
        results.append((H, sum(x % 2 for x in cc.values())))
print("Frobenius-symmetric oriented Petersens (H, #odd arcs):", sorted(Counter(results).items()))
check(set(h for h, _ in results) == {0, 3, 6} and all(no < 15 for _, no in results), "no all-odd symmetric orientation")

# =============================================================================================== E
section("E. Conway-Gordon-Sachs parity on random PL embeddings")
nrng = np.random.default_rng(12345)
for c in ORDER6:
    G = fam6[c]
    cyc = [tuple(x) for x in nx.simple_cycles(G) if len(x) >= 3]
    pairs = [(C1, C2) for C1, C2 in itertools.combinations(cyc, 2) if not set(C1) & set(C2)]
    nl_seen = set()
    for bends in (0, 1, 2):
        for trial in range(60):
            pos = {v: nrng.normal(size=3) for v in G.nodes()}
            ep = {}
            for (u, v) in G.edges():
                pts = [pos[u]] + [nrng.normal(size=3) for _ in range(bends)] + [pos[v]]
                ep[(u, v)] = pts
                ep[(v, u)] = pts[::-1]

            def poly(C):
                pts = []
                for i in range(len(C)):
                    pts.extend(ep[(C[i], C[(i + 1) % len(C)])][:-1])
                return np.array(pts)
            tot = 0
            nl = 0
            for C1, C2 in pairs:
                x = L.lk_polygons(poly(C1), poly(C2))
                check(abs(x - round(x)) < 1e-9, "integral linking number")
                tot += int(round(x))
                nl += (round(x) != 0)
            check(tot % 2 == 1, "CGS parity " + NAME6[c])
            if bends == 0:
                nl_seen.add(nl)
    print("  %-6s disjoint cycle pairs %2d: sum of lk odd in all 180 embeddings; linked pairs in straight ones: %s"
          % (NAME6[c], len(pairs), sorted(nl_seen)))
    if NAME6[c] == "K6":
        check(nl_seen <= {1, 3}, "straight K6: 1 or 3 linked pairs")

# =============================================================================================== F
section("F. Theorem L: linear K6 and the Gale tournament (exhaustive, exact)")
dirs = [(round(997 * math.cos((k + 0.37) * math.pi / 6)), round(997 * math.sin((k + 0.37) * math.pi / 6))) for k in range(6)]
recs = []
ntc = 0
for sigma in itertools.permutations(range(6)):
    inv = [0] * 6
    for i in range(6):
        inv[sigma[i]] = i
    for signs in itertools.product((1, -1), repeat=6):
        if signs[0] == -1:
            continue  # g -> -g is a rotation: same chirotope
        g = [(signs[i] * dirs[sigma[i]][0], signs[i] * dirs[sigma[i]][1]) for i in range(6)]
        angs = sorted(math.atan2(y, x) for x, y in g)
        gaps = [angs[i + 1] - angs[i] for i in range(5)] + [angs[0] + 2 * math.pi - angs[5]]
        if max(gaps) >= math.pi:
            continue
        ntc += 1
        w = L.positive_weights(g)
        P = L.affine_dual([(w[i] * g[i][0], w[i] * g[i][1]) for i in range(6)])
        chi = L.chirotope(P)
        check(all(chi.values()), "general position")
        T = [[1 if i != j and g[i][0] * g[j][1] - g[i][1] * g[j][0] > 0 else 0 for j in range(6)] for i in range(6)]
        check(L.is_locally_transitive(T) and not L.is_transitive(T), "Gale tournament circular, non-transitive")
        eps = [signs[inv[p]] for p in range(6)]

        def epsx(p):
            p %= 12
            return eps[p] if p < 6 else -eps[p - 6]
        walk_linked = set()
        for k in range(6):
            Wk = frozenset([inv[p] for p in range(k, 6) if eps[p] == -1] + [inv[p] for p in range(0, k) if eps[p] == 1])
            if len(Wk) == 3 and epsx(k - 1) == epsx(k):
                walk_linked.add(min(Wk, frozenset(range(6)) - Wk, key=lambda s: tuple(sorted(s))))
        linked = set()
        for A, B in L.TRIANGLE_PAIRS:
            sab = [L.seg_pierces_chi(chi, B[i], B[(i + 1) % 3], A[0], A[1], A[2]) for i in range(3)]
            sba = [L.seg_pierces_chi(chi, A[i], A[(i + 1) % 3], B[0], B[1], B[2]) for i in range(3)]
            lk = sum(sab)
            check(lk in (-1, 0, 1) and lk == sum(sba), "|lk| <= 1 and symmetric")
            if ntc % 997 == 1:  # spot-check the chirotope route against direct exact geometry
                check(lk == L.lk_triangles(P, A, B), "chirotope route = direct geometry")
            pab = sum(1 for x in sab if x)
            pba = sum(1 for x in sba if x)
            check(pab == L.pierce_rule(T, A, B) and pba == L.pierce_rule(T, B, A), "piercing rule")
            check((lk != 0) == (pab == 1) == (pba == 1), "linked iff exactly one piercing")
            nonsep = L.is_transitive(L.switch(T, set(A)))
            check(nonsep or lk == 0, "linked implies non-separable")
            if lk:
                linked.add(frozenset(A))
        check(linked == walk_linked, "walk criterion")
        check(len(linked) in (1, 3), "1 or 3 linked pairs")
        recs.append((L.digraph6(6, L.tour_arcs(T)), len(linked)))
cans = L.labelg([r[0] for r in recs])
bycls = defaultdict(Counter)
for (d6, nl), cn in zip(recs, cans):
    bycls[cn][nl] += 1
print("totally cyclic order/sign patterns:", ntc, "of", 720 * 32)
check(ntc == 18720, "18720 totally cyclic patterns")
tours = subprocess.run(["gentourng", "-q", "6"], capture_output=True, check=True).stdout.decode().split()
TT = {}
for s in tours:
    A = L.parse_tourn(s, 6)
    cn = L.labelg([L.digraph6(6, L.tour_arcs(A))])[0]
    H, cc = L.py_hp_counts(6, L.tour_arcs(A))
    TT[cn] = dict(A=A, H=H, odd=sum(x % 2 for x in cc.values()), circ=L.is_locally_transitive(A),
                  trans=L.is_transitive(A))
check(len(TT) == 56, "56 tournaments on 6 vertices")
circ_nontrans = {cn for cn, d in TT.items() if d['circ'] and not d['trans']}
check(set(bycls) == circ_nontrans, "Gale tournaments = circular non-transitive 6-tournaments")
C3TT2 = [[0] * 6 for _ in range(6)]
grp = [0, 0, 1, 1, 2, 2]
for i in range(6):
    for j in range(6):
        if i != j:
            C3TT2[i][j] = (1 if i < j else 0) if grp[i] == grp[j] else (1 if (grp[j] - grp[i]) % 3 == 1 else 0)
c3tt2 = L.labelg([L.digraph6(6, L.tour_arcs(C3TT2))])[0]
p7m = L.labelg([L.digraph6(6, L.tour_arcs(P7m))])[0]
RT7m = [row[1:] for row in [[1 if i != j and (j - i) % 7 in (1, 2, 3) else 0 for j in range(7)] for i in range(7)][1:]]
rt7m = L.labelg([L.digraph6(6, L.tour_arcs(RT7m))])[0]
for cn in sorted(bycls, key=lambda x: TT[x]['H']):
    tag = {c3tt2: "C3[TT2,TT2,TT2]", rt7m: "RT7 - v"}.get(cn, "")
    print("   Gale class %-10s H = %2d, odd arcs %2d, patterns %5d, linked pairs %s %s" % (
        cn, TT[cn]['H'], TT[cn]['odd'], sum(bycls[cn].values()), dict(bycls[cn]), tag))
check(all(set(bycls[cn]) == ({3} if cn == c3tt2 else {1}) for cn in bycls), "3 linked iff C3[TT2,TT2,TT2]")
check(sorted(TT[cn]['H'] for cn in bycls) == [17, 23, 23, 41, 45], "H values of the Gale classes")
maxH = max(d['H'] for d in TT.values())
maximizers = sorted(cn for cn, d in TT.items() if d['H'] == maxH)
check(maxH == 45 and sorted(maximizers) == sorted([c3tt2, p7m]), "the two H = 45 maximizers")
check(TT[p7m]['odd'] == 15 and not TT[p7m]['circ'] and TT[c3tt2]['odd'] == 9 and TT[c3tt2]['circ'], "maximizers' parities")
print("the two 6-vertex H-maximizers (H = 45): C3[TT2,TT2,TT2] (circular, 9 odd arcs, 3 linked pairs) and "
      "P7 - v (not circular, all 15 arcs odd, Gale dual of no linear K6)")
# the piercing-rule count on non-circular tournaments is not a linking count (asymmetric / even values occur)
lam = Counter()
for cn, d in TT.items():
    A = d['A']
    vals = [(L.pierce_rule(A, X, Y) == 1, L.pierce_rule(A, Y, X) == 1) for X, Y in L.TRIANGLE_PAIRS]
    lam[(d['circ'], all(a == b for a, b in vals), sum(a for a, _ in vals))] += 1
print("piercing-rule count over all 56 tournaments (circular?, symmetric?, count):", sorted(lam.items()))
check(all(sym and cnt in (1, 3) for (circ, sym, cnt) in lam if circ), "circular: symmetric, 1 or 3")

# F2: the same walk theorem for linear K_{n+3} in R^n, n odd (floating point, random configurations)
frng = np.random.default_rng(2026)
for nd, trials in ((3, 300), (5, 300), (7, 120), (9, 60)):
    Np = nd + 3
    r = Np // 2
    bound = r if r % 2 == 1 else r - 1
    pairsN = []
    for A in itertools.combinations(range(Np), r):
        B = tuple(sorted(set(range(Np)) - set(A)))
        if A < B:
            pairsN.append((A, B))
    seen = Counter()
    for t in range(trials):
        Pf = frng.normal(size=(Np, nd))
        linked = set()
        for A, B in pairsN:
            l1 = L.lk_simplices(Pf, A, B)
            l2 = L.lk_simplices(Pf, B, A)
            check(abs(l1) <= 1 and abs(l1) == abs(l2), "|lk| <= 1 and symmetric")
            if l1:
                linked.add(frozenset(A))
        check(linked == L.walk_linked(L.gale_float(Pf)), "walk criterion in R^%d" % nd)
        check(len(linked) % 2 == 1 and len(linked) <= bound, "odd and bounded")
        seen[len(linked)] += 1
    print("linear K_%d in R^%d: %d random configurations, numbers of linked pairs %s (bound %d); walk criterion exact"
          % (Np, nd, trials, sorted(seen.items()), bound))

# =============================================================================================== G
section("G. Collatz side")
# Theorem 5.1: Cayley digraphs of odd abelian groups are arc-even
crng = random.Random(7)
ncay = 0
for n in (3, 5, 7, 9, 11, 13):
    for trial in range(12):
        Sset = [s for s in range(1, n) if crng.random() < 0.45] or [1]
        arcsC = [(x, (x + s) % n) for x in range(n) for s in Sset]
        (H, cc), = L.engine_list(EXE, n, arcsC, [0])
        check(all(x % 2 == 0 for x in cc), "Cayley Z/%d arc-even" % n)
        ncay += 1
Z33 = [(a, b) for a in range(3) for b in range(3)]
idx = {v: i for i, v in enumerate(Z33)}
for trial in range(12):
    Sset = [v for v in Z33 if v != (0, 0) and crng.random() < 0.5] or [(0, 1)]
    arcsC = [(idx[x], idx[((x[0] + s[0]) % 3, (x[1] + s[1]) % 3)]) for x in Z33 for s in Sset]
    (H, cc), = L.engine_list(EXE, 9, arcsC, [0])
    check(all(x % 2 == 0 for x in cc), "Cayley Z3xZ3 arc-even")
    ncay += 1
print("Theorem 5.1 checked on %d random Cayley digraphs (Z/n, n odd <= 13, and Z3 x Z3)" % ncay)
# control: even order fails
(Hc, ccc), = L.engine_list(EXE, 6, [(x, (x + s) % 6) for x in range(6) for s in (1, 2)], [0])
check(any(x % 2 for x in ccc), "even-order control has an odd arc")
print("control Cay(Z/6, {1,2}): H = %d, odd arcs %d (odd order is needed)" % (Hc, sum(x % 2 for x in ccc)))

print("\ncode digraphs Cay(Z/(2^L-1), C(w)) and their punctures:")
for Lw in (2, 3, 4):
    M = 2 ** Lw - 1
    for w in L.necklaces(Lw):
        r = L.code_real(w)
        C = sorted(set((r * pow(2, i, M)) % M for i in range(Lw)))
        arcsC = [(x, (x + cc_) % M) for x in range(M) for cc_ in C]
        (H, cc), = L.engine_list(EXE, M, arcsC, [0])
        arcs2 = [(x - 1, y - 1) for (x, y) in arcsC if x and y]
        (H2, cc2), = L.engine_list(EXE, M - 1, arcs2, [0])
        print("   L=%d w=%s C=%s: Cayley H=%d, odd arcs %d/%d | punctured H=%d, odd arcs %d/%d" % (
            Lw, "".join(map(str, w)), C, H, sum(x % 2 for x in cc), len(cc), H2, sum(x % 2 for x in cc2), len(cc2)))
        check(all(x % 2 == 0 for x in cc), "code digraph arc-even")
        if Lw == 3:
            check(H == 189 and set(cc) == {54} and H2 == 45 and all(x % 2 for x in cc2), "L = 3: P7 and P7 - 0")
        if Lw == 4 and C == [1, 2, 4, 8]:
            check(H == 239160 and H2 == 47208 and all(x % 2 == 0 for x in cc2), "L = 4 puncture all-even")

print("\nMersenne-Legendre circulants T_L and code subtournaments of QR_M:")
for Lw in (3, 5, 7, 13, 17, 19, 31):
    M = 2 ** Lw - 1
    check(L.is_prime(M), "Mersenne prime")
    check(L.legendre(2, M) == 1, "2 is a QR mod M")
    SL = [d for d in range(1, Lw) if L.legendre(2 ** d - 1, M) == 1]
    check(len(SL) == (Lw - 1) // 2 and all((Lw - d) not in SL for d in SL), "T_L is a tournament")
    QRL = sorted({(x * x) % Lw for x in range(1, Lw)})
    paley_like = (Lw % 4 == 3) and any(sorted((a * d) % Lw for d in SL) == QRL for a in range(1, Lw))
    rot_like = any(sorted((a * d) % Lw for d in SL) == list(range(1, (Lw + 1) // 2)) for a in range(1, Lw))
    print("   L=%2d: S_L = %s; multiplier-equivalent to Paley: %s; to the rotational tournament: %s" % (
        Lw, SL, paley_like, rot_like))
    if Lw == 3:
        check(SL == [1], "T_3 = C3")
    if Lw in (5, 7):
        check(rot_like, "T_5, T_7 rotational")
    if Lw >= 5:
        check(not paley_like, "T_L not Paley for L >= 5")
    if Lw <= 7:
        # every unit code induces T_L (c in QR) or its reverse (c in NQR) inside QR_M
        QRM = {(x * x) % M for x in range(1, M)}
        for cu in range(1, M):
            if math.gcd(cu, M) != 1:
                continue
            pts = [(cu * pow(2, i, M)) % M for i in range(Lw)]
            for i in range(Lw):
                for j in range(Lw):
                    if i != j:
                        arc = ((pts[j] - pts[i]) % M) in QRM
                        want = ((j - i) % Lw in SL) if L.legendre(cu, M) == 1 else ((i - j) % Lw in SL)
                        check(arc == want, "code subtournament")
print("every unit cycle code induces T_L (or its reverse) in QR_M for L = 3, 5, 7: OK")

print("\nLegendre sign of the real code, and integral cycles (q = 3, 5; T = non-shortcut, T1 = shortcut):")
chi_int = defaultdict(set)
chi_nonint = defaultdict(set)
for Lw in (3, 5, 7):
    M = 2 ** Lw - 1
    for w in L.necklaces(Lw):
        r = L.code_real(w)
        chi = L.legendre(r, M)
        check(L.legendre(-L.code_2adic(w), M) == -chi, "2-adic sign opposite to real sign for L <= 7")
        ints = []
        for sc in (False, True):
            if not sc and any(w[i] == 1 and w[(i + 1) % Lw] == 1 for i in range(Lw)):
                continue
            for q in (3, 5):
                x = L.cycle_point(w, q, sc)
                check(L.follows_word(x, w, q, sc), "cycle point follows its word")
                if x.denominator == 1:
                    ints.append(("T1" if sc else "T", q, int(x)))
        (chi_int if ints else chi_nonint)[Lw].add(chi)
        if ints or Lw <= 5:
            print("   L=%d w=%s real code %3d chi=%+d  integral: %s" % (Lw, "".join(map(str, w)), r, chi, ints))
check(chi_int[5] == {1, -1} and chi_nonint[5] == {1, -1} and chi_int[7] == {1, -1}, "chi does not separate")
neq = sum(1 for w in L.necklaces(13) if L.legendre(L.code_real(w), 8191) == L.legendre(-L.code_2adic(w), 8191))
print("L = 13: real and 2-adic Legendre signs agree for %d of %d necklaces (opposite for all words when L <= 7)"
      % (neq, len(L.necklaces(13))))
check(neq == 280, "L = 13 sign relation fails")

print("\nq-interpolation:")
for w, sc, exp in [((1, 0, 0), False, [3, 5]), ((1, 0), True, [3, 5]), ((1, 0, 0), True, [7, 9]),
                   ((1, 1, 0, 0, 0), True, [5])]:
    integral = [q for q in range(1, 200, 2) if L.cycle_point(w, q, sc).denominator == 1]
    print("   word %s %-12s x(1..9) = %s ; integral for odd q < 200: %s" % (
        "".join(map(str, w)), "shortcut" if sc else "non-shortcut",
        [str(L.cycle_point(w, q, sc)) for q in (1, 3, 5, 7, 9)], integral))
    check(integral == exp, "integrality set")
check(L.cycle_point((1, 0, 0), 3, False) == 1 and L.cycle_point((1, 0, 0), 5, False) == -1, "{1,4,2} and {-1,-4,-2}")
check(L.cycle_point((1, 1, 0, 0, 0), 3, True) == Fraction(5, 23), "5/23")
qr = random.Random(99)
for trial in range(400):
    Lw = qr.randint(1, 12)
    w = tuple(qr.randint(0, 1) for _ in range(Lw))
    if sum(w) == 0:
        continue
    q = qr.randrange(1, 60, 2)
    k = sum(w)
    coeff = L.cw_poly(w)
    cq = sum(cf * q ** e for e, cf in enumerate(coeff))
    check(len(coeff) <= k and all(cf >= 0 for cf in coeff) and cq > 0, "c_w(q) shape")
    if 2 ** Lw != q ** k:
        check(L.cycle_point(w, q, True) == Fraction(cq, 2 ** Lw - q ** k), "x_w(q) = c_w(q)/(2^L - q^k)")
print("x_w(q) = c_w(q)/(2^L - q^k), deg c_w <= k-1, non-negative coefficients: 400 random checks OK")

print("\npositive cycles of x -> x/2, (qx+1)/2 for odd q < 200 (starts <= 20000):")
found = {}
for q in range(3, 200, 2):
    f = L.positive_cycles(q)
    if f:
        found[q] = sorted(f.items())
        print("   q = %3d (q mod 4 = %d): cycles (min, length) %s" % (q, q % 4, found[q]))
check(sorted(found) == [3, 5, 7, 15, 31, 63, 127, 181], "q with positive cycles")
check(all(f[0][0] == 1 for q, f in found.items() if q != 181) and found[5] == [(1, 5), (13, 7), (17, 7)]
      and found[181] == [(27, 15), (35, 15)], "cycle lists")

# =============================================================================================== H
section("H. Two-code tournaments; anti-circulant tournaments; an all-odd tournament on 14 vertices")
PAR = os.path.join(WORK, "procgen_petersen_par")
subprocess.run(["gcc", "-O2", "-o", PAR, os.path.join(HERE, "procgen_petersen_20261001_par.c")], check=True)
# H.1 the bitset parity engine against the exact engine
prng = random.Random(5)
for trial in range(200):
    n = prng.randint(2, 12)
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if prng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    if trial % 3 == 0:
        for i in range(n):
            for j in range(n):
                if i != j and prng.random() < 0.2:
                    A[i][j] = 1
    arcsA = L.tour_arcs(A)
    (Hx, ccx), = L.engine_list(EXE, n, arcsA, [0])
    H2, par = L.par_engine(PAR, A)
    check(H2 == Hx % 2 and all(par[a] == c % 2 for a, c in zip(arcsA, ccx)), "parity engine vs exact engine")
print("parity engine agrees with the exact engine on 200 random digraphs")

# H.2 two-code tournaments Theta(w) = QR_M on r(w)<2> u -c(w)<2>
print("two-code tournaments Theta(w) (M = 2^L - 1 prime):")
theta_can = {}
for Lw in (3, 5, 7):
    M = 2 ** Lw - 1
    QRM = {(x * x) % M for x in range(1, M)}
    for w in L.necklaces(Lw):
        r, c = L.code_real(w), L.code_2adic(w)
        Vt = sorted({(r * pow(2, i, M)) % M for i in range(Lw)} | {(-c * pow(2, i, M)) % M for i in range(Lw)})
        idx = {v: i for i, v in enumerate(Vt)}
        A = [[1 if x != y and (y - x) % M in QRM else 0 for y in Vt] for x in Vt]
        (Hx, ccx), = L.engine_list(EXE, len(Vt), L.tour_arcs(A), [0])
        nodd = sum(x % 2 for x in ccx)
        rev = tuple(reversed(w))
        palin = any(rev[i:] + rev[:i] == w for i in range(Lw))
        theta_can[(Lw, w)] = (L.labelg([L.digraph6(len(Vt), L.tour_arcs(A))])[0], nodd == len(ccx))
        if Lw < 7 or not palin or w == (0, 0, 0, 0, 0, 0, 1):
            print("   L=%d w=%s |V|=%2d H=%d odd arcs %d/%d palindromic=%s" % (Lw, "".join(map(str, w)), len(Vt), Hx, nodd,
                                                                         len(ccx), palin))
        if Lw == 3:
            check(nodd == 15 and Hx == 45, "Theta at L = 3 is P7 - 0")
        if Lw == 5:
            check(nodd == 35 and Hx == 15505, "Theta at L = 5")
        if Lw == 7:
            check((nodd == 91) == palin and Hx == (24540117 if palin else 24641505), "Theta at L = 7: all-odd iff palindromic")
print("   L = 7: all 14 palindromic necklaces give the same all-odd tournament (H = 24540117); the 4 chiral ones give 63/91")

# H.3 the 14-vertex all-odd tournament QR_127[mu_14], independent exact check
mu14 = sorted(x for x in range(1, 127) if pow(x, 14, 127) == 1)
check(mu14 == sorted({pow(2, i, 127) for i in range(7)} | {127 - pow(2, i, 127) for i in range(7)}), "mu14 = +-<2>")
QR127 = {(x * x) % 127 for x in range(1, 127)}
A14 = [[1 if x != y and (y - x) % 127 in QR127 else 0 for y in mu14] for x in mu14]
Hpy, cpy = L.exact_hp_python(A14)
(Hx, ccx), = L.engine_list(EXE, 14, L.tour_arcs(A14), [0])
check(Hpy == Hx == 24540117 and [cpy[a] for a in L.tour_arcs(A14)] == ccx and all(x % 2 for x in ccx), "QR127[mu14] all-odd")
print("QR_127 restricted to mu14 = %s: H = %d, all 91 arc counts odd, values %s (two independent engines)"
      % (mu14, Hpy, sorted(set(cpy.values()))))
can14 = L.labelg([L.digraph6(14, L.tour_arcs(A14))])[0]
check(all(can == can14 for (Lw, w), (can, ao) in theta_can.items() if Lw == 7 and ao), "Theta(palindromic) = QR127[mu14]")
H2b, parb = L.par_engine(PAR, A14)
check(H2b == 1 and all(v == 1 for v in parb.values()), "QR127[mu14] all-odd (bitset engine, a third method)")

# H.4 anti-circulant census
print("anti-circulant tournaments T_s on Z/2m (x -> y iff s(y - x) = (-1)^x):")
allodd_classes = {}
for m in (3, 5, 7, 9, 11):
    n = 2 * m
    reps = {}
    for free in itertools.product((1, -1), repeat=m):
        sx = L.anti_full_s(m, free)
        reps.setdefault(L.anti_pattern_key(m, sx), sx)
    found = []
    for key, sx in reps.items():
        A = L.anti_tournament(m, sx)
        check(all(A[i][j] + A[j][i] == 1 for i in range(n) for j in range(i + 1, n)), "tournament")
        # +2 is an automorphism, +1 an anti-automorphism
        check(all(A[(i + 2) % n][(j + 2) % n] == A[i][j] and A[(i + 1) % n][(j + 1) % n] == A[j][i]
                  for i in range(n) for j in range(n) if i != j), "anti-circulant symmetries")
        if m <= 7:
            (Hx, ccx), = L.engine_list(EXE, n, L.tour_arcs(A), [0])
            par = {a: c % 2 for a, c in zip(L.tour_arcs(A), ccx)}
            H2 = Hx % 2
        else:
            H2, par = L.par_engine(PAR, A)
        check(H2 == 1, "Redei")
        check(all(par[(x, (x + m) % n)] == 1 for x in range(n) if A[x][(x + m) % n]), "Theorem 6.1: antipodal arcs odd")
        if all(v == 1 for v in par.values()):
            found.append(A)
    cans = sorted(set(L.labelg([L.digraph6(n, L.tour_arcs(A)) for A in found])))
    allodd_classes[n] = cans
    print("   N = %2d: %4d patterns, %3d orbit representatives, all-odd isomorphism classes: %d" % (n, 2 ** m, len(reps), len(cans)))
check([len(allodd_classes[n]) for n in (6, 10, 14, 18, 22)] == [1, 1, 1, 2, 4], "all-odd anti-circulant class counts")
check(allodd_classes[14] == [can14], "N = 14 all-odd anti-circulant class = QR127[mu14]")
p7m_can = L.labelg([L.digraph6(6, L.tour_arcs(P7m))])[0]
check(allodd_classes[6] == [p7m_can], "N = 6 class = P7 - v")

# H.5 prime-derived anti-circulants QR_p[mu_2m]
print("QR_p restricted to mu_2m (p = 3 mod 4 prime, p < 2500): classes and all-odd status")
for m in (3, 5, 7, 9, 11):
    n = 2 * m
    byclass = defaultdict(list)
    for pp in range(7, 2500):
        if pp % 4 != 3 or not L.is_prime(pp) or ((pp - 1) // 2) % m:
            continue
        A = L.mu_tournament(pp, m)
        byclass[L.labelg([L.digraph6(n, L.tour_arcs(A))])[0]].append(pp)
    for cn, ps in sorted(byclass.items(), key=lambda kv: kv[1][0]):
        ao = cn in allodd_classes[n]
        print("   N = %2d: all-odd %-5s primes %s%s" % (n, ao, ps[:10], " ..." if len(ps) > 10 else ""))
        if m == 3:
            check(ao == all(pp % 8 == 7 for pp in ps) and len(set(pp % 8 for pp in ps)) == 1, "m = 3: all-odd iff p = 7 mod 8")
        if m == 5:
            check(ao == all(pp % 8 == 3 for pp in ps) and len(set(pp % 8 for pp in ps)) == 1, "m = 5: all-odd iff p = 3 mod 8 (p < 2500)")
    nao = sum(1 for cn in byclass if cn in allodd_classes[n])
    check(nao == {3: 1, 5: 1, 7: 1, 9: 1, 11: 4}[m], "number of prime-derived all-odd classes")
    if m == 11:
        check(set(allodd_classes[22]) <= set(byclass), "every N = 22 all-odd anti-circulant class is prime-derived")
    q_minus = 2 * m + 1
    if L.is_prime(q_minus) and q_minus % 4 == 3:
        cnq = L.labelg([L.digraph6(n, L.tour_arcs(L.mu_tournament(q_minus, m)))])[0]
        check(cnq in allodd_classes[n], "QR_q - 0 all-odd for q = %d" % q_minus)

# H.5b the m = 3 and m = 5 criteria from the algebra of the cyclotomic signs s(d) = chi(zeta^d - 1):
#   m = 3 (zeta = -omega): s = (1, a, -chi(2));  m = 5 (zeta = -eta): s = (-ab, -b, ab, a, -chi(2)),
#   a = chi(1 - eta), b = chi(1 - eta^2) (using chi(1 - x^-1) = -chi(1 - x) for x of odd order).
key3 = L.anti_pattern_key(3, L.anti_full_s(3, (1, 1, -1)))
key5 = L.anti_pattern_key(5, L.anti_full_s(5, (1, 1, -1, 1, 1)))
for aa in (1, -1):
    for c2 in (1, -1):
        k = L.anti_pattern_key(3, L.anti_full_s(3, (1, aa, -c2)))
        check((k == key3) == (c2 == 1), "m = 3 class depends only on chi(2)")
        for bb in (1, -1):
            k = L.anti_pattern_key(5, L.anti_full_s(5, (-aa * bb, -bb, aa * bb, aa, -c2)))
            check((k == key5) == (c2 == -1), "m = 5 class depends only on chi(2)")
# the all-odd orbit keys above are the all-odd classes found in H.4
check(L.labelg([L.digraph6(6, L.tour_arcs(L.anti_tournament(3, L.anti_full_s(3, (1, 1, -1)))))])[0] in allodd_classes[6] and
      L.labelg([L.digraph6(10, L.tour_arcs(L.anti_tournament(5, L.anti_full_s(5, (1, 1, -1, 1, 1)))))])[0] in allodd_classes[10],
      "reference all-odd patterns")
# the predicted signs agree with direct computation for every prime in range
for pp in range(7, 2500):
    if pp % 4 != 3 or not L.is_prime(pp):
        continue
    for m in (3, 5):
        if ((pp - 1) // 2) % m:
            continue
        A = L.mu_tournament(pp, m)
        key = L.anti_pattern_key(m, {d: (1 if A[0][d] else -1) for d in range(1, 2 * m)})
        want = (L.legendre(2, pp) == 1) if m == 3 else (L.legendre(2, pp) == -1)
        check((key == (key3 if m == 3 else key5)) == want, "prime criterion m = %d" % m)
print("QR_p[mu6] all-odd iff p = 7 (mod 8), QR_p[mu10] all-odd iff p = 3 (mod 8): PROVED by the sign algebra "
      "(checked for every a, b, chi(2)) and confirmed for all p < 2500")

# H.6 automorphism groups of the all-odd classes are the shifts by 2 (order m), N <= 18
for n, cans in allodd_classes.items():
    if n > 18:
        continue
    for cn in cans:
        # rebuild a representative and count automorphisms with networkx
        m = n // 2
        rep = None
        for free in itertools.product((1, -1), repeat=m):
            A = L.anti_tournament(m, L.anti_full_s(m, free))
            if L.labelg([L.digraph6(n, L.tour_arcs(A))])[0] == cn:
                rep = A
                break
        Gd = nx.DiGraph(L.tour_arcs(rep))
        na = sum(1 for _ in nx.algorithms.isomorphism.DiGraphMatcher(Gd, Gd).isomorphisms_iter())
        print("   all-odd class N = %2d: |Aut| = %d" % (n, na))
        check(na == m, "|Aut| = m")

# H.7 N = 26 with the one-array anti-circulant engine (256 MB)
ANTI = os.path.join(WORK, "procgen_petersen_anti")
subprocess.run(["gcc", "-O2", "-o", ANTI, os.path.join(HERE, "procgen_petersen_20261001_anti.c")], check=True)
arng = random.Random(3)
for m in (3, 5, 7, 9, 11):
    frees = list(itertools.product((1, -1), repeat=m))
    arng.shuffle(frees)
    for free in frees[:8]:
        A = L.anti_tournament(m, L.anti_full_s(m, free))
        H2, par = L.anti_engine(ANTI, A)
        H2b, parb = L.par_engine(PAR, A)
        n = 2 * m
        check(H2 == H2b and all(parb[(a, b)] == pa for d, (a, b, pa) in par.items()), "anti engine vs parity engine")
        check(all(par[min((b - a) % n, (a - b) % n)][2] == pb for (a, b), pb in parb.items()), "parity constant on sigma-orbits")
print("anti-circulant engine agrees with the parity engine (40 tournaments, N <= 22) and parities are constant on sigma-orbits")
A27 = L.f27_tournament()
check(all(A27[i][j] == A27[(j + 1) % 26][(i + 1) % 26] for i in range(26) for j in range(26) if i != j), "QR27 - 0 anti-circulant")
H2, par = L.anti_engine(ANTI, A27)
check(H2 == 1 and all(p3 == 1 for (_, _, p3) in par.values()), "QR27 - 0 all-odd (THM-4524, independent engine)")
# for palindromic words of length 13 the two-code vertex set is a coset of mu26 in F_8191^*
mu26 = {x for x in range(1, 8191) if pow(x, 26, 8191) == 1}
npal = 0
for w in L.necklaces(13):
    rev = tuple(reversed(w))
    if any(rev[i:] + rev[:i] == w for i in range(13)):
        r, c = L.code_real(w), L.code_2adic(w)
        Vt = {(r * pow(2, i, 8191)) % 8191 for i in range(13)} | {(-c * pow(2, i, 8191)) % 8191 for i in range(13)}
        check(Vt == {(r * x) % 8191 for x in mu26}, "two-code set = r * mu26 for palindromic words, L = 13")
        npal += 1
print("L = 13: %d palindromic necklaces, each two-code vertex set a coset of mu26" % npal)
A8191 = L.mu_tournament(8191, 13)
H2, par = L.anti_engine(ANTI, A8191)
pat8191 = [par[d][2] for d in range(1, 14)]
check(pat8191 == [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1], "QR8191[mu26] (L = 13 two-code tournament) not all-odd")
A2003 = L.mu_tournament(2003, 13)
H2, par = L.anti_engine(ANTI, A2003)
check(all(p3 == 1 for (_, _, p3) in par.values()), "QR2003[mu26] all-odd")
c27, c2003 = L.labelg([L.digraph6(26, L.tour_arcs(A27)), L.digraph6(26, L.tour_arcs(A2003))])
check(c27 != c2003, "two different all-odd classes at N = 26")
print("N = 26: QR27 - 0 all-odd (confirms THM-4524 with an independent engine); QR2003[mu26] all-odd and not isomorphic to it;")
print("        the L = 13 two-code tournament QR8191[mu26] is NOT all-odd (parities by difference d = 1..13: %s)" % pat8191)

print()
print("%d checks, %.0f s" % (NCHECK, time.time() - T0))
print("ALL CHECKS PASSED")
