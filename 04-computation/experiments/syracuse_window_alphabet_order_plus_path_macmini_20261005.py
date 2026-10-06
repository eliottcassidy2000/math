#!/usr/bin/env python3
"""
syracuse_window_alphabet_order_plus_path_macmini_20261005.py -- the general-m window alphabet, structurally.

Object (opus S12, syracuse_window_alphabet_20261005): a window (x, Sx, ..., S^(m-1)x) of the Syracuse map
read as a tournament on m vertices: time arcs x_i -> x_(i+1), every other pair oriented by numerical order
(DESC gauge: larger -> smaller; ASC gauge: smaller -> larger).  The "m-alphabet" is the set of
isomorphism classes that occur, with their Haar frequencies (a_i i.i.d., P(a = k) = 2^-k).

Structure proved here (elementary) and checked:
 (W1) In the generic (large-x / Haar) regime the window depends only on the valuation word through the
      order pattern of v_j = j log2(3) - A_j (A_j = a_1 + ... + a_j, A_0 = 0): x_j/x_0 = 2^(v_j) (1 + o(1)).
      A letter a_k > (m-1) log2 3 - (m-2) makes every interval sum containing k negative, so letters above
      the cap C_m = ceil((m-1) log2 3 - (m-2)) can be lumped: the Haar law is an exact dyadic rational.
 (W2) The time path is a Hamiltonian path of the window tournament and its complement is the restriction
      of a linear order (hence acyclic).  Define F_m = classes of m-tournaments having a Hamiltonian path
      whose complement is acyclic ("order + path").  Then W_m (either gauge) is contained in F_m.
 (W3) Reversing all arcs maps the DESC window of (x_0..x_{m-1}) to the ASC window of the reversed sequence
      (S12 section D); the reversed sequence is a backward orbit, so the two gauges give DIFFERENT class
      alphabets (m = 4: DESC 4 classes, ASC 2).  F_m is closed under converse, so converse(W_m) is in F_m.
Computed: W_m^DESC, W_m^ASC (Haar support and exact law) for m <= 7, F_m for m <= 7 (gentourng), and the
comparison; S12's m = 4 (11:5:5:11) and m = 5 (8 classes, masses 19/64, 13/64, 1/8, 1/8, 5/64, 5/64, 3/64, 3/64)
are the controls.
Reproduce: python3 syracuse_window_alphabet_order_plus_path_macmini_20261005.py [MMAX]
"""
import sys, math, itertools, subprocess
from fractions import Fraction
from collections import defaultdict, Counter
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher

MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 7
THETA = math.log2(3)
CHECKS = 0
def check(c, msg):
    global CHECKS
    CHECKS += 1
    if not c:
        print("CHECK FAILED:", msg); sys.exit(1)

def window_arcs(word, gauge):
    """word = (a_1..a_{m-1}); returns frozenset of arcs on vertices 0..m-1 (time order)."""
    m = len(word) + 1
    v = [0.0]; A = 0
    for a in word:
        A += a; v.append(len(v) * THETA - A)
    arcs = set()
    for i in range(m):
        for j in range(i + 1, m):
            if j == i + 1:
                arcs.add((i, j))                       # time arc
            else:
                big, small = (i, j) if v[i] > v[j] else (j, i)
                arcs.add((big, small) if gauge == "DESC" else (small, big))
    return frozenset(arcs)

def scores(arcs, m):
    return tuple(sorted(Counter(a for a, _ in arcs).get(i, 0) for i in range(m)))

class Classes:
    def __init__(self, m):
        self.m = m; self.reps = []; self.buckets = defaultdict(list)
    def key(self, arcs):
        G = nx.DiGraph(list(arcs))
        return (scores(arcs, self.m), nx.weisfeiler_lehman_graph_hash(G, iterations=3))
    def find(self, arcs):
        G = nx.DiGraph(list(arcs)); G.add_nodes_from(range(self.m))
        for idx in self.buckets[self.key(arcs)]:
            if DiGraphMatcher(self.reps[idx][1], G).is_isomorphic():
                return idx
        return None
    def add(self, arcs):
        idx = self.find(arcs)
        if idx is None:
            G = nx.DiGraph(list(arcs)); G.add_nodes_from(range(self.m))
            idx = len(self.reps); self.reps.append((arcs, G)); self.buckets[self.key(arcs)].append(idx)
        return idx

def haar_alphabet(m, gauge):
    cap = math.ceil((m - 1) * THETA - (m - 2))
    cls = Classes(m)
    mass = defaultdict(Fraction)
    labeled = {}
    for word in itertools.product(range(1, cap + 2), repeat=m - 1):
        p = Fraction(1)
        for a in word:
            p *= Fraction(1, 2 ** a) if a <= cap else Fraction(1, 2 ** cap)   # lump letters > cap
        arcs = window_arcs(word, gauge)
        if arcs not in labeled:
            labeled[arcs] = cls.add(arcs)
        mass[labeled[arcs]] += p
    check(sum(mass.values()) == 1, "Haar masses sum to 1")
    return cls, mass, len(labeled)

def gentourng(m):
    r = subprocess.run(["gentourng", "-q", str(m)], capture_output=True, text=True)
    out = []
    for line in r.stdout.split("\n"):
        line = line.strip()
        if len(line) == m * (m - 1) // 2 and set(line) <= {"0", "1"}:
            arcs = set(); k = 0
            for i in range(m):
                for j in range(i + 1, m):
                    arcs.add((i, j) if line[k] == "1" else (j, i)); k += 1
            out.append(frozenset(arcs))
    return out

def is_order_plus_path(arcs, m):
    """Does the tournament have a Hamiltonian path whose complement is acyclic?"""
    succ = defaultdict(set)
    for a, b in arcs: succ[a].add(b)
    for perm in itertools.permutations(range(m)):
        if all(perm[i + 1] in succ[perm[i]] for i in range(m - 1)):
            path = {(perm[i], perm[i + 1]) for i in range(m - 1)}
            rest = nx.DiGraph([a for a in arcs if a not in path]); rest.add_nodes_from(range(m))
            if nx.is_directed_acyclic_graph(rest):
                return perm
    return None

def describe(arcs, m):
    sc = scores(arcs, m)
    c3 = math.comb(m, 3) - sum(math.comb(s, 2) for s in sc)
    G = nx.DiGraph(list(arcs))
    L = max(len(c) for c in nx.strongly_connected_components(G))
    conv = frozenset((b, a) for a, b in arcs)
    return sc, c3, L, conv

print(f"theta = log2 3 = {THETA:.6f}")
A000568 = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}
summary = []
for m in range(3, MMAX + 1):
    cap = math.ceil((m - 1) * THETA - (m - 2))
    print("\n" + "=" * 78)
    print(f"m = {m}: letters lumped above cap C_m = {cap}; {(cap + 1) ** (m - 1)} word types")
    print("=" * 78)
    res = {}
    for gauge in ("DESC", "ASC"):
        cls, mass, nlab = haar_alphabet(m, gauge)
        res[gauge] = (cls, mass, nlab)
    clsD, massD, nlabD = res["DESC"]; clsA, massA, nlabA = res["ASC"]
    # F_m from gentourng
    hosts = gentourng(m)
    check(len(hosts) == A000568[m], f"gentourng count m={m}")
    F = [arcs for arcs in hosts if is_order_plus_path(arcs, m) is not None]
    # map window classes onto gentourng classes to compare sets
    allcls = Classes(m)
    host_idx = [allcls.add(h) for h in hosts]
    F_idx = set(allcls.add(h) for h in F)
    WD_idx = set(allcls.add(rep[0]) for rep in clsD.reps)
    WA_idx = set(allcls.add(rep[0]) for rep in clsA.reps)
    check(WD_idx <= F_idx and WA_idx <= F_idx, "W_m subset of F_m")
    # gauge relation (S12 section D): converse(DESC window of x_0..x_{m-1}) = ASC window of the REVERSED
    # sequence x_{m-1}..x_0, which is a backward orbit, not a forward window; so the two class alphabets
    # need not coincide (at m = 4 DESC has 4 classes, ASC has 2).  We only record both.
    WD_conv = set()
    for rep in clsD.reps:
        conv = frozenset((b, a) for a, b in rep[0]); WD_conv.add(allcls.add(conv))
    conv_in_F = WD_conv <= F_idx
    print(f"  |W_m^DESC| = {len(WD_idx)}  |W_m^ASC| = {len(WA_idx)}  |F_m| = {len(F_idx)}  |all classes| = {len(hosts)}"
          f"   distinct labeled DESC windows = {nlabD}, ASC = {nlabA}")
    print(f"  W_m^DESC == F_m ? {WD_idx == F_idx}      W_m^DESC == W_m^ASC (as sets)? {WD_idx == WA_idx}      converse(W_m^DESC) in F_m: {conv_in_F}   |W_m^DESC u W_m^ASC| = {len(WD_idx | WA_idx)}")
    missing = sorted(F_idx - WD_idx)
    if missing:
        print(f"  classes in F_m but not in the DESC Haar alphabet: {len(missing)}")
    # table of DESC classes with Haar mass
    print("  DESC alphabet (scores, c3, largest strong comp L, self-converse?, Haar mass, in ASC alphabet?):")
    rows = []
    for idx, rep in enumerate(clsD.reps):
        sc, c3, L, conv = describe(rep[0], m)
        selfconv = clsD.find(conv) == idx
        rows.append((massD[idx], sc, c3, L, selfconv, allcls.add(rep[0]) in WA_idx))
    for mass_, sc, c3, L, selfconv, inA in sorted(rows, reverse=True):
        print(f"     {str(sc):26s} c3={c3:2d} L={L}  self-conv={str(selfconv):5s}  mass={str(mass_):>12s} = {float(mass_):.5f}   in ASC: {inA}")
    # m = 4 and m = 5 controls (S12)
    if m == 4:
        got = sorted(float(v) for v in massD.values())
        check(sorted(massD.values()) == sorted([Fraction(11, 32), Fraction(5, 32), Fraction(5, 32), Fraction(11, 32)]), "S12 m=4 DESC law 11:5:5:11")
        check(sorted(massA.values()) == sorted([Fraction(11, 16), Fraction(5, 16)]), "S12 m=4 ASC law 11/16, 5/16")
    if m == 5:
        check(sorted(massD.values(), reverse=True) == [Fraction(19, 64), Fraction(13, 64), Fraction(1, 8), Fraction(1, 8), Fraction(5, 64), Fraction(5, 64), Fraction(3, 64), Fraction(3, 64)], "S12 m=5 DESC law")
    # regular tournaments (odd m) in F_m?
    if m % 2 == 1:
        reg = [allcls.add(h) for h in hosts if scores(h, m) == ((m - 1) // 2,) * m]
        print(f"  regular {m}-tournaments: {len(reg)} classes; in F_m: {sum(1 for r in reg if r in F_idx)}; in W_m^DESC: {sum(1 for r in reg if r in WD_idx)}")
    summary.append((m, len(WD_idx), len(F_idx), len(hosts), nlabD))

print("\n" + "=" * 78)
print("summary: m, |W_m^DESC|, |F_m| (order + Hamiltonian path), |classes|, labeled DESC windows")
for row in summary:
    print("  ", row)
print(f"\nALL {CHECKS} CHECKS PASSED")
