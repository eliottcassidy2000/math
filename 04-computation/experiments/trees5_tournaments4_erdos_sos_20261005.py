#!/usr/bin/env python3
"""
trees5_tournaments4_erdos_sos_20261005.py

Owner observation (2026-10-05): the 3 unlabeled trees on 5 vertices "correspond" to the
isomorphism classes of tournaments on 4 vertices once the two mutually converse classes
(scores (0,2,2,2) and (1,1,1,3)) are merged, giving 3.  Also: the Erdos--Sos threshold
(t-2)n/2 alternates between whole and half.

This script freezes the finite-exact facts:

  A. unlabeled trees on n vertices, n <= 10 (A000055), with the 5-vertex census:
     degree sequences, diameter, leaves, automorphism order, cycle lengths created by
     adding one edge (trees are the maximal acyclic graphs).
  B. tournaments on n vertices up to isomorphism (A000568), self-converse count (A002785),
     converse classes (iso + reversal) for n <= 6; comparison with trees on n+1 vertices.
  C. the arc-reversal metagraph on 4 vertices (nodes = iso classes, weight = number of arcs
     of a representative whose reversal lands in the class), and its converse quotient;
     the edge-rotation graph of the 3 trees on 5 vertices, for honest comparison.
  D. Erdos--Sos for t = 5: ex(n, T) for each of the three trees, all n <= 7 (graph atlas,
     exhaustive), the half/whole alternation of 3n/2, and the regular-graph (handshake)
     obstruction; which tree attains floor(3n/2).
  E. the forcing threshold E*(n) = floor(3n/2) + 1 written as a Collatz composite:
     (3n+1)/2 for odd n, 3(n/2)+1 for even n.

Universe/filters are explicit in each section.  No external data; OEIS values quoted for
comparison only.  Reproduce:  python3 trees5_tournaments4_erdos_sos_20261005.py
"""
import itertools, math, sys
from collections import Counter, defaultdict
import networkx as nx
from networkx.generators.atlas import graph_atlas_g

CHECKS = 0
def check(cond, msg):
    global CHECKS
    CHECKS += 1
    if not cond:
        print("CHECK FAILED:", msg)
        sys.exit(1)

print("=" * 78)
print("A. Unlabeled trees on n vertices (networkx.nonisomorphic_trees), n = 1..10")
print("=" * 78)
A000055 = {1: 1, 2: 1, 3: 1, 4: 2, 5: 3, 6: 6, 7: 11, 8: 23, 9: 47, 10: 106}
tree_counts = {}
for n in range(1, 11):
    if n == 1:
        c = 1
    elif n == 2:
        c = 1
    else:
        c = sum(1 for _ in nx.nonisomorphic_trees(n))
    tree_counts[n] = c
    check(c == A000055[n], f"trees on {n}")
print("n      :", " ".join(f"{n:4d}" for n in range(1, 11)))
print("trees  :", " ".join(f"{tree_counts[n]:4d}" for n in range(1, 11)), "  (= A000055)")

def tree_name(T):
    degs = sorted((d for _, d in T.degree()), reverse=True)
    if degs[0] == T.number_of_nodes() - 1:
        return "star K_{1,%d}" % (T.number_of_nodes() - 1)
    if degs[0] == 2:
        return "path P_%d" % T.number_of_nodes()
    return "spider/fork deg=" + str(degs)

print("\nThe three trees on 5 vertices:")
trees5 = list(nx.nonisomorphic_trees(5))
for T in trees5:
    degs = sorted((d for _, d in T.degree()), reverse=True)
    diam = nx.diameter(T)
    leaves = sum(1 for _, d in T.degree() if d == 1)
    aut = sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(T, T).isomorphisms_iter())
    # cycle lengths created by adding one non-edge: distance + 1
    created = sorted(Counter(nx.shortest_path_length(T, u, v) + 1
                             for u, v in itertools.combinations(T.nodes, 2)
                             if not T.has_edge(u, v)).items())
    print(f"  {tree_name(T):22s} degrees={degs} diameter={diam} leaves={leaves} |Aut|={aut}"
          f"  labeled copies={math.factorial(5)//aut}  added-edge cycle lengths {created}")

print("\n" + "=" * 78)
print("B. Tournaments on n vertices: iso classes, self-converse, converse classes (n <= 6)")
print("=" * 78)

def tournaments_iso_classes(n):
    """Return dict canon -> representative adjacency (frozenset of arcs) for all tournaments on n."""
    pairs = list(itertools.combinations(range(n), 2))
    perms = list(itertools.permutations(range(n)))
    classes = {}
    seen_canon = set()
    for bits in itertools.product((0, 1), repeat=len(pairs)):
        arcs = frozenset((i, j) if b else (j, i) for (i, j), b in zip(pairs, bits))
        # cheap invariant bucket: score sequence
        # canonical form: lexicographically least arc-set over relabelings
        best = None
        for p in perms:
            img = tuple(sorted((p[a], p[b]) for a, b in arcs))
            if best is None or img < best:
                best = img
        if best not in seen_canon:
            seen_canon.add(best)
            classes[best] = arcs
    return classes

def canon(arcs, n, perms):
    best = None
    for p in perms:
        img = tuple(sorted((p[a], p[b]) for a, b in arcs))
        if best is None or img < best:
            best = img
    return best

def converse(arcs):
    return frozenset((b, a) for a, b in arcs)

A000568 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456}
A002785 = {1: 1, 2: 1, 3: 2, 4: 2, 5: 8, 6: 12, 7: 88}   # self-converse tournaments (OEIS)
print("n  iso(A000568)  self-converse(A002785)  converse-classes  trees(n+1)  match?")
tour_classes = {}
for n in range(1, 7):
    perms = list(itertools.permutations(range(n)))
    cl = tournaments_iso_classes(n) if n >= 2 else {(): frozenset()}
    sc = 0
    conv_classes = set()
    for key, arcs in cl.items():
        ck = canon(converse(arcs), n, perms)
        if ck == key:
            sc += 1
        conv_classes.add(frozenset([key, ck]))
    tour_classes[n] = (cl, perms)
    iso = len(cl)
    ncc = len(conv_classes)
    check(iso == A000568[n], f"A000568({n})")
    check(sc == A002785[n], f"A002785({n})")
    check(ncc == (iso + sc) // 2, "orbit count")
    tr = tree_counts[n + 1]
    print(f"{n}  {iso:12d}  {sc:22d}  {ncc:16d}  {tr:10d}  {'YES' if ncc == tr else 'no'}")
n = 7
print(f"{n}  {A000568[7]:12d}  {A002785[7]:22d}  {(A000568[7]+A002785[7])//2:16d}  {tree_counts[8]:10d}  no   (OEIS values, not recomputed here)")
print("Converse-class counts 1,1,2,3,10,34,272 vs trees 1,1,2,3,6,11,23: equal exactly for n <= 4.")

print("\n" + "=" * 78)
print("C. Arc-reversal metagraph on 4 vertices and its converse quotient")
print("=" * 78)
cl4, perms4 = tour_classes[4]

def scores(arcs, n):
    return tuple(sorted(Counter(a for a, _ in arcs).get(v, 0) for v in range(n)))

def name4(arcs):
    s = scores(arcs, 4)
    return {(0, 1, 2, 3): "TT4", (1, 1, 2, 2): "S4 (strong)", (0, 2, 2, 2): "B=C3+sink", (1, 1, 1, 3): "A=C3+source"}[s]

def three_cycles(arcs, n):
    sc = Counter(a for a, _ in arcs)
    return math.comb(n, 3) - sum(math.comb(sc.get(v, 0), 2) for v in range(n))

def ham_paths(arcs, n):
    return sum(1 for p in itertools.permutations(range(n)) if all((p[i], p[i+1]) in arcs for i in range(n-1)))

def cycle_lengths(arcs, n):
    out = set()
    for k in range(3, n + 1):
        for sub in itertools.combinations(range(n), k):
            first = sub[0]
            for rest in itertools.permutations(sub[1:]):
                cyc = (first,) + rest
                if all((cyc[i], cyc[(i+1) % k]) in arcs for i in range(k)):
                    out.add(k)
                    break
            if k in out:
                break
    return sorted(out)

keys4 = list(cl4.keys())
label = {k: name4(cl4[k]) for k in keys4}
print("class            scores      labeled  |Aut|  3-cycles  Ham.paths  cycle lengths  converse")
for k in keys4:
    arcs = cl4[k]
    aut = sum(1 for p in perms4 if frozenset((p[a], p[b]) for a, b in arcs) == arcs)
    ck = canon(converse(arcs), 4, perms4)
    print(f"{label[k]:16s} {scores(arcs,4)}  {24//aut:7d}  {aut:5d}  {three_cycles(arcs,4):8d}  {ham_paths(arcs,4):9d}  {str(cycle_lengths(arcs,4)):13s}  {label[ck] if ck != k else 'self'}")
    check(ham_paths(arcs, 4) % 2 == 1, "Redei parity")

print("\nWeighted arc-reversal metagraph (row = class of a representative, entry = # of its 6 arcs")
print("whose reversal lands in the column class):")
W = {}
for k in keys4:
    arcs = cl4[k]
    row = Counter()
    for a in arcs:
        new = frozenset(arcs - {a} | {(a[1], a[0])})
        row[label[canon(new, 4, perms4)]] += 1
    W[label[k]] = row
order = ["TT4", "A=C3+source", "B=C3+sink", "S4 (strong)"]
print(" " * 16 + "".join(f"{o:>14s}" for o in order))
for r in order:
    print(f"{r:16s}" + "".join(f"{W[r][c]:14d}" for c in order))
# labeled consistency: labeled(r)*W[r][c] == labeled(c)*W[c][r]
lab = {"TT4": 24, "A=C3+source": 8, "B=C3+sink": 8, "S4 (strong)": 24}
for r in order:
    for c in order:
        check(lab[r] * W[r][c] == lab[c] * W[c][r], f"symmetric labeled flow {r}->{c}")
print("Labeled-flow symmetry labeled(r)*W[r][c] = labeled(c)*W[c][r] holds for all pairs.")
print("\nConverse quotient (merge A and B into M):")
Q = {"TT4": Counter(), "M": Counter(), "S4": Counter()}
m = {"TT4": "TT4", "A=C3+source": "M", "B=C3+sink": "M", "S4 (strong)": "S4"}
for r in order:
    for c in order:
        if r in ("TT4", "A=C3+source", "S4 (strong)"):  # one representative per quotient node
            Q[m[r]][m[c]] += W[r][c]
for r in ["TT4", "M", "S4"]:
    print(f"  {r:4s} -> " + ", ".join(f"{c}:{Q[r][c]}" for c in ["TT4", "M", "S4"]))
print("  Quotient is a TRIANGLE (every pair adjacent) with loops at TT4 and S4 and no loop at M.")

print("\nEdge-rotation graph of the 3 trees on 5 vertices (remove one edge, add one edge, stay a tree):")
name5 = {}
for T in trees5:
    name5[nx.weisfeiler_lehman_graph_hash(T)] = tree_name(T)
R = defaultdict(Counter)
for T in trees5:
    src = name5[nx.weisfeiler_lehman_graph_hash(T)]
    for e in list(T.edges):
        H = T.copy(); H.remove_edge(*e)
        for u, v in itertools.combinations(T.nodes, 2):
            if H.has_edge(u, v) or (u, v) == e or (v, u) == e:
                continue
            H2 = H.copy(); H2.add_edge(u, v)
            if nx.is_tree(H2):
                R[src][name5[nx.weisfeiler_lehman_graph_hash(H2)]] += 1
for src in R:
    print(f"  {src:22s} -> " + ", ".join(f"{t}:{c}" for t, c in sorted(R[src].items())))
check(R["path P_5"]["star K_{1,4}"] == 0 and R["star K_{1,4}"]["path P_5"] == 0, "path-star nonadjacent")
print("  Rotation graph is a PATH (path -- fork -- star) with loops at all three nodes: NOT the triangle above.")

print("\n" + "=" * 78)
print("D. Erdos--Sos at t = 5: ex(n, T) for the three trees, n <= 7 (exhaustive over the graph atlas)")
print("=" * 78)
atlas = graph_atlas_g()
def contains(G, T):
    GM = nx.algorithms.isomorphism.GraphMatcher(G, T)
    return GM.subgraph_is_monomorphic()
print("n   3n/2    floor  | ex(n,star)  ex(n,fork)  ex(n,path)  | 3-regular graph on n exists?  extremal graphs attaining floor(3n/2)")
for n in range(5, 8):
    graphs_n = [G for G in atlas if G.number_of_nodes() == n]
    ex = {}
    witnesses = {}
    for T in trees5:
        nm = tree_name(T)
        best = -1; wit = None
        for G in graphs_n:
            if G.number_of_edges() > best and not contains(G, T):
                best = G.number_of_edges(); wit = G
        ex[nm] = best; witnesses[nm] = wit
    bound = 3 * n / 2
    fl = math.floor(bound)
    for nm in ex:
        check(ex[nm] <= fl, f"ES bound at n={n} for {nm}")
    cubic = (3 * n) % 2 == 0 and n >= 4
    att = [nm for nm in ex if ex[nm] == fl]
    print(f"{n}   {bound:5.1f}   {fl:3d}   |  {ex['star K_{1,4}']:8d}   {ex[[k for k in ex if k.startswith('spider')][0]]:8d}   {ex['path P_5']:9d}   |  {'yes' if cubic else 'NO (handshake: 3n odd)':22s}  attained by: {att}")
    for nm in ex:
        w = witnesses[nm]
        degs = sorted((d for _, d in w.degree()), reverse=True)
        comps = sorted((len(c) for c in nx.connected_components(w)), reverse=True)
        print(f"      witness for {nm:22s}: {w.number_of_edges()} edges, degrees {degs}, components {comps}")
print("Erdos--Sos for all trees on 5 vertices is a THEOREM (paths: Erdos--Gallai 1959; stars: degree count;")
print("fork = diameter-3 spider: McLennan 2005 / spider results).  The bound 3n/2 is attained with equality")
print("only when 3n is even; the star attains floor(3n/2) for every n >= 5 (near-regular graphs of max degree 3);")
print("the path and fork attain 3n/2 only at 4 | n (disjoint K4's) in the range checked.")

print("\n" + "=" * 78)
print("E. The forcing threshold E*(n) = floor(3n/2)+1 as a Collatz composite (all n <= 10^6)")
print("=" * 78)
for n in range(1, 10 ** 6 + 1):
    E = (3 * n) // 2 + 1
    if n % 2 == 1:
        check(E == (3 * n + 1) // 2 and (3 * n + 1) % 2 == 0, "odd: (3n+1)/2")
    else:
        check(E == 3 * (n // 2) + 1, "even: 3(n/2)+1")
print("odd n : E*(n) = (3n+1)/2  = halve after 3x+1  (the Collatz odd step T(n))")
print("even n: E*(n) = 3(n/2)+1  = 3x+1 after halving (3x+1 applied to T(n) = n/2)")
print("Both are the two orders of composing {x -> 3x+1, x -> x/2}; this is shared syntax ((t-2)/2 = 3/2")
print("and rounding), not a transfer of structure.  Typed: ELEMENTARY IDENTITY / SYNTAX-LEVEL.")
print("\nmod-4 modes of 3n/2 at t = 5:")
for n in range(4, 13):
    E = (3 * n) // 2 + 1
    mode = {0: "0 mod 4: whole, (t-1)|n, disjoint K4 extremal for ALL three trees",
            1: "odd    : half, no 3-regular graph (handshake), star alone reaches the floor",
            2: "2 mod 4: whole, 3-regular graphs exist but (t-1) does not divide n: star-only equality",
            3: "odd    : half, no 3-regular graph (handshake), star alone reaches the floor"}[n % 4]
    print(f"  n={n:2d}  3n/2={3*n/2:5.1f}  E*={E:3d}  {mode}")

print(f"\nALL {CHECKS} CHECKS PASSED")

print("\n" + "=" * 78)
print("F. Otter's recursion: the parity-alternating half term (n <= 10)")
print("=" * 78)
# rooted trees r_n (A000081) via the Euler-transform recursion
N = 10
r = [0] * (N + 1); r[1] = 1
for n in range(2, N + 1):
    s = 0
    for k in range(1, n):
        s += sum(d * r[d] for d in range(1, k + 1) if k % d == 0) * r[n - k]
    r[n] = s // (n - 1)
A000081 = [0, 1, 1, 2, 4, 9, 20, 48, 115, 286, 719]
check(r == A000081, "A000081")
print("n  r_n   (1/2)sum r_i r_{n-i}   [n even] r_{n/2}/2   t_n = r_n - sum/2 + half-term   (A000055)")
for n in range(1, N + 1):
    s = sum(r[i] * r[n - i] for i in range(1, n))
    half = r[n // 2] if n % 2 == 0 else 0
    t2 = 2 * r[n] - s + half          # twice t_n, to stay in integers
    check(t2 % 2 == 0 and t2 // 2 == tree_counts[n], f"Otter at n={n}")
    print(f"{n:2d} {r[n]:4d}   {s/2:8.1f}              {half/2:5.1f}            {t2//2:4d}")
print("The half term r_{n/2}/2 is present exactly for even n (a bicentroid = symmetric central edge exists")
print("only when the two halves have n/2 vertices each).  Otter 1948, CITED; values VERIFIED n <= 10.")

print("\nCenter / centroid of the three trees on 5 vertices:")
for T in trees5:
    center = nx.center(T)
    # centroid: vertices minimizing the largest branch
    def maxbranch(v):
        H = T.copy(); H.remove_node(v)
        return max((len(c) for c in nx.connected_components(H)), default=0)
    mb = {v: maxbranch(v) for v in T.nodes}
    centroid = [v for v in T.nodes if mb[v] == min(mb.values())]
    print(f"  {tree_name(T):22s} diameter={nx.diameter(T)}  center size={len(center)} ({'bicentral' if len(center)==2 else 'central'})  centroid size={len(centroid)}")

print("\n" + "=" * 78)
print("G. Regular tournaments and the mean score (m-1)/2 (m <= 7)")
print("=" * 78)
for m in range(3, 8):
    mean = (m - 1) / 2
    note = "whole (regular tournaments exist)" if m % 2 == 1 else "half (no regular tournament; near-regular only)"
    print(f"  m={m}: mean score {mean:4.1f} = e(K_m)/m = Erdos--Sos density (t-2)/2 for t = m+1 = {m+1}:  {note}")
# near-regular 4-tournament is unique: S4
nr = [label[k] for k in keys4 if scores(cl4[k], 4) == (1, 1, 2, 2)]
check(nr == ["S4 (strong)"], "unique near-regular 4-tournament")
print("  The unique near-regular tournament on 4 vertices is S4 (scores 1,1,2,2); the unique")
print("  acyclic one is TT4; the converse pair A/B sits between them.")

print("\n" + "=" * 78)
print("H. Linear-order dictionary TT_m <-> P_{m+1}: cycles created by one reversal / one added edge")
print("=" * 78)
def cycles_through_lengths(arcs, n):
    out = Counter()
    for k in range(3, n + 1):
        for sub in itertools.combinations(range(n), k):
            first = sub[0]
            for rest in itertools.permutations(sub[1:]):
                cyc = (first,) + rest
                if all((cyc[i], cyc[(i+1) % k]) in arcs for i in range(k)):
                    out[k] += 1
    return out
for m in range(3, 7):
    TT = frozenset((i, j) for i in range(m) for j in range(i + 1, m))
    P = nx.path_graph(m + 1)
    rows = []
    for i in range(m):
        for j in range(i + 1, m):
            new = frozenset(TT - {(i, j)} | {(j, i)})
            lens = cycles_through_lengths(new, m)
            d = j - i
            expected = {L: 1 for L in range(3, d + 2)}  # one cycle of each length 3..d+1? check below
            created = sorted(lens.keys())
            check(created == list(range(3, d + 2)), f"TT_{m} reversal ({i},{j}) cycle lengths")
            # number of cycles of length L through the reversed arc: C(d-1, L-2)
            check(all(lens[L] == math.comb(d - 1, L - 2) for L in lens), f"TT_{m} reversal ({i},{j}) cycle counts")
            # path: adding edge (i, j+1)? use P_{m+1} with vertices 0..m, distance d' = d
    # path check
    for u in range(m + 1):
        for v in range(u + 2, m + 1):
            H = P.copy(); H.add_edge(u, v)
            cyc = nx.cycle_basis(H)
            check(len(cyc) == 1 and len(cyc[0]) == v - u + 1, f"P_{m+1} added edge ({u},{v})")
    print(f"  m={m}: reversing arc (i,j) of TT_{m} (distance d=j-i) creates exactly C(d-1, L-2) cycles of each")
    print(f"        length L = 3..d+1 (total 2^(d-1)-1, zero for d=1); adding edge (u,v) to P_{m+1} (distance d)")
    print(f"        creates exactly one cycle, of length d+1.  VERIFIED all arcs/non-edges.")

print("\n" + "=" * 78)
print("I. Collatz inverse step integrality: period-2 rhythm from ord_3(2) = 2 (m odd, m <= 10^5)")
print("=" * 78)
bad = 0
for mm in range(1, 10 ** 5, 2):
    adm = [k for k in range(1, 13) if (2 ** k * mm - 1) % 3 == 0]
    if mm % 3 == 0:
        check(adm == [], "multiples of 3 have no odd preimage")
    elif mm % 3 == 1:
        check(adm == [2, 4, 6, 8, 10, 12], "m = 1 mod 3: even k")
    else:
        check(adm == [1, 3, 5, 7, 9, 11], "m = 2 mod 3: odd k")
print("  n = (2^k m - 1)/3 is whole iff 2^k m = 1 (mod 3): k even for m = 1 (mod 3), k odd for m = 2 (mod 3),")
print("  never for 3 | m.  Period 2 = ord_3(2).  The Erdos--Sos alternation is 2 | 3n, i.e. 2 | n: period 2")
print("  from the prime 2 itself.  Same period, different mechanism (multiplicative order vs divisibility).")

print(f"\nALL {CHECKS} CHECKS PASSED (final)")

print("\n" + "=" * 78)
print("J. The Klein-group dictionary: end-folds of P_5 vs boundary arc flips of TT_4")
print("=" * 78)
# P_5 = 1-2-3-4-5.  f_L: reattach leaf 1 from 2 to 3.  f_R: reattach leaf 5 from 4 to 3.
def fold(T, leaf, old, new):
    H = T.copy(); H.remove_edge(leaf, old); H.add_edge(leaf, new); return H
P5 = nx.path_graph([1, 2, 3, 4, 5])
fL = lambda T: fold(T, 1, 2, 3) if T.has_edge(1, 2) else fold(T, 1, 3, 2)
fR = lambda T: fold(T, 5, 4, 3) if T.has_edge(5, 4) else fold(T, 5, 3, 4)
orbit_trees = {"id": P5, "L": fL(P5), "R": fR(P5), "LR": fL(fR(P5))}
check(nx.is_isomorphic(fR(fL(P5)), fL(fR(P5))), "folds commute")
check(nx.is_isomorphic(fL(fL(P5)), P5) and nx.is_isomorphic(fR(fR(P5)), P5), "folds are involutions")
for k, T in orbit_trees.items():
    check(nx.is_tree(T), f"{k} is a tree")
    print(f"  P5 after {k:2s}: {tree_name(T):22s} diameter {nx.diameter(T)}  (u = 4 - diameter = {4 - nx.diameter(T)})")
check(nx.is_isomorphic(orbit_trees["L"], orbit_trees["R"]), "fork_L ~ fork_R (reflection is an automorphism of P5)")
rho = {i: 6 - i for i in range(1, 6)}
check(nx.is_isomorphic(nx.relabel_nodes(P5, rho), P5) and set(map(frozenset, nx.relabel_nodes(P5, rho).edges)) == set(map(frozenset, P5.edges)), "rho is an automorphism of P5")
check(set(map(frozenset, nx.relabel_nodes(orbit_trees["L"], rho).edges)) == set(map(frozenset, orbit_trees["R"].edges)), "rho conjugates f_L to f_R")

# TT_4 on 1<2<3<4: i -> j for i < j.  g_L: reverse (1,3) [kills the source]; g_R: reverse (2,4) [kills the sink].
TT4 = frozenset((i, j) for i in range(1, 5) for j in range(i + 1, 5))
def flip(arcs, a):
    return frozenset(arcs - {a} | {(a[1], a[0])}) if a in arcs else frozenset(arcs - {(a[1], a[0])} | {a})
gL = lambda A: flip(A, (1, 3)); gR = lambda A: flip(A, (2, 4))
def sc(A): return tuple(sorted(Counter(x for x, _ in A).get(v, 0) for v in range(1, 5)))
orbit_tours = {"id": TT4, "L": gL(TT4), "R": gR(TT4), "LR": gL(gR(TT4))}
check(gR(gL(TT4)) == gL(gR(TT4)) and gL(gL(TT4)) == TT4, "flips commute / involutions")
names = {(0, 1, 2, 3): "TT4 (u=0)", (0, 2, 2, 2): "3-cycle over sink (u=1)", (1, 1, 1, 3): "source over 3-cycle (u=1)", (1, 1, 2, 2): "strong S4 (u=2)"}
for k, A in orbit_tours.items():
    print(f"  TT4 after {k:2s}: scores {sc(A)}  {names[sc(A)]}")
check(sorted(sc(A) for A in orbit_tours.values()) == sorted(names), "the four flips exhaust the four classes")
rhoT = {i: 5 - i for i in range(1, 5)}
img = frozenset((rhoT[a], rhoT[b]) for a, b in TT4)
check(img == converse(TT4), "rho is an ANTI-automorphism of TT4 (isomorphism onto the converse)")
imgL = frozenset((rhoT[a], rhoT[b]) for a, b in orbit_tours["L"])
check(imgL == converse(orbit_tours["R"]), "rho carries g_L(TT4) onto the converse of g_R(TT4)")
print("  rho(i) = 5-i is an automorphism of the path but an anti-automorphism of TT4: it identifies")
print("  fork_L with fork_R, and identifies the two diamonds only after converse-merging.  Merged label")
print("  u = #boundary defects (THM-584) <-> 4 - diameter;  H = 2u+1 = 1,3,5;  leaves = u+2 = 2,3,4.")

print("\n" + "=" * 78)
print("K. Fork-free structure lemma and exact ex(n, fork) = ex(n, P5) = 6 floor(n/4) + C(n mod 4, 2)")
print("=" * 78)
fork = [T for T in trees5 if tree_name(T).startswith("spider")][0]
def is_fork_free_by_structure(G):
    # claim: fork-free iff every component is a star, has max degree <= 2, or has <= 4 vertices
    for comp in nx.connected_components(G):
        H = G.subgraph(comp)
        if H.number_of_nodes() <= 4:
            continue
        degs = [d for _, d in H.degree()]
        if max(degs) <= 2:
            continue
        if max(degs) == H.number_of_nodes() - 1 and H.number_of_edges() == H.number_of_nodes() - 1:
            continue  # star
        return False
    return True
cnt = 0
for G in atlas:
    if G.number_of_nodes() == 0:
        continue
    check(is_fork_free_by_structure(G) == (not contains(G, fork)), f"fork-free structure on atlas graph {G.number_of_edges()}e")
    cnt += 1
print(f"  structure lemma agrees with direct subgraph search on all {cnt} atlas graphs (<= 7 vertices)")
import random
random.seed(20261005)
rc = 0
for trial in range(3000):
    n = random.randint(8, 12); p = random.choice([0.15, 0.25, 0.35, 0.5])
    G = nx.gnp_random_graph(n, p, seed=random.randint(0, 10**9))
    check(is_fork_free_by_structure(G) == (not contains(G, fork)), "fork-free structure random")
    rc += 1
print(f"  and on {rc} random graphs with 8..12 vertices")
print("  Proof: a vertex c of degree >= 3 whose neighbour a has a neighbour outside N[c] gives a fork; so N[c] is a")
print("  component; if deg c >= 4 any edge inside N(c) gives a fork, so the component is a star; if deg c = 3 the")
print("  component has 4 vertices.  Hence fork-free components are stars, paths, cycles or <= 4 vertices, and the")
print("  maximum edge count is 6 floor(n/4) + C(r,2): disjoint K4's plus a clique on the remainder (Erdos--Gallai/")
print("  Faudree--Schelp give the same value for P5).")
print("\n  n   3n/2   ex(star)  ex(path)=ex(fork)   deficit star   deficit path/fork")
for n in range(4, 17):
    exs = (3 * n) // 2
    exp = 6 * (n // 4) + math.comb(n % 4, 2)
    print(f"  {n:2d}  {3*n/2:5.1f}   {exs:5d}     {exp:8d}            {3*n/2-exs:4.1f}           {3*n/2-exp:4.1f}")
print("  deficit(star) has period 2 (handshake parity); deficit(path/fork) has period 4 (K4 packing residue).")

print(f"\nALL {CHECKS} CHECKS PASSED (final, with J and K)")

print("\n" + "=" * 78)
print("L. The end-defect Klein group for every m >= 4: TT_m with flips (1,3),(m-2,m); P_{m+1} with end folds")
print("=" * 78)
def canon_small(arcs, m):
    best = None
    for p in itertools.permutations(range(1, m + 1)):
        img = tuple(sorted((p[a-1], p[b-1]) for a, b in arcs))
        if best is None or img < best:
            best = img
    return best
print(" m  classes TT,B,A,S distinct?  A=conv(B)?  TT,S self-conv?  H(TT),H(B),H(A),H(S)   trees: path,forkL~forkR,double-fold  diameters  u=m-diam")
for m in range(4, 8):
    TTm = frozenset((i, j) for i in range(1, m + 1) for j in range(i + 1, m + 1))
    Bm = flip(TTm, (1, 3)); Am = flip(TTm, (m - 2, m)); Sm = flip(Bm, (m - 2, m))
    cls = [canon_small(X, m) for X in (TTm, Bm, Am, Sm)]
    distinct = len(set(cls)) == 4
    conv_ok = canon_small(converse(Bm), m) == cls[2]
    self_conv = canon_small(converse(TTm), m) == cls[0] and canon_small(converse(Sm), m) == cls[3]
    def hp(A):
        return sum(1 for p in itertools.permutations(range(1, m + 1)) if all((p[i], p[i+1]) in A for i in range(m - 1)))
    Hs = [hp(X) for X in (TTm, Bm, Am, Sm)]
    check(distinct and conv_ok and self_conv, f"Klein classes at m={m}")
    check(all(h % 2 == 1 for h in Hs), "Redei")
    P = nx.path_graph(range(1, m + 2))
    fl = fold(P, 1, 2, 3); fr = fold(P, m + 1, m, m - 1); both = fold(fl, m + 1, m, m - 1)
    check(nx.is_isomorphic(fl, fr) and not nx.is_isomorphic(fl, P) and not nx.is_isomorphic(both, fl) and not nx.is_isomorphic(both, P), f"tree orbits at m={m}")
    diams = [nx.diameter(X) for X in (P, fl, both)]
    check(diams == [m, m - 1, m - 2], "u = m - diameter")
    print(f" {m}   {str(distinct):5s}                   {str(conv_ok):5s}       {str(self_conv):5s}           {Hs}           3 trees (one pair merged by rho)    {diams}    [0,1,2]")
print(" For every m the end-defect group (Z_2)^2 gives 4 tournament classes (one converse pair) and 3 trees; the")
print(" merged label u = #end defects = m - diameter.  At m = 4 these exhaust all classes/trees; for m >= 5 they")
print(" are the boundary layer only (10 merged classes vs 6 trees at m = 5).  H(B_m) = H(A_m) = 3 for all m;")
print(" H(S_m) = 5 at m = 4 (the two flips interlock) and 9 for every m >= 5 (verified m <= 7).")

print("\nSpectral obstruction to any map between the merged 4-metagraph and the three-row automaton:")
import numpy as np
Qm = np.array([[3, 2, 1], [3, 0, 3], [1, 2, 3]], dtype=float)   # rows TT4, M, S4 (section C)
rot = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]], dtype=float)  # rows 1 -> 5 -> 3 -> 1 (collatz_mod6 Theorem 1.1)
ev_q = sorted(float(x) for x in np.round(np.linalg.eigvals(Qm).real, 6)); ev_r = [complex(round(z.real, 3), round(z.imag, 3)) for z in np.linalg.eigvals(rot)]
print("  merged metagraph (weighted, as above) eigenvalues:", ev_q, " (real; THM-584 V+ = {-2,2,6})")
print("  three-row rotation eigenvalues:", sorted(ev_r, key=lambda z: (z.real, z.imag)), " (cube roots of unity)")
check(ev_q == [-2.0, 2.0, 6.0], "THM-584 V+ spectrum reproduced")
print("  A symmetric 3-node graph with two loops cannot be carried onto a directed 3-cycle: shared datum is only '3'.")

print(f"\nALL {CHECKS} CHECKS PASSED (final, with L)")
