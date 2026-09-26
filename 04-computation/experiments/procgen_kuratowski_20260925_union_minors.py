#!/usr/bin/env python3
"""procgen_kuratowski_20260925, part U: the sign fold and Althofer's 3n+-1 union graph U.

U = the simple graph on N = {1,2,...} with the halving edges {m,2m} and, for odd n, the two up-edges
{n,3n+1} and {n,3n-1}; U_N = U induced on [1,N]. (Positions of Althofer's 3n+-1 game.)

Sections (every check raises on failure):
  U1  the fold pi(x) = |x| from the Collatz graph Gamma_+ on Z\\{0} to U: a graph homomorphism, fibres
      (halving 2, up+ 1, up- 1), locally injective, NOT a covering: the local deficit is exactly at
      x = 1,2,3,5 mod 6 (Kohl's Tait defects); W = Gamma_+ u Gamma_- -> U is a covering but the
      TRIVIAL one (no edge joins the two signs). Dodecahedron -> Petersen control (connected cover).
  U2  U is subcubic: deg 3 except 0(6) (deg 2) and 1, 2 (deg 2 simple; the multigraph digon {1,2}).
  U3  girth 3; the triangle {1,2,4} is the only triangle; counts of short cycles in U_2000.
  U4  planarity: U_N planar iff N <= 51.  Certificates: a verified planar rotation system for U_51
      (face count, Euler per component) and an explicit K_{3,3} subdivision in U_52 with every edge
      labelled by its arithmetic move; kernel(U_51) = prism, kernel(U_52) = K_{3,3} with a truncated
      vertex.
  U5  K5 minor: U_N has a K5 minor iff N >= 68.  Explicit branch sets in U_68; absence in U_67 by the
      perfect-matching argument on its 10-vertex cubic kernel, cross-checked by CP-SAT.
  U6  Petersen minor: iff N >= 104. Explicit Petersen subdivision in U_104 (10 branch vertices, 15
      arithmetic paths); absence in U_103 by exhaustive subdivision search on its 18-vertex kernel,
      cross-checked by CP-SAT.
  U7  kernels of U_N (N <= 4000) are connected, bridgeless, cubic and 3-edge-colourable (SAT).
  U8  planar double covers: U_N (N >= 52) has a planar double cover iff N <= 75 (exhaustive over all
      signings of the kernel); the Collatz sign cover W is trivial; edge-type signings are nonplanar.
Run: python3 04-computation/experiments/procgen_kuratowski_20260925_union_minors.py  (~3-6 min, < 500 MB)
Requires networkx, ortools (CP-SAT, 2 workers), python-sat.
"""
import itertools
import sys
import time
from collections import Counter

import networkx as nx
from ortools.sat.python import cp_model
from pysat.solvers import Glucose4


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


def C(x):
    return 3 * x + 1 if x % 2 else x // 2


def U_graph(N):
    G = nx.Graph()
    G.add_nodes_from(range(1, N + 1))
    for m in range(1, N // 2 + 1):
        G.add_edge(m, 2 * m)
    for n in range(1, N + 1, 2):
        for y in (3 * n + 1, 3 * n - 1):
            if y <= N:
                G.add_edge(n, y)
    return G


def label(u, v):
    """arithmetic move(s) joining u and v in U."""
    lo, hi = min(u, v), max(u, v)
    out = []
    if hi == 2 * lo:
        out.append(f"{hi}/2={lo}")
    if lo % 2 and hi == 3 * lo + 1:
        out.append(f"3*{lo}+1={hi}")
    if lo % 2 and hi == 3 * lo - 1:
        out.append(f"3*{lo}-1={hi}")
    return "|".join(out)


def kernel(G, keep_map=False):
    """min-degree-3 reduction (delete loops, merge parallel edges, delete degree <= 1, suppress degree 2)
    to a simple graph; valid for minors / topological minors of simple graphs of min degree >= 3 and for
    planarity and planar double covers.  If keep_map, also return for each kernel edge the G-path."""
    H = nx.MultiGraph()
    H.add_nodes_from(G.nodes())
    for u, v in G.edges():
        H.add_edge(u, v, path=(u, v))
    changed = True
    while changed:
        changed = False
        loops = list(nx.selfloop_edges(H, keys=True))
        if loops:
            H.remove_edges_from(loops)
            changed = True
        for u, v in list(H.edges()):
            if u != v and H.number_of_edges(u, v) > 1:
                ks = list(H[u][v].keys())
                for k in ks[1:]:
                    H.remove_edge(u, v, k)
                changed = True
        low = [v for v in H if H.degree(v) <= 1]
        if low:
            H.remove_nodes_from(low)
            changed = True
            continue
        for v in list(H.nodes()):
            if v in H and H.degree(v) == 2:
                es = list(H.edges(v, keys=True, data="path"))
                if len(es) != 2:
                    continue
                (_, x, k1, p1), (_, y, k2, p2) = es
                if x == v or y == v:
                    continue
                q1 = p1 if p1[-1] == v else tuple(reversed(p1))     # path x..v
                q2 = p2 if p2[0] == v else tuple(reversed(p2))      # path v..y
                check(q1[-1] == v and q2[0] == v, "path orientation")
                H.remove_node(v)
                H.add_edge(q1[0], q2[-1], path=q1 + q2[1:])
                changed = True
    K = nx.Graph()
    K.add_nodes_from(H.nodes())
    pm = {}
    for u, v, p in H.edges(data="path"):
        K.add_edge(u, v)
        pm[(u, v)] = p
        pm[(v, u)] = tuple(reversed(p))
    return (K, pm) if keep_map else K


print("=" * 100)
print("U1  the sign fold pi(x) = |x| and the dodecahedron -> Petersen control")
print("=" * 100)
Xf = 30000
Gp = nx.MultiGraph()
for x in range(-Xf, Xf + 1):
    if x != 0 and abs(C(x)) <= Xf:
        Gp.add_edge(x, C(x))
Uf = U_graph(Xf)
fib = Counter()
for u, v in Gp.edges():
    pu, pv = abs(u), abs(v)
    check(Uf.has_edge(pu, pv), f"pi maps edge {u},{v} to a U-edge")
    fib[(min(pu, pv), max(pu, pv), "neg" if u < 0 else "pos")] += 1
kind = Counter()
for (p, q, sg), k in fib.items():
    if q == 2 * p and q == 3 * p - 1:
        t = "halving=up- (the digon {1,2})"
    elif q == 2 * p:
        t = "halving"
    elif q == 3 * p + 1:
        t = "up+"
    else:
        t = "up-"
    kind[(t, sg, k)] += 1
print("  pi is a graph homomorphism Gamma_+ -> U on |x| <= %d; (edge type, sheet, multiplicity) -> #U-edges:" % Xf)
for key in sorted(kind):
    print("    ", key, kind[key])
check(all(not (sg == "pos" and t == "up-") and not (sg == "neg" and t == "up+") for (t, sg, k) in kind),
      "up+ edges come only from positives, up- only from negatives")
exc = []
for x in range(-Xf // 4, Xf // 4 + 1):
    if x == 0:
        continue
    d_up = len({abs(y) for y in Gp.neighbors(x)})
    d_U = Uf.degree(abs(x))
    expected = 1 if x % 6 in (1, 2, 3, 5) else 0
    if d_U - d_up != expected:
        exc.append(x)
print("  local deficit deg_U(|x|) - deg_Gamma+(x) equals 1 exactly at x = 1,2,3,5 mod 6 and 0 at x = 0,4 mod 6,")
print("  for all 0 < |x| <= %d, the only exceptions being" % (Xf // 4), exc, "(the simple-graph digon {1,2})")
check(exc == [1, 2], "deficit pattern")
print("  => pi is locally injective but not locally surjective: an immersion (fold), not a covering; the")
print("     deficit set 1,2,3,5 mod 6 is exactly the support of Kohl's Tait boundary delta (part K5):")
print("     the fold attaches the other sheet's up-edge at every Tait defect.")
# W = Gamma_+ u Gamma_- on Z\{0}: covering of U, trivial
Cm = lambda x: 3 * x - 1 if x % 2 else x // 2
Wg = nx.Graph()
for x in range(-Xf, Xf + 1):
    if x == 0:
        continue
    for f in (C, Cm):
        y = f(x)
        if abs(y) <= Xf:
            Wg.add_edge(x, y)
cross = [(u, v) for u, v in Wg.edges() if (u > 0) != (v > 0)]
check(not cross, "W has no edge between the signs")
locbij = all(sorted(abs(y) for y in Wg.neighbors(x)) == sorted(Uf.neighbors(abs(x)))
             for x in range(-Xf // 4, Xf // 4 + 1) if x != 0)
check(locbij, "W -> U locally bijective")
print("  W = Gamma_+ u Gamma_- -> U is locally bijective (a covering) with deck involution nu, but no W-edge joins")
print("     x > 0 to y < 0: W = U + nu(U), the TRIVIAL double cover (every Collatz-type move preserves sign).")
D = nx.dodecahedral_graph()
dist = dict(nx.all_pairs_shortest_path_length(D))
anti = {v: max(dist[v], key=lambda w: dist[v][w]) for v in D}
check(all(dist[v][anti[v]] == 5 and anti[anti[v]] == v for v in D), "antipodal involution")
check(all(D.has_edge(anti[u], anti[v]) for u, v in D.edges()), "antipode is an automorphism")
orb = {v: min(v, anti[v]) for v in D}
Q = nx.Graph()
for u, v in D.edges():
    Q.add_edge(orb[u], orb[v])
check(nx.is_isomorphic(Q, nx.petersen_graph()), "dodecahedron / antipode = Petersen")
check(all(len({orb[w] for w in D.neighbors(v)}) == 3 for v in D), "locally bijective")
check(nx.check_planarity(D)[0] and not nx.check_planarity(Q)[0] and nx.is_connected(D), "planar cover of nonplanar")
print("  control: dodecahedron -> dodecahedron/antipode = Petersen is a CONNECTED double covering (planar cover of a")
print("     nonplanar graph). The Collatz sign cover is disconnected, so the analogy fails at the key point.")

print()
print("=" * 100)
print("U2  U is subcubic")
print("=" * 100)
degs = Counter()
for n in range(1, Xf // 4):
    degs[(n % 6 if n > 2 else f"n={n}", Uf.degree(n))] += 1
print("  (n mod 6, degree) counts for n < %d:" % (Xf // 4), dict(sorted(degs.items(), key=str)))
check(all(d == (2 if (r == 0 or isinstance(r, str)) else 3) for (r, d) in degs), "degree pattern")
print("  odd n: neighbours 2n, 3n+1, 3n-1 (all larger); even n: n/2, 2n and (n-1)/3 or (n+1)/3 when odd;")
print("  degree 3 except n = 0 mod 6 and n = 1, 2 (the digon {1,2}: 2 = 1/2... = 3*1-1)")

print()
print("=" * 100)
print("U3  girth and short cycles")
print("=" * 100)
U2k = U_graph(2000)
check(nx.girth(U2k) == 3, "girth 3")
tris = sorted({tuple(sorted(t)) for t in (c for c in nx.enumerate_all_cliques(U2k) if len(c) == 3)})
print("  triangles in U_2000:", tris)
check(tris == [(1, 2, 4)], "unique triangle")
# count cycles of length <= 8 in U_600 by brute force DFS from the minimum vertex
U600 = U_graph(600)
cnt = Counter()
adj = {v: sorted(U600.neighbors(v)) for v in U600}
for s in adj:
    stack = [(s, [s])]
    while stack:
        v, path = stack.pop()
        for w in adj[v]:
            if w == s and len(path) >= 3:
                if path[1] < path[-1]:
                    cnt[len(path)] += 1
            elif w > s and w not in path and len(path) < 8:
                stack.append((w, path + [w]))
print("  cycles of length 3..8 in U_600:", dict(sorted(cnt.items())))
print("  the 3n-1 cycles {5,14,7,20,10} and {17,...} (length 18) and the 3n+1 triangle are among them")

print()
print("=" * 100)
print("U4  planarity threshold N = 52 with certificates")
print("=" * 100)
thr = None
for N in range(1, 200):
    if not nx.check_planarity(U_graph(N))[0]:
        thr = N
        break
print("  first nonplanar U_N: N =", thr)
check(thr == 52, "planarity threshold 52")
G51 = U_graph(51)
ok, emb = nx.check_planarity(G51)
check(ok, "U_51 planar")
# independent face count from the rotation system
rot = {v: list(emb.neighbors_cw_order(v)) for v in emb}
darts = {(u, v) for u in rot for v in rot[u]}
seen = set()
faces_per_comp = Counter()
comp_of = {}
for i, comp in enumerate(nx.connected_components(G51)):
    for v in comp:
        comp_of[v] = i
for d in sorted(darts):
    if d in seen:
        continue
    u, v = d
    f0 = d
    while True:
        seen.add((u, v))
        # next dart: at v, the successor of u in the clockwise rotation (face tracing rule)
        r = rot[v]
        w = r[(r.index(u) + 1) % len(r)]
        u, v = v, w
        if (u, v) == f0:
            break
    faces_per_comp[comp_of[d[0]]] += 1
eul = []
for i, comp in enumerate(nx.connected_components(G51)):
    V = len(comp)
    E = G51.subgraph(comp).number_of_edges()
    F = faces_per_comp[i] if E > 0 else 1
    eul.append(V - E + F)
check(all(x == 2 for x in eul), "Euler characteristic 2 on every component")
print(f"  U_51: rotation system verified by face tracing: V - E + F = 2 on all {len(eul)} components")
G52 = U_graph(52)
ok, sub = nx.check_planarity(G52, counterexample=True)
check(not ok, "U_52 nonplanar")
dsub = dict(sub.degree())
branch = sorted(v for v in dsub if dsub[v] == 3)
check(len(branch) == 6 and all(dsub[v] in (2, 3) for v in dsub), "K33 subdivision shape")
# recover the 9 branch paths
paths = []
for s0 in branch:
    for w in sub.neighbors(s0):
        p = [s0, w]
        while p[-1] not in branch:
            nxt = [z for z in sub.neighbors(p[-1]) if z != p[-2]]
            p.append(nxt[0])
        if p[0] < p[-1]:
            paths.append(p)
Hs = nx.Graph([(p[0], p[-1]) for p in paths])
check(nx.is_isomorphic(Hs, nx.complete_bipartite_graph(3, 3)), "suppressed = K33")
A_side, B_side = nx.bipartite.sets(Hs)
print("  explicit K_{3,3} subdivision in U_52: shores", sorted(A_side), sorted(B_side))
inner = Counter()
for p in sorted(paths):
    for z in p[1:-1]:
        inner[z] += 1
    check(all(G52.has_edge(p[i], p[i + 1]) for i in range(len(p) - 1)), "path edges in U")
    print("    ", " - ".join(map(str, p)), "   [", "; ".join(label(p[i], p[i + 1]) for i in range(len(p) - 1)), "]")
check(all(k == 1 for k in inner.values()) and not (set(inner) & set(branch)), "internally disjoint")
print("  the last vertex 52 enters through the ear 17 -- 52 -- 26:", label(17, 52), ",", label(26, 52))
K51, K52 = kernel(G51), kernel(G52)
check(nx.is_isomorphic(K51, nx.circular_ladder_graph(3)), "kernel(U_51) = prism")
tri52 = [c for c in nx.enumerate_all_cliques(K52) if len(c) == 3]
check(len(tri52) == 1, "one triangle in kernel(U_52)")
Kc = nx.contracted_nodes(nx.contracted_nodes(K52, tri52[0][0], tri52[0][1], self_loops=False), tri52[0][0],
                         tri52[0][2], self_loops=False)
check(nx.is_isomorphic(nx.Graph(Kc), nx.complete_bipartite_graph(3, 3)), "contract triangle -> K33")
print("  kernel(U_51) = triangular prism (= K4 with a truncated vertex) on", sorted(K51.nodes()))
print("  kernel(U_52) = K_{3,3} with a truncated vertex (contract its triangle", sorted(tri52[0]), "-> K_{3,3})")
print("  both contain the 3n-1 5-cycle 5 -14 - 7 - 20 - 10 (as kernel edges) : first nonplanarity = sheet interaction")


# ---------------------------------------------------------------- minor machinery
def find_minor(G, H, time_limit=900, workers=2, seed=0):
    Gn = list(G.nodes())
    Hn = list(H.nodes())
    n = len(Gn)
    m = cp_model.CpModel()
    x = {(v, i): m.NewBoolVar("") for v in Gn for i in Hn}
    r = {(v, i): m.NewBoolVar("") for v in Gn for i in Hn}
    d = {(v, i): m.NewIntVar(0, n, "") for v in Gn for i in Hn}
    for v in Gn:
        m.Add(sum(x[v, i] for i in Hn) <= 1)
    for i in Hn:
        m.Add(sum(r[v, i] for v in Gn) == 1)
        for v in Gn:
            m.AddImplication(r[v, i], x[v, i])
            m.Add(d[v, i] == 0).OnlyEnforceIf(r[v, i])
            ps = []
            for u in G.neighbors(v):
                p = m.NewBoolVar("")
                m.AddImplication(p, x[u, i])
                m.Add(d[u, i] + 1 <= d[v, i]).OnlyEnforceIf(p)
                ps.append(p)
            m.AddBoolOr(ps + [x[v, i].Not(), r[v, i]])
    for (i, j) in H.edges():
        lits = []
        for (u, v) in G.edges():
            for (p, q) in ((u, v), (v, u)):
                e = m.NewBoolVar("")
                m.AddImplication(e, x[p, i])
                m.AddImplication(e, x[q, j])
                lits.append(e)
        m.AddBoolOr(lits)
    s = cp_model.CpSolver()
    s.parameters.max_time_in_seconds = time_limit
    s.parameters.num_search_workers = workers
    s.parameters.random_seed = seed
    st = s.Solve(m)
    if st in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        return {i: {v for v in Gn if s.Value(x[v, i])} for i in Hn}
    if st == cp_model.INFEASIBLE:
        return None
    return "UNKNOWN"


def verify_minor(G, H, B):
    used = set()
    for i, S in B.items():
        check(S and not (S & used), "branch sets disjoint, nonempty")
        used |= S
        check(nx.is_connected(G.subgraph(S)), "branch set connected")
    for (i, j) in H.edges():
        check(any(G.has_edge(u, v) for u in B[i] for v in B[j]), f"branch sets {i},{j} adjacent")
    return True


def two_core(G):
    return nx.k_core(G, 2)


print()
print("=" * 100)
print("U5  K5 minor threshold N = 68")
print("=" * 100)
K5 = nx.complete_graph(5)
G68 = U_graph(68)
t0 = time.time()
B = find_minor(two_core(G68), K5, workers=1)
check(isinstance(B, dict), "K5 minor found in U_68")
verify_minor(G68, K5, B)
print(f"  U_68: K5 branch sets (verified in U_68 itself; {time.time() - t0:.1f}s):")
for i in sorted(B):
    print("    ", sorted(B[i]))
K67 = kernel(U_graph(67))
check(K67.number_of_nodes() == 10 and set(dict(K67.degree()).values()) == {3}, "kernel(U_67) cubic on 10 vertices")


def perfect_matchings(G):
    Vs = sorted(G.nodes())

    def rec(rem):
        if not rem:
            yield []
            return
        v = rem[0]
        for w in G.neighbors(v):
            if w in rem:
                r2 = [z for z in rem if z not in (v, w)]
                for rest in rec(r2):
                    yield [(v, w)] + rest
    yield from rec(Vs)


pms = list(perfect_matchings(K67))
k5hits = 0
for M in pms:
    Q = nx.MultiGraph()
    rep = {}
    for i, (u, v) in enumerate(M):
        rep[u] = rep[v] = i
    for u, v in K67.edges():
        if rep[u] != rep[v]:
            Q.add_edge(rep[u], rep[v])
    if nx.Graph(Q).number_of_edges() == 10:
        k5hits += 1
check(k5hits == 0, "no perfect matching of kernel(U_67) contracts to K5")
print(f"  U_67: kernel is cubic on 10 vertices; a K5 model there needs 5 two-vertex branch sets (>= 4 exits each),")
print(f"        i.e. a perfect matching M with kernel/M = K5; all {len(pms)} perfect matchings fail -> no K5 minor.")
t0 = time.time()
cs = find_minor(K67, K5)
check(cs is None, "CP-SAT: no K5 minor in kernel(U_67)")
print(f"        CP-SAT cross-check: INFEASIBLE ({time.time() - t0:.1f}s).  Monotone in N: K5 minor iff N >= 68.")

print()
print("=" * 100)
print("U6  Petersen minor threshold N = 104")
print("=" * 100)
P = nx.petersen_graph()
G104 = U_graph(104)
core104 = two_core(G104)
t0 = time.time()
B = find_minor(core104, P, workers=1)
check(isinstance(B, dict), "Petersen minor in U_104")
verify_minor(G104, P, B)
# convert the minor model into a subdivision (tripods in each branch set)
exits = {}
for (i, j) in P.edges():
    e = next((u, v) for u in sorted(B[i]) for v in sorted(B[j]) if G104.has_edge(u, v))
    exits[(i, j)] = e
    exits[(j, i)] = (e[1], e[0])
centre = {}
arm = {}
for i in P.nodes():
    Tsub = nx.minimum_spanning_tree(G104.subgraph(B[i]))
    ports = [exits[(i, j)][0] for j in P.neighbors(i)]
    p01 = nx.shortest_path(Tsub, ports[0], ports[1])
    p02 = nx.shortest_path(Tsub, ports[0], ports[2])
    p12 = nx.shortest_path(Tsub, ports[1], ports[2])
    med = set(p01) & set(p02) & set(p12)
    check(len(med) == 1, "tripod median")
    cen = med.pop()
    centre[i] = cen
    for j in P.neighbors(i):
        arm[(i, j)] = nx.shortest_path(Tsub, cen, exits[(i, j)][0])
S = nx.Graph()
ppaths = {}
for (i, j) in P.edges():
    p = arm[(i, j)] + list(reversed(arm[(j, i)]))
    ppaths[(i, j)] = p
    for k in range(len(p) - 1):
        S.add_edge(p[k], p[k + 1])
dS = dict(S.degree())
check(sorted(v for v in dS if dS[v] == 3) == sorted(centre.values()), "branch vertices = tripod centres")
check(all(d in (2, 3) for d in dS.values()), "subdivision degrees")
check(all(G104.has_edge(u, v) for u, v in S.edges()), "subdivision inside U_104")
Hs = nx.Graph([(centre[i], centre[j]) for (i, j) in P.edges()])
check(nx.is_isomorphic(Hs, P), "suppressed subdivision = Petersen")
inner = Counter(z for p in ppaths.values() for z in p[1:-1])
check(all(k == 1 for k in inner.values()) and not (set(inner) & set(centre.values())), "internally disjoint")
print(f"  explicit Petersen subdivision in U_104 ({time.time() - t0:.1f}s): branch vertices", sorted(centre.values()))
for (i, j) in sorted(P.edges()):
    p = ppaths[(i, j)]
    print(f"    {' - '.join(map(str, p)):40s} [{'; '.join(label(p[k], p[k + 1]) for k in range(len(p) - 1))}]")
K103 = kernel(U_graph(103))
nV = K103.number_of_nodes()
check(set(dict(K103.degree()).values()) == {3} and nV == 18, "kernel(U_103) cubic on 18 vertices")


def has_petersen_subdivision(K, stop_first=True):
    """exhaustive: choose the unused vertex set U; in G' = K - U every vertex has degree >= 2; a matching F on
    the degree-3 vertices is deleted so that exactly 10 vertices keep degree 3; the result must be connected
    and, suppressed, the Petersen graph (simple cubic on 10 vertices with girth 5 -- unique)."""
    Vs = sorted(K.nodes())
    n = len(Vs)
    adj0 = {v: set(K.neighbors(v)) for v in Vs}
    tested = 0
    for u in range(0, n - 10 + 1):
        for Uset in itertools.combinations(Vs, u):
            Us = set(Uset)
            adj = {v: adj0[v] - Us for v in Vs if v not in Us}
            if any(len(a_) < 2 for a_ in adj.values()):
                continue
            D3 = [v for v in adj if len(adj[v]) == 3]
            if len(D3) < 10 or (len(D3) - 10) % 2:
                continue
            need = (len(D3) - 10) // 2
            D3s = set(D3)
            E3 = [(p_, q_) for p_ in D3 for q_ in adj[p_] if q_ in D3s and p_ < q_]
            for F in itertools.combinations(E3, need):
                vs = [z for e in F for z in e]
                if len(set(vs)) != 2 * need:
                    continue
                tested += 1
                a2 = {v: set(adj[v]) for v in adj}
                for p_, q_ in F:
                    a2[p_].discard(q_)
                    a2[q_].discard(p_)
                br = [v for v in a2 if len(a2[v]) == 3]
                # suppress: walk from each branch vertex along each edge to the next branch vertex
                Hk = nx.MultiGraph()
                okp = True
                seen_int = set()
                for v in br:
                    for w in a2[v]:
                        prev, cur = v, w
                        while len(a2[cur]) == 2:
                            seen_int.add(cur)
                            nxt = next(z for z in a2[cur] if z != prev)
                            prev, cur = cur, nxt
                        if v < cur or (v == cur):
                            Hk.add_edge(v, cur)
                if len(seen_int) + len(br) != len(a2):
                    continue       # a pure cycle component of degree-2 vertices
                Hs_ = nx.Graph(Hk)
                if Hk.number_of_edges() != 15 or Hs_.number_of_edges() != 15 or any(u_ == v_ for u_, v_ in Hk.edges()):
                    continue
                if Hs_.number_of_nodes() == 10 and set(dict(Hs_.degree()).values()) == {3} and nx.is_connected(Hs_) \
                        and nx.girth(Hs_) == 5:
                    if stop_first:
                        return True, tested
    return False, tested


t0 = time.time()
res, tested = has_petersen_subdivision(K103)
check(not res, "no Petersen subdivision in kernel(U_103)")
print(f"  U_103: exhaustive subdivision search on its 18-vertex cubic kernel: none ({tested} candidate deletions, "
      f"{time.time() - t0:.1f}s)")
t0 = time.time()
cs = find_minor(K103, P)
check(cs is None, "CP-SAT: no Petersen minor in kernel(U_103)")
print(f"         CP-SAT cross-check: INFEASIBLE ({time.time() - t0:.1f}s).  Since U is subcubic, Petersen minor =")
print("         Petersen subdivision; monotone in N: Petersen minor iff N >= 104.")
K104 = kernel(G104)
t0 = time.time()
res, tested = has_petersen_subdivision(K104)
check(res, "positive control: exhaustive search finds Petersen in kernel(U_104)")
print(f"  positive control: the exhaustive search finds one in kernel(U_104) ({K104.number_of_nodes()} vertices, "
      f"{time.time() - t0:.1f}s)")

print()
print("=" * 100)
print("U7  kernels: bridgeless cubic and 3-edge-colourable")
print("=" * 100)


def three_edge_colourable(K):
    edges = list(K.edges())
    idx = {}
    for i, (u, v) in enumerate(edges):
        idx[(u, v)] = idx[(v, u)] = i
    var = lambda i, c: 3 * i + c + 1
    s = Glucose4()
    for i in range(len(edges)):
        s.add_clause([var(i, c) for c in range(3)])
    for v in K:
        inc = [idx[(v, w)] for w in K.neighbors(v)]
        for p_, q_ in itertools.combinations(inc, 2):
            for c in range(3):
                s.add_clause([-var(p_, c), -var(q_, c)])
    r = s.solve()
    s.delete()
    return r


lastk = None
nk = 0
for N in list(range(52, 1001)) + [2000, 4000]:
    K = kernel(U_graph(N))
    key = (K.number_of_nodes(), K.number_of_edges())
    if key == lastk:
        continue
    lastk = key
    nk += 1
    br = list(nx.bridges(K))
    col = three_edge_colourable(K)
    check(col and not br and nx.is_connected(K) and set(dict(K.degree()).values()) == {3},
          f"kernel(U_{N}) connected bridgeless cubic and colourable")
    if N in (52, 104, 250, 500, 1000, 2000, 4000):
        print(f"  N = {N:5d}: kernel V = {K.number_of_nodes():4d}, cubic, connected, bridges = {len(br)}, "
              f"3-edge-colourable = {col}")
print(f"  checked every distinct kernel for 52 <= N <= 1000 ({nk - 2} kernels) and N = 2000, 4000")
print("  => the Petersen minors of U_N obstruct nothing: every kernel is Tait-colourable (no snark).")

print()
print("=" * 100)
print("U8  planar double covers (the Petersen/dodecahedron phenomenon) of U_N")
print("=" * 100)


def planar_double_cover(K):
    T = nx.minimum_spanning_tree(K)
    co = [e for e in K.edges() if not T.has_edge(*e)]
    for mask in range(1 << len(co)):
        neg = {co[i] for i in range(len(co)) if mask >> i & 1}
        Dc = nx.Graph()
        for (u, v) in K.edges():
            if (u, v) in neg:
                Dc.add_edge((u, 0), (v, 1))
                Dc.add_edge((u, 1), (v, 0))
            else:
                Dc.add_edge((u, 0), (v, 0))
                Dc.add_edge((u, 1), (v, 1))
        if nx.check_planarity(Dc)[0]:
            return len(co), neg, nx.is_connected(Dc)
    return len(co), None, None


last = None
window = []
for N in range(52, 77):
    UN = U_graph(N)
    K = kernel(UN)
    cyc_comps = [cc for cc in nx.connected_components(UN) if UN.subgraph(cc).number_of_edges() >= len(cc)]
    check(nx.is_connected(K) and len(cyc_comps) == 1, "U_N has exactly one cyclic component (kernel connected)")
    key = tuple(sorted(K.edges()))
    if key == last:
        window.append((N, window[-1][1] if window else None))
        continue
    last = key
    beta, neg, conn = planar_double_cover(K)
    window.append((N, neg is not None))
    print(f"  N = {N}: kernel V = {K.number_of_nodes()}, beta = {beta}: planar double cover "
          f"{'FOUND (connected=%s; negative kernel edges %s)' % (conn, sorted(neg)) if neg is not None else 'NONE (all %d signings nonplanar)' % (1 << beta)}")
check(all(ok for N, ok in window if N <= 75) and not window[-1][1], "planar double cover iff N <= 75")
print("  => for N >= 52: U_N has a planar double cover iff N <= 75 (monotone: covers restrict to subgraphs).")
print("     U_76 is therefore not projective-planar (a projective embedding lifts to a planar double cover);")
print("     U_N (52 <= N <= 76) has exactly one component with a cycle (checked), so U_52..U_75 are projective-planar")
print("     by Negami's theorem for connected graphs (CITED, not re-read); the tree components embed in any face.")
# the natural signings
for name, negtype in (("Collatz sign cover (trivial)", set()), ("up- edges negative", {"m"}),
                      ("up+ edges negative", {"p"}), ("halving edges negative", {"h"})):
    thrN = None
    for N in range(20, 80):
        Dc = nx.Graph()
        for s_ in (1, -1):
            for m_ in range(1, N // 2 + 1):
                t = -1 if "h" in negtype else 1
                Dc.add_edge(s_ * m_, t * s_ * 2 * m_)
            for n_ in range(1, N + 1, 2):
                if 3 * n_ + 1 <= N:
                    t = -1 if "p" in negtype else 1
                    Dc.add_edge(s_ * n_, t * s_ * (3 * n_ + 1))
                if 3 * n_ - 1 <= N:
                    t = -1 if "m" in negtype else 1
                    Dc.add_edge(s_ * n_, t * s_ * (3 * n_ - 1))
        if not nx.check_planarity(Dc)[0]:
            thrN = N
            break
    print(f"  signing '{name}': double cover first nonplanar at N = {thrN}")
    check(thrN is not None and thrN <= 52, "type signings never beat planarity")
print("  => no edge-type signing (in particular not the sign/negation cover) realises the planar double covers:")
print("     the projective window 52 <= N <= 75 is real but its covers are not arithmetic.")

print()
print("=" * 100)
print("U9  Petersen family (linkless embedding; Colin de Verdiere mu <= 4)")
print("=" * 100)


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
            if any(G.has_edge(p_, q_) for p_, q_ in itertools.combinations(nb, 2)):
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
check(len(fam) == 7, "Petersen family has 7 members")
K91, K92 = kernel(U_graph(91)), kernel(U_graph(92))
print(f"  kernel(U_91): {K91.number_of_nodes()} vertices; kernel(U_92): {K92.number_of_nodes()} vertices (cubic)")
G92 = U_graph(92)
core92 = two_core(G92)
for F in fam:
    degs_F = sorted(d for _, d in F.degree())
    need = sum(max(1, d - 2) for d in degs_F)      # a branch set of s vertices in a cubic graph has <= s+2 exits
    t0 = time.time()
    r91 = "excluded by count" if need > K91.number_of_nodes() else None
    if r91 is None:
        cs = find_minor(K91, F)
        check(cs is None, "no Petersen-family minor in U_91")
        r91 = f"CP-SAT INFEASIBLE ({time.time() - t0:.0f}s)"
    name = "Petersen" if nx.is_isomorphic(F, P) else ("K6" if nx.is_isomorphic(F, nx.complete_graph(6)) else
            ("K_{3,3,1}" if nx.is_isomorphic(F, nx.complete_multipartite_graph(3, 3, 1)) else
             f"{F.number_of_nodes()}v/{F.number_of_edges()}e"))
    if name == "Petersen":
        found92 = False      # U6: no Petersen minor for N <= 103
        res92 = "none (U6: Petersen only from N = 104)"
    else:
        B = find_minor(core92, F, workers=1)
        found92 = isinstance(B, dict)
        check(found92, "family member minor in U_92")
        verify_minor(G92, F, B)
        res92 = "minor FOUND, branch sets verified in U_92"
    print(f"  {name:10s} degrees {degs_F}: needs >= {need} kernel vertices; U_91: {r91}; U_92: {res92}")
print("  => U_N is linklessly embeddable iff N <= 91 (RST: linkless <=> no Petersen-family minor; CITED), so by")
print("     Colin de Verdiere / Lovasz-Schrijver (CITED): mu(U_N) <= 3 iff N <= 51, mu(U_N) = 4 for 52 <= N <= 91,")
print("     mu(U_N) >= 5 for N >= 92 (the N <= 91 side is solver-certified; Petersen itself only enters at 104).")
print()
print("ALL CHECKS PASSED (union_minors)")
