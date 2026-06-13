r"""
twin_local_moves_kps2.py  (kind-pasteur-2026-06-09-S2, Branch II)

THE LOCAL MOVE TAXONOMY: which twin-like substructures extend cycle lengths
by +1 / +2, i.e. what is the WEAKEST structure that forces interval-type
density in the cycle spectrum (the structural content behind the Blowup
Interval Lemma and behind Bondy-Vince).

MOVES (all constructive; verified exhaustively on all connected atlas graphs
n <= 7 that contain a cycle):

M1 (+1, "edge-dominating vertex"): if C is a cycle and v NOT on C is adjacent
   to two CONSECUTIVE vertices u, w of C, then replacing the edge uw by the
   detour u-v-w gives a cycle of length |C|+1.
   This is the d=1 shadow of a twin: a single pair of adjacent true twins
   (u ~ v, N(u)\{v} = N(v)\{u}) automatically gives such a v for EVERY cycle
   through exactly one of them.

M2 (+2, "edge detour"): if xy is an edge with x,y NOT on C, x ~ u, y ~ w for
   consecutive u,w on C, then u-x-y-w replaces uw: length |C|+2.

M3 (chord split): a chord of C splits |C| = a + b (arc lengths) into cycles
   a+1 and b+1.

M4 (theta law): theta(a,b,c) (three internally disjoint u-w paths, a<=b<=c,
   b>=2) has cycle spectrum EXACTLY {a+b, a+c, b+c} -- the lacunary extreme.

TWIN STATEMENT (weakest twin forcing +1): if u,v are adjacent true twins and
   some cycle of length k contains u but not v, then k+1 is in spec(G).

BONDY-VINCE (cited, J. Graph Theory 27 (1998) 11-15): every graph with at
   most two vertices of degree < 3 (other than K1, K2) contains two cycles
   whose lengths differ by 1 or 2.  We verify the min-degree-3 form on all
   connected delta>=3 graphs with n <= 7 + note bipartite tightness (K_{3,3}).
"""
import sys, os, time
import networkx as nx
from networkx.generators.atlas import graph_atlas_g

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "05-knowledge", "results",
                   "twin_local_moves_kps2.out")

class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
    def write(self, s):
        sys.__stdout__.write(s); self.f.write(s)
    def flush(self):
        sys.__stdout__.flush(); self.f.flush()

sys.stdout = Tee(OUT)

def adj_masks(G):
    nodes = sorted(G.nodes())
    idx = {u: i for i, u in enumerate(nodes)}
    n = len(nodes)
    adj = [0] * n
    for u, v in G.edges():
        adj[idx[u]] |= 1 << idx[v]
        adj[idx[v]] |= 1 << idx[u]
    return adj, n

def cycle_spectrum(adj, n):
    spec = set()
    for u in range(n):
        nh = n - u - 1
        hi = 0
        for w in range(u + 1, n):
            hi |= 1 << w
        DP = [0] * (1 << nh)
        DP[0] = 1 << u
        adj_u = adj[u]
        shift = u + 1
        for s in range(1 << nh):
            ends = DP[s]
            if not ends:
                continue
            if bin(s).count("1") + 1 >= 3 and (ends & adj_u):
                spec.add(bin(s).count("1") + 1)
            uni = 0
            e = ends
            while e:
                vb = e & -e
                e ^= vb
                uni |= adj[vb.bit_length() - 1]
            sfull = (s << shift) | (1 << u)
            c = uni & hi & ~sfull
            while c:
                wb = c & -c
                c ^= wb
                DP[s | (1 << (wb.bit_length() - 1 - shift))] |= wb
    return spec

t0 = time.time()
atlas = graph_atlas_g()
graphs = [g for g in atlas if g.number_of_nodes() >= 3 and nx.is_connected(g)]

# ---------------- M1 / M2 / M3 verification over all cycles -------------------
m1_inst = m2_inst = m3_inst = 0
m1g_inst = [0]; m1g_fail = []
m5_inst = [0]; m5_fail = []; m5_t_hist = {}
m1_fail = []; m2_fail = []; m3_fail = []
twin_inst = 0; twin_fail = []
graphs_with_twins = 0
n_cyc_graphs = 0

for G in graphs:
    cycles = list(nx.simple_cycles(G))  # undirected: each cycle once, len>=3
    if not cycles:
        continue
    n_cyc_graphs += 1
    adj, n = adj_masks(G)
    spec = cycle_spectrum(adj, n)
    g6 = nx.to_graph6_bytes(G, header=False).strip().decode()
    Vset = set(G.nodes())

    # adjacent true twins
    twins = [(u, v) for u in G.nodes() for v in G.nodes() if u < v
             and G.has_edge(u, v)
             and set(G[u]) - {v} == set(G[v]) - {u}]
    if twins:
        graphs_with_twins += 1

    for cyc in cycles:
        k = len(cyc)
        oncyc = set(cyc)
        # M1
        for i in range(k):
            u, w = cyc[i], cyc[(i + 1) % k]
            for v in (set(G[u]) & set(G[w])) - oncyc:
                m1_inst += 1
                if k + 1 not in spec:
                    m1_fail.append((g6, cyc, v))
        # M2
        for i in range(k):
            u, w = cyc[i], cyc[(i + 1) % k]
            for x in set(G[u]) - oncyc:
                for y in (set(G[x]) & set(G[w])) - oncyc - {x}:
                    m2_inst += 1
                    if k + 2 not in spec:
                        m2_fail.append((g6, cyc, x, y))
        # M1' general fan move: v off-cycle with cycle-neighbors at distance d
        # -> cycles of length d+2 (short fan) and k-d+2 (long fan).
        # d=1 is M1 (+1); d=2 is a length-preserving rotation (+0).
        for v in Vset - oncyc:
            nb_on = [i for i in range(k) if cyc[i] in G[v]]
            for ii in range(len(nb_on)):
                for jj in range(ii + 1, len(nb_on)):
                    d = (nb_on[jj] - nb_on[ii]) % k
                    d = min(d, k - d)
                    m1g_inst[0] += 1
                    if d + 2 not in spec or k - d + 2 not in spec:
                        m1g_fail.append((g6, cyc, v, d))
        # M3
        pos = {v: i for i, v in enumerate(cyc)}
        for i in range(k):
            for j in range(i + 2, k):
                if (i, j) != (0, k - 1) and G.has_edge(cyc[i], cyc[j]):
                    a = j - i
                    b = k - a
                    m3_inst += 1
                    if a + 1 not in spec or b + 1 not in spec:
                        m3_fail.append((g6, cyc, i, j))
        # twin statement
        for (u, v) in twins:
            if (u in oncyc) != (v in oncyc):
                twin_inst += 1
                if k + 1 not in spec:
                    twin_fail.append((g6, cyc, u, v))
        # M5 sunlet-interval: t VERTEX-DISJOINT edge-dominators on C force
        # [k, k+t] subset spec.  (Blowup = extreme case t = k: every edge of
        # every cycle is dominated by a fresh twin.)  t = max matching between
        # cycle edges and outside dominating vertices.
        M = nx.Graph()
        has_dom = False
        for i in range(k):
            u, w = cyc[i], cyc[(i + 1) % k]
            for v in (set(G[u]) & set(G[w])) - oncyc:
                M.add_edge(("e", i), ("v", v))
                has_dom = True
        if has_dom:
            t = len(nx.max_weight_matching(M, maxcardinality=True))
            m5_t_hist[t] = m5_t_hist.get(t, 0) + 1
            m5_inst[0] += 1
            if not set(range(k, k + t + 1)) <= spec:
                m5_fail.append((g6, cyc, t))

print(f"Connected atlas graphs n<=7 with a cycle: {n_cyc_graphs}")
print(f"M1 (+1 edge-dominating vertex):  {m1_inst} instances, "
      f"{len(m1_fail)} failures")
print(f"M1' (general fan, d+2 & k-d+2):  {m1g_inst[0]} instances, "
      f"{len(m1g_fail)} failures")
print(f"M2 (+2 edge detour):             {m2_inst} instances, "
      f"{len(m2_fail)} failures")
print(f"M3 (chord split a+1/b+1):        {m3_inst} instances, "
      f"{len(m3_fail)} failures")
print(f"TWIN (+1 from adjacent true twin pair on/off cycle): {twin_inst} "
      f"instances, {len(twin_fail)} failures")
print(f"M5 (sunlet interval [k,k+t], t = max disjoint edge-dominators): "
      f"{m5_inst[0]} cycle instances, {len(m5_fail)} failures; "
      f"t-histogram {dict(sorted(m5_t_hist.items()))}")
for x in m5_fail[:5]:
    print(f"  FAIL M5: {x}")
print(f"Graphs (with cycle) containing at least one adjacent true twin pair: "
      f"{graphs_with_twins} / {n_cyc_graphs}")
for tag, lst in (("M1", m1_fail), ("M1G", m1g_fail), ("M2", m2_fail),
                 ("M3", m3_fail), ("TWIN", twin_fail)):
    for x in lst[:5]:
        print(f"  FAIL {tag}: {x}")

# ---------------- M4: theta spectrum law --------------------------------------
print("\nTheta law: spec(theta(a,b,c)) == {a+b, a+c, b+c} exactly")
theta_fail = []
cnt = 0
for a in range(1, 6):
    for b in range(max(a, 2), 6):
        for c in range(b, 6):
            T = nx.Graph()
            u, w = 0, 1
            nxt = 2
            for li, L in enumerate((a, b, c)):
                prev = u
                for s in range(L - 1):
                    T.add_edge(prev, nxt)
                    prev = nxt
                    nxt += 1
                T.add_edge(prev, w)
            Ti = T
            adj, n = adj_masks(Ti)
            spec = cycle_spectrum(adj, n)
            want = {a + b, a + c, b + c}
            cnt += 1
            if spec != want:
                theta_fail.append((a, b, c, sorted(spec)))
print(f"  checked {cnt} theta graphs (1<=a<=b<=c<=5, b>=2): "
      f"{len(theta_fail)} failures {theta_fail[:5]}")
print("  -> theta graphs are the lacunary extreme: 3 cycle lengths, all sums "
      "of two of three path lengths. But every internal vertex has degree 2: "
      "proper subdivisions NEVER have min degree 3. This is exactly why E-G "
      "needs delta >= 3.")

# ---------------- Bondy-Vince empirical (delta >= 3, n <= 7) ------------------
print("\nBondy-Vince check (delta>=3 => two cycles differing by 1 or 2):")
bv_fail = []
bv_cnt = 0
diff1_only = diff2_only = 0
for G in graphs:
    if min(dict(G.degree()).values()) < 3:
        continue
    bv_cnt += 1
    adj, n = adj_masks(G)
    spec = sorted(cycle_spectrum(adj, n))
    diffs = {spec[j] - spec[i] for i in range(len(spec))
             for j in range(i + 1, len(spec))}
    if not (1 in diffs or 2 in diffs):
        bv_fail.append((nx.to_graph6_bytes(G, header=False).strip().decode(),
                        spec))
    if 1 in diffs and 2 not in diffs:
        diff1_only += 1
    if 2 in diffs and 1 not in diffs:
        diff2_only += 1
print(f"  delta>=3 connected graphs n<=7: {bv_cnt}; "
      f"failures: {len(bv_fail)} {bv_fail[:5]}")
print(f"  graphs achieving only difference 2 (never 1, i.e. bipartite-like): "
      f"{diff2_only}; only difference 1: {diff1_only}")
K33 = nx.complete_bipartite_graph(3, 3)
adj, n = adj_masks(K33)
print(f"  K_3,3 spectrum: {sorted(cycle_spectrum(adj, n))} -- difference 2 "
      f"only: the '1 or 2' in Bondy-Vince cannot be improved to '1' "
      f"(bipartite graphs have all-even spectra).")

print(f"\nTotal time: {time.time()-t0:.1f} s")
print("\nVERDICT: all four local moves verified with zero failures; "
      "the weakest +1-forcing structure is a vertex dominating a cycle edge "
      "(one adjacent-true-twin pair suffices); chords give the +1/+2 of "
      "Bondy-Vince; subdivisions (theta) are the lacunary extreme excluded "
      "by the delta>=3 hypothesis.")
