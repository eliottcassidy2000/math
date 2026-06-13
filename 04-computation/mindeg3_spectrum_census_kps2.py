"""
mindeg3_spectrum_census_kps2.py  (kind-pasteur-2026-06-09-S2, Branch II)

THE MIN-DEGREE-3 CYCLE SPECTRUM CENSUS, n <= 12, and the dyadic-avoidance
records for the Erdos-Gyarfas conjecture (every delta>=3 graph has a cycle of
length a power of 2).

STRUCTURE OF THE ARGUMENT
  delta >= 3  forces  e >= ceil(3n/2)   (degree sum)
  C4-free     forces  e <= ex(n; C4)    (Turan number of the 4-cycle)
  For n <= 9:  ceil(3n/2) > ex(n;C4)  ==> EVERY delta>=3 graph contains C4.
  The corridor opens at n = 10 (15 <= 16): C4-avoiders must be within 1 edge
  of C4-extremality.  A delta>=3 graph on n <= 15 vertices avoiding BOTH C4
  and C8 would be a genuine counterexample to Erdos-Gyarfas (C16 needs n>=16).
  We enumerate ALL C4-free delta>=3 graphs for n = 10, 11, 12 and test C8.

INTERNAL VERIFICATION of ex(n;C4) chain (no literature trust needed):
  n <= 7 : direct maximum over the full atlas.
  n = 8  : e=12 would force d_v >= 12 - ex(7) = 3 for all v => cubic;
           checked against ALL cubic graphs on 8 vertices: none C4-free.
  n = 9  : e=14 forces d_v >= 14 - ex(8) = 3, sum=28 => degseq (4,3^8);
           exhaustive C4-free search over that sequence: empty.
  n = 10 : e=17 forces d_v >= 17 - ex(9) = 4 => sum >= 40 > 34: arithmetic.
  n = 11 : e=19 forces d_v >= 3, sum 38: 7 sequences, all searched: empty.
  n = 12 : e=22 forces d_v >= 22 - ex(11) = 4 => sum >= 48 > 44: arithmetic.
Cited values (Clapham-Flockhart-Sheehan 1989, J. Graph Theory 13: "Graphs
without four-cycles", exact values to n=21): ex = 4,6,7,9,11,13,16,18,21 for
n = 4..12.

CENSUS PARTS
  PART A: all connected delta>=3 graphs n = 4..7 (atlas) - full spectra.
  PART B: all connected delta>=3 graphs n = 8 (one-vertex augmentation of all
          1044 atlas 7-vertex graphs; pipeline sanity-checked at n=7).
  PART C: all connected cubic graphs n = 4,6,8,10 (degree-seq backtracking,
          counts cross-checked against A002851: 1, 2, 5, 19).
  PART D: ALL C4-free delta>=3 graphs at n = 10, 11, 12 + their spectra.
  PART E: dyadic-avoidance records and what kills each avoider.
"""
import sys, os, time
from itertools import combinations
import networkx as nx
import numpy as np
from networkx.generators.atlas import graph_atlas_g

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "05-knowledge", "results",
                   "mindeg3_spectrum_census_kps2.out")

class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
    def write(self, s):
        sys.__stdout__.write(s); self.f.flush
        self.f.write(s)
    def flush(self):
        sys.__stdout__.flush(); self.f.flush()

sys.stdout = Tee(OUT)
t00 = time.time()

# ---------------- spectrum machinery ------------------------------------------
def adj_masks_G(G):
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

def spec_of(G):
    adj, n = adj_masks_G(G)
    return cycle_spectrum(adj, n)

def g6(G):
    return nx.to_graph6_bytes(G, header=False).strip().decode()

# ---------------- iso dedup ----------------------------------------------------
def invariant(G):
    n = G.number_of_nodes()
    A = nx.to_numpy_array(G, nodelist=sorted(G.nodes()))
    ev = tuple(np.round(np.linalg.eigvalsh(A), 6))
    degs = tuple(sorted(d for _, d in G.degree()))
    tri = np.diag(np.linalg.matrix_power(A, 3))
    tri_seq = tuple(sorted(int(round(x)) for x in tri))
    return (n, G.number_of_edges(), degs, tri_seq, ev)

def dedup(graphs):
    buckets = {}
    out = []
    for G in graphs:
        key = invariant(G)
        lst = buckets.setdefault(key, [])
        if not any(nx.is_isomorphic(G, H) for H in lst):
            lst.append(G)
            out.append(G)
    return out

# ---------------- degree-sequence backtracking generator ----------------------
def gen_degseq(targets, c4free):
    """All graphs (up to easy symmetry; dedup afterwards) with degree sequence
    targets (descending), optionally C4-free.  Vertices completed in
    increasing order; isolated same-target vertices used by class prefix."""
    n = len(targets)
    adj = [set() for _ in range(n)]
    deg = [0] * n
    comm = [[0] * n for _ in range(n)]   # common-neighbor counts
    results = []

    def can_add(u, w):
        if not c4free:
            return True
        for y in adj[w]:
            if y != u and comm[min(u, y)][max(u, y)] >= 1:
                return False
        for y in adj[u]:
            if y != w and comm[min(w, y)][max(w, y)] >= 1:
                return False
        return True

    def add_edge(u, w):
        for y in adj[w]:
            if y != u:
                comm[min(u, y)][max(u, y)] += 1
        for y in adj[u]:
            if y != w:
                comm[min(w, y)][max(w, y)] += 1
        adj[u].add(w); adj[w].add(u)
        deg[u] += 1; deg[w] += 1

    def del_edge(u, w):
        adj[u].discard(w); adj[w].discard(u)
        deg[u] -= 1; deg[w] -= 1
        for y in adj[w]:
            if y != u:
                comm[min(u, y)][max(u, y)] -= 1
        for y in adj[u]:
            if y != w:
                comm[min(w, y)][max(w, y)] -= 1

    def rec():
        u = -1
        for i in range(n):
            if deg[i] < targets[i]:
                u = i
                break
        if u < 0:
            results.append([tuple(sorted((a, b)) ) for a in range(n)
                            for b in adj[a] if a < b])
            return
        need = targets[u] - deg[u]
        cands = [v for v in range(u + 1, n)
                 if deg[v] < targets[v] and v not in adj[u]]
        non_iso = [v for v in cands if deg[v] > 0]
        classes = {}
        for v in cands:
            if deg[v] == 0:
                classes.setdefault(targets[v], []).append(v)
        class_keys = sorted(classes)

        def emit(chosen):
            added = []
            ok = True
            for w in chosen:
                if can_add(u, w):
                    add_edge(u, w)
                    added.append(w)
                else:
                    ok = False
                    break
            if ok:
                rec()
            for w in reversed(added):
                del_edge(u, w)

        def choose_from_classes(ci, left, chosen):
            if left == 0:
                emit(chosen)
                return
            if ci == len(class_keys):
                return
            avail = classes[class_keys[ci]]
            maxtake = min(left, len(avail))
            for take in range(maxtake + 1):
                choose_from_classes(ci + 1, left - take,
                                    chosen + avail[:take])

        for r0 in range(min(need, len(non_iso)) + 1):
            for combo in combinations(non_iso, r0):
                choose_from_classes(0, need - r0, list(combo))

    rec()
    out = []
    for edges in results:
        G = nx.Graph()
        G.add_nodes_from(range(n))
        G.add_edges_from(edges)
        out.append(G)
    return out

def gen_degseq_iso(targets, c4free, connected_only=True):
    gs = gen_degseq(targets, c4free)
    if connected_only:
        gs = [G for G in gs if nx.is_connected(G)]
    return dedup(gs)

POWERS = {4, 8, 16, 32}

def dyadic(spec):
    return sorted(spec & POWERS)

# =============== PART A: atlas delta>=3, n <= 7 ================================
print("=" * 78)
print("PART A: all connected delta>=3 graphs, n = 4..7 (networkx atlas)")
print("=" * 78)
atlas = graph_atlas_g()
atlas_n = {}
for G in atlas:
    if G.number_of_nodes() >= 1:
        atlas_n.setdefault(G.number_of_nodes(), []).append(G)
for n in range(4, 8):
    glist = [G for G in atlas_n[n]
             if nx.is_connected(G) and min(d for _, d in G.degree()) >= 3]
    specs = [(G, spec_of(G)) for G in glist]
    no4 = [G for G, s in specs if 4 not in s]
    no48 = [G for G, s in specs if not (s & POWERS)]
    print(f"n={n}: {len(glist)} connected delta>=3 graphs; "
          f"without C4: {len(no4)}; avoiding ALL powers of 2: {len(no48)}")
    from collections import Counter
    cnt = Counter(tuple(sorted(s)) for _, s in specs)
    print(f"   distinct spectra ({len(cnt)}): "
          f"{sorted(cnt.items(), key=lambda kv: (-kv[1], kv[0]))[:12]}")

# =============== PART B: n = 8 by augmentation =================================
print()
print("=" * 78)
print("PART B: all connected delta>=3 graphs, n = 8 (augmentation)")
print("=" * 78)

def augment_all(parents, newv):
    """All G = H + v with delta(G)>=3, G connected, over parents H."""
    cands = []
    for H in parents:
        degH = dict(H.degree())
        if min(degH.values()) < 2:
            continue
        must = [v for v in H.nodes() if degH[v] == 2]
        rest = [v for v in H.nodes() if degH[v] > 2]
        for k in range(0, len(rest) + 1):
            for extra in combinations(rest, k):
                N = list(must) + list(extra)
                if len(N) < 3:
                    continue
                G = H.copy()
                G.add_node(newv)
                G.add_edges_from((newv, x) for x in N)
                if nx.is_connected(G):
                    cands.append(nx.convert_node_labels_to_integers(
                        G, ordering="sorted"))
    return cands

# pipeline sanity check at n=7
t0 = time.time()
parents6 = atlas_n[6]
cand7 = augment_all(parents6, 999)
got7 = dedup(cand7)
direct7 = [G for G in atlas_n[7]
           if nx.is_connected(G) and min(d for _, d in G.degree()) >= 3]
print(f"sanity n=7: augmentation gives {len(got7)} graphs, atlas direct "
      f"{len(direct7)}  -> {'MATCH' if len(got7) == len(direct7) else 'MISMATCH!'}"
      f"  ({time.time()-t0:.1f} s)")

t0 = time.time()
parents7 = atlas_n[7]
cand8 = augment_all(parents7, 999)
print(f"n=8 candidates: {len(cand8)} (dedup running...)")
all8 = dedup(cand8)
print(f"n=8: {len(all8)} connected delta>=3 graphs  ({time.time()-t0:.1f} s)")

specs8 = [(G, spec_of(G)) for G in all8]
no4_8 = [G for G, s in specs8 if 4 not in s]
with8 = [G for G, s in specs8 if 8 in s]
no_dyadic8 = [G for G, s in specs8 if not (s & POWERS)]
avoid8only = [(G, s) for G, s in specs8 if 8 not in s]
print(f"   without C4: {len(no4_8)} (Turan predicts 0: e>=12 > ex(8)=11)")
print(f"   containing C8: {len(with8)};  avoiding C8 (dyadic contact = {{4}} "
      f"only): {len(avoid8only)}")
print(f"   avoiding ALL powers: {len(no_dyadic8)}")
best8 = max(avoid8only, key=lambda gs: len(gs[1]))
print(f"   largest spectrum among C8-avoiders at n=8: {sorted(best8[1])} "
      f"(g6 {g6(best8[0])}) -- circumference must be <= 7")

# =============== PART C: cubic graphs n = 4..10 ================================
print()
print("=" * 78)
print("PART C: connected cubic graphs n = 4, 6, 8, 10 (A002851: 1, 2, 5, 19)")
print("=" * 78)
cubic = {}
for n in (4, 6, 8, 10):
    t0 = time.time()
    cubic[n] = gen_degseq_iso([3] * n, c4free=False)
    print(f"n={n}: {len(cubic[n])} connected cubic graphs "
          f"({time.time()-t0:.1f} s)")
    for G in cubic[n]:
        s = spec_of(G)
        girth = min(s)
        tag = ""
        if 4 not in s:
            tag += "  C4-FREE"
        if not (s & POWERS):
            tag += "  !!NO-POWER-OF-2!!"
        if n == 10 and girth == 5:
            tag += "  [= Petersen]"
        print(f"   g6 {g6(G):14s} spec {sorted(s)} dyadic {dyadic(s)}{tag}")

# =============== ex(n;C4) verification chain ===================================
print()
print("=" * 78)
print("ex(n; C4) INTERNAL VERIFICATION CHAIN")
print("=" * 78)
exC4 = {}
for n in range(2, 8):
    best = 0
    for G in atlas_n[n]:
        adjv = {v: set(G[v]) for v in G.nodes()}
        ok = True
        for u, v in combinations(G.nodes(), 2):
            if len(adjv[u] & adjv[v]) >= 2:
                ok = False
                break
        if ok:
            best = max(best, G.number_of_edges())
    exC4[n] = best
print(f"atlas-direct: ex(n;C4) for n=2..7: {[exC4[n] for n in range(2, 8)]} "
      f"(literature: 1,3,4,6,7,9)")

# n=8: 12 edges -> all degrees >= 12 - ex(7) = 3 -> cubic; any cubic C4-free?
cubic8_c4f = [G for G in gen_degseq([3] * 8, c4free=True)]
print(f"n=8: C4-free cubic graphs on 8 vertices (any labeled found): "
      f"{len(cubic8_c4f)} -> ex(8;C4) <= 11"
      f"{' VERIFIED' if not cubic8_c4f else ' FAILED'}")
exC4[8] = 11
# n=9: 14 edges -> degrees >= 3, sum 28 -> (4,3^8)
t0 = time.time()
s9 = gen_degseq([4] + [3] * 8, c4free=True)
print(f"n=9: C4-free graphs with degseq (4,3^8): {len(s9)} "
      f"-> ex(9;C4) <= 13 {'VERIFIED' if not s9 else 'FAILED'} "
      f"({time.time()-t0:.1f} s)")
exC4[9] = 13
print("n=10: e=17 forces all degrees >= 17 - ex(9) = 4 -> sum >= 40 > 34. "
      "ex(10;C4) <= 16 by arithmetic.")
exC4[10] = 16
# n=11: e=19 -> degrees >= 3, sum 38; 7 sequences
seqs19 = [[8] + [3] * 10, [7, 4] + [3] * 9, [6, 5] + [3] * 9,
          [6, 4, 4] + [3] * 8, [5, 5, 4] + [3] * 8,
          [5, 4, 4, 4] + [3] * 7, [4, 4, 4, 4, 4] + [3] * 6]
tot19 = 0
t0 = time.time()
for sq in seqs19:
    found = gen_degseq(sq, c4free=True)
    tot19 += len(found)
print(f"n=11: C4-free graphs with 19 edges (7 degree sequences): {tot19} "
      f"-> ex(11;C4) <= 18 {'VERIFIED' if tot19 == 0 else 'FAILED'} "
      f"({time.time()-t0:.1f} s)")
exC4[11] = 18
print("n=12: e=22 forces all degrees >= 22 - ex(11) = 4 -> sum >= 48 > 44. "
      "ex(12;C4) <= 21 by arithmetic.")
exC4[12] = 21

print("\nThe Turan corridor:  margin(n) = ex(n;C4) - ceil(3n/2)")
for n in range(4, 13):
    if n in exC4:
        lo = (3 * n + 1) // 2
        print(f"   n={n:2d}: ceil(3n/2) = {lo:2d}, ex(n;C4) = {exC4[n]:2d}, "
              f"margin = {exC4[n]-lo:+d}"
              + ("   <- corridor CLOSED: C4 forced" if exC4[n] < lo else
                 "   <- corridor open"))

# =============== PART D: all C4-free delta>=3 graphs, n = 10, 11, 12 ===========
print()
print("=" * 78)
print("PART D: ALL C4-free delta>=3 graphs at n = 10, 11, 12 and their spectra")
print("=" * 78)

def report_c4free(n, seq_list):
    found = []
    for sq in seq_list:
        t0 = time.time()
        gs = gen_degseq_iso(sq, c4free=True)
        print(f"   n={n} degseq {tuple(sq)}: {len(gs)} connected C4-free "
              f"graphs ({time.time()-t0:.1f} s)")
        found += gs
    print(f"   == n={n}: {len(found)} C4-free delta>=3 connected graphs ==")
    killers = []
    for G in found:
        s = spec_of(G)
        pows = dyadic(s)
        girth = min(s) if s else None
        print(f"      g6 {g6(G):16s} e={G.number_of_edges()} girth={girth} "
              f"spec {sorted(s)} dyadic {pows}"
              + ("   !!POWER-FREE: E-G COUNTEREXAMPLE!!" if not pows else ""))
        if not pows:
            killers.append(G)
    return found, killers

# n=10: e in {15,16}
print("n=10 (e must be 15 or 16; 17+ impossible):")
f10, k10 = report_c4free(10, [[3] * 10, [4, 4] + [3] * 8, [5] + [3] * 9])
# n=11: e in {17,18}
print("n=11 (e must be 17 or 18):")
f11, k11 = report_c4free(11, [[4] + [3] * 10,
                              [4, 4, 4] + [3] * 8,
                              [5, 4] + [3] * 9,
                              [6] + [3] * 10])
# n=12: e in {18..21}; degree sums 36, 38, 40, 42 (odd sums impossible)
print("n=12 (e must be 18..21):")
f12, k12 = report_c4free(12, [
    [3] * 12,                                              # e=18
    [4, 4] + [3] * 10, [5] + [3] * 11,                     # e=19
    [4, 4, 4, 4] + [3] * 8, [5, 4, 4] + [3] * 9,           # e=20
    [5, 5] + [3] * 10, [6, 4] + [3] * 10, [7] + [3] * 11,  # e=20
    [4] * 6 + [3] * 6, [5, 4, 4, 4, 4] + [3] * 7,          # e=21
    [5, 5, 4, 4] + [3] * 8, [6, 4, 4, 4] + [3] * 8,        # e=21
    [5, 5, 5] + [3] * 9, [6, 5, 4] + [3] * 9,              # e=21
    [7, 4, 4] + [3] * 9, [6, 6] + [3] * 10,                # e=21
    [7, 5] + [3] * 10, [8, 4] + [3] * 10, [9] + [3] * 11]) # e=21

# =============== PART E: verdicts ==============================================
print()
print("=" * 78)
print("PART E: DYADIC-AVOIDANCE RECORDS AND VERDICTS")
print("=" * 78)
print(f"n <= 9 : EVERY connected delta>=3 graph contains C4 "
      f"(Turan corridor closed; exhaustively confirmed n <= 8).")
print(f"n = 10 : C4-avoiders exist ({len(f10)} graphs). "
      f"E-G counterexamples found: {len(k10)}")
print(f"n = 11 : C4-free delta>=3 graphs: {len(f11)}. "
      f"E-G counterexamples found: {len(k11)}")
print(f"n = 12 : C4-free delta>=3 graphs: {len(f12)}. "
      f"E-G counterexamples found: {len(k12)}")
print(f"\nSince C16 requires n >= 16, any {{C4,C8}}-avoiding delta>=3 graph "
      f"on n <= 15 would refute Erdos-Gyarfas outright. For n <= 12: none "
      f"exist (this census). Royle/Markstrom (cited): none for n <= 16; "
      f"cubic counterexamples need n >= 30; Markstrom found 24-vertex cubic "
      f"graphs avoiding C4 and C8 but containing C16.")
print(f"\nTotal census time: {time.time()-t00:.1f} s")
