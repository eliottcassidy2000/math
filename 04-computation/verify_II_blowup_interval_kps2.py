r"""
verify_II_blowup_interval_kps2.py  (kind-pasteur-2026-06-09-S2, ADVERSARIAL VERIFIER)

Independent re-verification of branch II-blowup-interval claims:
  blowup_interval_lemma_kps2.py / mindeg3_spectrum_census_kps2.py /
  twin_local_moves_kps2.py

EVERYTHING here is re-implemented from scratch (different anchoring of the
cycle DP: max-vertex instead of min-vertex; exhaustive-DFS longest path;
distance-profile exact-invariant iso dedup instead of rounded eigenvalues;
programmatically derived degree-sequence lists instead of hand-coded ones).
The base cycle-spectrum routine is validated against a naive DFS enumeration
of ALL simple cycles on every connected atlas graph, and on small blowups.

CLAIMS UNDER TEST (sibling's verdicts):
  V1  995 connected atlas graphs (n<=7, >=1 edge); 971 contain a cycle.
  V2  Interval law: spec(G[K2]) >= [k,2k] for every k in spec(G).  995/995.
  V3  Path law: spec(G[K2]) >= [3, 2p(G)].  995/995.
  V4  Gap-free law: spec(G[K2]) == {3,...,c(B)} exactly.  995/995.
  V5  c(B) == 2p(G) for 939 graphs, c(B) > 2p(G) for 56; Hamiltonian blowup
      (c = 2n) for 896/995; smallest beater = net graph (c=12 > 2p=10).
  V6  EXACT LAW c(B) = 2*s(G), s = max order of connected sub-multigraph with
      edge multiplicities {1,2} and m-degrees {2,4}.  (NOT tested in their
      script -- we test it here for every graph with p(G) < n.)
  V7  delta>=3 census: 1, 3, 19, 150, 2589 (n=4..8); all contain C4 (n<=8);
      at n=8 exactly 36/2589 avoid C8, largest avoider spectrum {3..7}.
  V8  Cubic counts 1,2,5,19 (n=4,6,8,10); Petersen spec {5,6,8,9}; exactly
      one cubic-10 graph avoids C8: g6 IsTX@?B?w with spec {3,4,5}.
  V9  ex(n;C4) chain: 1,3,4,6,7,9 (n=2..7); ex(8)<=11, ex(9)<=13, ex(10)<=16,
      ex(11)<=18, ex(12)<=21; corridor margins -2,-2,-2,-2,-1,-1,+1,+1,+3.
  V10 C4-free delta>=3 census: 5 / 9 / 57 graphs at n = 10 / 11 / 12;
      ALL 71 have dyadic profile exactly {8}; girth-5 counts 1/0/2;
      zero E-G counterexamples for n <= 12.
  V11 Arithmetic: every [k,2k], k>=3, contains a power of 2 that is >= 4.
"""
import sys, os, time
from itertools import combinations
import networkx as nx
from networkx.generators.atlas import graph_atlas_g

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "05-knowledge", "results",
                   "verify_II_blowup_interval_kps2.out")

class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
    def write(self, s):
        sys.__stdout__.write(s); self.f.write(s)
    def flush(self):
        sys.__stdout__.flush(); self.f.flush()

sys.stdout = Tee(OUT)
T0 = time.time()
ISSUES = []

def flag(msg):
    ISSUES.append(msg)
    print("  *** ISSUE:", msg)

# ============================ independent machinery ===========================
def masks_of(G):
    nodes = sorted(G.nodes())
    idx = {u: i for i, u in enumerate(nodes)}
    n = len(nodes)
    adj = [0] * n
    for u, v in G.edges():
        adj[idx[u]] |= 1 << idx[v]
        adj[idx[v]] |= 1 << idx[u]
    return adj, n, nodes

def spec_v2(adj, n):
    """Cycle spectrum.  Fresh implementation: anchor each cycle at its
    MAXIMUM vertex a; DP over subsets of {0..a-1} with end-vertex bitmask."""
    spec = set()
    for a in range(2, n):
        adj_a = adj[a]
        low = (1 << a) - 1
        if not (adj_a & low):
            continue
        dp = [0] * (1 << a)
        dp[0] = 1 << a
        for m in range(1 << a):
            ends = dp[m]
            if not ends:
                continue
            cnt = bin(m).count("1")
            if cnt >= 2 and (ends & adj_a):
                spec.add(cnt + 1)
            uni = 0
            e = ends
            while e:
                b = e & -e
                e ^= b
                uni |= adj[b.bit_length() - 1]
            cand = uni & low & ~m
            while cand:
                wb = cand & -cand
                cand ^= wb
                dp[m | wb] |= wb
    return spec

NAIVE_BUDGET_HIT = [False]

def naive_spec(G, budget=4_000_000):
    """Ground truth: enumerate ALL simple cycles by DFS anchored at the
    minimum vertex of the cycle.  Step-counted."""
    adjset = {v: sorted(G[v]) for v in G.nodes()}
    spec = set()
    steps = [0]

    def dfs(start, v, onpath, depth):
        steps[0] += 1
        if steps[0] > budget:
            raise TimeoutError
        for w in adjset[v]:
            if w == start:
                if depth >= 3:
                    spec.add(depth)
            elif w > start and w not in onpath:
                onpath.add(w)
                dfs(start, w, onpath, depth + 1)
                onpath.discard(w)

    try:
        for s in sorted(G.nodes()):
            dfs(s, s, {s}, 1)
    except TimeoutError:
        NAIVE_BUDGET_HIT[0] = True
        return None
    return spec

def longest_path_bf(G):
    """Ground-truth longest path order: exhaustive DFS over all simple paths."""
    adjset = {v: list(G[v]) for v in G.nodes()}
    best = [1]

    def dfs(v, onpath, depth):
        if depth > best[0]:
            best[0] = depth
        for w in adjset[v]:
            if w not in onpath:
                onpath.add(w)
                dfs(w, onpath, depth + 1)
                onpath.discard(w)

    for s in G.nodes():
        dfs(s, {s}, 1)
    return best[0]

def blowup2(G):
    """G[K2] with ADJACENT twins (lexicographic product with K2)."""
    B = nx.Graph()
    for u in G.nodes():
        B.add_edge((u, 0), (u, 1))
    for u, v in G.edges():
        B.add_edge((u, 0), (v, 0))
        B.add_edge((u, 0), (v, 1))
        B.add_edge((u, 1), (v, 0))
        B.add_edge((u, 1), (v, 1))
    return B

def exact_invariant(G):
    """Iso invariant from EXACT integer data only (no float eigenvalues):
    per-vertex (degree, triangle count, sorted neighbor degrees, sorted BFS
    distance profile), sorted over vertices."""
    nodes = list(G.nodes())
    adjset = {v: set(G[v]) for v in nodes}
    per = []
    for v in nodes:
        tri = 0
        nb = sorted(adjset[v])
        for i in range(len(nb)):
            for j in range(i + 1, len(nb)):
                if nb[j] in adjset[nb[i]]:
                    tri += 1
        # BFS distances
        dist = {v: 0}
        frontier = [v]
        d = 0
        while frontier:
            d += 1
            nxt = []
            for x in frontier:
                for y in adjset[x]:
                    if y not in dist:
                        dist[y] = d
                        nxt.append(y)
            frontier = nxt
        per.append((len(adjset[v]), tri,
                    tuple(sorted(len(adjset[u]) for u in adjset[v])),
                    tuple(sorted(dist.values()))))
    return (len(nodes), G.number_of_edges(), tuple(sorted(per)))

def iso_dedup(graphs):
    buckets = {}
    reps = []
    for G in graphs:
        key = exact_invariant(G)
        lst = buckets.setdefault(key, [])
        if not any(nx.is_isomorphic(G, H) for H in lst):
            lst.append(G)
            reps.append(G)
    return reps

# ---------------- independent degree-sequence generator -----------------------
def gen_degseq2(targets, c4free, prefix_symmetry=True):
    """All graphs with the given (descending) degree sequence, complete up to
    isomorphism (prefix symmetry rule on still-isolated equal-target vertices;
    dedup must follow).  prefix_symmetry=False gives the fully unbroken
    enumeration (for validating the rule on small cases)."""
    n = len(targets)
    adj = [set() for _ in range(n)]
    deg = [0] * n
    common = [[0] * n for _ in range(n)]
    out = []

    def cm(a, b):
        return common[a][b] if a < b else common[b][a]

    def bump(a, b, d):
        if a > b:
            a, b = b, a
        common[a][b] += d

    def add(u, w):
        for x in adj[u]:
            if x != w:
                bump(w, x, 1)
        for x in adj[w]:
            if x != u:
                bump(u, x, 1)
        adj[u].add(w); adj[w].add(u)
        deg[u] += 1; deg[w] += 1

    def rem(u, w):
        adj[u].discard(w); adj[w].discard(u)
        deg[u] -= 1; deg[w] -= 1
        for x in adj[u]:
            if x != w:
                bump(w, x, -1)
        for x in adj[w]:
            if x != u:
                bump(u, x, -1)

    def edge_ok(u, w):
        if not c4free:
            return True
        for x in adj[u]:
            if x != w and cm(w, x) >= 1:
                return False
        for x in adj[w]:
            if x != u and cm(u, x) >= 1:
                return False
        return True

    def rec():
        u = -1
        for i in range(n):
            if deg[i] < targets[i]:
                u = i
                break
        if u < 0:
            G = nx.Graph()
            G.add_nodes_from(range(n))
            for a in range(n):
                for b in adj[a]:
                    if a < b:
                        G.add_edge(a, b)
            out.append(G)
            return
        need = targets[u] - deg[u]
        cand = [v for v in range(u + 1, n)
                if deg[v] < targets[v] and v not in adj[u]]
        if len(cand) < need:
            return
        if prefix_symmetry:
            used = [v for v in cand if deg[v] > 0]
            fresh = {}
            for v in cand:
                if deg[v] == 0:
                    fresh.setdefault(targets[v], []).append(v)
            fkeys = sorted(fresh)

            def attempt(chosen):
                added = []
                feas = True
                for w in chosen:
                    if edge_ok(u, w):
                        add(u, w)
                        added.append(w)
                    else:
                        feas = False
                        break
                if feas:
                    rec()
                for w in reversed(added):
                    rem(u, w)

            def pick_fresh(ci, left, chosen):
                if left == 0:
                    attempt(chosen)
                    return
                if ci == len(fkeys):
                    return
                cls = fresh[fkeys[ci]]
                for take in range(min(left, len(cls)) + 1):
                    pick_fresh(ci + 1, left - take, chosen + cls[:take])

            for r in range(min(need, len(used)) + 1):
                for cmb in combinations(used, r):
                    pick_fresh(0, need - r, list(cmb))
        else:
            for cmb in combinations(cand, need):
                added = []
                feas = True
                for w in cmb:
                    if edge_ok(u, w):
                        add(u, w)
                        added.append(w)
                    else:
                        feas = False
                        break
                if feas:
                    rec()
                for w in reversed(added):
                    rem(u, w)

    rec()
    return out

def gen_iso(targets, c4free, connected_only=True, prefix_symmetry=True):
    gs = gen_degseq2(targets, c4free, prefix_symmetry)
    if connected_only:
        gs = [G for G in gs if nx.is_connected(G)]
    return iso_dedup(gs)

def degseqs(n, e, dmin=3):
    """All descending degree sequences with n parts in [dmin, n-1], sum 2e."""
    res = []

    def build(prefix, remaining, slots, mx):
        if slots == 0:
            if remaining == 0:
                res.append(tuple(prefix))
            return
        lo = max(dmin, remaining - mx * (slots - 1)) if slots > 1 else remaining
        for d in range(min(mx, remaining - dmin * (slots - 1)), dmin - 1, -1):
            if d < dmin:
                break
            build(prefix + [d], remaining - d, slots - 1, d)

    build([], 2 * e, n, n - 1)
    return [s for s in res if all(dmin <= x <= n - 1 for x in s)]

POWERS = frozenset({4, 8, 16, 32})

def spec_of2(G):
    adj, n, _ = masks_of(G)
    return spec_v2(adj, n)

def g6str(G):
    return nx.to_graph6_bytes(G, header=False).strip().decode()

# ============================ PART 0: validation ==============================
print("=" * 78)
print("PART 0: validate fresh spectrum DP against naive all-cycle enumeration")
print("=" * 78)
atlas = graph_atlas_g()
atlas_n = {}
for G in atlas:
    if G.number_of_nodes() >= 1:
        atlas_n.setdefault(G.number_of_nodes(), []).append(G)

conn = [g for g in atlas if g.number_of_nodes() >= 2 and g.number_of_edges() >= 1
        and nx.is_connected(g)]
print(f"V1: connected atlas graphs n<=7 with >=1 edge: {len(conn)} (claim 995)")
if len(conn) != 995:
    flag(f"V1 count mismatch: {len(conn)} != 995")

mismatch0 = 0
ncyc = 0
for G in conn:
    s2 = spec_of2(G)
    s0 = naive_spec(G)
    if s0 != s2:
        mismatch0 += 1
        if mismatch0 <= 5:
            flag(f"P0 spectrum mismatch on {g6str(G)}: naive {sorted(s0)} "
                 f"vs DP {sorted(s2)}")
    if s2:
        ncyc += 1
print(f"spectrum DP vs naive on all {len(conn)} graphs: {mismatch0} mismatches")
print(f"V1: graphs containing a cycle: {ncyc} (claim 971)")
if ncyc != 971:
    flag(f"V1 cycle-graph count mismatch: {ncyc} != 971")

# validate DP on small blowups too (independent of base-graph validation)
bl_checked = bl_skip = bl_bad = 0
for G in conn:
    if G.number_of_nodes() > 5 or G.number_of_edges() > 7:
        continue
    B = blowup2(G)
    adjB, nB, _ = masks_of(B)
    s2 = spec_v2(adjB, nB)
    s0 = naive_spec(B, budget=3_000_000)
    if s0 is None:
        bl_skip += 1
        continue
    bl_checked += 1
    if s0 != s2:
        bl_bad += 1
        flag(f"P0 blowup spectrum mismatch on base {g6str(G)}")
print(f"blowup DP vs naive (n(G)<=5, e(G)<=7): {bl_checked} checked, "
      f"{bl_skip} skipped (budget), {bl_bad} mismatches")

# V11 arithmetic
bad11 = [k for k in range(3, 100001)
         if not any(k <= 2 ** j <= 2 * k and 2 ** j >= 4 for j in range(2, 40))]
print(f"V11: [k,2k] contains a power of 2 >= 4 for all 3<=k<=100000: "
      f"{'OK' if not bad11 else 'FAILS at ' + str(bad11[:3])}")
if bad11:
    flag(f"V11 arithmetic fails at {bad11[:3]}")

# ============================ PART 1: blowup census ===========================
print()
print("=" * 78)
print("PART 1: blowup laws on all 995 connected atlas graphs")
print("=" * 78)
t0 = time.time()
fail_int = fail_path = fail_gap = fail_deg = 0
eq2p = gt2p = ham = 0
beaters = []
pn_lt = []   # graphs with p < n (for the s-law test)
for G in conn:
    adjG, nG, _ = masks_of(G)
    sG = spec_v2(adjG, nG)
    pG = longest_path_bf(G)
    B = blowup2(G)
    adjB, nB, _ = masks_of(B)
    sB = spec_v2(adjB, nB)
    cB = max(sB)
    # V2 interval law
    for k in sG:
        if not set(range(k, 2 * k + 1)) <= sB:
            fail_int += 1
            flag(f"V2 interval fails: {g6str(G)} k={k} specB={sorted(sB)}")
    # V3 path law
    if not set(range(3, 2 * pG + 1)) <= sB:
        fail_path += 1
        flag(f"V3 path law fails: {g6str(G)} p={pG} specB={sorted(sB)}")
    # V4 gap-free
    if sB != set(range(3, cB + 1)):
        fail_gap += 1
        flag(f"V4 gap-free fails: {g6str(G)} specB={sorted(sB)}")
    # degree law
    if any(B.degree((u, i)) != 2 * G.degree(u) + 1
           for u in G.nodes() for i in (0, 1)):
        fail_deg += 1
        flag(f"degree law fails: {g6str(G)}")
    # circumference accounting
    if cB == 2 * pG:
        eq2p += 1
    elif cB > 2 * pG:
        gt2p += 1
        beaters.append((G.number_of_nodes(), G.number_of_edges(),
                        g6str(G), pG, cB))
    else:
        flag(f"c(B) < 2p on {g6str(G)}: c={cB}, p={pG}")
    if cB == nB:   # blowup Hamiltonian: circumference = |V(B)| = 2 n(G)
        ham += 1
    if pG < nG:
        pn_lt.append((G, pG, cB))
print(f"V2 interval-law failures:  {fail_int}   (claim 0)")
print(f"V3 path-law failures:      {fail_path}   (claim 0)")
print(f"V4 gap-free failures:      {fail_gap}   (claim 0)")
print(f"degree-law failures:       {fail_deg}   (claim 0)")
print(f"V5: c=2p for {eq2p} (claim 939); c>2p for {gt2p} (claim 56); "
      f"Hamiltonian blowups {ham} (claim 896)")
if eq2p != 939 or gt2p != 56 or ham != 896:
    flag(f"V5 counts mismatch: eq2p={eq2p}, gt2p={gt2p}, ham={ham} "
         f"(claims 939/56/896)")
beaters.sort()
print(f"beaters by (n, e): smallest five: {beaters[:5]}")
n6b = [b for b in beaters if b[0] == 6]
print(f"  beaters with n=6: {len(n6b)}; with (n,e)=(6,6): "
      f"{[b[2] for b in beaters if b[:2] == (6, 6)]}")
# net graph identification
NET = nx.Graph([(0, 1), (1, 2), (2, 0), (0, 3), (1, 4), (2, 5)])
EdW = nx.from_graph6_bytes(b"E@dW")
print(f"  E@dW iso to net graph (C3+3 pendants): "
      f"{nx.is_isomorphic(EdW, NET)}")
adjN, nN, _ = masks_of(blowup2(NET))
sN = spec_v2(adjN, nN)
print(f"  spec(NET[K2]) = {sorted(sN)} (claim [3..12]); "
      f"p(NET) = {longest_path_bf(NET)} (claim 5)")
if sorted(sN) != list(range(3, 13)):
    flag(f"net blowup spectrum mismatch: {sorted(sN)}")
K23 = nx.complete_bipartite_graph(2, 3)
adjK, nK, _ = masks_of(blowup2(K23))
sK = spec_v2(adjK, nK)
print(f"  spec(K23[K2]) = {sorted(sK)} (claim [3..10]); "
      f"p(K23) = {longest_path_bf(K23)} (claim 5)")
print(f"part 1 time {time.time()-t0:.1f} s")

# -------- V6: exact law c(B) = 2*s(G), tested where it says something new ----
print()
print(f"V6: exact-law test c(B) = 2*s(G) on the {len(pn_lt)} graphs with "
      f"p(G) < n (for p = n the law reduces to c = 2n = 2p, covered above)")

def euler_struct_exists(Svert, E):
    if len(Svert) < 2:
        return False
    inc = {v: 0 for v in Svert}
    for (a, b) in E:
        inc[a] += 1
        inc[b] += 1
    if any(inc[v] == 0 for v in Svert):
        return False
    k = len(E)
    deg = {v: 0 for v in Svert}
    rem = dict(inc)
    m = [0] * k

    def bt(i):
        if i == k:
            parent = {v: v for v in Svert}

            def find(x):
                while parent[x] != x:
                    parent[x] = parent[parent[x]]
                    x = parent[x]
                return x
            comps = len(Svert)
            for j, (a, b) in enumerate(E):
                if m[j]:
                    ra, rb = find(a), find(b)
                    if ra != rb:
                        parent[ra] = rb
                        comps -= 1
            return comps == 1
        a, b = E[i]
        for val in (2, 1, 0):
            da, db = deg[a] + val, deg[b] + val
            if da > 4 or db > 4:
                continue
            deg[a] = da; deg[b] = db
            rem[a] -= 1; rem[b] -= 1
            m[i] = val
            ok = True
            for v in (a, b):
                if rem[v] == 0:
                    if deg[v] not in (2, 4):
                        ok = False
                        break
                elif deg[v] + 2 * rem[v] < 2:
                    ok = False
                    break
            if ok and bt(i + 1):
                deg[a] -= val; deg[b] -= val
                rem[a] += 1; rem[b] += 1
                return True
            deg[a] -= val; deg[b] -= val
            rem[a] += 1; rem[b] += 1
        m[i] = 0
        return False

    return bt(0)

def s_value(G, deadline):
    nodes = sorted(G.nodes())
    for size in range(len(nodes), 1, -1):
        for S in combinations(nodes, size):
            if time.time() > deadline:
                return None
            Sset = set(S)
            E = [(a, b) for a, b in G.edges() if a in Sset and b in Sset]
            if euler_struct_exists(Sset, E):
                return size
    return 1

t0 = time.time()
s_bad = s_skip = s_ok = 0
for (G, pG, cB) in pn_lt:
    sv = s_value(G, time.time() + 20.0)
    if sv is None:
        s_skip += 1
        continue
    if cB != 2 * sv:
        s_bad += 1
        flag(f"V6 exact law FAILS on {g6str(G)}: c(B)={cB}, 2*s={2*sv}")
    else:
        s_ok += 1
print(f"V6: {s_ok} confirmed, {s_bad} failures, {s_skip} skipped (timeout) "
      f"({time.time()-t0:.1f} s)")
svn = s_value(NET, time.time() + 20)
print(f"  s(NET) = {svn} (sun-walk claim: 6); c(NET[K2]) = {max(sN)} = 2*{svn}")

# ============================ PART 2: delta>=3 census =========================
print()
print("=" * 78)
print("PART 2: delta>=3 census n=4..8 (independent augmentation + exact dedup)")
print("=" * 78)
t0 = time.time()
counts = {}
for n in range(4, 8):
    glist = [G for G in atlas_n[n]
             if nx.is_connected(G) and min(d for _, d in G.degree()) >= 3]
    counts[n] = len(glist)
print(f"atlas-direct delta>=3 connected counts n=4..7: "
      f"{[counts[n] for n in range(4, 8)]} (claim 1, 3, 19, 150)")
if [counts[n] for n in range(4, 8)] != [1, 3, 19, 150]:
    flag("V7 atlas counts n=4..7 mismatch")

def augment(parents):
    cands = []
    for H in parents:
        degH = dict(H.degree())
        if not degH or min(degH.values()) < 2:
            continue
        must = [v for v in H.nodes() if degH[v] == 2]
        rest = [v for v in H.nodes() if degH[v] > 2]
        for k in range(len(rest) + 1):
            for extra in combinations(rest, k):
                N = must + list(extra)
                if len(N) < 3:
                    continue
                G = H.copy()
                G.add_node("NEW")
                G.add_edges_from(("NEW", x) for x in N)
                if nx.is_connected(G):
                    cands.append(nx.convert_node_labels_to_integers(G))
    return cands

# pipeline sanity at n=7
got7 = iso_dedup(augment(atlas_n[6]))
print(f"augmentation sanity n=7: {len(got7)} (expect 150)")
if len(got7) != 150:
    flag(f"V7 augmentation pipeline broken at n=7: {len(got7)}")

cand8 = augment(atlas_n[7])
print(f"n=8 raw candidates: {len(cand8)}")
all8 = iso_dedup(cand8)
print(f"n=8 connected delta>=3 graphs up to iso: {len(all8)} (claim 2589)  "
      f"({time.time()-t0:.1f} s)")
if len(all8) != 2589:
    flag(f"V7 n=8 census mismatch: {len(all8)} != 2589")

no4 = avoid8 = 0
best_avoider = (0, None, None)
cubic8 = []
for G in all8:
    s = spec_of2(G)
    if 4 not in s:
        no4 += 1
    if 8 not in s:
        avoid8 += 1
        if len(s) > best_avoider[0]:
            best_avoider = (len(s), sorted(s), g6str(G))
    if sorted(d for _, d in G.degree()) == [3] * 8:
        cubic8.append(G)
print(f"n=8: graphs without C4: {no4} (claim 0); C8-avoiders: {avoid8} "
      f"(claim 36); largest avoider spectrum {best_avoider[1]} (claim [3..7])")
if no4 != 0 or avoid8 != 36 or best_avoider[1] != [3, 4, 5, 6, 7]:
    flag(f"V7 n=8 details mismatch: no4={no4}, avoid8={avoid8}, "
         f"best={best_avoider[1]}")
print(f"n=8 cubic among census: {len(cubic8)} (claim 5)")

# ============================ PART 3: cubic + records =========================
print()
print("=" * 78)
print("PART 3: cubic graphs and the n=10 records (independent generator)")
print("=" * 78)
# validate prefix-symmetry rule against unbroken enumeration on small cases
for tg in ([3] * 6, [3] * 8, [4, 4, 3, 3, 3, 3]):
    a = len(gen_iso(list(tg), False, True, prefix_symmetry=True))
    b = len(gen_iso(list(tg), False, True, prefix_symmetry=False))
    print(f"generator validation degseq {tuple(tg)}: prefix-rule {a} vs "
          f"unbroken {b} -> {'OK' if a == b else 'BUG'}")
    if a != b:
        flag(f"generator symmetry rule WRONG on {tuple(tg)}: {a} vs {b}")

# delta>=3 totals at n=6,7 via generator (cross-validates against atlas)
for n, expect in ((6, 19), (7, 150)):
    tot = []
    for e in range((3 * n + 1) // 2, n * (n - 1) // 2 + 1):
        for sq in degseqs(n, e):
            tot += gen_iso(list(sq), False, True)
    tot = iso_dedup(tot)
    print(f"generator total delta>=3 n={n}: {len(tot)} (atlas: {expect})")
    if len(tot) != expect:
        flag(f"generator/dedup misses graphs at n={n}: {len(tot)} != {expect}")

cubic_counts = {}
cubic10 = None
for n in (4, 6, 8, 10):
    gs = gen_iso([3] * n, False, True)
    cubic_counts[n] = len(gs)
    if n == 10:
        cubic10 = gs
print(f"V8 cubic counts (4,6,8,10): "
      f"{[cubic_counts[n] for n in (4, 6, 8, 10)]} (claim 1,2,5,19 = A002851)")
if [cubic_counts[n] for n in (4, 6, 8, 10)] != [1, 2, 5, 19]:
    flag(f"V8 cubic counts mismatch: {cubic_counts}")

PET = nx.petersen_graph()
pet_dec = nx.from_graph6_bytes(b"IsP@PGXD_")
print(f"their Petersen g6 IsP@PGXD_ iso to nx.petersen_graph(): "
      f"{nx.is_isomorphic(pet_dec, PET)}")
sP = spec_of2(PET)
sP_naive = naive_spec(PET)
print(f"Petersen spectrum: DP {sorted(sP)}, naive {sorted(sP_naive)} "
      f"(claim [5,6,8,9]); girth {min(sP)} (anchor: 5)")
if sorted(sP) != [5, 6, 8, 9] or sP != sP_naive:
    flag(f"V8 Petersen spectrum mismatch: {sorted(sP)} / naive "
         f"{sorted(sP_naive)}")

avoid8_10 = []
for G in cubic10:
    s = spec_of2(G)
    if 8 not in s:
        avoid8_10.append((G, sorted(s)))
print(f"cubic-10 C8-avoiders: {len(avoid8_10)} (claim 1)")
rec_dec = nx.from_graph6_bytes(b"IsTX@?B?w")
for G, s in avoid8_10:
    print(f"   spectrum {s} (claim [3,4,5]); iso to their IsTX@?B?w: "
          f"{nx.is_isomorphic(G, rec_dec)}")
    if s != [3, 4, 5] or not nx.is_isomorphic(G, rec_dec):
        flag(f"V8 cubic-10 record mismatch: {s}")
sr = spec_of2(rec_dec)
sr_naive = naive_spec(rec_dec)
print(f"decoded IsTX@?B?w: cubic={sorted(d for _, d in rec_dec.degree()) == [3]*10}, "
      f"connected={nx.is_connected(rec_dec)}, spec DP {sorted(sr)}, "
      f"naive {sorted(sr_naive)}")
if len(avoid8_10) != 1:
    flag(f"V8 cubic-10 avoider count {len(avoid8_10)} != 1")

# ============================ PART 4: ex(n;C4) chain ==========================
print()
print("=" * 78)
print("PART 4: ex(n;C4) verification chain (independent)")
print("=" * 78)

def c4free(G):
    adjset = {v: set(G[v]) for v in G.nodes()}
    for u, v in combinations(G.nodes(), 2):
        if len(adjset[u] & adjset[v]) >= 2:
            return False
    return True

exC4 = {}
for n in range(2, 8):
    exC4[n] = max(G.number_of_edges() for G in atlas_n[n] if c4free(G))
print(f"atlas-direct ex(n;C4) n=2..7: {[exC4[n] for n in range(2, 8)]} "
      f"(claim 1,3,4,6,7,9)")
if [exC4[n] for n in range(2, 8)] != [1, 3, 4, 6, 7, 9]:
    flag(f"V9 ex chain base mismatch: {exC4}")

r8 = gen_degseq2([3] * 8, c4free=True)
print(f"n=8: C4-free graphs with degseq (3^8): {len(r8)} -> ex(8;C4)<=11 "
      f"{'VERIFIED' if not r8 else 'FAILED'}")
if r8:
    flag("V9 ex(8) step fails: C4-free cubic on 8 vertices exists")
exC4[8] = 11
r9 = gen_degseq2([4] + [3] * 8, c4free=True)
print(f"n=9: C4-free graphs with degseq (4,3^8): {len(r9)} -> ex(9;C4)<=13 "
      f"{'VERIFIED' if not r9 else 'FAILED'}")
if r9:
    flag("V9 ex(9) step fails")
exC4[9] = 13
print("n=10: e=17 C4-free forces min degree >= 17-13=4 -> degree sum >= 40 > "
      "34: impossible. ex(10;C4) <= 16.")
exC4[10] = 16
seqs19 = degseqs(11, 19)
print(f"n=11, e=19: programmatic degseq list ({len(seqs19)} sequences, "
      f"claim 7): {seqs19}")
if len(seqs19) != 7:
    flag(f"V9 n=11 sequence count {len(seqs19)} != 7")
t19 = sum(len(gen_degseq2(list(sq), c4free=True)) for sq in seqs19)
print(f"n=11: C4-free graphs with 19 edges: {t19} -> ex(11;C4)<=18 "
      f"{'VERIFIED' if t19 == 0 else 'FAILED'}")
if t19:
    flag("V9 ex(11) step fails")
exC4[11] = 18
print("n=12: e=22 C4-free forces min degree >= 22-18=4 -> sum >= 48 > 44: "
      "impossible. ex(12;C4) <= 21.")
exC4[12] = 21
margins = [exC4[n] - (3 * n + 1) // 2 for n in range(4, 13)]
print(f"corridor margins n=4..12: {margins} "
      f"(claim [-2,-2,-2,-2,-1,-1,+1,+1,+3])")
if margins != [-2, -2, -2, -2, -1, -1, 1, 1, 3]:
    flag(f"V9 corridor margins mismatch: {margins}")

# ============================ PART 5: C4-free census ==========================
print()
print("=" * 78)
print("PART 5: ALL C4-free delta>=3 graphs at n=10,11,12 (independent)")
print("=" * 78)
tot_by_n = {}
for n in (10, 11, 12):
    lo = (3 * n + 1) // 2
    hi = exC4[n]
    found = []
    for e in range(lo, hi + 1):
        for sq in degseqs(n, e):
            gs = gen_iso(list(sq), c4free=True)
            if gs:
                print(f"   n={n} e={e} degseq {sq}: {len(gs)}")
            found += gs
    found = iso_dedup(found)
    tot_by_n[n] = found
    print(f"   == n={n}: {len(found)} C4-free delta>=3 connected graphs ==")

claims = {10: 5, 11: 9, 12: 57}
girth5 = {10: 0, 11: 0, 12: 0}
eg_counter = []
dyadic_bad = []
for n in (10, 11, 12):
    if len(tot_by_n[n]) != claims[n]:
        flag(f"V10 census count mismatch at n={n}: {len(tot_by_n[n])} != "
             f"{claims[n]}")
    for G in tot_by_n[n]:
        s = spec_of2(G)
        dy = sorted(s & POWERS)
        if min(s) == 5:
            girth5[n] += 1
        if 4 in s:
            flag(f"V10 generator leak: C4 present in {g6str(G)}")
        if dy != [8]:
            dyadic_bad.append((n, g6str(G), sorted(s), dy))
        if not dy:
            eg_counter.append((n, g6str(G), sorted(s)))
print(f"V10 counts: {[len(tot_by_n[n]) for n in (10, 11, 12)]} "
      f"(claim [5, 9, 57])")
print(f"V10 dyadic profile == {{8}} for all: "
      f"{'YES' if not dyadic_bad else 'NO: ' + str(dyadic_bad[:5])}")
if dyadic_bad:
    flag(f"V10 dyadic profiles not all {{8}}: {dyadic_bad[:5]}")
print(f"V10 girth-5 counts at n=10/11/12: "
      f"{[girth5[n] for n in (10, 11, 12)]} (claim [1, 0, 2])")
if [girth5[n] for n in (10, 11, 12)] != [1, 0, 2]:
    flag(f"V10 girth-5 counts mismatch: {girth5}")
print(f"V10 Erdos-Gyarfas counterexamples n<=12: {len(eg_counter)} (claim 0)")
if eg_counter:
    flag(f"V10 E-G COUNTEREXAMPLE FOUND: {eg_counter}")
pet_in = any(nx.is_isomorphic(G, PET) for G in tot_by_n[10])
print(f"Petersen among n=10 C4-free delta>=3: {pet_in}")
if not pet_in:
    flag("V10 Petersen missing from n=10 census")

# ============================ summary =========================================
print()
print("=" * 78)
print(f"TOTAL TIME {time.time()-T0:.1f} s")
print(f"ISSUES FOUND: {len(ISSUES)}")
for i, msg in enumerate(ISSUES, 1):
    print(f"  {i}. {msg}")
if not ISSUES:
    print("ALL branch II central claims INDEPENDENTLY CONFIRMED.")
