"""
verify2_II_blowup_adversarial_kps2.py  (kind-pasteur-2026-06-09-S2, SECOND adversarial
verifier of branch II-blowup-interval; independent of verify_II_blowup_interval_kps2.py)

Independently recomputes the load-bearing numbers of THM-456 with FRESH code:
  P0  fresh layered-frontier spectrum DP, validated against nx.simple_cycles
      (all connected n<=6 graphs + blowups of all connected n<=4 graphs) and
      closed-form anchors (C_k, K_n, K_{3,3}, Petersen {5,6,8,9}, theta(1,2,4)).
  P1  all 995 connected atlas graphs (2<=n<=7): LAW-1 interval, LAW-2 path floor,
      LAW-3 gap-free, degree law, c=2p counts (claim 939/56), Hamiltonian blowups
      (claim 896), NET graph beater (c=12 > 2p=10), and the exact-law upper object
      via an INDEPENDENT characterization: c(G[K2]) = 2*s'(G) where s'(G) = max
      number of distinct vertices over closed walks in G visiting each vertex at
      most twice (checked on ALL 995, not just the 144 with p<n).
  P2  n=8 delta>=3 census by my own augmentation + my own dedup (claim 2589),
      C4-free count (claim 0), C8-avoiders (claim 36), largest avoider spectrum
      (claim [3..7]), cubic among census (claim 5).
  P3  my own degree-sequence generator (fresh code), validated against atlas
      (n=6: 19, n=7: 150) and A002851 cubic counts (1,2,5,19); cubic-10 spectra:
      Petersen {5,6,8,9}, unique C8-avoider spec {3,4,5} iso to IsTX@?B?w.
  P4  ex(n;C4) chain n=2..12 re-derived (atlas direct n<=7; (3^8), (4,3^8),
      7 sequences at n=11 all empty; arithmetic steps at n=10, 12).
  P5  ALL C4-free delta>=3 graphs n=10,11,12 with my generator (claim 5/9/57),
      dyadic profile {8} for all 71, girth-5 counts [1,0,2], iso match with the
      five published n=10 g6 codes.
  P6  local-move instance recounts on all 971 cyclic graphs: M1 (claim 111358),
      M3 (claim 143523), TWIN (claim 23860), zero failures; Bondy-Vince on the
      173 delta>=3 graphs n<=7, diff-2-only graphs (claim 2, incl. K_{3,3}).
  P7  trivia: [k,2k] contains a power of 2 >= 4 for 3<=k<=100000.

Output saved to 05-knowledge/results/verify2_II_blowup_adversarial_kps2.out
"""
import sys, os, time
from itertools import combinations
import networkx as nx
import numpy as np
from networkx.generators.atlas import graph_atlas_g

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "05-knowledge", "results",
                   "verify2_II_blowup_adversarial_kps2.out")

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

def check(name, cond, detail=""):
    if cond:
        print(f"  [OK] {name}")
    else:
        ISSUES.append((name, detail))
        print(f"  [ISSUE] {name}  -- {detail}")

# ---------------- fresh machinery ----------------------------------------------
def masks(G):
    nodes = sorted(G.nodes())
    idx = {u: i for i, u in enumerate(nodes)}
    n = len(nodes)
    adj = [0] * n
    for u, v in G.edges():
        if u == v:
            continue
        adj[idx[u]] |= 1 << idx[v]
        adj[idx[v]] |= 1 << idx[u]
    return adj, n

def spectrum(adj, n):
    """Exact cycle spectrum. Layered frontier DP anchored at the cycle's
    minimum vertex; frontier = dict subset(full-space int) -> endpoint mask."""
    spec = set()
    full = (1 << n) - 1
    for a in range(n):
        if not (adj[a] & (full & ~((1 << (a + 1)) - 1))) and not (adj[a]):
            continue
        allowed = full & ~((1 << (a + 1)) - 1)      # only vertices > a
        adj_a = adj[a]
        cur = {1 << a: 1 << a}
        size = 1
        while cur:
            nxt = {}
            for S, ends in cur.items():
                if size >= 3 and (ends & adj_a):
                    spec.add(size)
                e = ends
                while e:
                    b = e & -e
                    e ^= b
                    u = b.bit_length() - 1
                    cand = adj[u] & allowed & ~S
                    while cand:
                        wb = cand & -cand
                        cand ^= wb
                        S2 = S | wb
                        nxt[S2] = nxt.get(S2, 0) | wb
            cur = nxt
            size += 1
        if spec == set(range(3, n + 1)):
            break                                   # spectrum already maximal
    return spec

def longest_path(adj, n):
    """Order of the longest path; unanchored layered frontier DP."""
    cur = {1 << v: 1 << v for v in range(n)}
    size = 1
    best = 1
    while cur:
        best = size
        nxt = {}
        for S, ends in cur.items():
            e = ends
            while e:
                b = e & -e
                e ^= b
                u = b.bit_length() - 1
                cand = adj[u] & ~S
                while cand:
                    wb = cand & -cand
                    cand ^= wb
                    S2 = S | wb
                    nxt[S2] = nxt.get(S2, 0) | wb
        cur = nxt
        size += 1
    return best

def s_walk(adj, n):
    """s'(G) = max #distinct vertices over closed walks in G that visit each
    vertex at most twice (>=2 distinct vertices required).  Independent
    characterization of the blowup circumference: c(G[K2]) = 2*s'(G).
    [cycle in G[K2] -> project (twin edges = stutters) -> closed <=2-visit walk;
    closed <=2-visit walk on q vertices -> lift, stuttering single visits -> 2q.]
    Walk anchored at its minimum vertex; counts kept 2 bits per vertex."""
    best = 0
    for a in range(n):
        if best >= n - a:           # cannot improve using only vertices >= a
            break
        if not adj[a]:
            continue
        start = (a, 1)              # count vector over positions a..n-1; a has 1
        seen = {start}
        stack = [start]
        adj_shift = [adj[v] >> a for v in range(n)]
        while stack:
            v, cnt = stack.pop()
            if (adj[v] >> a) & 1:   # can close back to anchor a
                q = 0
                c = cnt
                while c:
                    if c & 3:
                        q += 1
                    c >>= 2
                if q >= 2 and q > best:
                    best = q
            nb = adj_shift[v]
            pos = 0
            while nb:
                if nb & 1:
                    f = (cnt >> (2 * pos)) & 3
                    if f < 2:
                        st = (a + pos, cnt + (1 << (2 * pos)))
                        if st not in seen:
                            seen.add(st)
                            stack.append(st)
                nb >>= 1
                pos += 1
    return best

def blowup(G):
    B = nx.Graph()
    for u in G.nodes():
        B.add_edge((u, 0), (u, 1))
    for u, v in G.edges():
        for i in (0, 1):
            for j in (0, 1):
                B.add_edge((u, i), (v, j))
    return B

def naive_spec(G):
    """Ground truth via networkx's own cycle enumeration (independent algo)."""
    return {len(c) for c in nx.simple_cycles(G)}

def spec_of(G):
    adj, n = masks(G)
    return spectrum(adj, n)

def g6s(G):
    return nx.to_graph6_bytes(G, header=False).strip().decode()

# ---------------- iso dedup (own invariant + VF2) ------------------------------
def invariant(G):
    nodes = sorted(G.nodes())
    A = nx.to_numpy_array(G, nodelist=nodes)
    ev = tuple(np.round(np.sort(np.linalg.eigvalsh(A)), 5))
    degs = tuple(sorted(d for _, d in G.degree()))
    A3 = np.linalg.matrix_power(A, 3)
    tri = tuple(sorted(int(round(x)) for x in np.diag(A3)))
    return (G.number_of_nodes(), G.number_of_edges(), degs, tri, ev)

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

# =====================  P0: validate fresh machinery  ==========================
print("=" * 78)
print("P0: validate fresh spectrum/path/walk machinery")
print("=" * 78)
atlas = graph_atlas_g()
atlas_by_n = {}
for G in atlas:
    if G.number_of_nodes() >= 1:
        atlas_by_n.setdefault(G.number_of_nodes(), []).append(G)

conn = [G for n in range(2, 8) for G in atlas_by_n[n]
        if nx.is_connected(G) and G.number_of_edges() >= 1]
check("995 connected atlas graphs with an edge, 2<=n<=7",
      len(conn) == 995, f"got {len(conn)}")

mism = 0
for G in conn:
    if G.number_of_nodes() <= 6:
        if spec_of(G) != naive_spec(G):
            mism += 1
check("fresh DP == nx.simple_cycles on all connected graphs n<=6", mism == 0,
      f"{mism} mismatches")

mismB = 0
nb_checked = 0
for G in conn:
    if G.number_of_nodes() <= 4:
        B = blowup(G)
        Bi = nx.convert_node_labels_to_integers(B, ordering="sorted")
        if spec_of(Bi) != naive_spec(Bi):
            mismB += 1
        nb_checked += 1
check(f"fresh DP == nx.simple_cycles on blowups of all connected n<=4 "
      f"({nb_checked} blowups, up to 8 vertices)", mismB == 0,
      f"{mismB} mismatches")

# closed-form anchors
anch_ok = True
for k in range(3, 10):
    if spec_of(nx.cycle_graph(k)) != {k}:
        anch_ok = False
for m in range(4, 9):
    if spec_of(nx.complete_graph(m)) != set(range(3, m + 1)):
        anch_ok = False
if spec_of(nx.complete_bipartite_graph(3, 3)) != {4, 6}:
    anch_ok = False
pet = nx.petersen_graph()
pet_spec = spec_of(pet)
if pet_spec != {5, 6, 8, 9}:
    anch_ok = False
theta124 = nx.Graph([(0, 1), (0, 2), (2, 1), (0, 3), (3, 4), (4, 5), (5, 1)])
if spec_of(theta124) != {3, 5, 6}:
    anch_ok = False
check("closed-form anchors: C_k={k}, K_m=[3..m], K33={4,6}, "
      "Petersen={5,6,8,9}, theta(1,2,4)={3,5,6}", anch_ok,
      f"Petersen got {sorted(pet_spec)}")
check("Petersen girth 5 (anchor)", min(pet_spec) == 5, f"got {min(pet_spec)}")

# longest-path validation against brute force on small graphs
def brute_lp(G):
    best = 1
    nodes = list(G.nodes())
    def dfs(v, vis):
        nonlocal best
        best = max(best, len(vis))
        for w in G[v]:
            if w not in vis:
                dfs(w, vis | {w})
    for v in nodes:
        dfs(v, {v})
    return best
lp_mism = sum(1 for G in conn if G.number_of_nodes() <= 5
              and longest_path(*masks(G)) != brute_lp(G))
check("fresh longest-path DP == brute DFS on all connected n<=5", lp_mism == 0,
      f"{lp_mism} mismatches")
print(f"P0 time {time.time()-T0:.1f} s")

# =====================  P1: 995-graph blowup laws  =============================
print()
print("=" * 78)
print("P1: blowup laws on all 995 connected atlas graphs (fresh recompute)")
print("=" * 78)
t1 = time.time()
f_law1 = f_law2 = f_law3 = f_deg = f_2s = 0
n_eq2p = n_gt2p = n_lt2p = n_ham = 0
beaters = []
delta3_specs = []          # for P6 Bondy-Vince
for G in conn:
    n = G.number_of_nodes()
    adjG, nG = masks(G)
    sG = spectrum(adjG, nG)
    pG = longest_path(adjG, nG)
    sprime = s_walk(adjG, nG)
    B = nx.convert_node_labels_to_integers(blowup(G), ordering="sorted")
    adjB, nB = masks(B)
    sB = spectrum(adjB, nB)
    cB = max(sB)
    # LAW-1
    for k in sG:
        if not set(range(k, 2 * k + 1)) <= sB:
            f_law1 += 1
    # LAW-2
    if not set(range(3, 2 * pG + 1)) <= sB:
        f_law2 += 1
    # LAW-3 gap-free
    if sB != set(range(3, cB + 1)):
        f_law3 += 1
    # degree law
    degsG = dict(G.degree())
    for (u, i) in B.nodes() if False else []:
        pass
    Braw = blowup(G)
    if any(Braw.degree((u, i)) != 2 * degsG[u] + 1
           for u in G.nodes() for i in (0, 1)):
        f_deg += 1
    # exact upper object: c = 2*s'
    if cB != 2 * sprime:
        f_2s += 1
        print(f"    c != 2s' on {g6s(G)}: c={cB}, s'={sprime}")
    # circumference comparison
    if cB == 2 * pG:
        n_eq2p += 1
    elif cB > 2 * pG:
        n_gt2p += 1
        beaters.append((n, G.number_of_edges(), g6s(G), pG, cB))
    else:
        n_lt2p += 1
    if cB == 2 * n:
        n_ham += 1
    if n >= 4 and min(d for _, d in G.degree()) >= 3:
        delta3_specs.append((G, sorted(sG)))

check("LAW-1 interval [k,2k] subset spec(B): 0 failures", f_law1 == 0,
      f"{f_law1}")
check("LAW-2 path floor [3,2p] subset spec(B): 0 failures", f_law2 == 0,
      f"{f_law2}")
check("LAW-3 spec(B) gap-free = [3,c]: 0 failures", f_law3 == 0, f"{f_law3}")
check("degree law deg(u^i)=2deg(u)+1: 0 failures", f_deg == 0, f"{f_deg}")
check("EXACT LAW c(G[K2]) = 2*s'(G) on ALL 995 (independent walk-based s')",
      f_2s == 0, f"{f_2s} failures")
check("c = 2p for 939 graphs", n_eq2p == 939, f"got {n_eq2p}")
check("c > 2p for 56 graphs", n_gt2p == 56, f"got {n_gt2p}")
check("c < 2p never", n_lt2p == 0, f"got {n_lt2p}")
check("Hamiltonian blowups: 896", n_ham == 896, f"got {n_ham}")

NET = nx.Graph([(0, 1), (1, 2), (2, 0), (0, 3), (1, 4), (2, 5)])
adjN, nN = masks(NET)
sN = spectrum(adjN, nN)
pN = longest_path(adjN, nN)
spN = s_walk(adjN, nN)
BN = nx.convert_node_labels_to_integers(blowup(NET), ordering="sorted")
sBN = spectrum(*masks(BN))
check("NET: spec={3}, p=5, s'=6, spec(NET[K2])=[3..12]",
      sN == {3} and pN == 5 and spN == 6 and sBN == set(range(3, 13)),
      f"spec {sorted(sN)}, p {pN}, s' {spN}, blowup {sorted(sBN)}")
b6 = sorted(b for b in beaters if b[0] == 6)
check("smallest beaters have n=6 (none at n<=5)",
      min(b[0] for b in beaters) == 6, f"min n = {min(b[0] for b in beaters)}")
print(f"  beaters at n=6: {b6}")
net_iso = any(nx.is_isomorphic(nx.from_graph6_bytes(b[2].encode()), NET)
              for b in b6 if b[1] == 6)
check("NET graph is among the (n=6,e=6) beaters", net_iso, "not found")
print(f"P1 time {time.time()-t1:.1f} s")

# =====================  P2: n=8 census via own augmentation  ===================
print()
print("=" * 78)
print("P2: n=8 delta>=3 census (own augmentation + own dedup)")
print("=" * 78)
t2 = time.time()
direct = {}
for n in range(4, 8):
    direct[n] = [G for G in atlas_by_n[n]
                 if nx.is_connected(G) and min(d for _, d in G.degree()) >= 3]
check("atlas-direct delta>=3 counts n=4..7 = 1,3,19,150",
      [len(direct[n]) for n in range(4, 8)] == [1, 3, 19, 150],
      f"got {[len(direct[n]) for n in range(4, 8)]}")

def augment(parents, tag):
    cands = []
    for H in parents:
        deg = dict(H.degree())
        if not deg or min(deg.values()) < 2:
            continue
        must = [v for v in H.nodes() if deg[v] == 2]
        rest = [v for v in H.nodes() if deg[v] >= 3]
        for k in range(len(rest) + 1):
            for extra in combinations(rest, k):
                N = must + list(extra)
                if len(N) < 3:
                    continue
                G = H.copy()
                G.add_edges_from((tag, x) for x in N)
                if nx.is_connected(G):
                    cands.append(nx.convert_node_labels_to_integers(
                        G, ordering="sorted"))
    return cands

# sanity at n=7 first
got7 = dedup(augment(atlas_by_n[6], 999))
check("augmentation sanity at n=7: 150", len(got7) == 150, f"got {len(got7)}")

cand8 = augment(atlas_by_n[7], 999)
print(f"  raw n=8 candidates: {len(cand8)}")
all8 = dedup(cand8)
check("n=8 connected delta>=3 census = 2589", len(all8) == 2589,
      f"got {len(all8)}")
spec8 = [(G, spec_of(G)) for G in all8]
no4 = [G for G, s in spec8 if 4 not in s]
avoid8 = [(G, s) for G, s in spec8 if 8 not in s]
nopow = [G for G, s in spec8 if not (s & {4, 8})]
cubic8 = [G for G in all8 if all(d == 3 for _, d in G.degree())]
check("n=8: C4-free count = 0", len(no4) == 0, f"got {len(no4)}")
check("n=8: C8-avoiders = 36", len(avoid8) == 36, f"got {len(avoid8)}")
check("n=8: power-of-2 avoiders = 0", len(nopow) == 0, f"got {len(nopow)}")
best8 = max((sorted(s) for _, s in avoid8), key=len) if avoid8 else []
check("n=8: largest C8-avoider spectrum = [3,4,5,6,7]",
      best8 == [3, 4, 5, 6, 7], f"got {best8}")
check("n=8: cubic graphs in census = 5", len(cubic8) == 5,
      f"got {len(cubic8)}")
print(f"P2 time {time.time()-t2:.1f} s")

# =====================  P3: own degree-sequence generator  =====================
print()
print("=" * 78)
print("P3: own degree-sequence generator + cubic spectra")
print("=" * 78)
t3 = time.time()

def gen_seq(targets, c4free):
    """All labeled graphs with degree sequence `targets` (non-increasing),
    optionally C4-free; symmetry cut: untouched (degree-0) vertices of equal
    target are interchangeable -> only lowest-index prefixes chosen."""
    n = len(targets)
    adj = [set() for _ in range(n)]
    deg = [0] * n
    res = []

    def creates_c4(u, w):
        for y in adj[w]:
            if y == u:
                continue
            for z in adj[y]:
                if z != u and z != w and z in adj[u]:
                    return True
        return False

    def rec():
        u = -1
        for i in range(n):
            if deg[i] < targets[i]:
                u = i
                break
        if u < 0:
            res.append([(a, b) for a in range(n) for b in adj[a] if a < b])
            return
        need = targets[u] - deg[u]
        cands = [v for v in range(u + 1, n)
                 if deg[v] < targets[v] and v not in adj[u]]
        touched = [v for v in cands if deg[v] > 0]
        fresh = {}
        for v in cands:
            if deg[v] == 0:
                fresh.setdefault(targets[v], []).append(v)
        fkeys = sorted(fresh)

        def place(chosen):
            added = []
            good = True
            for w in chosen:
                if c4free and creates_c4(u, w):
                    good = False
                    break
                for y in adj[w]:
                    pass
                adj[u].add(w); adj[w].add(u)
                deg[u] += 1; deg[w] += 1
                added.append(w)
            if good:
                rec()
            for w in reversed(added):
                adj[u].discard(w); adj[w].discard(u)
                deg[u] -= 1; deg[w] -= 1

        def fill_fresh(ki, left, chosen):
            if left == 0:
                place(chosen)
                return
            if ki == len(fkeys):
                return
            avail = fresh[fkeys[ki]]
            for take in range(min(left, len(avail)) + 1):
                fill_fresh(ki + 1, left - take, chosen + avail[:take])

        for r in range(min(need, len(touched)) + 1):
            for combo in combinations(touched, r):
                fill_fresh(0, need - r, list(combo))

    rec()
    out = []
    for edges in res:
        G = nx.Graph()
        G.add_nodes_from(range(n))
        G.add_edges_from(edges)
        out.append(G)
    return out

def gen_iso(targets, c4free=False):
    gs = [G for G in gen_seq(list(targets), c4free) if nx.is_connected(G)]
    return dedup(gs)

def all_degseqs(n, parts_min, parts_max, esum):
    """non-increasing sequences of length n, values in [parts_min, parts_max],
    summing to esum"""
    out = []
    def rec(pos, prev, rem, seq):
        if pos == n:
            if rem == 0:
                out.append(tuple(seq))
            return
        lo = parts_min
        hi = min(prev, parts_max, rem - parts_min * (n - pos - 1))
        for v in range(hi, lo - 1, -1):
            rec(pos + 1, v, rem - v, seq + [v])
    rec(0, parts_max, esum, [])
    return out

# validation: full delta>=3 census at n=6 and n=7 via degree sequences
for nval, expect in ((6, 19), (7, 150)):
    found = []
    for ee in range(3 * nval, nval * (nval - 1) + 1):
        if ee % 2:
            continue
        for sq in all_degseqs(nval, 3, nval - 1, ee):
            found += gen_iso(sq)
    found = dedup(found)
    check(f"generator full delta>=3 census n={nval} = {expect}",
          len(found) == expect, f"got {len(found)}")

cubic_counts = []
cubic10 = None
for nval in (4, 6, 8, 10):
    cg = gen_iso([3] * nval)
    cubic_counts.append(len(cg))
    if nval == 10:
        cubic10 = cg
check("cubic counts n=4,6,8,10 = 1,2,5,19 (A002851)",
      cubic_counts == [1, 2, 5, 19], f"got {cubic_counts}")

pet_found = [G for G in cubic10 if nx.is_isomorphic(G, pet)]
check("Petersen among the 19 cubic-10 graphs", len(pet_found) == 1,
      f"got {len(pet_found)}")
c8avoid10 = [(G, spec_of(G)) for G in cubic10 if 8 not in spec_of(G)]
check("exactly 1 cubic-10 C8-avoider", len(c8avoid10) == 1,
      f"got {len(c8avoid10)}")
if len(c8avoid10) == 1:
    Gx, sx = c8avoid10[0]
    pub = nx.from_graph6_bytes(b"IsTX@?B?w")
    check("cubic-10 C8-avoider spectrum {3,4,5} and iso to IsTX@?B?w",
          sorted(sx) == [3, 4, 5] and nx.is_isomorphic(Gx, pub),
          f"spec {sorted(sx)}")
print(f"P3 time {time.time()-t3:.1f} s")

# =====================  P4: ex(n;C4) chain  ====================================
print()
print("=" * 78)
print("P4: ex(n;C4) chain re-derivation")
print("=" * 78)
t4 = time.time()
def is_c4_free(G):
    nodes = list(G.nodes())
    for u, v in combinations(nodes, 2):
        if len(set(G[u]) & set(G[v])) >= 2:
            return False
    return True

exC4 = {}
for n in range(2, 8):
    exC4[n] = max(G.number_of_edges() for G in atlas_by_n[n] if is_c4_free(G))
check("ex(n;C4) n=2..7 = 1,3,4,6,7,9 (atlas direct, own C4 test)",
      [exC4[n] for n in range(2, 8)] == [1, 3, 4, 6, 7, 9],
      f"got {[exC4[n] for n in range(2, 8)]}")
g8 = gen_seq([3] * 8, c4free=True)
check("no C4-free cubic graph on 8 vertices => ex(8;C4) <= 11", not g8,
      f"got {len(g8)} labeled graphs")
exC4[8] = 11
g9 = gen_seq([4] + [3] * 8, c4free=True)
check("no C4-free (4,3^8) graph on 9 vertices => ex(9;C4) <= 13", not g9,
      f"got {len(g9)}")
exC4[9] = 13
exC4[10] = 16   # e=17 forces min degree >= 17-13=4 -> sum >= 40 > 34
seqs11 = all_degseqs(11, 3, 8, 38)
tot11 = sum(len(gen_seq(list(sq), c4free=True)) for sq in seqs11)
check(f"n=11 e=19: {len(seqs11)} degree sequences (claim 7), all empty "
      f"=> ex(11;C4) <= 18", len(seqs11) == 7 and tot11 == 0,
      f"seqs {len(seqs11)}, found {tot11}")
exC4[11] = 18
exC4[12] = 21   # e=22 forces min degree >= 22-18=4 -> sum >= 48 > 44
margins = [exC4[n] - ((3 * n + 1) // 2) for n in range(4, 13)]
check("Turan corridor margins n=4..12 = [-2,-2,-2,-2,-1,-1,+1,+1,+3]",
      margins == [-2, -2, -2, -2, -1, -1, 1, 1, 3], f"got {margins}")
print(f"P4 time {time.time()-t4:.1f} s")

# =====================  P5: C4-free delta>=3 census n=10,11,12  ================
print()
print("=" * 78)
print("P5: ALL C4-free delta>=3 graphs n=10,11,12 (own generator)")
print("=" * 78)
t5 = time.time()
census = {}
for nval in (10, 11, 12):
    emin = (3 * nval + 1) // 2
    emax = exC4[nval]
    found = []
    for ee in range(emin, emax + 1):
        if (2 * ee) % 2:
            continue
        for sq in all_degseqs(nval, 3, nval - 1, 2 * ee):
            found += gen_iso(sq, c4free=True)
    census[nval] = dedup(found)
    print(f"  n={nval}: {len(census[nval])} C4-free delta>=3 connected graphs "
          f"(e in [{emin},{emax}])")
check("C4-free delta>=3 counts n=10,11,12 = 5, 9, 57",
      [len(census[n]) for n in (10, 11, 12)] == [5, 9, 57],
      f"got {[len(census[n]) for n in (10, 11, 12)]}")

all_dyadic_8 = True
girth5 = []
eg_cex = 0
for nval in (10, 11, 12):
    g5 = 0
    for G in census[nval]:
        s = spec_of(G)
        dy = sorted(s & {4, 8, 16})
        if dy != [8]:
            all_dyadic_8 = False
            print(f"    DRIFT: n={nval} {g6s(G)} dyadic {dy} spec {sorted(s)}")
        if not dy:
            eg_cex += 1
        if min(s) == 5:
            g5 += 1
    girth5.append(g5)
check("all 71 C4-free graphs have dyadic profile exactly {8}", all_dyadic_8)
check("girth-5 counts at n=10,11,12 = [1,0,2]", girth5 == [1, 0, 2],
      f"got {girth5}")
check("Erdos-Gyarfas counterexamples n<=12: 0", eg_cex == 0, f"got {eg_cex}")

pub10 = ["IsPHPGQCW", "IsPHPCSCW", "IsP@PGXD_", "I{dAH?TAo", "ITaJA`IPO"]
match = 0
for code in pub10:
    P = nx.from_graph6_bytes(code.encode())
    if any(nx.is_isomorphic(P, G) for G in census[10]):
        match += 1
check("their five published n=10 g6 codes match my census 1-1", match == 5,
      f"matched {match}/5")
pet_in = any(nx.is_isomorphic(pet, G) for G in census[10])
check("Petersen in n=10 census", pet_in)
ng3_10 = sum(1 for G in census[10] if min(spec_of(G)) == 3)
check("n=10: girth-3 count is 4 (THM-456 text says '3 of the 5' -- drift?)",
      ng3_10 == 4, f"got {ng3_10}")
e16 = [G for G in census[10] if G.number_of_edges() == 16]
check("two n=10 avoiders at e=16=ex(10;C4) with degseq (4,4,3^8)",
      len(e16) == 2 and all(sorted(d for _, d in G.degree()) ==
                            [3] * 8 + [4, 4] for G in e16),
      f"got {len(e16)}")
print(f"P5 time {time.time()-t5:.1f} s")

# =====================  P6: local-move instance recounts  ======================
print()
print("=" * 78)
print("P6: local-move instance recounts (M1, M3, TWIN) + Bondy-Vince")
print("=" * 78)
t6 = time.time()
m1 = m3 = tw = 0
m1f = m3f = twf = 0
ncyc = 0
for G in conn:
    if G.number_of_nodes() < 3:
        continue
    cycles = [c for c in nx.simple_cycles(G)]
    if not cycles:
        continue
    ncyc += 1
    s = spec_of(G)
    twins = [(u, v) for u, v in combinations(sorted(G.nodes()), 2)
             if G.has_edge(u, v) and set(G[u]) - {v} == set(G[v]) - {u}]
    for cyc in cycles:
        k = len(cyc)
        on = set(cyc)
        for i in range(k):
            u, w = cyc[i], cyc[(i + 1) % k]
            for v in (set(G[u]) & set(G[w])) - on:
                m1 += 1
                if k + 1 not in s:
                    m1f += 1
        for i in range(k):
            for j in range(i + 2, k):
                if (i, j) != (0, k - 1) and G.has_edge(cyc[i], cyc[j]):
                    a = j - i
                    b = k - a
                    m3 += 1
                    if a + 1 not in s or b + 1 not in s:
                        m3f += 1
        for (u, v) in twins:
            if (u in on) != (v in on):
                tw += 1
                if k + 1 not in s:
                    twf += 1
check("971 cyclic connected graphs", ncyc == 971, f"got {ncyc}")
check("M1 instances = 111358, failures 0", m1 == 111358 and m1f == 0,
      f"got {m1}, fails {m1f}")
check("M3 instances = 143523, failures 0", m3 == 143523 and m3f == 0,
      f"got {m3}, fails {m3f}")
check("TWIN instances = 23860, failures 0", tw == 23860 and twf == 0,
      f"got {tw}, fails {twf}")

bv_fail = 0
d2only = []
for G, s in delta3_specs:
    diffs = {b - a for a in s for b in s if b > a}
    if not (1 in diffs or 2 in diffs):
        bv_fail += 1
    if 2 in diffs and 1 not in diffs:
        d2only.append(G)
check("Bondy-Vince holds on all 173 delta>=3 graphs n<=7",
      len(delta3_specs) == 173 and bv_fail == 0,
      f"{len(delta3_specs)} graphs, {bv_fail} fails")
k33_in = any(nx.is_isomorphic(G, nx.complete_bipartite_graph(3, 3))
             for G in d2only)
check("exactly 2 graphs achieve only difference 2, incl. K_{3,3}",
      len(d2only) == 2 and k33_in,
      f"got {len(d2only)}, K33 {'in' if k33_in else 'NOT in'}")
print(f"P6 time {time.time()-t6:.1f} s")

# =====================  P7: trivia  ============================================
print()
bad = [k for k in range(3, 100001)
       if not any(k <= 2 ** j <= 2 * k for j in range(2, 19))]
check("[k,2k] contains a power of 2 >= 4 for all 3<=k<=100000", not bad,
      f"first bad {bad[:3]}")

print()
print("=" * 78)
print(f"TOTAL TIME {time.time()-T0:.1f} s")
print(f"ISSUES FOUND: {len(ISSUES)}")
for name, detail in ISSUES:
    print(f"  - {name}: {detail}")
if not ISSUES:
    print("ALL RECHECKED CLAIMS CONFIRMED by second independent implementation.")
