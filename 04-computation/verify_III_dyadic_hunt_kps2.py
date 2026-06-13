"""
kind-pasteur-2026-06-09-S2 : ADVERSARIAL VERIFICATION of branch III-dyadic-hunt.

Independent recomputation of the load-bearing claims using a DIFFERENT algorithm:
the sibling's checker counts cycles via canonical min-root DFS with direction dedup
(u > v1).  Here we count every cycle 2L times (all rotations x both directions) by
enumerating closed simple walks, assert divisibility by 2L, and divide.  Any
off-by-one in cycle length, any dedup error, any factor error in their checker
would show up as a disagreement.  A third leg uses networkx.simple_cycles driven
by fresh code.

Targets:
  T1. counter self-anchors: K4, K5, cube Q3, ring C8, Petersen (built from Kneser
      K(5,2), NOT from networkx.petersen_graph).
  T2. Heawood built as Fano-plane incidence graph (NOT from networkx builtin):
      c6=28, c8=21, total 213.
  T3. McGee from own LCF builder [12,7,-7]^8; identification as THE McGee graph
      via uniqueness of the (3,7)-cage (cubic, n=24, girth 7): c7=32, c8=34 (the
      S710 correction), c9=16, c16=492, total 5608, |Aut|=32.
  T4. Tutte-Coxeter from own LCF [-13,-9,7,-7,9,13]^5 (unique (3,8)-cage):
      c8=90, total 41400.
  T5. n=28 girth-5 C8-free specimens 0/1 from dyadic_gap_hunt_kps2.out:
      cubic, connected, girth 5, c8=0, c16=731/614, 3-connected, non-isomorphic.
  T6. n=30 specimen 3: c8=0, c16=726  (< the 736 reported as the stage-2 n=30
      minimum -- documenting the [:3] reporting flaw).
  T7. dihedral Cayley graphs built FROM SCRATCH (own group arithmetic):
      D_18 refl(0,1,9): n=36, connected, girth 4, c4=9, c8=0, c16=0, explicit C32.
      D_19 refl(0,1,8): n=38, connected, girth 6, c4=0, c8=0, c16=1995.
  T8. Exoo G78 class 1 (adjacency from dyadic_exoo_g78_kps2.out): cubic, connected,
      c3=22, c4=0, c8=0, c16=0, explicit C32 and C64 certificates.
  T9. Markstrom planar n=24 class (dyadic_markstrom_rediscovery class 1): planar,
      c3=7, c4=0, c8=0, c16=228, |Aut|=3.
  T10. n=26 census sample (class 18): cubic, connected, girth 3, c4=0, c8=0, c16=454.

Output -> 05-knowledge/results/verify_III_dyadic_hunt_kps2.out
"""
import sys, time, random
from collections import deque, Counter
from itertools import combinations

OUT = []
FAILS = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def check(name, got, expect):
    ok = (got == expect)
    if not ok:
        FAILS.append(f"{name}: expected {expect}, got {got}")
    log(f"  [{'PASS' if ok else 'FAIL'}] {name}: got {got} expected {expect}")
    return ok
def save():
    with open("05-knowledge/results/verify_III_dyadic_hunt_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

# ------------------------------------------------------------------ structures
def parse_adj(s):
    """parse '0:[1, 2]; 1:[0, 2]; ...' into adjacency lists"""
    adj = {}
    for part in s.split(";"):
        part = part.strip()
        v, lst = part.split(":")
        adj[int(v)] = sorted(int(x) for x in lst.strip()[1:-1].split(","))
    n = max(adj) + 1
    return [adj[i] for i in range(n)]

def bfs_dist(adj, s):
    n = len(adj)
    INF = n + 99
    dist = [INF] * n
    dist[s] = 0
    q = deque([s])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if dist[w] > dist[u] + 1:
                dist[w] = dist[u] + 1
                q.append(w)
    return dist

def is_connected(adj):
    return max(bfs_dist(adj, 0)) < len(adj)

# ------------------------------------------------- INDEPENDENT cycle counter
def cyc(adj, L):
    """Count simple cycles of length L by enumerating ALL closed simple walks
    of length L (every root, every direction) and dividing by 2L.
    Completely different dedup logic from the sibling's canonical-root counter."""
    n = len(adj)
    if L < 3 or L > n:
        return 0
    adjset = [set(a) for a in adj]
    total = 0
    for root in range(n):
        dist = bfs_dist(adj, root)
        # path positions 0..L-1; vertex at position p needs dist <= L - p
        on = [False] * n
        on[root] = True
        # stack frames: (vertex, position, iterator)
        stack = [(root, 0, iter(adj[root]))]
        while stack:
            u, p, it = stack[-1]
            moved = False
            for w in it:
                if p == L - 1:
                    continue  # closure handled when frame was pushed
                if on[w]:
                    continue
                if dist[w] > L - (p + 1):
                    continue
                if p + 1 == L - 1:
                    # w is the last vertex: count closure directly, do not recurse
                    if root in adjset[w]:
                        total += 1
                    continue
                on[w] = True
                stack.append((w, p + 1, iter(adj[w])))
                moved = True
                break
            if not moved:
                on[u] = False
                stack.pop()
    assert total % (2 * L) == 0, f"closed-walk total {total} not divisible by {2*L}"
    return total // (2 * L)

def min_cycle_len(adj, upto):
    for L in range(3, upto + 1):
        if cyc(adj, L) > 0:
            return L
    return None

def find_cycle_len(adj, L, seed=0, budget=60_000_000):
    """Return an explicit simple cycle of length L as a vertex list, or None/'budget'.
    Randomized-order DFS with distance pruning, early exit."""
    n = len(adj)
    if L < 3 or L > n:
        return None
    rng = random.Random(seed)
    adjset = [set(a) for a in adj]
    expansions = 0
    roots = list(range(n))
    rng.shuffle(roots)
    for root in roots:
        dist = bfs_dist(adj, root)
        on = [False] * n
        on[root] = True
        path = [root]
        def shuffled(u):
            a = list(adj[u])
            rng.shuffle(a)
            return iter(a)
        stack = [(root, 0, shuffled(root))]
        while stack:
            expansions += 1
            if expansions > budget:
                return "budget"
            u, p, it = stack[-1]
            moved = False
            for w in it:
                if on[w] or dist[w] > L - (p + 1):
                    continue
                if p + 1 == L - 1:
                    if root in adjset[w]:
                        return path + [w]
                    continue
                on[w] = True
                path.append(w)
                stack.append((w, p + 1, shuffled(w)))
                moved = True
                break
            if not moved:
                on[u] = False
                stack.pop()
                path.pop()
        # root exhausted; try next root (cycle through root absent)
    return None

def verify_cycle(adj, cycle, L):
    """certificate check: L distinct vertices, consecutive adjacency, closure"""
    if cycle is None or cycle == "budget" or len(cycle) != L or len(set(cycle)) != L:
        return False
    return all(cycle[(i + 1) % L] in adj[cycle[i]] for i in range(L))

# ------------------------------------------------------------- networkx leg
import networkx as nx
def to_nx(adj):
    G = nx.Graph()
    G.add_nodes_from(range(len(adj)))
    for u in range(len(adj)):
        for w in adj[u]:
            if u < w:
                G.add_edge(u, w)
    return G

def nx_spectrum(adj):
    return Counter(len(c) for c in nx.simple_cycles(to_nx(adj)))

def aut_order(adj):
    G = to_nx(adj)
    return sum(1 for _ in nx.vf2pp_all_isomorphisms(G, G, node_label=None))

# --------------------------------------------------------- own constructions
def lcf_graph(n, lcf, reps):
    adj = [[] for _ in range(n)]
    def add(a, b):
        if b not in adj[a]:
            adj[a].append(b)
            adj[b].append(a)
    for i in range(n):
        add(i, (i + 1) % n)
    seq = lcf * reps
    assert len(seq) == n
    for i, d in enumerate(seq):
        add(i, (i + d) % n)
    return adj

def kneser_petersen():
    verts = list(combinations(range(5), 2))
    idx = {v: i for i, v in enumerate(verts)}
    adj = [[] for _ in range(10)]
    for a, b in combinations(verts, 2):
        if not (set(a) & set(b)):
            adj[idx[a]].append(idx[b])
            adj[idx[b]].append(idx[a])
    return adj

def fano_heawood():
    lines = [tuple(sorted(((0 + i) % 7, (1 + i) % 7, (3 + i) % 7))) for i in range(7)]
    adj = [[] for _ in range(14)]
    for p in range(7):
        for li, line in enumerate(lines):
            if p in line:
                adj[p].append(7 + li)
                adj[7 + li].append(p)
    return adj

def dihedral_cayley(m, refls):
    """Cay(D_m, {s r^j : j in refls}). Elements: (0,a)=r^a, (1,a)=s r^a.
    Own group arithmetic: r^a * s r^b = s r^{b-a};  s r^a * s r^b = r^{b-a}."""
    def mul(g, h):
        t1, a = g
        t2, b = h
        if t1 == 0 and t2 == 0:
            return (0, (a + b) % m)
        if t1 == 0 and t2 == 1:
            return (1, (b - a) % m)
        if t1 == 1 and t2 == 0:
            return (1, (a + b) % m)
        return (0, (b - a) % m)
    elems = [(0, a) for a in range(m)] + [(1, a) for a in range(m)]
    idx = {e: i for i, e in enumerate(elems)}
    S = [(1, j) for j in refls]
    for s in S:  # involutions => S = S^{-1}
        assert mul(s, s) == (0, 0)
    adj = [[] for _ in range(2 * m)]
    for g in elems:
        for s in S:
            h = mul(g, s)
            if idx[h] not in adj[idx[g]]:
                adj[idx[g]].append(idx[h])
                adj[idx[h]].append(idx[g])
    return adj

# ================================================================== specimens
N28_SPEC0 = parse_adj("0:[2, 13, 23]; 1:[6, 7, 24]; 2:[0, 10, 19]; 3:[11, 14, 17]; 4:[13, 17, 27]; 5:[7, 9, 11]; 6:[1, 16, 18]; 7:[1, 5, 8]; 8:[7, 16, 25]; 9:[5, 14, 21]; 10:[2, 15, 17]; 11:[3, 5, 18]; 12:[19, 22, 25]; 13:[0, 4, 24]; 14:[3, 9, 26]; 15:[10, 20, 24]; 16:[6, 8, 20]; 17:[3, 4, 10]; 18:[6, 11, 27]; 19:[2, 12, 20]; 20:[15, 16, 19]; 21:[9, 23, 25]; 22:[12, 23, 26]; 23:[0, 21, 22]; 24:[1, 13, 15]; 25:[8, 12, 21]; 26:[14, 22, 27]; 27:[4, 18, 26]")
N28_SPEC1 = parse_adj("0:[1, 5, 9]; 1:[0, 12, 21]; 2:[6, 7, 10]; 3:[6, 8, 13]; 4:[11, 21, 24]; 5:[0, 20, 26]; 6:[2, 3, 14]; 7:[2, 8, 23]; 8:[3, 7, 18]; 9:[0, 22, 25]; 10:[2, 18, 24]; 11:[4, 16, 22]; 12:[1, 13, 19]; 13:[3, 12, 24]; 14:[6, 15, 21]; 15:[14, 16, 23]; 16:[11, 15, 17]; 17:[16, 25, 27]; 18:[8, 10, 26]; 19:[12, 20, 22]; 20:[5, 19, 25]; 21:[1, 4, 14]; 22:[9, 11, 19]; 23:[7, 15, 27]; 24:[4, 10, 13]; 25:[9, 17, 20]; 26:[5, 18, 27]; 27:[17, 23, 26]")
N30_SPEC3 = parse_adj("0:[1, 8, 13]; 1:[0, 15, 21]; 2:[9, 22, 23]; 3:[10, 15, 29]; 4:[8, 14, 23]; 5:[9, 27, 29]; 6:[9, 24, 26]; 7:[16, 26, 29]; 8:[0, 4, 12]; 9:[2, 5, 6]; 10:[3, 12, 21]; 11:[19, 22, 28]; 12:[8, 10, 20]; 13:[0, 18, 25]; 14:[4, 25, 28]; 15:[1, 3, 20]; 16:[7, 17, 21]; 17:[16, 19, 25]; 18:[13, 23, 28]; 19:[11, 17, 26]; 20:[12, 15, 27]; 21:[1, 10, 16]; 22:[2, 11, 24]; 23:[2, 4, 18]; 24:[6, 22, 27]; 25:[13, 14, 17]; 26:[6, 7, 19]; 27:[5, 20, 24]; 28:[11, 14, 18]; 29:[3, 5, 7]")
G78_CLASS1 = parse_adj("0:[1, 2, 63]; 1:[0, 3, 4]; 2:[0, 5, 6]; 3:[1, 4, 13]; 4:[1, 3, 5]; 5:[2, 4, 6]; 6:[2, 5, 41]; 7:[8, 9, 20]; 8:[7, 10, 11]; 9:[7, 12, 13]; 10:[8, 11, 42]; 11:[8, 10, 12]; 12:[9, 11, 13]; 13:[3, 9, 12]; 14:[15, 16, 55]; 15:[14, 17, 18]; 16:[14, 19, 20]; 17:[15, 18, 24]; 18:[15, 17, 19]; 19:[16, 18, 20]; 20:[7, 16, 19]; 21:[22, 23, 77]; 22:[21, 24, 25]; 23:[21, 26, 27]; 24:[17, 22, 25]; 25:[22, 24, 26]; 26:[23, 25, 27]; 27:[23, 26, 59]; 28:[29, 30, 70]; 29:[28, 31, 32]; 30:[28, 33, 34]; 31:[29, 32, 52]; 32:[29, 31, 33]; 33:[30, 32, 34]; 34:[30, 33, 45]; 35:[36, 37, 62]; 36:[35, 38, 39]; 37:[35, 40, 41]; 38:[36, 39, 49]; 39:[36, 38, 40]; 40:[37, 39, 41]; 41:[6, 37, 40]; 42:[10, 43, 44]; 43:[42, 45, 46]; 44:[42, 47, 48]; 45:[34, 43, 46]; 46:[43, 45, 47]; 47:[44, 46, 48]; 48:[44, 47, 56]; 49:[38, 50, 51]; 50:[49, 52, 53]; 51:[49, 54, 55]; 52:[31, 50, 53]; 53:[50, 52, 54]; 54:[51, 53, 55]; 55:[14, 51, 54]; 56:[48, 57, 58]; 57:[56, 59, 60]; 58:[56, 61, 62]; 59:[27, 57, 60]; 60:[57, 59, 61]; 61:[58, 60, 62]; 62:[35, 58, 61]; 63:[0, 64, 65]; 64:[63, 66, 67]; 65:[63, 68, 69]; 66:[64, 67, 77]; 67:[64, 66, 68]; 68:[65, 67, 69]; 69:[65, 68, 73]; 70:[28, 71, 72]; 71:[70, 73, 74]; 72:[70, 75, 76]; 73:[69, 71, 74]; 74:[71, 73, 75]; 75:[72, 74, 76]; 76:[72, 75, 77]; 77:[21, 66, 76]")
MARK24_PLANAR = parse_adj("0:[6, 17, 22]; 1:[2, 18, 23]; 2:[1, 13, 19]; 3:[8, 15, 16]; 4:[13, 14, 17]; 5:[9, 10, 11]; 6:[0, 10, 22]; 7:[12, 15, 20]; 8:[3, 9, 15]; 9:[5, 8, 18]; 10:[5, 6, 11]; 11:[5, 10, 23]; 12:[7, 16, 21]; 13:[2, 4, 17]; 14:[4, 19, 20]; 15:[3, 7, 8]; 16:[3, 12, 21]; 17:[0, 4, 13]; 18:[1, 9, 23]; 19:[2, 14, 20]; 20:[7, 14, 19]; 21:[12, 16, 22]; 22:[0, 6, 21]; 23:[1, 11, 18]")
N26_CLASS18 = parse_adj("0:[8, 13, 22]; 1:[5, 19, 24]; 2:[3, 10, 11]; 3:[2, 5, 16]; 4:[9, 14, 19]; 5:[1, 3, 6]; 6:[5, 9, 23]; 7:[8, 16, 17]; 8:[0, 7, 25]; 9:[4, 6, 15]; 10:[2, 15, 17]; 11:[2, 18, 24]; 12:[18, 20, 25]; 13:[0, 14, 15]; 14:[4, 13, 23]; 15:[9, 10, 13]; 16:[3, 7, 18]; 17:[7, 10, 23]; 18:[11, 12, 16]; 19:[1, 4, 22]; 20:[12, 21, 25]; 21:[20, 22, 24]; 22:[0, 19, 21]; 23:[6, 14, 17]; 24:[1, 11, 21]; 25:[8, 12, 20]")

def main():
    t0 = time.time()
    log("=" * 100)
    log("ADVERSARIAL VERIFICATION: branch III-dyadic-hunt  (verify_III_dyadic_hunt_kps2)")
    log("Independent closed-walk/2L cycle counter + fresh constructions + networkx leg")
    log("=" * 100)

    # ---- T1 anchors
    log("\nT1. counter self-anchors")
    K4 = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]
    K5 = [[j for j in range(5) if j != i] for i in range(5)]
    Q3 = [[i ^ 1, i ^ 2, i ^ 4] for i in range(8)]
    C8 = [[(i - 1) % 8, (i + 1) % 8] for i in range(8)]
    check("K4 c3", cyc(K4, 3), 4)
    check("K4 c4", cyc(K4, 4), 3)
    check("K5 c3", cyc(K5, 3), 10)
    check("K5 c4", cyc(K5, 4), 15)
    check("K5 c5", cyc(K5, 5), 12)
    check("cube c4", cyc(Q3, 4), 6)
    check("cube c6", cyc(Q3, 6), 16)
    check("cube c8 (Ham)", cyc(Q3, 8), 6)
    check("ring C8 c8", cyc(C8, 8), 1)
    PET = kneser_petersen()
    check("Petersen(Kneser) girth", min_cycle_len(PET, 10), 5)
    check("Petersen c5", cyc(PET, 5), 12)
    check("Petersen c6", cyc(PET, 6), 10)
    check("Petersen c7", cyc(PET, 7), 0)
    check("Petersen c8", cyc(PET, 8), 15)
    check("Petersen c9", cyc(PET, 9), 20)
    check("Petersen c10 (non-Ham)", cyc(PET, 10), 0)
    check("Petersen total (nx leg)", sum(nx_spectrum(PET).values()), 57)
    save()

    # ---- T2 Heawood
    log("\nT2. Heawood as Fano incidence graph (own construction)")
    HEA = fano_heawood()
    check("Heawood cubic", all(len(a) == 3 for a in HEA), True)
    check("Heawood girth", min_cycle_len(HEA, 14), 6)
    check("Heawood c6", cyc(HEA, 6), 28)
    check("Heawood c8", cyc(HEA, 8), 21)
    spec = nx_spectrum(HEA)
    check("Heawood total (nx leg)", sum(spec.values()), 213)
    check("Heawood spectrum (nx leg)", dict(sorted(spec.items())),
          {6: 28, 8: 21, 10: 84, 12: 56, 14: 24})
    save()

    # ---- T3 McGee
    log("\nT3. McGee via own LCF builder; identity = unique (3,7)-cage")
    MCG = lcf_graph(24, [12, 7, -7], 8)
    check("McGee cubic", all(len(a) == 3 for a in MCG), True)
    check("McGee girth (=> IS the (3,7)-cage)", min_cycle_len(MCG, 24), 7)
    check("McGee c7", cyc(MCG, 7), 32)
    check("McGee c8  [THE S710 CORRECTION]", cyc(MCG, 8), 34)
    check("McGee c9", cyc(MCG, 9), 16)
    check("McGee c16", cyc(MCG, 16), 492)
    check("McGee |Aut|", aut_order(MCG), 32)
    mspec = nx_spectrum(MCG)
    check("McGee total (nx leg)", sum(mspec.values()), 5608)
    check("McGee c8 (nx leg)", mspec[8], 34)
    save()

    # ---- T4 Tutte-Coxeter
    log("\nT4. Tutte-Coxeter via own LCF builder; identity = unique (3,8)-cage")
    TC = lcf_graph(30, [-13, -9, 7, -7, 9, 13], 5)
    check("TutteCoxeter cubic", all(len(a) == 3 for a in TC), True)
    check("TutteCoxeter girth (=> IS the (3,8)-cage)", min_cycle_len(TC, 8), 8)
    check("TutteCoxeter c8", cyc(TC, 8), 90)
    tspec = nx_spectrum(TC)
    check("TutteCoxeter total (nx leg)", sum(tspec.values()), 41400)
    save()

    # ---- T5 n=28 specimens
    log("\nT5. n=28 girth-5 C8-free specimens (the 'min order = 28' witnesses)")
    for tag, A, c5e, c6e, c7e, c16e in (("spec0", N28_SPEC0, 4, 12, 8, 731),
                                        ("spec1", N28_SPEC1, 5, 12, 12, 614)):
        check(f"n28 {tag} cubic", all(len(a) == 3 for a in A), True)
        check(f"n28 {tag} connected", is_connected(A), True)
        check(f"n28 {tag} girth", min_cycle_len(A, 8), 5)
        check(f"n28 {tag} c5", cyc(A, 5), c5e)
        check(f"n28 {tag} c6", cyc(A, 6), c6e)
        check(f"n28 {tag} c7", cyc(A, 7), c7e)
        check(f"n28 {tag} c8 == 0", cyc(A, 8), 0)
        check(f"n28 {tag} c16", cyc(A, 16), c16e)
        check(f"n28 {tag} 3-connected", nx.node_connectivity(to_nx(A)), 3)
    check("n28 spec0 !iso spec1 (>=2 classes)",
          nx.is_isomorphic(to_nx(N28_SPEC0), to_nx(N28_SPEC1)), False)
    save()

    # ---- T6 n=30 spec3 (frontier reporting flaw)
    log("\nT6. n=30 specimen 3 (excluded by specimens[n][:3] in stage 2)")
    check("n30 spec3 girth", min_cycle_len(N30_SPEC3, 8), 5)
    check("n30 spec3 c8 == 0", cyc(N30_SPEC3, 8), 0)
    c16_30 = cyc(N30_SPEC3, 16)
    check("n30 spec3 c16 (vs reported stage-2 'min' 736)", c16_30, 726)
    log(f"  => stage-2 frontier value 736 at n=30 is NOT the min over their own specimens"
        f" (spec3 has {c16_30}); monotonic-growth claim still holds (614 < {c16_30}).")
    save()

    # ---- T7 dihedral Cayley from scratch
    log("\nT7. dihedral 3-reflection Cayley graphs, own group arithmetic")
    D18 = dihedral_cayley(18, (0, 1, 9))
    check("D18(0,1,9) n", len(D18), 36)
    check("D18(0,1,9) cubic", all(len(a) == 3 for a in D18), True)
    check("D18(0,1,9) connected", is_connected(D18), True)
    check("D18(0,1,9) girth", min_cycle_len(D18, 8), 4)
    check("D18(0,1,9) c4", cyc(D18, 4), 9)
    check("D18(0,1,9) c8 == 0", cyc(D18, 8), 0)
    check("D18(0,1,9) c16 == 0", cyc(D18, 16), 0)
    cyc32 = find_cycle_len(D18, 32, seed=1)
    check("D18(0,1,9) C32 certificate", verify_cycle(D18, cyc32, 32), True)
    if verify_cycle(D18, cyc32, 32):
        log(f"    C32 witness: {cyc32}")
    log("  => dyadic spectrum of D18(0,1,9) is exactly {4,32} (n=36<64 excludes C64): CONFIRMED")
    D19 = dihedral_cayley(19, (0, 1, 8))
    check("D19(0,1,8) n", len(D19), 38)
    check("D19(0,1,8) connected", is_connected(D19), True)
    check("D19(0,1,8) girth", min_cycle_len(D19, 8), 6)
    check("D19(0,1,8) c4 == 0", cyc(D19, 4), 0)
    check("D19(0,1,8) c8 == 0", cyc(D19, 8), 0)
    check("D19(0,1,8) c16", cyc(D19, 16), 1995)
    save()

    # ---- T8 G78
    log("\nT8. Exoo G78 reconstruction, class 1")
    A = G78_CLASS1
    check("G78 n", len(A), 78)
    check("G78 cubic", all(len(a) == 3 for a in A), True)
    check("G78 connected", is_connected(A), True)
    check("G78 c3", cyc(A, 3), 22)
    check("G78 c4 == 0", cyc(A, 4), 0)
    check("G78 c8 == 0", cyc(A, 8), 0)
    t = time.time()
    c16_g78 = cyc(A, 16)
    check("G78 c16 == 0", c16_g78, 0)
    log(f"    (c16 recount took {time.time()-t:.1f}s)")
    w32 = find_cycle_len(A, 32, seed=2)
    check("G78 C32 certificate", verify_cycle(A, w32, 32), True)
    if verify_cycle(A, w32, 32):
        log(f"    C32 witness: {w32}")
    w64 = find_cycle_len(A, 64, seed=3)
    check("G78 C64 certificate", verify_cycle(A, w64, 64), True)
    if verify_cycle(A, w64, 64):
        log(f"    C64 witness: {w64}")
    save()

    # ---- T9 Markstrom planar n=24
    log("\nT9. Markstrom planar n=24 class (triple-lock graph)")
    A = MARK24_PLANAR
    check("M24 cubic", all(len(a) == 3 for a in A), True)
    check("M24 connected", is_connected(A), True)
    check("M24 planar", nx.check_planarity(to_nx(A))[0], True)
    check("M24 c3", cyc(A, 3), 7)
    check("M24 c4 == 0", cyc(A, 4), 0)
    check("M24 c8 == 0", cyc(A, 8), 0)
    check("M24 c16", cyc(A, 16), 228)
    check("M24 |Aut|", aut_order(A), 3)
    save()

    # ---- T10 n=26 census sample
    log("\nT10. n=26 census sample: class 18")
    A = N26_CLASS18
    check("n26c18 cubic", all(len(a) == 3 for a in A), True)
    check("n26c18 connected", is_connected(A), True)
    check("n26c18 girth", min_cycle_len(A, 8), 3)
    check("n26c18 c3", cyc(A, 3), 1)
    check("n26c18 c4 == 0", cyc(A, 4), 0)
    check("n26c18 c8 == 0", cyc(A, 8), 0)
    check("n26c18 c16", cyc(A, 16), 454)
    save()

    # ---- verdict
    log("\n" + "=" * 100)
    if FAILS:
        log(f"VERDICT: {len(FAILS)} FAILURES:")
        for f in FAILS:
            log("  !! " + f)
    else:
        log("VERDICT: ALL CHECKS PASSED -- every recomputed number matches the sibling's claims.")
    log(f"total time {time.time()-t0:.0f}s")
    save()

if __name__ == "__main__":
    main()
