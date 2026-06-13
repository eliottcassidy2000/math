"""
kind-pasteur-2026-06-09-S2 : BRANCH III (dyadic-gap hunt, HYP-2359 / Erdos-Gyarfas)
RECONSTRUCT EXOO'S G78 (arXiv:1403.5636 Sec. 3): the smallest KNOWN cubic graph with
no 2^m-cycle for m <= 4 (f(4) <= 78; lower bound 54 = unpublished Markstrom).

Construction: Petersen -> replace one vertex by K3 (= G12, 12 vertices, cubic)
-> replace 11 of the 12 vertices by Exoo's H7 gadget (one vertex left bare).
The u/v/w attachment wiring per copy is given only by a figure in the paper, so we
search: choice of bare vertex (12) x per-copy attachment permutation (6^11) via
simulated annealing minimizing (#C8, then #C16). c4=0 is automatic (H7 is C4-free
and any projected cycle has blowup length >= 7).

NEW DATA the paper does not give: the C32 / C64 status of the resulting graph(s).
At n=78 < 128 an Erdos-Gyarfas counterexample must avoid {4,8,16,32,64}; we check
whether 32- and 64-cycles exist in every {4,8,16}-free wiring found.

Output -> 05-knowledge/results/dyadic_exoo_g78_kps2.out
"""
import sys, os, time, random
from itertools import permutations
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from dyadic_cycle_checker_kps2 import count_cycles_len, girth, edges_to_adj
from dyadic_gap_hunt_kps2 import count_cycles_capped, is_connected

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def save():
    with open("05-knowledge/results/dyadic_exoo_g78_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

# ---- H7 (Exoo Fig. 2): u,a,b,v,c,d,w = 0..6 ; attachments u=0 (deg2), v=3, w=6 ----
H7_EDGES = [(0,1),(0,2),(1,3),(1,4),(3,4),(2,5),(2,6),(5,6),(4,5)]
ATT = [0, 3, 6]   # slot 0 = u, slot 1 = v, slot 2 = w
PERMS = list(permutations(range(3)))

def petersen_adj():
    # outer 0-4, inner 5-9, GP(5,2)
    edges = []
    for i in range(5):
        edges.append((i, (i + 1) % 5))
        edges.append((i, 5 + i))
        edges.append((5 + i, 5 + (i + 2) % 5))
    return edges

def build_G12():
    """Petersen with vertex 0 replaced by a triangle {10,11,12}->relabeled 0..11."""
    pe = petersen_adj()
    nbr0 = sorted({b for (a, b) in pe if a == 0} | {a for (a, b) in pe if b == 0})
    edges = [(a, b) for (a, b) in pe if a != 0 and b != 0]
    tri = [10, 11, 12]
    for t, x in zip(tri, nbr0):
        edges.append((t, x))
    edges += [(10, 11), (11, 12), (10, 12)]
    # relabel to 0..11 (old vertices 1..9 -> 0..8, 10,11,12 -> 9,10,11)
    relab = {}
    for v in range(1, 10):
        relab[v] = v - 1
    relab[10], relab[11], relab[12] = 9, 10, 11
    return [(relab[a], relab[b]) for (a, b) in edges]

G12_EDGES = build_G12()
G12_NBRS = [[] for _ in range(12)]
for (a, b) in G12_EDGES:
    G12_NBRS[a].append(b)
    G12_NBRS[b].append(a)
for x in G12_NBRS:
    x.sort()

def build_G78(bare, wire):
    """bare = vertex of G12 left as-is; wire[x] = perm index for each replaced x.
    Replaced vertex x occupies 7*pos[x] .. 7*pos[x]+6 ; bare vertex is vertex 77."""
    repl = [v for v in range(12) if v != bare]
    pos = {v: i for i, v in enumerate(repl)}
    N = 7 * 11 + 1
    edges = []
    for v in repl:
        base = 7 * pos[v]
        for (a, b) in H7_EDGES:
            edges.append((base + a, base + b))
    def port(v, y):
        """vertex of copy v facing neighbor y"""
        slot = G12_NBRS[v].index(y)
        p = PERMS[wire[pos[v]]]
        return 7 * pos[v] + ATT[p[slot]]
    done = set()
    for (a, b) in G12_EDGES:
        e = frozenset((a, b))
        if e in done:
            continue
        done.add(e)
        pa = 77 if a == bare else port(a, b)
        pb = 77 if b == bare else port(b, a)
        edges.append((pa, pb))
    return edges_to_adj(N, edges)

def main():
    t0 = time.time()
    log("=" * 100)
    log("EXOO G78 RECONSTRUCTION (arXiv:1403.5636): SA over bare-vertex x wiring space")
    log("=" * 100)
    log(f"G12 edges ({len(G12_EDGES)}): {sorted(tuple(sorted(e)) for e in G12_EDGES)}")
    adjT = edges_to_adj(12, G12_EDGES)
    log(f"G12: cubic={all(len(x)==3 for x in adjT)} girth={girth(adjT)} "
        f"c3={count_cycles_len(adjT,3)} c5={count_cycles_len(adjT,5)} c6={count_cycles_len(adjT,6)}")
    save()

    rng = random.Random(20260609)
    found = []
    BIG = 10 ** 6
    for bare in range(12):
        # two-stage SA on wiring: cost = BIG*c8 + c16 (c16 only evaluated when c8=0,
        # with cap-based early exit against the incumbent)
        best = None
        best_w = None
        for restart in range(4):
            wire = [rng.randrange(6) for _ in range(11)]
            def cost(w, incumbent):
                adj = build_G78(bare, w)
                c8 = count_cycles_len(adj, 8)
                if c8 > 0:
                    return BIG * c8
                cap = None if incumbent is None else max(0, incumbent) + 50
                c16 = count_cycles_capped(adj, 16, cap=cap)
                return c16
            cur = cost(wire, None)
            if best is None or cur < best:
                best, best_w = cur, wire[:]
            T = 30.0
            stall = 0
            for it in range(1500):
                w2 = wire[:]
                w2[rng.randrange(11)] = rng.randrange(6)
                s2 = cost(w2, best if best < BIG else None)
                if s2 <= cur or rng.random() < pow(2.718, -(s2 - cur) / max(T, 1e-9)):
                    wire, cur = w2, s2
                if cur < best:
                    best, best_w = cur, wire[:]
                    stall = 0
                else:
                    stall += 1
                T *= 0.997
                if best == 0 or stall > 500:
                    break
            if best == 0:
                break
        adj = build_G78(bare, best_w)
        c8 = count_cycles_len(adj, 8)
        c16 = count_cycles_capped(adj, 16)
        line = f"  bare={bare:2d}: best (c8,c16) = ({c8},{c16}) ({time.time()-t0:.0f}s)"
        if c8 == 0 and c16 == 0:
            found.append((bare, best_w[:], adj))
            line += "  *** {4,8,16}-FREE ***"
        log(line)
        save()

    log("")
    log("=" * 100)
    log(f"{{4,8,16}}-free wirings found: {len(found)}")
    log("=" * 100)
    import networkx as nx
    classes = []
    for (bare, w, adj) in found:
        G = nx.Graph((u, v) for u, vs in enumerate(adj) for v in vs if u < v)
        if not any(nx.is_isomorphic(G, H) for H in classes):
            classes.append(G)
    log(f"distinct iso classes among them: {len(classes)}")
    for i, G in enumerate(classes):
        adj = edges_to_adj(78, list(G.edges()))
        g = girth(adj)
        conn = is_connected(adj)
        c3 = count_cycles_len(adj, 3)
        c4 = count_cycles_len(adj, 4)
        c8 = count_cycles_len(adj, 8)
        c16 = count_cycles_capped(adj, 16)
        # existence checks for C32 / C64 (cap=0 -> early exit on first cycle found)
        t = time.time()
        c32x = count_cycles_capped(adj, 32, cap=0, budget=80_000_000)
        t32 = time.time() - t
        t = time.time()
        c64x = count_cycles_capped(adj, 64, cap=0, budget=80_000_000)
        t64 = time.time() - t
        def word(x, tt):
            return ("EXISTS" if x and x > 0 else ("inconclusive(budget)" if x == -1 else "NONE")) + f" ({tt:.0f}s)"
        log(f"class {i}: n=78 cubic conn={conn} girth={g} c3={c3} c4={c4} c8={c8} c16={c16}")
        log(f"   C32: {word(c32x, t32)}   C64: {word(c64x, t64)}")
        log(f"   NEW DATA vs Exoo paper: dyadic fate of G78-type graphs at 32/64.")
        log("   adj: " + "; ".join(f"{v}:{sorted(adj[v])}" for v in range(78)))
        save()
    if not found:
        log("NO {4,8,16}-free wiring found by SA -- report c8/c16 floors above (negative result).")
    log(f"total time {time.time()-t0:.0f}s")
    save()

if __name__ == "__main__":
    main()
