"""
kind-pasteur-2026-06-09-S2 : BRANCH III (dyadic-gap hunt, HYP-2359 / Erdos-Gyarfas)
Part 1: FAST exact cycle-length checker + verification on known graphs.

count_cycles_len(adj, L): exact number of simple cycles of length L.
  Canonical form: smallest label s on the cycle is the DFS root; only vertices > s
  may appear on the path; direction dedup via (second vertex) < (last vertex).
  Pruning: BFS distance back to s within the {v >= s} induced subgraph.

Verification targets (all independent literature values):
  - Petersen (n=10, girth 5): cycle spectrum c5=12, c6=10, c7=0, c8=15, c9=20;
    total cycles = 57 = psi(6) (Markstrom 2004 Table 1: Petersen is psi-extremal, psi=57).
  - Heawood (n=14, girth 6): total cycles = 213 = psi(8) (Markstrom Table 1).
  - McGee (n=24, girth 7): total cycles = 5608 = psi(13) lower bound attained by McGee
    (Markstrom Table 1). KEY QUESTION: does McGee have an 8-cycle? (S710 says McGee -> C16,
    implying c8=0; if c8(McGee)=0 then McGee is one of Markstrom's 4 cubic n=24 graphs
    with no C4 and no C8 -- Table 3 of Markstrom 2004.)
  - Tutte-Coxeter (n=30, girth 8): total cycles = 41400 = psi(16) bound (Markstrom Table 1).
  - Coxeter graph (n=28, girth 7): does it avoid C8? (would put it among the 251 at n=28).

Output saved to 05-knowledge/results/dyadic_cycle_checker_kps2.out
"""
import sys, time
from collections import deque

OUT = []
def log(s=""):
    print(s)
    OUT.append(str(s))

# ---------------------------------------------------------------- core checker
def bfs_dist_restricted(adj, s, n):
    """BFS distances from s in subgraph induced on {v >= s}. INF = n+1."""
    INF = n + 1
    dist = [INF] * n
    dist[s] = 0
    q = deque([s])
    while q:
        u = q.popleft()
        du = dist[u]
        for w in adj[u]:
            if w >= s and dist[w] == INF:
                dist[w] = du + 1
                q.append(w)
    return dist

def count_cycles_len(adj, L, early_exit=False, budget=None):
    """Exact count of simple cycles of length L in graph given by adjacency lists.
    If early_exit: return 1 as soon as one cycle is found.
    budget: max number of DFS node-expansions; returns -1 if exceeded (inconclusive)."""
    n = len(adj)
    if L < 3 or L > n:
        return 0
    count = 0
    expansions = 0
    for s in range(n):
        dist = bfs_dist_restricted(adj, s, n)
        # quick reject: need a cycle through s of length L => some neighbor pair usable
        on_path = [False] * n
        on_path[s] = True
        # iterative DFS: stack of (vertex, depth, neighbor-iterator-index, first_neighbor)
        # path[1] = first vertex after s; closing vertex must be > path[1]
        nbrs_s = [w for w in adj[s] if w > s and dist[w] <= L - 1]
        for v1 in nbrs_s:
            # DFS from v1 at depth 1
            on_path[v1] = True
            stack = [(v1, 1, iter(adj[v1]))]
            while stack:
                u, depth, it = stack[-1]
                if budget is not None:
                    expansions += 1
                    if expansions > budget:
                        # unwind
                        for (x, d, _) in stack:
                            on_path[x] = False
                        on_path[s] = False
                        return -1
                advanced = False
                for w in it:
                    if depth == L - 1:
                        # close the cycle: w must equal s, and u > v1 for dedup
                        # (we scan the iterator fully at the last level)
                        if w == s and u > v1:
                            count += 1
                            if early_exit:
                                return 1
                        continue
                    if w <= s or on_path[w]:
                        continue
                    if dist[w] > L - depth - 1:
                        continue
                    # descend
                    on_path[w] = True
                    stack.append((w, depth + 1, iter(adj[w])))
                    advanced = True
                    break
                if not advanced:
                    on_path[u] = False
                    stack.pop()
            on_path[v1] = False
        on_path[s] = False
    return count

def girth(adj):
    """Girth via BFS from every vertex."""
    n = len(adj)
    best = n + 1
    for s in range(n):
        dist = [-1] * n
        par = [-1] * n
        dist[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            if 2 * dist[u] + 1 >= best:
                continue
            for w in adj[u]:
                if dist[w] == -1:
                    dist[w] = dist[u] + 1
                    par[w] = u
                    q.append(w)
                elif par[u] != w and par[w] != u:
                    c = dist[u] + dist[w] + 1
                    if c < best:
                        best = c
    return best if best <= n else 0

def cycle_spectrum(adj, maxlen=None):
    n = len(adj)
    if maxlen is None:
        maxlen = n
    return {L: count_cycles_len(adj, L) for L in range(3, maxlen + 1)}

def edges_to_adj(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    return adj

# ---------------------------------------------------------------- known graphs
import networkx as nx
def nx_to_adj(G):
    G = nx.convert_node_labels_to_integers(G, ordering="sorted")
    n = G.number_of_nodes()
    return edges_to_adj(n, G.edges())

GRAPHS = {
    "Petersen      (n=10)": nx_to_adj(nx.petersen_graph()),
    "Heawood       (n=14)": nx_to_adj(nx.heawood_graph()),
    "MoebiusKantor (n=16)": nx_to_adj(nx.moebius_kantor_graph()),
    "Pappus        (n=18)": nx_to_adj(nx.pappus_graph()),
    "Desargues     (n=20)": nx_to_adj(nx.desargues_graph()),
    "Dodecahedron  (n=20)": nx_to_adj(nx.dodecahedral_graph()),
    "McGee         (n=24)": nx_to_adj(nx.LCF_graph(24, [12, 7, -7], 8)),
    "Nauru         (n=24)": nx_to_adj(nx.LCF_graph(24, [5, -9, 7, -7, 9, -5], 4)),
    "Coxeter       (n=28)": None,  # built below
    "TutteCoxeter  (n=30)": nx_to_adj(nx.LCF_graph(30, [-13, -9, 7, -7, 9, 13], 5)),
}

# Coxeter graph (28 vertices, girth 7): standard construction.
# Vertices: hub-free: 4 fans: a central 7-cycle? Use the standard description:
# 24 vertices in three 'rings' + 4? Simplest reliable source: networkx has no builtin;
# use the known construction: vertices Z7 x {0,1,2,3}; described by circulant layers:
#   layer 0: 7-cycle skips 1 (i,0)-(i+1,0)
#   layer 1: skips 2  (i,1)-(i+2,1)
#   layer 2: skips 3  (i,2)-(i+3,2)
#   spokes: (i,3)-(i,0),(i,3)-(i,1),(i,3)-(i,2)
cox_edges = []
for i in range(7):
    cox_edges.append(((i, 0), ((i + 1) % 7, 0)))
    cox_edges.append(((i, 1), ((i + 2) % 7, 1)))
    cox_edges.append(((i, 2), ((i + 3) % 7, 2)))
    cox_edges.append(((i, 3), (i, 0)))
    cox_edges.append(((i, 3), (i, 1)))
    cox_edges.append(((i, 3), (i, 2)))
idx = {}
for e in cox_edges:
    for v in e:
        if v not in idx:
            idx[v] = len(idx)
GRAPHS["Coxeter       (n=28)"] = edges_to_adj(28, [(idx[a], idx[b]) for a, b in cox_edges])

def main():
    log("=" * 100)
    log("DYADIC CYCLE CHECKER -- verification on known graphs (kind-pasteur-S2, Branch III)")
    log("=" * 100)
    log("")
    log(f"{'graph':24s} {'n':>3s} {'girth':>5s} {'c4':>4s} {'c5':>4s} {'c6':>4s} {'c7':>5s} "
        f"{'c8':>5s} {'c16':>7s} {'pow2-status'}")
    for name, adj in GRAPHS.items():
        n = len(adj)
        assert all(len(a) == 3 for a in adj), f"{name} not cubic!"
        g = girth(adj)
        c4 = count_cycles_len(adj, 4)
        c5 = count_cycles_len(adj, 5)
        c6 = count_cycles_len(adj, 6)
        c7 = count_cycles_len(adj, 7)
        c8 = count_cycles_len(adj, 8)
        c16 = count_cycles_len(adj, 16)
        c32 = count_cycles_len(adj, 32) if n >= 32 else 0
        pres = [L for L, c in [(4, c4), (8, c8), (16, c16), (32, c32)] if c > 0]
        status = f"has C{pres[0]}" if pres else "NO C4/C8/C16(/C32)"
        log(f"{name:24s} {n:3d} {g:5d} {c4:4d} {c5:4d} {c6:4d} {c7:5d} {c8:5d} {c16:7d} {status}")
    log("")

    # full spectra + total cycle counts for the psi-extremal verification graphs
    log("-" * 100)
    log("FULL CYCLE SPECTRA (verification against Markstrom 2004 Table 1 psi values)")
    log("-" * 100)
    targets = {
        "Petersen      (n=10)": 57,     # psi(6)
        "Heawood      (n=14)".replace("Heawood", "Heawood "): 213,  # psi(8)
        "McGee         (n=24)": 5608,   # psi(13) attained by McGee
        "TutteCoxeter  (n=30)": 41400,  # psi(16) attained by Tutte-Coxeter
    }
    targets = {"Petersen      (n=10)": 57, "Heawood       (n=14)": 213,
               "McGee         (n=24)": 5608, "TutteCoxeter  (n=30)": 41400}
    for name, psi in targets.items():
        adj = GRAPHS[name]
        n = len(adj)
        t0 = time.time()
        spec = cycle_spectrum(adj, n)
        tot = sum(spec.values())
        ok = "MATCH" if tot == psi else f"MISMATCH (expected {psi})"
        log(f"{name}: total cycles = {tot}  [psi target {psi}: {ok}]  ({time.time()-t0:.1f}s)")
        log("   spectrum: " + ", ".join(f"c{L}={c}" for L, c in spec.items() if c > 0))
    log("")
    log("Petersen literature check: c5=12 c6=10 c8=15 c9=20 expected.")
    log("McGee c8 == 0  <=>  McGee is one of Markstrom's 4 cubic n=24 graphs w/o C4,C8.")

    with open("05-knowledge/results/dyadic_cycle_checker_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")
    log("\nsaved -> 05-knowledge/results/dyadic_cycle_checker_kps2.out")

if __name__ == "__main__":
    main()
