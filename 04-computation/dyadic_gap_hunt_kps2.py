"""
kind-pasteur-2026-06-09-S2 : BRANCH III (dyadic-gap hunt, HYP-2359 / Erdos-Gyarfas / Erdos 64)
Part 2: THE HUNT. Simulated annealing over girth>=5 cubic graphs at n in [20,40],
minimizing #C8; then, from C8-free specimens, hill-climb minimizing #C16 inside the
{girth>=5, C8-free} subspace.

Context (Markstrom 2004, Congr. Numer. 171): all cubic graphs on <29 vertices contain a
cycle of length in {4,8,16}; the smallest cubic graphs with no C4 and no C8 have 24
vertices (exactly 4 of them; 23 at n=26, 251 at n=28). A cycle of length 32 needs >=32
vertices, so at n<=31 a cubic graph avoiding {4,8,16} would be a full counterexample to
Erdos-Gyarfas. We hunt in the HARDER girth>=5 regime (no C3 either).

Warm starts at n=24 (McGee), n=28 (Coxeter), n=30 (Tutte-Coxeter) -- cage modifications.

Output -> 05-knowledge/results/dyadic_gap_hunt_kps2.out
"""
import sys, os, time, random, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dyadic_cycle_checker_kps2 import (count_cycles_len, girth, cycle_spectrum,
                                       edges_to_adj, GRAPHS, bfs_dist_restricted)
from collections import deque

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def save():
    with open("05-knowledge/results/dyadic_gap_hunt_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

rng = random.Random(20260609)

# ---------------------------------------------------------- capped C16 counter
def count_cycles_capped(adj, L, cap=None, budget=None):
    """Like count_cycles_len but stops once count > cap (returns cap+1) and/or
    budget expansions exceeded (returns -1)."""
    n = len(adj)
    if L < 3 or L > n:
        return 0
    count = 0
    expansions = 0
    for s in range(n):
        dist = bfs_dist_restricted(adj, s, n)
        on_path = [False] * n
        on_path[s] = True
        nbrs_s = [w for w in adj[s] if w > s and dist[w] <= L - 1]
        for v1 in nbrs_s:
            on_path[v1] = True
            stack = [(v1, 1, iter(adj[v1]))]
            while stack:
                u, depth, it = stack[-1]
                if budget is not None:
                    expansions += 1
                    if expansions > budget:
                        return -1
                advanced = False
                for w in it:
                    if depth == L - 1:
                        if w == s and u > v1:
                            count += 1
                            if cap is not None and count > cap:
                                return count
                        continue
                    if w <= s or on_path[w]:
                        continue
                    if dist[w] > L - depth - 1:
                        continue
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

# ---------------------------------------------------------- graph plumbing
def adj_to_edgeset(adj):
    return set(frozenset((u, w)) for u in range(len(adj)) for w in adj[u] if u < w)

def edgeset_to_adj(n, es):
    adj = [[] for _ in range(n)]
    for e in es:
        a, b = tuple(e)
        adj[a].append(b)
        adj[b].append(a)
    return adj

def is_connected(adj):
    n = len(adj)
    seen = [False] * n
    seen[0] = True
    q = deque([0])
    c = 1
    while q:
        u = q.popleft()
        for w in adj[u]:
            if not seen[w]:
                seen[w] = True
                c += 1
                q.append(w)
    return c == n

def new_edge_short_cycle(adj, x, y, maxbad=4):
    """True if edge (x,y) lies on a cycle of length 3 or 4 in adj (girth violation)."""
    # C3: common neighbor
    Nx = set(adj[x])
    for w in adj[y]:
        if w != x and w in Nx:
            return True
    # C4 through (x,y): u in N(x)\{y}, v in N(y)\{x}, u~v, u!=v
    for u in adj[x]:
        if u == y:
            continue
        Nu = set(adj[u])
        for v in adj[y]:
            if v != x and v != u and v in Nu:
                return True
    return False

def random_cubic_girth5(n, rng, max_tries=200):
    """Random cubic graph via pairing model, then switch away C3/C4."""
    import networkx as nx
    for _ in range(max_tries):
        try:
            G = nx.random_regular_graph(3, n, seed=rng.randint(0, 10**9))
        except Exception:
            continue
        if not nx.is_connected(G):
            continue
        adj = [sorted(G[v]) for v in range(n)]
        # anneal away short cycles: objective = #C3 + #C4
        es = adj_to_edgeset(adj)
        cur = count_cycles_len(adj, 3) + count_cycles_len(adj, 4)
        tries = 0
        while cur > 0 and tries < 20000:
            tries += 1
            adj2, es2, ok = propose_switch(adj, es, rng, enforce_girth5=False)
            if not ok:
                continue
            v = count_cycles_len(adj2, 3) + count_cycles_len(adj2, 4)
            if v <= cur:
                adj, es, cur = adj2, es2, v
        if cur == 0 and is_connected(adj):
            return adj
    return None

def propose_switch(adj, es, rng, enforce_girth5=True):
    """Random double edge switch. Returns (new_adj, new_es, ok)."""
    n = len(adj)
    eslist = list(es)
    e1 = eslist[rng.randrange(len(eslist))]
    e2 = eslist[rng.randrange(len(eslist))]
    if e1 == e2:
        return adj, es, False
    a, b = tuple(e1)
    c, d = tuple(e2)
    if len({a, b, c, d}) < 4:
        return adj, es, False
    if rng.random() < 0.5:
        new1, new2 = (a, c), (b, d)
    else:
        new1, new2 = (a, d), (b, c)
    f1, f2 = frozenset(new1), frozenset(new2)
    if f1 in es or f2 in es or f1 == f2:
        return adj, es, False
    es2 = set(es)
    es2.discard(e1)
    es2.discard(e2)
    es2.add(f1)
    es2.add(f2)
    adj2 = edgeset_to_adj(n, es2)
    if enforce_girth5:
        if new_edge_short_cycle(adj2, *new1) or new_edge_short_cycle(adj2, *new2):
            return adj, es, False
    return adj2, es2, True

# ---------------------------------------------------------- stage 1: minimize #C8
def anneal_c8(adj0, rng, iters=9000, T0=4.0, T1=0.05, label=""):
    """SA minimizing #C8 over girth>=5 cubic graphs. Returns (best_c8, best_adj)."""
    adj = [list(a) for a in adj0]
    es = adj_to_edgeset(adj)
    cur = count_cycles_len(adj, 8)
    best, best_adj = cur, [list(a) for a in adj]
    for i in range(iters):
        if best == 0:
            break
        T = T0 * (T1 / T0) ** (i / iters)
        adj2, es2, ok = propose_switch(adj, es, rng, enforce_girth5=True)
        if not ok:
            continue
        v = count_cycles_len(adj2, 8)
        if v <= cur or rng.random() < math.exp((cur - v) / max(T, 1e-9)):
            adj, es, cur = adj2, es2, v
            if cur < best:
                best, best_adj = cur, [list(a) for a in adj]
    return best, best_adj

# ---------------------------------------------------------- stage 2: minimize #C16 in {c8=0}
def descend_c16(adj0, rng, max_evals=250, time_cap=300.0, label=""):
    """Hill-climb (with plateau walk) on #C16 within {girth>=5, c8=0}."""
    adj = [list(a) for a in adj0]
    es = adj_to_edgeset(adj)
    cur = count_cycles_capped(adj, 16)
    best, best_adj = cur, [list(a) for a in adj]
    t0 = time.time()
    evals = 0
    fails = 0
    while evals < max_evals and time.time() - t0 < time_cap and best > 0 and fails < 4000:
        adj2, es2, ok = propose_switch(adj, es, rng, enforce_girth5=True)
        if not ok:
            fails += 1
            continue
        if count_cycles_len(adj2, 8) != 0:
            fails += 1
            continue
        v = count_cycles_capped(adj2, 16, cap=cur + 5)
        evals += 1
        if v <= cur:
            adj, es, cur = adj2, es2, v
            fails = 0
            if cur < best:
                best, best_adj = cur, [list(a) for a in adj]
    return best, best_adj, evals

def dyadic_profile(adj):
    n = len(adj)
    p = {}
    for L in (3, 4, 5, 6, 7, 8):
        p[L] = count_cycles_len(adj, L)
    p[16] = count_cycles_len(adj, 16)
    if n >= 32:
        p[32] = count_cycles_capped(adj, 32, budget=80_000_000)
    return p

def fmt_adj(adj):
    return "; ".join(f"{u}:{sorted(adj[u])}" for u in range(len(adj)))

# ---------------------------------------------------------- main
def main():
    t_start = time.time()
    log("=" * 100)
    log("DYADIC-GAP HUNT: girth>=5 cubic graphs minimizing #C8 then #C16  (kind-pasteur-S2)")
    log("=" * 100)
    NS = list(range(20, 41, 2))
    RESTARTS = 8
    ITERS = 60000
    frontier = {}
    specimens = {}   # n -> list of (c8=0) adjacency lists
    warm = {24: GRAPHS["McGee         (n=24)"],
            28: GRAPHS["Coxeter       (n=28)"],
            30: GRAPHS["TutteCoxeter  (n=30)"]}

    for n in NS:
        log(f"\n--- n={n} " + "-" * 80)
        best_n = None
        best_adj_n = None
        found = []
        starts = []
        for r in range(RESTARTS):
            g0 = random_cubic_girth5(n, rng)
            if g0 is not None:
                starts.append((f"rand{r}", g0))
        if n in warm:
            starts.append(("warm-cage", [list(a) for a in warm[n]]))
        for tag, g0 in starts:
            t0 = time.time()
            c8_0 = count_cycles_len(g0, 8)
            b, badj = anneal_c8(g0, rng, iters=ITERS, label=f"n={n}/{tag}")
            dt = time.time() - t0
            log(f"  [{tag:9s}] start c8={c8_0:4d} -> best c8={b:4d}   ({dt:.0f}s)")
            if best_n is None or b < best_n:
                best_n, best_adj_n = b, badj
            if b == 0:
                found.append(badj)
        frontier[n] = best_n
        # independent verification of the best specimen
        g = girth(best_adj_n)
        conn = is_connected(best_adj_n)
        c8v = count_cycles_len(best_adj_n, 8)
        log(f"  BEST at n={n}: c8={best_n} (verify recount={c8v}), girth={g}, connected={conn}")
        if found:
            seen_keys = set()
            uniq = []
            for a in found:
                k = tuple(sorted(tuple(sorted(x)) for x in a))
                if k not in seen_keys:
                    seen_keys.add(k)
                    uniq.append(a)
            specimens[n] = uniq
            for i, a in enumerate(uniq):
                p = dyadic_profile(a)
                log(f"  C8-FREE specimen {i}: girth={girth(a)} connected={is_connected(a)} "
                    f"profile c3={p[3]} c4={p[4]} c5={p[5]} c6={p[6]} c7={p[7]} c8={p[8]} "
                    f"c16={p[16]}" + (f" c32={p.get(32)}" if 32 in p else ""))
                log(f"    adj: {fmt_adj(a)}")
        save()

    log("\n" + "=" * 100)
    log("FRONTIER TABLE (stage 1): n -> min #C8 found over girth>=5 cubic graphs")
    log("=" * 100)
    for n in NS:
        log(f"  n={n:2d}: min #C8 = {frontier[n]}")
    save()

    # ---------------- stage 2: kill C16 inside {c8=0}
    log("\n" + "=" * 100)
    log("STAGE 2: minimize #C16 within {girth>=5, C8-free}  (the {8,16}-avoider hunt)")
    log("=" * 100)
    frontier16 = {}
    for n in sorted(specimens):
        best16 = None
        best16_adj = None
        for i, a in enumerate(specimens[n][:3]):
            c16_0 = count_cycles_capped(a, 16)
            b, badj, ev = descend_c16(a, rng, max_evals=8000, time_cap=300.0)
            log(f"  n={n} spec{i}: c16 {c16_0} -> {b}  ({ev} c16-evals)")
            if best16 is None or b < best16:
                best16, best16_adj = b, badj
            save()
        frontier16[n] = best16
        if best16 is not None and best16_adj is not None:
            p = dyadic_profile(best16_adj)
            ok8 = count_cycles_len(best16_adj, 8)
            log(f"  n={n} BEST {{8,16}}-profile: c8={ok8} c16={best16} girth={girth(best16_adj)} "
                f"connected={is_connected(best16_adj)}")
            if best16 == 0:
                log("  *** {8,16}-AVOIDER FOUND *** full profile: " + str(p))
                log("    adj: " + fmt_adj(best16_adj))
                if n < 32:
                    log("  !!! n<32: NO C32 POSSIBLE -> THIS WOULD BE AN ERDOS-GYARFAS "
                        "COUNTEREXAMPLE. RE-VERIFY EVERYTHING. !!!")
    log("\nFRONTIER (stage 2): n -> min #C16 with c8=0, girth>=5:")
    for n in sorted(frontier16):
        log(f"  n={n:2d}: min #C16 = {frontier16[n]}")
    log(f"\ntotal time {time.time()-t_start:.0f}s")
    save()

if __name__ == "__main__":
    main()
