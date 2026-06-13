"""
kind-pasteur-2026-06-09-S2 : BRANCH III: the GIRTH LADDER of C8-freeness.
Known after this session: girth 3: min order 24 (Markstrom, exhaustive); girth 5: <=28
(our SA specimens); girth >=9: 58 (cages, trivial). Missing rungs: girth 6 and girth 7.
Hunt: SA minimizing #C8 over cubic graphs of girth >= 6 (and >=7) at n in 28..44.
S710's sharp sub-question (first dyadic window above girth): a girth-6 C8-free cubic
graph dodges window-8; a girth-7 one likewise (McGee does NOT: c8=34).

Output -> 05-knowledge/results/dyadic_girth67_hunt_kps2.out
"""
import sys, os, time, random, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from dyadic_cycle_checker_kps2 import count_cycles_len, girth
from dyadic_gap_hunt_kps2 import (count_cycles_capped, adj_to_edgeset, edgeset_to_adj,
                                  is_connected, fmt_adj)

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def save():
    with open("05-knowledge/results/dyadic_girth67_hunt_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

def short_cycle_through(adj, x, y, gmin):
    """True if edge (x,y) lies on a cycle shorter than gmin (cycle length <= gmin-1).
    BFS from x to y in G - edge(x,y); cycle len = dist+1 < gmin <=> dist <= gmin-2."""
    from collections import deque
    n = len(adj)
    dist = [-1] * n
    dist[x] = 0
    q = deque([x])
    lim = gmin - 2
    while q:
        u = q.popleft()
        if dist[u] >= lim:
            break
        for w in adj[u]:
            if u == x and w == y:
                continue
            if dist[w] == -1:
                dist[w] = dist[u] + 1
                if w == y:
                    return True
                q.append(w)
    return False

def propose_switch_g(adj, es, rng, gmin):
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
    es2.discard(e1); es2.discard(e2); es2.add(f1); es2.add(f2)
    adj2 = edgeset_to_adj(n, es2)
    if short_cycle_through(adj2, *new1, gmin) or short_cycle_through(adj2, *new2, gmin):
        return adj, es, False
    return adj2, es2, True

def random_cubic_girth_g(n, rng, gmin, max_tries=60):
    """Random cubic, then hill-climb away cycles shorter than gmin."""
    for _ in range(max_tries):
        try:
            G = nx.random_regular_graph(3, n, seed=rng.randint(0, 10**9))
        except Exception:
            continue
        if not nx.is_connected(G):
            continue
        adj = [sorted(G[v]) for v in range(n)]
        es = adj_to_edgeset(adj)
        def bad(a):
            return sum(count_cycles_len(a, L) for L in range(3, gmin))
        cur = bad(adj)
        tries = 0
        while cur > 0 and tries < 60000:
            tries += 1
            eslist = list(es)
            e1 = eslist[rng.randrange(len(eslist))]
            e2 = eslist[rng.randrange(len(eslist))]
            if e1 == e2:
                continue
            a, b = tuple(e1)
            c, d = tuple(e2)
            if len({a, b, c, d}) < 4:
                continue
            if rng.random() < 0.5:
                new1, new2 = (a, c), (b, d)
            else:
                new1, new2 = (a, d), (b, c)
            f1, f2 = frozenset(new1), frozenset(new2)
            if f1 in es or f2 in es or f1 == f2:
                continue
            es2 = set(es)
            es2.discard(e1); es2.discard(e2); es2.add(f1); es2.add(f2)
            adj2 = edgeset_to_adj(n, es2)
            v = bad(adj2)
            if v <= cur:
                adj, es, cur = adj2, es2, v
        if cur == 0 and is_connected(adj):
            return adj
    return None

def anneal_c8_g(adj0, rng, gmin, iters=60000, T0=4.0, T1=0.05):
    adj = [list(a) for a in adj0]
    es = adj_to_edgeset(adj)
    cur = count_cycles_len(adj, 8)
    best, best_adj = cur, [list(a) for a in adj]
    for i in range(iters):
        if best == 0:
            break
        T = T0 * (T1 / T0) ** (i / iters)
        adj2, es2, ok = propose_switch_g(adj, es, rng, gmin)
        if not ok:
            continue
        v = count_cycles_len(adj2, 8)
        if v <= cur or rng.random() < math.exp((cur - v) / max(T, 1e-9)):
            adj, es, cur = adj2, es2, v
            if cur < best:
                best, best_adj = cur, [list(a) for a in adj]
    return best, best_adj

def main():
    rng = random.Random(67)
    log("=" * 100)
    log("GIRTH-6 / GIRTH-7 C8-FREE HUNT (the girth ladder of C8-freeness)")
    log("=" * 100)
    for gmin, ns in [(6, [28, 30, 32, 34, 36, 40, 44]), (7, [30, 34, 38, 42, 46])]:
        log(f"\n### girth >= {gmin}:")
        for n in ns:
            best_n = None
            best_adj = None
            R = 6
            t0 = time.time()
            gen_fail = 0
            for r in range(R):
                g0 = random_cubic_girth_g(n, rng, gmin)
                if g0 is None:
                    gen_fail += 1
                    continue
                b, badj = anneal_c8_g(g0, rng, gmin, iters=60000)
                if best_n is None or b < best_n:
                    best_n, best_adj = b, badj
                if best_n == 0:
                    break
            if best_n is None:
                log(f"  n={n}: generation failed in all {R} restarts ({gen_fail} fails) "
                    f"({time.time()-t0:.0f}s)")
                continue
            line = (f"  n={n}: min c8 = {best_n} (girth={girth(best_adj)}, "
                    f"conn={is_connected(best_adj)}, genfail={gen_fail}) "
                    f"({time.time()-t0:.0f}s)")
            log(line)
            if best_n == 0:
                c16 = count_cycles_capped(best_adj, 16)
                c32 = count_cycles_capped(best_adj, 32, budget=60_000_000) if n >= 32 else 0
                log(f"    C8-FREE girth-{gmin} SPECIMEN! c16={c16} c32={c32}")
                log(f"    adj: {fmt_adj(best_adj)}")
            save()
    save()

if __name__ == "__main__":
    main()
