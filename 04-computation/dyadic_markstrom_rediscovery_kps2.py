"""
kind-pasteur-2026-06-09-S2 : BRANCH III (dyadic-gap hunt, HYP-2359 / Erdos-Gyarfas)
Part 4: C4-FREE hunt (triangles ALLOWED -- the true counterexample regime).

Markstrom 2004: the smallest cubic graphs with no C4 and no C8 have n=24 (exactly 4 iso
classes; one planar, built from K4 by repeated triangle-expansion); 23 classes at n=26,
251 at n=28. He verified ALL cubic graphs on <29 vertices contain a 2^k cycle.

This script:
 (1) REDISCOVERY/validation: anneal in {cubic, C4-free} at n=22 (expect min c8 > 0) and
     n=24 (expect c8=0 reachable, up to 4 iso classes). Verify found specimens' dyadic
     profile (only power-of-2 cycles should be C16 per Markstrom).
 (2) FRONTIER PUSH at n=30,32 (beyond Markstrom's exhaustive range): minimize c8 in the
     C4-free regime, then minimize c16 within {c4=0, c8=0}. At n=30: a {4,8,16}-avoider
     would be a GENUINE counterexample to Erdos-Gyarfas (C32 impossible at n<32).

Output -> 05-knowledge/results/dyadic_markstrom_rediscovery_kps2.out
"""
import sys, os, time, random, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dyadic_cycle_checker_kps2 import count_cycles_len, girth, cycle_spectrum, edges_to_adj
from dyadic_gap_hunt_kps2 import (count_cycles_capped, adj_to_edgeset, edgeset_to_adj,
                                  is_connected, fmt_adj, dyadic_profile)
from collections import deque
import networkx as nx

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def save():
    with open("05-knowledge/results/dyadic_markstrom_rediscovery_kps2.out", "w",
              encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

rng = random.Random(64)

def new_edge_c4(adj, x, y):
    """True if edge (x,y) lies on a 4-cycle."""
    for u in adj[x]:
        if u == y:
            continue
        Nu = set(adj[u])
        for v in adj[y]:
            if v != x and v != u and v in Nu:
                return True
    return False

def propose_switch_c4free(adj, es, rng):
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
    if new_edge_c4(adj2, *new1) or new_edge_c4(adj2, *new2):
        return adj, es, False
    return adj2, es2, True

def random_cubic_c4free(n, rng, max_tries=100):
    for _ in range(max_tries):
        try:
            G = nx.random_regular_graph(3, n, seed=rng.randint(0, 10**9))
        except Exception:
            continue
        if not nx.is_connected(G):
            continue
        adj = [sorted(G[v]) for v in range(n)]
        es = adj_to_edgeset(adj)
        cur = count_cycles_len(adj, 4)
        tries = 0
        while cur > 0 and tries < 30000:
            tries += 1
            # use unconstrained switch, hill-climb on c4
            eslist = list(es)
            e1 = eslist[rng.randrange(len(eslist))]
            e2 = eslist[rng.randrange(len(eslist))]
            if e1 == e2:
                continue
            a, b = tuple(e1); c, d = tuple(e2)
            if len({a, b, c, d}) < 4:
                continue
            if rng.random() < 0.5:
                new1, new2 = (a, c), (b, d)
            else:
                new1, new2 = (a, d), (b, c)
            f1, f2 = frozenset(new1), frozenset(new2)
            if f1 in es or f2 in es or f1 == f2:
                continue
            es2 = set(es); es2.discard(e1); es2.discard(e2); es2.add(f1); es2.add(f2)
            adj2 = edgeset_to_adj(n, es2)
            v = count_cycles_len(adj2, 4)
            if v <= cur:
                adj, es, cur = adj2, es2, v
        if cur == 0 and is_connected(adj):
            return adj
    return None

def anneal_c8_c4free(adj0, rng, iters=60000, T0=4.0, T1=0.05):
    adj = [list(a) for a in adj0]
    es = adj_to_edgeset(adj)
    cur = count_cycles_len(adj, 8)
    best, best_adj = cur, [list(a) for a in adj]
    for i in range(iters):
        if best == 0:
            break
        T = T0 * (T1 / T0) ** (i / iters)
        adj2, es2, ok = propose_switch_c4free(adj, es, rng)
        if not ok:
            continue
        v = count_cycles_len(adj2, 8)
        if v <= cur or rng.random() < math.exp((cur - v) / max(T, 1e-9)):
            adj, es, cur = adj2, es2, v
            if cur < best:
                best, best_adj = cur, [list(a) for a in adj]
    return best, best_adj

def descend_c16_c4free(adj0, rng, max_evals=8000, time_cap=300.0):
    adj = [list(a) for a in adj0]
    es = adj_to_edgeset(adj)
    cur = count_cycles_capped(adj, 16)
    best, best_adj = cur, [list(a) for a in adj]
    t0 = time.time()
    evals = 0
    fails = 0
    while evals < max_evals and time.time() - t0 < time_cap and best > 0 and fails < 5000:
        adj2, es2, ok = propose_switch_c4free(adj, es, rng)
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

def to_nx(adj):
    G = nx.Graph()
    G.add_nodes_from(range(len(adj)))
    for u in range(len(adj)):
        for w in adj[u]:
            if u < w:
                G.add_edge(u, w)
    return G

def main():
    t_start = time.time()
    log("=" * 100)
    log("C4-FREE HUNT (triangles allowed): rediscovery at n=22/24 + frontier push n=30/32")
    log("=" * 100)

    # ---------------- (1) rediscovery
    for n, expect in [(22, "expect min c8 > 0 (Markstrom: smallest C4+C8-free cubic is n=24)"),
                      (24, "expect c8=0 reachable (4 iso classes exist)")]:
        log(f"\n--- n={n}: {expect}")
        found = []
        best_overall = None
        for r in range(10):
            g0 = random_cubic_c4free(n, rng)
            if g0 is None:
                log(f"  [restart {r}] generation FAILED")
                continue
            b, badj = anneal_c8_c4free(g0, rng, iters=60000)
            if best_overall is None or b < best_overall:
                best_overall = b
            if b == 0:
                found.append(badj)
            log(f"  [restart {r}] best c8={b}")
        log(f"  n={n}: min c8 over restarts = {best_overall}")
        if found:
            # iso-dedupe
            classes = []
            for a in found:
                Ga = to_nx(a)
                if not any(nx.is_isomorphic(Ga, Gb) for Gb, _ in classes):
                    classes.append((Ga, a))
            log(f"  n={n}: C8-free specimens found: {len(found)}, distinct iso classes: "
                f"{len(classes)}" + ("  (Markstrom: exactly 4 exist)" if n == 24 else ""))
            for i, (Ga, a) in enumerate(classes):
                p = dyadic_profile(a)
                planar = nx.check_planarity(Ga)[0]
                log(f"   class {i}: girth={girth(a)} planar={planar} connected={is_connected(a)}"
                    f" c3={p[3]} c4={p[4]} c5={p[5]} c6={p[6]} c7={p[7]} c8={p[8]} c16={p[16]}")
                log(f"     adj: {fmt_adj(a)}")
        save()

    # ---------------- (2) frontier push at n=30, 32
    log("\n" + "=" * 100)
    log("FRONTIER PUSH: n=30 and n=32, C4-free regime (beyond Markstrom's exhaustive <29)")
    log("  At n=30: a {4,8,16}-avoider = GENUINE Erdos-Gyarfas counterexample (no C32 at n<32)")
    log("=" * 100)
    for n in (30, 32):
        log(f"\n--- n={n}")
        found = []
        best_overall = None
        for r in range(10):
            g0 = random_cubic_c4free(n, rng)
            if g0 is None:
                log(f"  [restart {r}] generation FAILED")
                continue
            b, badj = anneal_c8_c4free(g0, rng, iters=60000)
            if best_overall is None or b < best_overall:
                best_overall = b
            if b == 0:
                found.append(badj)
            log(f"  [restart {r}] best c8={b}")
        log(f"  n={n}: min c8 = {best_overall}; C8-free specimens: {len(found)}")
        save()
        # stage 2: kill C16
        best16 = None
        best16_adj = None
        for i, a in enumerate(found[:4]):
            c16_0 = count_cycles_capped(a, 16)
            b, badj, ev = descend_c16_c4free(a, rng, max_evals=8000, time_cap=300.0)
            log(f"  n={n} spec{i}: c16 {c16_0} -> {b}  ({ev} c16-evals)")
            if best16 is None or b < best16:
                best16, best16_adj = b, badj
            save()
        if best16 is not None:
            ok4 = count_cycles_len(best16_adj, 4)
            ok8 = count_cycles_len(best16_adj, 8)
            log(f"  n={n} BEST: c4={ok4} c8={ok8} c16={best16} girth={girth(best16_adj)} "
                f"connected={is_connected(best16_adj)}")
            log(f"    adj: {fmt_adj(best16_adj)}")
            if best16 == 0:
                if n < 32:
                    log("  !!! {4,8,16}-AVOIDER AT n<32 -> ERDOS-GYARFAS COUNTEREXAMPLE !!!")
                    log("  full spectrum: " + str(cycle_spectrum(best16_adj)))
                else:
                    c32 = count_cycles_capped(best16_adj, 32, budget=200_000_000)
                    log(f"  {{4,8,16}}-avoider at n=32! c32={c32} "
                        + ("<- COUNTEREXAMPLE!!!" if c32 == 0 else "(killed by C32)"))
    log(f"\ntotal time {time.time()-t_start:.0f}s")
    save()

if __name__ == "__main__":
    main()
