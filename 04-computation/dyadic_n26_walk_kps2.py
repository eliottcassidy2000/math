"""
kind-pasteur-2026-06-09-S2 : BRANCH III: complete the n=26 C4+C8-free class census
(23 classes exist, Markstrom Table 3) by RANDOM WALK inside the {cubic, C4-free, c8=0}
stratum: from each known class, do double-edge-switches that preserve C4-freeness and
C8-freeness; collect every iso class encountered. Reads seed classes from the collector
output. Reports the girth of every class found.

Output -> 05-knowledge/results/dyadic_n26_walk_kps2.out
"""
import sys, os, time, random, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from dyadic_cycle_checker_kps2 import count_cycles_len, girth
from dyadic_gap_hunt_kps2 import (count_cycles_capped, adj_to_edgeset, edgeset_to_adj,
                                  is_connected, fmt_adj)
from dyadic_markstrom_rediscovery_kps2 import propose_switch_c4free, to_nx

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def save():
    with open("05-knowledge/results/dyadic_n26_walk_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

def parse_classes(path, n_want):
    """Pull 'adj:' lines that follow 'n=<n_want> NEW class' lines."""
    txt = open(path, encoding="utf-8").read().splitlines()
    adjs = []
    take = False
    for line in txt:
        if f"n={n_want} NEW class" in line:
            take = True
            continue
        if take and "adj:" in line:
            s = line.split("adj:", 1)[1].strip()
            adj = {}
            for part in s.split(";"):
                part = part.strip()
                if not part:
                    continue
                v, lst = part.split(":")
                adj[int(v)] = eval(lst)
            adjs.append([sorted(adj[i]) for i in range(len(adj))])
            take = False
    return adjs

def main():
    rng = random.Random(2626)
    src = os.path.join("05-knowledge", "results", "dyadic_class_collector_kps2.out")
    seeds = parse_classes(src, 26)
    log(f"loaded {len(seeds)} seed classes at n=26 from collector output")
    classes = []
    for a in seeds:
        G = to_nx(a)
        if not any(nx.is_isomorphic(G, H) for H, _ in classes):
            classes.append((G, a))
    log(f"distinct seeds: {len(classes)} / 23 known to exist")
    t0 = time.time()
    TIME_CAP = 900.0
    walk_steps = 0
    while time.time() - t0 < TIME_CAP and len(classes) < 23:
        # pick a random known class, walk from it
        _, a = classes[rng.randrange(len(classes))]
        adj = [list(x) for x in a]
        es = adj_to_edgeset(adj)
        for step in range(400):
            adj2, es2, ok = propose_switch_c4free(adj, es, rng)
            if not ok:
                continue
            if count_cycles_len(adj2, 8) != 0:
                continue
            adj, es = adj2, es2
            walk_steps += 1
            G = to_nx(adj)
            if not any(nx.is_isomorphic(G, H) for H, _ in classes):
                classes.append((G, [list(x) for x in adj]))
                g = girth(adj)
                c3 = count_cycles_len(adj, 3)
                c16 = count_cycles_capped(adj, 16)
                log(f"  WALK class #{len(classes)}: girth={g} c3={c3} c16={c16} "
                    f"conn={is_connected(adj)}  ({time.time()-t0:.0f}s)")
                log(f"    adj: {fmt_adj(adj)}")
                save()
    log(f"\nFINAL: {len(classes)}/23 classes at n=26 after {walk_steps} successful walk "
        f"steps, {time.time()-t0:.0f}s")
    log(f"girth histogram: {sorted(girth(a) for _, a in classes)}")
    g5 = sum(1 for _, a in classes if girth(a) >= 5)
    log(f"classes with girth>=5: {g5}")
    save()

if __name__ == "__main__":
    main()
