"""
kind-pasteur-2026-06-09-S2 : BRANCH III: merge the n=26 C4+C8-free class censuses from
the collector (SA restarts) and the walk (in-stratum random walk); report distinct
classes, their girths, and whether the census of 23 (Markstrom Table 3) is complete.

Output -> 05-knowledge/results/dyadic_n26_census_merge_kps2.out
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from dyadic_cycle_checker_kps2 import count_cycles_len, girth
from dyadic_gap_hunt_kps2 import count_cycles_capped, is_connected, fmt_adj
from dyadic_markstrom_rediscovery_kps2 import to_nx

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))

def parse_adj_lines(path, marker):
    """Pull 'adj:' lines following lines containing the marker."""
    adjs = []
    take = False
    for line in open(path, encoding="utf-8"):
        if marker in line:
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
    base = os.path.join("05-knowledge", "results")
    A = parse_adj_lines(os.path.join(base, "dyadic_class_collector_kps2.out"),
                        "n=26 NEW class")
    B = parse_adj_lines(os.path.join(base, "dyadic_n26_walk_kps2.out"), "WALK class")
    log(f"collector n=26 classes: {len(A)}; walk classes: {len(B)}")
    classes = []
    for src, lst in (("collector", A), ("walk", B)):
        for a in lst:
            assert len(a) == 26 and all(len(x) == 3 for x in a)
            assert count_cycles_len(a, 4) == 0 and count_cycles_len(a, 8) == 0
            G = to_nx(a)
            if not any(nx.is_isomorphic(G, H) for H, _ in classes):
                classes.append((G, a))
    log(f"MERGED distinct classes: {len(classes)} / 23 (Markstrom Table 3)")
    gs = sorted(girth(a) for _, a in classes)
    log(f"girth histogram: {gs}")
    log(f"classes with girth >= 4: {sum(1 for _, a in classes if girth(a) >= 4)}")
    for i, (G, a) in enumerate(classes):
        log(f"  class {i+1:2d}: girth={girth(a)} c3={count_cycles_len(a,3)} "
            f"c16={count_cycles_capped(a,16)} planar={nx.check_planarity(G)[0]} "
            f"|Aut|={sum(1 for _ in nx.vf2pp_all_isomorphisms(G,G))} "
            f"conn={is_connected(a)}")
    with open(os.path.join(base, "dyadic_n26_census_merge_kps2.out"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

if __name__ == "__main__":
    main()
