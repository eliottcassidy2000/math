"""
kind-pasteur-2026-06-09-S2 : BRANCH III: collect ALL iso classes of C4+C8-free cubic
graphs at n=24 (Markstrom: exactly 4 exist) and as many as possible at n=26 (23 exist),
recording the GIRTH of each class.

THE SYLLOGISM: girth>=5 + C8-free  =>  C4-free + C8-free  =>  (by Markstrom's exhaustive
search) the graph is one of the 4 (n=24) / 23 (n=26) known classes. If every collected
class has girth 3, then NO girth>=5 cubic C8-free graph exists at that order -- modulo
classes we fail to find (we report coverage honestly). Combined with our SA hit at n=28
(girth-5 C8-free exists), this brackets the minimum order of a girth>=5 C8-free cubic
graph.

Output -> 05-knowledge/results/dyadic_class_collector_kps2.out
"""
import sys, os, time, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from dyadic_cycle_checker_kps2 import count_cycles_len, girth
from dyadic_gap_hunt_kps2 import count_cycles_capped, is_connected, fmt_adj
from dyadic_markstrom_rediscovery_kps2 import (random_cubic_c4free, anneal_c8_c4free,
                                               to_nx)

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def save():
    with open("05-knowledge/results/dyadic_class_collector_kps2.out", "w",
              encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

def collect(n, target, restarts, rng, iters=60000):
    classes = []  # list of (nxGraph, adj)
    t0 = time.time()
    for r in range(restarts):
        if len(classes) >= target:
            break
        g0 = random_cubic_c4free(n, rng)
        if g0 is None:
            continue
        b, badj = anneal_c8_c4free(g0, rng, iters=iters)
        if b != 0:
            continue
        # sanity: connected, cubic, c4=0, c8=0
        if not is_connected(badj):
            continue
        assert count_cycles_len(badj, 4) == 0 and count_cycles_len(badj, 8) == 0
        G = to_nx(badj)
        if not any(nx.is_isomorphic(G, H) for H, _ in classes):
            classes.append((G, badj))
            g = girth(badj)
            c3 = count_cycles_len(badj, 3)
            c16 = count_cycles_capped(badj, 16)
            planar = nx.check_planarity(G)[0]
            aut = sum(1 for _ in nx.vf2pp_all_isomorphisms(G, G))
            log(f"  n={n} NEW class #{len(classes)}: girth={g} c3={c3} c16={c16} "
                f"planar={planar} |Aut|={aut}  (restart {r}, {time.time()-t0:.0f}s)")
            log(f"    adj: {fmt_adj(badj)}")
            save()
    return classes

def main():
    rng = random.Random(13579)
    log("=" * 100)
    log("CLASS COLLECTOR: C4+C8-free cubic graphs at n=24 (4 exist) and n=26 (23 exist)")
    log("=" * 100)
    log("\n--- n=24, target 4 classes (Markstrom Table 3):")
    cl24 = collect(24, target=4, restarts=120, rng=rng)
    log(f"  n=24: collected {len(cl24)}/4 classes; girths: "
        f"{[girth(a) for _, a in cl24]}")
    save()
    log("\n--- n=26, target 23 classes (Markstrom Table 3):")
    cl26 = collect(26, target=23, restarts=300, rng=rng)
    log(f"  n=26: collected {len(cl26)}/23 classes; girths: "
        f"{sorted(girth(a) for _, a in cl26)}")
    g5 = [a for _, a in cl26 if girth(a) >= 5]
    log(f"  n=26 classes with girth>=5: {len(g5)}")
    log("")
    log("CONCLUSION INPUTS: if all 4 classes at n=24 have girth 3 (and Markstrom's search")
    log("was exhaustive), then no girth>=5 cubic C8-free graph exists on <=24 vertices;")
    log("girth>=5 C8-free found at n=28 by SA (dyadic_gap_hunt_kps2). n=26 status depends")
    log("on the 23-class girth histogram (coverage above).")
    save()

if __name__ == "__main__":
    main()
