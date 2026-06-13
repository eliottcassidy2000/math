"""
kind-pasteur-2026-06-09-S2 : BRANCH III (dyadic-gap hunt, HYP-2359 / Erdos-Gyarfas)
FOLLOW-UP to dyadic_cayley_b2b4_kps2: several dihedral 3-reflection Cayley graphs have
c8=0 AND c16=0 (e.g. D_18 refl(0,1,9), n=36). Their c32 entries were inconclusive
because dyadic_row counted ALL 32-cycles under a budget. Here we run EXISTENCE checks
(cap=0 -> early exit at first 32-cycle) with a large budget, on:
  * all even-m specimens with c8=c16=0 (dyadic contact possibly {4} only, since n<64
    means no C64 and cycle length <= n), and
  * the smallest girth-6 (C4-free) C8-free connected specimens D_19 refl(0,1,8) n=38
    and D_22 refl(0,1,5) n=44 (is their dyadic contact exactly {16}?-> no, c16>0; we
    check whether C32 exists too).
Output -> 05-knowledge/results/dyadic_dihedral_c32_kps2.out
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from dyadic_cycle_checker_kps2 import count_cycles_len, girth, edges_to_adj
from dyadic_gap_hunt_kps2 import count_cycles_capped, is_connected

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def save():
    with open("05-knowledge/results/dyadic_dihedral_c32_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

def dihedral_refl(m, j, k):
    edges = set()
    for t in (0, j, k):
        for i in range(m):
            edges.add(frozenset((i, m + (t - i) % m)))
    return edges_to_adj(2 * m, [tuple(e) for e in edges])

CASES = [
    # (m, j, k, note)
    (18, 1, 9,  "n=36 c4=9 c8=0 c16=0"),
    (18, 2, 9,  "n=36 c4=9 c8=0 c16=0"),
    (18, 4, 9,  "n=36 c4=9 c8=0 c16=0"),
    (20, 1, 10, "n=40 c4=10 c8=0 c16=0"),
    (20, 3, 10, "n=40 c4=10 c8=0 c16=0"),
    (22, 1, 11, "n=44 c4=11 c8=0 c16=0"),
    (24, 1, 12, "n=48 c4=12 c8=0 c16=0"),
    (24, 5, 12, "n=48 c4=12 c8=0 c16=0"),
    (26, 1, 13, "n=52 c4=13 c8=0 c16=0"),
    (28, 1, 14, "n=56 c4=14 c8=0 c16=0"),
    (30, 1, 15, "n=60 c4=15 c8=0 c16=0"),
    (19, 1, 8,  "n=38 girth=6 c4=0 c8=0 c16=1995 (smallest C4+C8-free in family)"),
    (22, 1, 5,  "n=44 girth=6 c4=0 c8=0 c16=2464"),
]

def main():
    t0 = time.time()
    log("=" * 100)
    log("C32 EXISTENCE in dihedral 3-reflection Cayley specimens (cap=0 early-exit)")
    log("=" * 100)
    for (m, j, k, note) in CASES:
        adj = dihedral_refl(m, j, k)
        n = 2 * m
        t = time.time()
        ex32 = count_cycles_capped(adj, 32, cap=0, budget=200_000_000)
        dt = time.time() - t
        status = "EXISTS" if ex32 and ex32 > 0 else ("inconclusive(budget)" if ex32 == -1 else "NONE")
        # recompute the basics for the record
        g = girth(adj)
        c4 = count_cycles_len(adj, 4)
        c8 = count_cycles_len(adj, 8)
        c16 = count_cycles_capped(adj, 16)
        full = ""
        if status == "NONE" and c8 == 0 and c16 == 0:
            full = "  *** dyadic spectrum = {4} ONLY (n<64: no longer 2-power possible) ***"
        elif status == "NONE" and c4 == 0 and c8 == 0:
            full = "  *** dyadic spectrum = {16} ONLY ***"
        log(f"  D_{m} refl(0,{j},{k}): n={n} girth={g} conn={is_connected(adj)} "
            f"c4={c4} c8={c8} c16={c16} C32:{status} ({dt:.0f}s)  [{note}]{full}")
        save()
    log(f"total time {time.time()-t0:.0f}s")
    save()

if __name__ == "__main__":
    main()
