"""
kind-pasteur-2026-06-09-S2 : BRANCH III: dyadic profiles of high-girth cubic cages.
Girth >= 9 kills C4 and C8 FOR FREE. The (3,9)-cages (n=58, 18 of them), (3,10)-cages
(n=70, 3), (3,11)-cage (n=112, 1), (3,12)-cage (n=126, 1 = Tutte 12-cage / Benson graph).
Data: A. E. Brouwer's cage pages, https://aeb.win.tue.nl/graphs/cages/cages.html
(.dr format: 'n<N>' then one neighbor-list line 'a,b,c;' per vertex 0..N-1).

For each cage: verify cubic, girth, c4=c8=0 (sanity), exact c16; c32 EXISTENCE via
early-exit DFS with expansion budget (existence only; -1 = budget exceeded/inconclusive).
A girth>=9 cubic graph on 32<=n<64 with c16=0=c32 would be an Erdos-Gyarfas
counterexample; n=58 cages are in that window for {16,32}.

Output -> 05-knowledge/results/dyadic_cage_profiles_kps2.out
"""
import sys, os, time, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dyadic_cycle_checker_kps2 import count_cycles_len, girth, edges_to_adj
from dyadic_gap_hunt_kps2 import count_cycles_capped, is_connected

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def save():
    with open("05-knowledge/results/dyadic_cage_profiles_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

def read_dr(path):
    txt = open(path).read().strip().splitlines()
    n = int(txt[0].lstrip("n"))
    adj = []
    for line in txt[1:]:
        line = line.strip().rstrip(".").rstrip(";")
        if not line:
            continue
        adj.append([int(x) for x in line.split(",")])
    assert len(adj) == n, (path, n, len(adj))
    # symmetrize/sanity: file lists each vertex's neighbors
    es = set()
    for u, nb in enumerate(adj):
        for w in nb:
            es.add(frozenset((u, w)))
    return edges_to_adj(n, [tuple(e) for e in es])

def has_cycle_len(adj, L, budget):
    r = count_cycles_capped(adj, L, cap=0, budget=budget)
    return r  # 0 = none, >=1 = exists, -1 = inconclusive

def main():
    base = os.path.join("05-knowledge", "results", "cages_kps2")
    files = sorted(glob.glob(os.path.join(base, "*.dr")))
    log("=" * 100)
    log("HIGH-GIRTH CUBIC CAGES: dyadic cycle profiles (girth>=9 => C8-free for free)")
    log("=" * 100)
    log(f"{'file':22s} {'n':>4s} {'girth':>5s} {'c8':>3s} {'c16':>6s} {'c32exists':>10s} {'time':>6s}")
    for path in files:
        name = os.path.basename(path)
        t0 = time.time()
        adj = read_dr(path)
        n = len(adj)
        assert all(len(a) == 3 for a in adj), name
        g = girth(adj)
        c8 = count_cycles_len(adj, 8)
        c16 = count_cycles_capped(adj, 16)
        c32 = has_cycle_len(adj, 32, budget=30_000_000)
        c32s = {0: "NO", -1: "inconcl."}.get(c32, "YES")
        log(f"{name:22s} {n:4d} {g:5d} {c8:3d} {c16:6d} {c32s:>10s} {time.time()-t0:5.1f}s")
        save()
    log("")
    log("Interpretation: every (3,9)-cage avoids C4+C8 by girth; the c16 column decides")
    log("whether any is a {4,8,16}-avoider on n=58 (then c32 would be the last gate <64).")
    save()

if __name__ == "__main__":
    main()
