"""
kind-pasteur-2026-06-09-S2 : BRANCH III addendum: structured candidates beyond GP/Cayley.
  - Flower snarks J_k (k odd, 4k vertices): J5 (n=20, girth 5), J7 (n=28, girth 6),
    J9 (n=36, girth 6). Dyadic profile c4/c8/c16(/c32).
  - The girth observation: a cubic graph of girth >= 9 contains NO C8 trivially; the
    smallest are the (3,9)-cages at n=58. At 32 <= n < 64 such a graph avoiding also
    C16 and C32 would be an Erdos-Gyarfas counterexample. We check c16 for one
    (3,9)-cage if obtainable; otherwise we report the reasoning only.
  - Deeper SA at n=26 (girth>=5): settle whether min c8 = 0 or 1 (stage-1 found 1).
    Also deeper runs at n=20,22,24 to confirm the floors 6,3,3.

Output -> 05-knowledge/results/dyadic_snark_cage_kps2.out
"""
import sys, os, time, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from dyadic_cycle_checker_kps2 import count_cycles_len, girth, cycle_spectrum, edges_to_adj
from dyadic_gap_hunt_kps2 import (count_cycles_capped, is_connected, fmt_adj,
                                  random_cubic_girth5, anneal_c8, dyadic_profile)

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def save():
    with open("05-knowledge/results/dyadic_snark_cage_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

def flower_snark(k):
    """J_k: k copies of K_{1,3} (center A_i, leaves B_i,C_i,D_i); B-cycle of length k;
    single 2k-cycle (C_1..C_k D_1..D_k). [Wikipedia construction]"""
    A = lambda i: 4 * (i % k)
    B = lambda i: 4 * (i % k) + 1
    C = lambda i: 4 * (i % k) + 2
    D = lambda i: 4 * (i % k) + 3
    edges = set()
    for i in range(k):
        edges.add(frozenset((A(i), B(i))))
        edges.add(frozenset((A(i), C(i))))
        edges.add(frozenset((A(i), D(i))))
        edges.add(frozenset((B(i), B(i + 1))))
    # 2k-cycle C_0..C_{k-1}, D_0..D_{k-1}
    cyc = [C(i) for i in range(k)] + [D(i) for i in range(k)]
    for j in range(2 * k):
        edges.add(frozenset((cyc[j], cyc[(j + 1) % (2 * k)])))
    return edges_to_adj(4 * k, [tuple(e) for e in edges])

def main():
    rng = random.Random(2026)
    log("=" * 100)
    log("FLOWER SNARKS J_k (4k vertices): dyadic profiles")
    log("=" * 100)
    for k in (5, 7, 9):
        adj = flower_snark(k)
        assert all(len(a) == 3 for a in adj)
        n = len(adj)
        g = girth(adj)
        c4 = count_cycles_len(adj, 4)
        c8 = count_cycles_len(adj, 8)
        c16 = count_cycles_capped(adj, 16)
        c32 = count_cycles_capped(adj, 32, budget=150_000_000) if n >= 32 else 0
        log(f"  J{k} (n={n}): girth={g} c4={c4} c8={c8} c16={c16}"
            + (f" c32={c32}" if n >= 32 else "") + f" connected={is_connected(adj)}")
    save()

    log("")
    log("=" * 100)
    log("GIRTH-9 OBSERVATION: cubic girth>=9 => NO C8 for free; smallest such = (3,9)-cages, n=58.")
    log("At 32<=n<64 an {8,16,32}-avoider suffices for a counterexample; girth>=9 handles the 8.")
    log("(3,9)-cage check deferred to a HoG download if available -- see session notes.")
    log("=" * 100)

    log("")
    log("=" * 100)
    log("DEEPER SA at n=20..26 girth>=5: confirm stage-1 floors (6,3,3,1) / try to reach 0 at 26")
    log("=" * 100)
    for n, floor in [(20, 6), (22, 3), (24, 3), (26, 1)]:
        best_n = None
        hits = {}
        t0 = time.time()
        R = 16
        for r in range(R):
            g0 = random_cubic_girth5(n, rng)
            if g0 is None:
                continue
            b, badj = anneal_c8(g0, rng, iters=250000, T0=5.0, T1=0.02)
            hits[b] = hits.get(b, 0) + 1
            if best_n is None or b < best_n:
                best_n, best_adj = b, badj
        log(f"  n={n}: min c8 = {best_n} over {R} deep restarts (250K iters each), "
            f"distribution {dict(sorted(hits.items()))}, prev floor {floor}  "
            f"({time.time()-t0:.0f}s)")
        if best_n == 0:
            p = dyadic_profile(best_adj)
            log(f"    NEW C8-FREE at n={n}! profile: {p}")
            log(f"    adj: {fmt_adj(best_adj)}")
        elif best_n is not None and best_n <= floor:
            # save the extremal witness for the frontier record
            log(f"    extremal witness (c8={best_n}): girth={girth(best_adj)} "
                f"connected={is_connected(best_adj)}")
            log(f"    adj: {fmt_adj(best_adj)}")
        save()

if __name__ == "__main__":
    main()
