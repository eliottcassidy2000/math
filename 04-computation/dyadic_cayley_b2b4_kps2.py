"""
kind-pasteur-2026-06-09-S2 : BRANCH III (dyadic-gap hunt, HYP-2359 / Erdos-Gyarfas)
COMPLETION RUN: the dyadic_gp_cayley_kps2.py run was interrupted during section (B2);
its .out contains (A) GP graphs and (B1) circulants but NO dihedral / Z2xZm data.
This script completes:
  (B0) independent networkx.simple_cycles recheck of the McGee c8=34 correction
  (B2) dihedral D_m, three reflections {s, s r^j, s r^k}, 2m<=80
  (B3) dihedral D_m, {s, r^a, r^-a}, 2m<=80
  (B4) Z_2 x Z_m, S={(1,0),(0,+-a)}, m even, 2m<=80
Output -> 05-knowledge/results/dyadic_cayley_b2b4_kps2.out
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
    with open("05-knowledge/results/dyadic_cayley_b2b4_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

def dyadic_row(adj):
    n = len(adj)
    g = girth(adj)
    c4 = count_cycles_len(adj, 4)
    c8 = count_cycles_len(adj, 8)
    c16 = c32 = None
    if c8 == 0 and is_connected(adj):
        c16 = count_cycles_capped(adj, 16)
        if n >= 32:
            c32 = count_cycles_capped(adj, 32, budget=5_000_000)
    return g, c4, c8, c16, c32

def main():
    t0 = time.time()

    # ---------- (B0) McGee independent recheck --------------------------------------
    log("=" * 100)
    log("(B0) INDEPENDENT RECHECK: McGee graph cycle counts via networkx.simple_cycles")
    log("=" * 100)
    import networkx as nx
    McGee = nx.LCF_graph(24, [12, 7, -7], 8)
    from collections import Counter
    hist = Counter()
    for cyc in nx.simple_cycles(McGee, length_bound=9):
        hist[len(cyc)] += 1
    log(f"  McGee (LCF [12,7,-7]^8, n=24): nx cycle histogram up to length 9: "
        + ", ".join(f"c{L}={hist[L]}" for L in sorted(hist)))
    log(f"  girth = {min(hist)}; |Aut| = {len(list(nx.algorithms.isomorphism.GraphMatcher(McGee, McGee).isomorphisms_iter()))}")
    adjM = edges_to_adj(24, list(McGee.edges()))
    log(f"  my checker on same graph: c7={count_cycles_len(adjM,7)} c8={count_cycles_len(adjM,8)} "
        f"c9={count_cycles_len(adjM,9)} girth={girth(adjM)}")
    log("  VERDICT: McGee HAS 8-cycles (c8=34) despite girth 7 -- S710's 'McGee -> C16' was an")
    log("  enumeration-order artifact, NOT a C8-free statement.")
    save()

    # ---------- (B2) dihedral, three reflections ------------------------------------
    log("")
    log("=" * 100)
    log("(B2) dihedral D_m, three reflections {s, s r^j, s r^k} (0<j<k<m), 2m<=80")
    log("     (bipartite cubic; canonical class = gap multiset {j, k-j, m-k} up to symmetry)")
    log("=" * 100)
    findings = []
    cnt = 0
    free_list = []
    for m in range(3, 41):
        seen = set()
        for j in range(1, m):
            for k in range(j + 1, m):
                gaps = tuple(sorted((j, k - j, m - k)))
                if gaps in seen:
                    continue
                seen.add(gaps)
                edges = set()
                for t in (0, j, k):
                    for i in range(m):
                        edges.add(frozenset((i, m + (t - i) % m)))
                adj = edges_to_adj(2 * m, [tuple(e) for e in edges])
                if any(len(x) != 3 for x in adj):
                    continue
                cnt += 1
                g, c4, c8, c16, c32 = dyadic_row(adj)
                if c8 == 0:
                    free_list.append((f"D_{m} refl(0,{j},{k})", 2 * m, g, c4, c16, c32,
                                      is_connected(adj)))
        save()
    log(f"  checked {cnt}; C8-free: {len(free_list)}")
    for t in free_list:
        log(f"    {t[0]}: n={t[1]} girth={t[2]} c4={t[3]} c16={t[4]} c32={t[5]} conn={t[6]}")
    findings += free_list
    save()

    # ---------- (B3) dihedral, reflection + rotation --------------------------------
    log("")
    log("(B3) dihedral D_m, {s, r^a, r^-a} (1<=a<m/2), 2m<=80:")
    cnt = 0
    free_list = []
    for m in range(3, 41):
        for a in range(1, (m + 1) // 2):
            edges = set()
            for i in range(m):
                edges.add(frozenset((i, (i + a) % m)))
                edges.add(frozenset((m + i, m + (i + a) % m)))
                edges.add(frozenset((i, m + (-i) % m)))
            adj = edges_to_adj(2 * m, [tuple(e) for e in edges])
            if any(len(x) != 3 for x in adj):
                continue
            cnt += 1
            g, c4, c8, c16, c32 = dyadic_row(adj)
            if c8 == 0:
                free_list.append((f"D_{m} {{s,r^+-{a}}}", 2 * m, g, c4, c16, c32,
                                  is_connected(adj)))
        save()
    log(f"  checked {cnt}; C8-free: {len(free_list)}")
    for t in free_list:
        log(f"    {t[0]}: n={t[1]} girth={t[2]} c4={t[3]} c16={t[4]} c32={t[5]} conn={t[6]}")
    findings += free_list
    save()

    # ---------- (B4) Z_2 x Z_m ------------------------------------------------------
    log("")
    log("(B4) Z_2 x Z_m, S={(1,0),(0,+-a)}, m even, 2m<=80:")
    cnt = 0
    free_list = []
    for m in range(4, 41, 2):
        for a in range(1, (m + 1) // 2):
            edges = set()
            for b in range(2):
                for i in range(m):
                    v = b * m + i
                    edges.add(frozenset((v, b * m + (i + a) % m)))
                    edges.add(frozenset((v, (1 - b) * m + i)))
            adj = edges_to_adj(2 * m, [tuple(e) for e in edges])
            if any(len(x) != 3 for x in adj):
                continue
            cnt += 1
            g, c4, c8, c16, c32 = dyadic_row(adj)
            if c8 == 0:
                free_list.append((f"Z2xZ_{m} a={a}", 2 * m, g, c4, c16, c32, is_connected(adj)))
        save()
    log(f"  checked {cnt}; C8-free: {len(free_list)}")
    for t in free_list:
        log(f"    {t[0]}: n={t[1]} girth={t[2]} c4={t[3]} c16={t[4]} c32={t[5]} conn={t[6]}")
    findings += free_list

    log("")
    log("=" * 100)
    log(f"SUMMARY (B2-B4): C8-free cubic Cayley specimens (incl. disconnected/small): {len(findings)}")
    log(f"  connected, n>=8, C8-free: "
        + str([t[0] for t in findings if t[6] and t[1] >= 8]))
    log(f"total time {time.time()-t0:.0f}s")
    save()

if __name__ == "__main__":
    main()
