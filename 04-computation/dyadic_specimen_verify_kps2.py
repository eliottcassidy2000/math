"""
kind-pasteur-2026-06-09-S2 : BRANCH III: ADVERSARIAL VERIFICATION of hunt specimens
with an INDEPENDENT cycle enumerator (networkx.simple_cycles, length_bound).
Cross-checks the n=28 girth-5 C8-free specimens and the n=24 C4-free (Markstrom-class)
rediscoveries; also isomorphism tests among specimens and against the McGee graph.

Output -> 05-knowledge/results/dyadic_specimen_verify_kps2.out
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import networkx as nx
from collections import Counter
from dyadic_cycle_checker_kps2 import count_cycles_len, girth

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))

def parse_adj(s):
    adj = {}
    for part in s.split(";"):
        part = part.strip()
        if not part:
            continue
        v, lst = part.split(":")
        adj[int(v)] = eval(lst)
    n = len(adj)
    out = [sorted(adj[i]) for i in range(n)]
    return out

SPEC28_0 = ("0:[2, 13, 23]; 1:[6, 7, 24]; 2:[0, 10, 19]; 3:[11, 14, 17]; 4:[13, 17, 27]; "
            "5:[7, 9, 11]; 6:[1, 16, 18]; 7:[1, 5, 8]; 8:[7, 16, 25]; 9:[5, 14, 21]; "
            "10:[2, 15, 17]; 11:[3, 5, 18]; 12:[19, 22, 25]; 13:[0, 4, 24]; 14:[3, 9, 26]; "
            "15:[10, 20, 24]; 16:[6, 8, 20]; 17:[3, 4, 10]; 18:[6, 11, 27]; 19:[2, 12, 20]; "
            "20:[15, 16, 19]; 21:[9, 23, 25]; 22:[12, 23, 26]; 23:[0, 21, 22]; 24:[1, 13, 15]; "
            "25:[8, 12, 21]; 26:[14, 22, 27]; 27:[4, 18, 26]")
SPEC28_1 = ("0:[1, 5, 9]; 1:[0, 12, 21]; 2:[6, 7, 10]; 3:[6, 8, 13]; 4:[11, 21, 24]; "
            "5:[0, 20, 26]; 6:[2, 3, 14]; 7:[2, 8, 23]; 8:[3, 7, 18]; 9:[0, 22, 25]; "
            "10:[2, 18, 24]; 11:[4, 16, 22]; 12:[1, 13, 19]; 13:[3, 12, 24]; 14:[6, 15, 21]; "
            "15:[14, 16, 23]; 16:[11, 15, 17]; 17:[16, 25, 27]; 18:[8, 10, 26]; 19:[12, 20, 22]; "
            "20:[5, 19, 25]; 21:[1, 4, 14]; 22:[9, 11, 19]; 23:[7, 15, 27]; 24:[4, 10, 13]; "
            "25:[9, 17, 20]; 26:[5, 18, 27]; 27:[17, 23, 26]")
SPEC28_2 = ("0:[6, 12, 13]; 1:[2, 8, 14]; 2:[1, 9, 12]; 3:[12, 21, 23]; 4:[5, 11, 26]; "
            "5:[4, 18, 22]; 6:[0, 9, 20]; 7:[10, 16, 26]; 8:[1, 10, 13]; 9:[2, 6, 23]; "
            "10:[7, 8, 27]; 11:[4, 13, 16]; 12:[0, 2, 3]; 13:[0, 8, 11]; 14:[1, 19, 27]; "
            "15:[24, 25, 26]; 16:[7, 11, 25]; 17:[18, 20, 24]; 18:[5, 17, 19]; 19:[14, 18, 23]; "
            "20:[6, 17, 21]; 21:[3, 20, 22]; 22:[5, 21, 25]; 23:[3, 9, 19]; 24:[15, 17, 27]; "
            "25:[15, 16, 22]; 26:[4, 7, 15]; 27:[10, 14, 24]")
# Markstrom-class rediscoveries at n=24 (C4-free, triangles allowed)
M24_0 = ("0:[3, 4, 5]; 1:[4, 9, 11]; 2:[3, 20, 21]; 3:[0, 2, 5]; 4:[0, 1, 17]; 5:[0, 3, 12]; "
         "6:[15, 19, 22]; 7:[8, 14, 20]; 8:[7, 11, 14]; 9:[1, 12, 19]; 10:[13, 16, 17]; "
         "11:[1, 8, 13]; 12:[5, 9, 22]; 13:[10, 11, 17]; 14:[7, 8, 19]; 15:[6, 22, 23]; "
         "16:[10, 18, 23]; 17:[4, 10, 13]; 18:[16, 21, 23]; 19:[6, 9, 14]; 20:[2, 7, 21]; "
         "21:[2, 18, 20]; 22:[6, 12, 15]; 23:[15, 16, 18]")
M24_1 = ("0:[6, 17, 22]; 1:[2, 18, 23]; 2:[1, 13, 19]; 3:[8, 15, 16]; 4:[13, 14, 17]; "
         "5:[9, 10, 11]; 6:[0, 10, 22]; 7:[12, 15, 20]; 8:[3, 9, 15]; 9:[5, 8, 18]; "
         "10:[5, 6, 11]; 11:[5, 10, 23]; 12:[7, 16, 21]; 13:[2, 4, 17]; 14:[4, 19, 20]; "
         "15:[3, 7, 8]; 16:[3, 12, 21]; 17:[0, 4, 13]; 18:[1, 9, 23]; 19:[2, 14, 20]; "
         "20:[7, 14, 19]; 21:[12, 16, 22]; 22:[0, 6, 21]; 23:[1, 11, 18]")
M24_2 = ("0:[1, 16, 21]; 1:[0, 2, 5]; 2:[1, 3, 10]; 3:[2, 17, 19]; 4:[5, 8, 17]; 5:[1, 4, 20]; "
         "6:[9, 11, 23]; 7:[19, 20, 21]; 8:[4, 10, 13]; 9:[6, 13, 23]; 10:[2, 8, 13]; "
         "11:[6, 18, 22]; 12:[15, 16, 20]; 13:[8, 9, 10]; 14:[17, 18, 21]; 15:[12, 16, 23]; "
         "16:[0, 12, 15]; 17:[3, 4, 14]; 18:[11, 14, 22]; 19:[3, 7, 22]; 20:[5, 7, 12]; "
         "21:[0, 7, 14]; 22:[11, 18, 19]; 23:[6, 9, 15]")

def nx_of(adj):
    G = nx.Graph()
    G.add_nodes_from(range(len(adj)))
    for u in range(len(adj)):
        for w in adj[u]:
            if u < w:
                G.add_edge(u, w)
    return G

def main():
    log("=" * 90)
    log("INDEPENDENT VERIFICATION (networkx.simple_cycles) of hunt specimens")
    log("=" * 90)
    specs = [("n28-spec0 (girth5)", SPEC28_0), ("n28-spec1 (girth5)", SPEC28_1),
             ("n28-spec2 (girth5)", SPEC28_2), ("n24-M0 (C4free)", M24_0),
             ("n24-M1 (C4free,planar)", M24_1), ("n24-M2 (C4free)", M24_2)]
    Gs = {}
    for name, s in specs:
        adj = parse_adj(s)
        G = nx_of(adj)
        Gs[name] = G
        cnt = Counter()
        for c in nx.simple_cycles(G, length_bound=9):
            cnt[len(c)] += 1
        ok_cubic = all(d == 3 for _, d in G.degree())
        log(f"{name}: cubic={ok_cubic} connected={nx.is_connected(G)} "
            f"nx-counts c3={cnt[3]} c4={cnt[4]} c5={cnt[5]} c6={cnt[6]} c7={cnt[7]} "
            f"c8={cnt[8]} c9={cnt[9]}")
        log(f"   my-checker  c3={count_cycles_len(adj,3)} c4={count_cycles_len(adj,4)} "
            f"c5={count_cycles_len(adj,5)} c6={count_cycles_len(adj,6)} "
            f"c7={count_cycles_len(adj,7)} c8={count_cycles_len(adj,8)} "
            f"c9={count_cycles_len(adj,9)}   girth={girth(adj)}")
    log("")
    log("isomorphism tests:")
    log(f"  n28-spec1 ~ n28-spec2 ? {nx.is_isomorphic(Gs['n28-spec1 (girth5)'], Gs['n28-spec2 (girth5)'])}")
    log(f"  n28-spec0 ~ n28-spec1 ? {nx.is_isomorphic(Gs['n28-spec0 (girth5)'], Gs['n28-spec1 (girth5)'])}")
    McGee = nx.LCF_graph(24, [12, 7, -7], 8)
    for nm in ("n24-M0 (C4free)", "n24-M1 (C4free,planar)", "n24-M2 (C4free)"):
        log(f"  {nm} ~ McGee ? {nx.is_isomorphic(Gs[nm], McGee)}")
    log(f"  n24-M1 planar? {nx.check_planarity(Gs['n24-M1 (C4free,planar)'])[0]} "
        f"(Markstrom graph is THE planar one: n=24, 36 edges, girth 3)")
    a1 = Gs['n24-M1 (C4free,planar)']
    log(f"  n24-M1: |Aut| = {sum(1 for _ in nx.vf2pp_all_isomorphisms(a1, a1))} "
        f"(MathWorld: Markstrom graph has 3 automorphisms)")
    with open("05-knowledge/results/dyadic_specimen_verify_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

if __name__ == "__main__":
    main()
