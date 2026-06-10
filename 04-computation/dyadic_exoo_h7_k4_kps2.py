"""
kind-pasteur-2026-06-09-S2 : BRANCH III (dyadic-gap hunt, HYP-2359 / Erdos-Gyarfas)
LITERATURE CROSS-VALIDATION of Exoo (arXiv:1403.5636, 'Three Graphs and the
Erdos-Gyarfas Conjecture'):
  * Exoo's H7 gadget (7 vertices, 9 edges, attachment vertices u,v,w of degree 2;
    no 2^m-cycle for any m; d(u,v)=d(u,w)=2, d(v,w)=3).
  * Exoo states one of Markstrom's four minimal n=24 cubic {C4,C8}-free graphs
    is obtained from K4 by replacing three vertices with H7 and the fourth by K3.
  * We enumerate ALL 6^3 = 216 attachment wirings, find which are {C4,C8}-free,
    and test isomorphism against the planar n=24 class our SA hunt rediscovered
    (dyadic_markstrom_rediscovery_kps2.out class 1: c3=7, planar, |Aut|=3).
  * Also: f(2)=10 cross-check (Exoo: exactly three cubic order-10 graphs with no C4,
    incl. Petersen) against our mindeg3 census (which found exactly 3).
Output -> 05-knowledge/results/dyadic_exoo_h7_k4_kps2.out
"""
import sys, os, time
from itertools import permutations
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from dyadic_cycle_checker_kps2 import count_cycles_len, girth, edges_to_adj
from dyadic_gap_hunt_kps2 import is_connected
import networkx as nx

OUT = []
def log(s=""):
    print(s, flush=True)
    OUT.append(str(s))
def save():
    with open("05-knowledge/results/dyadic_exoo_h7_k4_kps2.out", "w", encoding="utf-8") as f:
        f.write("\n".join(OUT) + "\n")

# ---- H7 (Exoo Fig. 2): u,a,b,v,c,d,w = 0..6 ---------------------------------------
H7_EDGES = [(0,1),(0,2),(1,3),(1,4),(3,4),(2,5),(2,6),(5,6),(4,5)]
ATTACH = {"u": 0, "v": 3, "w": 6}

def main():
    t0 = time.time()
    log("=" * 100)
    log("EXOO H7 GADGET VERIFICATION (arXiv:1403.5636 Fig. 2)")
    log("=" * 100)
    adj7 = edges_to_adj(7, H7_EDGES)
    G7 = nx.Graph(H7_EDGES)
    deg = sorted(d for _, d in G7.degree())
    log(f"  H7: 7 vertices, {G7.number_of_edges()} edges, degree seq {deg} (need three 2s, four 3s)")
    log(f"  cycles: c3={count_cycles_len(adj7,3)} c4={count_cycles_len(adj7,4)} "
        f"c5={count_cycles_len(adj7,5)} c6={count_cycles_len(adj7,6)} c7={count_cycles_len(adj7,7)}")
    log(f"  no 2^m-cycle for any m: {count_cycles_len(adj7,4)==0}")
    sp = dict(nx.all_pairs_shortest_path_length(G7))
    log(f"  d(u,v)={sp[0][3]} d(u,w)={sp[0][6]} d(v,w)={sp[3][6]}  (Exoo: 2, 2, 3)")
    save()

    log("")
    log("=" * 100)
    log("K4 CONSTRUCTION: three H7 copies + one K3, all 216 attachment wirings")
    log("=" * 100)
    # vertices: copy X in {0,1,2} occupies 7X..7X+6 ; K3 = 21,22,23
    # K4 'macro' edges: (A,B),(A,C),(B,C),(A,D),(B,D),(C,D); D = K3, its vertex 21+i
    # links to macro-neighbor i (i=0 -> A, 1 -> B, 2 -> C).
    # for H7 copy X the macro neighbors in fixed order:
    macro_nbrs = {0: [1, 2, "D"], 1: [0, 2, "D"], 2: [0, 1, "D"]}
    attach_names = ["u", "v", "w"]
    results = {}
    free_classes = []
    for p0 in permutations(range(3)):
        for p1 in permutations(range(3)):
            for p2 in permutations(range(3)):
                perms = [p0, p1, p2]
                # port[X][Y] = vertex of copy X that faces macro neighbor Y
                port = {}
                for X in range(3):
                    for slot, Y in enumerate(macro_nbrs[X]):
                        port[(X, Y)] = 7 * X + ATTACH[attach_names[perms[X][slot]]]
                edges = []
                for X in range(3):
                    for (a, b) in H7_EDGES:
                        edges.append((7 * X + a, 7 * X + b))
                edges += [(21, 22), (22, 23), (21, 23)]
                # macro edges between H7 copies (each pair once)
                for (X, Y) in [(0, 1), (0, 2), (1, 2)]:
                    edges.append((port[(X, Y)], port[(Y, X)]))
                # K3 vertex 21+i to copy i
                for i in range(3):
                    edges.append((21 + i, port[(i, "D")]))
                adj = edges_to_adj(24, edges)
                assert all(len(x) == 3 for x in adj), "not cubic"
                c4 = count_cycles_len(adj, 4)
                c8 = count_cycles_len(adj, 8)
                key = (c4, c8)
                results[key] = results.get(key, 0) + 1
                if c4 == 0 and c8 == 0:
                    G = nx.Graph(edges)
                    if not any(nx.is_isomorphic(G, H) for H in free_classes):
                        free_classes.append(G)
    log(f"  216 wirings -> (c4,c8) histogram: {dict(sorted(results.items()))}")
    log(f"  {{C4,C8}}-free wirings: {results.get((0,0),0)}; distinct iso classes among them: "
        f"{len(free_classes)}")
    save()

    # ---- compare against our SA-rediscovered planar class -------------------------
    planar_adj = {0:[6,17,22], 1:[2,18,23], 2:[1,13,19], 3:[8,15,16], 4:[13,14,17],
                  5:[9,10,11], 6:[0,10,22], 7:[12,15,20], 8:[3,9,15], 9:[5,8,18],
                  10:[5,6,11], 11:[5,10,23], 12:[7,16,21], 13:[2,4,17], 14:[4,19,20],
                  15:[3,7,8], 16:[3,12,21], 17:[0,4,13], 18:[1,9,23], 19:[2,14,20],
                  20:[7,14,19], 21:[12,16,22], 22:[0,6,21], 23:[1,11,18]}
    P = nx.Graph((u, v) for u, vs in planar_adj.items() for v in vs)
    log("")
    log("COMPARISON with SA-rediscovered planar n=24 class (markstrom_rediscovery class 1):")
    for i, G in enumerate(free_classes):
        pl, _ = nx.check_planarity(G)
        aut = sum(1 for _ in nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter())
        iso = nx.is_isomorphic(G, P)
        adjc = edges_to_adj(24, list(G.edges()))
        log(f"  construction class {i}: planar={pl} |Aut|={aut} girth={girth(adjc)} "
            f"c3={count_cycles_len(adjc,3)} c16={count_cycles_len(adjc,16)} "
            f"iso-to-rediscovered-planar={iso}")
    log("")
    log(f"f(3)=24 structural confirmation: Exoo's K4/H7/K3 description "
        f"{'REPRODUCES' if any(nx.is_isomorphic(G, P) for G in free_classes) else 'DOES NOT reproduce'} "
        f"the planar Markstrom graph our SA hunt rediscovered independently.")
    log(f"total time {time.time()-t0:.0f}s")
    save()

if __name__ == "__main__":
    main()
