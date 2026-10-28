"""Referee: HYP-9435(3).  (a) referee's n=10 example; (b) the author's proposed n=19 construction
(theta_3^2 with a directed triangle on the midpoints + 5 pairwise non-isomorphic rigid oriented trees)."""
from core import automorphisms, PermGroup
from decide import decide

def report(name, n, arcs):
    el = automorphisms(n, arcs)
    G = PermGroup(el)
    r = decide(n, arcs, el, want_counts=False)
    print(f"{name}: n={n} |A|={G.N} involutions={sum(1 for o in G.order_of if o == 2)} "
          f"defects(|K|,delta)={r['defects']} linear={r['linear']} weak={r['weak']}")

report("K4(C3+source)+K1+K2+dirP3", 10,
       [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3), (3, 1), (5, 6), (7, 8), (8, 9)])
arcs = [(0, 2), (0, 3), (0, 4), (2, 1), (3, 1), (4, 1), (2, 3), (3, 4), (4, 2)]   # theta + triangle, 5 vertices
arcs += [(6, 7)]                                   # K2 (vertex 5 = K1)
arcs += [(8, 9), (9, 10)]                          # directed P3
arcs += [(11, 12), (12, 13), (13, 14)]             # directed P4
arcs += [(15, 16), (16, 17), (18, 17)]             # P4 oriented ->,->,<- (rigid)
report("author's HYP-9435(3) construction", 19, arcs)
