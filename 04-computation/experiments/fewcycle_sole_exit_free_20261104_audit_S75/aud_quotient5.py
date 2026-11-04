"""aud_quotient5.py - which order-9 witnesses v => C3 => x => y => C3' => v have an almost transitive 5-vertex quotient
(transitive order with the arc between the two extreme vertices reversed; Achour-Boudabbous-Boussairi), and which have
Burzio-Demaria's bineutral A_5 (x_i -> x_j iff j < i - 1 or j = i + 1; THM-4868's run-length-2 class)?
Quotient vertices: c (= C3), x, y, c2 (= C3'), v."""
from itertools import permutations

names = ['c', 'x', 'y', 'c2', 'v']


def tour(vx, vy):
    arcs = {('c', 'x'), ('c', 'y'), ('c', 'c2'), ('x', 'y'), ('x', 'c2'), ('y', 'c2'), ('v', 'c'), ('c2', 'v')}
    arcs.add(('v', 'x') if vx else ('x', 'v'))
    arcs.add(('v', 'y') if vy else ('y', 'v'))
    return arcs


def almost_transitive(arcs):
    for p in permutations(names):
        if all(((p[i], p[j]) in arcs) != (i == 0 and j == 4) for i in range(5) for j in range(i + 1, 5)):
            return p
    return None


def bineutral(arcs):
    for p in permutations(names):
        if all(((p[i], p[j]) in arcs) == (j < i - 1 or j == i + 1) for i in range(5) for j in range(5) if i != j):
            return p
    return None


for vx in (1, 0):
    for vy in (1, 0):
        a = tour(vx, vy)
        print("v->x" if vx else "x->v", "v->y" if vy else "y->v",
              "| almost transitive order:", almost_transitive(a), "| bineutral A_5 labelling:", bineutral(a))
