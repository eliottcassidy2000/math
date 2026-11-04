# The four order-9 route-(b) witnesses (out-rows from aud_routes classes 9): vertex connectivity, arc connectivity,
# and the bad arc. Fraisse-Thomassen (1987): every 2-strong tournament minus any one arc is Hamiltonian, so each
# witness must have a cut vertex.
import itertools
W = [([34, 36, 33, 7, 15, 472, 159, 287, 95], (4, 3)), ([18, 20, 17, 7, 456, 31, 175, 303, 111], (5, 3)),
     ([10, 12, 9, 448, 15, 31, 183, 311, 119], (5, 4)), ([10, 12, 9, 464, 39, 15, 183, 311, 119], (4, 5))]
def strong(V, arcs):
    V = list(V)
    if len(V) <= 1: return True
    def reach(s, fwd):
        seen = {s}; st = [s]
        while st:
            x = st.pop()
            for y in V:
                if y not in seen and ((x, y) in arcs if fwd else (y, x) in arcs):
                    seen.add(y); st.append(y)
        return seen
    return len(reach(V[0], True)) == len(V) and len(reach(V[0], False)) == len(V)
def ham_cycle(V, arcs):
    V = list(V); s = V[0]
    for p in itertools.permutations(V[1:]):
        c = (s,) + p
        if all((c[i], c[(i + 1) % len(c)]) in arcs for i in range(len(c))): return True
    return False
for rows, bad in W:
    n = len(rows)
    arcs = {(u, v) for u in range(n) for v in range(n) if (rows[u] >> v) & 1}
    cut_vertices = [v for v in range(n) if not strong([u for u in range(n) if u != v], {a for a in arcs if v not in a})]
    arc_conn2 = all(strong(range(n), arcs - {a}) for a in arcs)
    nonham = [a for a in arcs if not ham_cycle(range(n), arcs - {a})] if True else None
    print('rows', rows, '| cut vertices', cut_vertices, '| 2-arc-strong', arc_conn2, '| arcs e with T - e non-Hamiltonian', nonham, '| reported', bad)
