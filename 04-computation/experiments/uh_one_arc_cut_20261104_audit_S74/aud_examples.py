# Small explicit checks: the note's order-4 route-(a) example, and an order-5 witness for route (c)
# (a uniquely Hamiltonian tournament with a non-critical vertex of in- and out-degree >= 2).
import itertools
def hcycles(n, arcs):
    cnt = 0
    for p in itertools.permutations(range(1, n)):
        c = (0,) + p
        if all((c[i], c[(i + 1) % n]) in arcs for i in range(n)): cnt += 1
    return cnt
def strong(V, arcs):
    V = list(V)
    def reach(s, fwd):
        seen = {s}; st = [s]
        while st:
            x = st.pop()
            for y in V:
                if y not in seen and (((x, y) if fwd else (y, x)) in arcs): seen.add(y); st.append(y)
        return seen
    return len(reach(V[0], True)) == len(V) and len(reach(V[0], False)) == len(V)
# (a) the note's example: 0->3, 3->1, 3->2, 1->0, 2->0, pair {1,2} deleted
D = {(0, 3), (3, 1), (3, 2), (1, 0), (2, 0)}
print('route (a) example: strong', strong(range(4), D), '| Hamiltonian cycles', hcycles(4, D),
      '| with 1->2 restored: hc', hcycles(4, D | {(1, 2)}), '| with 2->1 restored: hc', hcycles(4, D | {(2, 1)}))
# (c) S4 = cycle 0->1->2->3->0 with chords 0->2, 1->3; add v=4 with out-neighbourhood O, |O| = 2 (cyclic intervals)
S4 = {(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)}
for O in [(0, 1), (1, 2), (2, 3), (3, 0)]:
    arcs = set(S4) | {(4, o) for o in O} | {(i, 4) for i in range(4) if i not in O}
    hc = hcycles(5, arcs)
    nc = strong(range(4), {a for a in arcs if 4 not in a})
    print('S4 + v with N+(v) =', set(O), '| hc =', hc, '| T - v strong:', nc, '| out-degree of v = 2, in-degree = 2')
