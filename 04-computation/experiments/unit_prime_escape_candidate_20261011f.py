"""An actual t60 unit-graph bank beyond the tested small-prime field screens.

Exact Q(sqrt(-3),sqrt(-11),sqrt(-239)) coordinates. The solver searches only
for bounded colouring witnesses; its failure is never an UNSAT certificate.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import prod
import boundary_palette_compression_20261011 as p

RAD = (-3, -11, -239)
Z = (Q(0),)*8
ONE = (Q(1),)+Z[1:]
S = (Q(0), Q(1))+Z[2:]
T = (Q(0), Q(0), Q(1))+Z[3:]
R = Z[:4]+(Q(1),)+Z[5:]
FACTORS = {(i, j): prod(RAD[b] for b in range(3) if (i & j) >> b & 1)
           for i in range(8) for j in range(8)}
CHECKS = 0


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(a, c):
    return tuple(x*c for x in a)


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    ans = [Q(0)]*8
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            if x and y:
                ans[i ^ j] += x*y*FACTORS[i, j]
    return tuple(ans)


def conjugate(a):
    return tuple(x*(-1 if i.bit_count() % 2 else 1) for i, x in enumerate(a))


def norm(a):
    return mul(a, conjugate(a))


def unit_graph(points):
    p.need(type(points) is tuple and len(set(points)) == len(points)
           and all(type(z) is tuple and len(z) == 8
                   and all(type(x) is Q for x in z) for z in points), 'distinct exact field points')
    edges = frozenset((i, j) for i, j in combinations(range(len(points)), 2)
                      if norm(sub(points[i], points[j])) == ONE)
    return p.Graph(tuple(range(len(points))), edges)


def bent_chain():
    a = scale(add(ONE, S), Q(1, 2))
    a2 = sub(a, ONE)
    b = scale(add(scale(ONE, 5), T), Q(1, 6))
    pos, points, cells = Z, [Z], []
    for direction in [ONE]*3+[scale(b, -1)]*3+[conjugate(b)]*3:
        cell = (pos, add(pos, mul(direction, a)),
                add(pos, mul(direction, a2)), add(pos, mul(direction, S)))
        cells.append(cell)
        points.extend(cell[1:])
        pos = cell[3]
    points = tuple(dict.fromkeys(points))
    return points, tuple(cells), pos, a


def colour_witness(g, k, budget=100000):
    """Bounded DSATUR: return a verified witness or an explicitly unresolved result."""
    p.audit_graph(g)
    p.need(type(k) is int and k >= 1 and type(budget) is int and budget >= 1, 'exact palette and budget')
    adj = {v: set() for v in g.vertices}
    for a, b in g.edges:
        adj[a].add(b)
        adj[b].add(a)
    cs, nodes, cutoff = {}, 0, False

    def rec(used):
        nonlocal nodes, cutoff
        if len(cs) == len(g.vertices):
            return tuple(cs[v] for v in g.vertices)
        nodes += 1
        if nodes > budget:
            cutoff = True
            return None
        def key(v):
            return len({cs[u] for u in adj[v] if u in cs}), len(adj[v]), -v
        v = max((v for v in g.vertices if v not in cs), key=key)
        forbidden = {cs[u] for u in adj[v] if u in cs}
        for colour in range(min(k, used+1)):
            if colour not in forbidden:
                cs[v] = colour
                ans = rec(max(used, colour+1))
                if ans is not None:
                    return ans
                del cs[v]
                if cutoff:
                    return None
        return None
    result = rec(0)
    if result is not None:
        p.need(p.valid_colouring(g, k, result), 'independent literal witness check')
    return result, nodes, cutoff


def rejects(fn):
    try:
        fn()
    except ValueError:
        return True
    return False


def main():
    # The eight basis monomials are independent: the three square classes
    # -3,-11,-239 are independent in Q*/Q*2. All coordinate computations are exact.
    for i in range(8):
        ei = tuple(Q(int(i == j)) for j in range(8))
        check(mul(ei, ONE) == ei and conjugate(conjugate(ei)) == ei, 'basis operations')
        for j in range(8):
            ej = tuple(Q(int(j == q)) for q in range(8))
            check(conjugate(mul(ei, ej)) == mul(conjugate(ei), conjugate(ej)), 'conjugation product')
    base, cells, tip, a = bent_chain()
    omega = scale(add(scale(ONE, 119), R), Q(1, 120))
    check(norm(omega) == ONE and norm(tip) == scale(ONE, 60), 'rotation and endpoint norms')
    check(tip == sub(scale(S, 3), mul(S, T)), 'certified endpoint sqrt33+3sqrt-3')
    check(norm(sub(tip, mul(omega, tip))) == ONE, 'actual closing unit edge')
    cell_edges = ((0, 1), (0, 2), (1, 2), (1, 3), (2, 3))
    for cell in cells:
        check(len(set(cell)) == 4, 'injective diamond copy')
        check(all(norm(sub(cell[i], cell[j])) == ONE for i, j in cell_edges), 'all diamond unit edges')

    declared = {tuple(sorted((base.index(cell[i]), base.index(cell[j]))))
                for cell in cells for i, j in cell_edges}
    extra = unit_graph(base).edges-declared
    corners = {base.index(cell[0]) for cell in cells} | {base.index(tip)}
    check(extra == {(6, 12)} and all(set(e) <= corners for e in extra),
          'actual extra unit edge contradicts three-colour corner propagation')
    print('Bent base extra unit edge(6,12) joins forced-equal3-colour chain corners; base is not3-colourable.')

    unrotated = base
    both = tuple(dict.fromkeys(base+tuple(mul(omega, x) for x in base)))
    thick = tuple(dict.fromkeys(add(x, y) for x in both for y in (Z, ONE, a)))
    for name, points in (('bent9-diamond chain', unrotated), ('t60 joined copies', both),
                         ('t60 triangular thickening', thick)):
        g = unit_graph(points)
        best, used = None, None
        for k in (4, 5, 6):
            witness, nodes, cutoff = colour_witness(g, k, 20000)
            print(name, 'palette', k, 'witness', witness is not None, 'nodes', nodes, 'cutoff', cutoff)
            if witness is not None:
                best, used = witness, k
                break
        print('Candidate', name, 'vertices', len(points), 'unit edges', len(g.edges))
        if best is None:
            print('UNRESOLVED colouring search; no uncolourability conclusion.')
        else:
            check(used == 4, 'this declared finite candidate has a four-colour witness')
            check(p.valid_colouring(g, used, best), 'full unit graph witness')
            check(all(best[u] != best[v] for u, v in g.edges), 'every omitted unit edge remains separated')
            print('Colouring receipt:', best)
        if name != 'bent9-diamond chain':
            closing = tuple(sorted((points.index(tip), points.index(mul(omega, tip)))))
            check(closing in g.edges, 'marked exact unit closing edge')
            print('Retained closing edge:', closing)
        print('Coordinates basis RAD=', RAD, '; ordered rational coefficient tuples:')
        print([tuple(str(c) for c in z) for z in points])
    check(rejects(lambda: unit_graph((Z, Z))), 'coincident vertices rejected')
    check(rejects(lambda: unit_graph(((False,)+Z[1:], ONE))), 'boolean coordinate alias rejected')
    check(rejects(lambda: colour_witness(p.Graph((0,), frozenset()), True)), 'boolean palette alias rejected')
    check(rejects(lambda: colour_witness(p.Graph((0,), frozenset()), 4, True)), 'boolean budget alias rejected')
    print('Every joined candidate contains the nine-diamond propagation contradiction at palette3.')
    print('Prime-screen escape is inherited from unit_defect_frontier_20261011f; it is not a chromatic lower bound.')
    print(f'PASS: {CHECKS} exact checks. No six-equality fragment certified.')


if __name__ == '__main__':
    main()
