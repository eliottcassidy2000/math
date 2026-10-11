"""Exact marked two-port colouring counts and support-safe composition.

New elementary transfer algebra; no six-colour equality gadget is claimed.
Normal and optimized runs use the same always-active exact controls.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import combinations, product
from math import factorial
import boundary_palette_compression_20261011 as p
import boundary_context_weights_20261011d as old

CHECKS = 0


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def natural(n, why):
    p.need(type(n) is int and n >= 0, why)


@dataclass(frozen=True)
class Counts:
    palette: int
    equal: int
    different: int


def audit_counts(c):
    p.need(type(c) is Counts, 'typed count pair')
    p.need(type(c.palette) is int and c.palette >= 2, 'exact palette at least two')
    natural(c.equal, 'exact equal count')
    natural(c.different, 'exact different count')


def compatible(a, b):
    audit_counts(a)
    audit_counts(b)
    p.need(a.palette == b.palette, 'same palette')


def parallel_counts(a, b):
    compatible(a, b)
    return Counts(a.palette, a.equal*b.equal, a.different*b.different)


def series_counts(a, b):
    compatible(a, b)
    k, e, d, f, g = a.palette, a.equal, a.different, b.equal, b.different
    return Counts(k, e*f+(k-1)*d*g, e*g+d*f+(k-2)*d*g)


def matrix(c):
    audit_counts(c)
    k, e, d = c.palette, c.equal, c.different
    return ((e, (k-1)*d), (d, e+(k-2)*d))


def multiply(a, b):
    return tuple(tuple(sum(a[i][z]*b[z][j] for z in range(2))
                       for j in range(2)) for i in range(2))


def different_probability(c):
    audit_counts(c)
    z = c.equal+(c.palette-1)*c.different
    p.need(z > 0, 'nonempty colouring law')
    return F((c.palette-1)*c.different, z)


def spectral(c):
    audit_counts(c)
    z = c.equal+(c.palette-1)*c.different
    p.need(z > 0, 'nonempty colouring law')
    return F(c.equal-c.different, z)


@dataclass(frozen=True)
class TwoPort:
    graph: p.Graph
    ports: tuple


def audit_gadget(g):
    p.need(type(g) is TwoPort, 'typed two-port graph')
    p.check_domain(g.graph, g.ports, 2)
    p.need(len(g.ports) == 2, 'exactly two distinct retained ports')


def literal_counts(g, k):
    audit_gadget(g)
    p.need(type(k) is int and k >= 2, 'exact palette at least two')
    return Counts(k, old.conditional_count(g.graph, g.ports, k, (0, 0)),
                  old.conditional_count(g.graph, g.ports, k, (0, 1)))


def witness(g, k, equal):
    """Find and authenticate a literal extension; None means exhaustive failure."""
    audit_gadget(g)
    p.need(type(k) is int and k >= 2 and type(equal) is bool, 'exact query types')
    values = {g.ports[0]: 0, g.ports[1]: 0 if equal else 1}
    interior = tuple(v for v in g.graph.vertices if v not in values)
    for cs in product(range(k), repeat=len(interior)):
        values.update(zip(interior, cs))
        answer = tuple(values[v] for v in g.graph.vertices)
        if p.valid_colouring(g.graph, k, answer):
            return answer
    return None


def compose_gadgets(a, b, operation):
    """Fresh-copy interiors; a serial join shares just the intermediate port."""
    audit_gadget(a)
    audit_gadget(b)
    p.need(type(operation) is str and operation in ('S', 'P'), 'serial/parallel tag')
    if operation == 'P':
        ma, mb, first = dict(zip(a.ports, (0, 1))), dict(zip(b.ports, (0, 1))), 2
    else:
        ma, mb, first = dict(zip(a.ports, (0, 2))), dict(zip(b.ports, (2, 1))), 3
    for g, mapping in ((a, ma), (b, mb)):
        for v in g.graph.vertices:
            if v not in mapping:
                mapping[v] = first
                first += 1
    edges = frozenset(tuple(sorted((mapping[u], mapping[v])))
                      for g, mapping in ((a, ma), (b, mb)) for u, v in g.graph.edges)
    return TwoPort(p.Graph(tuple(range(first)), edges), (0, 1))


def compile_tree(tree, k):
    """Authenticate leaves and export actual marked graph plus exact counts.

    Grammar: ('leaf', TwoPort) or ('S'/'P', left, right). This composes
    abstract graphs; it does not authenticate simultaneous plane embeddings.
    """
    p.need(type(k) is int and k >= 2, 'exact palette')
    p.need(type(tree) is tuple and len(tree) >= 1 and type(tree[0]) is str,
           'typed finite expression')
    tag = tree[0]
    if tag == 'leaf':
        p.need(len(tree) == 2, 'leaf arity')
        return tree[1], literal_counts(tree[1], k)
    p.need(tag in ('S', 'P') and len(tree) == 3, 'composition arity')
    a, ca = compile_tree(tree[1], k)
    b, cb = compile_tree(tree[2], k)
    graph = compose_gadgets(a, b, tag)
    return graph, (series_counts(ca, cb) if tag == 'S' else parallel_counts(ca, cb))


def edge():
    return TwoPort(p.Graph((0, 1), frozenset({(0, 1)})), (0, 1))


def fan(r):
    natural(r, 'exact common-neighbour count')
    return TwoPort(p.Graph(tuple(range(r+2)),
                          frozenset((u, v) for u in (0, 1) for v in range(2, r+2))), (0, 1))


def fan_counts(k, r):
    p.need(type(k) is int and k >= 2, 'exact palette')
    natural(r, 'exact fan size')
    return Counts(k, (k-1)**r, (k-2)**r)


def double_cone(interior):
    p.audit_graph(interior)
    mapping = {v: i+2 for i, v in enumerate(interior.vertices)}
    edges = frozenset(tuple(sorted((mapping[u], mapping[v]))) for u, v in interior.edges)
    edges |= frozenset((u, v) for u in (0, 1) for v in mapping.values())
    return TwoPort(p.Graph(tuple(range(len(mapping)+2)), edges), (0, 1))


def colouring_count(graph, k):
    p.audit_graph(graph)
    natural(k, 'exact possibly empty palette')
    return sum(all(cs[graph.vertices.index(u)] != cs[graph.vertices.index(v)]
                   for u, v in graph.edges)
               for cs in product(range(k), repeat=len(graph.vertices)))


def integer_zero_test(k, row_mass_upper, probability_upper):
    """Sufficient implication only; callers owe truthful bounds for SAME law.

    Returns whether p_D <= probability_upper and Z <= row_mass_upper would
    imply d=0. It does not establish either premise from numerical inputs.
    """
    p.need(type(k) is int and k >= 2, 'exact palette')
    p.need(type(row_mass_upper) is F and row_mass_upper > 0, 'positive rational row bound')
    p.need(type(probability_upper) is F and 0 <= probability_upper <= 1,
           'rational probability upper bound')
    return probability_upper*row_mass_upper < k-1


def small_graphs(n):
    pairs = tuple(combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        yield p.Graph(tuple(range(n)), frozenset(e for i, e in enumerate(pairs) if mask >> i & 1))


def expressions(leaves):
    if leaves == 1:
        return (('leaf', edge()),)
    return tuple((tag, a, b) for left in range(1, leaves)
                 for a in expressions(left) for b in expressions(leaves-left) for tag in ('S', 'P'))


def rejects(call):
    try:
        call()
    except ValueError:
        return True
    return False


def main():
    global CHECKS
    CHECKS = 0
    pair_cases = 0
    for k in range(2, 9):
        for e, d, f, g in product(range(4), repeat=4):
            a, b = Counts(k, e, d), Counts(k, f, g)
            c = series_counts(a, b)
            check(matrix(c) == multiply(matrix(a), matrix(b)), 'serial matrix')
            independent_e = sum((e if z == 0 else d)*(f if z == 0 else g) for z in range(k))
            independent_d = sum((e if z == 0 else d)*(f if z == 1 else g) for z in range(k))
            check((c.equal, c.different) == (independent_e, independent_d), 'named intermediate colours')
            check((c.equal > 0) == ((e*f > 0) or (d*g > 0)), 'equal support')
            check((c.different > 0) == ((e*g > 0) or (d*f > 0) or (k >= 3 and d*g > 0)),
                  'different support')
            if (e+d)*(f+g):
                check(spectral(c) == spectral(a)*spectral(b), 'serial spectral multiplication')
            if k >= 3 and d*g:
                check(c.different > 0 and parallel_counts(a, b).different > 0, 'finite positive closure')
            if d*g and e >= d and f >= g:
                check(1 <= F(c.equal, c.different) <= min(F(e, d), F(f, g)), 'series cannot improve both ratios')
            pair_cases += 1
    expression_cases = 0
    for n in range(1, 5):
        for tree in expressions(n):
            for k in range(2, 7):
                graph, c = compile_tree(tree, k)
                check(c == literal_counts(graph, k), 'actual finite expression')
                if k >= 3:
                    got = witness(graph, k, False)
                    check(got is not None and p.valid_colouring(graph.graph, k, got), 'constructive unequal witness')
                expression_cases += 1
    check(series_counts(Counts(2, 0, 1), Counts(2, 0, 1)) == Counts(2, 1, 0), 'two-colour parity exception')
    for k in range(2, 9):
        for r in range(9):
            c = fan_counts(k, r)
            if r <= 4:
                check(c == literal_counts(fan(r), k), 'fan direct count')
            if k >= 3:
                check(c.different > 0, 'positive finite fan support')
            if k >= 3 and r:
                check(different_probability(c) < different_probability(fan_counts(k, r-1)), 'strict odds amplification')
        for r in (10, 100, 1000):
            c = fan_counts(k, r)
            check(c.equal == (k-1)**r and c.different == (k-2)**r, 'large exact symbolic fan')
    check([different_probability(fan_counts(6, r)) for r in range(3)] == [F(5, 6), F(4, 5), F(16, 21)],
          'unit common-neighbour bank')
    check(F(fan_counts(7, 1).different, fan_counts(7, 1).equal) == F(5, 6), 'seven-colour5/6')
    check(F(fan_counts(8, 1).different, fan_counts(8, 1).equal) == F(6, 7), 'eight-colour6/7')
    check(F(fan_counts(6, 1).different, fan_counts(6, 1).equal) == F(4, 5), 'six-colour4/5')
    # Attach leaves or isolated vertices to a retained centre; both count channels scale.
    base = fan(1)
    base_counts = literal_counts(base, 6)
    for r in range(4):
        new_vertices = tuple(range(3, 3+r))
        for leaves, factor in ((True, 5**r), (False, 6**r)):
            edges = base.graph.edges | (frozenset((0, v) for v in new_vertices) if leaves else frozenset())
            graph = TwoPort(p.Graph(base.graph.vertices+new_vertices, edges), (0, 1))
            got = literal_counts(graph, 6)
            check(got == Counts(6, factor*base_counts.equal, factor*base_counts.different), 'leaf/isolate common scalar')
            check(different_probability(got) == different_probability(base_counts), 'leaf odds unchanged')
    cone_cases = 0
    for n in range(5):
        for interior in small_graphs(n):
            for k in range(2, 7):
                got = literal_counts(double_cone(interior), k)
                check(got.equal == colouring_count(interior, k-1), 'cone equal chromatic count')
                check(got.different == colouring_count(interior, k-2), 'cone different chromatic count')
                cone_cases += 1
    diamond = double_cone(p.Graph((0, 1), frozenset({(0, 1)})))
    check(literal_counts(fan(2), 3) == Counts(3, 4, 1), 'independent interiors fail equality')
    check(literal_counts(diamond, 3) == Counts(3, 2, 0), 'one cross edge creates equality')
    check(literal_counts(TwoPort(diamond.graph, (2, 3)), 3) == Counts(3, 0, 1), 'same graph opposite terminal roles')
    for k in range(2, 7):
        interior = p.Graph(tuple(range(k-1)), frozenset(combinations(range(k-1), 2)))
        check(literal_counts(double_cone(interior), k) == Counts(k, factorial(k-1), 0), 'abstract equality selector')
    spindle = p.spindle()
    positive_edges = []
    for uv in sorted(spindle.edges):
        gadget = TwoPort(p.Graph(spindle.vertices, spindle.edges-{uv}), uv)
        c = literal_counts(gadget, 3)
        check(c.equal > 0 and c.different == 0, 'exact unit-separated three-colour control')
        got = witness(gadget, 3, True)
        check(got is not None and p.valid_colouring(gadget.graph, 3, got), 'positive unit witness')
        positive_edges.append((uv, c.equal))
    degree_cases = 0
    for n in range(2, 5):
        for graph in small_graphs(n):
            for uv in combinations(graph.vertices, 2):
                if uv in graph.edges:
                    continue
                for k in range(2, 5):
                    c = literal_counts(TwoPort(graph, uv), k)
                    if c.equal and not c.different:
                        degree = {v: sum(v in e for e in graph.edges) for v in graph.vertices}
                        check(all(degree[v] >= k-1 for v in uv), 'critical port degree')
                        check(sum(degree[v] >= k for v in graph.vertices if v not in uv) >= k-1,
                              'critical interior core')
                        degree_cases += 1
    # The strict absolute threshold is necessary: d=1 survives equality in the bound.
    check(not integer_zero_test(3, F(6), F(1, 3)), 'probability threshold equality hostile')
    check(integer_zero_test(3, F(6), F(1, 4)), 'conditional strict integer threshold')
    for r in range(1, 15):
        c = fan_counts(6, r)
        check(not integer_zero_test(6, F(c.equal+5*c.different), different_probability(c)),
              'truthful fan bound never proves missing support')
    bad = (
        lambda: audit_counts(Counts(True, 1, 1)), lambda: audit_counts(Counts(6, True, 1)),
        lambda: audit_counts(Counts(6, 1, 1.0)), lambda: audit_counts(Counts(6, -1, 1)),
        lambda: series_counts(Counts(6, 1, 1), Counts(7, 1, 1)),
        lambda: literal_counts(TwoPort(edge().graph, (0, True)), 6),
        lambda: literal_counts(TwoPort(edge().graph, (0, 0)), 6),
        lambda: literal_counts(TwoPort(edge().graph, (0, 2)), 6),
        lambda: compile_tree(('S', ('leaf', edge())), 6),
        lambda: compile_tree(('leaf', edge()), 6.0), lambda: fan(True),
        lambda: witness(fan(1), 6, 1), lambda: different_probability(Counts(6, 0, 0)),
        lambda: integer_zero_test(6, 1.0, F(0)),
    )
    for call in bad:
        check(rejects(call), 'malformed public input')
    print('existence_count_amplifiers_20261011e')
    print('scope: exact finite two-port counts; no six-colour plane tester')
    print(f'pair cases={pair_cases}; edge expression/palette cases={expression_cases}; cone cases={cone_cases}')
    print(f'critical degree positive controls={degree_cases}; malformed controls={len(bad)}')
    print(f'unit-separated k3 controls={positive_edges}')
    print('k6 common-neighbour bank p_D: 5/6,4/5,16/21; k7 ratio5/6; k8 ratio6/7')
    print('diamond opposite-port counts=(2,0), adjacent-port counts=(0,1); k6 abstract selector=(120,0)')
    print(f'exact checks={CHECKS}')


if __name__ == '__main__':
    main()
