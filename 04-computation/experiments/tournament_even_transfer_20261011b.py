"""Typed even-graph counts, exact finite-group colouring flows and CRT sidecars.

No six/seven-chromatic unit-distance witness is claimed.  The finite controls
test standard identities and exact reconstruction, not a new counting theorem.
"""
from collections import Counter, deque
from fractions import Fraction as F
from itertools import combinations, permutations, product
import json

CHECKS = 0


def need(ok, why):
    if not ok:
        raise ValueError(why)


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(why)


def graph(n, edges):
    need(type(n) is int and n >= 0, 'exact nonnegative graph order')
    need(type(edges) is tuple, 'ordered edge tuple')
    for e in edges:
        need(type(e) is tuple and len(e) == 2 and all(type(v) is int for v in e),
             'exact edge labels')
        need(0 <= e[0] < e[1] < n, 'simple edge in canonical orientation')
    need(tuple(sorted(set(edges))) == edges, 'sorted distinct edges')
    return edges


def modulus(q):
    need(type(q) is int and q >= 2, 'exact palette modulus at least two')


def vector(values, length, q):
    need(type(values) is tuple and len(values) == length and
         all(type(v) is int and 0 <= v < q for v in values), 'canonical residue vector')


def forest(n, edges):
    graph(n, edges)
    adj = [[] for _ in range(n)]
    for i, (u, v) in enumerate(edges):
        adj[u].append((v, i)); adj[v].append((u, i))
    parents = [None]*n
    order, roots = [], []
    seen = set()
    for root in range(n):
        if root in seen:
            continue
        roots.append(root); seen.add(root)
        queue = deque([root])
        while queue:
            u = queue.popleft(); order.append(u)
            for v, index in adj[u]:
                if v not in seen:
                    seen.add(v); parents[v] = (u, index); queue.append(v)
    tree = {pair[1] for pair in parents if pair is not None}
    chords = tuple(i for i in range(len(edges)) if i not in tree)
    return tuple(parents), tuple(order), tuple(roots), chords


def divergence(n, edges, q, values):
    graph(n, edges); modulus(q); vector(values, len(edges), q)
    result = [0]*n
    for (u, v), value in zip(edges, values):
        result[u] -= value; result[v] += value
    return tuple(v % q for v in result)


def complete_flow(n, edges, q, chord_values):
    """Bijection (Z/q)^beta -> Kirchhoff flows; works for composite q."""
    parents, order, roots, chords = forest(n, edges); modulus(q)
    vector(chord_values, len(chords), q)
    values = [0]*len(edges)
    for i, v in zip(chords, chord_values):
        values[i] = v
    debt = list(divergence(n, edges, q, tuple(values)))
    for vertex in reversed(order):
        if parents[vertex] is None:
            continue
        parent, index = parents[vertex]
        sign = 1 if edges[index][1] == vertex else -1
        values[index] = -sign*debt[vertex] % q
        debt[parent] = (debt[parent]+debt[vertex]) % q
        debt[vertex] = 0
    need(all(debt[r] == 0 for r in roots), 'component conservation')
    return tuple(values)


def all_flows(n, edges, q):
    _, _, _, chords = forest(n, edges); modulus(q)
    return tuple(complete_flow(n, edges, q, v)
                 for v in product(range(q), repeat=len(chords)))


def chromatic_from_flows(n, edges, q):
    """Exact signed finite Fourier sum; existence of a flow alone is insufficient."""
    graph(n, edges); modulus(q)
    counts = Counter(sum(x != 0 for x in flow) for flow in all_flows(n, edges, q))
    signed = sum(number*(-1)**s*(q-1)**(len(edges)-s) for s, number in counts.items())
    result = F(q**n*signed, q**len(edges))
    need(result.denominator == 1 and result >= 0, 'integer colouring count')
    return int(result), tuple(sorted(counts.items()))


def chromatic_direct(n, edges, q):
    graph(n, edges); modulus(q)
    return sum(all(colours[u] != colours[v] for u, v in edges)
               for colours in product(range(q), repeat=n))


def differences(n, edges, q, colours):
    graph(n, edges); modulus(q); vector(colours, n, q)
    return tuple((colours[v]-colours[u]) % q for u, v in edges)


def colour_from_cut(n, edges, q, edge_values, root_colours):
    """Recover a proper colouring iff the supplied nonzero edge vector is a cut.

    Retains one root colour for every connected component, including isolates.
    A Kirchhoff flow is not automatically a cut vector.
    """
    parents, order, roots, _ = forest(n, edges); modulus(q)
    vector(edge_values, len(edges), q); vector(root_colours, len(roots), q)
    need(all(edge_values), 'nowhere-zero edge differences required')
    colours = [None]*n
    for root, colour in zip(roots, root_colours):
        colours[root] = colour
    for v in order:
        if parents[v] is not None:
            u, index = parents[v]
            sign = 1 if edges[index][1] == v else -1
            colours[v] = (colours[u]+sign*edge_values[index]) % q
    result = tuple(colours)
    need(differences(n, edges, q, result) == edge_values, 'all chord cycle sums vanish')
    return result


def crt6(binary, ternary):
    need(type(binary) is tuple and type(ternary) is tuple and len(binary) == len(ternary),
         'same marked coordinates required')
    vector(binary, len(binary), 2); vector(ternary, len(ternary), 3)
    return tuple((3*a+4*b) % 6 for a, b in zip(binary, ternary))


def chromatic_six_crt(n, edges):
    graph(n, edges)
    signed = 0
    for a in all_flows(n, edges, 2):
        for b in all_flows(n, edges, 3):
            support = sum(x != 0 or y != 0 for x, y in zip(a, b))
            signed += (-1)**support*5**(len(edges)-support)
    result = F(6**n*signed, 6**len(edges))
    need(result.denominator == 1 and result >= 0, 'CRT colouring count')
    return int(result)


def renamed_graph(n, edges, perm):
    graph(n, edges)
    need(type(perm) is tuple and all(type(x) is int for x in perm) and
         sorted(perm) == list(range(n)), 'exact vertex permutation')
    return tuple(sorted(tuple(sorted((perm[u], perm[v]))) for u, v in edges))


def relabel_flow(n, edges, q, values, perm):
    renamed = renamed_graph(n, edges, perm); modulus(q); vector(values, len(edges), q)
    mapped = {}
    for (u, v), value in zip(edges, values):
        a, b = perm[u], perm[v]
        mapped[tuple(sorted((a, b)))] = value if a < b else -value % q
    return renamed, tuple(mapped[e] for e in renamed)


def automorphism_even(n, edges):
    graph(n, edges)
    return all(sum(perm[u] > perm[v] for u, v in edges) % 2 == 0
               for perm in permutations(range(n)) if renamed_graph(n, edges, perm) == edges)


def eulerian(n, edges):
    graph(n, edges)
    return not any(divergence(n, edges, 2, (1,)*len(edges)))


def graph_mask(n, edges):
    graph(n, edges); slots = {e: i for i, e in enumerate(combinations(range(n), 2))}
    return sum(1 << slots[e] for e in edges)


def graph_canonical(n, edges):
    return min(graph_mask(n, renamed_graph(n, edges, p)) for p in permutations(range(n)))


def tournament_relabel(n, mask, perm):
    # Internal finite reader: bit one reverses the canonical low-to-high arc.
    edges = tuple(combinations(range(n), 2)); slots = {e: i for i, e in enumerate(edges)}
    out = 0
    for i, (u, v) in enumerate(edges):
        a, b = perm[u], perm[v]
        bit = ((mask >> i) & 1) ^ (a > b)
        out |= bit << slots[tuple(sorted((a, b)))]
    return out


def small_count(n):
    slots = tuple(combinations(range(n), 2)); perms = tuple(permutations(range(n)))
    graphs = tuple(tuple(e for i, e in enumerate(slots) if mask >> i & 1)
                   for mask in range(1 << len(slots)))
    even = tuple(g for g in graphs if automorphism_even(n, g))
    euler = tuple(g for g in graphs if eulerian(n, g))
    tour_classes = {min(tournament_relabel(n, mask, p) for p in perms)
                    for mask in range(1 << len(slots))}
    even_classes = {graph_canonical(n, g) for g in even}
    for perm in perms:
        tour_fixed = sum(tournament_relabel(n, mask, perm) == mask
                         for mask in range(1 << len(slots)))
        signed_fixed = sum((-1)**sum(perm[u] > perm[v] for u, v in g)
                           for g in graphs if renamed_graph(n, g, perm) == g)
        check(tour_fixed == signed_fixed, 'signed fixed-point character identity')
    return {'n': n, 'labelled_tournaments': 1 << len(slots),
            'labelled_automorphism_even': len(even), 'labelled_eulerian': len(euler),
            'tournament_types': len(tour_classes), 'automorphism_even_types': len(even_classes),
            'eulerian_types': len({graph_canonical(n, g) for g in euler})}


def cycle_basis_image(n, cycles):
    """Finite source-marked control for THM4069; not a basis-free constructor."""
    span = {frozenset()}
    for cycle in cycles:
        span |= {g ^ frozenset(cycle) for g in tuple(span)}
    canonical = {g: graph_canonical(n, tuple(sorted(g))) for g in span}
    image = set()
    for g in span:
        for cycle in cycles:
            a, b = canonical[g], canonical[g ^ frozenset(cycle)]
            if a != b:
                image.add(tuple(sorted((a, b))))
    return tuple(sorted(set(canonical.values()))), tuple(sorted(image))


def reject(call):
    try:
        call()
    except ValueError:
        check(True, 'malformed or insufficient certificate rejected')
        return
    check(False, 'rejection required')


def main():
    table = [small_count(n) for n in range(1, 5)]
    check([r['tournament_types'] for r in table] == [1, 1, 2, 4], 'unlabelled tournaments')
    check(all(r['tournament_types'] == r['automorphism_even_types'] for r in table),
          'correct unlabelled count equality')
    check(table[2]['labelled_tournaments'] == 8 and table[2]['labelled_automorphism_even'] == 4,
          'labelled sizes differ')
    check(table[3]['eulerian_types'] == 3, 'Eulerian types are not the tournament sequence')
    triangle = ((0, 1), (0, 2), (1, 2))
    path = ((0, 1), (1, 2))
    check(eulerian(3, triangle) and not automorphism_even(3, triangle), 'Eulerian not sign-even')
    check(not eulerian(3, path) and automorphism_even(3, path), 'sign-even not Eulerian')
    check(not eulerian(3, ((0, 1),)) and not automorphism_even(3, ((0, 1),)), 'neither even')
    check(eulerian(3, ()) and automorphism_even(3, ()), 'both even')
    check(renamed_graph(2, (), (1, 0)) == (), 'empty graph has full stabilizer')
    check(all(tournament_relabel(2, mask, (1, 0)) != mask for mask in (0, 1)),
          'no equivariant image of the empty graph in tournaments')
    star = (triangle, ((0, 1), (0, 3), (1, 3)), ((0, 2), (0, 3), (2, 3)))
    path_basis = (triangle, ((0, 1), (0, 3), (1, 2), (2, 3)),
                  ((1, 2), (1, 3), (2, 3)))
    sv, se = cycle_basis_image(4, star); pv, pe = cycle_basis_image(4, path_basis)
    check(sv == pv and len(sv) == 3 and len(se) == 2 and len(pe) == 3,
          'same Eulerian quotient vertices, P3 versus K3 adjacency')
    # Every finite graph on <=4 vertices, q=2..7: two independent count paths.
    universes = 0
    for n in range(1, 5):
        possible = tuple(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            edges = tuple(e for i, e in enumerate(possible) if mask >> i & 1)
            _, _, roots, chords = forest(n, edges)
            check(len(chords) == len(edges)-n+len(roots), 'exact flow dimension')
            for q in range(2, 8):
                flows = all_flows(n, edges, q)
                check(len(set(flows)) == q**len(chords), 'bijective composite flow parametrization')
                for flow in flows:
                    check(not any(divergence(n, edges, q, flow)), 'independent Kirchhoff check')
                count, _ = chromatic_from_flows(n, edges, q)
                check(count == chromatic_direct(n, edges, q), 'signed sum versus literal colourings')
                if q == 6:
                    check(chromatic_six_crt(n, edges) == count, 'CRT support-union identity')
                universes += 1
    # Every proper colouring of two marked graphs, including a disconnected graph.
    lift_controls = 0
    for n, edges in ((3, triangle), (5, ((0, 1), (1, 2), (3, 4)))):
        _, _, roots, _ = forest(n, edges)
        for q in (2, 3, 6, 7):
            for colours in product(range(q), repeat=n):
                if any(colours[u] == colours[v] for u, v in edges):
                    continue
                values = differences(n, edges, q, colours)
                got = colour_from_cut(n, edges, q, values, tuple(colours[r] for r in roots))
                check(got == colours, 'exact colour reconstruction with root phases')
                lift_controls += 1
    # q=6 needs overlap, not just separate support enumerators.
    bowtie = ((0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (3, 4))
    a = (1, 1, 0, 0, 1, 0)
    b_same = (1, 2, 0, 0, 1, 0)
    b_other = (0, 0, 1, 2, 0, 1)
    for b in (b_same, b_other):
        check(not any(divergence(5, bowtie, 2, a)) and
              not any(divergence(5, bowtie, 3, b)), 'CRT hostile is genuinely native')
        joint = crt6(a, b)
        check(not any(divergence(5, bowtie, 6, joint)), 'CRT joined flow')
    check(sum(x != 0 for x in crt6(a, b_same)) == 3 and
          sum(x != 0 for x in crt6(a, b_other)) == 6, 'same marginals opposite-sign joint weights')
    check(F(-1, 5)**3 == F(-1, 125) and F(-1, 5)**6 == F(1, 15625), 'exact overlap hostile weights')
    check(chromatic_six_crt(3, triangle) == 120 and chromatic_direct(3, triangle, 2) == 0,
          'proper six-colouring is not a proper binary and proper ternary pair')
    flow = complete_flow(3, triangle, 6, (1,))
    check(all(flow) and not any(divergence(3, triangle, 6, flow)), 'nowhere-zero flow hostile')
    reject(lambda: colour_from_cut(3, triangle, 6, flow, (0,)))
    for p in permutations(range(3)):
        es, f = relabel_flow(3, triangle, 6, flow, p)
        check(not any(divergence(3, es, 6, f)), 'orientation-covariant transport')
    # Both even notions coexist with arbitrarily large abstract chromatic number.
    for k in (2, 3):
        n = 2*k
        edges = tuple((u, v) for u, v in combinations(range(n), 2) if u//2 != v//2)
        check(eulerian(n, edges) and automorphism_even(n, edges), 'double-part graph both even')
        colours = tuple(v//2 for v in range(n))
        vals = differences(n, edges, k, colours)
        check(colour_from_cut(n, edges, k, vals, (0,)) == colours, 'retained part colouring')
        check(all((2*u, 2*v) in edges for u, v in combinations(range(k), 2)), 'matching clique lower bound')
    # Exact colours 6 and 7 are reconstructed, not inferred from an even count.
    for q in (6, 7):
        complete = tuple(combinations(range(q), 2)); colours = tuple(range(q))
        check(colour_from_cut(q, complete, q, differences(q, complete, q, colours), (0,)) == colours,
              'marked q-clique reconstruction, no unit-distance claim')
    reject(lambda: graph(True, ()))
    reject(lambda: graph(3, ((False, 1),)))
    reject(lambda: complete_flow(3, triangle, 6.0, (1,)))
    reject(lambda: complete_flow(3, triangle, 6, (True,)))
    reject(lambda: colour_from_cut(3, triangle, 3, (1, 2, 1), (True,)))
    reject(lambda: colour_from_cut(3, triangle, 3, (1, 2, 1), ()))
    reject(lambda: colour_from_cut(3, triangle, 3, (1, 0, 1), (0,)))
    reject(lambda: crt6((1, 0), (1,)))
    reject(lambda: relabel_flow(3, triangle, 6, flow, (False, 1, 2)))
    print(json.dumps({'scope': 'counts versus marked finite colouring certificates; no plane six/seven witness',
                      'checks': CHECKS, 'small_count_table': table,
                      'graph_palette_universe': universes, 'proper_colour_lifts': lift_controls,
                      'triangle_six_colourings': 120,
                      'crt_same_support_weights': ['-1/125', '1/15625'],
                      'critical_flow_dimension': 'beta >= (k-3)n/2+1 for connected k-critical graphs'}, indent=2))


if __name__ == '__main__':
    main()
