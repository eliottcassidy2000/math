#!/usr/bin/env python3
"""Exact checks for endpoint quotient fibres and a faithful simple-graph carrier.
No external packages. Python 3.10+. All masks read pairs i<j lexicographically;
bit 1 means i -> j. No search/randomness is needed for the reported witnesses.
"""
from itertools import permutations, product
from collections import Counter, defaultdict
from fractions import Fraction
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def tournament(n, mask):
    a = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            a[i][j] = (mask >> k) & 1
            a[j][i] = 1 - a[i][j]
            k += 1
    return a


def circulant(n, steps):
    return [[int((j - i) % n in steps) for j in range(n)] for i in range(n)]


def ordinal_sum(a, b):
    p, q = len(a), len(b)
    c = [[0] * (p + q) for _ in range(p + q)]
    for i in range(p + q):
        for j in range(p + q):
            c[i][j] = a[i][j] if i < p and j < p else b[i-p][j-p] if i >= p and j >= p else int(i < p)
    return c


def terminal_counts(a):
    n = len(a)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for s in range(1, 1 << n):
        for v in range(n):
            if (s >> v) & 1:
                dp[s][v] += sum(dp[s ^ (1 << v)][u] for u in range(n) if a[u][v] and ((s >> u) & 1))
    return dp[-1]


def brute_paths(a):
    n = len(a)
    return [p for p in permutations(range(n)) if all(a[p[i]][p[i+1]] for i in range(n-1))]


def enumerate_paths_dfs(a):
    n = len(a)
    out = []
    def visit(p, remaining):
        if not remaining:
            out.append(tuple(p))
        else:
            for v in remaining:
                if not p or a[p[-1]][v]:
                    visit(p + [v], remaining - {v})
    visit([], set(range(n)))
    return out


def action_graph_components(a, terminal=None):
    group = tournament_automorphisms(a)
    paths = set(enumerate_paths_dfs(a))
    if terminal is not None:
        paths = {p for p in paths if p[-1] == terminal}
        group = [g for g in group if g[terminal] == terminal]
    histogram = Counter()
    while paths:
        p = min(paths)
        component = {tuple(g[v] for v in p) for g in group}
        require(component <= paths, 'check: component <= paths')
        paths -= component
        histogram[len(component)] += 1
    return dict(sorted(histogram.items()))


def tournament_automorphisms(a):
    n = len(a)
    classes = defaultdict(list)
    for v in range(n):
        classes[sum(a[v])].append(v)
    blocks = list(classes.values())
    out = []
    for images in product(*(permutations(b) for b in blocks)):
        p = [None] * n
        for block, image in zip(blocks, images):
            for v, w in zip(block, image):
                p[v] = w
        if all(a[i][j] == a[p[i]][p[j]] for i in range(n) for j in range(i+1, n)):
            out.append(tuple(p))
    return out


def orbit_representatives(group, n):
    todo, reps = set(range(n)), []
    while todo:
        v = min(todo)
        reps.append(v)
        todo -= {p[v] for p in group}
    return reps


def graph_carrier(a):
    """Core types are implicit in the graph, never fed to the automorphism test."""
    n = len(a)
    graph = [set() for _ in range(n)]
    markers = [0] * n
    def new(marker):
        graph.append(set())
        markers.append(marker)
        return len(graph) - 1
    def edge(v, w):
        graph[v].add(w)
        graph[w].add(v)
    for u in range(n):
        for v in range(n):
            if a[u][v]:
                l, r = new(1), new(2)
                edge(u, l); edge(l, r); edge(r, v)
    core_size = len(graph)
    for v in range(core_size):
        previous = v
        for _ in range(markers[v]):
            w = new(0)
            edge(previous, w)
            previous = w
    return graph


def graph_automorphism_count(g):
    """Independent generic uncoloured graph backtracking after 1-WL refinement."""
    n = len(g)
    color = [len(vs) for vs in g]
    while True:
        signatures = [(color[v], tuple(sorted(color[w] for w in g[v]))) for v in range(n)]
        label = {s: i for i, s in enumerate(sorted(set(signatures)))}
        new = [label[s] for s in signatures]
        if len(set(new)) == len(set(color)):
            color = new
            break
        color = new
    cells = defaultdict(set)
    for v, c in enumerate(color):
        cells[c].add(v)
    mapping, used = {}, set()
    def recurse():
        if len(mapping) == n:
            return 1
        chosen, candidates = None, None
        # Start at a constrained vertex, and propagate along mapped neighbours.
        order = sorted((v for v in range(n) if v not in mapping), key=lambda v: (len(cells[color[v]]), -sum(w in mapping for w in g[v])))
        for v in order:
            cand = [w for w in cells[color[v]] - used if all((x in g[v]) == (y in g[w]) for x, y in mapping.items())]
            if not cand:
                return 0
            if candidates is None or len(cand) < len(candidates):
                chosen, candidates = v, cand
                if len(cand) == 1:
                    break
        total = 0
        for w in candidates:
            mapping[chosen] = w; used.add(w)
            total += recurse()
            used.remove(w); del mapping[chosen]
        return total
    return recurse()


def test_tournament(a, expected_end=None, brute=False, graph=True):
    n = len(a)
    e = terminal_counts(a)
    if expected_end is not None:
        require(e == expected_end, 'check: e == expected_end')
    group = tournament_automorphisms(a)
    boundary = []
    for v in orbit_representatives(group, n):
        stab = sum(p[v] == v for p in group)
        require(e[v] % stab == 0, 'check: e[v] % stab == 0')
        boundary.append({'representative': v, 'stabilizer': stab, 'raw_terminal': e[v], 'quotient_fibre': e[v] // stab})
    require(sum(b['quotient_fibre'] for b in boundary) * len(group) == sum(e), "check: sum(b['quotient_fibre'] for b in boundary) * len(group) == sum(e)")
    if brute:
        paths = brute_paths(a)
        require([sum(p[-1] == v for p in paths) for v in range(n)] == e, 'check: [sum(p[-1] == v for p in paths) for v in range(n)] == e')
        for b in boundary:
            v = b['representative']
            stabgroup = [g for g in group if g[v] == v]
            canonical = {min(tuple(g[x] for x in p) for g in stabgroup) for p in paths if p[-1] == v}
            require(len(canonical) == b['quotient_fibre'], "check: len(canonical) == b['quotient_fibre']")
    row = {'n': n, 'H': sum(e), 'aut': len(group), 'tc': sum(e)//len(group), 'terminal': e, 'boundary_orbits': boundary}
    if graph:
        carrier = graph_carrier(a)
        ga = graph_automorphism_count(carrier)
        leaves = sum(len(vs) == 1 for vs in carrier)
        require(ga == len(group), 'check: ga == len(group)')
        require(len(carrier) == (5*n*n - 3*n)//2, 'check: len(carrier) == (5*n*n - 3*n)//2')
        require(sum(map(len, carrier))//2 == 3*n*(n-1), 'check: sum(map(len, carrier))//2 == 3*n*(n-1)')
        require(leaves == n*(n-1), 'check: leaves == n*(n-1)')
        row.update(graph_vertices=len(carrier), graph_edges=sum(map(len,carrier))//2, graph_aut=ga, graph_leaves=leaves)
    return row


def main():
    output = {}
    maxima = []
    for n in range(1, 6):
        maxima.append(max(sum(terminal_counts(tournament(n, m))) for m in range(1 << (n*(n-1)//2))))
    require(maxima == [1, 1, 3, 5, 15], 'check: maxima == [1, 1, 3, 5, 15]')
    output['exhaustive_H_maxima_n1_to5'] = maxima
    a, b = tournament(6, 26), tournament(7, 290)
    output['endpoint7_minimal'] = test_tournament(a, [1, 5, 7, 3, 2, 1], brute=True)
    output['endpoint21_minimal'] = test_tournament(b, [3, 9, 21, 9, 7, 3, 3], brute=True)
    c3 = circulant(3, {1})
    output['nontrivial_action_endpoint7'] = test_tournament(ordinal_sum(c3, a), [0, 0, 0, 3, 15, 21, 9, 6, 3])
    output['nontrivial_action_endpoint21'] = test_tournament(ordinal_sum(c3, b), [0, 0, 0, 9, 27, 63, 27, 21, 9, 9])
    output['endpoint7_action_graph'] = {'whole': action_graph_components(ordinal_sum(c3, a)), 'terminal5': action_graph_components(ordinal_sum(c3, a), 5)}
    output['endpoint21_action_graph'] = {'whole': action_graph_components(ordinal_sum(c3, b)), 'terminal5': action_graph_components(ordinal_sum(c3, b), 5)}
    require(output['endpoint7_action_graph'] == {'whole': {3: 19}, 'terminal5': {3: 7}}, "check: output['endpoint7_action_graph'] == {'whole': {3: 19}, 'terminal5': {3: 7}}")
    require(output['endpoint21_action_graph'] == {'whole': {3: 55}, 'terminal5': {3: 21}}, "check: output['endpoint21_action_graph'] == {'whole': {3: 55}, 'terminal5': {3: 21}}")
    output['circulant7'] = test_tournament(circulant(7, {1, 2, 3}), [25]*7, brute=True)
    output['paley7'] = test_tournament(circulant(7, {1, 2, 4}), [27]*7, brute=True)
    output['paley_naive_division'] = str(Fraction(27, 21))
    checks = 0
    for n in (3, 4):
        for mask in range(1 << (n*(n-1)//2)):
            test_tournament(tournament(n, mask), graph=True, brute=True)
            checks += 1
    output['exhaustive_graph_checks_orders3_4'] = checks
    print(json.dumps(output, indent=2))

if __name__ == '__main__':
    main()
