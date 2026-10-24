#!/usr/bin/env python3
"""Independent boundary audit. No imports from either production packet."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations
from pathlib import Path
from math import lcm
import json

ROOT = Path(__file__).resolve().parent


def need(test, message):
    if not test:
        raise RuntimeError(message)


def pairs(n):
    return [(a, b) for a in range(n) for b in range(a + 2, n)]


def matrix(n, x):
    a = [[False] * n for _ in range(n)]
    for u in range(n - 1):
        a[u][u + 1] = True
    for k, (u, v) in enumerate(pairs(n)):
        a[u][v] = bool(x & (1 << k))
        a[v][u] = not a[u][v]
    return a


def reflect(n, x):
    edges = pairs(n)
    position = {p: i for i, p in enumerate(edges)}
    return sum(((x >> i) & 1) << position[n - 1 - v, n - 1 - u]
               for i, (u, v) in enumerate(edges))


def path(a, p):
    return all(a[u][v] for u, v in zip(p, p[1:]))


def word(a, p):
    return sum(int(a[p[u]][p[v]]) << k for k, (u, v) in enumerate(pairs(len(a))))


def auto(a, p):
    return all(a[u][v] == a[p[u]][p[v]] for u in range(len(a)) for v in range(len(a)))


def moves(p, width=None):
    result = set()
    for old in range(len(p)):
        rest = p[:old] + p[old + 1:]
        for new in range(len(p)):
            if old != new and (width is None or abs(old - new) + 1 <= width):
                result.add(rest[:new] + (p[old],) + rest[new:])
    return result


def graph_components(vertices, neighbors):
    unseen = set(vertices)
    result = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        component, todo = {root}, [root]
        while todo:
            for y in neighbors(todo.pop()):
                if y in unseen:
                    unseen.remove(y)
                    component.add(y)
                    todo.append(y)
        result.append(component)
    return result


def native_census(n):
    """Recover classes by literal Hamiltonian permutations, not canonicalization."""
    N = 1 << len(pairs(n))
    all_perms = list(permutations(range(n)))
    unseen, classes, owner = set(range(N)), [], {}
    while unseen:
        x = min(unseen)
        a = matrix(n, x)
        paths = {p for p in all_perms if path(a, p)}
        fibres = Counter(word(a, p) for p in paths)
        autos = {p for p in all_perms if auto(a, p)}
        need(set(fibres.values()) == {len(autos)}, 'all observer fibres have full Aut size')
        need(set(fibres) <= unseen, 'isomorphism fibres overlap')
        unseen -= set(fibres)
        for y in fibres:
            owner[y] = len(classes)
        classes.append({'root': x, 'matrix': a, 'paths': paths, 'words': set(fibres), 'autos': autos})
    merged = Counter()
    blue = Counter()
    wordnode = {}
    for x in range(N):
        key = tuple(sorted({owner[x], owner[reflect(n, x)]}))
        wordnode[x] = key
        merged[key] += 1
        if reflect(n, x) == x:
            blue[key] += 1
    loops = Counter()
    for x in range(N):
        cx = x ^ (N - 1)
        if x < cx and wordnode[x] == wordnode[cx]:
            loops['blue' if reflect(n, x) == x else 'black'] += 1
    histogram = Counter((d, blue[u]) for u, d in merged.items())
    return classes, {'n': n, 'A': len(classes), 'V': len(merged), 'S': len(blue),
                     'N': N, 'B': sum(blue.values()), 'M2': sum(d*d for d in merged.values()),
                     'blue_loops': loops['blue'], 'black_loops': loops['black'],
                     'histogram': [[d, b, c] for (d, b), c in sorted(histogram.items())]}


def signature(degrees, blue, black):
    """Canonical coloured multigraph signature, preserving loops/multiplicities."""
    count = len(degrees)
    return min(tuple((blue[tuple(sorted((p[i], p[j])))], black[tuple(sorted((p[i], p[j])))])
                     for i in range(count) for j in range(i, count))
               for p in permutations(range(count)))


def partitions(size):
    """All set partitions encoded in restricted-growth order."""
    def visit(prefix, maximum):
        if len(prefix) == size:
            yield prefix
        else:
            for v in range(maximum + 2):
                yield from visit(prefix + (v,), max(maximum, v))
    yield from visit((0,), 0)


def tiny_realization_audit():
    """Compare every odd R-stable partition of X4 to every eligible graph."""
    n, N = 4, 8
    actual = set()
    ordinary_count = 0
    for ordinary in partitions(N):
        masses = Counter(ordinary)
        if any(v % 2 == 0 for v in masses.values()):
            continue
        image = {}
        for x, c in enumerate(ordinary):
            image.setdefault(c, set()).add(ordinary[reflect(n, x)])
        if any(len(v) != 1 for v in image.values()):
            continue
        ordinary_count += 1
        keys = sorted({tuple(sorted((c, next(iter(image[c]))))) for c in masses})
        idx = {key: i for i, key in enumerate(keys)}
        merged = [idx[tuple(sorted((c, next(iter(image[c])))))] for c in ordinary]
        blue, black, degrees = Counter(), Counter(), Counter(merged)
        for x in range(N // 2):
            edge = tuple(sorted((merged[x], merged[x ^ (N - 1)])))
            (blue if reflect(n, x) == x else black)[edge] += 1
        need(all(v % 2 == 0 for v in black.values()), 'tiny black twins')
        actual.add(signature(degrees, blue, black))
    candidates = set()
    loop_cases = 0
    # Two blue edges plus one K edge have six incidence slots.
    for endpoint_partition in partitions(6):
        blue, half, bdeg, kdeg = Counter(), Counter(), Counter(), Counter()
        for slot in (0, 2, 4):
            u, v = endpoint_partition[slot:slot + 2]
            target, degree = (blue, bdeg) if slot < 4 else (half, kdeg)
            target[tuple(sorted((u, v)))] += 1
            degree[u] += 1
            degree[v] += 1
        vertices = set(endpoint_partition)
        if any((bdeg[u] and bdeg[u] % 2 == 0) or
               (not bdeg[u] and kdeg[u] % 2 == 0) for u in vertices):
            continue
        black = Counter({e: 2*c for e, c in half.items()})
        candidates.add(signature(vertices, blue, black))
        loop_cases += any(u == v for u, v in blue) or any(u == v for u, v in half)
    need(actual == candidates, 'n4 realization theorem exhaustive equivalence')
    need(loop_cases > 0, 'tiny exhaustive audit did not cover loops')
    return {'n': 4, 'odd_R_stable_ordinary_partitions': ordinary_count,
            'coloured_graph_isomorphism_types': len(actual), 'eligible_loop_cases': loop_cases}


def verify_model(g):
    n = g['n']
    N = 1 << len(pairs(n))
    merged = g['merged_by_word']
    ordinary = [tuple(c) for c in g['ordinary_by_word']]
    need(len(merged) == len(ordinary) == N, 'model word carrier')
    mass, ordinary_mass = Counter(merged), Counter(ordinary)
    need(all(v % 2 for v in ordinary_mass.values()), 'odd ordinary fibres')
    image = defaultdict(set)
    bdeg, blue, black = Counter(), Counter(), Counter()
    scores = defaultdict(set)
    for x in range(N):
        rx, cx = reflect(n, x), x ^ (N - 1)
        need(rx != cx and merged[x] == merged[rx], 'actual free CR and merged R action')
        image[ordinary[x]].add(ordinary[rx])
        need(ordinary[x][0] == merged[x], 'ordinary node ownership')
        if rx == x:
            bdeg[merged[x]] += 1
        if x < cx:
            edge = tuple(sorted((merged[x], merged[cx])))
            (blue if rx == x else black)[edge] += 1
        scores[ordinary[x]].add(tuple(sorted(map(sum, matrix(n, x)))))
    need(all(len(v) == 1 for v in image.values()), 'ordinary R-stability')
    for u in mass:
        sheets = {c for c in ordinary_mass if c[0] == u}
        expected = {(u, 0)} if g['self_converse'][u] else {(u, 0), (u, 1)}
        need(sheets == expected, 'exact ordinary sheets at merged vertex')
        need(bool(bdeg[u]) == g['self_converse'][u], 'positive blue iff SC')
        for c in sheets:
            target = c if g['self_converse'][u] else (u, 1-c[1])
            need(image[c] == {target}, 'ordinary reflection fixes or exchanges sheets')
    need(any(len(v) > 1 for v in scores.values()), 'required nonnative score witness')
    need(all(v % 2 == 0 for v in black.values()), 'literal black multiplicities even')
    computed_degree = Counter()
    for (u, v), c in (blue + black).items():
        computed_degree[u] += c
        computed_degree[v] += c
    need(computed_degree == mass, 'loops twice and total degree equals merged mass')
    need(blue == Counter({(u, v): c for u, v, c in g['blue_edges']}), 'saved blue edges')
    need(black == Counter({(u, v): 2*c for u, v, c in g['black_half_edges']}), 'saved K edges')
    need([mass[u] for u in range(len(mass))] == g['degrees'], 'saved degree vector')
    for u, d in mass.items():
        need((bdeg[u] % 2 == 1 and d % 2 == 1) or (bdeg[u] == 0 and d % 4 == 2), 'node parity')
    V, S, A = len(mass), len(bdeg), len(ordinary_mass)
    need(A == 2*V-S, 'ordinary count identity')
    def connected(edges):
        vertices = {u for e in edges for u in e}
        adj = defaultdict(set)
        for u, v in edges:
            adj[u].add(v)
            adj[v].add(u)
        return len(graph_components(vertices, lambda u: adj[u])) == 1
    return {'V': V, 'S': S, 'A': A, 'N': N, 'B': sum(bdeg.values()),
            'M2': sum(d*d for d in mass.values()),
            'histogram': Counter((d, bdeg[u]) for u, d in mass.items()),
            'blue_loops': sum(c for (u, v), c in blue.items() if u == v),
            'black_loops': sum(c for (u, v), c in black.items() if u == v),
            'all_layers_connected': all(connected(e) for e in (blue, black, blue + black)),
            'support': set(blue) | set(black)}


def moment_controls(saved):
    result = []
    for r in range(1, 13):
        nodes = list(range(1, 2*r + 2, 2))
        weights = []
        for a in nodes:
            weight = Fraction(1)
            for b in nodes:
                if a != b:
                    weight *= Fraction(2-b, a-b)
            weights.append(weight)
        D = lcm(*(w.denominator for w in weights))
        relation = {2: D, **{a: -int(D*w) for a, w in zip(nodes, weights)}}
        need(D % 2 == 0, 'denominator even')
        need(all(sum(c*d**q for d, c in relation.items()) == 0 for q in range(r + 1)), 'moment annihilation')
        sides = [Counter({d: sign*c for d, c in relation.items() if sign*c > 0}) for sign in (1, -1)]
        N0 = sum(d*c for d, c in sides[0].items())
        if N0 % 2 == 0:
            for side in sides:
                side[1] += 1
            N0 += 1
        S = [sum(c for d, c in side.items() if d % 2) for side in sides]
        n = 3
        while True:
            m, f = (n-1)*(n-2)//2, (n-1)//2
            N, B = 1 << m, 1 << ((m+f)//2)
            if B > 2*max(S) and N > B+2*N0:
                break
            n += 1
        need((N-N0) % 2 == 1 and abs(S[0]-S[1]) == D, 'center parity and S separation')
        for s in S:
            need((B-2*s) > 0 and (B-2*s) % 2 == 0, 'positive integral blue loops')
            need((N-B-2*N0+2*s) > 0 and (N-B-2*N0+2*s) % 4 == 0, 'positive integral K loops')
        if r <= len(saved):
            row = saved[r - 1]
            need(row['signed_relation'] == {str(d): c for d, c in relation.items()}, 'saved exact moment relation')
            need(row['n'] == n and row['S'] == [s+1 for s in S], 'saved chosen order and SC counts')
        result.append({'r': r, 'D': D, 'n': n, 'all_padding_integral': True})
    # Direct hostile check of the arithmetic normalization.
    for d in range(1, 256):
        if d % 4 == 0:
            continue
        need(d*d % 8 == (1 if d % 2 else 4), 'merged second moment congruence')
        for j in range(2, 9):
            need(pow(d, 1 << j, 1 << (j + 2)) == d % 2, 'higher merged degree congruence')
    return result


def inertia_checks(classes, saved):
    n = len(classes[0]['matrix'])
    byroot = {c['root_word']: c for c in saved['components']}
    arrows, raw, lower, partial_count, unsaturated = 0, 0, Fraction(0), 0, 0
    for c in classes:
        a, paths, words, autos = c['matrix'], c['paths'], c['words'], c['autos']
        H, h, aut_order = len(paths), len(words), len(autos)
        pc = graph_components(paths, lambda p: moves(p) & paths)
        need(len(pc) == 1, 'unrestricted relocation path graph disconnected')
        row = byroot[c['root']]
        need((row['word_count'], row['holonomy_order'], row['H']) == (h, aut_order, H), 'producer component census')
        need(all(auto(a, tuple(g)) for g in row['holonomy_generators']), 'saved generators authenticated')
        generated, todo = {tuple(range(n))}, [tuple(range(n))]
        while todo:
            p = todo.pop()
            for q in row['holonomy_generators']:
                product = tuple(p[q[u]] for u in range(n))
                if product not in generated:
                    generated.add(product)
                    todo.append(product)
        need(generated == autos, 'saved generators close to full source Aut')
        # Count normalized directed arrows from source path graph, dividing by Aut.
        path_edges = sum(len(moves(p) & paths) for p in paths)
        need(path_edges % aut_order == 0, 'free action on marked arrows')
        arrows += path_edges // aut_order
        need(row['directed_relocation_edges'] == path_edges // aut_order, 'saved arrow count')
        raw += H
        # Completely independent partial holonomy order: a marked connected lift
        # above b words has b*g paths. Quotient these marked components by words.
        lifted = graph_components(paths, lambda p: moves(p, 3) & paths)
        seen_word_sets = set()
        for component in lifted:
            projected = frozenset(word(a, p) for p in component)
            if projected in seen_word_sets:
                continue
            need(all(not (projected & old) for old in seen_word_sets), 'projected components overlap unequally')
            seen_word_sets.add(projected)
            b = len(projected)
            need(len(component) % b == 0, 'lift degree integral')
            g = len(component) // b
            need(g <= aut_order and aut_order % g == 0, 'lift holonomy subgroup order')
            lower += Fraction(b*g, H)
            partial_count += 1
            unsaturated += b*g != H
        need(set().union(*seen_word_sets) == words, 'partial graph full coverage')
    need(saved['directed_relocation_edges'] == arrows, 'global relocation arrow census')
    if n >= 3:
        need(Fraction((1 << len(pairs(n)))) * Fraction(n*n+n-4, 8) == arrows, 'global arrow polynomial')
    need(saved['raw_closed_reobservation_count'] == raw, 'raw loop count is sum H')
    need(lower <= len(classes) <= partial_count, 'count sandwich')
    return {'n': n, 'classes': len(classes), 'directed_arrows': arrows,
            'literal_path_graphs_connected': True, 'span3_components': partial_count,
            'span3_lower_bound': str(lower), 'span3_unsaturated_components': unsaturated}


def locality_checks():
    result = []
    for n in range(3, 9):
        a, identity = matrix(n, 0), tuple(range(n))
        tried = 0
        for left in range(n):
            for right in range(left+1, n+1):
                if left == 0 and right == n:
                    continue
                for middle in permutations(identity[left:right]):
                    p = identity[:left]+middle+identity[right:]
                    if p != identity:
                        tried += 1
                        need(not path(a, p), 'backward source proper window not frozen')
        need(all(path(a, identity[k:]+identity[:k]) for k in range(n)), 'global rotations missing')
        aut_order = sum(auto(a, p) for p in permutations(range(n)))
        need(aut_order == (3 if n == 3 else 1), 'backward source rigidity')
        result.append({'family': 'T_star', 'n': n, 'window_permutations_checked': tried, 'aut_order': aut_order})
    for k in range(1, 5):
        n, width = 3*k, k+1
        a = [[False]*n for _ in range(n)]
        for u in range(n):
            for v in range(n):
                a[u][v] = (u < v if u//k == v//k else (v//k-u//k) % 3 == 1)
        identity = tuple(range(n))
        need(path(a, identity), 'concatenated inflation observer')
        need(auto(a, tuple((u+k) % n for u in range(n))), 'inflation rotation automorphism')
        tried = 0
        for left in range(n):
            for right in range(left+1, min(n, left+width)+1):
                for middle in permutations(identity[left:right]):
                    p = identity[:left]+middle+identity[right:]
                    if p != identity:
                        tried += 1
                        need(not path(a, p), 'inflation window not frozen')
        result.append({'family': 'C3_TT_k', 'k': k, 'n': n, 'width': width,
                       'window_permutations_checked': tried, 'certified_auto_order_at_least': 3})
    return result


def main():
    cert = json.loads((ROOT/'sources/identifiability/certificate.json').read_text())
    native_saved = json.loads((ROOT/'sources/identifiability/native_reference.json').read_text())
    inertia_saved = json.loads((ROOT/'sources/inertia/normal.json').read_text())
    native, inertia = [], []
    for n in range(1, 7):
        classes, summary = native_census(n)
        inertia.append(inertia_checks(classes, inertia_saved['orders'][n-1]))
        if n >= 3:
            row = native_saved[n-3]
            for key in ('A', 'V', 'S', 'N', 'B', 'M2', 'blue_loops', 'black_loops'):
                need(summary[key] == row[key], f'native reference {key}')
            expected_hist = Counter((v['degree'], v['blue_degree']) for v in row['nodes'])
            need(summary['histogram'] == [[d,b,c] for (d,b),c in sorted(expected_hist.items())], 'native signature histogram')
            native.append(summary)
        print(f'Native and relocation independent checks n={n}: PASS', flush=True)
    for row, independently in zip(inertia_saved['partial_graph_bounds'], inertia[2:]):
        need((row['components'], Fraction(row['inertia_lower_bound']), row['unsaturated_components']) ==
             (independently['span3_components'], Fraction(independently['span3_lower_bound']), independently['span3_unsaturated_components']), 'saved span3 bounds')
    bounded = [verify_model(g) for g in cert['bounded_models']]
    stars = [verify_model(g) for g in cert['models']]
    palette = {(d,b) for d,b,c in native[-1]['histogram']}
    for row in bounded:
        need(set(row['histogram']) <= palette and row['histogram'][1,1] == 2, 'bounded native palette')
        need((row['N'],row['B'],row['V'],row['M2'],row['blue_loops'],row['black_loops']) == (1024,64,34,48724,2,24), 'bounded aggregate and loop totals')
        need(row['all_layers_connected'], 'bounded layer connectivity')
    need([g['S'] for g in bounded] == [12,20], 'bounded S separation')
    need(stars[0]['support'] == stars[1]['support'], 'same star support')
    need([(g['V'],g['M2'],g['S'],g['A']) for g in stars] == [(10,1006066,2,18),(10,1006066,10,10)], 'star collision')
    ambiguity = [s for s in range(0,35,2) if 3*s-4 <= 64 <= 9*s-16 and (s-4*34+3*48724) % 8 == 0]
    need(ambiguity == [12,20], 'exhaustive aggregate ambiguity')
    def public_model(row):
        return {k:v for k,v in row.items() if k not in ('histogram','support')}
    report = {'status': 'PASS', 'native': native, 'inertia': inertia,
              'tiny_realization_exhaustive': tiny_realization_audit(),
              'bounded_models': list(map(public_model, bounded)),
              'same_support_models': list(map(public_model, stars)),
              'exact_relaxed_A_identification_set': [48,56],
              'finite_moment_padding': moment_controls(cert['hierarchy']),
              'locality': locality_checks(),
              'limits': ['Native independent census stops at n=6; no independent n=7 full rerun.',
                         'No claim that unbounded-span relocation connectivity holds at every order.',
                         'Countermodels are relaxed odd-fibre partitions, not native tournament censuses.']}
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
