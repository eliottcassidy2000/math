#!/usr/bin/env python3
"""Independent finite checks for equitable-cell holonomy and SCC error transfer.

No imports from the marked-growth engine. All assertions use explicit errors
and remain active under python -O. Adjacency rows are outgoing bitsets.
"""
from collections import Counter
from itertools import combinations, permutations, product
from math import factorial, gcd
from pathlib import Path
import hashlib
import json

COUNTS = Counter()


def need(value, label):
    COUNTS[label] += 1
    if not value:
        raise RuntimeError(label)


def arc(rows, v, w):
    return (rows[v] >> w) & 1


def validate(rows):
    n = len(rows)
    need(all(not arc(rows, v, v) for v in range(n)), 'no loops')
    need(all(arc(rows, v, w) + arc(rows, w, v) == 1
             for v in range(n) for w in range(v)), 'tournament')


def induce(rows, vertices):
    return tuple(sum(arc(rows, v, w) << j for j, w in enumerate(vertices))
                 for v in vertices)


def refine(rows):
    """All-equal start; retain old colors, directed out-counts into all colors."""
    n = len(rows)
    colors = (0,) * n
    while True:
        labels = sorted(set(colors))
        signatures = tuple((colors[v], tuple(sum(arc(rows, v, w)
                         for w in range(n) if colors[w] == label)
                         for label in labels)) for v in range(n))
        ids = {s: i for i, s in enumerate(sorted(set(signatures)))}
        new = tuple(ids[s] for s in signatures)
        if len(set(new)) == len(set(colors)):
            return tuple(tuple(v for v in range(n) if new[v] == c)
                         for c in sorted(set(new)))
        colors = new


def aut(rows, unrestricted=False):
    """Literal arc tests, optionally over all n! maps, otherwise degree-pruned.

    No stable-refinement, holonomy, or prescribed triangle block input is used.
    """
    n = len(rows)
    degrees = tuple(row.bit_count() for row in rows)
    cells = [tuple(range(n))] if unrestricted else [
        tuple(v for v in range(n) if degrees[v] == d)
        for d in sorted(set(degrees))]
    result = []
    if unrestricted:
        for p in permutations(range(n)):
            if all(arc(rows, v, w) == arc(rows, p[v], p[w])
                   for v in range(n) for w in range(v)):
                result.append(tuple(p))
    else:
        candidates = {v: tuple(w for w in range(n) if degrees[w] == degrees[v])
                      for v in range(n)}
        order = sorted(range(n), key=lambda v: (len(candidates[v]), v))
        mapping = [-1] * n
        used = set()
        def visit(index):
            if index == n:
                result.append(tuple(mapping))
                return
            v = order[index]
            for w in candidates[v]:
                if w in used or any(arc(rows, v, u) != arc(rows, w, mapping[u])
                                    for u in order[:index]):
                    continue
                mapping[v] = w
                used.add(w)
                visit(index + 1)
                used.remove(w)
                mapping[v] = -1
        visit(0)
    return tuple(result)


def orbits(group, n):
    return tuple(sorted({tuple(sorted({g[v] for g in group})) for v in range(n)}))


def error(rows, group=None):
    cells = refine(rows)
    sizes = {v: len(cell) for cell in cells for v in cell}
    if group is None:
        group = aut(rows)
    return sum(sizes[orbit[0]] - len(orbit)
               for orbit in orbits(group, len(rows)))


def scc(rows):
    """Literal reachability equivalence; ordered by descending external wins."""
    n = len(rows)
    reach = [rows[v] | (1 << v) for v in range(n)]
    for k in range(n):
        for v in range(n):
            if reach[v] >> k & 1:
                reach[v] |= reach[k]
    cells = {tuple(w for w in range(n) if (reach[v] >> w & 1)
                   and (reach[w] >> v & 1)) for v in range(n)}
    return tuple(sorted(cells, key=lambda c: -sum(arc(rows, c[0], w)
                  for w in range(n) if w not in c)))


def join(parts):
    rows = []
    offset = 0
    total = sum(map(len, parts))
    for part in parts:
        later = ((1 << total) - 1) ^ ((1 << (offset + len(part))) - 1)
        rows.extend((row << offset) | later for row in part)
        offset += len(part)
    return tuple(rows)


def check_equitable(rows):
    cells = refine(rows)
    for a in cells:
        need(len(a) % 2 == 1, 'stable cells odd')
        need({sum(arc(rows, v, w) for w in a) for v in a}
             == {(len(a) - 1) // 2}, 'stable cells internally regular')
        for b in cells:
            if a == b:
                continue
            degrees = {sum(arc(rows, v, w) for w in b) for v in a}
            need(len(degrees) == 1, 'equitable crossing rows')
            r = next(iter(degrees))
            incoming = {sum(arc(rows, v, w) for v in a) for w in b}
            need(len(incoming) == 1, 'equitable crossing columns')
            need(len(a) * r == len(b) * next(iter(incoming)), 'crossing balance')
            if gcd(len(a), len(b)) == 1:
                need(r in (0, len(b)), 'coprime crossing uniform')


def check_scc(rows):
    cells = refine(rows)
    components = scc(rows)
    rebuilt = []
    for component in components:
        local = refine(induce(rows, component))
        rebuilt.extend(tuple(component[v] for v in c) for c in local)
    need(set(cells) == set(rebuilt), 'SCC restriction of stable refinement')
    need(error(rows) == sum(error(induce(rows, c)) for c in components),
         'SCC error additivity')


def triangle_coordinates(rows, cell):
    start = min(cell)
    nxt = next(v for v in cell if arc(rows, start, v))
    last = next(v for v in cell if v not in (start, nxt))
    need(arc(rows, nxt, last) and arc(rows, last, start), 'directed triangle cell')
    return (start, nxt, last)


def holonomy(rows):
    cells = refine(rows)
    need(all(len(c) in (1, 3) for c in cells), 'holonomy domain')
    triangles = [triangle_coordinates(rows, c) for c in cells if len(c) == 3]
    edges = []
    for i, j in combinations(range(len(triangles)), 2):
        a, b = triangles[i], triangles[j]
        degree = sum(arc(rows, a[0], w) for w in b)
        if degree in (0, 3):
            continue
        target = 1 if degree == 1 else 0
        sigma = tuple(next(y for y, w in enumerate(b)
                           if arc(rows, v, w) == target) for v in a)
        need(len(set(sigma)) == 3, 'nonuniform block matching')
        epsilon = (sigma[1] - sigma[0]) % 3
        need(epsilon in (1, 2), 'nonzero matching multiplier')
        need(all((sigma[x] - sigma[0]) % 3 == epsilon * x % 3
                 for x in range(3)), 'matching affine coordinates')
        edges.append((i, j, 1 if epsilon == 1 else -1))
    adj = [[] for _ in triangles]
    for i, j, eps in edges:
        adj[i].append((j, eps))
        adj[j].append((i, eps))
    seen = set()
    balanced = []
    unbalanced = []
    for root in range(len(triangles)):
        if root in seen:
            continue
        signs = {root: 1}
        todo = [root]
        good = True
        while todo:
            v = todo.pop()
            for w, eps in adj[v]:
                proposed = signs[v] * eps
                if w not in signs:
                    signs[w] = proposed
                    todo.append(w)
                elif signs[w] != proposed:
                    good = False
        seen.update(signs)
        (balanced if good else unbalanced).append(tuple(sorted(signs)))
    predicted = []
    for values in product(range(3), repeat=len(triangles)):
        if all((values[j] - eps * values[i]) % 3 == 0 for i, j, eps in edges):
            mapping = list(range(len(rows)))
            for triangle, shift in zip(triangles, values):
                for x, vertex in enumerate(triangle):
                    mapping[vertex] = triangle[(x + shift) % 3]
            predicted.append(tuple(mapping))
    need(len(predicted) == 3 ** len(balanced), 'signed kernel dimension')
    return tuple(sorted(predicted)), edges, balanced, unbalanced


def triangle_blocks(specifications):
    """Three C3s; each pair gets uniform or permutation/complement crossing."""
    rows = [0] * 9
    for cell in range(3):
        for x in range(3):
            rows[3 * cell + x] |= 1 << (3 * cell + (x + 1) % 3)
    for (i, j), (degree, sigma) in zip(combinations(range(3), 2), specifications):
        for x in range(3):
            for y in range(3):
                forward = degree == 3 or (degree == 1 and y == sigma[x]) or (
                    degree == 2 and y != sigma[x])
                u, v = 3 * i + x, 3 * j + y
                if forward:
                    rows[u] |= 1 << v
                else:
                    rows[v] |= 1 << u
    return tuple(rows)


def all_tournaments(n):
    pairs = tuple(combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        rows = [0] * n
        for bit, (v, w) in enumerate(pairs):
            if mask >> bit & 1:
                rows[v] |= 1 << w
            else:
                rows[w] |= 1 << v
        yield tuple(rows)


def regular_tournaments(n):
    """Degree-pruned exhaustive labels; pruning only impossible remaining degrees."""
    pairs = tuple(combinations(range(n), 2))
    rows = [0] * n
    degrees = [0] * n
    remaining = [n - 1] * n
    target = (n - 1) // 2
    found = []
    def visit(index):
        if index == len(pairs):
            need(degrees == [target] * n, 'regular leaf degrees')
            found.append(tuple(rows))
            return
        a, b = pairs[index]
        remaining[a] -= 1
        remaining[b] -= 1
        for v, w in ((a, b), (b, a)):
            degrees[v] += 1
            rows[v] |= 1 << w
            if all(degrees[x] <= target <= degrees[x] + remaining[x]
                   for x in (a, b)):
                visit(index + 1)
            rows[v] ^= 1 << w
            degrees[v] -= 1
        remaining[a] += 1
        remaining[b] += 1
    visit(0)
    need(len(found) == len(set(found)), 'regular labels unique')
    return set(found)


def regular_classes(n):
    universe = regular_tournaments(n)
    full_count = len(universe)
    classes = []
    while universe:
        representative = min(universe)
        group = aut(representative, unrestricted=True)
        orbit = {induce(representative, p) for p in permutations(range(n))}
        need(orbit <= universe, 'regular class disjoint exhaustive partition')
        need(len(orbit) * len(group) == factorial(n), 'label orbit stabilizer')
        universe.difference_update(orbit)
        classes.append({'rows': representative, 'automorphism_order': len(group),
                        'labelled_orbit_size': len(orbit),
                        'vertex_orbits': orbits(group, n),
                        'refinement_error': error(representative, group),
                        'strong': len(scc(representative)) == 1})
    return full_count, classes


def main():
    small = Counter()
    for n in range(1, 6):
        for rows in all_tournaments(n):
            validate(rows)
            check_equitable(rows)
            check_scc(rows)
            need(error(rows) == 0, 'all labelled through five exact')
            if all(len(c) in (1, 3) for c in refine(rows)):
                predicted, _, _, unbalanced = holonomy(rows)
                need(predicted == tuple(sorted(aut(rows))), 'small holonomy literal group')
                need(not unbalanced, 'small holonomy no unbalanced component')
            small[n] += 1
    regular = {}
    for n in (1, 3, 5, 7):
        count, classes = regular_classes(n)
        regular[n] = {'labelled_count': count, 'classes': classes}
    need([regular[n]['labelled_count'] for n in (1, 3, 5, 7)] == [1, 2, 24, 2640],
         'regular labelled census sizes')
    need([len(regular[n]['classes']) for n in (1, 3, 5, 7)] == [1, 1, 1, 3],
         'regular class census sizes')
    need(sorted(c['automorphism_order'] for c in regular[7]['classes']) == [3, 7, 21],
         'regular seven group orders')
    need(sum(c['refinement_error'] for c in regular[7]['classes']) == 14,
         'regular seven total error')
    need(all(c['strong'] for c in regular[7]['classes']), 'regular seven strong')

    choices = [(0, (0, 1, 2)), (3, (0, 1, 2))] + [
        (d, p) for d in (1, 2) for p in permutations(range(3))]
    block_counts = Counter()
    group_distribution = Counter()
    error_distribution = Counter()
    for specs in product(choices, repeat=3):
        rows = triangle_blocks(specs)
        validate(rows)
        cells = refine(rows)
        block_counts['constructed'] += 1
        if not all(len(c) in (1, 3) for c in cells):
            need([len(c) for c in cells] == [9], 'excluded block family regular nine')
            block_counts['one_regular_nine_cell'] += 1
            continue
        block_counts['holonomy_domain'] += 1
        check_equitable(rows)
        predicted, edges, balanced, unbalanced = holonomy(rows)
        literal = tuple(sorted(aut(rows)))
        need(predicted == literal, 'block holonomy full literal group')
        expected_error = 6 * sum(map(len, unbalanced))
        need(error(rows, literal) == expected_error, 'holonomy error formula')
        group_distribution[len(literal)] += 1
        error_distribution[expected_error] += 1
        # Gauge covariance: independently flip all eight possible generator choices.
        for tau in product((1, -1), repeat=3):
            for values in product(range(3), repeat=3):
                original = all((values[j] - eps * values[i]) % 3 == 0
                               for i, j, eps in edges)
                changed = all((tau[j] * values[j] - tau[j] * eps * tau[i]
                               * tau[i] * values[i]) % 3 == 0 for i, j, eps in edges)
                need(original == changed, 'signed generator gauge covariance')
    matching_classes = {}
    for specs in product([(1, p) for p in permutations(range(3))], repeat=3):
        rows = triangle_blocks(specs)
        rotations = []
        for shifts in product(range(3), repeat=3):
            order = tuple(3 * cell + (x + shifts[cell]) % 3
                          for cell in range(3) for x in range(3))
            rotations.append(induce(rows, order))
        representative = min(rotations)
        matching_classes[representative] = len(aut(rows))
    need(len(matching_classes) == 16, '216 matchings give sixteen literal classes')
    need(Counter(matching_classes.values()) == Counter({3: 12, 1: 4}),
         'matching classes twelve symmetric four rigid')
    need(sum(error(rows) for rows in matching_classes) == 72,
         'matching family strong error seventy two')
    identity = (0, 1, 2)
    odd = (0, 2, 1)
    hostile = triangle_blocks(((1, identity), (1, identity), (1, odd)))
    balanced = triangle_blocks(((1, identity),) * 3)
    hostile_group = aut(hostile, unrestricted=True)
    balanced_group = aut(balanced, unrestricted=True)
    need(len(hostile_group) == 1, 'unrestricted nine factorial hostile rigidity')
    need(len(balanced_group) == 3, 'unrestricted nine factorial balanced symmetry')
    need(len(scc(hostile)) == 1, 'native hostile strong')
    need(error(hostile, hostile_group) == 18, 'native hostile error eighteen')
    need(sorted(map(len, refine(hostile))) == [3, 3, 3], 'native hostile cells three')
    need(sorted(row.bit_count() for row in hostile) == [3] * 3 + [4] * 3 + [5] * 3,
         'native hostile degree separation')

    old_hostile = (74, 84, 97, 22, 37, 11, 56)
    need(error(old_hostile) == 14, 'inherited seven hostile')
    scc_cases = []
    parts = [(), (0,), (0, 1), (2, 4, 1), old_hostile, hostile]
    for left, right in product(parts, repeat=2):
        rows = join((left, right))
        check_scc(rows)
        need(error(rows) == error(left) + error(right), 'ordered join additivity')
        scc_cases.append(len(rows))
    prime_bounds = {}
    for p in (5, 7):
        count = (factorial(p - 1) + (p - 1) ** 2) // p
        prime_bounds[p] = 2 * p * (p - 1) * (count - p + 1)
    need(prime_bounds == {5: 160, 7: 8568}, 'inherited prime family arithmetic')
    need(6 * 14 + 2 * 72 + 160 == 388, 'order ten lower bound arithmetic')
    report = {
        'schema': 'session50-cell-holonomy-v1',
        'small_labelled_counts': dict(small),
        'regular_census': regular,
        'three_triangle_family': dict(block_counts),
        'family_group_orders': dict(sorted(group_distribution.items())),
        'family_errors': dict(sorted(error_distribution.items())),
        'matching_only_family': {'labelled_presentations': 216, 'classes': 16,
            'class_group_orders': {1: 4, 3: 12}, 'strong_error_sum': 72},
        'derived_total_error_order_nine_lower_bound': 114,
        'inherited_prime_family_error_lower_bounds': prime_bounds,
        'derived_total_error_order_ten_lower_bound': 388,
        'native_hostile_nine': {'rows': hostile, 'stable_cells': refine(hostile),
            'literal_automorphisms': hostile_group, 'error': 18,
            'tested_unrestricted_permutations': factorial(9)},
        'balanced_control_nine': {'rows': balanced, 'automorphism_order': len(balanced_group),
            'tested_unrestricted_permutations': factorial(9)},
        'scc_join_controls': {'pairs': len(scc_cases), 'maximum_order': max(scc_cases)},
        'derived_strong_errors': {'through_order_six': 0, 'order_seven': 14, 'order_eight': 0},
        'derived_total_error_order_eight': 28,
        'checks': dict(sorted(COUNTS.items())),
    }
    payload = json.dumps(report, sort_keys=True, indent=2) + '\n'
    Path('receipt.json').write_text(payload)
    print(payload, end='')
    print('PASS all explicit finite controls; all-order claims depend on THEOREM.md proofs')


if __name__ == '__main__':
    main()
