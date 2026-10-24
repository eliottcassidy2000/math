#!/usr/bin/env python3
"""Independent finite certificate for the complete strong order-ten error.

Standard library only; no tournament atlas, target count, network, or external
source reads. Derives positive-error sectors from literal small presentations.
All checks are explicit exceptions and survive python -O. Outgoing bit rows.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb, factorial
import hashlib
import json

CHECKS = Counter()


def require(condition, name):
    CHECKS[name] += 1
    if not condition:
        raise RuntimeError(name)


def arc(rows, v, w):
    return (rows[v] >> w) & 1


def validate(rows):
    n = len(rows)
    require(all(0 <= row < 1 << n for row in rows), 'row bounds')
    require(all(not arc(rows, v, v) for v in range(n)), 'no loops')
    require(all(arc(rows, v, w) + arc(rows, w, v) == 1
                for v in range(n) for w in range(v)), 'tournament')


def move(rows, p):
    out = [0] * len(rows)
    for v, row in enumerate(rows):
        bits = row
        while bits:
            bit = bits & -bits
            out[p[v]] |= 1 << p[bit.bit_length() - 1]
            bits -= bit
    return tuple(out)


def strong(rows):
    n = len(rows)
    reach = [rows[v] | (1 << v) for v in range(n)]
    for k in range(n):
        for v in range(n):
            if (reach[v] >> k) & 1:
                reach[v] |= reach[k]
    return all(mask == (1 << n) - 1 for mask in reach)


def refine(rows):
    cells = (tuple(range(len(rows))),)
    while True:
        new = []
        for cell in cells:
            groups = {}
            for v in cell:
                key = tuple(sum(arc(rows, v, w) for w in block) for block in cells)
                groups.setdefault(key, []).append(v)
            new.extend(tuple(groups[key]) for key in sorted(groups))
        if len(new) == len(cells):
            return tuple(new)
        cells = tuple(new)


def automorphisms(rows):
    """Unrestricted literal arc backtracking, pruned only by vertex degree.

    Does not use stable cells, named blocks, or signed/matrix formulas.
    """
    n = len(rows)
    degrees = [r.bit_count() for r in rows]
    targets = [[w for w in range(n) if degrees[v] == degrees[w]] for v in range(n)]
    order = sorted(range(n), key=lambda v: (len(targets[v]), v))
    assignment = [-1] * n
    out = []
    def visit(k, used):
        if k == n:
            out.append(tuple(assignment))
            return
        v = order[k]
        for w in targets[v]:
            if (used >> w) & 1:
                continue
            if any(arc(rows, v, u) != arc(rows, w, assignment[u]) for u in order[:k]):
                continue
            assignment[v] = w
            visit(k + 1, used | (1 << w))
            assignment[v] = -1
    visit(0, 0)
    require(len(out) > 0, 'identity automorphism exists')
    return tuple(out)


def error(rows, actions=None):
    if actions is None:
        actions = automorphisms(rows)
    cells = refine(rows)
    size = {v: len(c) for c in cells for v in c}
    orbits = {tuple(sorted({p[v] for p in actions})) for v in range(len(rows))}
    return sum(size[o[0]] - len(o) for o in orbits)


def digest_rows(rows):
    return hashlib.sha256(json.dumps(sorted(rows), separators=(',', ':')).encode()).hexdigest()


def regular_tournaments(n):
    """Enumerate only the specified regular labelled sector, with exact pruning."""
    target = (n - 1) // 2
    pairs = tuple(combinations(range(n), 2))
    rows = [0] * n
    degree = [0] * n
    remaining = [n - 1] * n
    found = set()
    def visit(k):
        if k == len(pairs):
            require(all(x == target for x in degree), 'regular degree terminal')
            found.add(tuple(rows))
            return
        u, v = pairs[k]
        remaining[u] -= 1
        remaining[v] -= 1
        for a, b in ((u, v), (v, u)):
            degree[a] += 1
            rows[a] |= 1 << b
            if all(degree[w] <= target <= degree[w] + remaining[w] for w in (u, v)):
                visit(k + 1)
            rows[a] ^= 1 << b
            degree[a] -= 1
        remaining[u] += 1
        remaining[v] += 1
    visit(0)
    return found


def check_regular_seven():
    universe = regular_tournaments(7)
    require(len(universe) == 2640, 'regular seven labelled count')
    pending = set(universe)
    reps = []
    group = tuple(permutations(range(7)))
    while pending:
        rows = min(pending)
        orbit = {move(rows, p) for p in group}
        require(orbit <= pending, 'regular seven full permutation orbit disjoint')
        pending.difference_update(orbit)
        actions = automorphisms(rows)
        require(len(orbit) * len(actions) == factorial(7), 'regular seven orbit stabilizer')
        require(strong(rows), 'regular seven strong')
        e = error(rows, actions)
        reps.append({'rows': rows, 'automorphisms': len(actions), 'error': e,
                     'labelled_orbit': len(orbit)})
    require(sorted(r['automorphisms'] for r in reps) == [3, 7, 21], 'regular seven actions')
    require(sum(r['error'] for r in reps) == 14, 'regular seven total error')
    return reps, {'labelled': len(universe), 'classes': reps, 'sha256': digest_rows(universe)}


def carousel(n):
    return tuple(sum(1 << ((x + k) % n) for k in range(1, (n + 1) // 2)) for x in range(n))


def regular_matrices(n, d):
    row_options = tuple(sum(1 << x for x in c) for c in combinations(range(n), d))
    cols = [0] * n
    rows = []
    out = set()
    def visit(k):
        if k == n:
            require(all(c == d for c in cols), 'matrix column terminal')
            out.add(tuple(rows))
            return
        for row in row_options:
            for x in range(n):
                cols[x] += (row >> x) & 1
            if all(c <= d <= c + n - k - 1 for c in cols):
                rows.append(row)
                visit(k + 1)
                rows.pop()
            for x in range(n):
                cols[x] -= (row >> x) & 1
    visit(0)
    return out


def matrix_move(matrix, a, b):
    n = len(matrix)
    out = [0] * n
    for x in range(n):
        for y in range(n):
            if (matrix[x] >> y) & 1:
                out[(x + a) % n] |= 1 << ((y + b) % n)
    return tuple(out)


def two_cells(matrix):
    n = len(matrix)
    base = carousel(n)
    rows = list(base) + [r << n for r in base]
    for x in range(n):
        for y in range(n):
            u, v = x, n + y
            if (matrix[x] >> y) & 1:
                rows[u] |= 1 << v
            else:
                rows[v] |= 1 << u
    return tuple(rows)


def check_five_cells():
    require(len(automorphisms(carousel(5))) == 5, 'C5 full automorphism order')
    result = []
    all_errors = []
    for d in (1, 2):
        universe = regular_matrices(5, d)
        require(len(universe) == (120 if d == 1 else 2040), 'five crossing labelled count')
        pending = set(universe)
        orbit_histogram = Counter()
        error_histogram = Counter()
        witnesses = {}
        while pending:
            matrix = min(pending)
            orbit = {matrix_move(matrix, a, b) for a, b in product(range(5), repeat=2)}
            require(orbit <= pending, 'five crossing rotation orbit disjoint')
            pending.difference_update(orbit)
            rows = two_cells(matrix)
            validate(rows)
            require(strong(rows), 'five crossing strong')
            require(sorted(map(len, refine(rows))) == [5, 5], 'five actual stable cells')
            actions = automorphisms(rows)
            stabilizer = sum(matrix_move(matrix, a, b) == matrix for a, b in product(range(5), repeat=2))
            require(len(actions) == stabilizer, 'five full action matches independent stabilizer')
            require(len(orbit) * len(actions) == 25, 'five orbit stabilizer')
            e = error(rows, actions)
            require(e == (40 if len(actions) == 1 else 0), 'five actual orbit error')
            orbit_histogram[len(orbit)] += 1
            error_histogram[e] += 1
            all_errors.append(e)
            witnesses.setdefault(e, {'rows': rows, 'matrix': matrix, 'automorphisms': len(actions)})
        fixed_counts = Counter()
        for a, b in product(range(5), repeat=2):
            fixed = sum(matrix_move(matrix, a, b) == matrix for matrix in universe)
            expected = len(universe) if (a, b) == (0, 0) else (comb(5, d) if a and b else 0)
            require(fixed == expected, 'five independent fixed matrix counts')
            fixed_counts[fixed] += 1
        rigid_formula = (len(universe) - 4 * comb(5, d)) // 25
        require(error_histogram[40] == rigid_formula, 'five rigid formula')
        require(sum(f * multiplicity for f, multiplicity in fixed_counts.items()) == 25 * sum(orbit_histogram.values()),
                'five Burnside cross check')
        result.append({'crossing_degree': d, 'labelled_matrices': len(universe),
                       'orbit_histogram': sorted(orbit_histogram.items()),
                       'error_histogram': sorted(error_histogram.items()),
                       'fixed_count_histogram': sorted(fixed_counts.items()),
                       'witnesses': witnesses, 'sha256': digest_rows(universe)})
    require(factorial(5) ** 2 * (Fraction(1, 10) + Fraction(1, 6) * Fraction(1, 4)) == 2040,
            'bipartite cycle analytic matrix count')
    require(sum(all_errors) == 3360, 'five sector error')
    return {'sectors': result, 'total_error': sum(all_errors)}


def all_tournaments(n):
    pairs = tuple(combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        rows = [0] * n
        for k, (u, v) in enumerate(pairs):
            a, b = (u, v) if (bits >> k) & 1 else (v, u)
            rows[a] |= 1 << b
        yield tuple(rows)


def substitute(quotient, parts):
    offsets = []
    total = 0
    rows = []
    for part in parts:
        offsets.append(total)
        rows.extend(row << total for row in part)
        total += len(part)
    for i, part in enumerate(parts):
        for j, target in enumerate(parts):
            if arc(quotient, i, j):
                mask = ((1 << len(target)) - 1) << offsets[j]
                for v in range(offsets[i], offsets[i] + len(part)):
                    rows[v] |= mask
    return tuple(rows)


def check_seven_substitution(reps):
    universe = {q for q in all_tournaments(4) if strong(q)}
    group = tuple((0,) + p for p in permutations((1, 2, 3)))
    pending = set(universe)
    quotient_reps = []
    expansions = []
    while pending:
        q = min(pending)
        orbit = {move(q, p) for p in group}
        require(orbit <= pending, 'rooted quotient disjoint orbit')
        pending.difference_update(orbit)
        quotient_reps.append(q)
        for rep in reps:
            rows = substitute(q, [rep['rows'], (0,), (0,), (0,)])
            validate(rows)
            require(strong(rows), 'seven substitution strong')
            require(sorted(map(len, refine(rows))) == [1, 1, 1, 7], 'seven substitution actual cells')
            actions = automorphisms(rows)
            require(len(actions) == rep['automorphisms'], 'seven substitution full action')
            e = error(rows, actions)
            require(e == rep['error'], 'seven substitution inherited error')
            expansions.append({'rows': rows, 'error': e, 'automorphisms': len(actions), 'quotient': q})
    require(len(universe) == 24 and len(quotient_reps) == 4, 'strong rooted four quotient counts')
    require(sum(x['error'] for x in expansions) == 56, 'seven sector error')
    return {'labelled_quotients': len(universe), 'rooted_classes': quotient_reps,
            'expansions': expansions, 'total_error': sum(x['error'] for x in expansions)}


def triangle_presentation(majority, matches, attachment):
    rows = [0] * 10
    for i in range(3):
        for x in range(3):
            rows[3 * i + x] |= 1 << (3 * i + (x + 1) % 3)
    for (i, j), sigma in zip(combinations(range(3), 2), matches):
        for x in range(3):
            for y in range(3):
                forward = y != sigma[x] if arc(majority, i, j) else y == sigma[x]
                u, v = 3 * i + x, 3 * j + y
                a, b = (u, v) if forward else (v, u)
                rows[a] |= 1 << b
    for i in range(3):
        for x in range(3):
            a, b = (3 * i + x, 9) if attachment[i] else (9, 3 * i + x)
            rows[a] |= 1 << b
    return tuple(rows)


def triangle_block_group():
    result = []
    for pi in permutations(range(3)):
        for shifts in product(range(3), repeat=3):
            p = tuple(3 * pi[i] + (x + shifts[i]) % 3 for i in range(3) for x in range(3)) + (9,)
            result.append(p)
    return tuple(result)


def check_three_triangles():
    matchings = tuple(permutations(range(3)))
    group = triangle_block_group()
    universe = {}
    rejected = Counter()
    positive_controls = {}
    for q in all_tournaments(3):
        cyclic = strong(q)
        for t in product(range(2), repeat=3):
            for matches in product(matchings, repeat=3):
                rows = triangle_presentation(q, matches, t)
                unbalanced = 1
                for sigma in matches:
                    unbalanced *= 1 if (sigma[1] - sigma[0]) % 3 == 1 else -1
                is_strong = strong(rows)
                require(is_strong == (0 < sum(t) < 3), 'triangle singleton strong iff nonconstant')
                cells = sorted(map(len, refine(rows)))
                if not is_strong:
                    expected = [1, 9] if cyclic else [1, 3, 3, 3]
                    require(cells == expected, 'constant attachment hostile actual partition')
                    rejected[(cyclic, tuple(cells))] += 1
                    positive_controls.setdefault(('collapsed' if cyclic else 'nonstrong'), rows)
                    continue
                require(cells == [1, 3, 3, 3], 'triangle singleton actual stable partition')
                require(rows not in universe, 'triangle presentation injective')
                universe[rows] = (cyclic, unbalanced == -1)
    pending = set(universe)
    hist = Counter()
    witnesses = {}
    while pending:
        rows = min(pending)
        orbit = {move(rows, p) for p in group}
        require(orbit <= pending, 'triangle all block orbits disjoint complete')
        pending.difference_update(orbit)
        cyclic, unbalanced = universe[rows]
        require(all(universe[r] == (cyclic, unbalanced) for r in orbit), 'triangle orbit preserves sector and parity')
        actions = automorphisms(rows)
        expected = 1 if unbalanced else 3
        require(len(actions) == expected, 'triangle independent full automorphism')
        require(len(orbit) * len(actions) == len(group), 'triangle orbit stabilizer')
        e = error(rows, actions)
        require(e == (18 if unbalanced else 0), 'triangle actual error')
        hist[(cyclic, len(actions), e)] += 1
        witnesses.setdefault((cyclic, e), {'rows': rows, 'automorphisms': len(actions), 'orbit_size': len(orbit)})
    require(hist[(False, 1, 18)] == 24 and hist[(True, 1, 18)] == 8, 'triangle rigid class split')
    total = sum(key[2] * count for key, count in hist.items())
    require(total == 576, 'triangle total sector error')
    require(sum(hist.values()) == 128, 'triangle total strong class count')
    return {'strong_labelled_presentations': len(universe), 'strong_classes': sum(hist.values()),
            'class_histogram': [{'cyclic_majority': k[0], 'automorphisms': k[1], 'error': k[2], 'classes': v}
                                for k, v in sorted(hist.items())],
            'rejected_constant_presentations': [{'cyclic_majority': k[0], 'stable_sizes': k[1], 'presentations': v}
                                                for k, v in sorted(rejected.items())],
            'witnesses': [{'cyclic_majority': k[0], 'error': k[1], **v} for k, v in sorted(witnesses.items())],
            'hostile_rows': positive_controls, 'total_error': total, 'sha256': digest_rows(universe)}


def odd_partitions(n, floor=1):
    if not n:
        yield ()
    for x in range(floor, n + 1, 2):
        for rest in odd_partitions(n - x, x):
            yield (x,) + rest


def run():
    partitions = sorted(odd_partitions(10))
    require(len(partitions) == 10, 'all order ten odd partitions')
    expected = [(1,) * 10, (1,) * 7 + (3,), (1,) * 4 + (3, 3), (1, 3, 3, 3),
                (1,) * 5 + (5,), (1, 1, 3, 5), (5, 5), (1, 1, 1, 7), (3, 7), (1, 9)]
    require(partitions == sorted(expected), 'odd partition explicit classification')
    reps, seven = check_regular_seven()
    five = check_five_cells()
    substituted = check_seven_substitution(reps)
    triangles = check_three_triangles()
    es10 = five['total_error'] + substituted['total_error'] + triangles['total_error']
    require(es10 == 3992, 'complete strong error order ten')
    es7 = sum(r['error'] for r in reps)
    # Previously certified Session50 theorem input; not re-derived in this checker.
    es9_inherited = 716
    low_a = (1, 1, 1, 2)
    a_squared_3 = sum(low_a[i] * low_a[3 - i] for i in range(4))
    total = es10 + 2 * es9_inherited + a_squared_3 * es7
    require(total == 5508, 'SCC transfer order ten with inherited Es9')
    return {'status': 'PASS', 'scope': 'complete strong order-ten error from exhaustive positive sectors; total uses inherited Es9',
            'odd_partitions': partitions, 'regular_seven': seven, 'five_cells': five,
            'seven_cell_substitutions': substituted, 'three_triangle_cells': triangles,
            'strong_error_order_ten': es10, 'inherited_strong_error_order_nine': es9_inherited,
            'A_squared_coefficient_three': a_squared_3, 'total_error_order_ten': total,
            'check_counts': dict(sorted(CHECKS.items())),
            'boundary': {'full_order_ten_tournament_enumeration': False,
                         'target_tournament_count_input': False, 'external_runtime_inputs': False,
                         'full_target_automorphism_tests': 'only representatives of proved exceptional families'}}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
