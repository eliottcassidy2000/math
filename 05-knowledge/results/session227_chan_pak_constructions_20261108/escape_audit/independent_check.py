#!/usr/bin/env python3
"""Independent exact audit of Session227 boundary escape.

Standard library only; no producer imports, downloaded bank, randomized sample,
or disabled-under-optimization assertion.  Tournament rows encode out-neighbors.
The main engine deletes the FIRST vertex recursively.  A separate literal
permutation/split engine is the small-order cross-check.
"""
from functools import lru_cache
from itertools import combinations, permutations, product
import json


def need(value, label):
    if not value:
        raise RuntimeError(label)


def tournaments(n):
    edges = tuple(combinations(range(n), 2))
    for code in range(1 << len(edges)):
        rows = [0] * n
        for k, (a, b) in enumerate(edges):
            u, v = (a, b) if code & (1 << k) else (b, a)
            rows[u] |= 1 << v
        yield tuple(rows)


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask -= bit


def subsets(mask):
    sub = mask
    while True:
        yield sub
        if not sub:
            break
        sub = (sub - 1) & mask


class Counts:
    def __init__(self, rows):
        self.rows = tuple(rows)
        self.n = len(rows)
        self.full = (1 << self.n) - 1
        self.reverse = tuple(sum(1 << w for w in range(self.n)
                                 if rows[w] & (1 << v))
                             for v in range(self.n))
        need(all(not rows[i] & (1 << i) and rows[i] & ~self.full == 0
                 for i in range(self.n)), 'row bounds')
        need(all(bool(rows[i] & (1 << j)) != bool(rows[j] & (1 << i))
                 for i in range(self.n) for j in range(i)), 'tournament')

    @lru_cache(None)
    def first(self, mask, v):
        if not mask & (1 << v):
            return 0
        rest = mask ^ (1 << v)
        if not rest:
            return 1
        return sum(self.first(rest, w) for w in bits(rest & self.rows[v]))

    @lru_cache(None)
    def last(self, mask, v):
        if not mask & (1 << v):
            return 0
        rest = mask ^ (1 << v)
        if not rest:
            return 1
        return sum(self.last(rest, w) for w in bits(rest & self.reverse[v]))

    @lru_cache(None)
    def h(self, mask):
        return sum(self.first(mask, v) for v in bits(mask)) if mask else 1

    def state(self, r):
        return (self.h(self.full), self.last(self.full, r),
                sum(self.first(m, r) * self.h(self.full ^ m)
                    for m in range(self.full + 1)),
                self.h(self.full ^ (1 << r)))

    def cover(self, mask, w, v):
        # Both bins are nonempty: the first ends at w, second starts at v.
        return sum(self.last(m, w) * self.first(mask ^ m, v)
                   for m in subsets(mask))


def brute(rows, r):
    def path(p):
        return all(rows[a] & (1 << b) for a, b in zip(p, p[1:]))
    n = len(rows)
    D = E = U = s = 0
    for p in permutations(range(n)):
        if path(p):
            D += 1
            E += p[-1] == r
            s += p[0] == r
        for cut in range(n):
            U += p[cut] == r and path(p[:cut]) and path(p[cut:])
    d = sum(path(p) for p in permutations(v for v in range(n) if v != r))
    return (D, E, U, d), s


def graft(a, ar, b, br):
    n = len(a)
    result = list(a) + [row << n for row in b]
    for u in range(n):
        for v in range(len(b)):
            if u == ar:
                result[n + v] |= 1 << u
            else:
                result[u] |= 1 << (n + v)
    return tuple(result), n + br


def op(x, letter):
    D, E, U, d = x
    if letter == 'O':
        return D, 0, 2 * U, d
    if letter == 'Z':
        return D - E + U, d, 2 * U - d, d
    if letter == 'J':
        return D - E + U, D - E, D + U, D
    raise ValueError(letter)


def literal_op(a, root, letter):
    if letter == 'O':
        n = len(a)
        return tuple(row | (1 << n) for row in a) + (0,), root
    rows, newroot = graft(a, root, (0,), 0)
    return rows, root if letter == 'Z' else newroot


def f(x):
    D, E, U, d = x
    return U - D - d + E


def boundary(x):
    D, E, U, d = x
    return E == 0 or E == d or f(x) == 0


def strong(rows):
    full = (1 << len(rows)) - 1
    for v in range(len(rows)):
        seen = 1 << v
        while True:
            new = seen
            for u in bits(seen):
                new |= rows[u]
            if new == seen:
                break
            seen = new
        if seen != full:
            return False
    return True


def multiply(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def det(a):
    n = len(a)
    total = 0
    for p in permutations(range(n)):
        sign = (-1) ** sum(p[i] > p[j] for i in range(n) for j in range(i + 1, n))
        term = sign
        for i in range(n):
            term *= a[i][p[i]]
        total += term
    return total


def mv(matrix, x):
    return tuple(sum(a * b for a, b in zip(row, x)) for row in matrix)


def main():
    sources = []
    root_checks = brute_checks = pair_checks = 0
    root_degrees = {}
    for n in range(1, 6):
        for rows in tournaments(n):
            c = Counts(rows)
            for r in range(n):
                x = c.state(r)
                D, E, U, d = x
                out = list(bits(rows[r]))
                root_checks += 1
                root_degrees[str(len(out))] = root_degrees.get(str(len(out)), 0) + 1
                need(D >= 1 and 0 <= E <= min(D, d) and U >= D, 'physical cone')
                need((U == d) == (len(out) == 0), 'root singleton-suffix iff sink')
                need((D == E) == (len(out) == 0), 'all paths end at root iff sink')
                need((f(x) == 0) == (len(out) <= 1), 'f equality iff degree <= 1')
                if n <= 3:
                    need(boundary(x), 'no physical escape before order four')
                mask = c.full ^ (1 << r)
                exact_f = sum(c.cover(mask, w, v) for v in out for w in out if v != w)
                need(exact_f == f(x), 'unmatched split interpretation of f')
                lower = 0
                for v, w in combinations(out, 2):
                    remainder = mask ^ (1 << v) ^ (1 << w)
                    pair_bound = 2 * c.h(remainder)
                    need(c.cover(mask, w, v) + c.cover(mask, v, w) >= pair_bound,
                         'pairwise telescoping lower bound')
                    lower += pair_bound
                    pair_checks += 1
                need(f(x) >= lower >= len(out) * (len(out) - 1), 'inherited f bound')
                direct, start = brute(rows, r)
                need(x == direct and c.first(c.full, r) == start, 'permutation cross-check')
                brute_checks += 1
            if n <= 4:
                sources.append((rows, c.state(0), c.first(c.full, 0)))
            # Memoized instance methods retain keys until cleared.
            c.first.cache_clear(); c.last.cache_clear(); c.h.cache_clear()
    need(root_checks == 5405 and len(sources) == 75, 'explicit bank size')

    # The elementary telescoping inequality requires no tournament assumptions.
    # Exhaust all binary a_i,b_i sequences through length seven.
    binary_cases = 0
    for k in range(8):
        for middle in product(range(4), repeat=k):
            seq = [(1, 0)] + [(j >> 1, j & 1) for j in middle] + [(0, 1)]
            total = 0
            for (a, b), (ap, bp) in zip(seq, seq[1:]):
                ci = a * bp + (1 - b) * (1 - ap)
                need(ci - (a - b - ap + bp) == (1 - a) * (1 - bp) + b * ap,
                     'telescoping identity')
                total += ci
            need(total >= 2, 'telescoped lower bound including empty remainder')
            binary_cases += 1

    pairs = escape_pairs = nontrivial_pairs = singleton_pairs = generator_checks = 0
    degree_cases = {}
    for a, ax, astart in sources:
        for letter in 'OZJ':
            rows, r = literal_op(a, 0, letter)
            c = Counts(rows)
            observed = c.state(r)
            need(observed == op(ax, letter), 'literal generator formula')
            need(boundary(observed), 'one-letter boundary containment')
            if len(a) > 1:
                need(c.first(c.full, r) == astart, 'O/Z/J preserve root starts on nonsingletons')
            generator_checks += 1
            c.first.cache_clear(); c.last.cache_clear(); c.h.cache_clear()
        for b, bx, bstart in sources:
            rows, r = graft(a, 0, b, 0)
            c = Counts(rows)
            x = c.state(r)
            D, E, U, d = x
            need(rows[r].bit_count() == 1 + b[0].bit_count(), 'output root degree')
            if len(a) > 1:
                need(c.first(c.full, r) == astart * bx[2], 'root-start product, nonsingleton old graph')
            else:
                need(c.first(c.full, r) == bstart, 'singleton old graph root-start exception')
            if len(b) == 1:
                need(x == op(ax, 'J') and boundary(x), 'singleton gadget J exception')
                singleton_pairs += 1
            else:
                nontrivial_pairs += 1
                expected_escape = bool(a[0]) and bool(b[0])
                need(not boundary(x) == expected_escape, 'all-order escape iff, finite control')
                need(d - E >= ax[3] * bx[3] > 0, 'deleted endpoint injection')
                need((E > 0) == bool(a[0]), 'E positivity iff old root nonsink')
                need((f(x) > 0) == bool(b[0]), 'f positivity iff gadget root nonsink')
                degree_key = f'{bool(a[0])},{bool(b[0])}'
                degree_cases[degree_key] = degree_cases.get(degree_key, 0) + 1
                escape_pairs += int(expected_escape)
            pairs += 1
            c.first.cache_clear(); c.last.cache_clear(); c.h.cache_clear()
    need(pairs == 5625 and singleton_pairs == 75 and nontrivial_pairs == 5550,
         'graft universe cardinality')

    # Universal last-letter obstruction is algebraic even on nonphysical inputs.
    arbitrary_controls = 0
    for x in product(range(-2, 3), repeat=4):
        for letter in 'OZJ':
            need(boundary(op(x, letter)), 'boundary containment without physical assumption')
            arbitrary_controls += 1

    L = ((1, -1, 3, -2), (0, 0, 1, -1),
         (2, -1, 4, -1), (1, -1, 1, 0))
    linear_checks = 0
    for a, ax, _ in sources:
        rows, r = graft(a, 0, (2, 0), 0)
        c = Counts(rows)
        need(c.state(r) == mv(L, ax), 'two-vertex matrix')
        zzrows, zzroot = literal_op(a, 0, 'Z')
        zzrows, zzroot = literal_op(zzrows, zzroot, 'Z')
        need(zzrows == rows and r == len(a), 'L is two Z births with first-newcomer rerooting')
        linear_checks += 1
        sink_rows, sink_r = graft(a, 0, (2, 0), 1)
        sink_c = Counts(sink_rows)
        need(sink_c.state(sink_r) == op(op(ax, 'Z'), 'J'), 'sink-rooted edge is ZJ')
        need(not boundary(c.state(r)) == bool(a[0]), 'edge escape iff old root nonsink')
        c.first.cache_clear(); c.last.cache_clear(); c.h.cache_clear()
    t = (1, 0, 2, 1)
    rows, r = graft((2, 0), 0, (2, 0), 0)
    lx, _ = brute(rows, r)
    need(lx == mv(L, t) == (5, 1, 9, 3), 'minimal escaping witness')
    need(not boundary(lx), 'strict escape witness')
    closed_rows, closed_root = literal_op(rows, r, 'J')
    closed, _ = brute(closed_rows, closed_root)
    need(closed == (13, 4, 14, 5), 'literal H13 after J')
    scalar_hostile = op(op(t, 'J'), 'J')
    need(scalar_hostile == (5, 2, 6, 3), 'same scalar, different reachable state')
    need(scalar_hostile[0] == lx[0] and scalar_hostile != lx,
         'state separation does not imply scalar separation')
    need(op(op(t, 'O'), 'J') == (5, 1, 5, 1), 'producer scalar-hostile OJt')
    need(strong(rows) and strong(closed_rows), 'strong literal H5/H13 witnesses')
    M = ((0, 1, -1), (3, 4, 1), (1, 1, 1))
    Minv = ((-3, 2, -5), (2, -1, 3), (1, -1, 3))
    identity = tuple(tuple(int(i == j) for j in range(3)) for i in range(3))
    need(multiply(M, Minv) == multiply(Minv, M) == identity, 'two-sided integral inverse')
    need(det(M) == -1, 'M determinant')
    need(det(tuple(row[:3] for row in L[:3])) == -1, 'L nonzero third-order minor')
    need(mv(L, (-2, -1, 1, 1)) == (0, 0, 0, 0), 'L kernel')
    need(tuple(L[0][j] - 2 * L[1][j] - L[3][j] for j in range(4)) == (0, 0, 0, 0),
         'L image in plane')
    M2, M3 = multiply(M, M), multiply(multiply(M, M), M)
    need(all(M3[i][j] - 5 * M2[i][j] + M[i][j] + identity[i][j] == 0
             for i in range(3) for j in range(3)), 'cubic identity')
    plane_tests = 0
    for e0, u0, d0 in product(range(-2, 3), repeat=3):
        x0 = (2 * e0 + d0, e0, u0, d0)
        y0 = mv(L, x0)
        need((y0[1], y0[2], y0[3]) == mv(M, (e0, u0, d0)), 'restricted plane matrix')
        plane_tests += 1
    family, x, family_rows, family_root = [], t, (2, 0), 0
    for m in range(7):
        closedx = op(x, 'J')
        entry = {'m': m, 'h': x[0], 'closed_h': closedx[0]}
        if m <= 3:
            cc = Counts(family_rows)
            need(cc.state(family_root) == x, 'literal L iterate')
            need(cc.first(cc.full, family_root) == 2 ** m, 'literal root-start power')
            if m:
                need(strong(family_rows) and not boundary(x), 'strong escaping L iterate')
            cr, rr = literal_op(family_rows, family_root, 'J')
            need(Counts(cr).state(rr) == closedx and strong(cr), 'literal strong closing J')
            entry['literal_vertices'] = len(family_rows)
            family_rows, family_root = graft(family_rows, family_root, (2, 0), 0)
            cc.first.cache_clear(); cc.last.cache_clear(); cc.h.cache_clear()
        family.append(entry)
        x = mv(L, x)
    need([v['h'] for v in family] == [1, 5, 25, 119, 565, 2681, 12721], 'h sequence')
    need([v['closed_h'] for v in family] == [3, 13, 61, 289, 1371, 6505, 30865], 'closed sequence')
    for i in range(4):
        for key in ('h', 'closed_h'):
            need(family[i + 3][key] == 5 * family[i + 2][key] - family[i + 1][key] - family[i][key],
                 'scalar cubic recurrence')
    print(json.dumps({
        'status': 'PASS',
        'universal_bound_controls': {
            'universe': 'All labelled tournaments of orders 1..5, every root',
            'rooted_instances': root_checks,
            'permutation_profile_and_start_crosschecks': brute_checks,
            'outneighbor_pair_bounds': pair_checks,
            'root_outdegree_histogram': root_degrees,
            'binary_telescoping_cases_through_length7': binary_cases},
        'graft_controls': {
            'universe': 'All labelled A and B of orders 1..4, each root fixed at label 0',
            'source_count': len(sources), 'gadget_count': len(sources),
            'pairs': pairs, 'singleton_gadget_pairs': singleton_pairs,
            'nonsingleton_gadget_pairs': nontrivial_pairs,
            'escaping_pairs': escape_pairs,
            'old_and_gadget_nonsink_cases': degree_cases,
            'root_start_product_checks_including_singleton_exception': pairs},
        'literal_generator_controls': generator_checks,
        'arbitrary_integer_generator_controls': arbitrary_controls,
        'two_vertex_matrix': L, 'two_vertex_matrix_controls': linear_checks,
        'minimal_escape': {'source': t, 'gadget_rows': [2, 0], 'gadget_root': 0,
                           'output_rows': rows, 'output_root': r, 'state': lx,
                           'f': f(lx), 'J_output_rows': closed_rows,
                           'J_state': closed},
        'scalar_hostile_JJt': scalar_hostile,
        'restricted_plane_integer_controls': plane_tests,
        'restricted_matrix': M, 'restricted_inverse': Minv,
        'family': family,
        'scope': 'Exact finite controls, not an all-order proof or scalar-completeness claim; inherited universal bound is separately credited.'
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
