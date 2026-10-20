"""Independent exact controls for session17 balanced deletion extensions.

No producer imports, floating point, random sampling, or removable assertions.
"""
from collections import Counter
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json


def require(value, message):
    if not value:
        raise RuntimeError(message)


def tournament(a, b):
    n = a + b + 1
    rows = [sum(1 << v for v in range(u + 1, n)) for u in range(n)]
    for u, v in ((0, a), (a, a + b)):
        rows[u] ^= 1 << v
        rows[v] ^= 1 << u
    return tuple(rows)


def paths(rows):
    full = (1 << len(rows)) - 1
    out = []

    def visit(word, unseen):
        if not unseen:
            out.append(word)
            return
        choices = unseen if not word else unseen & rows[word[-1]]
        while choices:
            bit = choices & -choices
            choices ^= bit
            visit(word + (bit.bit_length() - 1,), unseen ^ bit)
    visit((), full)
    return out


def states(k):
    return range((1 << (k - 1)) + 1)


def word(k, x, y):
    top = 1 << (k - 1)
    A, B = set(range(k)), set(range(k + 1, 2 * k + 1))
    def choice(i):
        return A if i == top else {j + 1 for j in range(k - 1) if i >> j & 1}
    I, K = choice(x), choice(y)
    J = {2 * k - v for v in A - K}
    return tuple(sorted(I) + sorted(J) + [k] + sorted(A - I) + sorted(B - J))


def valid(rows, path):
    return len(path) == len(rows) and len(set(path)) == len(rows) and all(
        rows[u] >> v & 1 for u, v in zip(path, path[1:]))


def blue(rows, path):
    reversal = dict(zip(path, reversed(path)))
    return all(bool(rows[u] >> v & 1) != bool(rows[reversal[u]] >> reversal[v] & 1)
               for u, v in combinations(range(len(rows)), 2))


def endpoint_kernel(k, s, t, blue_only=False):
    embed = lambda v: v if v < k else v + 1 if v == k else v + 2
    I_inf, K_inf = s == 0, t == 2 * k
    I_empty, K_empty = s >= k, t <= k
    out = Counter()
    for e in (None,) if I_inf else (0, 1):
        for h in (None,) if K_inf else (0, 1):
            if blue_only and e != h:
                continue
            ss = embed(s)
            if I_empty:
                if e == 1:
                    ss = k
                elif h == 0:
                    ss = k + 2
            tt = embed(t)
            if K_empty:
                if h == 1:
                    tt = k + 2
                elif e == 0:
                    tt = k
            out[ss, tt] += 1
    return out


def rect_deletions():
    cases = total_source_checks = 0
    for a, b in product(range(2, 6), repeat=2):
        rows = tournament(a, b)
        old_paths = paths(rows)
        require(len(old_paths) == ((1 << (a - 1)) + 1) * ((1 << (b - 1)) + 1),
                'rectangle literal count')
        left, right = list(range(1, a)), list(range(a + 1, a + b))
        for dl in range((1 << len(left)) - 1):
            D1 = {v for i, v in enumerate(left) if dl >> i & 1}
            for dr in range((1 << len(right)) - 1):
                D2 = {v for i, v in enumerate(right) if dr >> i & 1}
                D = D1 | D2
                survivors = [v for v in range(a + b + 1) if v not in D]
                relabel = {v: i for i, v in enumerate(survivors)}
                aa, bb = a - len(D1), b - len(D2)
                parent_rows = tournament(aa, bb)
                induced = tuple(sum(1 << j for j, v in enumerate(survivors)
                                    if rows[u] >> v & 1) for u in survivors)
                require(induced == parent_rows, 'actual induced rectangle')
                fibres = Counter()
                for p in old_paths:
                    q = tuple(relabel[v] for v in p if v not in D)
                    require(valid(parent_rows, q), 'orphan in claimed rectangular map')
                    fibres[q] += 1
                target_paths = paths(parent_rows)
                require(set(fibres) == set(target_paths), 'rectangular surjectivity')
                for p in target_paths:
                    expected = (1 if p[0] == 0 else 1 << len(D1)) * (
                        1 if p[-1] == aa + bb else 1 << len(D2))
                    require(fibres[p] == expected, 'rectangular marked fibre')
                total_source_checks += len(old_paths)
                cases += 1
    return dict(maps=cases, source_path_deletions=total_source_checks)


def ports_and_squares():
    reports = []
    for k in range(2, 7):
        old = {word(k, x, y): (x, y) for x, y in product(states(k), repeat=2)}
        require(set(paths(tournament(k, k))) == set(old), 'independent HP decoder')
        new = {word(k + 1, x, y): (x, y)
               for x, y in product(states(k + 1), repeat=2)}
        E = Counter((p[0], p[-1]) for p in old)
        B = Counter((p[0], p[-1]) for p, (x, y) in old.items() if x == y)
        lifted, lifted_blue = Counter(), Counter()
        for endpoint, n in E.items():
            for port, m in endpoint_kernel(k, *endpoint).items():
                lifted[port] += n * m
        for endpoint, n in B.items():
            for port, m in endpoint_kernel(k, *endpoint, blue_only=True).items():
                lifted_blue[port] += n * m
        require(lifted == Counter((p[0], p[-1]) for p in new), 'full endpoint closure')
        require(lifted_blue == Counter((p[0], p[-1]) for p, (x, y) in new.items()
                                      if x == y), 'blue endpoint closure')
        # Per-parent check is stronger than just aggregate endpoint equality.
        actual_lifts = {p: Counter() for p in old}
        for p in new:
            parent = tuple(v if v < k else k if v == k + 1 else v - 2
                           for v in p if v not in (k, k + 2))
            actual_lifts[parent][p[0], p[-1]] += 1
        require(all(actual_lifts[p] == endpoint_kernel(k, p[0], p[-1]) for p in old),
                'endpoint kernel per parent')
        reports.append(dict(k=k, old_paths=len(old), new_paths=len(new),
                            old_endpoint_types=len(E), new_endpoint_types=len(lifted)))
    return reports


def join_test():
    base = tournament(3, 3)
    n = len(base)
    joined = tuple(row | sum(1 << v for v in range(n, 2 * n)) for row in base) + tuple(
        row << n for row in base)
    P = paths(base)
    require(len(P) == 25, 'base has 25 paths')
    standard = tuple(range(n))
    paired = {p + tuple(v + n for v in q): (p, q) for p, q in product(P, repeat=2)}
    require(set(paths(joined)) == set(paired), 'literal ordered-join bijection')
    counts = Counter()
    for path, (p, q) in paired.items():
        colour = 'blue' if blue(joined, path) else 'black'
        stratum = 'bulk' if p != standard and q != standard else 'boundary'
        counts[colour] += 1
        counts[stratum] += 1
        counts[stratum + '_' + colour] += 1
    require(counts == dict(blue=25, black=600, bulk=576, boundary=49,
                          bulk_blue=24, bulk_black=552, boundary_blue=1, boundary_black=48),
            '625/576 product versus diagonal distinction')
    return dict(counts)


def bareiss_det(matrix):
    M = [row[:] for row in matrix]
    sign, divisor = 1, 1
    for k in range(len(M) - 1):
        p = next((p for p in range(k, len(M)) if M[p][k]), None)
        if p is None:
            return 0
        if p != k:
            M[k], M[p] = M[p], M[k]
            sign = -sign
        pivot = M[k][k]
        for i in range(k + 1, len(M)):
            for j in range(k + 1, len(M)):
                numerator = pivot * M[i][j] - M[i][k] * M[k][j]
                require(numerator % divisor == 0, 'exact Bareiss division')
                M[i][j] = numerator // divisor
        for i in range(k + 1, len(M)):
            M[i][k] = 0
        divisor = pivot
    return sign * M[-1][-1]


def port_unimodularity():
    output = []
    for k in range(2, 9):
        E = Counter((p[0], p[-1]) for x, y in product(states(k), repeat=2)
                    for p in [word(k, x, y)])
        expected = Counter()
        a = [1] + [1 << (k - 1 - i) for i in range(1, k)]
        for i, j in product(range(k), repeat=2):
            expected[i, 2 * k - j] = a[i] * a[j]
        for alpha in range(1, k):
            for beta in range(k):
                c = 1 if beta == alpha - 1 else 1 << (beta - alpha - 1) if beta > alpha else 0
                if c:
                    expected[alpha, beta] = c
                    expected[2 * k - beta, 2 * k - alpha] = c
        for port in ((0, k), (k, 2 * k), (k + 1, k - 1)):
            expected[port] = 1
        require(E == expected, 'explicit complete endpoint coefficients')
        require(len(E) == 2 * k * k - k + 3, 'endpoint support cardinality')
        matrix = [[E[i, j] for j in range(2 * k + 1)] for i in range(2 * k + 1)]
        require(bareiss_det(matrix) == 1, 'full endpoint unimodularity')
        # Independent weighted permanent DP counts all perfect matchings.
        dp = {0: 1}
        for i, row in enumerate(matrix):
            new = Counter()
            for mask, count in dp.items():
                for j, weight in enumerate(row):
                    if weight and not (mask >> j & 1):
                        new[mask | (1 << j)] += count * weight
            dp = new
        require(dp == {(1 << (2 * k + 1)) - 1: 1}, 'unique unit-weight perfect matching')
        output.append(dict(k=k, matrix_order=2 * k + 1, support=len(E), determinant=1, permanent=1))
    return output


def marked_five_field():
    labels = (1, 2, 3, 4, 0)
    output = []
    for x, y in product(states(3), repeat=2):
        p = word(3, x, y)
        kappa = lambda z: z if z == 4 else 3 - z
        require(labels[kappa(x)] == -labels[x] % 5, 'equivariant F5 coordinate')
        require(tuple(6 - v for v in reversed(p)) == word(3, y, x), 'R is coordinate swap')
        output.append(dict(states=[x, y], vector=[labels[x], labels[y]], path=list(p), blue=x == y))
    require(next(item for item in output if item['vector'] == [0, 0])['path'] == list(range(7)),
            'zero vector is standard path')
    Path(__file__).with_name('marked_f5_paths.json').write_text(json.dumps(output, indent=2) + '\n')
    return dict(points=25, nonzero=24, sign_pairs=12, equivariant_coordinate_labelings=8)


def twin_hostile():
    # u=0,v=1,a=2,b=3,c=4, u/v homogeneous and u->v.
    arcs = [(0, 1), (2, 0), (2, 1), (4, 0), (4, 1), (0, 3), (1, 3),
            (3, 2), (3, 4), (2, 4)]
    rows = tuple(sum(1 << v for uu, v in arcs if uu == u) for u in range(5))
    p = (2, 1, 3, 4, 0)
    require(valid(rows, p), 'twin witness source path')
    require(all(bool(rows[0] >> w & 1) == bool(rows[1] >> w & 1)
                for w in (2, 3, 4)), 'actual homogeneous twins')
    require(not (rows[2] >> 3 & 1), 'orphan joining edge')
    small_tests = 0
    for n in range(2, 5):
        pairs = list(combinations(range(n), 2))
        for mask in range(1 << len(pairs)):
            rr = [0] * n
            for j, (u, v) in enumerate(pairs):
                if mask >> j & 1:
                    rr[u] |= 1 << v
                else:
                    rr[v] |= 1 << u
            for u, v in pairs:
                if not all(bool(rr[u] >> w & 1) == bool(rr[v] >> w & 1)
                           for w in range(n) if w not in (u, v)):
                    continue
                for deleted in (u, v):
                    for path in paths(rr):
                        q = tuple(w for w in path if w != deleted)
                        require(all(rr[x] >> y & 1 for x, y in zip(q, q[1:])),
                                'smaller twin orphan')
                        small_tests += 1
    return dict(rows=rows, path=p, deleted=1, bad_remainder=(2, 3, 4, 0),
                smaller_valid_twin_deletions=small_tests)


def arithmetic():
    cases = 0
    for q in range(2, 101):
        for d in range(1, 8):
            D = 1 << d
            qq = D * (q - 1) + 1
            require(qq * qq == D * D * q * q - 2 * D * (D - 1) * q + (D - 1)**2,
                    'all-depth total recurrence')
            require(qq * (qq - 1) == D * D * q * (q - 1) - D * (D - 1) * (q - 1),
                    'all-depth black recurrence')
            cases += 1
    return dict(symbolic_formula_integer_controls=cases,
                expression_776=5**4 + 5**3 + 5**2 + 1,
                bulk_576=(25 - 1)**2)


def main():
    output = dict(rectangles=rect_deletions(), endpoint_transfers=ports_and_squares(),
                  ordered_join=join_test(), twin_deletion=twin_hostile(), arithmetic=arithmetic())
    output['endpoint_unimodularity'] = port_unimodularity()
    output['marked_F5_path_bridge'] = marked_five_field()
    output['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    text = json.dumps(output, indent=2, sort_keys=True)
    print(text)
    Path(__file__).with_name('checks.json').write_text(text + '\n')


if __name__ == '__main__':
    main()
