#!/usr/bin/env python3
"""Independent literal controls. All failures remain active under Python -O."""
from itertools import combinations, product
from math import comb, gcd
from collections import Counter
from core import endpoints, h, covers, ext, induced, ordinal, strong, inflate


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def tournament(n, code):
    rows = [0] * n
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        if (code >> bit) & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return tuple(rows)


def legal(rows):
    n = len(rows)
    return all(not (rows[i] >> i & 1) for i in range(n)) and all(
        ((rows[i] >> j & 1) + (rows[j] >> i & 1)) == 1
        for i, j in combinations(range(n), 2))


def H(rows):
    return h(rows) if rows else 1


def scc_count(rows):
    n = len(rows)
    reach = []
    for i in range(n):
        seen = 1 << i
        while True:
            nxt = seen
            for v in range(n):
                if seen >> v & 1:
                    nxt |= rows[v]
            if nxt == seen:
                break
            seen = nxt
        reach.append(seen)
    classes = {tuple(j for j in range(n) if (reach[i] >> j & 1) and (reach[j] >> i & 1))
               for i in range(n)}
    return len(classes)


def carrier_test(rows, all_children=True):
    require(legal(rows), ('not tournament', rows))
    n = len(rows)
    old = H(rows)
    q = covers(rows)
    c = [[0] * n for _ in range(n)]
    for u, v in combinations(range(n), 2):
        k = q[u][v] + q[v][u]
        require(k % 2 == 0, ('odd K', rows, u, v, k))
        floor = H(induced(rows, [i for i in range(n) if i not in (u, v)]))
        require(k >= 2 * floor, ('deletion floor', rows, u, v, k, floor))
        a, b = (u, v) if rows[u] >> v & 1 else (v, u)
        c[a][b] = k // 2
    masks = range(1 << n) if all_children else sorted(set([0, (1 << n) - 1] +
                 [1 << i for i in range(n)] + [(1 << i) | (1 << j) for i, j in combinations(range(n), 2)]))
    equals = 0
    for mask in masks:
        cut = sum(c[u][v] for u in range(n) for v in range(n)
                  if mask >> u & 1 and not mask >> v & 1)
        got = H(ext(rows, mask))
        require(got == old + 2 * cut, ('canonical carrier', rows, mask, got, cut))
        no_old_outcut = all(not (rows[u] >> v & 1) for u in range(n) for v in range(n)
                            if mask >> u & 1 and not mask >> v & 1)
        require((got == old) == no_old_outcut, ('equality geometry', rows, mask))
        equals += got == old
    if all_children:
        require(equals == scc_count(rows) + 1, ('SCC count', rows, equals))
    return old, c


def observer(rows, port):
    ha_rows = ext(rows, 1 << port)
    haa_rows = ext(ha_rows, 1 << len(rows))
    a, b, c = H(rows), H(ha_rows), H(haa_rows)
    require((b - a) % 2 == 0 and (c - b) % 2 == 0, 'observer parity')
    return a, (b - a) // 2, (c - b) // 2


def gate(x, y):
    return ext((8, 1, 3, 6), 2 | x | (y << 2))


def wrapper(rows):
    n = len(rows)
    return tuple(r | (1 << n) for r in rows) + (1 << (n + 1), (1 << n) - 1)


def F2(rows):
    if not rows: return 1
    return 2 * H(rows) + sum(map(sum, covers(rows)))


def valuation3(x):
    ans = 0
    while x % 3 == 0:
        x //= 3
        ans += 1
    return ans


def mobius(values):
    ans = list(values)
    n = (len(values) - 1).bit_length()
    require(1 << n == len(values), 'cube size')
    for i in range(n):
        for mask in range(1 << n):
            if mask >> i & 1:
                ans[mask] -= ans[mask ^ (1 << i)]
    return ans


def linear_forest(n, edges):
    parent = list(range(n))
    degree = [0] * n
    def find(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for u, v in edges:
        degree[u] += 1
        degree[v] += 1
        if degree[u] > 2 or degree[v] > 2:
            return False
        a, b = find(u), find(v)
        if a == b:
            return False
        parent[a] = b
    return True


def main():
    counts = []
    for n in range(1, 5):
        count = 0
        for code in range(1 << comb(n, 2)):
            carrier_test(tournament(n, code))
            count += 1
        counts.append((n, count))
    for n, code in [(5, 8), (5, 759), (5, 341)]:
        carrier_test(tournament(n, code))
    U = (112, 97, 115, 103, 74, 16, 32)
    R4 = (996, 993, 994, 999, 975, 912, 800, 576, 128, 256)
    carrier_test(U)
    carrier_test(R4, all_children=False)
    print('Canonical carrier, deletion floors, SCC equality:', counts,
          '+ 3 order-five cases, U, and R4 controls')

    for x, y in product(range(2), repeat=2):
        g = gate(x, y)
        t = y * (1 - x)
        require(legal(g) and strong(g), ('gate connectivity', x, y))
        require(H(g) == 9 + 2 * t, ('gate H', x, y))
        require(F2(g) == 84 + 16 * t, ('gate F2', x, y))
        expected = (9 + 2 * y - 2 * x * y, 4 + 3 * x + 5 * y, 7 + 2 * x + 4 * y)
        require(observer(g, 4) == expected, ('gate observer', x, y, observer(g, 4)))
        print('Gate', (x, y), 'rows', g, 'H,F2', H(g), F2(g), 'observer', expected)

    for s in [1, 2]:
        n = s % 2
        by_count = Counter()
        by_valuation = Counter()
        values = []
        for mask in range(1 << (2 * s)):
            bits = [(mask >> i) & 1 for i in range(2 * s)]
            gs = [gate(bits[2 * i], bits[2 * i + 1]) for i in range(s)]
            k = sum(bits[2 * i + 1] * (1 - bits[2 * i]) for i in range(s))
            y = ordinal(*gs, tuple(0 for _ in range(n)))
            z = wrapper(y)
            require(legal(z) and strong(z), ('wrapper connectivity', s, mask))
            got = H(z)
            expected = 9 ** (s - k) * 11 ** k + 2 ** n * 84 ** (s - k) * 100 ** k
            require(got == expected, ('literal wrapper H', s, mask, got, expected))
            require(valuation3(got) == s - k, ('valuation', s, mask))
            v = observer(z, len(y))
            require(valuation3(gcd(*v)) == s - k, ('full observer valuation', s, mask, v))
            by_count[got] += 1
            by_valuation[s - k] += 1
            values.append(got)
        require(by_valuation == Counter({r: comb(s, r) * 3 ** r for r in range(s + 1)}),
                ('valuation distribution', s, by_valuation))
        require(sum(v % 3 != 0 for v in values) == 1, ('sharp support', s))
        print('Strong composition', s, 'counts', dict(sorted(by_count.items())),
              'valuations', dict(sorted(by_valuation.items())))
        if s == 2:
            coefficients = mobius(values)
            require(coefficients[-1] == 260, ('degree-four sharpness', coefficients[-1]))

    quotient = (8, 1, 3, 6)
    base_v = observer(quotient, 3)
    for bits in product(range(2), repeat=4):
        g1, g2 = gate(*bits[:2]), gate(*bits[2:])
        z = inflate(quotient, [g1, g2, (0,), (0,)])
        factor = 4 * bits[1] * (1 - bits[0]) * bits[3] * (1 - bits[2])
        require(legal(z) and strong(z), ('modular substitution connectivity', bits))
        v = observer(z, 11)
        require(all((a - factor * b) % 3 == 0 for a, b in zip(v, base_v)),
                ('modular substitution observer', bits, v, base_v, factor))
    print('Uniform substitution: all 16 strong two-gate marked observers match modulo3')

    edges = list(combinations(range(4), 2))
    values = [H(tournament(4, code)) for code in range(64)]
    coefficients = mobius(values)
    checked_zero = 0
    for mask, coefficient in enumerate(coefficients):
        subset = [edge for i, edge in enumerate(edges) if mask >> i & 1]
        if not linear_forest(4, subset):
            require(coefficient == 0, ('forbidden monomial', subset, coefficient))
            checked_zero += 1
    print('Full orientation cube: all', checked_zero, 'non-linear-forest coefficients vanish')

    # Independent finite minimum-support controls over composite groups.
    for modulus in [3, 5, 9, 15]:
        for c in range(modulus):
            nonzero = sum((v - c) % modulus != 0 for v in values)
            require(nonzero == 0 or nonzero >= 2 ** (6 - 3), ('minimum support', modulus, c))
    print('PASS: all literal, structural, polynomial, and valuation checks')

if __name__ == '__main__':
    main()
