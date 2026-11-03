#!/usr/bin/env python3
"""Independent finite control for joint tournament/SCC/converse count formulas.
Universe: all ordinary tournament classes, orders 0..6, by full augmentation.
No imported repository code/catalogue; only Python standard library.
"""
from collections import Counter
from itertools import permutations, product
from math import factorial
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def encode(rows, order):
    return sum(((rows[u] >> v) & 1) << k
               for k, (u, v) in enumerate((order[i], order[j])
               for i in range(len(order)) for j in range(i + 1, len(order))))


def decode(code, n):
    rows = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if (code >> k) & 1:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
            k += 1
    return tuple(rows)


def orders(rows):
    cells = {}
    for i, row in enumerate(rows):
        cells.setdefault(row.bit_count(), []).append(i)
    for blocks in product(*(permutations(cells[d]) for d in sorted(cells))):
        yield tuple(x for block in blocks for x in block)


def canonical(rows):
    return min(encode(rows, p) for p in orders(rows))


def hdp(rows):
    n = len(rows)
    if not n:
        return 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v]:
                avail = rows[v] & ~mask
                while avail:
                    bit = avail & -avail
                    w = bit.bit_length() - 1
                    dp[mask | bit][w] += dp[mask][v]
                    avail -= bit
    return sum(dp[-1])


def isstrong(rows):
    n = len(rows)
    if not n:
        return False
    reach = [r | (1 << i) for i, r in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] & (1 << k):
                reach[i] |= reach[k]
    return all(r == (1 << n) - 1 for r in reach)


def converse(rows):
    allv = (1 << len(rows)) - 1
    return tuple(allv ^ r ^ (1 << i) for i, r in enumerate(rows))


def finite_inventory(bound=6):
    classes = [{0}]
    data = []
    for n in range(bound + 1):
        if n:
            reps = set()
            for code in classes[n - 1]:
                old = decode(code, n - 1)
                for cut in range(1 << (n - 1)):
                    rows = tuple(r | ((1 << (n - 1)) if (cut >> i) & 1 else 0)
                                 for i, r in enumerate(old))
                    rows += (((1 << (n - 1)) - 1) ^ cut,)
                    reps.add(canonical(rows))
            classes.append(reps)
        layer = []
        for code in sorted(classes[n]):
            rows = decode(code, n)
            h = hdp(rows)
            literal = sum(all(rows[p[i]] & (1 << p[i + 1])
                              for i in range(n - 1))
                          for p in permutations(range(n)))
            check(h == literal, ('Hamilton DP/permutation mismatch', n, code))
            aut = sum(encode(rows, p) == encode(rows, tuple(range(n)))
                      for p in orders(rows))
            check(aut % 2 == 1 and h % aut == 0, ('orbit failure', n, code))
            # An acyclic tournament has no directed triangle, hence is transitive.
            triangles = sum(bool(rows[i] & (1 << j) and rows[j] & (1 << k)
                                 and rows[k] & (1 << i))
                            for i in range(n) for j in range(n) for k in range(n)) // 3
            check((h == 1) == (triangles == 0), ('H=1 boundary', n, code))
            layer.append(dict(code=code, h=h, aut=aut, strong=isstrong(rows),
                              self_converse=canonical(converse(rows)) == code))
        check(sum(factorial(n) // r['aut'] for r in layer) == 2 ** (n * (n - 1) // 2),
              ('labelled mass mismatch', n))
        data.append(layer)
    check([len(x) for x in data] == [1, 1, 1, 2, 4, 12, 56], 'class total')
    return data


def run():
    data = finite_inventory()
    ordinary = [Counter(r['h'] for r in layer) for layer in data]
    strong = [Counter(r['h'] for r in layer if r['strong']) for layer in data]
    fixed = [Counter(r['h'] for r in layer if r['self_converse']) for layer in data]
    center = [Counter(r['h'] for r in layer if r['strong'] and r['self_converse'])
              for layer in data]
    center[0] = Counter({1: 1})
    for n in range(1, len(data)):
        expected = Counter()
        for k in range(1, n + 1):
            for h, a in strong[k].items():
                for g, b in ordinary[n - k].items():
                    expected[h * g] += a * b
        check(expected == ordinary[n], ('SCC joint recurrence', n))
    for n in range(len(data)):
        expected = Counter()
        for m in range(n // 2 + 1):
            for u, a in ordinary[m].items():
                for v, b in center[n - 2 * m].items():
                    expected[u * u * v] += a * b
        check(expected == fixed[n], ('square-central recurrence', n))
    # Busch's f(7)=25 and monotonicity make the following h<=21 inventory global.
    low_strong = {}
    for h in range(3, 22, 2):
        low_strong[h] = [[n, strong[n][h], center[n][h]]
                         for n in range(3, len(data)) if strong[n][h]]
    expected_low = {3:[[3,1,1]],5:[[4,1,1]],7:[],9:[[5,2,2]],11:[[5,1,1]],
                    13:[[5,1,1]],15:[[5,2,2],[6,2,0]],17:[[6,2,2]],
                    19:[[6,2,0]],21:[]}
    check(low_strong == expected_low, 'complete low strong inventory')
    thresholds = {1:[0,1],3:[None,3],5:[4,None],7:[None,None],9:[6,5],
                  11:[None,5],13:[None,5],15:[None,5],17:[6,None],
                  19:[None,None],21:[None,None]}
    for n in range(len(data)):
        for h, limits in thresholds.items():
            limit = limits[n % 2]
            check(bool(fixed[n][h]) == (limit is not None and n >= limit),
                  ('self-converse threshold', n, h))
    return {'status':'FINITE-EXACT: all ordinary classes 0..6, no source filters',
            'classes':[len(x) for x in data],
            'labelled_mass':[sum(factorial(n)//r['aut'] for r in layer)
                             for n, layer in enumerate(data)],
            'F_n':[dict(sorted(c.items())) for c in ordinary],
            'self_converse_F_n':[dict(sorted(c.items())) for c in fixed],
            'strong_H_le_21':low_strong,
            'self_converse_min_order_even_odd_H_le_21':thresholds,
            'controls':['Hamilton subset DP = literal vertex permutations',
                        'ordinary-class labelled mass = 2^binom(n,2)',
                        'ordered-SCC joint coefficient recurrence',
                        'self-converse square-central joint recurrence',
                        'H=1 iff transitive',
                        'n=3 fixed classes have H=1 and H=3',
                        'strong H15 self-converse at 5, not at 6',
                        'H19 attained but never self-converse (using cited floor)']}


if __name__ == '__main__':
    print(json.dumps(run(), indent=2, sort_keys=True))
