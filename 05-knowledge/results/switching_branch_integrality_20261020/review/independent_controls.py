#!/usr/bin/env python3
"""Independent review controls: no producer or repository imports; explicit checks survive -O."""
from collections import Counter, deque
from fractions import Fraction
from itertools import permutations, product
from math import factorial, gcd, lcm
import json


def require(truth, detail):
    if not truth:
        raise RuntimeError(detail)


def cycle_order(p):
    seen = set()
    result = 1
    for i in range(len(p)):
        if i in seen:
            continue
        j, length = i, 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = p[j]
        result = lcm(result, length)
    return result


def odd_parts(n, smallest=1):
    if not n:
        yield ()
    for x in range(smallest, n + 1, 2):
        for suffix in odd_parts(n-x, x):
            yield (x,) + suffix


def branch_sum(n):
    answer = Fraction(0)
    for parts in odd_parts(n):
        k = len(parts)
        z = 1
        for length, count in Counter(parts).items():
            z *= length**count * factorial(count)
        pair_orbits = sum(x // 2 for x in parts)
        pair_orbits += sum(gcd(parts[i], parts[j]) for i in range(k) for j in range(i))
        answer += Fraction(2**(pair_orbits-k+1), z)
    return answer


def switching_control(n):
    # The normalized class is encoded through its ordered-triple invariant.
    # Work with pullbacks p: the new vertex i is the old vertex p[i].
    edges = [(i, j) for i in range(1, n) for j in range(i+1, n)]
    index = {edge: k for k, edge in enumerate(edges)}

    def direction(q, a, b):
        sign = a > b
        a, b = sorted((a, b))
        return int(sign) ^ (((q >> index[a, b]) & 1) if a else 0)

    def pullback(q, p):
        result = 0
        for k, (i, j) in enumerate(edges):
            a, b, c = p[0], p[i], p[j]
            bit = direction(q, a, b) ^ direction(q, a, c) ^ direction(q, b, c)
            result |= bit << k
        return result

    def matching_flip(q, g):
        # Flipping M changes each normalized triple by its number of matching edges.
        result = q
        for k, (i, j) in enumerate(edges):
            delta = int(g[0] == i) ^ int(g[0] == j) ^ int(g[i] == j)
            result ^= delta << k
        return result

    generators = []
    for i in range(n-1):
        p = list(range(n))
        p[i], p[i+1] = p[i+1], p[i]
        generators.append(p)
    canonical = {}
    representatives = []
    # Discover orbits by adjacent-transposition BFS, not by a full-permutation orbit scan.
    for q in range(1 << len(edges)):
        if q in canonical:
            continue
        representatives.append(q)
        queue = deque([q])
        canonical[q] = q
        while queue:
            x = queue.popleft()
            for p in generators:
                y = pullback(x, p)
                if y not in canonical:
                    canonical[y] = q
                    queue.append(y)
    all_permutations = list(permutations(range(n)))
    orders = {p: cycle_order(p) for p in all_permutations}
    records = []
    weighted = Fraction(0)
    for q in representatives:
        stabilizer = [p for p in all_permutations if pullback(q, p) == q]
        rho = Fraction(sum(orders[p] % 2 == 1 for p in stabilizer), len(stabilizer))
        weighted += rho
        flips = []
        for g in stabilizer:
            if orders[g] != 2:
                continue
            require(all(g[i] != i for i in range(n)), ('nonsemiregular involution', n, q, g))
            target = matching_flip(q, g)
            commuting_witnesses = [
                h for h in all_permutations
                if pullback(q, h) == target and all(h[g[i]] == g[h[i]] for i in range(n))
            ]
            flips.append({'g': g, 'target_mask': target, 'target_type': canonical[target],
                          'commuting_witness': min(commuting_witnesses) if commuting_witnesses else None})
        records.append({'mask': q, 'aut_order': len(stabilizer), 'rho': str(rho), 'flips': flips})
    require(weighted == branch_sum(n), ('stabilizer vs partition', n, weighted, branch_sum(n)))
    # A further independent Burnside count sums direct fixed labelled classes over all permutations.
    direct_odd = Fraction(sum(
        sum(pullback(q, p) == q for q in range(1 << len(edges)))
        for p in all_permutations if orders[p] % 2), factorial(n))
    require(direct_odd == weighted, ('restricted direct Burnside', n))
    return {'n': n, 'labelled_classes': len(canonical), 'type_count': len(records),
            'N_odd': str(weighted), 'N_lev': str(len(records)-weighted), 'types': records}


def group_controls():
    s4 = list(permutations(range(4)))
    a4 = [p for p in s4 if sum(p[i] > p[j] for i in range(4) for j in range(i+1, 4)) % 2 == 0]
    # GL(3,2) on its seven nonzero vectors realizes PSL(2,7).
    gl = []
    for columns in product(range(1, 8), repeat=3):
        image = []
        for v in range(1, 8):
            y = 0
            for i in range(3):
                if v & (1 << i):
                    y ^= columns[i]
            image.append(y)
        if len(set(image)) == 7 and 0 not in image:
            gl.append(tuple(y-1 for y in image))
    groups = {'C2': [(0, 1), (1, 0)], 'A4': a4, 'GL(3,2)=PSL(2,7)': gl}
    expected = {'C2': (2, Fraction(1, 2)), 'A4': (12, Fraction(3, 4)),
                'GL(3,2)=PSL(2,7)': (168, Fraction(5, 8))}
    result = {}
    for name, group in groups.items():
        order_histogram = Counter(map(cycle_order, group))
        odd_count = sum(count for order, count in order_histogram.items() if order % 2)
        rho = Fraction(odd_count, len(group))
        require((len(group), rho) == expected[name], (name, len(group), rho))
        for p in (2, 3, 5, 7):
            prime_part = len(group)
            while prime_part % p == 0:
                prime_part //= p
            root_count = sum(count for order, count in order_histogram.items() if prime_part % order == 0)
            regular_count = sum(count for order, count in order_histogram.items() if order % p)
            require(root_count == regular_count and root_count % prime_part == 0, (name, p))
        result[name] = {'order': len(group), 'element_orders': dict(sorted(order_histogram.items())),
                        'odd_order_elements': odd_count, 'rho': str(rho)}
    # The stated sufficient criterion really is not necessary.
    # C6 acts on a singleton and three points through C6 -> C3.
    # The 2-regular elements (powers 0,2,4) have fixed-point counts 4,1,1.
    terms = [Fraction(4, 6), Fraction(1, 6), Fraction(1, 6)]
    require(sum(terms) == 1 and any(x.denominator % 2 == 0 for x in terms), 'criterion nonnecessity')
    result['criterion_not_necessary'] = {'group': 'C6', 'set': 'one fixed point plus C6/C2',
                                         'prime': 2, 'terms': list(map(str, terms)), 'B_p': '1'}
    return result


def main():
    records = [switching_control(n) for n in (1, 2, 3, 4, 6)]
    by_n = {record['n']: record for record in records}
    require((by_n[2]['N_odd'], by_n[2]['N_lev']) == ('1/2', '1/2'), 'n2 exception')
    require([(c['mask'], c['aut_order'], c['rho']) for c in by_n[4]['types']] ==
            [(0, 4, '1/4'), (2, 12, '3/4')], 'n4 types')
    require(all(f['target_type'] == c['mask'] for c in by_n[4]['types'] for f in c['flips']), 'n4 fixed types')
    expected_pairs = {0: 40, 4: 41, 40: 0, 41: 4}
    require({c['mask']: set(f['target_type'] for f in c['flips'])
             for c in by_n[6]['types'] if c['flips']} == {a: {b} for a, b in expected_pairs.items()}, 'n6 pairing')
    inherited_rho = [Fraction(1, 2)]*19 + [Fraction(1, 4)]*3 + [Fraction(3, 4)] + [Fraction(1, 8)]*3 + [Fraction(5, 8)]
    require(len(inherited_rho) == 27 and sum(1-rho for rho in inherited_rho) == 15, 'n8 inherited arithmetic')
    require(branch_sum(8) == 64, 'n8 odd partition formula')
    output = {'status': 'ALL INDEPENDENT CONTROLS PASSED', 'switching': records, 'groups': group_controls(),
              'n8': {'scope': 'inherited census; arithmetic checked independently',
                     'even_aut_types': len(inherited_rho), 'N_lev': str(sum(1-rho for rho in inherited_rho))}}
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
