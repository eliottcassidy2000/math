"""Exact controls for Sylow tower classification and composition deficits.

Run from any directory; stdlib only. Checks remain active under python -O.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from itertools import combinations, permutations
import json
from math import factorial
from pathlib import Path


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def digits(n, p):
    out = []
    while n:
        out.append(n % p)
        n //= p
    return out


def val(n, p):
    answer = 0
    while n and n % p == 0:
        answer += 1
        n //= p
    return answer


def factorial_val(n, p):
    answer = 0
    while n:
        n //= p
        answer += n
    return answer


def tower_generators(p, k):
    n = p**k
    if not k:
        return []
    m = n // p
    answer = []
    for i in range(p):
        for g in tower_generators(p, k - 1):
            h = list(range(n))
            for j in range(m):
                h[i*m+j] = i*m+g[j]
            answer.append(tuple(h))
    answer.append(tuple(((i//m+1) % p)*m+i % m for i in range(n)))
    return answer


def standard_sylow(n, p):
    answer, blocks, offset = [], [], 0
    for k, multiplicity in enumerate(digits(n, p)):
        for _ in range(multiplicity):
            m = p**k
            blocks.append((offset, k))
            for g in tower_generators(p, k):
                h = list(range(n))
                for j in range(m):
                    h[offset+j] = offset+g[j]
                answer.append(tuple(h))
            offset += m
    need(offset == n, 'digit blocks partition the vertices')
    return answer, blocks


def signed_pair_orbits(n, generators):
    """Independent signed BFS on actual vertex-pair action, no tower formula."""
    remaining = set(combinations(range(n), 2))
    orbits = []
    while remaining:
        root = min(remaining)
        signs, queue = {root: 0}, deque([root])
        while queue:
            a, b = queue.popleft()
            for g in generators:
                x, y = g[a], g[b]
                pair = (min(x, y), max(x, y))
                sign = signs[(a, b)] ^ int(x > y)
                if pair in signs:
                    need(signs[pair] == sign, 'odd group cannot reverse a pair')
                else:
                    signs[pair] = sign
                    queue.append(pair)
        remaining.difference_update(signs)
        orbits.append(signs)
    return orbits


def from_orbits(n, orbits, mask):
    A = [[False]*n for _ in range(n)]
    for i, orbit in enumerate(orbits):
        bit = (mask >> i) & 1
        for (a, b), sign in orbit.items():
            A[a][b] = bool(bit ^ sign)
            A[b][a] = not A[a][b]
    return A


def tower_shape(A, offset, p, k):
    if not k:
        return
    m = p**(k-1)
    first = [[A[offset+i][offset+j] for j in range(m)] for i in range(m)]
    for block in range(p):
        o = offset+block*m
        need([[A[o+i][o+j] for j in range(m)] for i in range(m)] == first,
             'children are identical in standard coordinates')
    for i in range(p):
        for j in range(i+1, p):
            bit = A[offset+i*m][offset+j*m]
            need(all(A[offset+i*m+a][offset+j*m+b] == bit
                     for a in range(m) for b in range(m)), 'cross arcs uniform')
            need(bit == A[offset][offset+((j-i) % p)*m], 'quotient is circulant')
    tower_shape(A, offset, p, k-1)


def validate_shape(A, p, blocks):
    for offset, k in blocks:
        tower_shape(A, offset, p, k)
    for (o, k), (r, ell) in combinations(blocks, 2):
        need(all(A[o+i][r+j] == A[o][r]
                 for i in range(p**k) for j in range(p**ell)),
             'different Sylow orbits are modules')


def automorphisms(A):
    """Enumerate all arc-preserving bijections, using only scores and arcs."""
    n = len(A)
    scores = [sum(row) for row in A]
    candidates = [[j for j in range(n) if scores[i] == scores[j]] for i in range(n)]
    order = sorted(range(n), key=lambda i: (len(candidates[i]), i))
    mapping, used, result = [-1]*n, set(), []

    def search(depth):
        if depth == n:
            result.append(tuple(mapping))
            return
        x = order[depth]
        for y in candidates[x]:
            if y in used:
                continue
            if all(A[x][z] == A[y][mapping[z]] for z in order[:depth]):
                mapping[x] = y
                used.add(y)
                search(depth+1)
                used.remove(y)
                mapping[x] = -1
    search(0)
    return result


def transitive(n):
    return [[i < j for j in range(n)] for i in range(n)]


def cycle(p):
    return [[0 < (j-i) % p <= (p-1)//2 for j in range(p)] for i in range(p)]


def strong_quotient(n):
    need(n >= 3, 'strong quotient needs at least three vertices')
    A = transitive(n)
    A[0][n-1], A[n-1][0] = False, True
    return A


def strong(A):
    for source in range(len(A)):
        reached, queue = {source}, deque([source])
        while queue:
            i = queue.popleft()
            for j in range(len(A)):
                if A[i][j] and j not in reached:
                    reached.add(j)
                    queue.append(j)
        if len(reached) != len(A):
            return False
    return True


def tower(p, k):
    if k == 0:
        return transitive(1)
    return substitute(cycle(p), [tower(p, k-1)]*p)


def substitute(Q, children):
    locations = [(i, j) for i, child in enumerate(children) for j in range(len(child))]
    return [[children[i][a][b] if i == j else Q[i][j]
             for j, b in locations] for i, a in locations]


def compose(g, h):
    return tuple(g[x] for x in h)


def element_order(g):
    h, count = tuple(range(len(g))), 0
    while True:
        h = compose(g, h)
        count += 1
        if h == tuple(range(len(g))):
            return count


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    report = {'status': 'FINITE-EXACT', 'universe': {}, 'pair_orbit_checks': [],
              'exhaustive_fixed_tournaments': [], 'examples': {}}
    # No inherited graph/census filtering. Actual standard Sylow generators.
    for p in (3, 5, 7, 11):
        for n in range(1, 41):
            generators, blocks = standard_sylow(n, p)
            orbits = signed_pair_orbits(n, generators)
            ds = digits(n, p)
            s = sum(ds)
            predicted = s*(s-1)//2 + (p-1)//2*sum(k*a for k, a in enumerate(ds))
            need(len(orbits) == predicted, f'pair orbit formula p={p},n={n}')
            need(factorial_val(n, p) == (n-s)//(p-1), 'digit factorial identity')
            if n <= 12 and len(orbits) <= 12:
                for mask in range(1 << len(orbits)):
                    validate_shape(from_orbits(n, orbits, mask), p, blocks)
                report['exhaustive_fixed_tournaments'].append(
                    {'p': p, 'n': n, 'dimension': predicted, 'checked': 1 << predicted})
            else:
                for mask in (0, (1 << len(orbits))-1, sum(1 << j for j in range(0, len(orbits), 2))):
                    validate_shape(from_orbits(n, orbits, mask), p, blocks)
            report['pair_orbit_checks'].append({'p': p, 'n': n, 'dimension': predicted})
    # Independent complete tournament universe through n=5, all vertex permutations.
    maxima = []
    for n in range(1, 6):
        pairs = [{pair: 0} for pair in combinations(range(n), 2)]
        peak = {p: 0 for p in (2, 3, 5)}
        strong_peak = {p: -1 for p in (2, 3, 5)}
        for mask in range(1 << len(pairs)):
            A = from_orbits(n, pairs, mask)
            group = automorphisms(A)
            brute = [g for g in permutations(range(n))
                     if all(A[g[i]][g[j]] == A[i][j] for i, j in combinations(range(n), 2))]
            need(set(group) == set(brute), 'independent full-permutation automorphism census')
            for p in peak:
                peak[p] = max(peak[p], val(len(group), p))
                if strong(A):
                    strong_peak[p] = max(strong_peak[p], val(len(group), p))
        need(peak[2] == 0, 'no even automorphism groups')
        for p in (3, 5):
            need(peak[p] == factorial_val(n, p), 'sharp odd-prime maximum')
            if n >= 3:
                need(strong_peak[p] == factorial_val(n, p)-int(sum(digits(n, p)) == 2),
                     'sharp strong odd-prime maximum')
        maxima.append({'n': n, 'max_valuations': peak, 'strong_max_valuations': strong_peak})
    report['all_tournaments_through_five'] = maxima
    c3, one = cycle(3), transitive(1)
    r9 = substitute(c3, [c3]*3)
    examples = {
        'saturated_12': substitute(transitive(2), [r9, c3]),
        'order81_nonabelian_12': substitute(transitive(4), [r9, one, one, one]),
        'order81_abelian_12': substitute(transitive(4), [c3]*4),
        'carry_loss_6': substitute(c3, [transitive(2)]*3),
        'decoration_loss_9': substitute(c3, [c3, transitive(3), transitive(3)]),
        'saturated_9': r9,
        'strong_extremal_12': substitute(c3, [substitute(transitive(2), [c3, one])]*3),
        'paley_7': [[(j-i) % 7 in (1, 2, 4) for j in range(7)] for i in range(7)],
    }
    expected = {'saturated_12': 243, 'order81_nonabelian_12': 81,
                'order81_abelian_12': 81, 'carry_loss_6': 3,
                'decoration_loss_9': 3, 'saturated_9': 81,
                'strong_extremal_12': 81, 'paley_7': 21}
    groups = {}
    for name, A in examples.items():
        group = automorphisms(A)
        groups[name] = group
        need(len(group) == expected[name], f'automorphism order {name}')
        abelian = all(compose(g, h) == compose(h, g) for g in group for h in group)
        report['examples'][name] = {
            'n': len(A), 'automorphism_order': len(group), 'abelian': abelian,
            'strong': strong(A),
            'element_order_histogram': dict(sorted(Counter(element_order(g) for g in group).items())),
            'delta3': factorial_val(len(A), 3)-val(len(group), 3),
            'upper_triangle': ''.join(str(int(A[i][j])) for i, j in combinations(range(len(A)), 2))}
    need(not report['examples']['order81_nonabelian_12']['abelian'], 'nonabelian hostile')
    need(report['examples']['order81_abelian_12']['abelian'], 'abelian hostile')
    need(report['examples']['strong_extremal_12']['strong'], 'strong extremal positive')
    need(not report['examples']['saturated_12']['strong'], 'saturation at digit sum two hostile')
    paley = examples['paley_7']
    proper_modules = []
    for mask in range(1, 1 << 7):
        inside = [i for i in range(7) if mask & (1 << i)]
        if 2 <= len(inside) <= 6 and all(len({paley[j][i] for i in inside}) == 1
                                        for j in range(7) if j not in inside):
            proper_modules.append(inside)
    need(not proper_modules, 'Paley7 is a prime tournament')
    report['examples']['paley_7']['nontrivial_proper_modules'] = proper_modules
    strong_constructions = []
    for p in (3, 5, 7, 11):
        for n in range(3, 41):
            exponents = [k for k, a in enumerate(digits(n, p)) for _ in range(a)]
            defect = int(len(exponents) == 2)
            if defect:
                k = max(exponents)
                need(k > 0, 'n>=3 excludes the impossible strong size two')
                exponents.remove(k)
                exponents.extend([k-1]*p)
            children = [tower(p, k) for k in exponents]
            A = children[0] if len(children) == 1 else substitute(strong_quotient(len(children)), children)
            need(len(A) == n and strong(A), 'strong attaining construction')
            lower = sum((p**k-1)//(p-1) for k in exponents)
            need(lower == factorial_val(n, p)-defect, 'strong subgroup exact exponent')
            offset = 0
            for k in exponents:
                for g in tower_generators(p, k):
                    h = list(range(n))
                    for j in range(p**k):
                        h[offset+j] = offset+g[j]
                    need(all(A[h[i]][h[j]] == A[i][j] for i, j in combinations(range(n), 2)),
                         'every claimed strong subgroup generator preserves actual arcs')
                offset += p**k
            strong_constructions.append({'n': n, 'p': p, 'defect': defect,
                                         'subgroup_valuation': lower})
    report['strong_constructions'] = strong_constructions
    # Arithmetic carry recurrence independently checked with integer factorials.
    count = 0
    for p in (3, 5, 7, 11):
        for m in range(1, 101):
            quotient = factorial(p*m)//factorial(m)**p
            need(val(quotient, p)-1 == sum(digits(m, p))-1, 'homogeneous carry recurrence')
            count += 1
    report['homogeneous_carry_controls'] = count
    report['universe'] = {
        'pair_orbits': 'p in {3,5,7,11}; every n=1..40; no filters',
        'fixed_tournaments': 'all signed-orbit assignments if n<=12 and dimension<=12; three explicit controls otherwise',
        'maxima': 'all labelled tournaments n=1..5; all permutations independently enumerated',
        'examples': 'eight explicit arc matrices; all automorphisms enumerated from scores and arcs',
        'strong_constructions': 'p in {3,5,7,11}; every n=3..40; actual strong reachability and every claimed subgroup generator',
        'arithmetic': 'integer-factorial vs base-p digit identity, four primes and m=1..100'}
    with args.out.open('w', encoding='utf-8', newline='\n') as output:
        output.write(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'output': str(args.out), 'pair_orbit_cases': len(report['pair_orbit_checks']),
                      'fixed_objects': sum(r['checked'] for r in report['exhaustive_fixed_tournaments']),
                      'example_orders': {k: r['automorphism_order'] for k, r in report['examples'].items()}}, sort_keys=True))


if __name__ == '__main__':
    main()
