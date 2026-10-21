#!/usr/bin/env python3
"""Independent reviewer: prime-cycle inclusion/exclusion and signed paths.

No producer modules or receipts are imported. Output is deterministic.
All checks use explicit raises, including under python -O.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
from math import lcm
from pathlib import Path
import argparse
import json


def require(test, explanation):
    if not test:
        raise RuntimeError(explanation)


def compose(p, q):
    return tuple(p[x] for x in q)


def generated(generators):
    n = len(generators[0])
    identity = tuple(range(n))
    seen = {identity}
    queue = [identity]
    for p in queue:
        for a in generators:
            q = compose(a, p)
            if q not in seen:
                seen.add(q)
                queue.append(q)
    return sorted(seen)


def cycles(p):
    unused = set(range(len(p)))
    result = []
    while unused:
        x = min(unused)
        block = []
        while x in unused:
            unused.remove(x)
            block.append(x)
            x = p[x]
        result.append(tuple(sorted(block)))
    return tuple(sorted(result))


def join(P, Q):
    blocks = [set(block) for block in P]
    for block in Q:
        hit = [a for a in blocks if a.intersection(block)]
        merged = set(block).union(*hit)
        blocks = [a for a in blocks if not a.intersection(block)] + [merged]
    return tuple(sorted(tuple(sorted(a)) for a in blocks))


def prime(n):
    return n >= 2 and all(n % d for d in range(2, int(n ** .5) + 1))


def regular_color_polynomial(G):
    """Avoid named words and set-partition enumeration.

    A coloring is nonregular iff fixed by an element of prime order.
    Inclusion-exclusion on their cycle equivalences counts the complement;
    intersections of fixed sets are the joins of equivalence partitions.
    """
    n = len(G[0])
    constraints = set()
    for p in G:
        P = cycles(p)
        if prime(lcm(*(len(c) for c in P))):
            constraints.add(P)
    terms = Counter({tuple((i,) for i in range(n)): 1})
    for Q in sorted(constraints):
        update = Counter()
        for P, coefficient in terms.items():
            update[join(P, Q)] -= coefficient
        terms.update(update)
        terms = Counter({P: c for P, c in terms.items() if c})
    polynomial = Counter()
    for P, c in terms.items():
        polynomial[len(P)] += c
    polynomial = {k: c for k, c in sorted(polynomial.items()) if c}
    evaluate = lambda q: sum(c * q ** k for k, c in polynomial.items())
    return {
        'group_order': len(G),
        'distinct_prime_cycle_constraints': len(constraints),
        'regular_color_polynomial_coefficients': polynomial,
        'regular_named_colorings_q_1_to_4': [evaluate(q) for q in range(1, 5)]
    }


def fixes_color(p, c):
    return all(c[i] == c[p[i]] for i in range(len(p)))


def signed_path_audit():
    records = []
    for n in range(1, 6):
        count = 0
        kernel_objects = Counter()
        negative_objects = Counter()
        pairs = list(combinations(range(n), 2))
        all_perms = list(permutations(range(n)))
        for encoding in range(1 << len(pairs)):
            arcs = {(j, i) if encoding >> b & 1 else (i, j)
                    for b, (i, j) in enumerate(pairs)}
            positive = [p for p in all_perms
                        if all((p[i], p[j]) in arcs for i, j in arcs)]
            negative = [p for p in all_perms
                        if all((p[j], p[i]) in arcs for i, j in arcs)]
            for k in range(1, n + 1):
                for path in permutations(range(n), k):
                    if any((x, y) not in arcs for x, y in zip(path, path[1:])):
                        continue
                    U = [x for x in range(n) if x not in path]
                    index = {x: i for i, x in enumerate(U)}
                    H = [p for p in positive if all(p[x] == x for x in path)]
                    A = [p for p in negative
                         if all(p[x] == y for x, y in zip(path, reversed(path)))]
                    hu = [tuple(index[p[x]] for x in U) for p in H]
                    au = [tuple(index[p[x]] for x in U) for p in A]
                    require(len(hu) == len(set(hu)), 'ordinary residual faithful')
                    require(len(H) % 2 == 1, 'ordinary residual odd')
                    binary = [c for c in product(range(2), repeat=len(U))
                              if sum(fixes_color(p, c) for p in hu) == 1]
                    require(binary, 'ordinary residual binary existence')
                    require(min(sum(c) for c in binary) <= len(U) // 2,
                            'ordinary residual minority bound')
                    kernel = len(hu) + len(au) - len(set(hu + au))
                    if A:
                        negative_objects[len(U)] += 1
                    if kernel:
                        kernel_objects[len(U)] += 1
                    if len(U) >= 2:
                        require(not kernel, 'signed omitted action faithful')
                        for c in binary:
                            remainder = [p for p in au if fixes_color(p, c)]
                            require(len(remainder) <= 1, 'at most one anti symmetry')
                            if remainder:
                                a = remainder[0]
                                require(compose(a, a) == tuple(range(len(U))),
                                        'remaining anti map involutory')
                                v = next(i for i in range(len(U)) if a[i] != i)
                                refined = tuple(2 if i == v else c[i]
                                                for i in range(len(U)))
                                require(sum(fixes_color(p, refined) for p in hu + au) == 1,
                                        'third-color refinement regular')
                    if n >= 2 and len(U) == 1 and k % 2:
                        require(not A, 'odd marked path and singleton omission forbid anti')
                    count += 1
        records.append({'n': n, 'ordered_path_objects': count,
                        'with_negative_map_by_omitted_size': dict(sorted(negative_objects.items())),
                        'with_signed_kernel_by_omitted_size': dict(sorted(kernel_objects.items()))})
    return records


def main():
    translation = tuple((i + 1) % 7 for i in range(7))
    multiplier2 = tuple(2 * i % 7 for i in range(7))
    multiplier3 = tuple(3 * i % 7 for i in range(7))
    ordinary = generated([translation, multiplier2])
    signed = generated([translation, multiplier3])
    # GL(3,2) generated by swapping first two coordinates, cycling all three,
    # and adding coordinate 2 to coordinate 1. Carrier is integers 1,...,7.
    def vperm(function):
        return tuple(function(x) - 1 for x in range(1, 8))
    swap = vperm(lambda x: ((x & 1) << 1) | ((x & 2) >> 1) | (x & 4))
    rotate = vperm(lambda x: ((x << 1) & 7) | (x >> 2))
    shear = vperm(lambda x: x ^ ((x & 2) >> 1))
    fano = generated([swap, rotate, shear])
    require([len(ordinary), len(signed), len(fano)] == [21, 42, 168], 'group orders')
    reports = {name: regular_color_polynomial(G) for name, G in
               [('ordinary_paley', ordinary), ('signed_paley', signed), ('fano', fano)]}
    require(reports['ordinary_paley']['regular_named_colorings_q_1_to_4'] == [0, 84, 2016, 15960],
            'ordinary named-color census')
    require(reports['signed_paley']['regular_named_colorings_q_1_to_4'] == [0, 0, 1512, 14280],
            'signed named-color census')
    require(reports['fano']['regular_named_colorings_q_1_to_4'] == [0, 0, 0, 4032],
            'Fano named-color census')
    mapping = [1, 2, 4, 3, 6, 7, 5]
    plus = {frozenset((b + q) % 7 for q in (1, 2, 4)) for b in range(7)}
    minus = {frozenset((b - q) % 7 for q in (1, 2, 4)) for b in range(7)}
    require(not plus.intersection(minus), 'two distinct Fano systems')
    require(all(mapping[a] ^ mapping[b] ^ mapping[c] == 0 for a, b, c in plus),
            'published cyclic-to-vector identification')
    for p in ordinary:
        require({frozenset(p[x] for x in line) for line in plus} == plus,
                'ordinary preserves selected Fano system')
    for p in set(signed).difference(ordinary):
        require({frozenset(p[x] for x in line) for line in plus} == minus,
                'negative swaps selected Fano systems')
    reports['signed_partial_paths'] = signed_path_audit()
    reports['source_sha256'] = sha256(Path(__file__).read_bytes()).hexdigest()
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    Path(args.output).write_text(json.dumps(reports, indent=2) + '\n')
    print(json.dumps(reports, indent=2))


if __name__ == '__main__':
    main()
