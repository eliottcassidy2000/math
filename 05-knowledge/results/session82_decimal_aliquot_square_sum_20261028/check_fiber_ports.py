#!/usr/bin/env python3
"""Complete finite aliquot-fiber square-sum graph audit; Python stdlib only."""
import argparse
import itertools
import json
from math import isqrt
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def divisor_sum(n):
    require(type(n) is int and n >= 1, 'positive integer source required')
    if n == 1:
        return 0
    total = 1
    for d in range(2, isqrt(n) + 1):
        if n % d == 0:
            total += d
            if d * d != n:
                total += n // d
    return total


def complete_fibers(N):
    require(type(N) is int and 37 <= N <= 500, 'audit bound must be 37..500')
    bound = (N - 1) ** 2
    sums = [0] * (bound + 1)
    for d in range(1, bound // 2 + 1):
        for n in range(2 * d, bound + 1, d):
            sums[n] += d
    fibers = {t: [] for t in range(2, N + 1)}
    independent = {t: [] for t in fibers}
    for n in range(1, bound + 1):
        if sums[n] in fibers:
            fibers[sums[n]].append(n)
        t = divisor_sum(n)
        if t in independent:
            independent[t].append(n)
    require(fibers == independent, 'independent complete fiber comparison failed')
    return fibers


def relations(fibers):
    label = {n: t for t, values in fibers.items() for n in values}
    R = {(a, b): set() for a in fibers for b in fibers if a != b}
    for x in label:
        for k in range(isqrt(x) + 1, isqrt(x + max(label, default=0)) + 1):
            y = k * k - x
            if y in label and label[x] != label[y]:
                R[label[x], label[y]].add((x, y))
    for (a, b), pairs in R.items():
        direct = {(x, y) for x in fibers[a] for y in fibers[b]
                  if isqrt(x + y) ** 2 == x + y}
        require(pairs == direct, 'square-generation / pair-testing mismatch')
    return R


def lift(fibers, R, word):
    require(len(word) >= 1 and all(t in fibers for t in word), 'invalid fiber word')
    require(all(a != b for a, b in zip(word, word[1:])), 'adjacent labels must differ')
    reachable = set(fibers[word[0]])
    for a, b in zip(word, word[1:]):
        reachable = {y for x, y in R[a, b] if x in reachable}
    return reachable


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-label', type=int, default=100)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_name('fiber_ports_results.json'))
    args = parser.parse_args()
    F = complete_fibers(args.max_label)
    R = relations(F)
    ports = {b: {a: {y for x, y in R[a, b]} for a in F
                 if a != b and R[a, b]} for b in F}
    bad2 = []
    for b in F:
        for a, c in itertools.combinations(sorted(ports[b]), 2):
            feasible = bool(ports[b][a] & ports[b][c])
            require(feasible == bool(lift(F, R, [a, b, c])), 'port / DP mismatch')
            if not feasible:
                bad2.append((max(a, b, c), a, b, c))
    bad2.sort()
    require([v for v in bad2 if v[0] == 13] == [(13, 4, 13, 10)]
            and bad2[0][0] == 13, 'minimum length-2 obstruction changed')
    require(F[4] == [9] and F[10] == [14] and F[13] == [27, 35], 'seed fibers')
    seed = sorted(F[4] + F[10] + F[13])
    edges = [(x, y) for x, y in itertools.combinations(seed, 2)
             if isqrt(x + y) ** 2 == x + y]
    require(edges == [(9, 27), (14, 35)], 'seed induced edges')
    require(lift(F, R, [4, 13, 14]) == {22, 169}, 'positive lift control')
    # All actual distinct-label length-3 paths, counted once modulo reversal.
    # Orient by b<c (the two middle labels), not by endpoint order.
    bad3 = []
    for b in F:
        for c in ports[b]:
            if b >= c:
                continue
            for a, bp in ports[b].items():
                if a == c:
                    continue
                for d, cp in ports[c].items():
                    if d in (a, b):
                        continue
                    locally = (any(y in bp for y, z in R[b, c]) and
                               any(z in cp for y, z in R[b, c]))
                    globally = any(y in bp and z in cp for y, z in R[b, c])
                    require(globally == bool(lift(F, R, [a, b, c, d])),
                            'length-3 relation / DP mismatch')
                    if locally and not globally:
                        bad3.append((max(a, b, c, d), a, b, c, d))
    bad3.sort()
    require(bad3[0] == (37, 33, 14, 37, 12), 'length-3 minimum changed')
    require(sum(v[0] == 37 for v in bad3) == 1, 'length-3 minimum multiplicity')
    word = [33, 14, 37, 12]
    expected = [{(87, 169)}, {(22, 203), (169, 155)}, {(203, 121)}]
    require([R[a, b] for a, b in zip(word, word[1:])] == expected, 'long witness')
    # Brute Cartesian assignments, independent of relation composition, for
    # every distinct-label word of length 4 through label 13, and the long seed.
    brute_words = list(itertools.permutations(range(2, 14), 4)) + [tuple(word)]
    for w in brute_words:
        brute = {ns[-1] for ns in itertools.product(*(F[t] for t in w))
                 if all(isqrt(x + y) ** 2 == x + y for x, y in zip(ns, ns[1:]))}
        require(brute == lift(F, R, w), 'Cartesian / DP mismatch')
    if args.max_label == 100:
        require(len(bad2) == 4777 and len(bad3) == 1415, 'default census changed')
    result = {
        'scope': 'complete fibers 2..N; distinct-label paths; reversal identified',
        'max_label': args.max_label, 'source_bound': (args.max_label - 1) ** 2,
        'fibers_through_13': {t: F[t] for t in range(2, 14)},
        'nonliftable_length2_count': len(bad2), 'first_length2': bad2[0],
        'locally_liftable_nonliftable_length3_count': len(bad3),
        'first_length3': bad3[0],
        'length3_witness_fibers': {t: F[t] for t in word},
        'length3_witness_relations': [sorted(v) for v in expected],
        'cartesian_word_checks': len(brute_words),
        'independent_checks': ['full divisor-pair scan', 'square-pair adjacency',
                               'Cartesian lifts through label13 plus long seed'],
        'status': 'PASS'}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
