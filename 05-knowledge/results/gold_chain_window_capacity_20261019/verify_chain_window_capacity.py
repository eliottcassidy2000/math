#!/usr/bin/env python3
"""Bounded, exact replay for the equal-chain window-capacity continuation.

The note supplies the all-size proofs. No assertion statements are used, so
running under python -O preserves every guard.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import comb
import json

checks = 0


def need(condition, message):
    global checks
    if not condition:
        raise RuntimeError(message)
    checks += 1


def shuffle_words(m):
    for positions in combinations(range(2 * m), m):
        aset = set(positions)
        a = b = 0
        word = []
        for i in range(2 * m):
            if i in aset:
                word.append(a)
                a += 1
            else:
                word.append(m + b)
                b += 1
        yield tuple(word)


def binary(word, m):
    return ''.join('A' if x < m else 'B' for x in word)


def distance(u, v):
    return sum(x != y for x, y in zip(u, v))


def fiber_capacity(words, prefix_length):
    fibers = defaultdict(lambda: [0, 0])
    for w in words:
        fibers[w[prefix_length:]][0 if w[0] == 0 else 1] += 1
    return sum(min(a, b) for a, b in fibers.values()), fibers


def matching_size(adjacency):
    owner = {}

    def augment(a, seen):
        for b in adjacency[a]:
            if b in seen:
                continue
            seen.add(b)
            if b not in owner or augment(owner[b], seen):
                owner[b] = a
                return True
        return False

    for a in range(len(adjacency)):
        augment(a, set())
    return len(owner)


def main():
    rows = []
    prefix_obstructions = []
    for m in range(2, 8):
        words = list(shuffle_words(m))
        C = comb(2 * m, m)
        need(len(words) == C, ('word count', m))
        q = m // 2
        delta = comb(m - 1, q) * comb(m - 1, q - 1)
        actual, fibers = fiber_capacity(words, m + 1)
        symbolic = sum(
            comb(m - 1, a - 1) * min(comb(m, a - 1), comb(m, a))
            for a in range(1, m + 1)
        )
        need(actual == symbolic == C // 2 - delta, ('capacity', m))
        for suffix, (aa, bb) in fibers.items():
            a = m - sum(x < m for x in suffix)
            need(aa == comb(m, a - 1), ('A fiber count', m, suffix))
            need(bb == comb(m, a), ('B fiber count', m, suffix))
        threshold = Fraction(2 * C, 3 * C + 2 * delta)
        rows.append({
            'm': m, 'side': C // 2, 'window_capacity': actual,
            'deficiency': delta, 'gold_threshold': str(threshold),
        })

        extreme = tuple(range(2 * m))
        for k in range(1, 2 * m):
            _, fs = fiber_capacity(words, k)
            aa, bb = fs[extreme[k:]]
            if k <= m:
                need((aa, bb) == (1, 0), ('short prefix', m, k))
            else:
                r = k - m
                need(aa == comb(k - 1, m - 1), ('long prefix A', m, k))
                need(bb == comb(k - 1, m), ('long prefix B', m, k))
                need((aa - bb) * k == (m - r) * comb(k, m),
                     ('prefix imbalance formula', m, k))
            need(aa > bb, ('strict prefix imbalance', m, k))
            prefix_obstructions.append({'m': m, 'k': k, 'A': aa, 'B': bb})

    for m in range(2, 129):
        q = m // 2
        delta = comb(m - 1, q) * comb(m - 1, q - 1)
        positive_sum = sum(
            comb(m - 1, j) * max(0, comb(m, j) - comb(m, j + 1))
            for j in range(m)
        )
        need(positive_sum == delta, ('telescoping identity', m))

    m = 3
    words = list(shuffle_words(m))
    by_name = {binary(w, m): w for w in words}
    left = sorted((w for w in words if w[0] == 0), key=lambda w: binary(w, m))
    right = sorted((w for w in words if w[0] == m), key=lambda w: binary(w, m))
    u, v = by_name['ABBBAA'], by_name['ABBABA']
    expected = [
        ['BAAABB', 6, 6], ['BAABAB', 6, 6], ['BAABBA', 5, 4],
        ['BABAAB', 5, 4], ['BABABA', 4, 2], ['BABBAA', 2, 4],
        ['BBAAAB', 6, 5], ['BBAABA', 5, 3], ['BBABAA', 3, 5],
        ['BBBAAA', 4, 5],
    ]
    actual = [[binary(w, m), distance(u, w), distance(v, w)] for w in right]
    need(actual == expected, 'compression distance table')
    old_n = sum(distance(u, w) <= 5 for w in right)
    new_n = sum(distance(v, w) <= 5 for w in right)
    need((old_n, new_n) == (7, 8), 'compression increases neighborhood')

    def neighborhood(family):
        return {binary(w, m) for w in right
                if any(distance(x, w) <= 5 for x in family)}

    special_names = ['AAABBB', 'AABABB', 'AABBAB', 'ABAABB']
    special_n = neighborhood([by_name[name] for name in special_names])
    lex_n = neighborhood(left[:4])
    need((len(special_n), len(lex_n)) == (5, 10), 'lex minimization fails')
    adjacency = [[j for j, b in enumerate(right) if distance(a, b) <= 4]
                 for a in left]
    need(matching_size(adjacency) == 9, 'complete support-four capacity at m=3')

    def complement(w):
        return tuple(x + m if x < m else x - m for x in w)

    def is_dyck(w):
        height = 0
        for x in w:
            height += 1 if x < m else -1
            if height < 0:
                return False
        return True

    graph = [{j for j, b in enumerate(left) if distance(a, complement(b)) <= 5}
             for a in left]
    for i in range(len(left)):
        need(i not in graph[i], ('no complement loop', i))
        for j in range(len(left)):
            need((j in graph[i]) == (i in graph[j]), ('symmetric graph', i, j))
    independent = {i for i, w in enumerate(left) if is_dyck(w)}
    independent_n = set().union(*(graph[i] for i in independent))
    need(not independent.intersection(independent_n), 'Dyck family independent')
    need((len(independent), len(independent_n)) == (5, 5), 'Dyck Hall tightness')

    receipt = {
        'status': 'passed', 'guard_count': checks,
        'enumerated_chain_lengths': list(range(2, 8)),
        'binomial_identity_range': [2, 128],
        'capacity_rows': rows,
        'proper_prefix_obstructions': prefix_obstructions,
        'compression_distances': actual,
        'compression_neighborhood_sizes': [old_n, new_n],
        'lex_initial_four': [binary(w, m) for w in left[:4]],
        'lex_initial_neighborhood_size': len(lex_n),
        'alternative_four': special_names,
        'alternative_neighborhood': sorted(special_n),
        'm3_complete_support4_capacity': 9,
        'm3_independent_dyck_family_size': len(independent),
        'm3_independent_dyck_neighborhood_size': len(independent_n),
        'all_m_support_m_plus_2': 'open in this work',
    }
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
