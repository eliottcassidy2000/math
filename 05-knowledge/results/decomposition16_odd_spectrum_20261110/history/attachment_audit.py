#!/usr/bin/env python3
"""Backward recursion and permutation checks of the repaired H3 compiler."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def need(ok, description):
    if not ok:
        raise RuntimeError(description)


def source(a, b):
    n = a + b + 3
    matrix = [[int(i < j) for j in range(n)] for i in range(n)]
    matrix[a][a + 2], matrix[a + 2][a] = 0, 1
    return matrix


def adjoin(matrix, bits):
    return [row + [1 - bits[i]] for i, row in enumerate(matrix)] + [list(bits) + [0]]


def path_count(matrix):
    n = len(matrix)

    @lru_cache(None)
    def suffix(v, remaining):
        if not remaining:
            return 1
        return sum(suffix(j, remaining ^ (1 << j)) for j in range(n)
                   if remaining >> j & 1 and matrix[v][j])

    return sum(suffix(v, ((1 << n) - 1) ^ (1 << v)) for v in range(n))


def literal(matrix):
    return sum(all(matrix[u][v] for u, v in zip(p, p[1:]))
               for p in permutations(range(len(matrix))))


def chain_count(integer):
    # Endpoint-gap sum, independent of the producer's integer recurrence.
    word = [int(s) for s in bin(integer)[2:]] if integer else []
    nonadjacent = sum(word[i] * (1 - word[j]) * 2 ** (j - i - 1)
                      for i in range(len(word)) for j in range(i + 1, len(word)))
    adjacent = sum(word[i] * (1 - word[i + 1]) for i in range(len(word) - 1))
    return 1 + nonadjacent + adjacent


def canonical(A, B, k):
    left = [int(s) for s in bin(A)[2:]] if A else []
    right = [1 - int(s) for s in bin(B)[2:][::-1]] if B else []
    return adjoin(source(len(left), len(right)), left + [int(i < k) for i in range(3)] + right)


def main():
    count = literal_checks = uniform_checks = 0
    histograms = {}
    for total in range(5):
        for a in range(total + 1):
            b = total - a
            old = source(a, b)
            hist = Counter()
            for bits in product(range(2), repeat=len(old)):
                A = sum(bits[i] << (a - 1 - i) for i in range(a))
                B = sum((1 - bits[a + 3 + i]) << i for i in range(b))
                k = sum(bits[a:a + 3])
                matrix = adjoin(old, bits)
                got = path_count(matrix)
                expected = (12 * A * B + 3 * chain_count(A) + 3 * chain_count(B) - 3
                            + 3 * (3 - k) * A + 3 * k * B
                            + (3 - k) * (A % 2) + k * (B % 2)
                            + 2 * (k in (1, 2)))
                need(got == expected, ('all-contact formula', a, b, bits, got, expected))
                if len(matrix) <= 6:
                    need(got == literal(matrix), 'literal independent control')
                    literal_checks += 1
                if k in (0, 3):
                    word = list(bits[:a]) + [k // 3] + list(bits[a + 3:])
                    chain = [[int(i < j) for j in range(len(word))] for i in range(len(word))]
                    collapsed = adjoin(chain, word)
                    correction = 6 * A * (B + 1) if k == 0 else 6 * B * (A + 1)
                    need(got == 3 * path_count(collapsed) + correction, 'uniform reentry correction')
                    need(got % 3 == 0, 'uniform divisible by three')
                    uniform_checks += 1
                hist[got] += 1
                count += 1
            histograms[f'{a},{b}'] = dict(sorted(hist.items()))

    hostile_rows = [30, 20, 24, 18, 32, 15]
    hostile = [[r >> j & 1 for j in range(6)] for r in hostile_rows]
    split_path = (1, 2, 4, 5, 0, 3)
    need(all(hostile[u][v] for u, v in zip(split_path, split_path[1:])), 'split C3 witness')
    need(literal(hostile) == 27, 'hostile full count')
    collapsed_rows = [6, 4, 8, 3]
    collapsed = [[r >> j & 1 for j in range(4)] for r in collapsed_rows]
    need(literal(collapsed) == 5, 'hostile collapsed count')

    cap = 55
    apex = set()
    target_graph_checks = 0
    for q in range(cap):
        bits = [int(s) for s in bin(2 * q)[2:]] if q else []
        chain = [[int(i < j) for j in range(len(bits))] for i in range(len(bits))]
        value = path_count(adjoin(chain, bits))
        need(value == chain_count(2 * q), 'one-apex target graph')
        if value <= cap:
            apex.add(value)
        target_graph_checks += 1
    mixed = set()
    for A in range((cap - 5) // 6 + 1):
        for B in range((cap - 5 - 6 * A) // (12 * A + 3) + 1):
            value = path_count(canonical(A, B, 1))
            expected = 12*A*B + 6*A + 3*B + 3*chain_count(A) + 3*chain_count(B) + 2*(A%2) + B%2 - 1
            need(value == expected, 'mixed target graph')
            if value <= cap:
                mixed.add(value)
            target_graph_checks += 1
    uniform = {3*v for v in apex if 3*v <= cap}
    for A in range(1, (cap - 3) // 9 + 1):
        for B in range((cap - 3 - 9*A) // (12*A) + 1):
            value = path_count(canonical(A, B, 0))
            expected = 12*A*B + 9*A + 3*chain_count(A) + 3*chain_count(B) + 3*(A%2) - 3
            need(value == expected, 'uniform target graph')
            if value <= cap:
                uniform.add(value)
            target_graph_checks += 1
    misses = [h for h in range(1, cap + 1, 2) if h not in apex | mixed | uniform]
    need(misses == [7, 21, 55], 'first extra parent-lane obstruction')
    witness55_rows = [54, 124, 56, 113, 96, 64, 5]
    witness55 = [[r >> j & 1 for j in range(7)] for r in witness55_rows]
    need(literal(witness55) == path_count(witness55) == 55, 'H55 witness')
    deletions55 = []
    for v in range(7):
        keep = [i for i in range(7) if i != v]
        matrix = [[witness55[i][j] for j in keep] for i in keep]
        value = literal(matrix)
        need(value == path_count(matrix), 'H55 deletion')
        deletions55.append(value)
    need(deletions55 == [9, 25, 15, 19, 29, 29, 5], 'H55 minimum is sharp')
    print(json.dumps({
        'status': 'PASS',
        'universe': 'every cut on every TT_a -> C3 -> TT_b with a+b<=4',
        'direct_recursion_controls': count,
        'permutation_controls': literal_checks,
        'uniform_reentry_controls': uniform_checks,
        'target_inverse_graph_controls_through55': target_graph_checks,
        'parent1_or3_missing_odds_through55': misses,
        'sharp_H55_witness': {'rows': witness55_rows, 'H': 55, 'deletions': deletions55},
        'histograms': histograms,
        'refuted_draft_uniform_module_claim': {
            'rows': hostile_rows, 'H': 27, 'collapsed_rows': collapsed_rows,
            'collapsed_H': 5, 'split_triangle_path': split_path,
            'repaired_excess': 12,
        },
        'source_sha256': {
            str(path.relative_to(HERE.parent)): sha256(path.read_bytes()).hexdigest()
            for path in [HERE.parent / 'attachment/PROOF.md', HERE.parent / 'attachment/check.py']
        },
        'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
