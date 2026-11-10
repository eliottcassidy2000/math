#!/usr/bin/env python3
"""Small independent controls for recovered numerical decompositions.

No imports from the historical producers. Literal permutation counts and
subset endpoint DP agree on the complete 64 order-five fixed-path cube.
"""
from collections import Counter
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SOURCES = [
    '01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md',
    '05-knowledge/results/trijection_reset_20261020/arithmetic.md',
    '05-knowledge/results/tournament_joint_parity_20261102/parity.md',
    '05-knowledge/results/tournament_joint_parity_20261102/growth_obstruction.md',
    '05-knowledge/results/groups_of_five_primegame_merged_metagraph_20261031.md',
    '05-knowledge/results/boundary_carry_returns_20261103/README.md',
    '05-knowledge/results/six_roots_hidden_bits_20261103/NUMBER_ATLAS.md',
]


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def tournament(bits):
    # Fixed directed path 0 -> 1 -> 2 -> 3 -> 4.
    # The first four coordinates are the two mirror-moved pairs in THM-4733.
    edges = [(0, 2), (2, 4), (0, 3), (1, 4), (0, 4), (1, 3)]
    rows = [sum(1 << j for j in range(i + 1, 5)) for i in range(5)]
    for bit, (i, j) in zip(bits, edges):
        if bit:
            rows[i] ^= 1 << j
            rows[j] ^= 1 << i
    return rows


def literal(rows):
    return sum(all(rows[u] >> v & 1 for u, v in zip(p, p[1:]))
               for p in permutations(range(len(rows))))


def dynamic(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last, value in enumerate(dp[mask]):
            if value:
                for v in range(n):
                    if not mask >> v & 1 and rows[last] >> v & 1:
                        dp[mask | 1 << v][v] += value
    return sum(dp[-1])


def pairs(counter):
    return [[a, b] for a, b in sorted(counter.items())]


def main():
    all_counts = {}
    blue = Counter()
    for bits in product(range(2), repeat=6):
        rows = tournament(bits)
        h = literal(rows)
        require(h == dynamic(rows), ('two counters', bits))
        require(h % 2 == 1, ('tournament parity', bits))
        all_counts[bits] = h
        if bits[0] == bits[1] and bits[2] == bits[3]:
            blue[h] += 1
    require(blue == {1: 1, 3: 1, 9: 4, 11: 3, 13: 3, 15: 4}, 'blue histogram')

    fold_failures = []
    for fixed in product(range(2), repeat=2):
        for q, r in product(range(2), repeat=2):
            left = (1, 0, q, r) + fixed
            right = (0, 1, q, r) + fixed
            if all_counts[left] != all_counts[right]:
                fold_failures.append({'fixed': fixed, 'q_r': (q, r),
                    'left_rows': tournament(left), 'right_rows': tournament(right),
                    'H': [all_counts[left], all_counts[right]]})
    require(len(fold_failures) == 8, 'complete fold-failure count')

    divisors = [(a, b, c) for a in range(3) for b in range(2) for c in range(2)
                if (a, b, c) not in [(0, 0, 0), (2, 1, 1)]]
    ranks = Counter(a + b + c for a, b, c in divisors)
    heights = Counter(a - b - c for a, b, c in divisors)
    require(ranks == {1: 3, 2: 4, 3: 3}, 'proper-divisor ranks')
    require(heights == {-2: 1, -1: 3, 0: 2, 1: 3, 2: 1}, 'proper-divisor heights')
    fold = Counter((x[0] + x[1], x[2], x[3]) for x in product(range(2), repeat=4))
    require(len(fold) == 12 and Counter(fold.values()) == {1: 8, 2: 4}, 'fold fibres')

    # An explicitly NEW example of the requested histogram, not a recovered
    # definition: move one middle complement pair to ranks 1 and 3.
    example_pair = ((1, 1, 0, 0), (0, 0, 1, 1))
    regroupings = []
    for a in product(range(2), repeat=4):
        if sum(a) != 2:
            continue
        b = tuple(1 - v for v in a)
        hist = Counter()
        for x in product(range(2), repeat=4):
            g = 1 if x == a else 3 if x == b else sum(x)
            xc = tuple(1 - v for v in x)
            gc = 1 if xc == a else 3 if xc == b else sum(xc)
            require(gc == 4 - g, 'complement equivariance')
            hist[g] += 1
        require([hist[k] for k in range(5)] == [1, 5, 4, 5, 1], 'regrading')
        regroupings.append([a, b])

    result = {
        'status': 'PASS',
        'scope': '64 fixed-path order-five tournaments; 16 blue states; 16 four-bit words; 10 proper divisor exponent vectors',
        'literal_permutation_controls': 64,
        'blue_H_histogram': pairs(blue),
        'fold_fibre_sizes': pairs(Counter(fold.values())),
        'fold_H_failures': fold_failures,
        'divisor_rank_histogram': pairs(ranks),
        'divisor_height_histogram': pairs(heights),
        'new_regrouping_example': example_pair,
        'complement_symmetric_regroupings_checked': len(regroupings),
        'regrouped_histogram': [1, 5, 4, 5, 1],
        'provenance_limit': 'No original definition of the requested 1-5-4-5-1 grading was found in the bounded source recovery. The explicit regrouping is new and is not attributed to prior work.',
        'source_sha256': {path: sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCES},
        'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
