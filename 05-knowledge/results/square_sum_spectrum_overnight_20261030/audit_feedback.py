#!/usr/bin/env python3
"""Consumer for every positive witness, plus independent subset-DP controls.

No importing the C++ scanner or the matrix-construction script.
"""
import argparse
from collections import defaultdict
from itertools import permutations, product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def decode(length, number):
    require(0 <= number < 4 ** length, 'encoded word out of range')
    return tuple((number >> (2 * i)) & 3 for i in range(length - 1, -1, -1))


def gaps(word):
    total = 0
    for left, right in ((0, 1), (1, 0)):
        current = {(0, 0, 0): 1}
        for letter in word:
            following = defaultdict(int)
            for state, count in current.items():
                for j in range(3):
                    if state[j] == 0:
                        if j == 1 and not (letter >> left & 1):
                            continue
                        if j == 2 and not (letter >> right & 1):
                            continue
                    new = list(state)
                    new[j] = 1 if j == 2 or not (letter >> (left if j == 0 else right) & 1) else 2
                    following[tuple(new)] += count
            current = following
        total += sum(count for state, count in current.items()
                     if 2 not in state and (state[1] or left == 0))
    return total


def tournament(word):
    n = len(word) + 2
    rows = [0] * n
    rows[0] |= 2  # marked arc x -> y
    for i, a in enumerate(word, 2):
        for j in range(2):
            if a >> j & 1:
                rows[j] |= 1 << i
            else:
                rows[i] |= 1 << j
        for j in range(i + 1, n):
            rows[i] |= 1 << j
    return rows


def subset_dp(rows):
    n = len(rows)
    # Dictionary per used set; independent of feedback gaps or chain structure.
    dp = [{} for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for last, count in dp[mask].items():
            bits = rows[last] & ~mask
            while bits:
                bit = bits & -bits
                bits -= bit
                j = bit.bit_length() - 1
                d = dp[mask | bit]
                d[j] = d.get(j, 0) + count
    return sum(dp[-1].values())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs='?', default='two_feedback_100000.json')
    args = parser.parse_args()
    data = json.loads((HERE / args.file).read_text())
    seen = set()
    for h, length, number in data['witnesses']:
        require(h not in seen and h % 2 == 1 and h <= data['limit'], 'target schema')
        seen.add(h)
        require(gaps(decode(length, number)) == h, f'wrong witness {h}')
    expected = set(range(1, data['limit'] + 1, 2)) - {7, 21}
    require(seen == expected, 'coverage claim failed')
    controls = 0
    for length in range(5):
        for word in product(range(4), repeat=length):
            rows = tournament(word)
            expected = gaps(word)
            require(subset_dp(rows) == expected, 'subset-DP disagreement')
            if length <= 3:
                literal = sum(all(rows[a] >> b & 1 for a, b in zip(p, p[1:]))
                              for p in permutations(range(len(rows))))
                require(literal == expected, 'permutation disagreement')
            controls += 1
    # The new literature-prefix target receives a different exact algorithm too.
    special = next((w for w in data['witnesses'] if w[0] == 80557), None)
    if special:
        h, length, number = special
        require(subset_dp(tournament(decode(length, number))) == h, '80557 audit')
    require(gaps((1, 2, 0)) == 13, 'one-apex hostile lost')
    print(json.dumps({'status': 'passed', 'positive_witnesses': len(seen),
                      'covered_through': data['limit'] - (data['limit'] % 2 == 0),
                      'excluded_values_from_inherited_theorem': [7, 21],
                      'independent_small_controls': controls,
                      '80557_subset_dp_checked': special is not None,
                      'scope': 'All saved positive witnesses, not global completeness.'},
                     sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
