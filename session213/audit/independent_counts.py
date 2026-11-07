"""Independent Session213 exact count controls, no author imports.

Adjacency matrices, a backwards subset recurrence, and a permutation oracle
are deliberately separate from author bit-row forward implementations.
"""
import itertools
import json
from pathlib import Path


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def apex(bits):
    m = len(bits)
    a = [[False] * (m + 1) for _ in range(m + 1)]
    for i in range(m):
        for j in range(i + 1, m):
            a[i][j] = True
        a[i][m] = bool(bits[i])
        a[m][i] = not bits[i]
    return a


def count_paths(a, start=None, close=False):
    n = len(a)
    if not n:
        return 1
    counts = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        if start is None or v == start:
            counts[1 << v][v] = 1
    for mask in range(1, 1 << n):
        if not mask & (mask - 1):
            continue
        for last in range(n):
            if mask & (1 << last):
                previous = mask ^ (1 << last)
                counts[mask][last] = sum(
                    counts[previous][v] for v in range(n)
                    if (previous & (1 << v)) and a[v][last]
                )
    if close:
        check(start is not None and n >= 3, "cycle convention requires fixed root and n>=3")
        return sum(counts[-1][v] for v in range(n) if a[v][start])
    return sum(counts[-1])


def brute_paths(a, start=None):
    return sum(
        (start is None or word[0] == start)
        and all(a[u][v] for u, v in zip(word, word[1:]))
        for word in itertools.permutations(range(len(a)))
    )


def defect(bits):
    return 1 + 2 * sum(
        1 if j - i == 1 else 2 ** (j - i - 2)
        for i in range(len(bits)) for j in range(i + 1, len(bits))
        if bits[i] == 0 and bits[j] == 1
    )


def bits(q, m):
    return [(q >> i) & 1 for i in range(m)]


def main():
    receipt = {"status": "PASS", "implementation": "backwards subset DP; dense adjacency; no author imports"}
    apex_cases = sentinel_cases = permutation_cases = 0
    examples = []
    for m in range(0, 9):
        for q in range(1 << m):
            b = bits(q, m)
            d = count_paths(apex(b))
            check(d == defect(b), ("one-apex pair formula", m, q, d, defect(b)))
            check((d == 1) == ((q & (q + 1)) == 0), ("exact equality boundary", m, q, d))
            check(defect(bits(2 * q + 1, m + 1)) == d, ("odd recurrence", m, q))
            check(defect(bits(2 * q, m + 1)) == d + q + q % 2, ("even recurrence", m, q))
            check(defect(b + [0, 0]) == d, ("padding invariance", m, q))
            apex_cases += 1
            if m <= 6:
                check(brute_paths(apex(b)) == d, ("permutation oracle", m, q))
                permutation_cases += 1
            if m <= 6:
                a = apex([0, 1] + b)
                h = count_paths(a)
                root = count_paths(a, start=1)
                check(root == q + 1, ("sentinel endpoint", m, q, root))
                check(h == 2 * (q + 1) + d, ("sentinel unrestricted", m, q, h, d))
                check((h == 2 * (q + 1) + 1) == ((q & (q + 1)) == 0), ("sentinel equality", m, q))
                if m <= 3:
                    check(brute_paths(a) == h, ("sentinel permutation", m, q))
                    check(brute_paths(a, start=1) == root, ("sentinel rooted permutation", m, q))
                    permutation_cases += 2
                sentinel_cases += 1
                if q < 16 and m == 4:
                    examples.append({"q": q, "cycle_target": q + 1, "defect": d, "total_paths": h})
    receipt.update({"apex_cases_orders_1_to_9": apex_cases,
                    "sentinel_cases_orders_3_to_9": sentinel_cases,
                    "permutation_checks_orders_at_most_7": permutation_cases,
                    "sentinel_examples": examples})
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
