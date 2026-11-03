#!/usr/bin/env python3
"""No producer imports: permutation counts, rational rank, native port replay."""
import itertools
import json
from fractions import Fraction
from math import isqrt
from pathlib import Path


def check(condition, message):
    if not condition:
        raise ValueError(message)


def count(rows):
    return sum(all((rows[a] >> b) & 1 for a, b in zip(p, p[1:]))
               for p in itertools.permutations(range(len(rows))))


def extend(rows, incoming):
    n = len(rows)
    a = list(rows) + [0]
    for v in range(n):
        if v in incoming:
            a[v] += 2**n
        else:
            a[n] += 2**v
    return a


def run_commands(rows, ports, word):
    a, p = list(rows), list(ports)
    for cmd in word:
        mask, discard = divmod(cmd, 3)
        n = len(a)
        a = extend(a, {p[i] for i in range(2) if mask & 2**i})
        if discard < 2:
            p[discard] = n
    return a, p


def row_for_code(n, code):
    a = [0] * n
    for bit, (i, j) in enumerate(itertools.combinations(range(n), 2)):
        if code & 2**bit:
            a[i] += 2**j
        else:
            a[j] += 2**i
    h = count(a)
    singles = [count(extend(a, {i})) for i in range(n)]
    linear = [x - h for x in singles]
    pairs = [singles[i] + singles[j] - h - count(extend(a, {i, j}))
             for i, j in itertools.combinations(range(n), 2)]
    return [h] + linear + pairs


def rational_rank(matrix):
    a = [[Fraction(v) for v in row] for row in matrix]
    rank = 0
    for j in range(len(a[0])):
        i = next((i for i in range(rank, len(a)) if a[i][j]), None)
        if i is None:
            continue
        a[rank], a[i] = a[i], a[rank]
        value = a[rank][j]
        a[rank] = [v / value for v in a[rank]]
        for i in range(rank + 1, len(a)):
            value = a[i][j]
            a[i] = [x - value * y for x, y in zip(a[i], a[rank])]
        rank += 1
    return rank


def main():
    data = json.loads(Path(__file__).with_name("receipt.json").read_text())
    report = {"audit": "same agent; independent implementation; no producer imports"}
    ranks = []
    for row in data["exhaustive_cut_census"]:
        n = row["order"]
        matrix = [row_for_code(n, c) for c in row["basis_source_codes"]]
        r = rational_rank(matrix)
        check(r == row["linear_rank_mod_1000003"], "rank disagreement")
        ranks.append([n, r])
    report["exact_rational_rank_of_supplied_minors"] = ranks

    source = data["H31_shared_source"]
    b, ports = run_commands((2, 4, 1), (0, 1), (4, 8, 3, 5))
    check(b == source["rows"] and ports == [5, 3], "source command mismatch")
    old = [int(x) for x in source["shared_labels"]]
    words = {93: (3, 3), 95: (10,), 97: (6,)}
    targets = []
    for witness in source["witnesses"]:
        labels = [int(x) for x in witness["labels"]]
        rows, terminal = run_commands(b, ports, words[witness["target"]])
        check(rows == witness["rows"] and terminal == witness["ports"], "branch or ports")
        check(labels[:7] == old and len(set(labels)) == len(labels), "shared labels")
        for i, j in itertools.combinations(range(len(labels)), 2):
            z = labels[i] + labels[j]
            check(bool(rows[i] & 2**j) == (isqrt(z)**2 == z), "square arc mismatch")
        h = count(rows)
        check(h == witness["target"], "direct path count mismatch")
        reach = [r | 2**i for i, r in enumerate(rows)]
        for k in range(len(rows)):
            for i in range(len(rows)):
                if reach[i] & 2**k:
                    reach[i] |= reach[k]
        check(all(r == 2**len(rows) - 1 for r in reach), "branch not strong")
        targets.append(h)
    report["permutation_recounted_strong_targets"] = targets
    report["source_H_by_permutations"] = count(b)

    roots = {}
    for key, interval in data["MRB_orientation_separator"]["root_intervals"].items():
        n = int(key)
        lo, hi = map(Fraction, interval)
        check(lo**n <= n < hi**n, "invalid rational root bound")
        roots[n] = lo, hi
    upper = 2 * sum(roots[k + 1][1] - roots[k][0] for k in (2, 6, 10))
    check(upper == Fraction(data["MRB_orientation_separator"]["beta_minus_alpha_upper_bound"]),
          "root sum mismatch")
    check(22 * upper < Fraction(-319, 500), "source orientation sign")
    report["MRB_source_phase_gap_strict_upper_bound"] = str(22 * upper)
    # Independent four-character arithmetic, keeping the declared query clock.
    responses = [count(extend(b, s)) for s in (set(), {5}, {3}, {5, 3})]
    chars = ([1, 1, 1, 1], [-1, 1, -1, 1], [-1, -1, 1, 1], [1, -1, -1, 1])
    coefficients = [sum(a*b for a, b in zip(responses, c)) / 4 for c in chars]
    check(coefficients == [69, 5, 27, -6], "source Walsh coefficients")
    report["four_cube_coefficients"] = [int(v) for v in coefficients]
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
