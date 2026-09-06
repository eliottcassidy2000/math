#!/usr/bin/env python3
"""Independent n=4 no-three-in-line audit via ordered permutation pairs.

No repository imports; unlike the producer this enumerates permutations,
then explicitly groups their uncolored unions. Collinearity uses determinants.
"""
from collections import defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations
import sys


def need(test, label):
    if not test:
        raise RuntimeError(label)


def cycles(p, q):
    inv = {value: i for i, value in enumerate(p)}
    rho = [inv[value] for value in q]
    seen = set()
    lengths = []
    for i in range(len(p)):
        if i in seen:
            continue
        j, length = i, 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = rho[j]
        lengths.append(2*length)
    return tuple(sorted(lengths))


def collinear(points):
    p, q, *rest = points
    return all((q[0]-p[0])*(r[1]-p[1]) == (q[1]-p[1])*(r[0]-p[0])
               for r in rest)


def main():
    sys.stdout.reconfigure(newline="\n")
    n = 4
    perms = list(permutations(range(n)))
    unions = {}
    pair_count = 0
    for p in perms:
        for q in perms:
            if any(a == b for a, b in zip(p, q)):
                continue
            pair_count += 1
            board = frozenset((i, p[i]) for i in range(n)) | frozenset((i, q[i]) for i in range(n))
            typ = cycles(p, q)
            old = unions.setdefault(board, [typ, 0])
            need(old[0] == typ, "cycle type independent of coloring")
            old[1] += 1
    need(pair_count == 216 and len(unions) == 90, "complete permutation/uncolored universes")
    groups = defaultdict(list)
    zero_boards = zero_pairs = 0
    for board, (typ, weight) in unions.items():
        need(weight == 2**len(typ), "exact ordered-decomposition multiplicity")
        x3 = sum(collinear(t) for t in combinations(sorted(board), 3))
        x4 = sum(collinear(t) for t in combinations(sorted(board), 4))
        groups[typ].append((x3, x4))
        zero_boards += x3 == 0
        zero_pairs += weight*(x3 == 0)
    expected = {
        (4, 4): (18, Q(2), Q(56, 9), Q(1, 3), Q(1, 2)),
        (8,): (72, Q(2), Q(25, 18), Q(1, 6), Q(1, 36)),
    }
    actual = {}
    for typ, data in sorted(groups.items()):
        count = len(data)
        mean = Q(sum(x for x, _ in data), count)
        variance = Q(sum(x*x for x, _ in data), count)-mean*mean
        fourth = Q(sum(y for _, y in data), count)
        zero = Q(sum(x == 0 for x, _ in data), count)
        actual[typ] = (count, mean, variance, fourth, zero)
    need(actual == expected, "full n4 distribution table")
    need(Q(zero_boards, len(unions)) == Q(11, 90), "uniform uncolored zero probability")
    need(Q(zero_pairs, pair_count) == Q(5, 27), "uniform ordered-pair zero probability")
    print("status=INDEPENDENT_EXACT_AUDIT_PASS; producer/repository imports=none")
    print("universe=n4 ordered disjoint permutation pairs216; distinct uncolored boards90")
    print("class_tables="+str(actual))
    print("uniform_board_zero=11/90; uniform_ordered_pair_zero=5/27")
    print("decomposition_weight=2^(number of even-cycle components), checked on every board")


if __name__ == "__main__":
    main()
