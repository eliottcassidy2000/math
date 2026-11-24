"""Recover all line partners at small source heights and test Euclid targets.

Uses literal anchored boards only for the initial three entries on each axis.
All remaining singleton lengths run through the exact tribonacci recurrence.
"""
from collections import defaultdict
from itertools import combinations, product
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "chanpak_ternary_buckets_20261122"))
from target_flip import height

BOUND = 23
MAX_EXPONENT = 61
PERFECTS = {p: (1 << (p-1))*((1 << p)-1) for p in range(3, MAX_EXPONENT+1)}
TARGETS = {value+sign for value in PERFECTS.values() for sign in (-1, 1)}
BIG = PERFECTS[MAX_EXPONENT]
CAP = 2*BIG
# lambda > 9/5. Every tensor of d positive seeds on each ternary grid
# has F >= lambda^(sum(r)-2d), and d <= 3 for the source bank here.
MAX_GAPS = next(k for k in range(6, 300) if 9**(k-6) > CAP*5**(k-6))


def board(n, reversed_arcs, partner):
    rows = [0]*n
    for i, j in combinations(range(n), 2):
        backwards = ((i, j) in reversed_arcs) ^ (partner and j-i >= 2)
        u, v = (j, i) if backwards else (i, j)
        rows[u] |= 1 << v
    return rows


def panels():
    result = defaultdict(list)
    for n in range(3, 7):
        tiles = [(i, j) for i, j in combinations(range(n), 2) if j-i >= 2]
        for mask in range(1 << len(tiles)):
            arcs = tuple(t for i, t in enumerate(tiles) if mask >> i & 1)
            # A forward anchored board is strong iff every cut is crossed backwards.
            if any(not any(a <= cut < b for a, b in arcs) for cut in range(n-1)):
                continue
            h = height(board(n, arcs, False))
            if h <= BOUND:
                result[h].append((n, arcs))
    return result


def all_words(catalogue, remaining):
    yield ()
    for h, types in catalogue.items():
        if h <= remaining:
            for panel in types:
                for suffix in all_words(catalogue, remaining//h):
                    yield ((h, panel),) + suffix


def partner_height(word, gaps):
    n = gaps[0]
    reversed_arcs = set()
    for (_, (size, arcs)), gap in zip(word, gaps[1:]):
        reversed_arcs.update((n+i, n+j) for i, j in arcs)
        n += size+gap
    return height(board(n, reversed_arcs, True)) if n else 1


def contract(values, weights):
    span = len(values)//3
    a, b, c = weights
    return tuple(a*values[i] + b*values[i+span] + c*values[i+2*span]
                 for i in range(span))


def scan(word, weights, visit):
    d = 1+len(word)
    table = tuple(partner_height(word, gaps) for gaps in product(range(3), repeat=d))

    def extend(values, gaps, budget):
        if len(values) == 3:
            a, b, c = values
            for r in range(budget+1):
                if r >= 3 and a > CAP:
                    break
                if word or r >= 1:
                    visit(a, gaps+(r,))
                a, b, c = b, c, a+b+c
            return
        for r in range(budget+1):
            extend(contract(values, weights[r]), gaps+(r,), budget-r)
    extend(table, (), MAX_GAPS-1)


def triangle_checks(weights):
    word = ((3, (3, ((0, 2),))),)
    table = tuple(partner_height(word, gaps) for gaps in product(range(3), repeat=2))
    assert table == (1,5,9,5,13,23,9,23,41)

    def f(r, s):
        return sum(a*b for a, b in zip(contract(table, weights[r]), weights[s]))

    for r, s in ((0, 5), (1, 4), (3, 3), (3, 5)):
        assert f(r, s) == partner_height(word, (r, s))
    assert f(8,14) == f(11,11) == 2295169
    for pos in (8, 11):
        scores = [row.bit_count() for row in board(25, {(pos, pos+2)}, True)]
        assert scores.count(pos+1) == 3
    blue = [f(i, i) for i in range(20)]
    assert blue[:6] == [1,13,41,129,457,1533]
    for i in range(6, len(blue)):
        assert blue[i] == 2*blue[i-1]+3*blue[i-2]+6*blue[i-3]-blue[i-4]-blue[i-6]


def main():
    catalogue = panels()
    counts = {h: len(types) for h, types in sorted(catalogue.items())}
    print("strong panels by source height", counts)
    assert counts == {3: 1, 5: 5, 9: 18, 11: 11, 13: 13, 15: 38, 17: 34, 19: 38, 23: 92}
    weights = [(1,0,0), (0,1,0), (0,0,1)]
    while len(weights) < MAX_GAPS:
        weights.append(tuple(sum(axis) for axis in zip(*weights[-3:])))
    triangle_checks(weights)
    hits = defaultdict(set)
    lower = {1: None}
    upper = {1: None}
    for word in all_words(catalogue, BOUND):
        h = 1
        for factor, _ in word:
            h *= factor
        def visit(value, gaps):
            if value in TARGETS:
                hits[value].add(h)
            if value <= BIG and (h not in lower or lower[h] is None or value > lower[h][0]):
                lower[h] = value, word, gaps
            if value >= BIG and (h not in upper or upper[h] is None or value < upper[h][0]):
                upper[h] = value, word, gaps
        scan(word, weights, visit)
    for p, e in PERFECTS.items():
        if hits[e-1] or hits[e+1]:
            print("exponent", p, "lower sources", sorted(hits[e-1]), "upper sources", sorted(hits[e+1]))
    lo = max(lower.items(), key=lambda item: item[1][0])
    hi = min(upper.items(), key=lambda item: item[1][0])
    assert BIG//2 < lo[1][0] < BIG < hi[1][0] < CAP



    assert {p for p, e in PERFECTS.items() if hits[e-1] or hits[e+1]} == {3,4,5,7}
    assert lo[1][0]-BIG == -2213982461498407132682959208946669
    assert hi[1][0]-BIG == 1485690430693459244528868155770979
    print("nearest below", lo, "offset", lo[1][0]-BIG)
    print("nearest above", hi, "offset", hi[1][0]-BIG)
    print("all source heights <=", BOUND, "all admissible orders, all gap sums below", MAX_GAPS, "for outputs <=", CAP)


if __name__ == "__main__":
    main()
