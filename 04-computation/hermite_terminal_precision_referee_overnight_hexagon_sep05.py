#!/usr/bin/env python3
"""Independent exact referee for terminal-cluster two-jet precision.

Universe: all translated-to-zero integer sets of 1..6 nodes, diameter <=12,
p=2,3,5,7; 1,600 seeded irregular/deep signed sets; and 44 named/seeded
literal integer Smith forms, each tested at the four small primes.

No imports from the primary compiler. The first path computes rational
reciprocal sums, the second uses only valuation-tree terminal clusters,
and the third constructs the full ordinary-value/first-Hasse-derivative
integer matrix and calls SymPy's integer Smith implementation.
"""

from fractions import Fraction as Q
from itertools import combinations
from random import Random

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


def need(condition, payload):
    if not condition:
        raise AssertionError(payload)


def vp(a, p):
    if not a:
        return 10**9
    a = abs(int(a))
    v = 0
    while a % p == 0:
        a //= p
        v += 1
    return v


def rational_v(q, p):
    return vp(q.numerator, p) - vp(q.denominator, p) if q else 10**9


def local(nodes, p):
    values = []
    for x in nodes:
        total = sum(vp(x - y, p) for y in nodes if x != y)
        reciprocal = sum((Q(1, x - y) for y in nodes if x != y), Q())
        values.append(max(2 * total, 2 * total - rational_v(2 * reciprocal, p)))
    return max(values, default=0)


def tree(nodes, p):
    if len(nodes) < 2:
        return 0
    terminal = {}
    for x in nodes:
        depth = max(vp(x - y, p) for y in nodes if x != y)
        cluster = tuple(y for y in nodes if y == x or vp(x - y, p) >= depth)
        if all(vp(y - z, p) == depth for y, z in combinations(cluster, 2)):
            total = sum(vp(x - y, p) for y in nodes if x != y)
            terminal[cluster] = 2 * total + max(0, depth - (len(cluster) == p))
    need(bool(terminal), (nodes, p, "no terminal cluster"))
    return max(terminal.values())


def literal_snf(nodes):
    n = len(nodes)
    rows = []
    for x in nodes:
        rows += [[x**j for j in range(2 * n)],
                 [0] + [j * x**(j - 1) for j in range(1, 2 * n)]]
    diagonal = smith_normal_form(Matrix(rows), domain=ZZ)
    return abs(int(diagonal[-1, -1]))


def main():
    primes = (2, 3, 5, 7)
    count = 0
    for n in range(1, 7):
        for tail in combinations(range(1, 13), n - 1):
            nodes = (0,) + tail
            for p in primes:
                lhs, rhs = local(nodes, p), tree(nodes, p)
                need(lhs == rhs, (nodes, p, lhs, rhs))
                count += 1
    rng = Random(4436)
    for index in range(1600):
        n = rng.randrange(2, 17)
        p = rng.choice((2, 3, 5, 7, 11, 13))
        nodes = tuple(sorted(rng.sample(range(-500, 501), n)))
        if index % 2 == 0:
            nodes = tuple(sorted(set(p**rng.randrange(1, 5) * x for x in nodes)))
        need(local(nodes, p) == tree(nodes, p), (nodes, p))
        count += 1
    controls = [
        (0,), (0, 1), (0, 2), (0, 8), (0, 9),
        (0, 3, 6, 1), (0, 9, 18, 3, 1), tuple(range(5)),
        tuple(5 * x for x in range(5)) + (1, 6),
        tuple(7 * x for x in range(7)) + (1, 8),
    ]
    for exponent in range(5):
        controls += [tuple(2**exponent * x for x in row)
                     for row in ((0, 1, 2, 5), (0, 1, 3, 4))]
    for _ in range(24):
        controls.append(tuple(sorted(rng.sample(range(-10, 21), rng.randrange(1, 7)))))
    snf_count = 0
    for nodes in controls:
        largest = literal_snf(nodes)
        for p in primes:
            lhs, rhs = vp(largest, p), tree(nodes, p)
            need(lhs == rhs, (nodes, p, lhs, rhs))
            snf_count += 1
    need((count, snf_count) == (7944, 176), (count, snf_count))
    print("PASS exact reciprocal/tree comparisons", count)
    print("PASS independent literal integer Smith comparisons", snf_count)
    for nodes in controls[:10]:
        print(nodes, [(p, tree(nodes, p)) for p in primes])


if __name__ == "__main__":
    main()
