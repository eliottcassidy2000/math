#!/usr/bin/env python3
"""Exact constructive controls for marked root-character Laurent lifts.

Standard library only. No assertions: all gates survive python -O.
"""
import sys as _output_sys
_output_sys.stdout.reconfigure(newline="\n")
from collections import Counter, deque
from itertools import product
from math import gcd
import json
import random
import sys


def require(ok, msg):
    if not ok:
        raise ValueError(msg)


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def matmul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def det(a):
    """Bareiss determinant, exact over Z."""
    a = [row[:] for row in a]
    n, sign, prev = len(a), 1, 1
    for k in range(n - 1):
        if a[k][k] == 0:
            r = next((r for r in range(k + 1, n) if a[r][k]), None)
            if r is None:
                return 0
            a[k], a[r] = a[r], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                num = pivot * a[i][j] - a[i][k] * a[k][j]
                require(num % prev == 0, "Bareiss nonexact division")
                a[i][j] = num // prev
            a[i][k] = 0
        prev = pivot
    return sign * a[-1][-1]


def primes(m):
    return [p for p in range(2, m + 1) if m % p == 0 and all(p % d for d in range(2, p))]


def lift(c, m):
    """Return integral unimodular Q=C mod m, or None for determinant obstruction.

    Modular elementary row/column elimination records exact inverse operations.
    Local prime choices plus CRT create a unit pivot without a field assumption.
    """
    n = len(c)
    delta = det(c) % m
    if delta not in {1, m - 1}:
        return None
    a = [[v % m for v in row] for row in c]
    ri, si = identity(n), identity(n)

    def rowadd(i, j, t):
        a[i] = [(x + t * y) % m for x, y in zip(a[i], a[j])]
        # R -> E_ij(t) R, so R^-1 -> R^-1 E_ij(-t).
        for row in ri:
            row[j] -= t * row[i]

    def coladd(i, j, t):
        # Add t times column i to column j.
        for row in a:
            row[j] = (row[j] + t * row[i]) % m
        si[i] = [x - t * y for x, y in zip(si[i], si[j])]

    for k in range(n - 1):
        # For each prime, if a_kk vanishes, add one nonvanishing row.
        local = {p: [0] * (n - k - 1) for p in primes(m)}
        for p, choices in local.items():
            if a[k][k] % p == 0:
                r = next(r for r in range(k + 1, n) if a[r][k] % p)
                choices[r - k - 1] = 1
        for r in range(k + 1, n):
            t = next(t for t in range(m) if all(t % p == choices[r - k - 1] for p, choices in local.items()))
            if t:
                rowadd(k, r, t)
        v = a[k][k]
        require(gcd(v, m) == 1, "pivot not a unit")
        for r in range(k + 1, n):
            rowadd(r, k, -a[r][k] * pow(v, -1, m) % m)

        # D(u)=E12(u) E21(-u^-1) E12(u) E12(-1) E21(1) E12(-1).
        # Apply reversed order, since these multiply on the left.
        u, ui = pow(v, -1, m), v
        ops = [(k, k + 1, u), (k + 1, k, -ui), (k, k + 1, u),
               (k, k + 1, -1), (k + 1, k, 1), (k, k + 1, -1)]
        for i, j, t in reversed(ops):
            rowadd(i, j, t)
        require(a[k][k] == 1, "diagonal rescale failed")
        for j in range(k + 1, n):
            coladd(k, j, -a[k][j])
    require(a[-1][-1] in {1, m - 1}, "last determinant mismatch")
    if a[-1][-1] != 1:
        a[-1] = [-x % m for x in a[-1]]
        for row in ri:
            row[-1] *= -1
    require(a == identity(n), "elimination incomplete")
    q = matmul(ri, si)
    require(abs(det(q)) == 1, "lift not unimodular")
    require(all((q[i][j] - c[i][j]) % m == 0 for i in range(n) for j in range(n)), "lift changes a marked character")
    return q


def reachable_rank2(m):
    """Independent BFS using only integral shear/sign generators modulo m."""
    start = (1, 0, 0, 1)
    seen, todo = {start}, deque([start])
    while todo:
        a, b, c, d = todo.popleft()
        for z in (((a + c) % m, (b + d) % m, c, d),
                  (a, b, (a + c) % m, (b + d) % m),
                  (-a % m, -b % m, c, d)):
            if z not in seen:
                seen.add(z)
                todo.append(z)
    return seen


def main():
    receipt = {"status": "FINITE-EXACT controls; all-rank theorem has a separate proof", "rank2": []}
    total_lifts = 0
    for m in range(2, 13):
        seen = reachable_rank2(m)
        counts = Counter()
        lifted = 0
        for z in product(range(m), repeat=4):
            c = [list(z[:2]), list(z[2:])]
            delta = (z[0] * z[3] - z[1] * z[2]) % m
            if gcd(delta, m) != 1:
                continue
            counts[delta] += 1
            q = lift(c, m)
            require((q is not None) == (z in seen), "BFS disagrees with constructive lift")
            lifted += q is not None
        require(lifted == len(seen), "BFS includes noninvertible matrix")
        receipt["rank2"].append({"m": m, "invertible": sum(counts.values()), "liftable": lifted, "determinant_fibres": dict(sorted(counts.items()))})
        total_lifts += lifted

    rank3 = []
    for m in (2, 3):
        invertible = lifted = 0
        for z in product(range(m), repeat=9):
            c = [list(z[i:i+3]) for i in (0, 3, 6)]
            if gcd(det(c), m) != 1:
                continue
            invertible += 1
            require(lift(c, m) is not None, "universal small modulus failed")
            lifted += 1
        rank3.append({"m": m, "invertible": invertible, "lifted": lifted})
        total_lifts += lifted
    receipt["rank3_exhaustive"] = rank3

    rng = random.Random(261026)
    accepted = augmented = 0
    for n in range(1, 7):
        for m in range(2, 31):
            for _ in range(8):
                while True:
                    c = [[rng.randrange(m) for _ in range(n)] for _ in range(n)]
                    delta = det(c) % m
                    if gcd(delta, m) == 1:
                        break
                accepted += 1
                q = lift(c, m)
                require((q is not None) == (delta in {1, m - 1}), "determinant gate failed")
                if q is None:
                    d = [row + [0] for row in c] + [[0] * n + [pow(delta, -1, m)]]
                    require(lift(d, m) is not None, "one compensating character failed")
                    augmented += 1
    receipt["seeded_constructive_controls"] = {"seed": 261026, "ranks": [1, 2, 3, 4, 5, 6], "moduli": [2, 30], "invertible_inputs": accepted, "one_character_repairs": augmented}
    receipt["exhaustive_successful_lifts"] = total_lifts

    # Explicit smallest-modulus and six/seven controls; these are symbolic matrices.
    q5, q7 = [[2, 5], [5, 13]], [[2, 7], [7, 25]]
    for m, q in ((5, q5), (7, q7)):
        require(det(q) == 1, "displayed stabilized matrix determinant")
        require(q[0][0] % m == 2 and q[1][1] % m == pow(2, -1, m), "displayed characters")
        require(q[0][1] % m == q[1][0] % m == 0, "displayed cross characters")
        require(lift([[2]], m) is None, "scalar hostile disappeared")
    require(lift([[5]], 6) == [[-1]], "sixth-root inversion control")
    receipt["explicit_controls"] = {"m5_one_auxiliary_lift": q5, "m7_one_auxiliary_lift": q7, "m6_exponent5_lift": [[-1]], "affine_hostile": {"matrix": [[1, 1], [0, 1]], "collapsed_points": [[0, 0], [1, 0]], "common_image": [0, 0]}}
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
