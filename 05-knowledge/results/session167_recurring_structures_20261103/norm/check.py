#!/usr/bin/env python3
"""Exact certificates for the tribonacci guard norm shell; no assert statements."""
from fractions import Fraction
from functools import reduce
from math import gcd
import json
from pathlib import Path


def require(test, message):
    if not test:
        raise RuntimeError(message)


def C(v):
    s, u, d = v
    return u, s + u + d, s


def inverse_C(v):
    s, u, d = v
    return d, s, u - s - d


def norm(v):
    s, u, d = v
    return (d**3 + 2*d*d*s + d*d*u + 2*d*s*s - 2*d*s*u
            - d*u*u + 2*s**3 - 2*s*u*u + u**3)


def multiply(a, b):
    """Product mod x^k-(1+x+...+x^(k-1)), ascending coefficients."""
    k = len(a)
    c = [0]*(2*k-1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            c[i+j] += ai*bj
    for i in range(2*k-2, k-1, -1):
        for j in range(k):
            c[i-k+j] += c[i]
    return c[:k]


def determinant(a):
    n = len(a)
    b = [[Fraction(x) for x in row] for row in a]
    answer = Fraction(1)
    for j in range(n):
        pivot = next((i for i in range(j, n) if b[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            b[j], b[pivot] = b[pivot], b[j]
            answer = -answer
        t = b[j][j]
        answer *= t
        for i in range(j+1, n):
            factor = b[i][j]/t
            for ell in range(j, n):
                b[i][ell] -= factor*b[j][ell]
    require(answer.denominator == 1, "nonintegral determinant")
    return answer.numerator


def algebra_norm(a):
    k = len(a)
    columns = []
    for i in range(k):
        e = [0]*k
        e[i] = 1
        columns.append(multiply(a, e))
    return determinant([[columns[j][i] for j in range(k)] for i in range(k)])


def root_catalogue(K):
    B = 0
    while (B+1)**3 <= 1792*K:
        B += 1
    return [v for s in range(B+1) for u in range(B-s+1)
            for d in range(B-s-u+1) if u < s+d
            and norm(v := (s, u, d)) == K]


def path_start_counts(out):
    """Independent Hamilton-path subset DP; result indexed by first vertex."""
    n = len(out)
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        if mask & (mask-1) == 0:
            continue
        vertices = mask
        while vertices:
            bit = vertices & -vertices
            v = bit.bit_length()-1
            rest = mask ^ bit
            nxt = out[v] & rest
            value = 0
            while nxt:
                low = nxt & -nxt
                value += dp[rest][low.bit_length()-1]
                nxt ^= low
            dp[mask][v] = value
            vertices ^= bit
    return dp[-1]


def full_tiling(n):
    out = [0]*n
    for i in range(n):
        for j in range(i+1, n):
            a, b = (i, j) if j == i+1 else (j, i)
            out[a] |= 1 << b
    return out


def matmul(a, b, m):
    return [[sum(a[i][k]*b[k][j] for k in range(3)) % m
             for j in range(3)] for i in range(3)]


def matpow(a, n, m):
    b = [[int(i == j) for j in range(3)] for i in range(3)]
    while n:
        if n & 1:
            b = matmul(b, a, m)
        a = matmul(a, a, m)
        n //= 2
    return b


def prime_factors(n):
    factors = []
    d = 2
    while d*d <= n:
        if n % d == 0:
            factors.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors


def main():
    # An independent determinant implementation verifies every coefficient of
    # the cubic identity on a degree-determining 4x4x4 grid.
    for s in range(4):
        for u in range(4):
            for d in range(4):
                v = (s, u, d)
                require(norm(v) == algebra_norm([d, u-s, s]), "norm identity")
                require(norm(C(v)) == norm(v), "C norm identity")
                require(inverse_C(C(v)) == v, "inverse identity")
    # A degree <=3 polynomial in each variable is determined by these grids.
    f = lambda x: x**3-x*x-x-1
    require(f(Fraction(11, 6)) < 0 < f(Fraction(15, 8)), "root bracket")
    require(Fraction(8, 15)-Fraction(49, 256) > Fraction(1, 3), "imaginary bound")
    require(15**3 < 1792*2 < 16**3, "norm-two size bound")
    all_small = [(s, u, d) for s in range(16) for u in range(16-s)
                 for d in range(16-s-u)]
    boundary = [v for v in all_small if v[1] < v[0]+v[2]]
    require(len(all_small) == 816 and len(boundary) == 576, "certificate size")
    roots = [v for v in boundary if norm(v) == 2]
    require(roots == [(1, 0, 0)], "norm-two root catalogue")
    minima = [min(norm(v) for v in boundary if sum(v) == ell) for ell in range(1, 16)]
    require(minima == [1, 7, 4, 19, 11, 32, 17, 53, 26, 73, 44, 98, 53, 136, 76],
            "strict-cone minima certificate")
    require(root_catalogue(1) == [(0, 0, 1)], "norm-one root catalogue")
    small_catalogues = {K: root_catalogue(K) for K in range(1, 11)}
    expected = {1: [(0, 0, 1)], 2: [(1, 0, 0)], 3: [],
                4: [(1, 1, 1)], 5: [], 6: [], 7: [(1, 0, 1)],
                8: [(0, 0, 2)], 9: [], 10: []}
    require(small_catalogues == expected, "small norm catalogues")
    v = (1, 0, 0)
    for _ in range(100):
        require(norm(v) == 2 and reduce(gcd, v) == 1, "orbit invariant")
        if v != (1, 0, 0):
            require(sum(inverse_C(v)) < sum(v), "descent")
        v = C(v)
    physical = []
    v = C((1, 0, 0))
    for n in range(2, 14):
        out = full_tiling(n)
        starts = path_start_counts(out)
        smaller = path_start_counts(full_tiling(n-1))
        actual = (starts[-1], sum(starts)-starts[-1], smaller[-1])
        require(actual == v, f"physical endpoint state n={n}")
        physical.append({"n": n, "state": v, "H": sum(starts)})
        v = C(v)
    # Exact modular periods and first p-adic lift witnesses.
    matrix = [[0, 1, 0], [1, 1, 1], [1, 0, 0]]
    ident = [[int(i == j) for j in range(3)] for i in range(3)]
    periods = {}
    for p, order in [(3, 13), (5, 31), (47, 46)]:
        require(matpow(matrix, order, p) == ident, "period upper bound")
        for q in prime_factors(order):
            require(matpow(matrix, order//q, p) != ident, "period minimality")
        pp = matpow(matrix, order, p*p)
        witness = [[((pp[i][j]-ident[i][j])//p) % p for j in range(3)]
                   for i in range(3)]
        require(any(any(row) for row in witness), "nonzero lift witness")
        for k in range(1, 4):
            m, t = p**k, order*p**(k-1)
            require(matpow(matrix, t, m) == ident, "lift period")
            for q in prime_factors(t):
                require(matpow(matrix, t//q, m) != ident, "lift minimality")
        shell = sum(norm((s, u, d)) % p == 2 for s in range(p)
                    for u in range(p) for d in range(p))
        expected_shell = {3: 13, 5: 31, 47: 2116}[p]
        require(shell == expected_shell, "modular norm-shell cardinality")
        periods[p] = {"period": order, "lift_witness": witness,
                      "norm_two_points": shell, "C_orbits": shell//order}
    require([x for x in range(47) if (x**3-x*x-x-1) % 47 == 0] == [5, 17, 26],
            "split prime roots")
    require(determinant([[1, 1, 0], [1, 2, 1], [3, 3, 2]]) == 2,
            "scalar observation invertible modulo odd moduli")
    # Sharp boundary: higher-degree norm values do not determine one orbit.
    counterexamples = {
        4: {"first": [1, 1, 2, 2], "second": [1, 2, 2, 1], "norm": 3,
            "unit": [1, 1, 0, 0], "inverse": [2, -1, 2, -1]},
        5: {"first": [1, 1, 1, 1, 2], "second": [1, 1, 1, 2, 1], "norm": 2,
            "unit": [1, 0, 1, 1, 0], "inverse": [1, 3, 4, 1, -2]},
    }
    for k, c in counterexamples.items():
        require(algebra_norm(c["first"]) == c["norm"], "first higher norm")
        require(algebra_norm(c["second"]) == c["norm"], "second higher norm")
        require(algebra_norm(c["unit"]) == 1, "higher unit norm")
        base = [1, 1, 1, 0] if k == 4 else [1, 1, 0, 0, 0]
        require(algebra_norm(base) == c["norm"], "higher base norm")
        require(multiply(base, c["unit"]) == c["second"], "second orbit expression")
        x = [0, 1]+[0]*(k-2)
        first = base
        for _ in range(2 if k == 4 else 4):
            first = multiply(first, x)
        require(first == c["first"], "first orbit expression")
        require(multiply(c["unit"], c["inverse"]) == [1]+[0]*(k-1), "unit inverse")
        fk = lambda x: x**k-sum(x**j for j in range(k))
        require(fk(Fraction(19, 10)) < 0 < fk(2), "higher root bracket")
    summary = {"status": "PASS", "norm_identity_grid": 64,
               "complete_boundary_candidates": len(boundary),
               "norm_two_roots": roots, "small_norm_roots": small_catalogues,
               "physical_DP": physical, "prime_periods": periods,
               "higher_degree_counterexamples": counterexamples}
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
