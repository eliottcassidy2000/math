#!/usr/bin/env python3
"""Independent standard-library audit for proposed THM-4362.

This script starts from the displayed generator law in THM-4308(18)--(19).
It deliberately does not import any repository computation.
"""

from __future__ import annotations

import hashlib
import math
import pathlib
import sys


sys.stdout.reconfigure(newline="\n")


def require(condition: bool, detail=None) -> None:
    """Optimization-stable assertion: Python -O must run every check."""

    if not condition:
        raise RuntimeError(f"audit check failed: {detail!r}")


def C(n: int, k: int) -> int:
    """Ordinary combinatorial binomial coefficient (zero outside 0<=k<=n)."""

    return math.comb(n, k) if 0 <= k <= n else 0


def generator_terms(a: int, b: int, c: int, e: int, m: int):
    """Terms of x^(a+2b+e)t^(b+c+2e)(1+x^2t)^(c+e), projected to n<=m."""

    r0 = a + 2 * b + e
    n0 = b + c + 2 * e
    N = c + e
    return [
        (n0 + k, r0 + 2 * k, C(N, k))
        for k in range(N + 1)
        if n0 + k <= m
    ]


def weights(m: int):
    return [((-1) ** (n - 5)) * C(m + 3 - n, 3) for n in range(5, m + 1)]


def L_value(m: int, a: int, b: int, c: int, e: int) -> int:
    total = 0
    for n, r, coeff in generator_terms(a, b, c, e, m):
        if 5 <= n <= m and r == 2 * n - 10:
            total += ((-1) ** (n - 5)) * C(m + 3 - n, 3) * coeff
    return total


def hostile(m: int):
    d = m - 6
    eps = d % 2
    a = d
    b = 0
    e = eps
    c = (10 + d - 3 * eps) // 2
    return d, a, b, c, e


def generalized_L_value(
    m: int, ell: int, q: int, a: int, b: int, c: int, e: int
) -> int:
    """Diagonal-offset extension discovered during the clean-room audit."""

    nstar = (ell + 1) // 2
    total = 0
    for n, r, coeff in generator_terms(a, b, c, e, m):
        if nstar <= n <= m and r == 2 * n - ell:
            total += ((-1) ** (n - nstar)) * C(m + q - n, q) * coeff
    return total


def projected_generators(d: int, m: int):
    """Every (18) generator with a nonzero row-<=m projection."""

    for a in range(d + 1):
        for b in range(d - a + 1):
            # A projected term exists exactly when n0=b+c+2e<=m.
            for e in range((m - b) // 2 + 1):
                for c in range(m - b - 2 * e + 1):
                    yield a, b, c, e


checks = 0
columns = 0
annihilating_pairs = 0
generalized_columns = 0

# Exact stencils used by the row-nine P2 and proposed row-ten P3 arguments.
require(weights(9) == [35, -20, 10, -4, 1], "m=9 stencil")
require(weights(10) == [56, -35, 20, -10, 4, -1], "m=10 stencil")
checks += 2

# Exhaust every nonzero projected generator in a broad finite grid on the
# asserted side of the threshold.  This is a check, not the universal proof.
for m in range(7, 19):
    for d in range(0, m - 6):  # exactly 0 <= d <= m-7
        annihilating_pairs += 1
        for a, b, c, e in projected_generators(d, m):
            columns += 1
            value = L_value(m, a, b, c, e)
            require(value == 0, (m, d, a, b, c, e, value))
            checks += 1

            # If the diagonal is hit, the cancellation hypotheses used in
            # the proof really hold, including both truncation endpoints.
            support = a == 2 * c + 3 * e - 10
            hits = any(r == 2 * n - 10 and 5 <= n <= m
                       for n, r, _ in generator_terms(a, b, c, e, m))
            require(hits == support, (m, d, a, b, c, e, hits, support))
            checks += 1
            if support:
                N = c + e
                n0 = b + c + 2 * e
                nmax = n0 + N
                require(N >= 4, ("N", m, d, a, b, c, e, N))
                require(5 <= n0 <= m - 1, ("n0", m, n0))
                require(nmax == b + a + 10 <= m + 3, ("nmax", m, d, nmax))
                R = m + 3 - n0
                alternating = sum(
                    ((-1) ** k) * C(N, k) * C(R - k, 3)
                    for k in range(N + 1)
                )
                require(alternating == 0, ("finite difference", m, d, a, b, c, e))
                checks += 6

# Verify the exact boundary witness, its sign, and containment in every later
# filtration level.  We go far beyond the two stencils used downstream.
for m in range(7, 101):
    d0, a, b, c, e = hostile(m)
    N = c + e
    n0 = b + c + 2 * e
    nmax = n0 + N
    require(a >= 0 and b == 0 and c >= 0 and e in (0, 1), ("hostile indices", m))
    require(a + b == d0, ("hostile depth", m))
    require(a == 2 * c + 3 * e - 10, ("hostile support", m))
    require(N >= 5 and n0 >= 5 and nmax == m + 4, ("hostile endpoints", m))
    value = L_value(m, a, b, c, e)
    require(value == (-1) ** (d0 + 7) == (-1) ** (m - 1), ("hostile sign", m, value))
    checks += 7
    for d in range(d0, d0 + 13):
        require(a + b <= d and value != 0, ("nestedness", m, d))
        checks += 2

# Independently verify the partial alternating sum that controls the hostile
# sign, using only ordinary (nonnegative-top) binomial coefficients.
for N in range(5, 101):
    partial = sum(
        ((-1) ** k) * C(N, k) * C(N - 1 - k, 3)
        for k in range(N - 3)
    )
    require(partial == (-1) ** N, ("partial identity", N, partial))
    checks += 1

# Figurate-number claims: C(r,3) is the cumulative triangular sequence and a
# cubic, hence has identically zero fourth forward difference.
T = lambda z: z * (z + 1) // 2
for r in range(3, 101):
    require(C(r, 3) == sum(T(q) for q in range(1, r - 1)), ("tetrahedral", r))
    checks += 1
for r in range(3, 97):
    vals = [C(r + j, 3) for j in range(5)]
    fourth = vals[4] - 4 * vals[3] + 6 * vals[2] - 4 * vals[1] + vals[0]
    require(fourth == 0, ("fourth difference", r, fourth))
    checks += 1
for z in range(-100, 101):
    require(T(z + 2) - T(z - 2) == 2 + 4 * z, ("centered T", z))
    checks += 1

# Safe extension: the same proof is a diagonal-offset/order hierarchy.  For
# ell>=1, n*=ceil(ell/2), and every 0<=q<ceil(ell/3), the functional with
# weights (-1)^(n-n*) C(m+q-n,q) has sharp threshold d<=m+q-ell whenever
# m>=ell-q.  The maximal q gives the deepest threshold.
for ell in range(1, 19):
    nstar = (ell + 1) // 2
    for q in range((ell - 1) // 3 + 1):
        for m in range(ell - q, ell - q + 9):
            threshold = m + q - ell
            for d in range(threshold + 1):
                for a, b, c, e in projected_generators(d, m):
                    generalized_columns += 1
                    value = generalized_L_value(m, ell, q, a, b, c, e)
                    require(
                        value == 0,
                        ("generalized zero", ell, q, m, d, a, b, c, e, value),
                    )
                    checks += 1

            d0 = threshold + 1
            S = ell + d0  # = m+q+1
            eps = S % 2
            a, b, e, c = d0, 0, eps, (S - 3 * eps) // 2
            N = c + e
            n0 = c + 2 * e
            require(c >= 0 and a == 2 * c + 3 * e - ell, ("generalized support", ell, q, m))
            require(N >= q + 1 and n0 >= nstar, ("generalized minima", ell, q, m))
            require(n0 + N == m + q + 1, ("generalized endpoint", ell, q, m))
            value = generalized_L_value(m, ell, q, a, b, c, e)
            expected = (-1) ** (n0 - nstar + N + q + 1)
            require(
                value == expected and abs(value) == 1,
                ("generalized hostile", ell, q, m, value),
            )
            checks += 4

script_hash = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
print("THM-4362 CLEAN-ROOM REFEREE: PASS")
print(f"annihilating (m,d) pairs exhaustively checked: {annihilating_pairs}")
print(f"projected generator columns checked: {columns}")
print("boundary witnesses checked: m=7..100, plus 13 nested depths each")
print("partial hostile-sign identities checked: N=5..100")
print("stencil m=9:  35,-20,10,-4,1")
print("stencil m=10: 56,-35,20,-10,4,-1")
print("figurate checks: cumulative triangular, fourth difference, centered T identity")
print(f"general diagonal/order hierarchy columns checked: {generalized_columns}")
print(f"total assertions grouped as checks: {checks}")
print(f"script_sha256={script_hash}")
