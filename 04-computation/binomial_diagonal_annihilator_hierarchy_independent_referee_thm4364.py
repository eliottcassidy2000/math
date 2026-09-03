#!/usr/bin/env python3
"""Independent boundary audit for the proposed THM-4364 hierarchy.

This file uses only the source-normal monomial law printed in THM-4308.
It does not import any repository computation or scout implementation.
"""

from __future__ import annotations

import hashlib
import math
import pathlib
import sys


sys.stdout.reconfigure(newline="\n")


def require(condition: bool, detail=None) -> None:
    """Optimization-stable assertion: ``python -O`` runs every audit."""

    if not condition:
        raise RuntimeError(f"audit check failed: {detail!r}")


def C(n: int, k: int) -> int:
    """Ordinary binomial coefficient; zero unless 0 <= k <= n."""

    return math.comb(n, k) if 0 <= k <= n else 0


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def row_start(ell: int) -> int:
    return ceil_div(ell, 2)


def min_packet(ell: int) -> int:
    return ceil_div(ell, 3)


def generator_terms(a: int, b: int, c: int, e: int, m: int):
    """Projected terms of x^(a+2b+e)t^(b+c+2e)(1+x^2t)^(c+e)."""

    r0 = a + 2 * b + e
    n0 = b + c + 2 * e
    N = c + e
    return [
        (n0 + k, r0 + 2 * k, C(N, k))
        for k in range(N + 1)
        if n0 + k <= m
    ]


def L_direct(
    m: int, ell: int, q: int, a: int, b: int, c: int, e: int
) -> int:
    """Evaluate the proposed functional directly on one projected generator."""

    s = row_start(ell)
    total = 0
    for n, xdeg, coeff in generator_terms(a, b, c, e, m):
        if s <= n <= m and xdeg == 2 * n - ell:
            total += (-1) ** (n - s) * C(m + q - n, q) * coeff
    return total


def L_closed(
    m: int, ell: int, q: int, a: int, b: int, c: int, e: int
) -> int:
    """Closed evaluation using full or truncated finite differences."""

    s = row_start(ell)
    if a != 2 * c + 3 * e - ell:
        return 0
    N = c + e
    n0 = b + c + 2 * e
    nmax = n0 + N
    if n0 > m:
        return 0
    if nmax <= m + q:
        # Here R=m+q-n0 >= N, so only ordinary nonnegative-top
        # binomial coefficients occur in the complete difference.
        return (-1) ** (n0 - s) * C(m + q - nmax, q - N)

    delta = nmax - (m + q)
    K = m - n0
    require(K == N - q - delta >= 0, (m, ell, q, a, b, c, e, K, delta))
    return (-1) ** (n0 - s + K) * C(N - q - 1, delta - 1)


def projected_generators(d: int, m: int):
    """All THM-4308 generators with a+b<=d and a nonzero row<=m term."""

    for a in range(d + 1):
        for b in range(d - a + 1):
            for e in range((m - b) // 2 + 1):
                for c in range(m - b - 2 * e + 1):
                    yield a, b, c, e


def p0_representations(ell: int):
    """All P_0 packets that meet r=2n-ell, indexed by packet length N."""

    lo = min_packet(ell)
    hi = ell // 2
    for N in range(lo, hi + 1):
        c = 3 * N - ell
        e = ell - 2 * N
        require(c >= 0 and e >= 0 and 2 * c + 3 * e == ell, (ell, N, c, e))
        yield N, c, e


checks = 0
closed_form_columns = 0
p0_parameter_cases = 0
sufficient_columns = 0


# The equality 2c+3e=ell gives exactly the consecutive packet lengths
# ceil(ell/3),...,floor(ell/2).  This range drives both necessities.
for ell in range(2, 81):
    reps = list(p0_representations(ell))
    expected = list(range(min_packet(ell), ell // 2 + 1))
    require([N for N, _, _ in reps] == expected, (ell, reps, expected))
    checks += 1 + 3 * len(reps)


# Audit the exact one-generator closed form across all support, truncation,
# and order regimes.  This is independent of the universal proof in REPORT.md.
for ell in range(2, 19):
    s = row_start(ell)
    for q in range(0, min_packet(ell) + 5):
        for m in range(s, ell + 5):
            for d in range(0, 7):
                for a, b, c, e in projected_generators(d, m):
                    direct = L_direct(m, ell, q, a, b, c, e)
                    closed = L_closed(m, ell, q, a, b, c, e)
                    require(
                        direct == closed,
                        ("closed form", ell, q, m, d, a, b, c, e, direct, closed),
                    )
                    closed_form_columns += 1
                    checks += 1


# Exact P_0 classification on a broad grid:
#   L annihilates pi_m(P_0) iff q<ceil(ell/3) and m+q>=ell.
for ell in range(2, 41):
    s = row_start(ell)
    rho = min_packet(ell)
    reps = list(p0_representations(ell))
    for q in range(0, rho + 7):
        for m in range(s, ell + 8):
            values = [L_direct(m, ell, q, 0, 0, c, e) for _, c, e in reps]
            observed = all(value == 0 for value in values)
            predicted = q < rho and m + q >= ell
            require(observed == predicted, ("P0 iff", ell, q, m, values))
            p0_parameter_cases += 1
            checks += 1

            if not predicted:
                # Construct a deterministic witness for every excluded pair.
                if m + q < ell:
                    N = max(rho, ell - m)
                    c, e = 3 * N - ell, ell - 2 * N
                    delta = ell - (m + q)
                    expected_abs = C(N - q - 1, delta - 1)
                    expected_value = (-1) ** (m - s) * expected_abs
                else:
                    N = max(rho, ell - m)
                    c, e = 3 * N - ell, ell - 2 * N
                    expected_abs = C(m + q - ell, q - N)
                    expected_value = (-1) ** (ell - N - s) * expected_abs
                value = L_direct(m, ell, q, 0, 0, c, e)
                require(c >= 0 and e >= 0 and value == expected_value != 0,
                        ("P0 witness", ell, q, m, N, c, e, value, expected_value))
                checks += 1


# Unit hostile at the first excluded finite-difference order q=rho.  Put it
# at the first row satisfying the endpoint inequality, m=ell-q.
for ell in range(2, 121):
    s = row_start(ell)
    q = min_packet(ell)
    m = ell - q
    N = q
    c, e = 3 * N - ell, ell - 2 * N
    value = L_direct(m, ell, q, 0, 0, c, e)
    require(m >= s and c >= 0 and e >= 0, ("order indices", ell, q, m, c, e))
    require(value == (-1) ** (m - s) and abs(value) == 1,
            ("order unit", ell, q, m, value))
    checks += 3


# Unit hostile at the first excluded row m=ell-q-1 for every admissible q.
# The minimal packet N=rho always begins no later than this row.
for ell in range(2, 121):
    s = row_start(ell)
    rho = min_packet(ell)
    for q in range(rho):
        m = ell - q - 1
        N = rho
        c, e = 3 * N - ell, ell - 2 * N
        value = L_direct(m, ell, q, 0, 0, c, e)
        require(m >= s, ("row domain", ell, q, m, s))
        require(value == (-1) ** (m - s) and abs(value) == 1,
                ("row unit", ell, q, m, value))
        checks += 2


# In the admissible region, exhaust every projected generator through a broad
# finite range on the asserted side of the sharp depth threshold.
for ell in range(2, 21):
    rho = min_packet(ell)
    for q in range(rho):
        for m in range(ell - q, ell - q + 8):
            threshold = m + q - ell
            for d in range(threshold + 1):
                for a, b, c, e in projected_generators(d, m):
                    value = L_direct(m, ell, q, a, b, c, e)
                    require(value == 0,
                            ("depth zero", ell, q, m, d, a, b, c, e, value))
                    sufficient_columns += 1
                    checks += 1


# At the first excluded depth d0=m+q-ell+1, a unit hostile exists.  The parity
# correction is epsilon=(ell+d0) mod 2 (not d0 mod 2 when ell is odd).
for ell in range(2, 121):
    s = row_start(ell)
    rho = min_packet(ell)
    for q in range(rho):
        for m in range(ell - q, ell - q + 17):
            d0 = m + q - ell + 1
            S = ell + d0
            epsilon = S % 2
            a, b, e = d0, 0, epsilon
            c = (S - 3 * epsilon) // 2
            N = c + e
            n0 = c + 2 * e
            value = L_direct(m, ell, q, a, b, c, e)
            require(a + b == d0 and c >= 0 and e in (0, 1),
                    ("depth indices", ell, q, m, d0, a, b, c, e))
            require(a == 2 * c + 3 * e - ell,
                    ("depth support", ell, q, m, d0, a, c, e))
            require(N >= rho > q and n0 >= s and n0 + N == m + q + 1,
                    ("depth endpoints", ell, q, m, d0, N, n0))
            require(value == (-1) ** (m - s) and abs(value) == 1,
                    ("depth unit", ell, q, m, d0, value))
            checks += 5


# The first low-order and familiar tetrahedral specializations.
def stencil(m: int, ell: int, q: int):
    s = row_start(ell)
    return [(-1) ** (n - s) * C(m + q - n, q) for n in range(s, m + 1)]


require(stencil(8, 8, 2) == [15, -10, 6, -3, 1], "THM-4308 triangular stencil")
require(stencil(9, 10, 3) == [35, -20, 10, -4, 1], "tetrahedral m=9")
require(stencil(10, 10, 3) == [56, -35, 20, -10, 4, -1], "tetrahedral m=10")
require(L_direct(8, 8, 2, 3, 0, 4, 1) == 1, "THM-4308 triangular hostile")
require(L_direct(9, 10, 3, 3, 0, 5, 1) == 1, "THM-4358 tetrahedral hostile")
require(L_direct(10, 10, 3, 4, 0, 7, 0) == -1, "THM-4361 tetrahedral hostile")
checks += 6


# Audit the two binomial identities used in the proof with ordinary integers.
for N in range(1, 81):
    for q in range(0, N):
        for R in range(N, N + 20):
            full = sum((-1) ** k * C(N, k) * C(R - k, q) for k in range(N + 1))
            require(full == 0, ("full difference", N, q, R, full))
            checks += 1
        for delta in range(1, N - q + 1):
            K = N - q - delta
            partial = sum(
                (-1) ** k * C(N, k) * C(q + K - k, q)
                for k in range(K + 1)
            )
            expected = (-1) ** K * C(N - q - 1, delta - 1)
            require(partial == expected, ("partial difference", N, q, delta, partial))
            checks += 1


script_hash = hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest()
print("THM-4364 INDEPENDENT BOUNDARY REFEREE: PASS")
print("exact condition: q<ceil(ell/3) and m+q>=ell")
print("admissible-depth condition: d<=m+q-ell")
print(f"P0 parameter triples checked: {p0_parameter_cases}")
print(f"one-generator closed-form columns checked: {closed_form_columns}")
print(f"admissible projected columns checked: {sufficient_columns}")
print("unit boundaries checked: first excluded order, row, and depth through ell=120")
print("special stencils: q=2 triangular and q=3 tetrahedral (m=9,10)")
print(f"total grouped assertions: {checks}")
print(f"script_sha256={script_hash}")
