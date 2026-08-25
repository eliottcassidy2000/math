#!/usr/bin/env python3
"""Exact audit for THM-4031's endpoint-owner tail on THM-4025 rays.

The audit uses two representations: Fraction evaluation of the inherited
gate and an integer numerator for the endpoint certificate.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")


def require(condition: bool, *data: object) -> None:
    """Optimization-safe audit gate."""
    if not condition:
        raise AssertionError(data)


def posres(a: int, m: int) -> int:
    r = a % m
    require(r != 0, "zero residue", a, m)
    return r


def gaps(t: int, u: int) -> tuple[F, F]:
    g = gcd(t, u)
    return (
        F(posres(3 * t - 4 * u, 42 * g), 42 * u),
        F(posres(9 * t + 16 * u, 126 * g), 126 * u),
    )


def gate_margin(t: int, U: int) -> F:
    e = [gaps(t, u) for u in range(1, U + 1)]
    eps1 = min(x for x, _ in e)
    eps2 = min(y for _, y in e)
    B = max(F(0), F(2, 63) - eps1 - eps2)
    D = F(t * (2 * U - 1), 84 * U * (U - 1))
    return B - D


def endpoint_packet(t: int, U: int, k: int) -> tuple[int, int, int, int, F]:
    require(t % 2 == k % 2 == 1 and 1 < U <= t, "domain", t, U, k)
    u = k * U - 1
    g_big = gcd(k * t, u)
    g_small = gcd(t, u)
    require(g_big == g_small, "gcd identity", t, U, k, g_big, g_small)
    r1 = posres(k * (3 * t - 4 * U) + 4, 42 * g_small)
    r2 = posres(k * (9 * t + 16 * U) - 16, 126 * g_small)
    h = 3 * r1 + r2
    e1, e2 = gaps(k * t, u)
    require(e1 + e2 == F(h, 126 * u), "endpoint gap sum", t, U, k)
    a = 4 * U - 3 * t
    numerator = 2 * a * u - 2 * U * h - 3 * t
    lower = F(numerator, 252 * U * u)
    D = F(k * t * (2 * k * U - 1), 84 * k * U * (k * U - 1))
    require(lower == F(2, 63) - e1 - e2 - D, "margin numerator", t, U, k)
    return g_small, r1, r2, numerator, lower


def ceil_div(a: int, b: int) -> int:
    require(b > 0, "ceil denominator", a, b)
    return -((-a) // b)


def phase_tail(t: int, U: int) -> tuple[int, int, tuple[int, int, int, int]]:
    """Return exact endpoint-certificate tail and its last hostile phase."""
    require(4 * U > 3 * t, "positive slope", t, U)
    P = 126 * t
    last_fail = -1
    worst = None
    for s in range(1, P, 2):
        g, r1, r2, n0, _ = endpoint_packet(t, U, s)
        g2, R1, R2, _, _ = endpoint_packet(t, U, s + P)
        require((g, r1, r2) == (g2, R1, R2), "period", t, U, s)
        slope = 2 * (4 * U - 3 * t) * P * U
        jumps = max(0, ceil_div(-n0, slope))
        first_good = s + jumps * P
        require(endpoint_packet(t, U, first_good)[3] >= 0, "phase threshold", t, U, s)
        if first_good >= P:
            lf = first_good - P
            require(endpoint_packet(t, U, lf)[3] < 0, "last phase failure", t, U, s)
            if lf > last_fail:
                last_fail = lf
                worst = (s, g, r1, r2)
    tail = 1 if last_fail < 1 else last_fail + 2
    require(tail % 2 == 1, "tail parity", t, U, tail)
    return tail, last_fail, worst or (1, 0, 0, 0)


def smallest_odd_at_least(num: int, den: int) -> int:
    q = ceil_div(num, den)
    if q % 2 == 0:
        q += 1
    return max(q, 1)


def crude_tail(t: int, U: int) -> int:
    """All-phase bound from h <= 252t-4."""
    a = 4 * U - 3 * t
    hmax = 252 * t - 4
    # 2a(kU-1) >= 2Uhmax+3t.
    return smallest_odd_at_least(2 * a + 2 * U * hmax + 3 * t, 2 * a * U)


def main() -> None:
    checks = 0
    # Fraction-vs-integer and actual-margin domination controls.
    for t in range(3, 22, 2):
        for U in range(2, t + 1):
            if 4 * U <= 3 * t:
                continue
            for k in range(1, 16, 2):
                *_, lower = endpoint_packet(t, U, k)
                actual = gate_margin(k * t, k * U)
                require(actual >= lower, "actual dominates endpoint", t, U, k)
                checks += 1

    rays = ((5, 4), (11, 9), (11, 11), (13, 10), (17, 13), (101, 76))
    print("LRC14_OWNER_ENDPOINT_RATIONAL_TAIL_THM4031_EXACT_AUDIT")
    print(f"fraction_integer_checks={checks}")
    for t, U in rays:
        tail, last_fail, phase = phase_tail(t, U)
        crude = crude_tail(t, U)
        require(tail <= crude, "phase tail versus universal tail", t, U)
        # Sample the exact theorem at and beyond the tail without minimizing
        # huge owner universes for the largest ray.
        for k in (tail, tail + 2, tail + 126 * t):
            require(endpoint_packet(t, U, k)[3] >= 0, "tail sample", t, U, k)
        actual_last_closure = None
        if t <= 17:
            for k in range(1, min(tail + 2, 103), 2):
                if gate_margin(k * t, k * U) < 0:
                    actual_last_closure = k
        print(
            f"ray={t},{U};a={4*U-3*t};period={126*t};"
            f"endpoint_tail={tail};endpoint_last_fail={last_fail};"
            f"worst_phase={phase};crude_tail={crude};"
            f"actual_last_closure_le_101={actual_last_closure}"
        )

    # Complete hostile/head audit for the canonical revival ray: brute force
    # every odd scale below the proved endpoint tail, then use the theorem.
    canonical_tail = phase_tail(5, 4)[0]
    canonical_closures = tuple(
        k for k in range(1, canonical_tail, 2) if gate_margin(5 * k, 4 * k) < 0
    )
    require(
        canonical_closures == (1, 3, 5, 7, 11, 13),
        "canonical complete head",
        canonical_closures,
    )
    print(f"ray=5,4;complete_closures={canonical_closures};proved_tail={canonical_tail}")

    # Hostile boundary: on the other side the demand already exceeds the
    # largest possible B, at every odd dilation.
    for k in range(1, 100, 2):
        t, U = 3, 2
        D = F(k * t * (2 * k * U - 1), 84 * k * U * (k * U - 1))
        require(D > F(2, 63), "other-strip hostile", k, D)

    # The scaled endpoint u=kU does not decay: this is why the -1 sidecar is
    # load-bearing rather than cosmetic.
    for t, U in ((5, 4), (11, 9), (17, 13)):
        base = gaps(t, U)
        for k in range(1, 40, 2):
            require(gaps(k * t, k * U) == base, "scaled endpoint hostile", t, U, k)

    print("hostile_other_strip=(t,U)=(3,2):D>2/63>=B at every odd scale")
    print("hostile_scaled_endpoint=u=kU:gaps are invariant, not O(1/k)")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
