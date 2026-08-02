#!/usr/bin/env python3
"""Exact companion for THM-3152 on the resonant quadratic factorial pair.

The multi-place Newton degree-barcode lemma is proved abstractly; the declared
censuses and symbolic checks are finite-exact.  All polynomial arithmetic is
over Z (python-flint), and all valuations are exact.
"""

from fractions import Fraction
from math import comb, factorial, isqrt

import gmpy2
import sympy as sp
from flint import fmpz_poly


X = fmpz_poly([0, 1])


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % q for q in range(3, isqrt(n) + 1, 2))


def factor(n):
    ans = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            e = 0
            while n % p == 0:
                n //= p
                e += 1
            ans.append((p, e))
        p += 1 if p == 2 else 2
    if n > 1:
        ans.append((n, 1))
    return ans


def is_prime_power(n):
    return len(factor(n)) == 1


def pair(d):
    """Return A_(d-2)^d and A_(d-1)^d from THM-3124's recurrence."""
    old = fmpz_poly([1])
    current = fmpz_poly([d - 1, 2])
    d_power = d
    for n in range(1, d - 1):
        following = (
            fmpz_poly([d_power * (d - n - 1)])
            + X * current * (2 * (n + 1) * (2 * n + 1))
            + (old - X * old * (4 * d)) * (n * (n + 1))
        )
        old, current = current, following
        d_power *= d
    require(old.degree() == d - 2 and current.degree() == d - 1, (d, "degree"))
    return old, current


def direct_coefficients(n, d):
    facts = [factorial(k) for k in range(2 * n + 1)]
    return [
        comb(n, j)
        * sum(
            comb(n - j, ell)
            * d ** (n - j - ell)
            * (-1) ** ell
            * facts[2 * j + ell]
            for ell in range(n - j + 1)
        )
        for j in range(n + 1)
    ]


def vp(a, p):
    a = gmpy2.mpz(a)
    if not a:
        return None
    return int(gmpy2.remove(a, p)[1])


def lower_hull(poly, p):
    hull = []
    for x in range(poly.degree() + 1):
        y = vp(poly[x], p)
        if y is None:
            continue
        point = (x, y)
        while len(hull) >= 2:
            x0, y0 = hull[-2]
            x1, y1 = hull[-1]
            x2, y2 = point
            if Fraction(y1 - y0, x1 - x0) >= Fraction(y2 - y1, x2 - x1):
                hull.pop()
            else:
                break
        hull.append(point)
    return hull


def slope_lengths(poly, p):
    hull = lower_hull(poly, p)
    return {
        Fraction(y1 - y0, x1 - x0): x1 - x0
        for (x0, y0), (x1, y1) in zip(hull, hull[1:])
    }


def degree_barcode(P, Q, p):
    """Necessary degrees of a common Q-factor, viewed over Q_p."""
    return degree_barcode_many((P, Q), p)


def zero_order(poly):
    """Multiplicity of the coordinate root v=0."""
    for j in range(poly.degree() + 1):
        if poly[j] != 0:
            return j
    raise RuntimeError("zero polynomial has no Newton barcode")


def degree_barcode_many(polynomials, p):
    """Necessary common-factor degrees after retaining several exact rows."""
    ledgers = [slope_lengths(poly, p) for poly in polynomials]
    common_slopes = set(ledgers[0])
    for ledger in ledgers[1:]:
        common_slopes &= set(ledger)
    # Finite lower slopes omit the coordinate root v=0.  Retain it as a
    # separate denominator-one block; this makes the abstract lemma valid
    # without assuming nonzero constant coefficients.
    zero_capacity = min(zero_order(poly) for poly in polynomials)
    degrees = set(range(zero_capacity + 1))
    common = []
    for slope in sorted(common_slopes):
        capacity = min(ledger[slope] for ledger in ledgers)
        denominator = slope.denominator
        common.append((slope, capacity, denominator))
        degrees = {
            old + use
            for old in degrees
            for use in range(0, capacity + 1, denominator)
        }
    return degrees, tuple(common)


def first_full_pseudoremainder(P, Q, d):
    """Integral first full Euclidean row for n=d-2."""
    n = d - 2
    leading_quotient = 2 * (n + 1) * (2 * n + 1)
    remainder = (
        (2 * n - 1) * (Q - leading_quotient * X * P)
        + 2 * d * (n + 1) * P
    )
    require(remainder.degree() <= n - 1, (d, "pseudo-remainder degree"))
    return remainder


def intersect_trace_many(polynomials, primes):
    possible = set(range(1, min(poly.degree() for poly in polynomials) + 1))
    trace = []
    for p in primes:
        local, barcode = degree_barcode_many(polynomials, p)
        narrowed = possible & local
        if narrowed != possible:
            trace.append((p, tuple(sorted(narrowed)), barcode))
        possible = narrowed
        if not possible:
            break
    return possible, trace


def intersect_trace(P, Q, primes):
    return intersect_trace_many((P, Q), primes)


def polynomial_identity_checks():
    """Check the exact first full remainder and differential syzygy."""
    v = sp.symbols("v")
    checked = []
    for d in range(3, 31):
        n = d - 2
        direct_P = direct_coefficients(n, d)
        direct_Q = direct_coefficients(n + 1, d)
        P = sp.Poly.from_list(direct_P[::-1], gens=v).as_expr()
        Q = sp.Poly.from_list(direct_Q[::-1], gens=v).as_expr()
        if d <= 15:
            recurrence_P, recurrence_Q = pair(d)
            require(
                [int(recurrence_P[j]) for j in range(n + 1)] == direct_P,
                (d, "direct/recurrence P"),
            )
            require(
                [int(recurrence_Q[j]) for j in range(n + 2)] == direct_Q,
                (d, "direct/recurrence Q"),
            )
        C = 2 * (n + 1) * (2 * n + 1)
        D = -sp.Rational(2 * (n + 1) * d, 2 * n - 1)
        full_remainder = sp.expand(Q - (C * v + D) * P)
        require(sp.degree(full_remainder, v) <= n - 1, (d, "full remainder"))

        Delta = 1 - 4 * d * v
        G = 2 * (n + 1) * Delta * v**2
        H = d ** (n + 1) * ((4 * n + 6) * v - 1)
        K = (n + 1) * (4 * d * v**2 - (6 * d - 2) * v + 1)
        differential = sp.expand((2 * d * v - 1) * Q - (H - G * sp.diff(P, v) + K * P))
        require(differential == 0, (d, "differential syzygy"))

        x = sp.symbols("x")
        transformed = sp.Poly(
            sp.expand((4 * d) ** n * P.subs(v, (x + 1) / (4 * d))), x
        )
        square_coeffs = []
        for j in range(n + 1):
            integral = sum(
                comb(2 * (n - j), ell)
                * (2 * d) ** (2 * (n - j) - ell)
                * (-1) ** ell
                * factorial(2 * j + ell)
                for ell in range(2 * (n - j) + 1)
            )
            require(integral > 0, (d, j, "square coefficient positivity"))
            square_coeffs.append(comb(n, j) * integral)
        require(
            transformed == sp.Poly(sum(c * x**j for j, c in enumerate(square_coeffs)), x),
            (d, "discriminant-square transform"),
        )
        checked.append(d)
    return checked


def bounded_barcode_census():
    tested = 0
    max_primes_used = 0
    worst = None
    for d in range(4, 251):
        if is_prime_power(d - 1):
            continue
        P, Q = pair(d)
        primes = [p for p in range(2, 2 * d + 1) if is_prime(p)]
        possible, trace = intersect_trace(P, Q, primes)
        require(not possible, (d, possible, "bounded barcode survivor"))
        tested += 1
        if len(trace) > max_primes_used:
            max_primes_used = len(trace)
            worst = (d, trace)
    return tested, max_primes_used, worst


def three_exit_residual(d):
    return (
        not is_prime(d)
        and not is_prime(d - 2)
        and not is_prime_power(d - 1)
    )


def high_three_exit_residual_census():
    bank = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
    rows = []
    for d in range(1001, 1101):
        if not three_exit_residual(d):
            continue
        P, Q = pair(d)
        R = first_full_pseudoremainder(P, Q, d)
        require(R.degree() == d - 3, (d, "high three-exit residual remainder degree"))
        possible, trace = intersect_trace_many((P, Q, R), bank)
        require(not possible, (d, possible, "high three-exit residual barcode survivor"))
        rows.append((d, factor(d), factor(d - 1), factor(d - 2), trace))
    return bank, rows


def planted_control():
    # Both polynomials contain X+1.  Every local barcode intersection must
    # retain degree one; this attacks an accidentally overstrong criterion.
    P = (X + 1) * (X**3 + 2 * X + 3)
    Q = (X + 1) * (2 * X**4 + X**2 + 5)
    possible, trace = intersect_trace(P, Q, (2, 3, 5, 7, 11, 13))
    require(1 in possible, (possible, "planted factor lost"))
    P_zero = X * (X**3 + 2 * X + 3)
    Q_zero = X * (2 * X**4 + X**2 + 5)
    possible_zero, trace_zero = intersect_trace(P_zero, Q_zero, (2, 3, 5, 7, 11, 13))
    require(1 in possible_zero, (possible_zero, "coordinate factor lost"))
    return tuple(sorted(possible)), trace, tuple(sorted(possible_zero)), trace_zero


def main():
    print("FACTORIAL MULTI-PLACE NEWTON DEGREE BARCODE EXACT COMPANION FOR THM-3152")
    print("status=PROVED abstract lemma + FINITE-EXACT census; current PROVED THM-3146 is not used as a dependency")
    checked = polynomial_identity_checks()
    print("symbolic identities checked d=%d..%d count=%d" % (checked[0], checked[-1], len(checked)))
    print("identity: full quotient Q/P begins C*v-2(n+1)(n+2)/(2n-1)")
    print("identity: (2dv-1)Q=H-GP'+KP with H=d^(n+1)((4n+6)v-1)")
    print("identity: (4d)^n P((x+1)/(4d))=L((((2d-t)^2+x*t^2)^n))")

    tested, depth, worst = bounded_barcode_census()
    print("bounded non-prime-power predecessors d<=250 tested=%d survivors=0" % tested)
    print("max narrowing events=%d first_worst_d=%d" % (depth, worst[0]))

    bank, rows = high_three_exit_residual_census()
    print("high three-exit residual universe=1001..1100 rows=%d survivors=0 bank=%s" % (len(rows), bank))
    print("high census carrier=(P,Q,first full integral Euclidean pseudoremainder)")
    print("high three-exit residual d-list=%s" % ([d for d, _, _, _, _ in rows],))
    print(
        "high closure trace (d, narrowing primes, last nonempty size)=%s"
        % (
            [
                (
                    d,
                    tuple(p for p, _, _ in trace),
                    len(trace[-2][1]) if len(trace) >= 2 else d - 3,
                )
                for d, _, _, _, trace in rows
            ],
        )
    )

    P, Q = pair(1001)
    possible, trace = intersect_trace(P, Q, (2, 3))
    require(not possible, "d=1001 two-place witness")
    print(
        "d=1001 p=2 common blocks=%s"
        % ([(str(s), c, b) for s, c, b in trace[0][2]],)
    )
    print("d=1001 D_2 positive=%s" % (trace[0][1],))
    print(
        "d=1001 p=3 common blocks=%s"
        % ([(str(s), c, b) for s, c, b in trace[1][2]],)
    )
    local3, _ = degree_barcode(P, Q, 3)
    print("d=1001 D_3 positive=%s" % (tuple(sorted(local3 - {0})),))
    print("d=1001 D_2 intersection D_3 positive=() => gcd_Q(P,Q)=1")

    R = first_full_pseudoremainder(P, Q, 1001)
    possible_euclidean, trace_euclidean = intersect_trace_many((P, Q, R), (2, 3))
    require(not possible_euclidean, "d=1001 Euclidean-flag witness")
    local2e, blocks2e = degree_barcode_many((P, Q, R), 2)
    local3e, blocks3e = degree_barcode_many((P, Q, R), 3)
    print("d=1001 Euclidean p=2 common blocks=%s" % ([(str(s), c, b) for s, c, b in blocks2e],))
    print("d=1001 Euclidean D_2 positive=%s" % (tuple(sorted(local2e - {0})),))
    print("d=1001 Euclidean p=3 common blocks=%s" % ([(str(s), c, b) for s, c, b in blocks3e],))
    print("d=1001 Euclidean D_3 positive=%s" % (tuple(sorted(local3e - {0})),))
    print("d=1001 Euclidean D_2 intersection D_3 positive=() => gcd_Q(P,Q)=1")

    planted, planted_trace, planted_zero, planted_zero_trace = planted_control()
    print("planted degree-one common factor retained degrees=%s trace_events=%d" % (planted, len(planted_trace)))
    print("planted coordinate factor v retained degrees=%s trace_events=%d" % (planted_zero, len(planted_zero_trace)))
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
