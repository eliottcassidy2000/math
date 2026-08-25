#!/usr/bin/env python3
"""Exact structural probe for the THM-4029 sixty-phase tail."""

from fractions import Fraction as Q
from math import comb
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "04-computation"))

from lrc14_ap_cover_finite_owner_formula_thm4029 import (  # noqa: E402
    owner_formula_deficit,
    prove_phase_rational_law,
)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def divisors(n):
    return tuple(d for d in range(1, n + 1) if n % d == 0)


def minimal_period(values):
    n = len(values)
    return next(
        d for d in divisors(n)
        if all(values[r] == values[(r + d) % n] for r in range(n))
    )


def poly_mul(a, b):
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return tuple(out)


def quotient_falling_polynomial(skip):
    """Ascending coefficients of product_{0<=d<=5,d!=skip}(n-d)."""
    out = (Q(1),)
    for d in range(6):
        if d != skip:
            out = poly_mul(out, (Q(-d), Q(1)))
    return out


def main():
    laws, constants = prove_phase_rational_law(period=60, minimum_n=12)
    require(constants == {Q(127, 35)}, "wrong leading deficit constant")

    # A_c(r) is the aggregate coefficient of 1/(n-c) at n mod 60 = r.
    A = [[Q(0) for _ in range(60)] for _ in range(6)]
    for r, (_n0, terms, _selectors) in laws.items():
        for coefficient, shift in terms:
            A[shift][r] += coefficient

    print("AGGREGATE SHIFT COEFFICIENTS")
    for c in range(6):
        support = sum(value != 0 for value in A[c])
        print(
            f"c={c};minimal_period={minimal_period(A[c])};"
            f"support_phases={support};distinct={len(set(A[c]))}"
        )
    require(
        tuple(minimal_period(A[c]) for c in range(6)) == (60, 60, 60, 60, 30, 6),
        "wrong aggregate shift periods",
    )
    vectors = tuple(tuple(A[c][r] for c in range(6)) for r in range(60))
    print(
        f"full_A_vector_minimal_period={minimal_period(vectors)};"
        f"distinct_phase_vectors={len(set(vectors))}"
    )
    require(minimal_period(vectors) == 60 and len(set(vectors)) == 60,
            "full phase vector is not exactly 60-periodic")

    # Resolve the phase source by owner denominator.  Fibonacci/Ostrowski
    # convergent denominators <=6 are 1,2,3,5; q=4 is the extra 2-adic clock.
    by_q = {}
    for q in range(1, 7):
        qvectors = []
        for r in range(60):
            row = [Q(0)] * 6
            for owner, _side, _missing, _track, (coefficient, shift) in laws[r][2]:
                if owner.denominator == q:
                    row[shift] += coefficient
            qvectors.append(tuple(row))
        by_q[q] = tuple(qvectors)
        print(
            f"owner_denominator_q={q};coefficient_vector_minimal_period="
            f"{minimal_period(by_q[q])};distinct={len(set(by_q[q]))}"
        )
        require(minimal_period(by_q[q]) == q, f"owner denominator {q} has wrong period")
    fib_q = (1, 2, 3, 5)
    fib_vectors = tuple(
        tuple(sum(by_q[q][r][c] for q in fib_q) for c in range(6))
        for r in range(60)
    )
    without_q4 = tuple(
        tuple(sum(by_q[q][r][c] for q in (1, 2, 3, 5, 6)) for c in range(6))
        for r in range(60)
    )
    print(f"fibonacci_denominator_skeleton_period={minimal_period(fib_vectors)}")
    print(f"all_owners_except_q4_period={minimal_period(without_q4)}")
    require(minimal_period(fib_vectors) == 30, "Fibonacci denominator skeleton period is not 30")
    require(minimal_period(without_q4) == 30, "q=4 is not the required extra clock")

    # Clear Q_6(n)=n(n-1)...(n-5).  The numerator is a phase-polynomial.
    quotients = tuple(quotient_falling_polynomial(c) for c in range(6))
    numerator_polys = []
    for r in range(60):
        coefficients = [Q(0)] * 6
        for c in range(6):
            for degree, value in enumerate(quotients[c]):
                # The owner formula is D=(1/7) sum_c A_c(r)/(n-c).
                coefficients[degree] += A[c][r] * value / 7
        numerator_polys.append(tuple(coefficients))
    print("CLEARED-DENOMINATOR PHASE POLYNOMIAL")
    for degree in range(6):
        values = tuple(poly[degree] for poly in numerator_polys)
        print(
            f"degree={degree};minimal_period={minimal_period(values)};"
            f"distinct={len(set(values))};first={values[0]}"
        )
    print(
        f"phase_polynomial_minimal_period={minimal_period(tuple(numerator_polys))};"
        f"distinct_phase_polynomials={len(set(numerator_polys))}"
    )
    coefficient_periods = tuple(
        minimal_period(tuple(poly[degree] for poly in numerator_polys))
        for degree in range(6)
    )
    require(coefficient_periods == (60, 60, 60, 60, 60, 1),
            "wrong cleared-polynomial coefficient periods")
    require(minimal_period(tuple(numerator_polys)) == 60
            and len(set(numerator_polys)) == 60,
            "cleared phase polynomials are not exactly 60 distinct phases")
    require({poly[5] for poly in numerator_polys} == {Q(127, 35)},
            "wrong normalized leading coefficient")

    # Laurent-at-infinity phase moments: sum_c A_c(r)c^k.
    print("ASYMPTOTIC PHASE MOMENTS")
    for k in range(8):
        values = tuple(sum(A[c][r] * Q(c) ** k for c in range(6)) for r in range(60))
        print(
            f"k={k};minimal_period={minimal_period(values)};"
            f"distinct={len(set(values))};range=({min(values)},{max(values)})"
        )
        if k == 0:
            require(set(values) == {Q(127, 5)}, "wrong raw zeroth phase moment")
        if k == 1:
            require(minimal_period(values) == 60 and len(set(values)) == 25,
                    "wrong second-order phase moment")
            require((min(values), max(values)) == (Q(131, 10), Q(39, 2)),
                    "wrong second-order phase range")

    # Exact sixth 60-step finite-difference recurrence after clearing Q_6.
    def falling6(n):
        out = 1
        for d in range(6):
            out *= n - d
        return out

    failures = []
    for n in range(12, 500):
        lhs = Q(0)
        for j in range(7):
            nj = n + 60 * j
            b = owner_formula_deficit(nj + 1)
            lhs += Q((-1) ** (6 - j) * comb(6, j)) * falling6(nj) * b
        if lhs:
            failures.append((n, lhs))
            break
    print(f"cleared_tail_recurrence_n12_499={not failures};failures={failures}")
    require(not failures, "cleared recurrence failed")

    # Triangular/quasipolynomial track clock: E_s(n) is largest <=n congruent s mod q.
    # For N=qK+r, summing from n=q avoids the pre-arrival range n<s and gives
    # q^2*T_{K-1}+(r+1)(qK+s)-q*min(r+1,s).
    print("TRACK CLOCK TRIANGULAR DECOMPOSITION")
    for q in range(1, 7):
        for s in range(q):
            for r in range(q):
                for K in range(1, 30):
                    N = q * K + r
                    total = sum(n - ((n - s) % q) for n in range(q, N + 1))
                    triangular = K * (K - 1) // 2
                    formula = q * q * triangular + (r + 1) * (q * K + s) - q * min(r + 1, s)
                    require(total == formula, f"track clock identity failed at {(q, s, K, r)}")
            print(f"q={q};s={s};period={q};triangular_formula=PASS_K1_29")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
