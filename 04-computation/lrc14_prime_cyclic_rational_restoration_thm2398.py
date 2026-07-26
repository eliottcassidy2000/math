#!/usr/bin/env python3
"""Exact finite controls for THM-2398.

The proof itself is cyclotomic algebra.  This companion independently
checks finite coefficient banks, disjoint cross-correlations, the
prime/composite boundary, and the quantitative constants used in the
LRC interface.  Vanishing, counting, determinant, and rational-constant
checks are exact.  Trigonometric magnitude inequalities and the explicitly
labelled analytic hostiles are secondary floating controls after their
exact proofs or formulas are stated in the theorem.
"""

from __future__ import annotations

import cmath
import itertools
import math
from fractions import Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def cyclotomic_vector(coefficients: tuple[int, ...], k: int) -> tuple[int, ...]:
    """Coordinates of sum_s coefficients[s] zeta^(k*s) in Q(zeta_p).

    The basis is 1,zeta,...,zeta^(p-2), and
    zeta^(p-1)=-(1+...+zeta^(p-2)).
    """

    p = len(coefficients)
    require(p >= 3, "prime order must be at least three")
    out = [0] * (p - 1)
    for s, coefficient in enumerate(coefficients):
        exponent = (k * s) % p
        if exponent == p - 1:
            for j in range(p - 1):
                out[j] -= coefficient
        else:
            out[exponent] += coefficient
    return tuple(out)


def vanishes(coefficients: tuple[int, ...], k: int) -> bool:
    return not any(cyclotomic_vector(coefficients, k))


def coefficient_sweep(p: int, alphabet: tuple[int, ...]) -> tuple[int, int]:
    total = 0
    flat = 0
    for coefficients in itertools.product(alphabet, repeat=p):
        total += 1
        is_flat = len(set(coefficients)) == 1
        flat += int(is_flat)
        for k in range(1, p):
            require(
                vanishes(coefficients, k) == is_flat,
                f"prime-cyclic dichotomy failed at p={p}, k={k}, "
                f"coefficients={coefficients}",
            )
    return total, flat


def cyclic_cross_correlation(
    left: tuple[int, ...], right: tuple[int, ...]
) -> tuple[int, ...]:
    p = len(left)
    require(len(right) == p, "cross-correlation size mismatch")
    return tuple(
        sum(left[(r + shift) % p] * right[r] for r in range(p))
        for shift in range(p)
    )


def disjoint_pair_sweep(p: int) -> int:
    checked = 0
    for assignment in itertools.product(range(3), repeat=p):
        left = tuple(int(x == 1) for x in assignment)
        right = tuple(int(x == 2) for x in assignment)
        if not any(left) or not any(right):
            continue
        checked += 1
        correlation = cyclic_cross_correlation(left, right)
        require(correlation[0] == 0, "disjoint anchor was not zero")
        require(
            sum(correlation) == sum(left) * sum(right),
            "cross-correlation mass identity failed",
        )
        for k in range(1, p):
            require(
                not vanishes(correlation, k),
                "a disjoint positive prime-fibre correlation lost a colour",
            )
    expected = 3**p - 2 * 2**p + 1
    require(checked == expected, "disjoint-pair count mismatch")
    return checked


def fixed_two_root_sweep() -> tuple[int, Fraction]:
    """All singleton/adjacent A and disjoint two-root F controls at p=13."""

    p = 13
    checked = 0
    rational_floor = Fraction(4, 13**4)
    numerical_min = float("inf")
    zeta = cmath.exp(2j * math.pi / p)

    for left_support in ((0,), (0, 1)):
        available = [r for r in range(p) if r not in left_support]
        left = tuple(int(r in left_support) for r in range(p))
        for right_support in itertools.combinations(available, 2):
            right = tuple(int(r in right_support) for r in range(p))
            correlation = cyclic_cross_correlation(left, right)
            checked += 1
            require(correlation[0] == 0, "fixed-pattern anchor failed")
            for k in range(1, p):
                require(
                    not vanishes(correlation, k),
                    "fixed two-root pattern lost a prime colour",
                )
                a = sum(left[r] * zeta ** (-k * r) for r in range(p)) / p
                b = sum(right[r] * zeta ** (-k * r) for r in range(p)) / p
                product_abs = abs(a * b)
                numerical_min = min(numerical_min, product_abs)
                require(
                    product_abs > float(rational_floor),
                    "strict two-root rational floor failed",
                )

    require(checked == math.comb(12, 2) + math.comb(11, 2), "pattern count")
    return checked, rational_floor


def bareiss_determinant(matrix: list[list[int]]) -> int:
    """Fraction-free exact determinant."""

    n = len(matrix)
    require(all(len(row) == n for row in matrix), "matrix is not square")
    a = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if a[i][k]), None)
            if pivot is None:
                return 0
            a[k], a[pivot] = a[pivot], a[k]
            sign = -sign
        pivot_value = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot_value - a[i][k] * a[k][j]
                require(numerator % previous == 0, "Bareiss division failed")
                a[i][j] = numerator // previous
        previous = pivot_value
        for i in range(k + 1, n):
            a[i][k] = 0
    return sign * a[n - 1][n - 1]


def multiply_cyclotomic_vectors(
    left: tuple[int, ...], right: tuple[int, ...], p: int
) -> tuple[int, ...]:
    require(len(left) == len(right) == p - 1, "basis length mismatch")
    raw = [0] * p
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            raw[(i + j) % p] += a * b
    top = raw[p - 1]
    return tuple(raw[j] - top for j in range(p - 1))


def cyclotomic_norm(coefficients: tuple[int, ...]) -> int:
    p = len(coefficients)
    value = cyclotomic_vector(coefficients, 1)
    columns: list[tuple[int, ...]] = []
    for power in range(p - 1):
        basis = tuple(int(j == power) for j in range(p - 1))
        columns.append(multiply_cyclotomic_vectors(value, basis, p))
    matrix = [[columns[j][i] for j in range(p - 1)] for i in range(p - 1)]
    return bareiss_determinant(matrix)


def norm_controls() -> tuple[int, int]:
    samples = (
        (1,) + (0,) * 12,
        (1, 1) + (0,) * 11,
        (1, 2, 3) + (0,) * 10,
        tuple(range(13)),
        (0, 2, 0, 1, 3, 0, 4, 1, 0, 2, 5, 0, 1),
    )
    minimum_abs_norm = None
    for coefficients in samples:
        norm = cyclotomic_norm(coefficients)
        require(norm != 0, f"sample norm vanished: {coefficients}")
        minimum_abs_norm = (
            abs(norm)
            if minimum_abs_norm is None
            else min(minimum_abs_norm, abs(norm))
        )
        total = sum(coefficients)
        coefficient_gcd = 0
        for coefficient in coefficients:
            coefficient_gcd = math.gcd(coefficient_gcd, coefficient)
        variance = 13 * sum(a * a for a in coefficients) - total * total
        require(coefficient_gcd > 0 and variance > 0, "sample norm typing")
        require(
            abs(norm) % coefficient_gcd**12 == 0,
            "cyclotomic norm missed the coefficient gcd",
        )
        zeta = cmath.exp(2j * math.pi / 13)
        squared_sum = 0.0
        for k in range(1, 13):
            value = sum(
                coefficients[s] * zeta ** (k * s) for s in range(13)
            )
            squared_sum += abs(value) ** 2
            require(
                abs(value) + 1e-12 >= coefficient_gcd**6 / total**5,
                "paired cyclotomic norm floor failed",
            )
            require(
                abs(value) ** 2 + 1e-11
                >= coefficient_gcd**12 * (10 / variance) ** 5,
                "variance-sharpened norm floor failed",
            )
        require(
            abs(squared_sum - variance) < 1e-8,
            "nonzero-character Parseval variance failed",
        )
    require(minimum_abs_norm is not None, "no norm controls")
    return len(samples), minimum_abs_norm


def composite_hostile() -> tuple[Fraction, Fraction]:
    # C_4, kernel (0,1/2,0,1/2): m_1=(-i+i)/2=0.
    real = Fraction(0)
    imag = -Fraction(1, 2) + Fraction(1, 2)
    require(real == 0 and imag == 0, "composite hostile did not vanish")
    return real, imag


def analytic_hostile_residuals() -> tuple[float, float, float]:
    # Rationality hostile on C_13:
    # kappa_s=(1+eps*cos(4*pi*s/13))/13 is positive and m_1=0.
    p = 13
    epsilon = 0.5
    zeta = cmath.exp(2j * math.pi / p)
    kernel = tuple(
        (1.0 + epsilon * math.cos(4.0 * math.pi * s / p)) / p
        for s in range(p)
    )
    require(min(kernel) > 0.0, "analytic hostile was not strictly positive")
    require(max(kernel) - min(kernel) > 0.01, "analytic hostile was flat")
    m1 = sum(kernel[s] * zeta ** (-s) for s in range(p))

    # Anchored two-sided irrational hostile.  The s=0 weight is zero;
    # for s!=0 the exact weights are
    # (11+2*cos(2*pi*s/13))/130.  Their k=1 transform is
    # 11*(-1)+(12-1)=0.
    anchored = [0.0] + [
        (11.0 + 2.0 * math.cos(2.0 * math.pi * s / p)) / 130.0
        for s in range(1, p)
    ]
    require(min(anchored[1:]) > 0.0 and anchored[0] == 0.0, "anchor typing")
    require(abs(sum(anchored) - 1.0) < 1e-14, "anchor hostile mass")
    anchored_m1 = sum(anchored[s] * zeta ** (-s) for s in range(p))

    # Strictly positive rational nonuniform hostile on C_91:
    # half full-uniform plus half uniform on the order-seven subgroup.
    n = 91
    root = cmath.exp(2j * math.pi / n)
    kernel91 = [1.0 / (2 * n)] * n
    for j in range(7):
        kernel91[13 * j] += 1.0 / 14
    m91 = sum(kernel91[s] * root ** (-s) for s in range(n))
    require(max(kernel91) > min(kernel91) > 0.0, "C_91 hostile typing")
    require(abs(m1) < 1e-13, "irrational C_13 hostile residual too large")
    require(
        abs(anchored_m1) < 1e-13,
        "anchored irrational C_13 hostile residual too large",
    )
    require(abs(m91) < 1e-12, "rational C_91 hostile residual too large")
    return abs(m1), abs(anchored_m1), abs(m91)


def terminal_alignment_hostiles() -> tuple[int, int]:
    p = 13

    # Separate aggregate spectra but no common base co-support.
    left_rows = (
        (tuple(int(r == 0) for r in range(p)), Fraction(1, 2)),
        ((0,) * p, Fraction(1, 2)),
    )
    right_rows = (
        ((0,) * p, Fraction(1, 2)),
        (tuple(int(r == 1) for r in range(p)), Fraction(1, 2)),
    )
    correlation = [Fraction(0) for _ in range(p)]
    for row_index in range(2):
        left, weight = left_rows[row_index]
        right, right_weight = right_rows[row_index]
        require(weight == right_weight, "hostile base weights")
        row_correlation = cyclic_cross_correlation(left, right)
        for shift in range(p):
            correlation[shift] += weight * row_correlation[shift]
    require(not any(correlation), "co-support hostile correlation was nonzero")

    left_aggregate = tuple(Fraction(int(r == 0), 2) for r in range(p))
    right_aggregate = tuple(Fraction(int(r == 1), 2) for r in range(p))
    require(
        not vanishes(tuple(int(2 * x) for x in left_aggregate), 1),
        "left aggregate spectrum vanished",
    )
    require(
        not vanishes(tuple(int(2 * x) for x in right_aggregate), 1),
        "right aggregate spectrum vanished",
    )

    # Both packets coexist on every fibre, but without a pinned anchor the
    # rotating singleton bank gives the exactly uniform correlation.
    flat_correlation = [Fraction(0) for _ in range(p)]
    right = tuple(int(r == 0) for r in range(p))
    for translate in range(p):
        left = tuple(int(r == translate) for r in range(p))
        row_correlation = cyclic_cross_correlation(left, right)
        for shift in range(p):
            flat_correlation[shift] += Fraction(row_correlation[shift], p)
    require(
        len(set(flat_correlation)) == 1
        and flat_correlation[0] == Fraction(1, p),
        "flat restoration hostile failed",
    )
    for k in range(1, p):
        require(
            vanishes(tuple(1 for _ in range(p)), k),
            "flat restoration charged mode survived",
        )
    return 2, p


def lrc_constants() -> tuple[Fraction, Fraction, Fraction, Fraction]:
    # A fixed 13-unit two-root word has only thirteen cyclic translates.
    pattern_multiplier = Fraction(4, 13 * 13**4)
    require(pattern_multiplier == Fraction(4, 371_293), "pattern floor")
    last_lane = pattern_multiplier * Fraction(1, 1_391_208)
    common_core = pattern_multiplier * Fraction(33, 115_934)
    owner_cell = pattern_multiplier * Fraction(33, 753_571)
    owner_septimal = owner_cell / 7
    require(last_lane == Fraction(1, 129_136_447_986), "last-lane floor")
    require(common_core == Fraction(66, 21_522_741_331), "common-core floor")
    require(owner_cell == Fraction(132, 279_795_637_303), "owner-cell floor")
    require(
        owner_septimal == Fraction(132, 1_958_569_461_121),
        "septimal floor",
    )
    return last_lane, common_core, owner_cell, owner_septimal


def main() -> None:
    print("THM-2398 prime-cyclic rational restoration exact companion")

    for p, alphabet in (
        (3, (0, 1, 2)),
        (5, (0, 1, 2)),
        (7, (0, 1, 2)),
        (13, (0, 1)),
    ):
        total, flat = coefficient_sweep(p, alphabet)
        print(f"prime sweep p={p}: kernels={total}, flat={flat}, PASS")

    disjoint = disjoint_pair_sweep(7)
    print(f"disjoint two-sided sweep p=7: pairs={disjoint}, PASS")

    patterns, rational_floor = fixed_two_root_sweep()
    print(
        "p=13 singleton/adjacent x disjoint-two-root patterns: "
        f"{patterns}, strict product floor>{rational_floor}, PASS"
    )

    norm_count, minimum_norm = norm_controls()
    print(
        f"paired cyclotomic norm controls: {norm_count}, "
        f"min_abs_norm={minimum_norm}, exponent=5, PASS"
    )

    real, imag = composite_hostile()
    print(f"composite C4 anchored hostile: m1=({real},{imag}), PASS")

    irrational_residual, anchored_residual, c91_residual = (
        analytic_hostile_residuals()
    )
    print(
        "analytic hostile residuals: "
        f"C13<{irrational_residual:.3e}, "
        f"anchored<{anchored_residual:.3e}, "
        f"C91<{c91_residual:.3e}, PASS"
    )

    hostile_count, flat_cells = terminal_alignment_hostiles()
    print(
        f"terminal alignment hostiles: {hostile_count}, "
        f"flat_cells={flat_cells}, PASS"
    )

    last_lane, common_core, owner_cell, owner_septimal = lrc_constants()
    print(f"last-lane all-colour floor>{last_lane}")
    print(f"common-core all-colour floor>{common_core}")
    print(f"owner-cell all-colour floor>{owner_cell}")
    print(f"owner C7xC13 floor>{owner_septimal}")
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
