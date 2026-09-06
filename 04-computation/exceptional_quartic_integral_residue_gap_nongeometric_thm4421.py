#!/usr/bin/env python3
"""Exact verifier for the THM-4421 integral residue-gap no-go.

The computation uses only Python's standard library.  It audits the integral
tangent lattices at the retained exceptional-quartic triple from THM-4381 /
THM-4411, the polynomial evaluation formula for the collision period, and
sharp positive and hostile controls.
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd
import sys


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def det_bareiss(matrix: tuple[tuple[int, ...], ...]) -> int:
    """Exact determinant of a square integer matrix."""
    n = len(matrix)
    if n == 0:
        return 1
    require(all(len(row) == n for row in matrix), ("square-matrix", matrix))
    a = [list(row) for row in matrix]
    sign = 1
    previous = 1
    for k in range(n - 1):
        pivot_row = next((r for r in range(k, n) if a[r][k] != 0), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                require(numerator % previous == 0, ("bareiss-division", numerator, previous))
                a[i][j] = numerator // previous
        for i in range(k + 1, n):
            a[i][k] = 0
        previous = pivot
    return sign * a[n - 1][n - 1]


def minor_gcds(matrix: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    """Return determinantal divisors through the rank of an integer matrix."""
    nrows = len(matrix)
    ncols = len(matrix[0])
    divisors: list[int] = []
    for size in range(1, min(nrows, ncols) + 1):
        values: list[int] = []
        for rows in combinations(range(nrows), size):
            for cols in combinations(range(ncols), size):
                sub = tuple(tuple(matrix[i][j] for j in cols) for i in rows)
                value = abs(det_bareiss(sub))
                if value:
                    values.append(value)
        if not values:
            break
        divisors.append(reduce(gcd, values))
    return tuple(divisors)


def smith_nonzero_invariants(matrix: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    divisors = minor_gcds(matrix)
    previous = 1
    invariants: list[int] = []
    for divisor in divisors:
        require(divisor % previous == 0, ("smith-divisibility", divisor, previous))
        invariants.append(divisor // previous)
        previous = divisor
    return tuple(invariants)


def mat_vec(matrix: tuple[tuple[int, ...], ...], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a * b for a, b in zip(row, vector)) for row in matrix)


def left_vec_matrix(vector: tuple[int, ...], matrix: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    return tuple(sum(vector[i] * matrix[i][j] for i in range(len(vector))) for j in range(len(matrix[0])))


def rank_mod_prime(matrix: tuple[tuple[int, ...], ...], prime: int) -> int:
    a = [[entry % prime for entry in row] for row in matrix]
    nrows = len(a)
    ncols = len(a[0])
    rank = 0
    for col in range(ncols):
        pivot = next((r for r in range(rank, nrows) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        inv = pow(a[rank][col], -1, prime)
        a[rank] = [(inv * x) % prime for x in a[rank]]
        for row in range(nrows):
            if row != rank and a[row][col]:
                factor = a[row][col]
                a[row] = [(x - factor * y) % prime for x, y in zip(a[row], a[rank])]
        rank += 1
        if rank == nrows:
            break
    return rank


TANGENTS = ((3, -9), (3, 4), (3, 9))
# L(w)_i = det(t_i,w), with w=(c,e).
WEDGE_MOTION = ((9, 3), (-4, 3), (-9, 3))
ELL = (5, -18, 13)

# n_i + v_i t_i = (c,e).  Columns are c,e,v_-,v_0,v_+ and rows are
# C_-,E_-,C_0,E_0,C_+,E_+.
FULL_MOTION = (
    (1, 0, -3, 0, 0),
    (0, 1, 9, 0, 0),
    (1, 0, 0, -3, 0),
    (0, 1, 0, -4, 0),
    (1, 0, 0, 0, -3),
    (0, 1, 0, 0, -9),
)
FULL_LEFT_RELATION = (15, 5, 24, -18, -39, 13)

# h-values -> (delta C,delta E) on the three branches.
NORMAL_FROM_HVALUES = (
    (2, 0, 0),
    (-2, 0, 0),
    (0, 0, 0),
    (0, 4, 0),
    (0, 0, -2),
    (0, 0, -2),
)


def ell(values: tuple[int, int, int]) -> int:
    return sum(a * b for a, b in zip(ELL, values))


def wedge_image_solution(vector: tuple[int, int, int]) -> tuple[int, int] | None:
    """Exact integral solution to WEDGE_MOTION*(c,e)=vector, if it exists."""
    first, middle, last = vector
    if ell(vector) != 0 or first % 3:
        return None
    c_num = first - last
    e_num = first + last
    if c_num % 18 or e_num % 6:
        return None
    c = c_num // 18
    e = e_num // 6
    if mat_vec(WEDGE_MOTION, (c, e)) != vector:
        return None
    return c, e


def polynomial_values(coefficients: tuple[int, ...]) -> tuple[int, int, int]:
    minus = sum(coef * ((-1) ** degree) for degree, coef in enumerate(coefficients))
    zero = coefficients[0] if coefficients else 0
    plus = sum(coefficients)
    return minus, zero, plus


def parity_aggregates(coefficients: tuple[int, ...]) -> tuple[int, int]:
    odd = sum(coefficients[degree] for degree in range(1, len(coefficients), 2))
    positive_even = sum(coefficients[degree] for degree in range(2, len(coefficients), 2))
    return odd, positive_even


def rational_collision_solution(b: int, k: int) -> tuple[Fraction, ...]:
    """Solution for h-values (b-13k,b,b+5k), assuming ell=0."""
    c = Fraction(-12 * k)
    e = Fraction(4 * b - 16 * k)
    v_minus = Fraction(2 * (7 * k - b), 3)
    v_zero = Fraction(-4 * k)
    v_plus = Fraction(2 * (b - k), 3)
    return c, e, v_minus, v_zero, v_plus


def normal_vector_from_hvalues(values: tuple[int, int, int]) -> tuple[int, ...]:
    return mat_vec(NORMAL_FROM_HVALUES, values)


def full_motion_image_rational(solution: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(sum(Fraction(a) * b for a, b in zip(row, solution)) for row in FULL_MOTION)


def main() -> None:
    sys.stdout.reconfigure(newline="\n")
    checks = 0

    # The primitive tangent relation and its integral saturation index.
    tangent_column_1 = tuple(row[0] for row in WEDGE_MOTION)
    tangent_column_2 = tuple(row[1] for row in WEDGE_MOTION)
    cross = (
        tangent_column_1[1] * tangent_column_2[2] - tangent_column_1[2] * tangent_column_2[1],
        tangent_column_1[2] * tangent_column_2[0] - tangent_column_1[0] * tangent_column_2[2],
        tangent_column_1[0] * tangent_column_2[1] - tangent_column_1[1] * tangent_column_2[0],
    )
    require(cross == tuple(3 * x for x in ELL), ("primitive-cross", cross))
    require(reduce(gcd, map(abs, ELL)) == 1, "primitive-relation")
    require(left_vec_matrix(ELL, WEDGE_MOTION) == (0, 0), "wedge-left-kernel")
    checks += 3

    wedge_divisors = minor_gcds(WEDGE_MOTION)
    wedge_smith = smith_nonzero_invariants(WEDGE_MOTION)
    require(wedge_divisors == (1, 3), ("wedge-divisors", wedge_divisors))
    require(wedge_smith == (1, 3), ("wedge-smith", wedge_smith))
    checks += 2

    # Exact characterization im(L)=ker(ell) intersect {first=0 mod 3}.
    image_box_count = 0
    saturated_box_count = 0
    for first in range(-24, 25):
        for middle in range(-24, 25):
            for last in range(-24, 25):
                vector = (first, middle, last)
                in_saturation = ell(vector) == 0
                solution = wedge_image_solution(vector)
                expected_image = in_saturation and first % 3 == 0
                require((solution is not None) == expected_image,
                        ("wedge-image-membership", vector, solution, expected_image))
                if in_saturation:
                    saturated_box_count += 1
                if solution is not None:
                    image_box_count += 1
                checks += 1

    # The full collision-motion lattice has two 3-torsion factors.
    full_divisors = minor_gcds(FULL_MOTION)
    full_smith = smith_nonzero_invariants(FULL_MOTION)
    require(full_divisors == (1, 1, 1, 3, 9), ("full-divisors", full_divisors))
    require(full_smith == (1, 1, 1, 3, 3), ("full-smith", full_smith))
    require(left_vec_matrix(FULL_LEFT_RELATION, FULL_MOTION) == (0, 0, 0, 0, 0),
            "full-left-kernel")
    require(reduce(gcd, map(abs, FULL_LEFT_RELATION)) == 1, "full-relation-primitive")
    require(left_vec_matrix(FULL_LEFT_RELATION, NORMAL_FROM_HVALUES) == (20, -72, 52),
            "normal-period-row")
    require((20, -72, 52) == tuple(4 * x for x in ELL), "normal-period-scale")
    checks += 6

    # Characteristic three is exactly a bad tangent fibre.
    require(rank_mod_prime(TANGENTS, 3) == 1, "mod-three-tangent-rank")
    require(tuple(tuple(x % 3 for x in row) for row in TANGENTS)
            == ((0, 0), (0, 1), (0, 0)), "mod-three-tangent-rows")
    checks += 2

    # A legal rational target rescaling C -> C/3 saturates the wedge lattice.
    rescaled_wedge_motion = ((9, 1), (-4, 1), (-9, 1))
    require(minor_gcds(rescaled_wedge_motion) == (1, 1), "rescaled-wedge-divisors")
    require(smith_nonzero_invariants(rescaled_wedge_motion) == (1, 1),
            "rescaled-wedge-smith")
    require(left_vec_matrix(ELL, rescaled_wedge_motion) == (0, 0),
            "rescaled-wedge-kernel")
    rescaled_full_motion = (
        (1, 0, -1, 0, 0),
        (0, 1, 9, 0, 0),
        (1, 0, 0, -1, 0),
        (0, 1, 0, -4, 0),
        (1, 0, 0, 0, -1),
        (0, 1, 0, 0, -9),
    )
    require(minor_gcds(rescaled_full_motion) == (1, 1, 1, 1, 1),
            "rescaled-full-divisors")
    require(smith_nonzero_invariants(rescaled_full_motion) == (1, 1, 1, 1, 1),
            "rescaled-full-smith")
    checks += 5

    # Polynomial evaluation formula ell(h)=8*O+18*E.
    formula_count = 0
    for a0 in range(-3, 4):
        for a1 in range(-3, 4):
            for a2 in range(-3, 4):
                for a3 in range(-3, 4):
                    for a4 in range(-3, 4):
                        coefficients = (a0, a1, a2, a3, a4)
                        values = polynomial_values(coefficients)
                        odd, even = parity_aggregates(coefficients)
                        require(ell(values) == 8 * odd + 18 * even,
                                ("parity-formula", coefficients, values, odd, even))
                        formula_count += 1
                        checks += 1

    # Sharp conditional residue gap and the minimal visible cancellation.
    best_gap: int | None = None
    zero_pairs: list[tuple[int, int]] = []
    for odd in range(-12, 13):
        for even in range(-12, 13):
            period = 8 * odd + 18 * even
            if odd % 3:
                require(period != 0, ("residue-nonzero", odd, even))
                require(period % 2 == 0 and period % 3 != 0,
                        ("residue-alphabet", odd, even, period))
                candidate = abs(period)
                if best_gap is None or candidate < best_gap:
                    best_gap = candidate
            if period == 0 and (odd, even) != (0, 0):
                zero_pairs.append((odd, even))
                require(odd % 9 == 0 and even % 4 == 0,
                        ("zero-period-divisibility", odd, even))
                require(odd // 9 == -(even // 4),
                        ("zero-period-parameter", odd, even))
            checks += 1
    require(best_gap == 2, ("sharp-gap", best_gap))
    require(8 * (-2) + 18 * 1 == 2, "sharp-gap-witness")
    primitive_zero_pairs = [
        pair for pair in zero_pairs if gcd(abs(pair[0]), abs(pair[1])) == 1
    ]
    require(primitive_zero_pairs == [(-9, 4), (9, -4)],
            ("primitive-zero-pairs", primitive_zero_pairs))
    require(min(abs(odd) + abs(even) for odd, even in zero_pairs) == 13,
            "minimal-zero-pair-norm")
    checks += 3

    sharp_h = (0, -2, 1)  # x^2-2x
    visible_zero_h = (0, 9, -4)  # 9x-4x^2
    invisible_h = (0, -1, 0, 1)  # x^3-x
    require(polynomial_values(sharp_h) == (3, 0, -1), "sharp-h-values")
    require(ell(polynomial_values(sharp_h)) == 2, "sharp-h-period")
    require(polynomial_values(visible_zero_h) == (-13, 0, 5), "visible-zero-values")
    require(ell(polynomial_values(visible_zero_h)) == 0, "visible-zero-period")
    require(polynomial_values(invisible_h) == (0, 0, 0), "invisible-h-values")
    require(ell(polynomial_values(invisible_h)) == 0, "invisible-h-period")
    checks += 6

    # Every actual wedge vector W_h=12(h_-,h_0,h_+) with ell(h)=0
    # kills the wedge-lattice Z/3 torsion and has an integral common target
    # velocity.  This is checked on a bounded value universe.
    actual_wedge_count = 0
    for h_minus in range(-30, 31):
        for h_zero in range(-30, 31):
            for h_plus in range(-30, 31):
                values = (h_minus, h_zero, h_plus)
                if ell(values):
                    continue
                wedge_vector = tuple(12 * x for x in values)
                solution = wedge_image_solution(wedge_vector)
                require(solution is not None, ("actual-wedge-solution", values))
                require(solution == (
                    2 * (h_minus - h_plus) // 3,
                    2 * (h_minus + h_plus),
                ), ("actual-wedge-formula", values, solution))
                actual_wedge_count += 1
                checks += 1

    # On the polynomial collision kernel O=9k,E=-4k, the full integral
    # collision-motion obstruction is exactly rho=b-k mod 3.  Rational
    # collision always persists, and opposite torsion controls cancel.
    polynomial_kernel_count = 0
    integral_motion_count = 0
    for b in range(-12, 13):
        for k in range(-8, 9):
            coefficients = (b, 9 * k, -4 * k)
            values = polynomial_values(coefficients)
            require(values == (b - 13 * k, b, b + 5 * k),
                    ("kernel-values", b, k, values))
            require(ell(values) == 0, ("kernel-period", b, k))
            solution = rational_collision_solution(b, k)
            image = full_motion_image_rational(solution)
            normal = tuple(Fraction(x) for x in normal_vector_from_hvalues(values))
            require(image == normal, ("rational-motion", b, k, image, normal))
            is_integral = all(x.denominator == 1 for x in solution)
            require(is_integral == ((b - k) % 3 == 0),
                    ("integral-motion-residue", b, k, solution))
            polynomial_kernel_count += 1
            integral_motion_count += int(is_integral)
            checks += 4

    controls = {
        "constant_h=1": (1, 0),
        "visible_mixed_h=9x-4x^2": (0, 1),
        "sum_h=1+9x-4x^2": (1, 1),
    }
    control_rows: list[str] = []
    for label, (b, k) in controls.items():
        residue = (b - k) % 3
        solution = rational_collision_solution(b, k)
        is_integral = all(x.denominator == 1 for x in solution)
        require(is_integral == (residue == 0),
                ("control-integrality", label, residue, solution))
        control_rows.append(
            f"{label}: rho={residue}, integral_motion={is_integral}, solution={solution}"
        )
        checks += 1
    require((controls["constant_h=1"][0] - controls["constant_h=1"][1]) % 3 == 1,
            "constant-control-residue")
    require((controls["visible_mixed_h=9x-4x^2"][0]
             - controls["visible_mixed_h=9x-4x^2"][1]) % 3 == 2,
            "mixed-control-residue")
    require((controls["sum_h=1+9x-4x^2"][0]
             - controls["sum_h=1+9x-4x^2"][1]) % 3 == 0,
            "sum-control-residue")
    checks += 3

    # Pure Newton rays: every nonconstant monomial has period 8 or 18.
    monomial_periods = tuple(8 if degree % 2 else 18 for degree in range(1, 17))
    for degree, expected in enumerate(monomial_periods, 1):
        coefficients = tuple(0 for _ in range(degree)) + (1,)
        require(ell(polynomial_values(coefficients)) == expected,
                ("monomial-period", degree, expected))
        checks += 1

    print("status=PASS")
    print(f"checks={checks}")
    print(f"tangents={TANGENTS}")
    print(f"primitive_relation={ELL}")
    print(f"wedge_determinantal_divisors={wedge_divisors}")
    print(f"wedge_smith_nonzero={wedge_smith}; coker=Z+(Z/3)")
    print(f"full_determinantal_divisors={full_divisors}")
    print(f"full_smith_nonzero={full_smith}; coker=Z+(Z/3)^2")
    print(f"wedge_box_saturation={saturated_box_count}; wedge_box_image={image_box_count}")
    print(f"formula_polynomials_checked={formula_count}")
    print(f"sharp_residue_gap={best_gap}")
    print(f"primitive_visible_zero_pairs={primitive_zero_pairs}")
    print("sharp_split_h=x^2-2x; values=(3,0,-1); ell=2")
    print("visible_collision_h=9x-4x^2; values=(-13,0,5); ell=0")
    print(f"actual_wedge_kernel_vectors_checked={actual_wedge_count}")
    print(
        "polynomial_collision_kernel_checked="
        f"{polynomial_kernel_count}; integral_motion={integral_motion_count}"
    )
    for row in control_rows:
        print(row)
    print(f"monomial_periods_degrees_1_to_16={monomial_periods}")
    print("mod3_tangent_rank=1; endpoint_tangents_zero_mod3")
    print("after_C_to_C_over_3_wedge_smith_nonzero=(1,1)")
    print("after_C_to_C_over_3_full_smith_nonzero=(1,1,1,1,1)")
    print("verdict=conditional_residue_gap_plus_unconditional_transfer_no_go")


if __name__ == "__main__":
    main()
