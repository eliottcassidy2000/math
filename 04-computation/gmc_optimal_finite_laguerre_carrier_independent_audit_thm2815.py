"""Independent hostile audit for THM-2815.

This companion deliberately avoids SymPy and every other CAS.  It rebuilds
the Laguerre relation, factorial Hankel inverse, determinant, quotient
reduction, unique optimal relation, and MISTAKE-211 selector recovery using
only Python integers and fractions.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def trim(poly: list[Fraction]) -> list[Fraction]:
    answer = poly[:]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def multiply(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    answer = [Fraction(0) for _ in range(len(left) + len(right) - 1)]
    for i, left_coefficient in enumerate(left):
        for j, right_coefficient in enumerate(right):
            answer[i + j] += left_coefficient * right_coefficient
    return trim(answer)


def factorial_functional(poly: list[Fraction]) -> Fraction:
    return sum(
        (coefficient * factorial(degree) for degree, coefficient in enumerate(poly)),
        Fraction(0),
    )


def monic_laguerre(degree: int) -> list[Fraction]:
    return [
        Fraction(
            (-1) ** (degree + k) * factorial(degree) * comb(degree, k),
            factorial(k),
        )
        for k in range(degree + 1)
    ]


def monic_remainder(
    dividend: list[Fraction], divisor: list[Fraction]
) -> list[Fraction]:
    """Return the remainder on division by a monic polynomial."""
    require(divisor[-1] == 1, "divisor is not monic")
    work = trim(dividend)
    divisor_degree = len(divisor) - 1
    while len(work) - 1 >= divisor_degree:
        shift = len(work) - 1 - divisor_degree
        leading = work[-1]
        for j, coefficient in enumerate(divisor):
            work[shift + j] -= leading * coefficient
        work = trim(work)
    return work


def hankel(size: int) -> list[list[Fraction]]:
    return [
        [Fraction(factorial(i + j)) for j in range(size)]
        for i in range(size)
    ]


def inverse_entry(last_degree: int, i: int, j: int) -> Fraction:
    return Fraction(
        (-1) ** (i + j)
        * sum(
            comb(k, i) * comb(k, j)
            for k in range(max(i, j), last_degree + 1)
        ),
        factorial(i) * factorial(j),
    )


def inverse_hankel(last_degree: int) -> list[list[Fraction]]:
    return [
        [
            inverse_entry(last_degree, i, j)
            for j in range(last_degree + 1)
        ]
        for i in range(last_degree + 1)
    ]


def matrix_product(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    rows = len(left)
    middle = len(right)
    columns = len(right[0])
    require(len(left[0]) == middle, "matrix dimensions do not match")
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(middle)),
                Fraction(0),
            )
            for j in range(columns)
        ]
        for i in range(rows)
    ]


def identity(size: int) -> list[list[Fraction]]:
    return [
        [Fraction(int(i == j)) for j in range(size)]
        for i in range(size)
    ]


def bareiss_determinant(matrix: list[list[int]]) -> int:
    """Fraction-free determinant; the factorial Hankel pivots are nonzero."""
    size = len(matrix)
    if size == 0:
        return 1
    work = [row[:] for row in matrix]
    previous = 1
    for pivot_index in range(size - 1):
        pivot = work[pivot_index][pivot_index]
        require(pivot != 0, f"zero Bareiss pivot at {pivot_index}")
        for i in range(pivot_index + 1, size):
            for j in range(pivot_index + 1, size):
                numerator = (
                    work[i][j] * pivot
                    - work[i][pivot_index] * work[pivot_index][j]
                )
                require(
                    numerator % previous == 0,
                    f"nonexact Bareiss division at ({i},{j})",
                )
                work[i][j] = numerator // previous
        previous = pivot
    return work[-1][-1]


def coefficient_selector(
    last_degree: int, selected_degree: int
) -> list[Fraction]:
    inverse = inverse_hankel(last_degree)
    return [inverse[i][selected_degree] for i in range(last_degree + 1)]


laguerre_count = 0
for degree in range(1, 15):
    ell = monic_laguerre(degree)
    require(ell[-1] == 1, f"monicity failed at D={degree}")

    orthogonality = [
        factorial_functional(
            [Fraction(0)] * shift + ell
        )
        for shift in range(degree)
    ]
    require(
        orthogonality == [Fraction(0)] * degree,
        f"orthogonality failed at D={degree}",
    )

    first_detected_multiplier = factorial_functional(
        [Fraction(0)] * degree + ell
    )
    require(
        first_detected_multiplier == factorial(degree) ** 2,
        f"first multiplier failed at D={degree}",
    )
    require(
        factorial_functional(multiply(ell, ell))
        == factorial(degree) ** 2,
        f"squared norm failed at D={degree}",
    )

    # Monomials are a basis, so these checks audit the full linear horizon.
    for exponent in range(2 * degree):
        monomial = [Fraction(0)] * exponent + [Fraction(1)]
        remainder = monic_remainder(monomial, ell)
        require(
            len(remainder) <= degree,
            f"remainder degree failed at D={degree}, n={exponent}",
        )
        require(
            factorial_functional(remainder) == factorial(exponent),
            f"quotient readout failed at D={degree}, n={exponent}",
        )
    first_alias_remainder = monic_remainder(multiply(ell, ell), ell)
    require(
        first_alias_remainder == [Fraction(0)],
        f"square does not vanish in quotient at D={degree}",
    )
    for leading_coefficient in (Fraction(-3), Fraction(2)):
        multiplier = [
            leading_coefficient * coefficient for coefficient in ell
        ]
        for i in range(degree):
            multiplier[i] += Fraction((i + 1) * (degree + 2))
        first_layer_alias = multiply(multiplier, ell)
        require(
            len(first_layer_alias) - 1 == 2 * degree,
            f"first alias layer degree failed at D={degree}",
        )
        require(
            factorial_functional(first_layer_alias)
            == leading_coefficient * factorial(degree) ** 2,
            f"first alias layer readout failed at D={degree}",
        )

    # Solve the monic orthogonality relation independently with H^{-1}.
    H_inverse = inverse_hankel(degree - 1)
    right_side = [
        Fraction(-factorial(degree + row))
        for row in range(degree)
    ]
    lower_coefficients = [
        sum(
            (
                H_inverse[i][row] * right_side[row]
                for row in range(degree)
            ),
            Fraction(0),
        )
        for i in range(degree)
    ]
    require(
        lower_coefficients + [Fraction(1)] == ell,
        f"unique monic relation failed at D={degree}",
    )

    # One moment earlier, uniqueness fails by exactly one Laguerre direction.
    ell_previous = monic_laguerre(degree - 1)
    for parameter in (Fraction(-3), Fraction(-1), Fraction(0), Fraction(2), Fraction(5)):
        shortened_relation = ell[:]
        for i, coefficient in enumerate(ell_previous):
            shortened_relation[i] += parameter * coefficient
        require(
            shortened_relation[-1] == 1,
            f"short-horizon relation lost monicity at D={degree}",
        )
        for exponent in range(2 * degree - 1):
            monomial = [Fraction(0)] * exponent + [Fraction(1)]
            remainder = monic_remainder(monomial, shortened_relation)
            require(
                factorial_functional(remainder) == factorial(exponent),
                (
                    "short-horizon quotient failed at "
                    f"D={degree}, c={parameter}, n={exponent}"
                ),
            )
        final_probe = factorial_functional(
            [Fraction(0)] * (degree - 1) + shortened_relation
        )
        require(
            final_probe == parameter * factorial(degree - 1) ** 2,
            f"uniqueness-lock moment failed at D={degree}, c={parameter}",
        )
    laguerre_count += 1


selector_count = 0
for last_degree in range(0, 13):
    size = last_degree + 1
    H = hankel(size)
    H_inverse = inverse_hankel(last_degree)
    require(
        matrix_product(H, H_inverse) == identity(size),
        f"closed inverse formula failed at d={last_degree}",
    )
    require(
        matrix_product(H_inverse, H) == identity(size),
        f"reverse inverse formula failed at d={last_degree}",
    )

    determinant = bareiss_determinant(
        [[factorial(i + j) for j in range(size)] for i in range(size)]
    )
    expected_determinant = 1
    for j in range(size):
        expected_determinant *= factorial(j) ** 2
    require(
        determinant == expected_determinant,
        f"factorial Hankel determinant failed at d={last_degree}",
    )

    for selected_degree in range(size):
        phi = coefficient_selector(last_degree, selected_degree)
        laguerre_form = [Fraction(0) for _ in range(size)]
        for k in range(selected_degree, size):
            weight = Fraction(
                (-1) ** selected_degree * comb(k, selected_degree),
                factorial(selected_degree),
            )
            for i in range(k + 1):
                standard_laguerre_coefficient = Fraction(
                    (-1) ** i * comb(k, i), factorial(i)
                )
                laguerre_form[i] += weight * standard_laguerre_coefficient
        require(
            laguerre_form == phi,
            (
                "Laguerre-basis selector formula failed at "
                f"d={last_degree}, j={selected_degree}"
            ),
        )
        dual_values = [
            factorial_functional([Fraction(0)] * row + phi)
            for row in range(size)
        ]
        require(
            dual_values
            == [
                Fraction(int(row == selected_degree))
                for row in range(size)
            ],
            (
                "selector duality failed at "
                f"d={last_degree}, j={selected_degree}"
            ),
        )
        norm = factorial_functional(multiply(phi, phi))
        expected_norm = Fraction(
            sum(
                comb(k, selected_degree) ** 2
                for k in range(selected_degree, size)
            ),
            factorial(selected_degree) ** 2,
        )
        require(
            norm == expected_norm,
            f"selector norm failed at d={last_degree}, j={selected_degree}",
        )
        selector_count += 1


# Rebuild the MISTAKE-211 hostile without importing the primary companion.
mistake_coefficient_6 = Fraction(4)
mistake_coefficient_18 = Fraction(-4 * factorial(6), factorial(18))
G4 = [Fraction(0)] * 19
G4[6] = mistake_coefficient_6
G4[18] = mistake_coefficient_18
require(
    factorial_functional(G4) == 0,
    "MISTAKE-211 scalar cancellation failed",
)
recovered = []
for selected_degree in (6, 18):
    phi = coefficient_selector(18, selected_degree)
    recovered.append(factorial_functional(multiply(G4, phi)))
require(
    recovered == [mistake_coefficient_6, mistake_coefficient_18],
    "MISTAKE-211 multiplier selector failed",
)
require(
    all(value != 0 for value in recovered),
    "MISTAKE-211 hostile lost a private coefficient",
)


print("THM-2815 INDEPENDENT HOSTILE AUDIT")
print("engine=python-stdlib integers+fractions; no CAS")
print(f"laguerre_D_range=1..{laguerre_count}")
print("exact_horizon=2D-1; first_readout_alias_degree=2D")
print("first_alias_value=(D!)^2")
print("first_alias_layer=L(q*ell_D)=lc(q)*(D!)^2 for deg(q)=D")
print("optimal_relation=unique monic factorial-orthogonal degree-D relation")
print("uniqueness_lock=ell_D+c*ell_(D-1) survives 2D-2; moment 2D-1 kills c")
print("all_horizons=min_dimension=floor(H/2)+1; even=affine-line; odd=unique")
print(f"selector_d_range=0..12; selector_columns={selector_count}")
print("hankel_checks=two-sided_inverse+bareiss_determinant")
print(f"mistake211_scalar={factorial_functional(G4)}")
print(f"mistake211_private_coefficients={tuple(recovered)}")
print(
    "scope=selectors require extra L(G*phi_j) observations; "
    "scalar nullity supplies only L(G)"
)
print("ALL INDEPENDENT HOSTILE CHECKS PASSED")
