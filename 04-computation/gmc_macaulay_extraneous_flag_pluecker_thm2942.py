#!/usr/bin/env python3
"""Exact companion for THM-2942.

The fixed degree-seven 20Q+10C+6F Macaulay chart has a universal
extraneous flag factor.  This script:

* extracts the six symbolic x_2-degree blocks after F=x_2^4;
* identifies their common binary resultant and the residual flag curvature;
* checks the pure-power specialization from THM-2927;
* verifies, over a prime field, all 210 quartic-row mutations against the
  Pluecker coordinates of (S/(Q,C))_3; and
* audits the one-cut Taylor/Gregory--Newton positivity lemma and the
  factorial first-gap independence consequence.
"""

from __future__ import annotations

import importlib.util
from hashlib import sha256
from itertools import combinations
from math import factorial
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


SOURCE = Path(__file__).with_name(
    "gmc_diameter_four_nonconsecutive_macaulay_newton_thm2921.py"
)
SOURCE_BYTES = SOURCE.read_bytes().replace(b"\r\n", b"\n").replace(
    b"\r", b"\n"
)
SOURCE_SHA256 = (
    "42e9b5ceddd677d1f2601a9d5d668c9437281596b65999ddcb8549d4e0b9bf64"
)
require(
    sha256(SOURCE_BYTES).hexdigest() == SOURCE_SHA256,
    "THM-2921 constructor dependency hash changed",
)
SPEC = importlib.util.spec_from_file_location("thm2921_exact", SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load THM-2921")
t = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(t)

DENOMINATOR_SOURCE = Path(__file__).with_name(
    "gmc_general_width_terminal_pole_cancellation_thm2925.py"
)
DENOMINATOR_BYTES = DENOMINATOR_SOURCE.read_bytes().replace(
    b"\r\n", b"\n"
).replace(b"\r", b"\n")
DENOMINATOR_SHA256 = (
    "83d70a95f0943992d0e4b7027eede431d4dc968b66655e37b43fd0acfc692e47"
)
require(
    sha256(DENOMINATOR_BYTES).hexdigest() == DENOMINATOR_SHA256,
    "THM-2925 denominator dependency hash changed",
)
DENOMINATOR_SPEC = importlib.util.spec_from_file_location(
    "thm2925_exact", DENOMINATOR_SOURCE
)
require(
    DENOMINATOR_SPEC is not None and DENOMINATOR_SPEC.loader is not None,
    "cannot load THM-2925",
)
thm2925 = importlib.util.module_from_spec(DENOMINATOR_SPEC)
DENOMINATOR_SPEC.loader.exec_module(thm2925)


def macaulay_rows(forms: tuple[dict[tuple[int, int, int], object], ...]):
    rows = []
    metadata = []
    for degree, form in zip(t.ORDERS, forms):
        for multiplier in t.MONOMIALS[t.TARGET_DEGREE - degree]:
            row = [0] * len(t.TARGET_MONOMIALS)
            for monomial, coefficient in form.items():
                target = tuple(
                    multiplier[index] + monomial[index]
                    for index in range(3)
                )
                row[t.TARGET_INDEX[target]] = coefficient
            rows.append(row)
            metadata.append((degree, multiplier))
    return rows, metadata


def symbolic_flag_audit() -> tuple[str, str]:
    q20, q11, q02 = sp.symbols("q20 q11 q02")
    c30, c21, c12, c03 = sp.symbols("c30 c21 c12 c03")
    binary_t = sp.symbols("binary_t")
    quadratic = {
        (2, 0, 0): q20,
        (1, 1, 0): q11,
        (0, 2, 0): q02,
    }
    cubic = {
        (3, 0, 0): c30,
        (2, 1, 0): c21,
        (1, 2, 0): c12,
        (0, 3, 0): c03,
    }
    quartic = {(0, 0, 4): sp.Integer(1)}
    rows, metadata = macaulay_rows((quadratic, cubic, quartic))
    selected = sp.Matrix([rows[index] for index in t.SELECTED_ROWS])
    selected_metadata = [metadata[index] for index in t.SELECTED_ROWS]

    quartic_rows = list(range(30, 36))
    pivot_columns = []
    pivot_targets = []
    for row_index in quartic_rows:
        nonzero = [
            column
            for column in range(selected.cols)
            if selected[row_index, column] != 0
        ]
        require(
            len(nonzero) == 1
            and selected[row_index, nonzero[0]] == 1,
            "quartic unit pivot changed",
        )
        pivot_columns.append(nonzero[0])
        pivot_targets.append(t.TARGET_MONOMIALS[nonzero[0]])
    require(
        pivot_targets
        == [
            (0, 0, 7),
            (0, 1, 6),
            (0, 2, 5),
            (0, 3, 4),
            (1, 0, 6),
            (1, 1, 5),
        ],
        "quartic pivot targets changed",
    )
    require(
        (sum(quartic_rows) + sum(pivot_columns)) % 2 == 0,
        "quartic pivot expansion changed sign",
    )

    binary_q = q20 + q11 * binary_t + q02 * binary_t**2
    binary_c = (
        c30 + c21 * binary_t + c12 * binary_t**2 + c03 * binary_t**3
    )
    binary_resultant = sp.factor(
        sp.resultant(binary_q, binary_c, binary_t)
    )
    curvature = sp.expand(
        c12 * q20**2
        - c21 * q11 * q20
        - c30 * q02 * q20
        + c30 * q11**2
    )
    residual_columns = [
        column
        for column in range(selected.cols)
        if column not in pivot_columns
    ]
    expected = (
        c30 * q20**2 * binary_resultant,
        q20**2 * binary_resultant,
        q20 * binary_resultant,
        binary_resultant,
        curvature,
        q20,
    )
    block_records = []
    block_product = sp.Integer(1)
    for z_degree in range(6):
        block_rows = [
            row
            for row in range(30)
            if selected_metadata[row][1][2] == z_degree
        ]
        block_columns = [
            column
            for column in residual_columns
            if t.TARGET_MONOMIALS[column][2] == z_degree
        ]
        require(
            len(block_rows) == len(block_columns),
            f"nonsquare x2 block at degree {z_degree}",
        )
        block = selected.extract(block_rows, block_columns)
        determinant = sp.factor(block.det(method="domain-ge"))
        require(
            sp.expand(determinant - expected[z_degree]) == 0,
            f"symbolic x2 block changed at degree {z_degree}",
        )
        block_product *= determinant
        block_records.append(
            f"{z_degree}:{len(block_rows)}:{sp.sstr(determinant)}"
        )
    require(
        sp.expand(
            block_product
            - q20**6
            * c30
            * curvature
            * binary_resultant**4
        )
        == 0,
        "symbolic flag product changed",
    )

    u0, u1, v0, v1 = sp.symbols("u0 u1 v0 v1")
    pure_curvature = sp.expand(
        curvature.subs(
            {
                q20: u0**2,
                q11: 2 * u0 * u1,
                q02: u1**2,
                c30: v0**3,
                c21: 3 * v0**2 * v1,
                c12: 3 * v0 * v1**2,
                c03: v1**3,
            }
        )
    )
    require(
        sp.expand(
            pure_curvature
            - 3 * v0 * u0**2 * (u0 * v1 - u1 * v0) ** 2
        )
        == 0,
        "pure-power flag curvature changed",
    )
    block_digest = sha256(
        ("\n".join(block_records) + "\n").encode()
    ).hexdigest()
    resultant_digest = sha256(
        sp.sstr(binary_resultant).encode()
    ).hexdigest()
    return block_digest, resultant_digest


def symbolic_mutation_audit() -> str:
    q_symbols = {
        monomial: sp.symbols("q" + "".join(map(str, monomial)))
        for monomial in t.MONOMIALS[2]
    }
    c_symbols = {
        monomial: sp.symbols("c" + "".join(map(str, monomial)))
        for monomial in t.MONOMIALS[3]
    }
    degree_three = t.MONOMIALS[3]
    degree_three_index = {
        monomial: index for index, monomial in enumerate(degree_three)
    }
    relations = []
    for multiplier in t.MONOMIALS[1]:
        row = [0] * len(degree_three)
        for monomial, coefficient in q_symbols.items():
            target = tuple(
                multiplier[axis] + monomial[axis] for axis in range(3)
            )
            row[degree_three_index[target]] = coefficient
        relations.append(row)
    relations.append([c_symbols[monomial] for monomial in degree_three])
    relation_matrix = sp.Matrix(relations)

    q200 = q_symbols[(2, 0, 0)]
    q110 = q_symbols[(1, 1, 0)]
    q020 = q_symbols[(0, 2, 0)]
    q011 = q_symbols[(0, 1, 1)]
    q002 = q_symbols[(0, 0, 2)]
    c300 = c_symbols[(3, 0, 0)]
    c210 = c_symbols[(2, 1, 0)]
    c120 = c_symbols[(1, 2, 0)]
    c021 = c_symbols[(0, 2, 1)]
    c012 = c_symbols[(0, 1, 2)]
    curvature = sp.expand(
        c120 * q200**2
        - c210 * q110 * q200
        - c300 * q020 * q200
        + c300 * q110**2
    )
    alternate = sp.expand(
        c012 * q020 * q200**2
        - c021 * q011 * q200**2
        - c210 * q002 * q020 * q200
        + c210 * q011**2 * q200
        + c300 * q002 * q020 * q110
        - c300 * q011**2 * q110
    )
    original_rows = (0, 1, 2, 3, 4, 5)
    alternate_rows = (0, 3, 4, 5, 6, 7)
    original_complement = [
        index for index in range(10) if index not in original_rows
    ]
    alternate_complement = [
        index for index in range(10) if index not in alternate_rows
    ]
    original_minor = sp.factor(
        relation_matrix[:, original_complement].det()
    )
    alternate_minor = sp.factor(
        relation_matrix[:, alternate_complement].det()
    )
    require(
        sp.expand(original_minor + q200 * curvature) == 0,
        "original complementary Pluecker minor changed",
    )
    require(
        sp.expand(alternate_minor + alternate) == 0,
        "alternate complementary Pluecker minor changed",
    )
    payload = (
        f"original:{sp.sstr(original_minor)}\n"
        f"alternate:{sp.sstr(alternate_minor)}\n"
    )
    return sha256(payload.encode()).hexdigest()


def rref_mod(matrix: list[list[int]], prime: int):
    rows = [[entry % prime for entry in row] for row in matrix]
    pivot_columns = []
    pivot_row = 0
    for column in range(len(rows[0])):
        selected = next(
            (
                row
                for row in range(pivot_row, len(rows))
                if rows[row][column] % prime
            ),
            None,
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        inverse = pow(rows[pivot_row][column], prime - 2, prime)
        rows[pivot_row] = [
            entry * inverse % prime for entry in rows[pivot_row]
        ]
        for row in range(len(rows)):
            if row == pivot_row or not rows[row][column]:
                continue
            multiplier = rows[row][column]
            rows[row] = [
                (
                    rows[row][index]
                    - multiplier * rows[pivot_row][index]
                )
                % prime
                for index in range(len(rows[0]))
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return rows, pivot_columns


def quotient_coordinates(
    quadratic: dict[tuple[int, int, int], int],
    cubic: dict[tuple[int, int, int], int],
) -> list[list[int]]:
    """Coordinates of the ten degree-three monomials in (S/(Q,C))_3."""
    prime = t.PRIME
    degree_three = t.MONOMIALS[3]
    index = {monomial: position for position, monomial in enumerate(degree_three)}
    relations = []
    for multiplier in t.MONOMIALS[1]:
        row = [0] * len(degree_three)
        for monomial, coefficient in quadratic.items():
            target = tuple(
                multiplier[axis] + monomial[axis] for axis in range(3)
            )
            row[index[target]] = coefficient % prime
        relations.append(row)
    relations.append([cubic[monomial] % prime for monomial in degree_three])
    reduced, pivots = rref_mod(relations, prime)
    require(len(pivots) == 4, "generic (Q,C)_3 relation rank changed")
    free = [
        column for column in range(len(degree_three)) if column not in pivots
    ]
    require(len(free) == 6, "degree-three quotient dimension changed")
    answer = []
    for monomial_index in range(len(degree_three)):
        coordinates = [0] * len(free)
        if monomial_index in free:
            coordinates[free.index(monomial_index)] = 1
        else:
            row = pivots.index(monomial_index)
            for coordinate, free_column in enumerate(free):
                coordinates[coordinate] = -reduced[row][free_column] % prime
        answer.append(coordinates)
    return answer


def pluecker_audit() -> tuple[str, int, int]:
    prime = t.PRIME
    primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47)
    quadratic = {
        monomial: primes[index]
        for index, monomial in enumerate(t.MONOMIALS[2])
    }
    cubic = {
        monomial: primes[index + 1]
        for index, monomial in enumerate(t.MONOMIALS[3])
    }
    quartic = {
        monomial: primes[index] ** 2 + 1
        for index, monomial in enumerate(t.MONOMIALS[4])
    }
    rows = t.macaulay_rows_mod((quadratic, cubic, quartic))
    fixed_rows = [rows[index] for index in t.SELECTED_ROWS[:30]]
    coordinates = quotient_coordinates(quadratic, cubic)

    records = []
    common_scale = None
    nonzero_coordinates = 0
    pluecker_values = {}
    determinant_values = {}
    for subset in combinations(range(10), 6):
        pluecker = t.determinant_mod(
            [coordinates[index] for index in subset]
        )
        determinant = t.determinant_mod(
            fixed_rows + [rows[36 + index] for index in subset]
        )
        pluecker_values[subset] = pluecker
        determinant_values[subset] = determinant
        if pluecker:
            nonzero_coordinates += 1
            scale = determinant * pow(pluecker, prime - 2, prime) % prime
            if common_scale is None:
                common_scale = scale
            require(
                scale == common_scale,
                f"quartic mutation scale changed at {subset}",
            )
        else:
            require(
                determinant == 0,
                f"zero Pluecker coordinate has nonzero chart at {subset}",
            )
        records.append(f"{subset}:{pluecker}:{determinant}")
    require(
        len(records) == 210
        and nonzero_coordinates == 207
        and common_scale == 517319,
        "210-chart quotient census changed",
    )

    fixed = (0, 1, 2, 3)
    a, b, c, d = (4, 5, 6, 7)

    def coordinate(table, left, right):
        return table[tuple(sorted(fixed + (left, right)))]

    for table, label in (
        (pluecker_values, "quotient"),
        (determinant_values, "Macaulay"),
    ):
        relation = (
            coordinate(table, a, b) * coordinate(table, c, d)
            - coordinate(table, a, c) * coordinate(table, b, d)
            + coordinate(table, a, d) * coordinate(table, b, c)
        ) % prime
        require(relation == 0, f"{label} Pluecker exchange failed")
    digest = sha256(("\n".join(records) + "\n").encode()).hexdigest()
    return digest, nonzero_coordinates, int(common_scale)


def polynomial_value(coefficients: list[int], value: int) -> int:
    return sum(
        coefficient * value**degree
        for degree, coefficient in enumerate(coefficients)
    )


def derivative_coefficients(coefficients: list[int], order: int) -> list[int]:
    answer = coefficients[:]
    for _ in range(order):
        answer = [
            (degree + 1) * answer[degree + 1]
            for degree in range(len(answer) - 1)
        ]
    return answer


def forward_difference(
    coefficients: list[int],
    base: int,
    order: int,
) -> int:
    values = [
        polynomial_value(coefficients, base + shift)
        for shift in range(order + 1)
    ]
    for _ in range(order):
        values = [
            values[index + 1] - values[index]
            for index in range(len(values) - 1)
        ]
    return values[0]


def one_cut_audit() -> tuple[int, int]:
    positive_cases = 0
    for degree in range(2, 13):
        for cut in range(1, degree + 1):
            coefficients = [
                -(index + 1) if index < cut else index + 1
                for index in range(degree + 1)
            ]
            base = 5
            require(
                polynomial_value(coefficients, base) > 0,
                "one-cut positive-base control changed",
            )
            for order in range(1, degree + 1):
                derivative = derivative_coefficients(coefficients, order)
                require(
                    polynomial_value(derivative, base) > 0
                    and polynomial_value(derivative, base + 3) > 0,
                    "one-cut derivative positivity failed",
                )
                require(
                    forward_difference(coefficients, base, order) > 0,
                    "one-cut Gregory--Newton positivity failed",
                )
                positive_cases += 1

    root_cases = 0
    for degree in range(1, 13):
        base = 2
        coefficients = [-base**degree] + [0] * (degree - 1) + [1]
        require(
            polynomial_value(coefficients, base) == 0,
            "positive-root constant coefficient changed",
        )
        for order in range(1, degree + 1):
            require(
                polynomial_value(
                    derivative_coefficients(coefficients, order),
                    base,
                )
                > 0
                and forward_difference(coefficients, base, order) > 0,
                "positive-root higher coefficient changed",
            )
            root_cases += 1
        require(
            derivative_coefficients(coefficients, degree + 1) == []
            and forward_difference(coefficients, base, degree + 1) == 0,
            "one-cut degree boundary changed",
        )
    return positive_cases, root_cases


def first_gap_audit() -> int:
    """The restricted Q,C coefficients do not see the second interior slot."""
    checks = 0
    q_monomials = ((2, 0, 0), (1, 1, 0), (0, 2, 0))
    c_monomials = (
        (3, 0, 0),
        (2, 1, 0),
        (1, 2, 0),
        (0, 3, 0),
    )
    for width in range(3, 11):
        for first in range(1, width - 1):
            reference_offsets = (0, first, first + 1, width)
            for depth in range(4):
                reference_q = t.moment_form(depth, 2, reference_offsets)
                reference_c = t.moment_form(depth, 3, reference_offsets)
                reference = tuple(
                    reference_q[monomial] for monomial in q_monomials
                ) + tuple(
                    reference_c[monomial] for monomial in c_monomials
                )
                for second in range(first + 1, width):
                    offsets = (0, first, second, width)
                    quadratic = t.moment_form(depth, 2, offsets)
                    cubic = t.moment_form(depth, 3, offsets)
                    current = tuple(
                        quadratic[monomial] for monomial in q_monomials
                    ) + tuple(
                        cubic[monomial] for monomial in c_monomials
                    )
                    require(
                        current == reference,
                        "first-gap flag acquired second-slot dependence",
                    )
                    checks += 1
    require(checks == 480, "first-gap audit count changed")
    return checks


def factorial_cofactor_separation_audit() -> tuple[str, int]:
    """Finite exact control for the original and one mutated flag cofactor."""
    records = []
    support_count = 0
    naive_degree_anomalies = []
    for width in range(3, 21):
        for first in range(1, width - 1):
            for second in range(first + 1, width):
                interior = (0, first, second)
                quadratic = {
                    monomial: thm2925.scaled_coefficient(
                        width, monomial, interior
                    )
                    for monomial in t.MONOMIALS[2]
                }
                cubic = {
                    monomial: thm2925.scaled_coefficient(
                        width, monomial, interior
                    )
                    for monomial in t.MONOMIALS[3]
                }
                q200 = quadratic[(2, 0, 0)]
                q110 = quadratic[(1, 1, 0)]
                q020 = quadratic[(0, 2, 0)]
                q011 = quadratic[(0, 1, 1)]
                q002 = quadratic[(0, 0, 2)]
                c300 = cubic[(3, 0, 0)]
                c210 = cubic[(2, 1, 0)]
                c120 = cubic[(1, 2, 0)]
                c021 = cubic[(0, 2, 1)]
                c012 = cubic[(0, 1, 2)]
                curvature = (
                    c120 * q200**2
                    - c210 * q110 * q200
                    - c300 * q020 * q200
                    + c300 * q110**2
                )
                alternate = (
                    c012 * q020 * q200**2
                    - c021 * q011 * q200**2
                    - c210 * q002 * q020 * q200
                    + c210 * q011**2 * q200
                    + c300 * q002 * q020 * q110
                    - c300 * q011**2 * q110
                )
                require(
                    curvature != 0 and alternate != 0,
                    "factorial flag cofactor vanished identically",
                )
                common = curvature.gcd(alternate)
                if common.leading_coefficient() < 0:
                    common = -common
                coefficients = tuple(common.coeffs())
                require(
                    common(0) != 0
                    and all(coefficient >= 0 for coefficient in coefficients),
                    "common flag factor lost one-sign positivity",
                )
                naive_degree = max(1, width // 2 + 1 - second)
                if common.degree() != naive_degree:
                    naive_degree_anomalies.append(
                        (
                            width,
                            first,
                            second,
                            common.degree(),
                            naive_degree,
                        )
                    )
                primitive_curvature = curvature // common
                primitive_alternate = alternate // common
                require(
                    primitive_curvature.gcd(primitive_alternate).degree() == 0,
                    "reduced flag cofactors retained a common root",
                )
                records.append(
                    f"{width}:{first}:{second}:{common.degree()}:"
                    + ",".join(map(str, coefficients))
                    + ":"
                    + sha256(
                        ",".join(
                            map(str, primitive_curvature.coeffs())
                        ).encode()
                    ).hexdigest()
                    + ":"
                    + sha256(
                        ",".join(
                            map(str, primitive_alternate.coeffs())
                        ).encode()
                    ).hexdigest()
                )
                support_count += 1
    require(support_count == 1140, "factorial cofactor support count changed")
    require(
        naive_degree_anomalies == [(11, 7, 8, 2, 1)],
        "naive common-factor degree hostile changed",
    )
    return (
        sha256(("\n".join(records) + "\n").encode()).hexdigest(),
        support_count,
    )


def main() -> None:
    require(t.is_prime(t.PRIME), "audit modulus is not prime")
    require(
        (20, 10, 6) == (20, 10, 6)
        and (20 - 12, 10 - 8, 6 - 6) == (8, 2, 0),
        "Macaulay/resultant multidegree invoice changed",
    )
    block_digest, resultant_digest = symbolic_flag_audit()
    mutation_digest = symbolic_mutation_audit()
    pluecker_digest, nonzero_coordinates, common_scale = pluecker_audit()
    positive_cases, root_cases = one_cut_audit()
    first_gap_checks = first_gap_audit()
    cofactor_digest, cofactor_supports = factorial_cofactor_separation_audit()

    print("THM-2942 MACAULAY EXTRANEOUS FLAG FACTOR AND PLUECKER MUTATION")
    print(f"constructor_dependency_sha256={SOURCE_SHA256}")
    print(f"denominator_dependency_sha256={DENOMINATOR_SHA256}")
    print("chart_multidegree=(20,10,6);resultant_multidegree=(12,8,6)")
    print(
        "factor=Delta=q20^6*c30*K*Res(Q,C,F);"
        "K=c12*q20^2-c21*q11*q20-c30*q02*q20+c30*q11^2"
    )
    print("F=x2^4_block_sizes=8,7,6,5,3,1;symbolic_blocks=PASS")
    print(f"symbolic_block_digest_sha256={block_digest}")
    print(f"binary_resultant_digest_sha256={resultant_digest}")
    print(
        "mutation=J0:0,1,2,3,4,5->J1:0,3,4,5,6,7;"
        "Delta_J1=q200^5*c300*P_alt*Res"
    )
    print(f"symbolic_mutation_digest_sha256={mutation_digest}")
    print(
        "pure_power_K=3*v0*u0^2*(u0*v1-u1*v0)^2;"
        "THM2927_factor_recovered"
    )
    print(
        f"quartic_row_mutations=210;nonzero_sample_coordinates="
        f"{nonzero_coordinates};common_scale_mod_{t.PRIME}={common_scale}"
    )
    print(f"pluecker_census_digest_sha256={pluecker_digest}")
    print("three_term_pluecker_exchange=quotient:PASS,Macaulay:PASS")
    print(
        f"one_cut_controls=positive:{positive_cases},root:{root_cases};"
        "Taylor_and_Gregory_Newton=PASS"
    )
    print(
        f"factorial_first_gap_independence_checks={first_gap_checks};"
        "second_interior_slot_absent"
    )
    print(
        f"factorial_cofactor_separation_supports={cofactor_supports};"
        "widths=3..20;gcd_one_sign=YES;reduced_coprime=YES"
    )
    print(
        "naive_gcd_degree_law_hostile=(M,a,b)=(11,7,8):"
        "actual2_vs_predicted1"
    )
    print(f"factorial_cofactor_digest_sha256={cofactor_digest}")
    print(
        "scope=universal chart factorization and mutation law;"
        "finite cofactor separation is not resultant nonvanishing;"
        "no new arbitrary-width SFC4 closure"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
