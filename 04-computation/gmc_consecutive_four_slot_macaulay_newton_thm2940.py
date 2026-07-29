#!/usr/bin/env python3
"""Exact companion for THM-2940.

This script gives a short fixed-chart proof of first-window SFC(4) on
the translated consecutive support (n,n+1,n+2,n+3).  It hash-pins the
audited ordinary-monomial constructor of THM-2921, replaces its
diameter-four denominators by the exact width-three denominators, and
checks one degree-seven Macaulay minor on the whole depth ray.
"""

from __future__ import annotations

import importlib.util
from fractions import Fraction
from hashlib import sha256
from pathlib import Path


SOURCE_NAME = "gmc_diameter_four_nonconsecutive_macaulay_newton_thm2921.py"
SOURCE_SHA256 = (
    "42e9b5ceddd677d1f2601a9d5d668c9437281596b65999ddcb8549d4e0b9bf64"
)
OFFSETS = (0, 1, 2, 3)
ORDERS = (2, 3, 4)
AUDIT_DEPTHS = (0, 1, 2, 7, 31, 97, 197)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path: Path) -> bytes:
    return path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")


def load_constructor():
    path = Path(__file__).with_name(SOURCE_NAME)
    require(path.is_file(), f"missing constructor dependency: {path}")
    require(
        sha256(lf_bytes(path)).hexdigest() == SOURCE_SHA256,
        "THM-2921 constructor dependency hash changed",
    )
    spec = importlib.util.spec_from_file_location("thm2921_constructor", path)
    require(spec is not None and spec.loader is not None, "import spec failed")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


T = load_constructor()
fmpz_mat = T.fmpz_mat
fmpz_poly = T.fmpz_poly
POLY_X = T.POLY_X


def denominator(depth: int, order: int) -> int:
    answer = 1
    for root in (1, 2):
        answer *= (depth + root) ** (order - 1)
    answer *= (depth + 3) ** (order - 2)
    return answer


def denominator_poly(order: int):
    answer = fmpz_poly(1)
    for root in (1, 2):
        answer *= (POLY_X + root) ** (order - 1)
    answer *= (POLY_X + 3) ** (order - 2)
    return answer


def symbolic_scaled_form(order: int):
    """Construct the exact width-three denominator times one moment form."""
    answer = {}
    common_denominator = denominator_poly(order)
    for monomial in T.MONOMIALS[order]:
        directions = tuple(
            direction
            for direction, count in enumerate(monomial)
            for _ in range(count)
        )
        total = (fmpz_poly(0), fmpz_poly(1))
        for mask in range(1 << order):
            selected_offsets = tuple(
                OFFSETS[3]
                if mask & (1 << position)
                else OFFSETS[direction]
                for position, direction in enumerate(directions)
            )
            numerator, tensor_denominator = T.normalized_tensor_poly(
                selected_offsets
            )
            sign = -1 if mask.bit_count() % 2 else 1
            total = T.add_fraction(
                total,
                (sign * numerator, tensor_denominator),
            )
        numerator, coefficient_denominator = T.reduced_fraction(
            T.multinomial(monomial) * total[0],
            total[1],
        )
        require(
            common_denominator % coefficient_denominator == 0,
            f"width-three denominator escape: order={order}, {monomial}",
        )
        answer[monomial] = (
            common_denominator // coefficient_denominator
        ) * numerator
    return answer


def macaulay_rows(depth: int, *, scaled: bool) -> list[list[int]]:
    rows: list[list[int]] = []
    for order in ORDERS:
        form = T.moment_form(depth, order, OFFSETS)
        scale = denominator(depth, order) if scaled else 1
        coefficients = {
            monomial: coefficient * scale
            for monomial, coefficient in form.items()
        }
        require(
            all(
                coefficient.denominator == 1
                for coefficient in coefficients.values()
            ),
            f"nonintegral coefficient: n={depth}, order={order}",
        )
        for multiplier in T.MONOMIALS[T.TARGET_DEGREE - order]:
            row = [0] * len(T.TARGET_MONOMIALS)
            for monomial, coefficient in coefficients.items():
                target = tuple(
                    multiplier[index] + monomial[index]
                    for index in range(3)
                )
                row[T.TARGET_INDEX[target]] = coefficient.numerator
            rows.append(row)
    return rows


def selected_determinant(depth: int, *, scaled: bool = True) -> int:
    rows = macaulay_rows(depth, scaled=scaled)
    return int(
        fmpz_mat([rows[index] for index in T.SELECTED_ROWS]).det()
    )


def direct_modular_audit() -> int:
    checks = 0
    for depth in AUDIT_DEPTHS:
        direct_forms = tuple(
            T.direct_four_variable_form_mod(depth, order, OFFSETS)
            for order in ORDERS
        )
        for order, direct_form in zip(ORDERS, direct_forms):
            normalized_form = T.moment_form(depth, order, OFFSETS)
            reduced = {
                monomial: (
                    coefficient.numerator
                    * pow(
                        coefficient.denominator,
                        T.PRIME - 2,
                        T.PRIME,
                    )
                )
                % T.PRIME
                for monomial, coefficient in normalized_form.items()
            }
            require(
                reduced == direct_form,
                f"direct form mismatch: n={depth}, order={order}",
            )
        direct_rows = T.macaulay_rows_mod(direct_forms)
        direct_minor = T.determinant_mod(
            [direct_rows[index] for index in T.SELECTED_ROWS]
        )
        scale = (
            pow(denominator(depth, 2), 20, T.PRIME)
            * pow(denominator(depth, 3), 10, T.PRIME)
            * pow(denominator(depth, 4), 6, T.PRIME)
        ) % T.PRIME
        require(
            selected_determinant(depth) % T.PRIME
            == direct_minor * scale % T.PRIME,
            f"direct minor mismatch: n={depth}",
        )
        checks += 1
    return checks


def finite_differences(values: list[int]) -> tuple[int, ...]:
    first_values = []
    row = values
    while row:
        first_values.append(row[0])
        row = [
            row[index + 1] - row[index]
            for index in range(len(row) - 1)
        ]
    return tuple(first_values)


def sequence_digest(values: tuple[int, ...]) -> str:
    return sha256(",".join(map(str, values)).encode()).hexdigest()


def form_digest(forms) -> str:
    payload = []
    for order, form in zip(ORDERS, forms):
        for monomial in T.MONOMIALS[order]:
            payload.append(
                f"{order}:{monomial}:"
                + ",".join(map(str, form[monomial].coeffs()))
            )
    return sha256("|".join(payload).encode()).hexdigest()


def main() -> None:
    require(T.is_prime(T.PRIME), "audit modulus is not prime")
    require(len(T.SELECTED_ROWS) == 36, "selected row count changed")

    degree_bound = 20 * 2 + 10 * 5 + 6 * 8
    require(degree_bound == 138, "determinant degree invoice changed")

    mixed_tensor = T.difference_tensor(0, (0, 1), OFFSETS)
    mixed_coefficient = T.moment_form(0, 2, OFFSETS)[(1, 1, 0)]
    require(
        mixed_tensor != 0 and mixed_coefficient == 2 * mixed_tensor,
        "ordinary-monomial multinomial hostile failed",
    )

    symbolic_forms = tuple(
        symbolic_scaled_form(order)
        for order in ORDERS
    )
    symbolic_degrees = tuple(
        max(polynomial.degree() for polynomial in form.values())
        for form in symbolic_forms
    )
    require(
        symbolic_degrees == (2, 5, 8),
        "scaled row degrees changed",
    )

    specialization_checks = 0
    for depth in range(9):
        for order, symbolic_form in zip(ORDERS, symbolic_forms):
            numeric_form = T.moment_form(depth, order, OFFSETS)
            scale = denominator(depth, order)
            for monomial in T.MONOMIALS[order]:
                expected = numeric_form[monomial] * scale
                require(
                    expected.denominator == 1,
                    "symbolic cross-check lost integrality",
                )
                require(
                    int(symbolic_form[monomial](depth))
                    == expected.numerator,
                    (
                        "symbolic/numeric mismatch: "
                        f"n={depth}, order={order}, {monomial}"
                    ),
                )
                specialization_checks += 1

    values = tuple(
        selected_determinant(depth)
        for depth in range(degree_bound + 1)
    )
    newton = finite_differences(list(values))
    require(len(newton) == degree_bound + 1, "Newton count changed")
    require(all(value > 0 for value in newton), "Newton positivity failed")

    raw_zero = selected_determinant(0, scaled=False)
    scale_zero = (
        denominator(0, 2) ** 20
        * denominator(0, 3) ** 10
        * denominator(0, 4) ** 6
    )
    require(
        values[0] == raw_zero * scale_zero,
        "depth-zero original/scaled determinant mismatch",
    )

    modular_checks = direct_modular_audit()
    require(modular_checks == 7, "modular audit count changed")

    print("THM-2940 CONSECUTIVE FOUR-SLOT MACAULAY--NEWTON SHORT CLOSURE")
    print(f"constructor_dependency_sha256={SOURCE_SHA256}")
    print("offsets=(0,1,2,3);ordinary_multinomial_constructor=YES")
    print(
        "selected_rows="
        + ",".join(map(str, T.SELECTED_ROWS))
        + ";allocation=(20Q,10C,6F)"
    )
    print("denominators=D2:(1,1,0);D3:(2,2,1);D4:(3,3,2)")
    print("scaled_row_degrees=(2,5,8);determinant_degree_bound=138")
    print("multinomial_hostile=mixed_quadratic_ratio_2:PASS")
    print(
        f"form_digest={form_digest(symbolic_forms)};"
        f"specialization_checks={specialization_checks}"
    )
    print(
        "newton_base=0;"
        f"newton_count={len(newton)};"
        f"newton_digest={sequence_digest(newton)}"
    )
    print(
        f"n0_sign={'+' if values[0] > 0 else '-'};"
        f"n0_raw_digits={len(str(abs(raw_zero)))};"
        f"first_digits={len(str(newton[0]))};"
        f"last_digits={len(str(newton[-1]))};"
        "all_positive=YES"
    )
    print(
        f"independent_direct_constructor_mod_p={T.PRIME};"
        f"minor_checks={modular_checks};PASS"
    )
    print("symbolic_denominator_division=PASS")
    print("consequence=fixed maximal minor nonzero for every integer n>=0")
    print("scope=translated consecutive four-slot first-window SFC(4)")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
