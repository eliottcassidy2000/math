#!/usr/bin/env python3
"""Exact companion for THM-2924, the diameter-six SFC(4) atlas.

The correctly typed moment constructors are inherited from the hash-pinned
THM-2922 companion.  This script independently defines the width-six
denominators, proves their exact polynomial divisibility, and checks the same
fixed degree-seven Macaulay minor for all ten diameter-six support types.
"""

from __future__ import annotations

import importlib.util
from hashlib import sha256
from pathlib import Path

from flint import fmpz_mat, fmpz_poly


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


SOURCE = Path(__file__).with_name(
    "gmc_diameter_five_macaulay_newton_atlas_thm2922.py"
)
SOURCE_BYTES = SOURCE.read_bytes().replace(b"\r\n", b"\n")
require(
    sha256(SOURCE_BYTES).hexdigest()
    == "9f146073cf2f953bbbfcaaf20f1f58dbb25a105316d97a486e80db2be62b9b2e",
    "THM-2922 constructor dependency hash changed",
)
SPEC = importlib.util.spec_from_file_location("thm2922_exact", SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load THM-2922")
w5 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(w5)
t = w5.t


WIDTH = 6
FAMILIES = tuple(
    (0, first, second, WIDTH)
    for first in range(1, WIDTH)
    for second in range(first + 1, WIDTH)
)
NEWTON_BASES = tuple(
    5 if offsets[1] == 1 else 2 if offsets[1] == 2 else 0
    for offsets in FAMILIES
)
DEGREE_BOUND = 312
AUDIT_DEPTHS = (0, 1, 2, 11, 37, 113, 251)


def denominator_poly(order: int) -> fmpz_poly:
    answer = fmpz_poly(1)
    for shift in range(1, WIDTH):
        answer *= (t.POLY_X + shift) ** (order - 1)
    answer *= (t.POLY_X + WIDTH) ** (order - 2)
    return answer


def denominator(depth: int, order: int) -> int:
    answer = 1
    for shift in range(1, WIDTH):
        answer *= (depth + shift) ** (order - 1)
    answer *= (depth + WIDTH) ** (order - 2)
    return answer


def symbolic_scaled_form(
    order: int,
    offsets: tuple[int, int, int, int],
) -> dict[tuple[int, int, int], fmpz_poly]:
    """Construct D_order times every ordinary width-six coefficient."""
    common_denominator = denominator_poly(order)
    answer: dict[tuple[int, int, int], fmpz_poly] = {}
    for monomial in t.MONOMIALS[order]:
        directions = tuple(
            direction
            for direction, count in enumerate(monomial)
            for _ in range(count)
        )
        total = (fmpz_poly(0), fmpz_poly(1))
        for mask in range(1 << order):
            selected_offsets = tuple(
                offsets[3] if mask & (1 << position) else offsets[direction]
                for position, direction in enumerate(directions)
            )
            numerator, tensor_denominator = t.normalized_tensor_poly(
                selected_offsets
            )
            sign = -1 if mask.bit_count() % 2 else 1
            total = t.add_fraction(
                total,
                (sign * numerator, tensor_denominator),
            )
        numerator, coefficient_denominator = t.reduced_fraction(
            t.multinomial(monomial) * total[0],
            total[1],
        )
        require(
            common_denominator % coefficient_denominator == 0,
            f"denominator escape: {offsets}, m={order}, {monomial}",
        )
        answer[monomial] = (
            common_denominator // coefficient_denominator
        ) * numerator
    return answer


def selected_determinant(
    depth: int,
    offsets: tuple[int, int, int, int],
    *,
    scaled: bool = True,
) -> int:
    rows: list[list[int]] = []
    for order in t.ORDERS:
        scale = denominator(depth, order) if scaled else 1
        form = t.moment_form(depth, order, offsets)
        for multiplier in t.MONOMIALS[t.TARGET_DEGREE - order]:
            row = [0] * len(t.TARGET_MONOMIALS)
            for monomial, coefficient in form.items():
                coefficient *= scale
                require(
                    coefficient.denominator == 1,
                    f"nonintegral row: {offsets}, n={depth}, m={order}",
                )
                target = tuple(
                    multiplier[index] + monomial[index]
                    for index in range(3)
                )
                row[t.TARGET_INDEX[target]] = coefficient.numerator
            rows.append(row)
    return int(
        fmpz_mat([rows[index] for index in t.SELECTED_ROWS]).det()
    )


def direct_minor_mod(
    depth: int,
    offsets: tuple[int, int, int, int],
) -> int:
    forms = tuple(
        t.direct_four_variable_form_mod(depth, order, offsets)
        for order in t.ORDERS
    )
    rows = t.macaulay_rows_mod(forms)
    return t.determinant_mod(
        [rows[index] for index in t.SELECTED_ROWS]
    )


def finite_differences(values: list[int]) -> tuple[int, ...]:
    answer = []
    while values:
        answer.append(values[0])
        values = [
            values[index + 1] - values[index]
            for index in range(len(values) - 1)
        ]
    return tuple(answer)


def integer_sequence_digest(values: tuple[int, ...]) -> str:
    return sha256(",".join(map(str, values)).encode()).hexdigest()


def form_digest(
    forms: tuple[dict[tuple[int, int, int], fmpz_poly], ...],
) -> str:
    payload = []
    for order, form in zip(t.ORDERS, forms):
        for monomial in t.MONOMIALS[order]:
            payload.append(
                f"{order}:{monomial}:"
                + ",".join(map(str, form[monomial].coeffs()))
            )
    return sha256("|".join(payload).encode()).hexdigest()


def main() -> None:
    require(
        20 * 5 + 10 * 11 + 6 * 17 == DEGREE_BOUND,
        "determinant degree invoice changed",
    )
    require(
        len(FAMILIES) == 10
        and NEWTON_BASES.count(5) == 4
        and NEWTON_BASES.count(2) == 3
        and NEWTON_BASES.count(0) == 3,
        "family/base atlas changed",
    )

    print("THM-2924 DIAMETER-SIX MACAULAY--NEWTON ATLAS")
    print(
        "constructor=hash-pinned THM-2922/2921 ordinary multinomial forms;"
        "width-six denominators rebuilt"
    )
    print(
        "selected_rows="
        + ",".join(map(str, t.SELECTED_ROWS))
        + ";allocation=(20Q,10C,6F)"
    )
    print("scaled_row_degrees=(5,11,17);determinant_degree_bound=312")

    modular_checks = 0
    positive_newton_coefficients = 0
    for offsets, newton_base in zip(FAMILIES, NEWTON_BASES):
        forms = tuple(
            symbolic_scaled_form(order, offsets)
            for order in t.ORDERS
        )
        degrees = tuple(
            max(polynomial.degree() for polynomial in form.values())
            for form in forms
        )
        require(
            degrees == (5, 11, 17),
            f"scaled row degrees changed for {offsets}",
        )

        for depth in range(9):
            for order, form in zip(t.ORDERS, forms):
                numeric = t.moment_form(depth, order, offsets)
                scale = denominator(depth, order)
                for monomial, polynomial in form.items():
                    expected = numeric[monomial] * scale
                    require(
                        expected.denominator == 1
                        and int(polynomial(depth)) == expected.numerator,
                        f"symbolic/numeric mismatch: {offsets}, n={depth}",
                    )

        values = tuple(
            selected_determinant(depth, offsets)
            for depth in range(
                newton_base + DEGREE_BOUND + 2
            )
        )
        newton_with_exhaustion = finite_differences(
            list(
                values[
                    newton_base:newton_base + DEGREE_BOUND + 2
                ]
            )
        )
        newton = newton_with_exhaustion[:-1]
        require(
            len(newton) == DEGREE_BOUND + 1
            and all(value > 0 for value in newton)
            and newton_with_exhaustion[-1] == 0,
            f"Newton positivity/degree exhaustion failed for {offsets}",
        )
        positive_newton_coefficients += len(newton)
        require(
            all(values[depth] != 0 for depth in range(newton_base)),
            f"exceptional low depth vanished for {offsets}",
        )

        raw_zero = selected_determinant(0, offsets, scaled=False)
        scale_zero = (
            denominator(0, 2) ** 20
            * denominator(0, 3) ** 10
            * denominator(0, 4) ** 6
        )
        require(
            values[0] == raw_zero * scale_zero,
            f"original/scaled depth-zero mismatch for {offsets}",
        )

        for depth in AUDIT_DEPTHS:
            direct_minor = direct_minor_mod(depth, offsets)
            scale = (
                pow(denominator(depth, 2), 20, t.PRIME)
                * pow(denominator(depth, 3), 10, t.PRIME)
                * pow(denominator(depth, 4), 6, t.PRIME)
            ) % t.PRIME
            require(
                selected_determinant(depth, offsets) % t.PRIME
                == direct_minor * scale % t.PRIME,
                f"independent modular minor mismatch: {offsets}, n={depth}",
            )
            modular_checks += 1

        low_signs = "".join(
            "+" if values[depth] > 0 else "-"
            for depth in range(newton_base)
        ) or "none"
        print(
            f"family={offsets};"
            f"newton_base={newton_base};"
            f"low_signs={low_signs};"
            f"n0_raw_digits={len(str(abs(raw_zero)))};"
            f"form_digest={form_digest(forms)};"
            f"newton_digest={integer_sequence_digest(newton)}"
        )
        print(
            f"family={offsets};"
            f"newton_count={len(newton)};"
            f"first_digits={len(str(newton[0]))};"
            f"last_digits={len(str(newton[-1]))};"
            "all_positive=YES;degree_exhaustion=YES"
        )

    require(modular_checks == 70, "modular audit count changed")
    require(
        positive_newton_coefficients == 3130,
        "Newton coefficient census changed",
    )
    print("symbolic_denominator_division=PASS;families=10")
    print(
        f"positive_newton_coefficients={positive_newton_coefficients};"
        "degree_exhaustions=10;PASS"
    )
    print(
        f"independent_direct_constructor_mod_p={t.PRIME};"
        f"minor_checks={modular_checks};PASS"
    )
    print("consequence=fixed maximal minor nonzero for every integer n>=0")
    print(
        "scope=first-window SFC(4) for all ten translated "
        "four-slot supports of diameter exactly six"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
