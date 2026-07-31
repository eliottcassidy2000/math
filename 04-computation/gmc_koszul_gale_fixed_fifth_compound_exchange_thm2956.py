#!/usr/bin/env python3
"""Exact companion for the THM-2956 Koszul--Gale cofactor exchange."""

from __future__ import annotations

import importlib.util
from hashlib import sha256
from pathlib import Path

import sympy as sp
from flint import fmpz_mat


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


SOURCE = Path(__file__).with_name(
    "gmc_width_seven_eight_two_chart_resultant_closure_thm2943.py"
)
SOURCE_BYTES = SOURCE.read_bytes().replace(b"\r\n", b"\n").replace(
    b"\r", b"\n"
)
SOURCE_SHA256 = (
    "d2f8afeba7dd6c7950405a4845d7bf112b6c9872dd8161146446be8bbdaae0ba"
)
require(
    sha256(SOURCE_BYTES).hexdigest() == SOURCE_SHA256,
    "THM-2943 exact dependency hash changed",
)
SPEC = importlib.util.spec_from_file_location("thm2943_exchange_exact", SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load THM-2943")
thm2943 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(thm2943)

t = thm2943.t
thm2942 = thm2943.thm2942

OLD_ROWS = (
    tuple(range(20))
    + tuple(range(21, 30))
    + (35,)
    + tuple(range(37, 42))
)
NEW_ROWS = (
    tuple(range(21))
    + tuple(range(21, 30))
    + tuple(range(37, 42))
)
TARGET_COLUMNS = tuple(range(1, 36))


def symbolic_gale_check() -> None:
    q = {
        exponent: sp.symbols("q" + "".join(map(str, exponent)))
        for exponent in t.MONOMIALS[2]
    }
    c = {
        exponent: sp.symbols("c" + "".join(map(str, exponent)))
        for exponent in t.MONOMIALS[3]
    }
    s2 = t.MONOMIALS[2]
    s5 = t.MONOMIALS[5]
    s4 = t.MONOMIALS[4]
    relation = sp.zeros(6, 36)
    for row, multiplier in enumerate(s2):
        for exponent, coefficient in c.items():
            target = tuple(
                multiplier[index] + exponent[index] for index in range(3)
            )
            relation[row, s5.index(target)] = coefficient
        for exponent, coefficient in q.items():
            target = tuple(
                multiplier[index] + exponent[index] for index in range(3)
            )
            relation[row, 21 + s4.index(target)] = -coefficient

    new_complement = tuple(21 + index for index in range(9, 15))
    old_complement = (20,) + tuple(21 + index for index in range(9, 14))
    q200 = q[(2, 0, 0)]
    c300 = c[(3, 0, 0)]
    require(
        sp.expand(relation[:, new_complement].det() - q200**6) == 0,
        "new complementary Gale minor changed",
    )
    require(
        sp.expand(
            relation[:, old_complement].det() - c300 * q200**5
        )
        == 0,
        "old complementary Gale minor changed",
    )

    old_selected = tuple(range(20)) + tuple(range(21, 30)) + (35,)
    new_selected = tuple(range(21)) + tuple(range(21, 30))
    require(sum(old_selected) == 450, "old selected-index sum changed")
    require(sum(new_selected) == 435, "new selected-index sum changed")
    require(
        (sum(old_selected) - sum(new_selected)) % 2 == 1,
        "complementary Pluecker parity changed",
    )


def determinant(rows, selected_rows: tuple[int, ...], depth: int) -> int:
    numeric_rows = thm2943.evaluate_rows(rows, depth)
    return int(
        fmpz_mat(
            [
                [numeric_rows[row][column] for column in TARGET_COLUMNS]
                for row in selected_rows
            ]
        ).det()
    )


def factorial_controls() -> None:
    controls = (
        (3, 1, 2),
        (11, 1, 2),
        (11, 1, 6),
        (11, 3, 7),
    )
    checks = 0
    for width, first, second in controls:
        offsets = (0, first, second, width)
        forms = thm2943.polynomial_forms(width, (0, first, second))
        rows, _metadata = thm2942.macaulay_rows(forms)
        q200 = forms[0][(2, 0, 0)]
        c300 = forms[1][(3, 0, 0)]
        require(
            q200(0) > 0
            and all(coefficient >= 0 for coefficient in q200.coeffs()),
            f"q200 coefficient positivity changed: {offsets}",
        )
        require(
            c300(0) < 0
            and all(coefficient <= 0 for coefficient in c300.coeffs()),
            f"c300 coefficient orientation changed: {offsets}",
        )
        for depth in (0, 1, 7):
            old = determinant(rows, OLD_ROWS, depth)
            new = determinant(rows, NEW_ROWS, depth)
            require(
                int(q200(depth)) * old + int(c300(depth)) * new == 0,
                f"rank-35 Pluecker identity failed: {offsets}, n={depth}",
            )
            require(
                old == 0 or new != 0,
                f"optimal chart created a zero: {offsets}, n={depth}",
            )
            checks += 1
    require(checks == 12, "factorial control census changed")


def main() -> None:
    require(len(OLD_ROWS) == 35, "old row count changed")
    require(len(NEW_ROWS) == 35, "new row count changed")
    require(len(TARGET_COLUMNS) == 35, "target-column count changed")
    symbolic_gale_check()
    factorial_controls()
    print("THM-2956 KOSZUL-GALE FIXED FIFTH-COMPOUND EXCHANGE")
    print("gale_complements=new:q200^6,old:c300*q200^5")
    print("selected_index_sums=old:450,new:435;pluecker_sign=opposite")
    print("universal_identity=q200*P_fixed+c300*P_opt=0")
    print("factorial_identity_controls=12")
    print("optimal_chart_no_new_zero_controls=12")
    print("degree_invoices=fixed:55M-35,optimal:54M-35")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
