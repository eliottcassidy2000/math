#!/usr/bin/env python3
"""Exact width-20 modular gate for THM-2955's fixed rank-35 cofactor.

For support (0,1,2,20), the fixed THM-2949 cofactor has two positive-real
seams and no one-sign Gregory--Newton certificate within the old 4M search.
This companion proves nonvanishing at every nonnegative integral depth by
exact division and a complete root-free residue table modulo 97.
"""

from __future__ import annotations

import importlib.util
import sys
from hashlib import sha256
from math import gcd
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "gmc_fixed_rank_thirty_five_cofactor_newton_atlas_thm2949.py"
EXPECTED_SOURCE_SHA256 = (
    "9a1c7068e079e232dc97fd6eb925621aa74b3d636380a85995f8e0db8b30aa54"
)
OFFSETS = (0, 1, 2, 20)
DEGREE = 55 * 20 - 35
PRIME = 97


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_digest(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


require(
    lf_digest(SOURCE) == EXPECTED_SOURCE_SHA256,
    "THM-2949 source changed",
)
SPEC = importlib.util.spec_from_file_location("thm2949_for_2955", SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load THM-2949")
thm2949 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = thm2949
SPEC.loader.exec_module(thm2949)


def primitive_coefficients(polynomial: object) -> tuple[int, ...]:
    coefficients = tuple(int(value) for value in polynomial.coeffs())
    content = 0
    for coefficient in coefficients:
        content = gcd(content, abs(coefficient))
    require(content != 0, "zero crossing factor")
    primitive = tuple(coefficient // content for coefficient in coefficients)
    if primitive[-1] < 0:
        primitive = tuple(-coefficient for coefficient in primitive)
    require(
        gcd(*map(abs, primitive)) == 1,
        "primitive content normalization failed",
    )
    return primitive


def evaluate_mod(coefficients: tuple[int, ...], residue: int) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * residue + coefficient) % PRIME
    return answer


def negative_linear_factor() -> object:
    """Construct the explicit degree-425 negative-root factor."""
    variable = thm2949.thm2943.X
    answer = (2 * variable + 1) ** 5
    answer *= (variable + 1) ** 6
    answer *= (variable + 2) ** 21
    for root in range(3, 7):
        answer *= (variable + root) ** 26
    for root in range(7, 11):
        answer *= (variable + root) ** 24
    for root in range(11, 14):
        answer *= (variable + root) ** 20
    for root in range(14, 21):
        answer *= (variable + root) ** 19
    require(answer.degree() == 425, "negative-root factor degree changed")
    return answer


def main() -> None:
    rows = thm2949.forms_and_rows(20, 1, 2)
    values = [
        thm2949.fixed_cofactor_value(rows, depth)
        for depth in range(DEGREE + 1)
    ]
    cofactor = thm2949.thm2943.interpolate(values)
    require(cofactor.degree() == DEGREE, "cofactor degree changed")
    for depth in (DEGREE + 1, DEGREE + 7, DEGREE + 80):
        require(
            int(cofactor(depth)) == thm2949.fixed_cofactor_value(rows, depth),
            f"outside-grid determinant mismatch at depth {depth}",
        )

    quadratic, cubic, _quartic = thm2949.thm2943.polynomial_forms(
        20,
        (0, 1, 2),
    )
    q200 = quadratic[(2, 0, 0)]
    c300 = cubic[(3, 0, 0)]
    require(
        q200.degree() == 19
        and all(value > 0 for value in q200.coeffs()),
        "q200 positivity/degree changed",
    )
    require(
        c300.degree() == 39
        and all(value < 0 for value in c300.coeffs()),
        "c300 orientation/degree changed",
    )

    negative = negative_linear_factor()
    forced = c300 * q200**5 * negative
    require(cofactor % forced == 0, "explicit degree-425 division failed")
    crossing = cofactor // forced
    require(crossing.degree() == 506, "crossing quotient degree changed")
    primitive = primitive_coefficients(crossing)
    require(len(primitive) == 507, "primitive quotient degree changed")

    residues = tuple(evaluate_mod(primitive, residue) for residue in range(PRIME))
    require(all(residues), "degree-506 quotient acquired a root mod 97")
    residue_digest = sha256(",".join(map(str, residues)).encode()).hexdigest()
    primitive_digest = sha256(str(primitive).encode()).hexdigest()

    integer_signs = tuple(
        "+" if cofactor(depth) > 0 else "-" if cofactor(depth) < 0 else "0"
        for depth in range(201)
    )
    require("0" not in integer_signs, "low integral cofactor zero")
    require(
        integer_signs[:62] == ("+",) * 62
        and integer_signs[62:117] == ("-",) * 55
        and integer_signs[117:] == ("+",) * 84,
        "low sign-run profile changed",
    )

    print("WIDTH-20 FIXED RANK-35 MOD-97 INTEGER-ROOT GATE")
    print(f"support={OFFSETS};cofactor_degree={cofactor.degree()}")
    print(
        f"q200_degree={q200.degree()};q200_coefficients_positive=YES;"
        f"minus_c300_degree={c300.degree()};"
        "minus_c300_coefficients_positive=YES"
    )
    print(
        "quotient_degree=931;crossing_degree=506;"
        "crossing_defined_by_exact_B425_division=YES;"
        "negative_linear_degree=425"
    )
    print(
        "negative_linear_profile="
        "-1/2^5,-1^6,-2^21,-3..-6^26,-7..-10^24,"
        "-11..-13^20,-14..-20^19"
    )
    print(f"crossing_primitive_sha256={primitive_digest}")
    print(f"modulus={PRIME};residue_nonzeros={sum(bool(x) for x in residues)}/97")
    print("mod97_values=" + ",".join(map(str, residues)))
    print(f"mod97_values_sha256={residue_digest}")
    print("outside_grid_direct_determinants=3")
    print("integer_sign_runs_0_to_200=0..61:+,62..116:-,117..200:+")
    print("all_nonnegative_integer_depths_nonzero=YES")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
