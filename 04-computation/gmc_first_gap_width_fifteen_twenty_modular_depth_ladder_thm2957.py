#!/usr/bin/env python3
"""Exact THM-2957 modular ladder for (0,1,2,M), 15<=M<=20.

This research companion imports the canonical THM-2949 constructor and no
scratch mathematics.  For each width it exact-divides the fixed cofactor by

    c300*q200^5*B_M,

where B_M is the explicit negative-root floor-law factor, and proves the
remaining primitive core has no integral root using a frozen root-free prime.
It also checks three determinants outside interpolation, the complete low
integer sign runs, and a one-sign Newton vector at the final positive-tail
base.  Width 17 is the first member whose tail base exceeds 4M.
"""

from __future__ import annotations

import importlib.util
import sys
from hashlib import sha256
from math import gcd
from pathlib import Path


ROOT = Path(__file__).parents[1]
SOURCE = ROOT / (
    "04-computation/"
    "gmc_fixed_rank_thirty_five_cofactor_newton_atlas_thm2949.py"
)
SOURCE_BYTES = SOURCE.read_bytes().replace(b"\r\n", b"\n").replace(
    b"\r", b"\n"
)
SOURCE_SHA256 = (
    "9a1c7068e079e232dc97fd6eb925621aa74b3d636380a85995f8e0db8b30aa54"
)
if sha256(SOURCE_BYTES).hexdigest() != SOURCE_SHA256:
    raise RuntimeError("THM-2949 exact dependency hash changed")
SPEC = importlib.util.spec_from_file_location("thm2949_modular_ladder", SOURCE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load THM-2949")
thm2949 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = thm2949
SPEC.loader.exec_module(thm2949)

thm2943 = thm2949.thm2943
thm2942 = thm2949.thm2942
require = thm2949.require
x = thm2949.t.POLY_X


WIDTHS = tuple(range(15, 21))
GATE_PRIMES = {
    15: 83,
    16: 71,
    17: 107,
    18: 73,
    19: 79,
    20: 97,
}
CORE_DEGREES = {
    15: 378,
    16: 403,
    17: 431,
    18: 453,
    19: 482,
    20: 506,
}
TAIL_BASES = {
    15: 52,
    16: 63,
    17: 75,
    18: 88,
    19: 102,
    20: 117,
}
SIGN_RUNS = {
    15: "0..28:+,29..51:-,52..225:+",
    16: "0..33:+,34..62:-,63..256:+",
    17: "0..39:+,40..74:-,75..289:+",
    18: "0..46:+,47..87:-,88..324:+",
    19: "0..53:+,54..101:-,102..361:+",
    20: "0..61:+,62..116:-,117..400:+",
}
CORE_DIGESTS = {
    15: "b6fea710768bfa6be0c7c4c529c1698f1a0ddff0dbf2e95e568272a020823423",
    16: "505e5c935a11ec2beff10d22d16768a0bfe1878b37c50ac19f8ada852b32950b",
    17: "447b5b527c3b230aebb088a9733ff0480e1c8822f3ffe32deecbcaec463669bc",
    18: "75a1081c87b3322fa3f81cd6bd35520cf963223e141a1dc84ca020395f801148",
    19: "86c1a256e878afa7b9890818430b77d87134bcdb68654699cb0bb385fc5a6acf",
    20: "a0176d5c76c67c931de4217da581d1599e29d573514917f416638c1d2497355a",
}
GATE_DIGESTS = {
    15: "5c3f2d2af0bf23808b9a7b0fc023a38eb3d18d57c49e86e58f6c5bdd074dc26e",
    16: "a1f03010ddecd4296d9ef80ccc6cb9503c3d4a7dfc33fb153dc23cd5ce9a5208",
    17: "bce5b1c385b75d3dc865c89672f64a120ca926d3a5bb4885f2f49d37b3f75242",
    18: "7ed43f14e61981c1c0cec8fd100ffe811639fc09c6e2247d65840d61429007b2",
    19: "0738bb31e61f4d2b2d77bc91a2c14599beb3f64ecfc1dfdb28cdc03bf515f59c",
    20: "81743f0f0a65a8c3d3b92e313b8ee73b3376cd656b9c9943d93be34f07181978",
}


def primitive_coefficients(poly) -> tuple[int, ...]:
    coefficients = tuple(int(value) for value in poly.coeffs())
    content = 0
    for coefficient in coefficients:
        content = gcd(content, abs(coefficient))
    require(content != 0, "zero modular core")
    coefficients = tuple(value // content for value in coefficients)
    if coefficients[-1] < 0:
        coefficients = tuple(-value for value in coefficients)
    require(gcd(*map(abs, coefficients)) == 1, "primitive content failure")
    return coefficients


def evaluate_mod(
    coefficients: tuple[int, ...],
    residue: int,
    prime: int,
) -> int:
    answer = 0
    for coefficient in reversed(coefficients):
        answer = (answer * residue + coefficient) % prime
    return answer


def negative_factor(width: int):
    """The stable first-gap B_M floor-law factor for 15<=M<=20."""
    answer = (2 * x + 1) ** 5 * (x + 1) ** 6 * (x + 2) ** 21
    for root in range(3, width + 1):
        exponent = (
            26
            if 3 * root <= width
            else 24
            if 2 * root <= width
            else 20
            if 3 * root <= 2 * width
            else 19
        )
        answer *= (x + root) ** exponent
    return answer


def sign_runs(polynomial, limit: int) -> str:
    signs = tuple(
        "+" if polynomial(depth) > 0
        else "-"
        if polynomial(depth) < 0
        else "0"
        for depth in range(limit + 1)
    )
    require("0" not in signs, "integral zero in sign window")
    records = []
    start = 0
    for index in range(1, len(signs) + 1):
        if index == len(signs) or signs[index] != signs[start]:
            records.append(f"{start}..{index - 1}:{signs[start]}")
            start = index
    return ",".join(records)


def audit_width(width: int) -> str:
    offsets = (0, 1, 2, width)
    forms = thm2943.polynomial_forms(width, (0, 1, 2))
    rows = thm2949.forms_and_rows(width, 1, 2)
    degree = 55 * width - 35
    values = [
        thm2949.fixed_cofactor_value(rows, depth)
        for depth in range(degree + 1)
    ]
    polynomial = thm2943.interpolate(values)
    require(polynomial.degree() == degree, f"degree changed: M={width}")

    outside = (degree + 4 * width - 2, degree + 4 * width - 1, degree + 4 * width)
    for depth in outside:
        require(
            int(polynomial(depth)) == thm2949.fixed_cofactor_value(rows, depth),
            f"outside determinant mismatch: M={width},n={depth}",
        )

    q200 = forms[0][(2, 0, 0)]
    c300 = forms[1][(3, 0, 0)]
    require(
        q200.degree() == width - 1
        and q200(0) > 0
        and all(value > 0 for value in q200.coeffs()),
        f"q200 positivity changed: M={width}",
    )
    require(
        c300.degree() == 2 * width - 1
        and c300(0) < 0
        and all(value < 0 for value in c300.coeffs()),
        f"c300 orientation changed: M={width}",
    )
    b_factor = negative_factor(width)
    forced = c300 * q200**5 * b_factor
    require(polynomial % forced == 0, f"explicit B division failed: M={width}")
    core = polynomial // forced
    require(
        core.degree() == CORE_DEGREES[width],
        f"core degree changed: M={width}",
    )
    core_coefficients = primitive_coefficients(core)
    core_digest = sha256(str(core_coefficients).encode()).hexdigest()
    require(core_digest == CORE_DIGESTS[width], f"core digest changed: M={width}")

    prime = GATE_PRIMES[width]
    gate_values = tuple(
        evaluate_mod(core_coefficients, residue, prime)
        for residue in range(prime)
    )
    require(all(gate_values), f"modular root appeared: M={width}")
    gate_digest = sha256(",".join(map(str, gate_values)).encode()).hexdigest()
    require(gate_digest == GATE_DIGESTS[width], f"gate digest changed: M={width}")

    runs = sign_runs(polynomial, width * width)
    require(runs == SIGN_RUNS[width], f"sign runs changed: M={width}")
    base = TAIL_BASES[width]
    newton = thm2949.newton_coefficients(
        [int(polynomial(base + offset)) for offset in range(degree + 1)]
    )
    require(
        newton[0] > 0 and all(value >= 0 for value in newton),
        f"tail Newton cone changed: M={width}",
    )
    newton_digest = sha256(",".join(map(str, newton)).encode()).hexdigest()

    record = (
        f"M={width};degree={degree};B_degree={b_factor.degree()};"
        f"core_degree={core.degree()};root_free_prime={prime};"
        f"core_digest={core_digest};gate_digest={gate_digest};"
        f"sign_runs={runs};GN_tail_base={base};"
        f"GN_tail_digest={newton_digest};"
        f"fourM_status={'PASS' if base <= 4 * width else 'FAIL'};"
        "outside_determinants=3;integer_nonzero=YES"
    )
    print(record, flush=True)
    return record


def main() -> None:
    print("THM-2957 FIRST-GAP WIDTH 15-20 MODULAR DEPTH LADDER")
    print(f"thm2949_dependency_sha256={SOURCE_SHA256}")
    print("supports=(0,1,2,M),15<=M<=20;fixed_rank35_cofactor=YES")
    records = [audit_width(width) for width in WIDTHS]
    require(TAIL_BASES[16] <= 4 * 16, "last 4M pass changed")
    require(TAIL_BASES[17] > 4 * 17, "first 4M failure changed")
    print(
        "record_digest="
        + sha256("\n".join(records).encode()).hexdigest()
    )
    print("first_4M_failure=M17;sharp_hostile=M20")
    print("no_monotone_prime_law_claimed=YES")
    print("consequence=first-window_SFC4_on_six_first-gap_supports")
    print(
        "scope=integer_depths;not_complete_width15-20;"
        "not_arbitrary_width;not_real-ray positivity;not NC2 or GMC2"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
