#!/usr/bin/env python3
"""Exact THM-2959 modular continuation for (0,1,2,M), 21<=M<=24.

The companion imports the promoted THM-2949 constructor.  For each width it
interpolates the fixed rank-35 cofactor, exact-divides

    c300*q200^5*B_M,

where B_M is the THM-2957 floor-law factor.  At width 21 the quotient has the
additional harmless negative-root factor n+17.  After removing precisely that
surplus, each primitive core has a frozen finite-field root-free witness.

The program also checks three direct determinants outside interpolation, every
integer sign through M^2, and the complete Gregory--Newton vector at the final
positive-tail base.  All truth-bearing arithmetic is exact.
"""

from __future__ import annotations

import importlib.util
import sys
from hashlib import sha256
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / (
    "04-computation/"
    "gmc_fixed_rank_thirty_five_cofactor_newton_atlas_thm2949.py"
)
SOURCE_BYTES = SOURCE.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
SOURCE_SHA256 = sha256(SOURCE_BYTES).hexdigest()
EXPECTED_SOURCE_SHA256 = (
    "9a1c7068e079e232dc97fd6eb925621aa74b3d636380a85995f8e0db8b30aa54"
)
if SOURCE_SHA256 != EXPECTED_SOURCE_SHA256:
    raise RuntimeError("THM-2949 exact dependency hash changed")
SPEC = importlib.util.spec_from_file_location("thm2949_for_thm2959", SOURCE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("cannot load THM-2949")
thm2949 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = thm2949
SPEC.loader.exec_module(thm2949)

thm2943 = thm2949.thm2943
require = thm2949.require
x = thm2949.t.POLY_X
WIDTHS = tuple(range(21, 25))

EXPECTED = {
    21: {
        "core_degree": 531,
        "prime": 127,
        "core_digest": (
            "3826e18fe51b3f8f3f7359257355823fc107cdb414b54447728ec813dac0e1fe"
        ),
        "gate_digest": (
            "365eeeafc09b137555bae4711574baa3b14c7d2b36c73bab9bf027c41535b0e7"
        ),
        "runs": "0..70:+,71..133:-,134..441:+",
        "base": 134,
        "newton_digest": (
            "d8d073a344a7649afee6d4373c65d6bf70d52662585aa6816850ae5ed41b8fca"
        ),
    },
    22: {
        "core_degree": 557,
        "prime": 107,
        "core_digest": (
            "5fa6bb6af27d164d4d69e40f3da1be5a7f7885968068701f25a1e29c8f7d8bd7"
        ),
        "gate_digest": (
            "9432915e09a30033831510b7832068ec47a2635147962af06dd5a5425f36b488"
        ),
        "runs": "0..79:+,80..150:-,151..484:+",
        "base": 151,
        "newton_digest": (
            "13494fc018c276ea8fedb0fc83e9aca50e66d533696f84ee982818d7a848fee7"
        ),
    },
    23: {
        "core_degree": 585,
        "prime": 113,
        "core_digest": (
            "4aa61de9c1793b40354540dbdc06bb771aa49795ddd2e5d694a5d90743bc0c8b"
        ),
        "gate_digest": (
            "deb2e21c89809291ca19d728dc6a8682cbfd04b689d44936c1bd88006cd4c31c"
        ),
        "runs": "0..89:+,90..169:-,170..529:+",
        "base": 170,
        "newton_digest": (
            "36c93dd5e27fb71436ceae3de415fe2b7e01d895bcfe743bc30498af008e8e7c"
        ),
    },
    24: {
        "core_degree": 607,
        "prime": 137,
        "core_digest": (
            "0379f6fdc92f0380698f7e48fc3b8b4a9a6aeba4cbdd52f3765ca23a710992a2"
        ),
        "gate_digest": (
            "315d7cea1008f6546286f580889b1d65efe4d144140bda83c5629482e8982129"
        ),
        "runs": "0..99:+,100..189:-,190..576:+",
        "base": 190,
        "newton_digest": (
            "081beb5f76430d613afabe15c1a98bf17b37163c937bab9ecfb37df4f5320945"
        ),
    },
}

SURPLUS_NEGATIVE_FACTORS = {
    21: ((17, 1),),
}


def primitive_coefficients(poly) -> tuple[int, ...]:
    raw_coefficients = tuple(poly.coeffs())
    require(
        all(value == int(value) for value in raw_coefficients),
        "modular core acquired a nonintegral coefficient",
    )
    coefficients = tuple(int(value) for value in raw_coefficients)
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
    """The THM-2957 floor-law factor, exactly tested on the new bank."""
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


def sign_runs(poly, limit: int) -> tuple[str, int]:
    signs = tuple(
        "+" if poly(depth) > 0
        else "-"
        if poly(depth) < 0
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
    last_negative = max(
        (index for index, sign in enumerate(signs) if sign == "-"),
        default=-1,
    )
    return ",".join(records), last_negative + 1


def newton_at_base(poly, degree: int, base: int) -> str:
    newton = thm2949.newton_coefficients(
        [int(poly(base + offset)) for offset in range(degree + 1)]
    )
    require(
        newton[0] > 0 and all(value >= 0 for value in newton),
        "one-sign Newton cone failed",
    )
    return sha256(",".join(map(str, newton)).encode()).hexdigest()


def audit_width(width: int) -> str:
    expected = EXPECTED[width]
    forms = thm2943.polynomial_forms(width, (0, 1, 2))
    rows = thm2949.forms_and_rows(width, 1, 2)
    degree = 55 * width - 35
    values = [
        thm2949.fixed_cofactor_value(rows, depth)
        for depth in range(degree + 1)
    ]
    polynomial = thm2943.interpolate(values)
    require(polynomial.degree() == degree, f"degree changed: M={width}")

    outside = (
        degree + 4 * width - 2,
        degree + 4 * width - 1,
        degree + 4 * width,
    )
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
    require(polynomial % forced == 0, f"floor-law division failed: M={width}")
    core = polynomial // forced
    surplus_profile = SURPLUS_NEGATIVE_FACTORS.get(width, ())
    for root, exponent in surplus_profile:
        factor = (x + root) ** exponent
        require(core % factor == 0, f"surplus negative factor missing: M={width}")
        require(core(-root) == 0, f"surplus root hostile changed: M={width}")
        core //= factor

    coefficients = primitive_coefficients(core)
    core_digest = sha256(str(coefficients).encode()).hexdigest()
    require(
        core.degree() == expected["core_degree"]
        and core_digest == expected["core_digest"],
        f"primitive core census changed: M={width}",
    )

    gate_prime = expected["prime"]
    gate_values = tuple(
        evaluate_mod(coefficients, residue, gate_prime)
        for residue in range(gate_prime)
    )
    require(all(gate_values), f"root appeared in frozen gate: M={width}")
    gate_digest = sha256(",".join(map(str, gate_values)).encode()).hexdigest()
    require(gate_digest == expected["gate_digest"], f"gate digest changed: M={width}")

    runs, candidate_base = sign_runs(polynomial, width * width)
    require(
        runs == expected["runs"] and candidate_base == expected["base"],
        f"integer sign census changed: M={width}",
    )
    newton_digest = newton_at_base(polynomial, degree, candidate_base)
    require(
        newton_digest == expected["newton_digest"],
        f"Newton digest changed: M={width}",
    )

    record = (
        f"M={width};degree={degree};B_floor_degree={b_factor.degree()};"
        f"surplus_negative_factors={surplus_profile or 'NONE'};"
        f"B_effective_degree="
        f"{b_factor.degree() + sum(e for _, e in surplus_profile)};"
        f"core_degree={core.degree()};root_free_prime={gate_prime};"
        f"core_digest={core_digest};gate_digest={gate_digest};"
        f"sign_runs={runs};GN_tail_base={candidate_base};"
        f"GN_tail_digest={newton_digest};"
        f"fourM_status={'PASS' if candidate_base <= 4 * width else 'FAIL'};"
        "outside_determinants=3;integer_nonzero=YES"
    )
    print(record, flush=True)
    return record


def main() -> None:
    print("THM-2959 FIRST-GAP WIDTH 21-24 MODULAR DEPTH CONTINUATION")
    print(f"thm2949_dependency_sha256={SOURCE_SHA256}")
    print("supports=(0,1,2,M),21<=M<=24;fixed_rank35_cofactor=YES")
    records = [audit_width(width) for width in WIDTHS]
    print("record_digest=" + sha256("\n".join(records).encode()).hexdigest())
    print("width21_surplus=(n+17)^1;negative_root_only=YES")
    print("all_GN_tail_bases_exceed_4M=YES")
    print("consequence=first-window_SFC4_on_four_first-gap_supports")
    print(
        "scope=integer_depths;not_complete_width21-24;"
        "not_arbitrary_width;not_real-ray positivity;"
        "not_nextprime_rank_shortcut;not NC2 or GMC2"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
