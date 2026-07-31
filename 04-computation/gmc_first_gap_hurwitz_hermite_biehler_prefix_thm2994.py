#!/usr/bin/env python3
"""Exact Hurwitz/Hermite--Biehler prefix for first-gap norm cores.

Canonical THM-2994 companion.  For the genuine THM-2969 wall-stripped cores
C_M, 6 <= M <= 16, it certifies Hurwitz stability by a primitive
fraction-free Routh array for the reciprocal core.  It also supplies rigorous
Arb corroboration through M=12 and directly isolates the
Hermite--Biehler legs through M=8.

The script deliberately does *not* infer the Newton-circuit no-return property
from stability.  Its final exact control is THM-2991's three-cluster global
return witness P_(2,20).
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

from flint import fmpz, fmpz_poly


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "04-computation"
    / "gmc_first_gap_wall_stripped_norm_core_atlas_thm2969.py"
)
SOURCE_BYTES = SOURCE.read_bytes().replace(bytes((13, 10)), bytes((10,))).replace(
    bytes((13,)), bytes((10,))
)
SOURCE_SHA256 = sha256(SOURCE_BYTES).hexdigest()
EXPECTED_SOURCE_SHA256 = (
    "5be5df0fd058f436593339f78cd6240a47975082d22929c135ba811458753bf5"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


require(
    SOURCE_SHA256 == EXPECTED_SOURCE_SHA256,
    "THM-2969 exact dependency hash changed",
)


EXPECTED_ROUTH = {
    6: (121, 38376, "bf5e565214f4cf3257d1a6a1055b8bf6b58b9ae68bd55c1aa32236dbd814f3bb"),
    7: (144, 55189, "7e1ed405d784ce7c07664867ae9c3fc32547bd512057357881a18c41aa2e5f61"),
    8: (164, 74026, "60844d75371190769d5db5435ebe5da049c693afd29ea48b1ddf261d11b436ab"),
    9: (184, 97503, "c9bea1dc9b8ac06c50d46d928e34bc1caaeb26fd2fd24debdcb94ffcc2863cb5"),
    10: (205, 119890, "7151ba9b535619b9006d6f72abeeef6a8ff979c75a5476b3b04024d968dfcd9f"),
    11: (226, 148547, "ce7eddf2ecee0dda531b1bfc7cc8123a1402c62c08cffddaaecdf49238103fa3"),
    12: (244, 174313, "0082b706509af1d93a5afe5f0882d0395e19def94b80105a50a5da3cc7f53ae8"),
    13: (268, 209343, "eb765ac99562c13db36d8a47dfd4142f04bb196164b153b2d1faa4ad4d28bf68"),
    14: (288, 242649, "474e601046f515353c6a6aa2f10cdb67328c9c9c39cf61bbec9d53249c1d7339"),
    15: (308, 283228, "e70693e4413446420592c976921aeddd7d6c73529fa0380c427b3a95de226b68"),
    16: (329, 321582, "5d3065df61595d91a34815b597d702fba8ca2c1e8e13ee2448369fb3634d48fc"),
}

EXPECTED_ROOT_COUNTS = {
    6: (37, 84),
    7: (44, 100),
    8: (46, 118),
    9: (54, 130),
    10: (59, 146),
    11: (66, 160),
    12: (64, 180),
}

EXPECTED_HB = {
    6: (60, 60, "E<O...<O"),
    7: (72, 71, "E<O...<E"),
    8: (82, 81, "E<O...<E"),
}

EXPECTED_RECORD_DIGEST = (
    "614dd91698c171c177f029c8fe733561cbe47001cc2d3a61312dbf81ddc3f35f"
)


def load_base(width: int):
    name = f"thm2969_hurwitz_hb_worker_{width}"
    spec = importlib.util.spec_from_file_location(name, SOURCE)
    require(spec is not None and spec.loader is not None, "cannot load THM-2969")
    base = importlib.util.module_from_spec(spec)
    sys.modules[name] = base
    spec.loader.exec_module(base)
    return base


def primitive_integer_row(row: list[fmpz]) -> list[fmpz]:
    content = fmpz(0)
    for value in row:
        content = content.gcd(abs(value))
    require(content > 0, "Routh row vanished")
    if content > 1:
        row = [value // content for value in row]
    return row


def fraction_free_routh_first_column(
    coefficients_descending: list[int],
) -> list[fmpz]:
    """Return a positive-row-scaled exact Routh first column.

    If stored rows U,L are positive multiples alpha,beta of the ordinary
    rational rows, then

        L_0 U_(j+1) - U_0 L_(j+1)

    is alpha*beta times the positive ordinary lower pivot times the next
    ordinary row.  Removing positive primitive content is therefore lawful.
    """

    degree = len(coefficients_descending) - 1
    width = (degree + 2) // 2
    upper = [fmpz(value) for value in coefficients_descending[0::2]]
    lower = [fmpz(value) for value in coefficients_descending[1::2]]
    upper += [fmpz(0)] * (width - len(upper))
    lower += [fmpz(0)] * (width - len(lower))
    upper = primitive_integer_row(upper)
    lower = primitive_integer_row(lower)
    first = [upper[0], lower[0]]
    require(upper[0] > 0 and lower[0] > 0, "initial Routh pivot is nonpositive")

    for _index in range(2, degree + 1):
        require(lower[0] > 0, "Routh pivot is nonpositive")
        row = [
            lower[0] * upper[j + 1] - upper[0] * lower[j + 1]
            for j in range(width - 1)
        ] + [fmpz(0)]
        row = primitive_integer_row(row)
        first.append(row[0])
        upper, lower = lower, row

    require(len(first) == degree + 1, "Routh first-column length changed")
    require(all(value > 0 for value in first), "nonpositive Routh pivot")
    return first


def expanded_complex_roots(poly) -> list:
    return [
        root
        for root, multiplicity in poly.complex_roots()
        for _index in range(multiplicity)
    ]


def root_record(core, width: int) -> dict[str, object] | None:
    if width not in EXPECTED_ROOT_COUNTS:
        return None
    roots = expanded_complex_roots(core)
    require(len(roots) == core.degree(), f"root count changed: M={width}")
    real_count = 0
    nonreal_count = 0
    for root in roots:
        require(root.real < 0, f"left-half-plane isolation failed: M={width}")
        if root.imag == 0:
            real_count += 1
        else:
            require(
                not root.imag.contains(0),
                f"root reality undecided: M={width}",
            )
            nonreal_count += 1
    require(
        (real_count, nonreal_count) == EXPECTED_ROOT_COUNTS[width],
        f"real/nonreal root census changed: M={width}",
    )
    return {
        "certified_left": len(roots),
        "nonreal": nonreal_count,
        "real": real_count,
    }


def hb_record(core, width: int) -> dict[str, object] | None:
    if width not in EXPECTED_HB:
        return None
    coefficients = [int(value) for value in core.coeffs()]
    even_leg = fmpz_poly(
        [
            value if index % 2 == 0 else -value
            for index, value in enumerate(coefficients[0::2])
        ]
    )
    odd_leg = fmpz_poly(
        [
            value if index % 2 == 0 else -value
            for index, value in enumerate(coefficients[1::2])
        ]
    )
    roots_by_label: list[tuple[object, str]] = []
    for poly, label in ((even_leg, "E"), (odd_leg, "O")):
        roots = expanded_complex_roots(poly)
        require(len(roots) == poly.degree(), f"HB root count changed: M={width}")
        for root in roots:
            require(root.imag == 0, f"HB leg acquired nonreal root: M={width}")
            require(root.real > 0, f"HB leg acquired nonpositive root: M={width}")
            roots_by_label.append((root.real, label))

    roots_by_label.sort(key=lambda item: float(item[0].mid()))
    for left, right in zip(roots_by_label, roots_by_label[1:]):
        require(left[0] < right[0], f"HB roots overlap: M={width}")
        require(left[1] != right[1], f"HB alternation failed: M={width}")
    pattern = (
        f"{roots_by_label[0][1]}<{roots_by_label[1][1]}"
        f"...<{roots_by_label[-1][1]}"
    )
    expected_even, expected_odd, expected_pattern = EXPECTED_HB[width]
    require(
        even_leg.degree() == expected_even
        and odd_leg.degree() == expected_odd
        and pattern == expected_pattern,
        f"HB degree/pattern changed: M={width}",
    )
    return {
        "even_degree": even_leg.degree(),
        "odd_degree": odd_leg.degree(),
        "pattern": pattern,
    }


def audit_width(width: int) -> dict[str, object]:
    base = load_base(width)
    capture = io.StringIO()
    with contextlib.redirect_stdout(capture):
        core = base.audit_width(width)[0]
    degree, expected_height, expected_digest = EXPECTED_ROUTH[width]
    require(core.degree() == degree, f"core degree changed: M={width}")

    # Native coefficients are constant-first.  Reading them as a descending
    # vector certifies C^vee(s)=s^d C(1/s), not C itself.  C(0)>0, and
    # Re(1/z)=Re(z)/|z|^2, so the two Hurwitz statements are equivalent.
    coefficients = [int(value) for value in core.coeffs()]
    require(coefficients[0] > 0, f"core constant term vanished: M={width}")
    first = fraction_free_routh_first_column(coefficients)
    first_digest = sha256(
        ("\n".join(str(value) for value in first) + "\n").encode()
    ).hexdigest()
    max_height = max(abs(value).bit_length() for value in first)
    require(
        max_height == expected_height and first_digest == expected_digest,
        f"Routh record changed: M={width}",
    )

    return {
        "M": width,
        "core_digest": base.polynomial_digest(core),
        "degree": degree,
        "hb": hb_record(core, width),
        "max_first_column_pivot_height_bits": max_height,
        "orientation": "reciprocal_constant_first_as_descending",
        "positive_pivots": len(first),
        "root_census": root_record(core, width),
        "routh_digest": first_digest,
    }


def separation_controls() -> str:
    # Hurwitz alone need not imply ULC.
    quadratic_ratio = Fraction(1, 4)
    require(quadratic_ratio < 1, "quadratic Hurwitz/ULC hostile changed")

    # THM-2991's smallest three-cluster PF-infinity global return P_(2,20).
    coefficients = (1, 45, 607, 2283, 2920, 1200)
    normalized = tuple(
        Fraction(value, choose)
        for value, choose in zip(coefficients, (1, 5, 10, 10, 5, 1))
    )
    ratios = tuple(
        normalized[index] ** 2
        / (normalized[index - 1] * normalized[index + 1])
        for index in range(1, 5)
    )
    expected = (
        Fraction(810, 607),
        Fraction(368449, 205470),
        Fraction(5212089, 3544880),
        Fraction(42632, 34245),
    )
    require(
        ratios == expected
        and all(ratio > 1 for ratio in ratios)
        and ratios[0] < ratios[1]
        and ratios[3] < ratios[0],
        "THM-2991 three-cluster global-return control changed",
    )
    return (
        "quadratic_hurwitz_ratio=1/4;"
        "pf_infinity_global_return="
        "810/607,368449/205470,5212089/3544880,42632/34245;"
        "R4_lt_R1=1"
    )


def run_checks(workers: int) -> str:
    widths = tuple(range(6, 17))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        records = list(pool.map(audit_width, widths))
    records.sort(key=lambda record: int(record["M"]))
    record_lines = [
        json.dumps(record, sort_keys=True, separators=(",", ":"))
        for record in records
    ]
    record_digest = sha256(("\n".join(record_lines) + "\n").encode()).hexdigest()
    require(record_digest == EXPECTED_RECORD_DIGEST, "global record digest changed")
    lines = [
        "FIRST-GAP HURWITZ/HERMITE-BIEHLER PREFIX",
        f"dependency_THM2969_sha256={SOURCE_SHA256}",
        "scope=M6..16;exact_fraction_free_reciprocal_Routh=1;all_pivots_positive=1",
        *record_lines,
        "direct_HB_Arb_scope=M6..8;full_core_Arb_scope=M6..12",
        separation_controls(),
        f"record_digest={record_digest}",
        "all_exact_checks=PASS",
        "scope_boundary=no_all_width_no_return_or_GMC_consequence",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(1 <= args.workers <= 8, "worker count outside 1..8")
    transcript = run_checks(args.workers)
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
