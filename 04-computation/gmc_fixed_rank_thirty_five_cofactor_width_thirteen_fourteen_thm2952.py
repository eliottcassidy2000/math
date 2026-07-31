#!/usr/bin/env python3
"""Exact width-thirteen/fourteen continuation of THM-2949.

For every normalized four-slot support

    (0, a, b, M),  1 <= a < b < M,  M in {13, 14},

this companion reconstructs the same fixed 35-by-35 Macaulay cofactor
used by THM-2949.  It proves integer-depth nonvanishing by an exact
finite prefix and a shifted Gregory--Newton tail, with three direct
determinant controls outside each interpolation grid.

Only integer/rational FLINT arithmetic is truth-bearing.  The imported
THM-2949 constructor is LF-hash pinned.
"""

from __future__ import annotations

import importlib.util
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path


WIDTHS = (13, 14)
EXPECTED_WIDTH_DIGESTS = {
    13: "c1e7eaf3a794773573a664f9e83c439817efa6795f4c690b2849b79b7acad923",
    14: "2158531e878df93c158a85d83de9b4f7b10a6d25b3c1143a71d9b2ff3721c3d8",
}
EXPECTED_GLOBAL_DIGEST = (
    "02566dbe576277c3c1582ca778505353496c81d5d5530f9c430f43cc51d0f515"
)
EXPECTED_MAX_BASES = {13: 35, 14: 43}
EXPECTED_BASE_ZERO_COUNTS = {13: 22, 14: 23}
EXPECTED_MIXED_PREFIX_COUNTS = {13: 12, 14: 17}
EXPECTED_MAX_BASE_SUPPORTS = {
    13: ((1, 12),),
    14: ((1, 2), (1, 13)),
}
OUTSIDE_DEPTHS = ("degree+1", "degree+2", "2*degree+3")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


SOURCE = (
    Path(__file__).resolve().parents[1]
    / "04-computation"
    / "gmc_fixed_rank_thirty_five_cofactor_newton_atlas_thm2949.py"
)
SOURCE_BYTES = SOURCE.read_bytes().replace(b"\r\n", b"\n").replace(
    b"\r", b"\n"
)
SOURCE_SHA256 = (
    "9a1c7068e079e232dc97fd6eb925621aa74b3d636380a85995f8e0db8b30aa54"
)
require(
    sha256(SOURCE_BYTES).hexdigest() == SOURCE_SHA256,
    "THM-2949 exact dependency hash changed",
)
SPEC = importlib.util.spec_from_file_location("thm2949_width13_14", SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load THM-2949")
thm2949 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(thm2949)

require(
    thm2949.DELETE_ROW_POSITION == 30
    and thm2949.DELETE_GLOBAL_ROW == 36
    and thm2949.DELETE_COLUMN == 0
    and thm2949.DELETE_TARGET == (0, 0, 7),
    "THM-2949 fixed cofactor changed",
)


def integer_digest(values: list[int]) -> str:
    return sha256(",".join(map(str, values)).encode()).hexdigest()


def record_digest(records: list[str]) -> str:
    return sha256("\n".join(records).encode()).hexdigest()


def audit_family(task: tuple[int, int, int]) -> dict[str, object]:
    width, first, second = task
    offsets = (0, first, second, width)
    degree = 55 * width - 35
    rows = thm2949.forms_and_rows(width, first, second)

    values = [
        thm2949.fixed_cofactor_value(rows, depth)
        for depth in range(degree + 1)
    ]
    polynomial = thm2949.thm2943.interpolate(values)
    polynomial = thm2949.normalized_positive_leading(polynomial)
    require(
        polynomial.degree() == degree,
        f"fixed-cofactor degree changed: {offsets}",
    )

    outside = (degree + 1, degree + 2, 2 * degree + 3)
    for depth in outside:
        require(
            int(polynomial(depth))
            == thm2949.fixed_cofactor_value(rows, depth),
            f"outside-grid determinant mismatch: {offsets}, n={depth}",
        )

    newton = thm2949.newton_coefficients(
        [int(polynomial(depth)) for depth in range(degree + 1)]
    )
    base = None
    for candidate in range(4 * width + 1):
        if newton[0] > 0 and all(value >= 0 for value in newton):
            base = candidate
            break
        newton = thm2949.shift_newton(newton)
    require(base is not None, f"no Newton tail by base 4M: {offsets}")
    require(
        int(polynomial(base)) > 0,
        f"Newton-tail initial value is not positive: {offsets}",
    )

    low_signs = "".join(
        "+"
        if polynomial(depth) > 0
        else "-"
        if polynomial(depth) < 0
        else "0"
        for depth in range(base)
    )
    require("0" not in low_signs, f"integer-prefix zero: {offsets}")

    record = (
        f"family={offsets};degree={degree};"
        f"GN_base={base};orientation=+1;"
        f"low_signs={low_signs};"
        f"constant_sign={1 if polynomial(0)>0 else -1:+d};"
        f"newton_digest={integer_digest(newton)};"
        "outside_controls=3"
    )
    return {
        "width": width,
        "support": (first, second),
        "base": base,
        "low_signs": low_signs,
        "record": record,
    }


def main() -> None:
    tasks = [
        (width, first, second)
        for width in WIDTHS
        for first, second in combinations(range(1, width), 2)
    ]
    with ProcessPoolExecutor(max_workers=4) as executor:
        audited = list(executor.map(audit_family, tasks))

    all_records: list[str] = []
    total_outside_controls = 0
    summaries = []
    for width in WIDTHS:
        family_rows = [row for row in audited if row["width"] == width]
        records = [str(row["record"]) for row in family_rows]
        all_records.extend(records)
        digest = record_digest(records)
        require(
            digest == EXPECTED_WIDTH_DIGESTS[width],
            f"width-{width} record digest changed",
        )

        family_count = comb(width - 1, 2)
        require(len(records) == family_count, "support census changed")
        bases = [int(row["base"]) for row in family_rows]
        max_base = max(bases)
        base_zero = sum(base == 0 for base in bases)
        mixed_prefixes = sum(
            "+" in str(row["low_signs"]) and "-" in str(row["low_signs"])
            for row in family_rows
        )
        max_supports = tuple(
            row["support"]
            for row in family_rows
            if int(row["base"]) == max_base
        )
        require(
            max_base == EXPECTED_MAX_BASES[width],
            f"width-{width} maximum Newton base changed",
        )
        require(
            base_zero == EXPECTED_BASE_ZERO_COUNTS[width],
            f"width-{width} base-zero count changed",
        )
        require(
            mixed_prefixes == EXPECTED_MIXED_PREFIX_COUNTS[width],
            f"width-{width} mixed-prefix count changed",
        )
        require(
            max_supports == EXPECTED_MAX_BASE_SUPPORTS[width],
            f"width-{width} maximum-base supports changed",
        )
        total_outside_controls += 3 * family_count
        summaries.append(
            (
                width,
                family_count,
                55 * width - 35,
                max_base,
                base_zero,
                mixed_prefixes,
                max_supports,
                digest,
            )
        )

    require(len(audited) == 144, "total family count changed")
    require(
        total_outside_controls == 432,
        "outside-grid control count changed",
    )
    require(
        record_digest(all_records) == EXPECTED_GLOBAL_DIGEST,
        "global record digest changed",
    )

    hostile = next(
        row
        for row in audited
        if row["width"] == 14 and row["support"] == (1, 2)
    )
    require(
        hostile["base"] == 43
        and hostile["low_signs"]
        == "++++++++++++++++++++++++-------------------",
        "width-fourteen real-ray hostile changed",
    )

    print("THM-2952 FIXED RANK-35 COFACTOR WIDTH-13/14 CONTINUATION")
    print(f"thm2949_dependency_sha256={SOURCE_SHA256}")
    print(
        "fixed_delete=rowpos30/global36/column0/target(0,0,7);"
        f"outside_grid={','.join(OUTSIDE_DEPTHS)}"
    )
    for (
        width,
        families,
        degree,
        max_base,
        base_zero,
        mixed_prefixes,
        max_supports,
        digest,
    ) in summaries:
        print(
            f"width={width};families={families};degree={degree};"
            f"max_GN_base={max_base};base0={base_zero};"
            f"mixed_low_prefixes={mixed_prefixes};"
            f"max_base_supports={max_supports};"
            f"family_digest={digest}"
        )
    print("total_families=144;outside_grid_determinant_controls=432")
    print(f"global_record_digest={EXPECTED_GLOBAL_DIGEST}")
    print(
        "width14_real_ray_hostile=support:(0,1,2,14);"
        "integer_sign_word:+24,-19,+tail_from_43;"
        "two_positive_nonintegral_crossings"
    )
    print("all_fixed_cofactors_nonzero_at_every_integer_depth=PASS")
    print("consequence=first-window_SFC4_widths_13_14")
    print(
        "scope=four distinct slots;integer depths;no width15,"
        "real-ray positivity,arbitrary-width SFC4,NC2,or GMC2"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
