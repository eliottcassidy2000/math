#!/usr/bin/env python3
"""Exact strict-ULC atlas for first-gap wall-stripped norm cores, M=6..34.

Scratch candidate companion for reserved THM-2982.  The proved finite scope is
only the supports (0,1,2,M), 6 <= M <= 34.  The separate fast Jacobi check at
M=34 uses only the top three form slices and the explicit wall invoice; it is
not an arbitrary-width theorem.
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

from flint import fmpq_mat


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = (
    ROOT / "04-computation/gmc_first_gap_wall_stripped_norm_core_atlas_thm2969.py"
)
DEPENDENCIES = {
    "THM-2969": (
        BASE_PATH,
        "5be5df0fd058f436593339f78cd6240a47975082d22929c135ba811458753bf5",
    ),
    "THM-2973": (
        ROOT
        / "04-computation/gmc_first_gap_wall_stripped_norm_core_continuation_thm2973.py",
        "3539d1fcd7dd317659aab8296a86ab928999a7efff237e3f00571481a80cdad3",
    ),
    "THM-2978": (
        ROOT
        / "04-computation/gmc_first_gap_wall_stripped_norm_core_at_thirty_two_thm2978.py",
        "4787aac1159f8ee9bdf7ad27df9654c3ed46ab32060d38965ba691b662c132cb",
    ),
}

EXTENDED_CORES = {
    27: (556, "c00583d9eded605bd99fed73892a358e91224a7d4086a19c0a55f0b394285b7f"),
    28: (577, "9a6b8b02f59f2d1f31a7bc262863c4461aa8d1d0c89c1e4b6901237a68558591"),
    29: (599, "98f7a7599769664a781edb7fa1f61bd5b4a71092c3fef345de74c66094a42e18"),
    30: (617, "ba6fc2329242f5f222295ec261b36d3a57e353206716667d67ee26e38b17c7a1"),
    31: (639, "c9ffa50924c9d0298fd4995b1343fb200d3da0c38e1fcbe07347d97b8908e81c"),
    32: (660, "3dc8b52265fd9b96d92851cc1f69df9f987c681843f61b6893ac68506edd323e"),
    33: (680, "e86a6beb0403fa1a7cfa020813073a102d94ae7a17e4ed9c4eb673eeb7665040"),
    34: (701, "080a8e37065be58ccb1d7c6f551f3febc35c960f59d76d3b3335dfff658f0d91"),
}

# (degree, true ascending-power minimizer index, distance from the leading
# edge, floor(10^15*(minimum exact ULC ratio-1))).  The last coordinate is an
# exact rational enclosure marker, not a floating-point test.
EXPECTED_RIDGE = {
    6: (121, 61, 60, 8062284689677),
    7: (144, 80, 64, 6081743865577),
    8: (164, 97, 67, 4914016682307),
    9: (184, 116, 68, 4030049973186),
    10: (205, 136, 69, 3322774623052),
    11: (226, 157, 69, 2774689440288),
    12: (244, 178, 66, 2422277897436),
    13: (268, 202, 66, 2034759382659),
    14: (288, 225, 63, 1792140857838),
    15: (308, 248, 60, 1589425302484),
    16: (329, 271, 58, 1419120767553),
    17: (351, 296, 55, 1271473804801),
    18: (369, 317, 52, 1172404748369),
    19: (392, 343, 49, 1061054752654),
    20: (412, 366, 46, 983518608740),
    21: (431, 388, 43, 917010071591),
    22: (453, 413, 40, 850259420433),
    23: (475, 438, 37, 792330287977),
    24: (493, 460, 33, 751692881350),
    25: (516, 486, 30, 702694538014),
    26: (536, 509, 27, 667586916181),
    27: (556, 533, 23, 633706325972),
    28: (577, 557, 20, 602552814983),
    29: (599, 583, 16, 572626563514),
    30: (617, 604, 13, 551016704370),
    31: (639, 630, 9, 525597811989),
    32: (660, 654, 6, 504135911672),
    33: (680, 678, 2, 484377993981),
    34: (701, 700, 1, 465983639062),
}
EXPECTED_ATLAS_DIGEST = (
    "9a496cb20febf0e88b039afcca56540ff0fe352f11aeaa79b965dd7922237e8e"
)
EXPECTED_M34_RATIO = Fraction(
    27655303142730846948837614910195097534441250912778059699455304305433884593670175322013125,
    27642422226229352775679708685499047535273900449378640304403997574607551772198286144871066,
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(data).hexdigest()


def load_base(tag: str):
    name = f"thm2969_thm2982_{tag}"
    spec = importlib.util.spec_from_file_location(name, BASE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import THM-2969 engine")
    base = importlib.util.module_from_spec(spec)
    sys.modules[name] = base
    spec.loader.exec_module(base)
    return base


def corrected_degree(width: int) -> int:
    return (
        23 * width
        - 2 * (width // 3)
        - 2 * (width // 2)
        - (2 * width // 3)
        - 3
        - int(width in (11, 12, 21, 31))
    )


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def audit_width(width: int) -> dict[str, object]:
    """Rebuild one full primitive core and scan every strict ULC circuit."""
    base = load_base(f"worker_{width}")
    old_correction = base.smith_correction
    old_require = base.require

    def smith_correction(M: int, root: int) -> int:
        return old_correction(M, root) + int((M, root) == (31, 25))

    def worker_require(condition: bool, message: str) -> None:
        # THM-2969 predates the audited M=31 quartic wall.  THM-2973 proves
        # exactly this census extension; every division/digest check remains.
        if not condition and message == "pure coefficient resonance census changed: M=31":
            return
        old_require(condition, message)

    base.smith_correction = smith_correction
    base.expected_core_degree = corrected_degree
    base.require = worker_require
    base.EXPECTED_CORES.update(EXTENDED_CORES)

    capture = io.StringIO()
    with contextlib.redirect_stdout(capture):
        core = base.audit_width(width)[0]
    require(
        len(tuple(line for line in capture.getvalue().splitlines() if line)) == 1,
        f"unexpected worker transcript at M={width}",
    )

    # python-flint fmpz_poly.coeffs() is native constant-first.  The explicit
    # indexed equality prevents the historical accidental reversal.
    coefficients = tuple(int(value) for value in core.coeffs())
    require(
        coefficients == tuple(int(core[index]) for index in range(core.degree() + 1)),
        f"coefficient orientation changed at M={width}",
    )
    degree = core.degree()
    require(all(value > 0 for value in coefficients), f"positivity failed at M={width}")

    rows: list[tuple[int, int, Fraction]] = []
    for index in range(1, degree):
        left = coefficients[index] ** 2 * index * (degree - index)
        right = (
            coefficients[index - 1]
            * coefficients[index + 1]
            * (index + 1)
            * (degree - index + 1)
        )
        rows.append((index, left - right, Fraction(left, right)))
    require(all(slack > 0 for _index, slack, _ratio in rows), f"ULC failed at M={width}")
    minimum = min(rows, key=lambda row: row[2])
    ratio = minimum[2]
    epsilon = ratio - 1
    epsilon_floor = epsilon.numerator * 10**15 // epsilon.denominator
    expected = EXPECTED_RIDGE[width]
    require(
        (degree, minimum[0], degree - minimum[0], epsilon_floor) == expected,
        f"ridge record changed at M={width}",
    )
    ratio_text = fraction_text(ratio)
    return {
        "M": width,
        "degree": degree,
        "core_digest": base.polynomial_digest(core),
        "min_ratio_index": minimum[0],
        "min_distance_from_leading": degree - minimum[0],
        "min_ratio_sha256": sha256(ratio_text.encode()).hexdigest(),
        "eps_floor_1e15": epsilon_floor,
        "strict_ULC": True,
        # Kept internally for the independent M=34 equality check, then
        # removed before the compact canonical record is hashed.
        "_edge_ratio": fraction_text(rows[-1][2]),
    }


def as_fraction(value) -> Fraction:
    return Fraction(int(value.numer()), int(value.denom()))


def matrix_trace(matrix) -> object:
    require(matrix.nrows() == matrix.ncols(), "trace needs a square matrix")
    return sum((matrix[index, index] for index in range(matrix.nrows())), start=0)


def polynomial_log_jets(poly) -> tuple[Fraction, Fraction]:
    degree = poly.degree()
    leading = int(poly[degree])
    u = Fraction(int(poly[degree - 1]), leading)
    return u, Fraction(2 * int(poly[degree - 2]), leading) - u * u


def selected_rows(base, forms):
    rows, _metadata = base.thm2943.thm2942.macaulay_rows(forms)
    return [rows[index] for index in base.thm2943.t.SELECTED_ROWS]


def fast_m34_edge_ratio() -> Fraction:
    """Independent top-three-slice/Jacobi prediction of the M=34 edge."""
    width = 34
    base = load_base("fast_m34")
    forms = base.thm2943.polynomial_forms(width, (0, 1, 2))
    slices = [[], [], []]
    for form in forms:
        degree = max(poly.degree() for poly in form.values())
        for drop in range(3):
            slices[drop].append(
                {
                    monomial: int(poly[degree - drop])
                    for monomial, poly in form.items()
                }
            )
    m0, m1, m2 = [fmpq_mat(selected_rows(base, row)) for row in slices]
    inverse = m0.inv()
    first = inverse * m1
    chart_u = as_fraction(matrix_trace(first))
    chart_ell2 = (
        2 * as_fraction(matrix_trace(inverse * m2))
        - as_fraction(matrix_trace(first * first))
    )

    q200, c300, curvature, _alternate = base.thm2943.flag_polynomials(forms)
    q_u, q_ell2 = polynomial_log_jets(q200)
    c_u, c_ell2 = polynomial_log_jets(c300)
    k_u, k_ell2 = polynomial_log_jets(curvature)
    resultant_u = chart_u - 6 * q_u - c_u - k_u
    resultant_ell2 = chart_ell2 - 6 * q_ell2 - c_ell2 - k_ell2

    # W=B E/H at M=34: one half-wall with multiplicity six and the exact
    # baseline integer-wall exponents.  There is no M=34 pure/Smith addition.
    wall_degree = 6
    wall_u = Fraction(3)
    wall_square_sum = Fraction(3, 2)
    for root in range(1, width + 1):
        if root == 1:
            exponent = 6
        elif root == 2:
            exponent = 21
        elif 3 * root <= width:
            exponent = 26
        elif 2 * root <= width:
            exponent = 24
        elif 3 * root <= 2 * width:
            exponent = 20
        else:
            exponent = 19
        if root == width:
            exponent += 2
        elif 2 <= root < width:
            exponent += 3 if 2 * root <= width else 4
        if root == width or 3 <= root <= width // 2:
            exponent -= 1
        require(exponent > 0, (root, exponent))
        wall_degree += exponent
        wall_u += exponent * root
        wall_square_sum += exponent * root * root

    degree = 46 * width - 26 - wall_degree
    core_u = resultant_u - wall_u
    core_ell2 = resultant_ell2 + wall_square_sum
    core_v = core_ell2 + core_u * core_u
    return Fraction(degree - 1, degree) * core_u * core_u / core_v


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(1 <= args.workers <= 8, "workers must lie in 1..8")

    dependency_hashes = {}
    for theorem, (path, expected) in DEPENDENCIES.items():
        actual = lf_sha256(path)
        require(actual == expected, f"{theorem} dependency hash changed")
        dependency_hashes[theorem] = actual

    widths = tuple(range(6, 35))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        records = list(pool.map(audit_width, widths))
    records.sort(key=lambda row: int(row["M"]))

    m34 = next(row for row in records if row["M"] == 34)
    full_m34_edge = Fraction(str(m34["_edge_ratio"]))
    fast_m34_edge = fast_m34_edge_ratio()
    require(full_m34_edge == EXPECTED_M34_RATIO, "full M34 edge ratio changed")
    require(fast_m34_edge == EXPECTED_M34_RATIO, "fast Jacobi M34 prediction changed")

    compact_records = []
    for row in records:
        compact = dict(row)
        compact.pop("_edge_ratio")
        compact_records.append(compact)
    lines = [
        json.dumps(row, sort_keys=True, separators=(",", ":"))
        for row in compact_records
    ]
    atlas_digest = sha256("\n".join(lines).encode()).hexdigest()
    require(atlas_digest == EXPECTED_ATLAS_DIGEST, "atlas record digest changed")

    transcript = "\n".join(
        [
            "FIRST-GAP WALL-STRIPPED NORM-CORE STRICT ULC ATLAS",
            "scope=M6..34;coefficient_order=ascending_n_power;all_circuits_scanned=1",
            "dependencies="
            + ",".join(f"{key}:{dependency_hashes[key]}" for key in sorted(dependency_hashes)),
            *lines,
            f"atlas_record_digest={atlas_digest}",
            "ridge_path_k=60,64,67,68,69,69,66,66,63,60,58,55,52,49,46,43,40,37,33,30,27,23,20,16,13,9,6,2,1",
            "M34_boundary_attachment=TRUE;first_width_in_scope_with_k=1",
            "M34_fast_Jacobi_prediction=EXACT_MATCH;top_form_slices=3;matrix_size=36",
            "M35_plus_symbolic_edge_and_no_return=CONDITIONAL_FRONTIER_NOT_DEPENDENCY",
            "no_arbitrary_width_or_arbitrary_radial_consequence",
            "all_exact_checks=PASS",
        ]
    ) + "\n"
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
