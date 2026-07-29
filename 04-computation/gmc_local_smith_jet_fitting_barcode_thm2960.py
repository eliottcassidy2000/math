#!/usr/bin/env python3
"""Exact companion for THM-2960.

For the first-gap support (0,1,2,M), 6 <= M <= 20, this script:

* recovers the local Smith bars of THM-2949's fixed 35-by-35 cofactor
  from lower block-Toeplitz jet nullities;
* separates the exact q200^5*c300 pure-coefficient contribution;
* proves divisibility by the corrected negative-depth factor B^Smith_M;
* verifies the two genuine matrix-level sporadic corrections; and
* compares both full 36-by-36 charts locally, proving on this finite bank
  that their common negative-depth content after q200^5*c300 is
  B^Smith_M*E_M.

No determinant is interpolated, expanded, or reduced modulo a prime.
"""

from __future__ import annotations

import importlib.util
import sys
from fractions import Fraction
from hashlib import sha256
from math import comb
from pathlib import Path

from flint import fmpq, fmpq_mat, fmpz_mat


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
SOURCE = HERE / "gmc_fixed_rank_thirty_five_cofactor_newton_atlas_thm2949.py"
SOURCE_BYTES = SOURCE.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
SOURCE_SHA256 = (
    "9a1c7068e079e232dc97fd6eb925621aa74b3d636380a85995f8e0db8b30aa54"
)
require(
    sha256(SOURCE_BYTES).hexdigest() == SOURCE_SHA256,
    "THM-2949 dependency changed",
)
SPEC = importlib.util.spec_from_file_location("thm2949_for_thm2960", SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot import THM-2949")
thm2949 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = thm2949
SPEC.loader.exec_module(thm2949)

t = thm2949.t
x = t.POLY_X
WIDTHS = tuple(range(6, 21))
FIXED_TRUNCATIONS = 5
FULL_TRUNCATIONS = 7

ORIGINAL_ROWS = tuple(thm2949.SELECTED)
ALTERNATE_ROWS = thm2949.thm2943.selection(thm2949.thm2943.MUTATED_F)
require(
    ALTERNATE_ROWS
    == tuple(sorted((set(ORIGINAL_ROWS) - {37, 38}) | {42, 43})),
    "stable row mutation changed",
)


# The tuple is (d_1,...,d_5), with d_k the nullity of multiplication
# by the local matrix in the kth truncated power-series module.
EXPECTED_PROFILES = {
    "half": (5, 5, 5, 5, 5),
    "one": (2, 4, 5, 6, 6),
    "two": (9, 17, 20, 21, 21),
    "early": (14, 23, 26, 26, 26),
    "middle_left": (14, 23, 24, 24, 24),
    "middle_right": (9, 17, 20, 20, 20),
    "high": (8, 15, 19, 19, 19),
    "terminal": (14, 19, 19, 19, 19),
    "q_resonance": (10, 18, 22, 24, 24),
    "c_resonance": (9, 16, 20, 20, 20),
    "sporadic_11_9": (9, 16, 20, 20, 20),
    "sporadic_12_5": (14, 23, 25, 25, 25),
}

EXPECTED_FIXED_GLOBAL_DIGEST = (
    "8be6356fe0ad6683de9433ff9d4bdce7c813534aaed6becabf70e9b0da0283ea"
)
EXPECTED_BRIDGE_GLOBAL_DIGEST = (
    "035a164a768deb837d2f483c707eb124519709ae45d5fdef0ec08e3547e41f2a"
)


def valuation(poly, factor) -> int:
    answer = 0
    while poly and poly % factor == 0:
        poly //= factor
        answer += 1
    return answer


def taylor_coefficient(poly, center, order: int):
    if not callable(poly):
        return int(poly) if order == 0 else 0
    coefficients = tuple(int(value) for value in poly.coeffs())
    return sum(
        (
            coefficients[degree]
            * comb(degree, order)
            * center ** (degree - order)
            for degree in range(order, len(coefficients))
        ),
        Fraction(0, 1) if isinstance(center, Fraction) else 0,
    )


def as_fmpq(value):
    value = Fraction(value)
    return fmpq(value.numerator, value.denominator)


def jet_nullities(rows, center, truncations: int) -> tuple[int, ...]:
    row_count = len(rows)
    column_count = len(rows[0])
    require(row_count == column_count, "square local chart required")
    size = row_count
    jets = [
        [
            [
                taylor_coefficient(rows[row][column], center, order)
                for column in range(size)
            ]
            for row in range(size)
        ]
        for order in range(truncations)
    ]
    answer = []
    for truncation in range(1, truncations + 1):
        block = [[0] * (size * truncation) for _ in range(size * truncation)]
        for out_degree in range(truncation):
            for in_degree in range(out_degree + 1):
                jet = jets[out_degree - in_degree]
                for row in range(size):
                    target = block[out_degree * size + row]
                    start = in_degree * size
                    for column in range(size):
                        target[start + column] = jet[row][column]
        matrix = (
            fmpq_mat([[as_fmpq(value) for value in row] for row in block])
            if isinstance(center, Fraction)
            else fmpz_mat(block)
        )
        answer.append(size * truncation - matrix.rank())
    return tuple(answer)


def smith_multiset(nullities: tuple[int, ...]) -> str:
    at_least = tuple(
        nullities[index] - (nullities[index - 1] if index else 0)
        for index in range(len(nullities))
    )
    exact = []
    for exponent in range(1, len(at_least)):
        count = at_least[exponent - 1] - at_least[exponent]
        if count:
            exact.append(f"{exponent}^{count}")
    require(at_least[-1] == 0, "Smith exponent exceeds jet truncation")
    return ",".join(exact)


def q_resonance(width: int, root: int) -> bool:
    return width % 6 == 1 and 3 * root == 2 * width + 1


def c_resonance(width: int, root: int) -> bool:
    return width % 4 == 1 and 4 * root == 3 * width + 1


def chamber(width: int, root: int) -> str:
    if root == 1:
        return "one"
    if root == 2:
        return "two"
    if root == width:
        return "terminal"
    if (width, root) == (11, 9):
        return "sporadic_11_9"
    if (width, root) == (12, 5):
        return "sporadic_12_5"
    if q_resonance(width, root):
        return "q_resonance"
    if c_resonance(width, root):
        return "c_resonance"
    if 3 * root <= width:
        return "early"
    if 2 * root <= width:
        return "middle_left"
    if 3 * root <= 2 * width:
        return "middle_right"
    return "high"


def floor_b_multiplicity(width: int, root: int) -> int:
    if root == 1:
        return 6
    if root == 2:
        return 21
    if 3 * root <= width:
        return 26
    if 2 * root <= width:
        return 24
    if 3 * root <= 2 * width:
        return 20
    return 19


def b_multiplicity(width: int, root: int) -> int:
    return floor_b_multiplicity(width, root) + int(
        (width, root) in ((11, 9), (12, 5))
    )


def e_multiplicity(width: int, root: int) -> int:
    if root == 1:
        return 0
    if root == width:
        return 2
    return 3 if 2 * root <= width else 4


def expected_q_valuation(width: int, root: int) -> int:
    return int(q_resonance(width, root))


def expected_c_valuation(width: int, root: int) -> int:
    return int(c_resonance(width, root))


def polynomial_data(width: int):
    source_rows = thm2949.forms_and_rows(width, 1, 2)
    cofactor = tuple(
        tuple(
            source_rows[row][column]
            for column in thm2949.COFACTOR_COLUMNS
        )
        for row in thm2949.COFACTOR_ROWS
    )
    original = tuple(tuple(source_rows[row]) for row in ORIGINAL_ROWS)
    alternate = tuple(tuple(source_rows[row]) for row in ALTERNATE_ROWS)
    forms = thm2949.thm2943.polynomial_forms(width, (0, 1, 2))
    q200 = forms[0][(2, 0, 0)]
    c300 = forms[1][(3, 0, 0)]
    return cofactor, original, alternate, q200, c300


def audit_fixed_width(width: int, data) -> str:
    cofactor, _original, _alternate, q200, c300 = data
    half_profile = jet_nullities(cofactor, Fraction(-1, 2), FIXED_TRUNCATIONS)
    require(half_profile == EXPECTED_PROFILES["half"], "half profile changed")
    require(
        valuation(q200, 2 * x + 1) == 0
        and valuation(c300, 2 * x + 1) == 0,
        "pure coefficient acquired a half-root",
    )
    records = [
        f"half:{half_profile}:{smith_multiset(half_profile)}:B=5:q=0:c=0"
    ]
    chamber_counts = {name: 0 for name in EXPECTED_PROFILES}
    chamber_counts["half"] = 1
    b_degree = 5

    for root in range(1, width + 1):
        name = chamber(width, root)
        profile = jet_nullities(cofactor, -root, FIXED_TRUNCATIONS)
        require(
            profile == EXPECTED_PROFILES[name],
            f"Smith profile changed: M={width},r={root},got={profile}",
        )
        factor = x + root
        q_order = valuation(q200, factor)
        c_order = valuation(c300, factor)
        require(
            q_order == expected_q_valuation(width, root)
            and c_order == expected_c_valuation(width, root),
            f"pure root law changed: M={width},r={root}",
        )
        beta_floor = floor_b_multiplicity(width, root)
        beta = b_multiplicity(width, root)
        require(
            profile[-1] == beta + 5 * q_order + c_order,
            f"local order decomposition failed: M={width},r={root}",
        )
        b_degree += beta
        chamber_counts[name] += 1
        records.append(
            f"{root}:{name}:{profile}:{smith_multiset(profile)}:"
            f"B={beta}:Bfloor={beta_floor}:q={q_order}:c={c_order}"
        )

    floor_degree = (
        19 * width
        + 2 * (width // 3)
        + 4 * (width // 2)
        + (2 * width) // 3
        - 20
    )
    require(
        b_degree == floor_degree + int(width in (11, 12)),
        f"corrected B degree law changed: M={width}",
    )
    digest = sha256("\n".join(records).encode()).hexdigest()
    nonzero_counts = ",".join(
        f"{name}:{count}"
        for name, count in chamber_counts.items()
        if count
    )
    q_root = next(
        (root for root in range(1, width + 1) if q_resonance(width, root)),
        0,
    )
    c_root = next(
        (root for root in range(1, width + 1) if c_resonance(width, root)),
        0,
    )
    return (
        f"M={width};centres={width + 1};B_degree={b_degree};"
        f"B_floor_degree={floor_degree};"
        f"q_resonance={q_root};c_resonance={c_root};"
        f"chambers={nonzero_counts};record_digest={digest};"
        f"sporadic_correction={int(width in (11, 12))};"
        "local_order_decomposition=PASS"
    )


def certifies_order(profile: tuple[int, ...], target: int) -> bool:
    return profile[-1] == target and profile[-1] == profile[-2]


def profile_to_target(rows, center, target: int) -> tuple[int, ...]:
    profile = jet_nullities(rows, center, FIXED_TRUNCATIONS)
    if profile[-1] < target or (
        profile[-1] == target and profile[-1] != profile[-2]
    ):
        profile = jet_nullities(rows, center, FULL_TRUNCATIONS)
    return profile


def audit_bridge_width(width: int, data) -> str:
    _cofactor, original, alternate, q200, c300 = data
    records = []

    original_half = profile_to_target(original, Fraction(-1, 2), 6)
    alternate_half = profile_to_target(alternate, Fraction(-1, 2), 6)
    require(
        original_half[-1] >= 6
        and alternate_half[-1] >= 6
        and (
            certifies_order(original_half, 6)
            or certifies_order(alternate_half, 6)
        ),
        f"half common-order certificate changed: M={width}",
    )
    records.append(
        f"half:P{original_half}:A{alternate_half}:common=6:B=5:E=1"
    )

    for root in range(1, width + 1):
        q_order = valuation(q200, x + root)
        c_order = valuation(c300, x + root)
        beta = b_multiplicity(width, root)
        epsilon = e_multiplicity(width, root)
        common_order = 5 * q_order + c_order + beta + epsilon
        original_profile = profile_to_target(original, -root, common_order)
        alternate_profile = profile_to_target(alternate, -root, common_order)
        require(
            original_profile[-1] >= common_order
            and alternate_profile[-1] >= common_order
            and (
                certifies_order(original_profile, common_order)
                or certifies_order(alternate_profile, common_order)
            ),
            "two-chart/35-minor seam mismatch: "
            f"M={width},r={root},P={original_profile},A={alternate_profile}",
        )
        records.append(
            f"{root}:P{original_profile}:A{alternate_profile}:"
            f"common={common_order}:q={q_order}:c={c_order}:"
            f"B={beta}:E={epsilon}"
        )

    digest = sha256("\n".join(records).encode()).hexdigest()
    return (
        f"M={width};centres={width + 1};half_common=6;"
        f"integer_common_after_q5c=B+E;record_digest={digest}"
    )


def barcode_control() -> str:
    exponents = (0, 1, 2, 4)
    nullities = tuple(
        sum(min(truncation, exponent) for exponent in exponents)
        for truncation in range(1, 6)
    )
    require(nullities == (3, 5, 6, 7, 7), "barcode control changed")
    at_least = tuple(
        nullities[index] - (nullities[index - 1] if index else 0)
        for index in range(len(nullities))
    )
    require(at_least == (3, 2, 1, 1, 0), "barcode recovery changed")
    return "barcode_control=alpha(0,1,2,4);d(3,5,6,7,7);recovery=PASS"


def row_rank_hostile(data_by_width) -> str:
    cofactor, _original, _alternate, _q200, _c300 = data_by_width[20]
    forms = thm2949.thm2943.polynomial_forms(20, (0, 1, 2))
    factor = x + 3
    form_minima = tuple(
        min(valuation(poly, factor) for poly in form.values())
        for form in forms
    )
    require(form_minima == (0, 0, 0), "row-gcd hostile changed")
    profile = jet_nullities(cofactor, -3, FIXED_TRUNCATIONS)
    require(profile[:3] == (14, 23, 26), "rank hostile changed")
    return (
        "row_rank_hostile=M20,r3;form_row_minima=0,0,0;"
        "ordinary_corank=14;det_order=26;higher_jets_are_essential=YES"
    )


def main() -> None:
    print("THM-2960 LOCAL SMITH-JET FITTING BARCODE")
    print(f"thm2949_dependency_sha256={SOURCE_SHA256}")
    print(
        "widths=6..20;fixed_matrix_size=35;full_chart_size=36;"
        "characteristic=0"
    )
    print(barcode_control())

    data_by_width = {width: polynomial_data(width) for width in WIDTHS}
    fixed_records = []
    for width in WIDTHS:
        record = audit_fixed_width(width, data_by_width[width])
        fixed_records.append(record)
        print("fixed;" + record, flush=True)
    fixed_digest = sha256("\n".join(fixed_records).encode()).hexdigest()
    require(
        fixed_digest == EXPECTED_FIXED_GLOBAL_DIGEST,
        "fixed-cofactor global record digest changed",
    )
    print("fixed_global_record_digest=" + fixed_digest)

    bridge_records = []
    for width in WIDTHS:
        record = audit_bridge_width(width, data_by_width[width])
        bridge_records.append(record)
        print("bridge;" + record, flush=True)
    bridge_digest = sha256("\n".join(bridge_records).encode()).hexdigest()
    require(
        bridge_digest == EXPECTED_BRIDGE_GLOBAL_DIGEST,
        "two-chart bridge global record digest changed",
    )
    print("bridge_global_record_digest=" + bridge_digest)

    print(
        "q_hostile=M7,r5;q_order=1;c_order=0;"
        "floor_B=19;smith_order=24"
    )
    print(
        "c_hostile=M9,r7;q_order=0;c_order=1;"
        "floor_B=19;smith_order=20"
    )
    print("matrix_sporadics=(11,9):+1,(12,5):+1")
    print(row_rank_hostile(data_by_width))
    print(
        "E_M=(2n+1)*prod(r=2..floor(M/2),(n+r)^3)*"
        "prod(r=floor(M/2)+1..M-1,(n+r)^4)*(n+M)^2"
    )
    print("determinant_interpolation_used=NO;finite_field_inference_used=NO")
    print(
        "scope=first-gap M6..20 local factors only;"
        "no arbitrary-width core nonvanishing"
    )
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
