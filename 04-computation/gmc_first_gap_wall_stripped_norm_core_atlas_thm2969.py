#!/usr/bin/env python3
"""Finite exact atlas for the first-gap wall-stripped norm core.

Canonical finite-exact scope:

    support (0,1,2,M),  6 <= M <= 26.

For the two stable full Macaulay charts P_M,A_M, the script verifies

    gcd(P_M,A_M)
      ~ q200^5*c300*R_M*H_M,

where R_M is the genuine specialized ternary resultant and

    H_M=(n+M)*prod_{r=3..floor(M/2)}(n+r).

It then removes the local Smith factor B_M and the full-chart seam E_M.
The resulting primitive core is checked to equal the pure resultant
after its complementary wall factor is removed.  All assertions are
finite exact assertions over ZZ/QQ; no finite-field inference is made.

The script also records a sharp stopping hostile to a naive M -> M+6
shift/divisibility recurrence, and a calibrated-chart resultant control.
The latter bounds common prime divisors of the *calibrated* chart values;
it does not produce a fixed prime bank for the original minors.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import sys
from hashlib import sha256
from pathlib import Path


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "04-computation"
    / "gmc_width_seven_eight_two_chart_resultant_closure_thm2943.py"
)
SOURCE_BYTES = SOURCE.read_bytes().replace(b"\r\n", b"\n").replace(
    b"\r", b"\n"
)
SOURCE_SHA256 = sha256(SOURCE_BYTES).hexdigest()
EXPECTED_SOURCE_SHA256 = (
    "d2f8afeba7dd6c7950405a4845d7bf112b6c9872dd8161146446be8bbdaae0ba"
)
require(
    SOURCE_SHA256 == EXPECTED_SOURCE_SHA256,
    "THM-2943 exact dependency hash changed",
)
SPEC = importlib.util.spec_from_file_location("thm2943_norm_core_atlas", SOURCE)
require(SPEC is not None and SPEC.loader is not None, "cannot load THM-2943")
thm2943 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = thm2943
SPEC.loader.exec_module(thm2943)

X = thm2943.X


EXPECTED_CORES = {
    6: (
        121,
        "c984248b9099d5a83637ab8487122c453505b181e38c79750328243a216ecfc3",
    ),
    7: (
        144,
        "cbdb0d648d29917073b4590b526fa3d981f6a15fc906a4ee4cca263dfd31fb26",
    ),
    8: (
        164,
        "6dc83f2406f9d575973baacaeb34fe23dd2b5efc772a7cc02a490cd90b12e749",
    ),
    9: (
        184,
        "2be4e9988411aeb920a276f2e78303d83b839161cbaa783da0ba482891115e93",
    ),
    10: (
        205,
        "f9c6698e52e6df4f88e0e4c79655d91f9f4e0483a505c7203e7d0b3174c79835",
    ),
    11: (
        226,
        "00744ab8a8830097754c40c470e1deadab2c88a028081aac1182fc651a344bdc",
    ),
    12: (
        244,
        "25d309de2b271884fd825a0b9042b546e5aee8b85d030912e775c6f1c0a983cd",
    ),
    13: (
        268,
        "deb7f4951bb11fc7953077536c909900837f3b97d5389443002b9c9553526b98",
    ),
    14: (
        288,
        "5ab8a419af6a3088dcbef1edbc6625b88d791ffdadfb632d3ceb0be4d4053d23",
    ),
    15: (
        308,
        "3695f4d1db417f8c4f3f5ea939741889d568d37acf1313b965c79a85884a4bee",
    ),
    16: (
        329,
        "f2fdb558e32aa1625c35a4494eb07de261ead0f0490aad2172c912b85827d444",
    ),
    17: (
        351,
        "07fe2a1065d2572f0fcff835832d37d3cc8523d1c9fa509f78342082f32debc4",
    ),
    18: (
        369,
        "71ce4e997650c189ec332a2b61b33e836456cd72c6c5f00f7a8a6b9ff9144bd7",
    ),
    19: (
        392,
        "59c7a8314c96b6325d7f2280ff8c96dd6e035bf78139045b710bfa958b8bfb93",
    ),
    20: (
        412,
        "5302b208b689d6dd4378494ffce5f3c8d2c60225006ff3a7c5d4380852737fc2",
    ),
    21: (
        431,
        "c188ab377c593c1a1035446fddfe0a8443fd96a46f58b51241d9af89f1f29475",
    ),
    22: (
        453,
        "44d0182663549969d3ff17223a06e76d94d2c8ff732f70a63db5be9480c9d153",
    ),
    23: (
        475,
        "75c487966c9c32eccb90b2d711ec591425087260ea2fffe485b45ee0b71758cb",
    ),
    24: (
        493,
        "b9bd1d69db049097fb20ce3b6f77f8163b68d3e0036f328bb4a0ca2ceaea7c1a",
    ),
    25: (
        516,
        "6833241cbec29b9bc930a0573dc963ce143341890d159bd1ab11b2eacf4b2383",
    ),
    26: (
        536,
        "fba227bd0bb3d3e5621b6bab15dd44b35a21de990cc828dcab2960b3ea24e068",
    ),
}

EXPECTED_CORE_RECORD_DIGEST = (
    "29bcdeade4601f818730957a6a326fa83c1cc1cb17df6189cda6c1082bbab255"
)
EXPECTED_FLAG_RECORD_DIGEST = (
    "a73fcfc2ed9a41ba86176c8e6409ac3f86d2f69df2f17f7584f74ae055e4859b"
)
EXPECTED_M6_CALIBRATED_RESULTANT = {
    "bits": 2380,
    "digits": 717,
    "digest": (
        "68785d6e7d64e3a73786911729c1adad14650349b5672e69818b309a9902d3a1"
    ),
}


def primitive_positive(poly):
    if hasattr(poly, "denom") and poly.denom() != 1:
        poly = poly.numer()
    elif hasattr(poly, "numer"):
        poly = poly.numer()
    require(poly != 0, "cannot normalize zero polynomial")
    if poly.leading_coefficient() < 0:
        poly = -poly
    content = poly.content()
    require(content > 0, "nonpositive polynomial content")
    return poly // content


def polynomial_digest(poly) -> str:
    coefficients = tuple(int(value) for value in primitive_positive(poly).coeffs())
    return sha256(str(coefficients).encode()).hexdigest()


def associated(left, right) -> bool:
    return primitive_positive(left) == primitive_positive(right)


def valuation(poly, factor) -> int:
    answer = 0
    while poly % factor == 0:
        poly //= factor
        answer += 1
    return answer


def floor_beta(width: int, root: int) -> int:
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


def smith_correction(width: int, root: int) -> int:
    # (11,9) and (21,17) are the pure f400 resonances.  (12,5) is
    # the separate matrix-level Smith sporadic.
    return int((width, root) in ((11, 9), (12, 5), (21, 17)))


def smith_factor(width: int):
    answer = (2 * X + 1) ** 5
    for root in range(1, width + 1):
        exponent = floor_beta(width, root) + smith_correction(width, root)
        answer *= (X + root) ** exponent
    return answer


def seam_factor(width: int):
    half = width // 2
    answer = 2 * X + 1
    for root in range(2, width):
        answer *= (X + root) ** (3 if root <= half else 4)
    answer *= (X + width) ** 2
    return answer


def expected_flag(width: int):
    answer = X + width
    for root in range(3, width // 2 + 1):
        answer *= X + root
    return answer


def expected_core_degree(width: int) -> int:
    correction = int(width in (11, 12, 21))
    return (
        23 * width
        - 2 * (width // 3)
        - 2 * (width // 2)
        - (2 * width // 3)
        - 3
        - correction
    )


def is_pf2(coefficients: tuple[int, ...]) -> bool:
    return all(
        coefficients[index] ** 2
        >= coefficients[index - 1] * coefficients[index + 1]
        for index in range(1, len(coefficients) - 1)
    )


def audit_width(width: int):
    forms = thm2943.polynomial_forms(width, (0, 1, 2))
    rows, _metadata = thm2943.thm2942.macaulay_rows(forms)
    degree_bound = 58 * width - 36
    original, mutated = thm2943.interpolate_pair(rows, degree_bound)
    require(
        original.degree() == degree_bound and mutated.degree() == degree_bound,
        f"full chart degree changed: M={width}",
    )

    for depth in range(degree_bound + 1, degree_bound + 4):
        numeric = thm2943.evaluate_rows(rows, depth)
        require(
            thm2943.determinant(numeric, thm2943.ORIGINAL_F) == original(depth)
            and thm2943.determinant(numeric, thm2943.MUTATED_F)
            == mutated(depth),
            f"outside interpolation control failed: M={width},n={depth}",
        )

    q200, c300, curvature, alternate = thm2943.flag_polynomials(forms)
    require(
        original % (q200**6 * c300 * curvature) == 0
        and mutated % (q200**5 * c300 * alternate) == 0,
        f"THM-2942 chart factorization failed: M={width}",
    )
    resultant = original // (q200**6 * c300 * curvature)
    require(
        resultant == mutated // (q200**5 * c300 * alternate),
        f"specialized resultants disagree: M={width}",
    )
    require(
        resultant.degree() == 46 * width - 26,
        f"specialized resultant degree changed: M={width}",
    )

    flag = primitive_positive((q200 * curvature).gcd(alternate))
    flag_expected = expected_flag(width)
    require(
        flag == primitive_positive(flag_expected),
        f"first-gap common flag formula changed: M={width}",
    )

    raw_common = original.gcd(mutated)
    require(
        associated(raw_common, q200**5 * c300 * resultant * flag),
        f"raw common gcd factorization changed: M={width}",
    )
    expected_common_degree = degree_bound - ((9 * width - 5) // 2)
    require(
        raw_common.degree() == expected_common_degree,
        f"raw common gcd degree changed: M={width}",
    )

    b_factor = smith_factor(width)
    e_factor = seam_factor(width)
    divisor = q200**5 * c300 * b_factor * e_factor
    require(
        raw_common % divisor == 0,
        f"q5*c*B*E division failed: M={width}",
    )
    core = primitive_positive(raw_common // divisor)

    require(
        (b_factor * e_factor) % flag_expected == 0,
        f"flag is not absorbed by the local walls: M={width}",
    )
    resultant_wall = (b_factor * e_factor) // flag_expected
    require(
        resultant % resultant_wall == 0
        and associated(core, resultant // resultant_wall),
        f"core is not the wall-stripped pure resultant: M={width}",
    )

    coefficients = tuple(int(value) for value in core.coeffs())
    require(
        coefficients
        and all(value > 0 for value in coefficients)
        and is_pf2(coefficients),
        f"strict coefficient positivity/PF2 failed: M={width}",
    )
    expected_degree, expected_digest = EXPECTED_CORES[width]
    require(
        core.degree() == expected_degree == expected_core_degree(width)
        and polynomial_digest(core) == expected_digest,
        f"core degree/digest changed: M={width}",
    )

    f400 = forms[2][(4, 0, 0)]
    q200_roots = tuple(
        root for root in range(1, width + 1) if q200 % (X + root) == 0
    )
    c300_roots = tuple(
        root for root in range(1, width + 1) if c300 % (X + root) == 0
    )
    f400_roots = tuple(
        root for root in range(1, width + 1) if f400 % (X + root) == 0
    )
    expected_q200_roots = (
        ((2 * width + 1) // 3,) if width % 6 == 1 else ()
    )
    expected_c300_roots = (
        ((3 * width + 1) // 4,) if width % 4 == 1 else ()
    )
    expected_f400_roots = (
        (9,) if width == 11 else (17,) if width == 21 else ()
    )
    require(
        q200_roots == expected_q200_roots
        and c300_roots == expected_c300_roots
        and f400_roots == expected_f400_roots,
        f"pure coefficient resonance census changed: M={width}",
    )

    if width == 25:
        root17 = X + 17
        root19 = X + 19
        require(
            (
                valuation(q200, root17),
                valuation(c300, root17),
                valuation(b_factor, root17),
                valuation(e_factor, root17),
                valuation(flag, root17),
                valuation(resultant_wall, root17),
                valuation(raw_common, root17),
                valuation(divisor, root17),
                valuation(core, root17),
            )
            == (1, 0, 19, 4, 0, 23, 28, 28, 0)
            and (
                valuation(q200, root19),
                valuation(c300, root19),
                valuation(b_factor, root19),
                valuation(e_factor, root19),
                valuation(flag, root19),
                valuation(resultant_wall, root19),
                valuation(raw_common, root19),
                valuation(divisor, root19),
                valuation(core, root19),
            )
            == (0, 1, 19, 4, 0, 23, 24, 24, 0),
            "M=25 simultaneous q/c wall invoice changed",
        )

    calibrated_original = primitive_positive(original // raw_common)
    calibrated_mutated = primitive_positive(mutated // raw_common)
    require(
        calibrated_original.gcd(calibrated_mutated).degree() == 0,
        f"calibrated charts are not coprime: M={width}",
    )

    core_record = (
        f"M={width};degree={core.degree()};digest={polynomial_digest(core)};"
        "strict_positive=1;PF2=1"
    )
    flag_roots = tuple(range(3, width // 2 + 1)) + (width,)
    flag_record = (
        f"M={width};degree={flag.degree()};"
        f"roots={','.join(map(str, flag_roots))}"
    )
    print(
        f"core;{core_record};resultant_degree={resultant.degree()};"
        f"flag_degree={flag.degree()};calibrated_degrees="
        f"{calibrated_original.degree()},{calibrated_mutated.degree()};"
        f"q200_roots={','.join(map(str, q200_roots)) or 'NONE'};"
        f"c300_roots={','.join(map(str, c300_roots)) or 'NONE'};"
        f"f400_roots={','.join(map(str, f400_roots)) or 'NONE'}",
        flush=True,
    )
    return (
        core,
        calibrated_original,
        calibrated_mutated,
        core_record,
        flag_record,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output_handle = None
    if args.output is not None:
        output_handle = args.output.open("w", encoding="utf-8", newline="\n")
        sys.stdout = output_handle

    print("FIRST-GAP WALL-STRIPPED NORM-CORE ATLAS")
    print(f"thm2943_dependency_sha256={SOURCE_SHA256}")
    print("scope=first_gap_(0,1,2,M);M=6..26;characteristic=0")
    print(
        "H_M=(n+M)*prod(r=3..floor(M/2),(n+r));"
        "C_M=primitive(gcd(P,A)/(q200^5*c300*B_M*E_M))"
    )
    results = {}
    core_records = []
    flag_records = []
    for width in range(6, 27):
        result = audit_width(width)
        results[width] = result
        core_records.append(result[3])
        flag_records.append(result[4])

    core_digest = sha256("\n".join(core_records).encode()).hexdigest()
    flag_digest = sha256("\n".join(flag_records).encode()).hexdigest()
    require(
        core_digest == EXPECTED_CORE_RECORD_DIGEST,
        "global core record digest changed",
    )
    require(
        flag_digest == EXPECTED_FLAG_RECORD_DIGEST,
        "global flag record digest changed",
    )
    print("core_record_digest=" + core_digest)
    print("flag_record_digest=" + flag_digest)

    baseline_six_step = all(
        (
            expected_core_degree(width + 6)
            + int(width + 6 in (11, 12, 21))
        )
        - (
            expected_core_degree(width)
            + int(width in (11, 12, 21))
        )
        == 124
        for width in range(6, 21)
    )
    require(baseline_six_step, "baseline six-width degree step changed")
    shift_hits_7_13 = tuple(
        shift
        for shift in range(-13, 14)
        if results[13][0]
        .gcd(primitive_positive(results[7][0](X + shift)))
        .degree()
        > 0
    )
    shift_hits_6_12 = tuple(
        shift
        for shift in range(-12, 13)
        if results[12][0]
        .gcd(primitive_positive(results[6][0](X + shift)))
        .degree()
        > 0
    )
    require(
        not shift_hits_7_13 and not shift_hits_6_12,
        "simple six-width shift hostile disappeared",
    )
    print(
        "degree_law=23M-2floor(M/3)-2floor(M/2)-floor(2M/3)-3-epsilon;"
        "epsilon_support=11,12,21;baseline_d(M+6)-d(M)=124"
    )
    print(
        "shift_hostile=C13_gcd_C7(n+s)=1_for_s=-13..13;"
        "C12_gcd_C6(n+s)=1_for_s=-12..12"
    )

    u6, v6 = results[6][1], results[6][2]
    calibrated_resultant = abs(int(u6.resultant(v6)))
    resultant_digest = sha256(str(calibrated_resultant).encode()).hexdigest()
    require(
        calibrated_resultant.bit_length()
        == EXPECTED_M6_CALIBRATED_RESULTANT["bits"]
        and len(str(calibrated_resultant))
        == EXPECTED_M6_CALIBRATED_RESULTANT["digits"]
        and resultant_digest == EXPECTED_M6_CALIBRATED_RESULTANT["digest"],
        "M6 calibrated-chart resultant changed",
    )
    require(
        all(
            calibrated_resultant
            % math.gcd(abs(int(u6(depth))), abs(int(v6(depth))))
            == 0
            for depth in range(51)
        ),
        "M6 sampled resultant divisibility failed",
    )
    print(
        "calibrated_divisor_lemma="
        "gcd(U_M(n),V_M(n))_divides_abs(Res_n(U_M,V_M))"
    )
    print(
        f"M6_calibrated_resultant_bits={calibrated_resultant.bit_length()};"
        f"digits={len(str(calibrated_resultant))};"
        f"digest={resultant_digest};sample_depths=0..50"
    )
    print(
        "scope_boundary=coefficient_positivity_is_archimedean;"
        "no_fixed_prime_bank_for_uncalibrated_minors"
    )
    print("outside_paired_depth_controls=63")
    print("outside_determinant_evaluations=126")
    print(
        "M25_wall_invoice=q17:(1,0,19,4,0,23,28,28,0);"
        "c19:(0,1,19,4,0,23,24,24,0)"
    )
    print("determinant_interpolation_used=YES;finite_field_inference_used=NO")
    print("all_exact_checks=PASS")
    if output_handle is not None:
        output_handle.flush()
        output_handle.close()


if __name__ == "__main__":
    main()
