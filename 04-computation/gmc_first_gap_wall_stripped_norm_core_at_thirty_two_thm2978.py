#!/usr/bin/env python3
"""Exact first-gap wall-stripped norm core at width 32.

Finite-exact scope:

    support (0,1,2,32).

This companion deliberately proves no width-33 or arbitrary-width statement.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import sys
from hashlib import sha256
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = (
    ROOT
    / "04-computation"
    / "gmc_first_gap_wall_stripped_norm_core_atlas_thm2969.py"
)
BASE_BYTES = BASE_PATH.read_bytes().replace(b"\r\n", b"\n").replace(
    b"\r", b"\n"
)
BASE_SHA256 = sha256(BASE_BYTES).hexdigest()
EXPECTED_BASE_SHA256 = (
    "5be5df0fd058f436593339f78cd6240a47975082d22929c135ba811458753bf5"
)
if BASE_SHA256 != EXPECTED_BASE_SHA256:
    raise RuntimeError("THM-2969 dependency hash changed")


WIDTH = 32
EXPECTED_CORE = (
    660,
    "3dc8b52265fd9b96d92851cc1f69df9f987c681843f61b6893ac68506edd323e",
)
EXPECTED_ROOTS = ((), (), ())
EXPECTED_W_ORDERS = tuple(
    6
    if root == 1
    else 24
    if root == 2
    else 28
    if 3 <= root <= 10
    else 26
    if 11 <= root <= 16
    else 24
    if 17 <= root <= 21
    else 23
    if 22 <= root <= 31
    else 20
    for root in range(1, 33)
)


def load_base():
    name = f"thm2969_width32_worker_{id(object())}"
    spec = importlib.util.spec_from_file_location(name, BASE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import THM-2969 companion")
    base = importlib.util.module_from_spec(spec)
    sys.modules[name] = base
    spec.loader.exec_module(base)
    return base


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    output_handle = None
    if args.output is not None:
        output_handle = args.output.open("w", encoding="utf-8", newline="\n")
        sys.stdout = output_handle

    base = load_base()
    base.EXPECTED_CORES[WIDTH] = EXPECTED_CORE

    capture = io.StringIO()
    with contextlib.redirect_stdout(capture):
        core, calibrated_original, calibrated_mutated, core_record, flag_record = (
            base.audit_width(WIDTH)
        )
    lines = tuple(line for line in capture.getvalue().splitlines() if line)
    require(len(lines) == 1, "unexpected THM-2969 worker transcript")

    coefficients = tuple(int(value) for value in core.coeffs())
    require(core.degree() == EXPECTED_CORE[0], "core degree changed")
    require(base.polynomial_digest(core) == EXPECTED_CORE[1], "core digest changed")
    require(all(value > 0 for value in coefficients), "core positivity changed")
    require(base.is_pf2(coefficients), "core PF2 changed")
    require(
        calibrated_original.gcd(calibrated_mutated).degree() == 0,
        "calibrated chart gcd changed",
    )

    forms = base.thm2943.polynomial_forms(WIDTH, (0, 1, 2))
    q200, c300, _curvature, _alternate = base.thm2943.flag_polynomials(forms)
    f400 = forms[2][(4, 0, 0)]
    roots = tuple(
        tuple(
            root
            for root in range(1, WIDTH + 1)
            if polynomial % (base.X + root) == 0
        )
        for polynomial in (q200, c300, f400)
    )
    require(roots == EXPECTED_ROOTS, "integer pure-wall census changed")
    core_integer_roots = tuple(
        root
        for root in range(1, WIDTH + 1)
        if core % (base.X + root) == 0
    )
    require(not core_integer_roots, "wall-stripped core gained an integer root")
    require(
        all(base.smith_correction(WIDTH, root) == 0 for root in range(1, 33)),
        "unexpected Smith correction at width 32",
    )

    b_factor = base.smith_factor(WIDTH)
    e_factor = base.seam_factor(WIDTH)
    h_factor = base.expected_flag(WIDTH)
    require((b_factor * e_factor) % h_factor == 0, "flag absorption changed")
    w_factor = (b_factor * e_factor) // h_factor
    require(
        (
            b_factor.degree(),
            e_factor.degree(),
            h_factor.degree(),
            w_factor.degree(),
        )
        == (693, 108, 15, 786),
        "wall-factor degree invoice changed",
    )
    require(
        base.expected_core_degree(WIDTH) == EXPECTED_CORE[0],
        "uncorrected degree law changed",
    )
    w_orders = tuple(
        base.valuation(w_factor, base.X + root) for root in range(1, 33)
    )
    require(w_orders == EXPECTED_W_ORDERS, "wall-complement orders changed")
    require(
        base.valuation(w_factor, 2 * base.X + 1) == 6,
        "half-wall complement order changed",
    )
    require(
        (
            46 * WIDTH - 26,
            58 * WIDTH - 36 - ((9 * WIDTH - 5) // 2),
            5 * q200.degree()
            + c300.degree()
            + b_factor.degree()
            + e_factor.degree(),
            core.degree(),
        )
        == (1446, 1679, 1019, 660),
        "resultant/common/divisor/core degree invoice changed",
    )

    print("FIRST-GAP NORM-CORE WIDTH-32 CONTINUATION")
    print(f"thm2969_dependency_sha256={BASE_SHA256}")
    print("scope=first_gap_(0,1,2,32);characteristic=0")
    print(lines[0])
    print(flag_record)
    print(
        "pure_integer_walls=q:NONE,c:NONE,f:NONE;"
        "exceptional_M32_smith_correction=NONE;core_integer_roots=NONE"
    )
    print(
        "degree_law=23M-2floor(M/3)-2floor(M/2)-floor(2M/3)-3;"
        "M32_degree=660"
    )
    print("calibrated_chart_gcd_degree=0;wall_complement_integral=1")
    print(
        "degree_invoice=B693,E108,H15,W786,resultant1446,"
        "raw_gcd1679,removed1019,core660"
    )
    print(
        "W_integer_orders=r1:6,r2:24,r3..10:28,r11..16:26,"
        "r17..21:24,r22..31:23,r32:20;half_wall_order=6"
    )
    print(
        "associated_over_Q=G~q^5*c*R*H~q^5*c*B*E*N;"
        "R~W*N;W=(B*E)/H"
    )
    print(
        "positivity_chain=n>=0=>N(n)>0_and_W(n)>0=>R(n)!=0=>"
        "no_common_ternary_zero=>first_four_factorial_moments_not_all_zero"
    )
    print("determinant_interpolants=2;paired_depth_controls=3")
    print("individual_outside_determinant_evaluations=6")
    print("M33_boundary=NOT_COMPUTED;no_arbitrary_width_or_radial_extrapolation")
    print("all_exact_checks=PASS")

    if output_handle is not None:
        output_handle.flush()
        output_handle.close()


if __name__ == "__main__":
    main()
