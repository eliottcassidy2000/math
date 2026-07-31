#!/usr/bin/env python3
"""Canonical exact continuation of the first-gap norm-core atlas.

Finite-exact scope:

    support (0,1,2,M),  27 <= M <= 31.

The M=31 calculation adds the simple quartic wall (M,r)=(31,25)
classified by THM-2964.  Nothing in this companion is extrapolated to
M=32 or to arbitrary width.
"""

from __future__ import annotations

import argparse
import contextlib
import importlib.util
import io
import sys
from concurrent.futures import ProcessPoolExecutor
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


EXPECTED_CORES = {
    27: (
        556,
        "c00583d9eded605bd99fed73892a358e91224a7d4086a19c0a55f0b394285b7f",
    ),
    28: (
        577,
        "9a6b8b02f59f2d1f31a7bc262863c4461aa8d1d0c89c1e4b6901237a68558591",
    ),
    29: (
        599,
        "98f7a7599769664a781edb7fa1f61bd5b4a71092c3fef345de74c66094a42e18",
    ),
    30: (
        617,
        "ba6fc2329242f5f222295ec261b36d3a57e353206716667d67ee26e38b17c7a1",
    ),
    31: (
        639,
        "c9ffa50924c9d0298fd4995b1343fb200d3da0c38e1fcbe07347d97b8908e81c",
    ),
}
EXPECTED_ROOTS = {
    27: ((), (), ()),
    28: ((), (), ()),
    29: ((), (22,), ()),
    30: ((), (), ()),
    31: ((21,), (), (25,)),
}
EXPECTED_INVOICES = {
    (29, 22): (0, 1, 19, 4, 0, 23, 24, 24, 0),
    (31, 21): (1, 0, 19, 4, 0, 23, 28, 28, 0),
    (31, 25): (0, 0, 20, 4, 0, 24, 24, 24, 0),
}
EXPECTED_CORE_RECORD_DIGEST = (
    "cf5513fec8a47aea1b090dca3fa1d931d2a07704cdd576cd20626b197f409fb3"
)
EXPECTED_FLAG_RECORD_DIGEST = (
    "92a1d2eb0a13263f16884fdd796d4e910e048a74522ecb9f985121b72bceead3"
)


def load_base():
    name = f"thm2969_extension_worker_{id(object())}"
    spec = importlib.util.spec_from_file_location(name, BASE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot import THM-2969 companion")
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


def audit_width(width: int) -> dict[str, object]:
    base = load_base()
    old_correction = base.smith_correction
    old_require = base.require

    def smith_correction(M: int, root: int) -> int:
        return old_correction(M, root) + int((M, root) == (31, 25))

    def require(condition: bool, message: str) -> None:
        if (
            not condition
            and message == "pure coefficient resonance census changed: M=31"
        ):
            return
        old_require(condition, message)

    base.smith_correction = smith_correction
    base.expected_core_degree = corrected_degree
    base.require = require
    base.EXPECTED_CORES.update(EXPECTED_CORES)

    capture = io.StringIO()
    with contextlib.redirect_stdout(capture):
        core, calibrated_original, calibrated_mutated, core_record, flag_record = (
            base.audit_width(width)
        )

    forms = base.thm2943.polynomial_forms(width, (0, 1, 2))
    q200, c300, _curvature, _alternate = base.thm2943.flag_polynomials(
        forms
    )
    f400 = forms[2][(4, 0, 0)]
    roots = tuple(
        tuple(
            root
            for root in range(1, width + 1)
            if polynomial % (base.X + root) == 0
        )
        for polynomial in (q200, c300, f400)
    )
    old_require(roots == EXPECTED_ROOTS[width], f"root census changed: M={width}")
    for polynomial, polynomial_roots in zip((q200, c300, f400), roots):
        for root in polynomial_roots:
            old_require(
                base.valuation(polynomial, base.X + root) == 1,
                f"pure wall ceased to be simple: {(width, root)}",
            )
    old_require(
        calibrated_original.gcd(calibrated_mutated).degree() == 0,
        f"calibrated charts ceased to be coprime: M={width}",
    )

    b_factor = base.smith_factor(width)
    e_factor = base.seam_factor(width)
    h_factor = base.expected_flag(width)
    old_require(
        (b_factor * e_factor) % h_factor == 0,
        f"complementary wall division changed: M={width}",
    )
    w_factor = (b_factor * e_factor) // h_factor
    invoices = {}
    for (invoice_width, root), expected in EXPECTED_INVOICES.items():
        if invoice_width != width:
            continue
        factor = base.X + root
        q_order = base.valuation(q200, factor)
        c_order = base.valuation(c300, factor)
        b_order = base.valuation(b_factor, factor)
        e_order = base.valuation(e_factor, factor)
        h_order = base.valuation(h_factor, factor)
        w_order = base.valuation(w_factor, factor)
        core_order = base.valuation(core, factor)
        raw_order = 5 * q_order + c_order + w_order + core_order + h_order
        removed_order = 5 * q_order + c_order + b_order + e_order
        invoice = (
            q_order,
            c_order,
            b_order,
            e_order,
            h_order,
            w_order,
            raw_order,
            removed_order,
            core_order,
        )
        old_require(invoice == expected, f"wall invoice changed: {(width, root)}")
        invoices[root] = invoice

    lines = tuple(line for line in capture.getvalue().splitlines() if line)
    old_require(len(lines) == 1, f"unexpected worker transcript: M={width}")
    return {
        "width": width,
        "line": lines[0],
        "core_record": core_record,
        "flag_record": flag_record,
        "roots": roots,
        "invoices": invoices,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    parser.add_argument("--workers", type=int, default=5)
    args = parser.parse_args()
    if not 1 <= args.workers <= 5:
        raise RuntimeError("workers must lie in 1..5")

    output_handle = None
    if args.output is not None:
        output_handle = args.output.open("w", encoding="utf-8", newline="\n")
        sys.stdout = output_handle

    print("FIRST-GAP NORM-CORE WIDTH-27-31 CONTINUATION")
    print(f"thm2969_dependency_sha256={BASE_SHA256}")
    print("scope=first_gap_(0,1,2,M);M=27..31;characteristic=0")
    print("epsilon_support=11,12,21,31;new_quartic_wall=(31,25)")

    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(audit_width, range(27, 32)))
    results.sort(key=lambda result: int(result["width"]))
    for result in results:
        print(result["line"])

    core_records = [str(result["core_record"]) for result in results]
    flag_records = [str(result["flag_record"]) for result in results]
    core_digest = sha256("\n".join(core_records).encode()).hexdigest()
    flag_digest = sha256("\n".join(flag_records).encode()).hexdigest()
    if core_digest != EXPECTED_CORE_RECORD_DIGEST:
        raise RuntimeError("extension core record digest changed")
    if flag_digest != EXPECTED_FLAG_RECORD_DIGEST:
        raise RuntimeError("extension flag record digest changed")
    print("core_record_digest=" + core_digest)
    print("flag_record_digest=" + flag_digest)

    if not all(
        corrected_degree(width + 6)
        + int(width + 6 in (11, 12, 21, 31))
        - corrected_degree(width)
        - int(width in (11, 12, 21, 31))
        == 124
        for width in range(21, 26)
    ):
        raise RuntimeError("corrected six-width degree step changed")
    print(
        "degree_law=23M-2floor(M/3)-2floor(M/2)-floor(2M/3)-3-epsilon;"
        "corrected_baseline_d(M+6)-d(M)=124_for_M=21..25"
    )
    print(
        "M29_c22_invoice=(0,1,19,4,0,23,24,24,0);"
        "order=(q,c,B,E,H,W,raw,removed,core)"
    )
    print(
        "M31_q21_invoice=(1,0,19,4,0,23,28,28,0);"
        "M31_f25_invoice=(0,0,20,4,0,24,24,24,0)"
    )
    print("M31_pure_roots=q21,f25;both_simple;corrected_core_degree=639")
    print("widths=5;determinant_interpolants=10;paired_depth_controls=15")
    print("individual_outside_determinant_evaluations=30")
    print("M32_boundary=NOT_COMPUTED;no_arbitrary_width_extrapolation")
    print("all_exact_checks=PASS")

    if output_handle is not None:
        output_handle.flush()
        output_handle.close()


if __name__ == "__main__":
    main()
