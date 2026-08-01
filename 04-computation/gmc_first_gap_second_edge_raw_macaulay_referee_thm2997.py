#!/usr/bin/env python3
"""Independent raw-chart referee for THM-2997.

This companion does not reconstruct the P_j tables by interpolation.  At each
audited physical width it instead builds the actual factorial Q,C,F forms,
extracts their top four n-slices, forms the selected 36 by 36 Macaulay matrices,
and differentiates the raw determinant through order three.  It then removes
the universal THM-2942 flag factor q200^6*c300*K and compares the resulting
logarithmic coefficients with the frozen P_1,P_2,P_3 table of THM-2997.

The distinction is deliberate: the raw leading selected matrix is invertible
on the audited physical locus but is not the identity.  Identity appears only
after local response/leading-matrix normalization.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

import sympy as sp
from flint import fmpq_mat


ROOT = Path(__file__).resolve().parents[1]
ATLAS_PATH = (
    ROOT
    / "04-computation"
    / "gmc_first_gap_wall_stripped_norm_core_strict_ulc_through_thirty_four_thm2982.py"
)
TARGET_PATH = (
    ROOT
    / "04-computation"
    / "gmc_first_gap_wall_stripped_all_width_second_edge_circuit_positivity_thm2997.py"
)
EXPECTED_ATLAS_SHA256 = (
    "645353f04d2143f91b7b213e7223761c65c72a81c0ebb522184643fdd97ca24b"
)
EXPECTED_TARGET_SHA256 = (
    "40959bd9e47fb9ea7bfd9b0ac98b6a303aa1b4b840047b9eef5b2ee214997d5e"
)
WIDTHS = (6, 7, 8, 9, 10, 11, 12, 20, 33, 34)
EXPECTED_RECORD_DIGEST = (
    "b11188341b9f7387e021e7ae3c08e541b724dc6bd9e91af9061ddf9bbed85015"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(data).hexdigest()


def load_module(name: str, path: Path, expected_hash: str):
    require(lf_sha256(path) == expected_hash, f"dependency hash changed: {path.name}")
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


ATLAS = load_module("thm2997_raw_referee_atlas", ATLAS_PATH, EXPECTED_ATLAS_SHA256)
TARGET = load_module("thm2997_raw_referee_target", TARGET_PATH, EXPECTED_TARGET_SHA256)
BASE = ATLAS.load_base("thm2997_raw_referee")
P1, P2, P3 = TARGET.load_jets()


def as_fraction(value) -> Fraction:
    return Fraction(int(value.numer()), int(value.denom()))


def sympy_fraction(value) -> Fraction:
    value = sp.cancel(value)
    require(bool(value.is_Rational), f"expected rational value, got {value}")
    return Fraction(int(value.p), int(value.q))


def trace(matrix) -> object:
    require(matrix.nrows() == matrix.ncols(), "trace needs a square matrix")
    return sum((matrix[index, index] for index in range(matrix.nrows())), start=0)


def polynomial_log_coefficients(poly) -> tuple[Fraction, Fraction, Fraction]:
    """Return [t],[t^2],[t^3] of log(a0+a1*t+a2*t^2+a3*t^3+...)."""
    degree = poly.degree()
    require(degree >= 3, "polynomial needs four leading slices")
    a0 = Fraction(int(poly[degree]))
    a1 = Fraction(int(poly[degree - 1]))
    a2 = Fraction(int(poly[degree - 2]))
    a3 = Fraction(int(poly[degree - 3]))
    require(a0 != 0, "zero leading coefficient")
    x1, x2, x3 = a1 / a0, a2 / a0, a3 / a0
    return (
        x1,
        x2 - x1 * x1 / 2,
        x3 - x1 * x2 + x1 * x1 * x1 / 3,
    )


def top_four_slices(forms):
    slices = [[], [], [], []]
    for form in forms:
        degree = max(poly.degree() for poly in form.values())
        for drop in range(4):
            slices[drop].append(
                {
                    monomial: int(poly[degree - drop])
                    for monomial, poly in form.items()
                }
            )
    return slices


def selected_chart_log_coefficients(slices):
    matrices = [
        fmpq_mat(ATLAS.selected_rows(BASE, forms_at_drop))
        for forms_at_drop in slices
    ]
    m0, m1, m2, m3 = matrices
    identity = fmpq_mat(m0.nrows(), m0.ncols())
    for index in range(identity.nrows()):
        identity[index, index] = 1
    require(m0 != identity, "raw selected leading matrix unexpectedly identity")
    inverse = m0.inv()
    x1, x2, x3 = inverse * m1, inverse * m2, inverse * m3
    return (
        as_fraction(trace(x1)),
        as_fraction(trace(x2)) - as_fraction(trace(x1 * x1)) / 2,
        as_fraction(trace(x3))
        - as_fraction(trace(x1 * x2))
        + as_fraction(trace(x1 * x1 * x1)) / 3,
    )


def expected_log_coefficients(width: int):
    substitution = {
        TARGET.M: width,
        TARGET.U: 2**width,
        TARGET.V: 3**width,
    }
    d_value = 2 ** (2 * width) + 3 * 2**width - 3 * 3**width - 1
    values = [sympy_fraction(poly.subs(substitution)) for poly in (P1, P2, P3)]
    return (
        values[0] / d_value**2,
        values[1] / (16 * d_value**4),
        -values[2] / (128 * d_value**6),
    )


def audit_width(width: int):
    forms = BASE.thm2943.polynomial_forms(width, (0, 1, 2))
    slices = top_four_slices(forms)
    raw = selected_chart_log_coefficients(slices)
    q200, c300, curvature, _alternate = BASE.thm2943.flag_polynomials(forms)
    q_log = polynomial_log_coefficients(q200)
    c_log = polynomial_log_coefficients(c300)
    k_log = polynomial_log_coefficients(curvature)
    resultant = tuple(
        raw[index] - 6 * q_log[index] - c_log[index] - k_log[index]
        for index in range(3)
    )
    expected = expected_log_coefficients(width)
    require(resultant == expected, f"raw-minus-flag jet mismatch at M={width}")
    if width == 6:
        require(
            resultant[2]
            == Fraction(958351870363086969113, 6204146484375000),
            "M=6 third-log-coefficient anchor changed",
        )
    serialized = json.dumps(
        [[value.numerator, value.denominator] for value in resultant],
        separators=(",", ":"),
    )
    return {
        "M": width,
        "raw_chart_identity": False,
        "flag_leaders_nonzero": all(
            poly[poly.degree()] != 0 for poly in (q200, c300, curvature)
        ),
        "raw_minus_flag_matches_P123": True,
        "jet_digest": sha256(serialized.encode("ascii")).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = [audit_width(width) for width in WIDTHS]
    lines = [json.dumps(row, sort_keys=True, separators=(",", ":")) for row in records]
    record_digest = sha256("\n".join(lines).encode("ascii")).hexdigest()
    require(record_digest == EXPECTED_RECORD_DIGEST, "raw-chart record digest changed")
    transcript = "\n".join(
        [
            "THM-2997 INDEPENDENT RAW MACAULAY REFEREE",
            f"atlas_dependency_sha256={lf_sha256(ATLAS_PATH)}",
            f"target_dependency_sha256={lf_sha256(TARGET_PATH)}",
            "method=actual factorial top-four slices;raw selected 36x36 log determinant;subtract q200^6*c300*K",
            "normalization=raw leading selected matrix invertible and nonidentity;identity only after local inverse/response normalization",
            "audited_widths=" + ",".join(str(width) for width in WIDTHS),
            *lines,
            f"record_digest={record_digest}",
            "all_exact_checks=PASS",
        ]
    ) + "\n"
    if args.output is None:
        print(transcript, end="")
    else:
        args.output.write_text(transcript, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
