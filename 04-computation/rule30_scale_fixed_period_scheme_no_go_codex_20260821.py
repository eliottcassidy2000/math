#!/usr/bin/env python3
"""Exact linear-operator controls for scale-fixed Rule 30 profile schemes.

THM-3512's projective recurrence linearizes after restoring cyclic amplitudes:
fixed period-n profiles are projectivized eigenspaces of
(T_n A)_j=A_(2j)+A_(2j+1).  This companion freezes the period 3, 5, and 7
spectra and the 2-adic eigenvalue filter used by the proof reflection.  It is
not a finite Rule 30 graph and constructs no physical profile.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import sympy as sp


THEOREM = Path(__file__).resolve().parents[1] / "01-canon/theorems/THM-3512-rule30-van-der-put-haar-cocycle-and-profinite-automaton-boundary.md"
THEOREM_SHA256 = "82c177c9db3c112f0d9bd7e493d8cfd6db53926ba46151b0d9d2d2ef61199d81"


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def operator(period: int) -> sp.Matrix:
    matrix = sp.zeros(period)
    for row in range(period):
        matrix[row, 2 * row % period] += 1
        matrix[row, (2 * row + 1) % period] += 1
    return matrix


def main() -> None:
    require(hashlib.sha256(THEOREM.read_bytes()).hexdigest() == THEOREM_SHA256,
            "THM-3512 theorem hash")
    variable = sp.symbols("lambda")
    expected = {
        3: (variable - 2) * (variable - 1) * (variable + 1),
        5: ((variable - 2) * (variable - 1) * (variable + 1)
            * (variable**2 + 1)),
        7: ((variable - 2) * (variable - 1)**2
            * (variable**2 + variable + 1)**2),
    }
    records = []
    for period in (3, 5, 7):
        matrix = operator(period)
        characteristic = sp.factor(matrix.charpoly(variable).as_expr())
        require(sp.expand(characteristic - expected[period]) == 0,
                (period, "characteristic polynomial"))
        constant = sp.ones(period, 1)
        require(matrix * constant == 2 * constant, (period, "constant eigenline"))
        eig2_dimension = period - (matrix - 2 * sp.eye(period)).rank()
        eig1_dimension = period - (matrix - sp.eye(period)).rank()
        eigminus1_dimension = period - (matrix + sp.eye(period)).rank()
        cyclo3_dimension = period - (matrix**2 + matrix + sp.eye(period)).rank()
        records.append((
            period,
            str(characteristic),
            eig2_dimension,
            eig1_dimension,
            eigminus1_dimension,
            cyclo3_dimension,
        ))

    require(records == [
        (3, "(lambda - 2)*(lambda - 1)*(lambda + 1)", 1, 1, 1, 0),
        (5, "(lambda - 2)*(lambda - 1)*(lambda + 1)*(lambda**2 + 1)",
         1, 1, 1, 0),
        (7, "(lambda - 2)*(lambda - 1)**2*(lambda**2 + lambda + 1)**2",
         1, 2, 0, 4),
    ], "spectral records")

    # Every factor other than lambda-2 is monic with odd constant term, so
    # all of its algebraic roots are 2-adic units.  The only eigenvalue with
    # positive 2-adic valuation is therefore 2, on the constant amplitude
    # line; its projective ratios are g_j=-A_(j+1)/A_j=-1 and have gap 1.
    non_two_factors = (
        variable - 1,
        variable + 1,
        variable**2 + 1,
        variable**2 + variable + 1,
    )
    require(all(int(sp.Poly(item, variable).TC()) % 2 for item in non_two_factors),
            "non-two factors have unit constant term")

    semantic = (THEOREM_SHA256, tuple(records), tuple(map(str, non_two_factors)))
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("RULE30_SCALE_FIXED_PERIOD_SCHEME_NO_GO_20260821")
    print("status=PROVED_LINEARIZATION_CONSEQUENCE+FINITE-EXACT_SPECTRA;no_prize")
    print(f"dependency_sha256={THEOREM_SHA256}")
    print("operator=(T_n*A)_j=A_(2j)+A_(2j+1);odd_indices_mod_n")
    print("spectral_records=" + repr(tuple(records)).replace(" ", ""))
    print("period5_saturated_geometry=projective_points_only")
    print("period7_saturated_geometry=projective_points_and_projective_lines;genus_zero")
    print("physical_filter=only_lambda_2;constant_amplitude;G_j=-1;gap=1")
    print("obstruction=scale_fixed_repeats_gap_1;THM-3512_no_consecutive_gap_1")
    print(f"semantic_sha256={semantic_sha256}")
    print("scope=periods_3_5_7_exact_scale_fixed_profiles_only;no_multiscale_cycle")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
