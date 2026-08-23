#!/usr/bin/env python3
"""Exact algebraic companion for THM-3823's compact-Euler obstruction."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


q, A, C = sp.symbols("q A C")
R = (q - 3) * (q + 1) * (q + 2)
den = 3 * q**2 + 7
Aq = -2 * q**2 * R / den**2
Cq = -q * R / den

Delta = (
    A * (C + 5 * A) * (4 * C + 19 * A) * (3 * C - 17 * A)
    + C**2 * (162 * A**3 + 126 * A**2 * C - 4 * C**3)
    - 27 * A**2 * C**4
)
zero(Delta.subs({A: Aq, C: Cq}), "branch normalization lies on Delta")

# The affine normalization is P1 minus the two denominator roots and
# infinity.  No numerator cancels either finite puncture, and both A,C have
# a pole at infinity.
A_num, A_den = map(sp.Poly, sp.cancel(Aq).as_numer_denom())
C_num, C_den = map(sp.Poly, sp.cancel(Cq).as_numer_denom())
gate(sp.gcd(A_num, sp.Poly(den, q)).degree() == 0,
     "A has no cancellation at finite normalization punctures")
gate(sp.gcd(C_num, sp.Poly(den, q)).degree() == 0,
     "C has no cancellation at finite normalization punctures")
gate(A_num.degree() - A_den.degree() == 1,
     "A has a pole at normalization infinity")
gate(C_num.degree() - C_den.degree() == 2,
     "C has a pole at normalization infinity")
gate(sp.discriminant(den, q) == -84,
     "two distinct finite normalization punctures")

# Exactly the four distinct values q=0,-2,-1,3 map to the square-zero
# origin.  THM-3811 supplies uniqueness away from this fibre via the unique
# double (or triple) root.
origin_preimage = sp.Poly(q * R, q)
gate(origin_preimage.degree() == 4, "origin has four normalization preimages")
gate(sp.gcd(origin_preimage, origin_preimage.diff()).degree() == 0,
     "origin normalization preimages are distinct")
gate(tuple(sorted(sp.solve(origin_preimage.as_expr(), q))) == (-2, -1, 0, 3),
     "origin preimage values")

# The two triple-root parameters q^2=7/3 are distinct, affine, and neither
# lies over the origin.  These are the other two omitted target values.
triple = sp.Poly(3 * q**2 - 7, q)
gate(sp.gcd(triple, sp.Poly(den, q)).degree() == 0,
     "triple parameters avoid normalization punctures")
gate(sp.gcd(triple, origin_preimage).degree() == 0,
     "triple parameters avoid the origin fibre")
gate(sp.discriminant(triple.as_expr(), q) == 84,
     "two distinct triple-root parameters")

# Compact-support Euler ledger.  The affine normalization has chi=-1;
# identifying four points to one lowers chi by three.  Off the branch there
# are three sheets; on the branch minus the three omitted values there is
# exactly the simple companion sheet.
chi_normalization = 2 - 3
chi_branch = chi_normalization - 4 + 1
chi_off_branch = 1 - chi_branch
chi_simple_stratum = chi_branch - 3
chi_U = 3 * chi_off_branch + chi_simple_stratum
gate((chi_normalization, chi_branch, chi_off_branch, chi_simple_stratum, chi_U)
     == (-1, -4, 5, -7, 8), "compact Euler ledger")

# A degree-d plane atlas has exact signed Euler sheet debt 8d-1.  It cannot
# be finite etale, since then compact Euler multiplicativity would say
# 1=8d for a positive integer d.
d = sp.symbols("d", integer=True, positive=True)
zero(8 * d - 1 - (d * chi_U - 1), "Euler sheet-debt formula")
gate(all(1 != 8 * value for value in range(1, 20)),
     "finite-cover contradiction positive-degree controls")

semantic = {
    "normalization": "P1 minus roots(3q^2+7) and infinity; chi_c=-1",
    "singularity": "q=-2,-1,0,3 are identified at the origin; chi_c(branch)=-4",
    "omitted": "origin plus two nonzero triple-root values",
    "strata": "off branch: 3 sheets over chi=5; branch-open: 1 sheet over chi=-7",
    "surface": "chi_c(U)=3*5-7=8",
    "finite": "no finite etale A2_C->U_C",
    "atlas_debt": "generic degree d forces signed Euler debt 8d-1",
    "later": "THM-3841 excludes the remaining nonfinite polynomial atlas",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "inactive Python assert")

print("theorem=THM-3823-nonlinear-cubic-compact-euler-finite-atlas-obstruction")
print("normalization=P1_minus_3_points;chi_c=-1")
print("origin_preimage=q=-2,-1,0,3;four_to_one_correction=-3")
print("branch=chi_c=-4")
print("projection_strata=off_branch:chi5_times3;simple_branch:chi-7_times1")
print("surface=chi_c(U_C)=8")
print("finite_atlas=no_finite_etale_A2_C_to_U_C")
print("quasifinite_atlas=generic_degree_d_has_signed_Euler_debt_8d-1")
print("later=nonfinite_plane_atlas_excluded_by_THM3841")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
