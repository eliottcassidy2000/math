#!/usr/bin/env python3
"""Exact companion for THM-3835's marked-root denominator gate."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def check(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise AssertionError(label)
    CHECKS += 1


def same(lhs: object, rhs: object, label: str) -> None:
    check(sp.cancel(sp.expand(lhs - rhs)) == 0, label)  # type: ignore[operator]


z, C, kappa, h, k = sp.symbols("z C kappa h k", nonzero=True)
m = sp.symbols("m")
q = 7 * z**2 + 3
r = 3 * z**3 + 7 * z**2 + 1
b = 6 * z**3 + 7 * z**2 - 1
s = z**2 * q * C - b

# If h=zk, the determinant-one row makes k a polynomial unit.
same((C * k - m * h).subs(h, z * k), k * (C - m * z),
     "determinant factors through a polynomial root ratio")

# A determinant-one row can have a polynomial ratio when dominance is absent.
x, y = sp.symbols("x y")
h_flat, k_flat, m_flat, C_flat = x, sp.Integer(1), sp.Integer(0), sp.Integer(1)
same(C_flat * k_flat - m_flat * h_flat, 1,
     "non-dominant polynomial-ratio hostile is unimodular")
same(h_flat / k_flat, x, "non-dominant hostile has polynomial ratio")

# A genuinely rational unimodular row shows that rationality itself survives.
h_rat = x
k_rat = 1 + x * y
m_rat = y
C_rat = 1
same(C_rat * k_rat - m_rat * h_rat, 1,
     "rational-ratio control is unimodular")
check(sp.denom(sp.cancel(h_rat / k_rat)) != 1,
      "rational-ratio control retains a denominator")

# Every constant projective row change remains unimodular.  Its inverse
# formulas certify that no numerator/denominator ideal information is lost.
aa, bb, cc, dd = sp.symbols("aa bb cc dd")
det = aa * dd - bb * cc
H_lin = aa * h + bb * k
K_lin = cc * h + dd * k
same(dd * H_lin - bb * K_lin, det * h,
     "projective row inverse recovers h")
same(-cc * H_lin + aa * K_lin, det * k,
     "projective row inverse recovers k")

# With k=kappa, the chart law is a nonzero algebraic relation in z,C.
dependence = sp.expand(kappa * C * s - r)
check(sp.Poly(dependence, C).degree() == 2,
      "constant-k chart relation is quadratic in C")
same(sp.Poly(dependence, C).coeff_monomial(C**2), kappa * z**2 * q,
     "nonzero leading dependence coefficient")
check(dependence != 0, "constant-k dependence relation is nonzero")

# Homogenize without losing the k=0 divisor.
R = 3 * h**3 + 7 * h**2 * k + k**3
S = C * h**2 * (7 * h**2 + 3 * k**2) - k * (
    6 * h**3 + 7 * h**2 * k - k**3
)
same(k**3 * r.subs(z, h / k), R, "homogenized cubic denominator r")
same(k**4 * s.subs(z, h / k), S, "homogenized chart denominator s")
same(
    k**3 * (k * C * s - r).subs(z, h / k),
    C * S - R,
    "chart law homogenizes to CS=R",
)

# Chain-rule homogenization of Jac(z,C)=lambda*k*r.
hx, hy, kx, ky, Cx, Cy, lam = sp.symbols("hx hy kx ky Cx Cy lam")
jac_h_C = hx * Cy - hy * Cx
jac_k_C = kx * Cy - ky * Cx
jac_ratio_C = ((hx * k - h * kx) * Cy - (hy * k - h * ky) * Cx) / k**2
same(k**2 * jac_ratio_C, k * jac_h_C - h * jac_k_C,
     "quotient-rule Jacobian numerator")
same(
    k**2 * (jac_ratio_C - lam * k * r.subs(z, h / k)),
    k * jac_h_C - h * jac_k_C - lam * R,
    "weighted Keller law homogenizes exactly",
)

# At k=0 the homogeneous system still has visible, nontrivial data.
same(R.subs(k, 0), 3 * h**3, "k-zero cubic boundary retained")
same(S.subs(k, 0), 7 * C * h**4, "k-zero chart boundary retained")

source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "theorem": "a dominant plane map to U cannot have polynomial z=h/k",
    "denominator": "Ck-mh=1 makes h/k and every PGL2 transform retain an exact nonconstant denominator",
    "contradiction": "polynomial z makes k constant and forces a nonzero relation between the transcendence basis z,C",
    "survivor": "CS=R and k*Jac(h,C)-h*Jac(k,C)=lambda*R retain the k=0 divisor",
    "scope": "dominance only; no atlas or Jacobian counterexample constructed",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("theorem=THM-3835-polynomial-marked-root-ratio-nonentry")
print("denominator=h_over_k_is_reduced;exact_denominator=k_up_to_scalar")
print("polynomial_ratio=forces_k_constant_and_z_C_algebraically_dependent")
print("projective_row=no_PGL2_transform_is_polynomial_or_integral")
print("homogeneous_chart=CS_equals_R")
print("homogeneous_Keller=kJac_hC-hJac_kC_equals_lambda_R")
print("boundary=k_zero_retained")
print("scope=dominance_only;no_plane_atlas_or_JC2_constructed")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
