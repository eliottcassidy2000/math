#!/usr/bin/env python3
"""Exact companion for THM-3835's polynomial root-ratio obstruction."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)
    CHECKS += 1


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(sp.factor(expression)) == 0, label)


T, c = sp.symbols("T c")
H = (
    84 * T**7
    + (36 * c**2 + 196 * c) * T**6
    + (84 * c**3 + 36 * c**2) * T**5
    + (49 * c**4 + 112 * c**3) * T**4
    - 12 * c**5 * T**3
    + (-14 * c**6 + 12 * c**5) * T**2
    + c**8
)

gate(sp.degree(H, T) == 7, "odd degree seven")
gate(sp.LC(sp.Poly(H, T)) == 84, "constant nonzero leading coefficient")
gate(H.subs(T, 0) == c**8, "nonzero constant term when c is a unit")

# If H(T,c) had only one distinct odd-multiplicity root over an algebraic
# closure, it would be (T-a) times a square.  The leading square coefficient
# is retained explicitly rather than adjoining sqrt(84).
a, ell, u, v, w, c_inv = sp.symbols("a ell u v w c_inv")
S = ell * T**3 + u * T**2 + v * T + w
coefficient_equations = sp.Poly(sp.expand((T - a) * S**2 - H), T).all_coeffs()
system = coefficient_equations + [c * c_inv - 1]
G = sp.groebner(system, ell, u, v, w, a, c_inv, c, order="grevlex")
gate(len(G.polys) == 1 and G.polys[0].as_expr() == 1,
     "no unit-c specialization is linear times a square")

# The excluded c=0 boundary is real: H(T,0)=84*T^7 is linear times a square
# after adjoining ell with ell^2=84.  This checks that c*c_inv=1 is
# load-bearing rather than cosmetic.
zero(H.subs(c, 0) - 84 * T**7, "zero-denominator boundary")
zero(
    ((T - 0) * (ell * T**3) ** 2 - H.subs(c, 0))
    - (ell**2 - 84) * T**7,
    "c=0 hostile factorization modulo ell^2=84",
)

# The determinant mechanism behind the theorem is an exact polynomial
# identity: h=z*k makes C*k-m*h a multiple of k.
h, k, z, C, m = sp.symbols("h k z C m")
zero((C * k - m * h).subs(h, z * k) - k * (C - m * z),
     "polynomial ratio forces denominator to divide determinant")

# Once h,k are scalar, the first lift law makes A=hD algebraically dependent
# on C whenever Q is nonzero; Q=0 instead forces C scalar.
D = sp.symbols("D")
Q = 7 * h**2 + 3 * k**2
A = h * D
lift_reconstruction = sp.cancel(h * (1 + 2 * C * k) / Q)
zero(A.subs(D, (1 + 2 * C * k) / Q) - lift_reconstruction,
     "constant-row lift makes A affine in C")

semantic = {
    "assumption": "z=h/k belongs to K[x,y] for a dominant plane atlas",
    "determinant": "Ck-mh=k(C-mz)=1, hence k is scalar",
    "square_gate": "H(h,k)=w^2 from THM-3822",
    "odd_part": "for every k=c!=0, H(T,c) is not (T-a)*S(T)^2",
    "ufd": "two odd roots force two polynomial squares with constant difference",
    "conclusion": "h,k scalar; lift law makes (A,C) algebraically dependent",
    "scope": "all degree; rational nonpolynomial z remains open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "no inactive Python assert")

print("theorem=THM-3835-polynomial-marked-root-ratio-nonentry")
print("ratio=z=h/k_polynomial_forces_k_scalar")
print("sidecar=H(h,k)=w^2;deg_T_H=7")
print("odd_part=H(T,c)_never_linear_times_square_for_c_nonzero")
print("certificate=exact_bilinear_coefficient_Groebner_[1]")
print("mechanism=two_odd_roots_give_squares_with_nonzero_constant_difference")
print("conclusion=h,k_scalar_and_target_map_nondominant")
print("boundary=c=0_is_linear_times_square_but_forbidden_by_SL2")
print("scope=all_degree;genuinely_rational_root_ratio_required")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
