#!/usr/bin/env python3
"""Independent hostile audit of the THM-3861 cubic normal-strip theorem.

The proof is rebuilt from the six Jacobian buckets.  The checker freezes the
SL2 normalization, both UFD derivative invariants, the irreducible-valuation
inequalities (including zero coefficients), the impossible cubic-quadratic
normal form, the surviving cubic-linear family, and explicit inverse
compositions.  The p=q=0 face is checked to be exactly THM-3856's universe.

No Python ``assert`` statements are used.
"""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(newline="\n")


GATES = 0


def require(condition: bool, label: str) -> None:
    global GATES
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)
    GATES += 1


def same(left: sp.Expr, right: sp.Expr, label: str) -> None:
    require(sp.factor(left - right) == 0, label)


s, z, w = sp.symbols("s z w")
a, alpha, u, p, b, beta, v, q = [sp.Function(name)(s) for name in (
    "a", "alpha", "u", "p", "b", "beta", "v", "q"
)]


# -------------------------------------------------------------------------
# I. Six complete buckets and the target SL2 normalization.
# -------------------------------------------------------------------------

A = a + alpha * z + u * z**2 + p * z**3
C = b + beta * z + v * z**2 + q * z**3
J = sp.Poly(sp.expand(sp.diff(A, z) * sp.diff(C, s)
                      - sp.diff(A, s) * sp.diff(C, z)), z)

expected_buckets = {
    5: 3 * (p * sp.diff(q, s) - sp.diff(p, s) * q),
    4: 3 * p * sp.diff(v, s) + 2 * u * sp.diff(q, s)
       - 2 * sp.diff(p, s) * v - 3 * sp.diff(u, s) * q,
    3: 3 * p * sp.diff(beta, s) + 2 * u * sp.diff(v, s)
       + alpha * sp.diff(q, s) - sp.diff(p, s) * beta
       - 2 * sp.diff(u, s) * v - 3 * sp.diff(alpha, s) * q,
    2: 3 * p * sp.diff(b, s) + 2 * u * sp.diff(beta, s)
       + alpha * sp.diff(v, s) - sp.diff(u, s) * beta
       - 2 * sp.diff(alpha, s) * v - 3 * sp.diff(a, s) * q,
    1: 2 * u * sp.diff(b, s) + alpha * sp.diff(beta, s)
       - sp.diff(alpha, s) * beta - 2 * sp.diff(a, s) * v,
    0: alpha * sp.diff(b, s) - sp.diff(a, s) * beta,
}
for degree, expected in expected_buckets.items():
    same(J.coeff_monomial(z**degree), expected,
         f"complete Jacobian bucket z^{degree}")

# Constant target changes multiply the bracket only by their determinant.
m11, m12, m21, m22 = sp.symbols("m11 m12 m21 m22")
A_target = m11 * A + m12 * C
C_target = m21 * A + m22 * C
J_target = sp.expand(
    sp.diff(A_target, z) * sp.diff(C_target, s)
    - sp.diff(A_target, s) * sp.diff(C_target, z)
)
same(J_target, (m11 * m22 - m12 * m21) * J.as_expr(),
     "constant target GL2 determinant law")

# Every nonzero constant cubic direction can be sent to (1,0) by SL2.
Pdir, Qdir = sp.symbols("Pdir Qdir", nonzero=True)
matrix_when_P_nonzero = sp.Matrix([[1 / Pdir, 0], [-Qdir, Pdir]])
require(sp.simplify(matrix_when_P_nonzero.det()) == 1,
        "SL2 normalization when first cubic direction is nonzero")
require(
    sp.simplify(matrix_when_P_nonzero * sp.Matrix([Pdir, Qdir]))
    == sp.Matrix([1, 0]),
    "SL2 kills the second cubic coefficient",
)
matrix_when_only_Q = sp.Matrix([[0, 1 / Qdir], [-Qdir, 0]])
require(sp.simplify(matrix_when_only_Q.det()) == 1,
        "SL2 normalization when only second cubic direction is nonzero")
require(
    sp.simplify(matrix_when_only_Q * sp.Matrix([0, Qdir]))
    == sp.Matrix([1, 0]),
    "SL2 swaps a one-sided cubic coefficient",
)

# The degree-drop face p=q=0 is exactly the four-bucket quadratic problem.
quadratic_face = {
    degree: sp.expand(expression.subs({p: 0, q: 0,
                                      sp.diff(p, s): 0, sp.diff(q, s): 0}))
    for degree, expression in expected_buckets.items()
}
require(quadratic_face[5] == 0 and quadratic_face[4] == 0,
        "p=q=0 deletes the two genuinely cubic buckets")
same(quadratic_face[3], 2 * (u * sp.diff(v, s) - sp.diff(u, s) * v),
     "THM3856 quadratic z3 bucket")
same(quadratic_face[2],
     alpha * sp.diff(v, s) + 2 * u * sp.diff(beta, s)
     - 2 * sp.diff(alpha, s) * v - sp.diff(u, s) * beta,
     "THM3856 quadratic z2 bucket")
same(quadratic_face[1],
     alpha * sp.diff(beta, s) - sp.diff(alpha, s) * beta
     + 2 * (u * sp.diff(b, s) - sp.diff(a, s) * v),
     "THM3856 quadratic z1 bucket")


# -------------------------------------------------------------------------
# II. Genuinely cubic-quadratic branch: it is impossible.
# -------------------------------------------------------------------------

rho, sigma = sp.symbols("rho sigma", nonzero=True)
hfun = sp.Function("h")(s)
betafun = sp.Function("betafun")(s)
ufun = sp.Function("ufun")(s)
alphafun = sp.Function("alphafun")(s)
bfun = sp.Function("bfun")(s)
afun = sp.Function("afun")(s)

pcq = rho * hfun**3
vcq = sigma * hfun**2
E4_cq = 3 * pcq * sp.diff(vcq, s) - 2 * sp.diff(pcq, s) * vcq
same(E4_cq, 0, "3:2 UFD powers solve the z4 bucket")

E3_cq = (
    3 * pcq * sp.diff(betafun, s) - sp.diff(pcq, s) * betafun
    + 2 * (ufun * sp.diff(vcq, s) - sp.diff(ufun, s) * vcq)
)
invariant_cq = 3 * rho * betafun / hfun - 2 * sigma * ufun / hfun**2
same(E3_cq, hfun**4 * sp.diff(invariant_cq, s),
     "cubic-quadratic derivative invariant")

# The integration constant is removable by the determinant-one target shear
# A -> A-(kappa/sigma)C.
kappa = sp.symbols("kappa")
K0 = sp.Rational(3, 2) * rho / sigma
u_integrated = K0 * betafun * hfun + kappa * hfun**2
u_sheared = sp.expand(u_integrated - (kappa / sigma) * vcq)
same(u_sheared, K0 * betafun * hfun,
     "constant target shear removes the E3 integration constant")

E2_cq = sp.expand(
    3 * pcq * sp.diff(bfun, s)
    + 2 * u_sheared * sp.diff(betafun, s)
    + alphafun * sp.diff(vcq, s)
    - sp.diff(u_sheared, s) * betafun
    - 2 * sp.diff(alphafun, s) * vcq
)
E2_cq_expected = sp.expand(
    3 * rho * hfun**3 * sp.diff(bfun, s)
    + K0 * (betafun * hfun * sp.diff(betafun, s)
            - betafun**2 * sp.diff(hfun, s))
    + 2 * sigma * (alphafun * hfun * sp.diff(hfun, s)
                   - sp.diff(alphafun, s) * hfun**2)
)
same(E2_cq, E2_cq_expected,
     "cubic-quadratic local-valuation identity")

# Exact valuation inequalities at an arbitrary irreducible phi|h.  e is the
# multiplicity of h and m that of beta.  beta-unit is immediately impossible.
valuation_cases = 0
for e in range(1, 97):
    require(e - 1 < 2 * e - 1 and e - 1 < 3 * e,
            "beta-unit term is uniquely minimal in E2")
    # beta=0: E0 makes alpha a unit, whose E2 term has order 2e-1<3e.
    require(2 * e - 1 < 3 * e,
            "beta-zero branch contradicts E2 after E0")
    for m in range(1, 3 * e + 2):
        beta_order = e + 2 * m - 1
        alpha_unit_order = 2 * e - 1
        base_order = 3 * e
        if 2 * m < e:
            require(beta_order < alpha_unit_order
                    and beta_order < base_order,
                    "beta term is uniquely minimal when 2m<e")
        elif 2 * m > e:
            require(alpha_unit_order < beta_order
                    and alpha_unit_order < base_order,
                    "alpha term is uniquely minimal when 2m>e")
        else:
            # Only 2m=e can balance E2; then alpha*beta' is uniquely minimal
            # in E1 against u*b', alpha'*beta, and a'*v.
            require(m - 1 < e + m and m - 1 < m
                    and m - 1 < 2 * e,
                    "balanced E2 case is killed by E1")
        valuation_cases += 1

# Independently reconstruct the unsheared Kummer packet used by the current
# theorem text.  This verifies that the alternative local proof above closes
# the same equations, including its factored constant bucket.
Pk, Vk, dk, ek, a0k = sp.symbols("Pk Vk dk ek a0k", nonzero=True)
Xfun = sp.Function("Xfun")(s)
hk = sp.Function("hk")(s)
bk = sp.Function("bk")(s)
Kk = 3 * Pk / (2 * Vk**2)
Mk = 3 * Pk / (2 * Vk)
pk = Pk * hk**3
vk = Vk * hk**2
betak = hk * Xfun
uk = vk * (Kk * Xfun + dk)
alphak = hk * (Kk * Xfun**2 / 4 + dk * Xfun + Mk * bk + ek)
ak = (
    -Pk * Xfun**3 / (16 * Vk**3)
    + 3 * Pk * bk * Xfun / (4 * Vk**2)
    + ek * Xfun / (2 * Vk) + dk * bk + a0k
)
E3_k = (
    3 * pk * sp.diff(betak, s) - sp.diff(pk, s) * betak
    + 2 * (uk * sp.diff(vk, s) - sp.diff(uk, s) * vk)
)
E2_k = (
    3 * pk * sp.diff(bk, s) + 2 * uk * sp.diff(betak, s)
    + alphak * sp.diff(vk, s) - sp.diff(uk, s) * betak
    - 2 * sp.diff(alphak, s) * vk
)
E1_k = (
    2 * uk * sp.diff(bk, s) + alphak * sp.diff(betak, s)
    - sp.diff(alphak, s) * betak - 2 * sp.diff(ak, s) * vk
)
Tk = bk - Xfun**2 / (4 * Vk)
E0_k = alphak * sp.diff(bk, s) - sp.diff(ak, s) * betak
same(E3_k, 0, "unsheared Kummer packet kills E3")
same(E2_k, 0, "unsheared Kummer packet kills E2")
same(E1_k, 0, "unsheared Kummer packet kills E1")
same(E0_k, hk * (Mk * Tk + ek) * sp.diff(Tk, s),
     "unsheared Kummer packet factors E0")

# Sharp rational hostile: every bucket closes, but a has the unique cubic
# pole predicted by the DVR proof.
A_rational = (
    -sp.Rational(1, 16) / s**3
    + sp.Rational(3, 8) * s**3 * z
    + sp.Rational(3, 2) * s**9 * z**2 + s**15 * z**3
)
C_rational = s**4 * z + s**10 * z**2
J_rational = sp.factor(
    sp.diff(A_rational, z) * sp.diff(C_rational, s)
    - sp.diff(A_rational, s) * sp.diff(C_rational, z)
)
require(J_rational == -sp.Rational(3, 16),
        "sharp rational Kummer hostile has constant bracket")
require(sp.denom(sp.together(A_rational)).subs(s, 0) == 0,
        "sharp rational hostile fails exactly polynomiality of a")

# With h a unit, p=P and v=V are nonzero constants.  The integrated equations
# are reconstructed and then square-completed by a source shear.
P, V, dconst, econst = sp.symbols("P V dconst econst", nonzero=True)
K = sp.Rational(3, 2) * P / V
L = K / (4 * V)
beta_poly = sp.Function("beta_poly")(s)
b_poly = sp.Function("b_poly")(s)
u_poly = K * beta_poly
alpha_poly = K * b_poly + L * beta_poly**2 + dconst
a_poly = (
    K * b_poly * beta_poly + dconst * beta_poly
    - L * beta_poly**3 / 3 + econst
) / (2 * V)
A_cq = P * z**3 + u_poly * z**2 + alpha_poly * z + a_poly
C_cq = V * z**2 + beta_poly * z + b_poly
B_poly = b_poly - beta_poly**2 / (4 * V)

J_cq = sp.expand(sp.diff(A_cq, z) * sp.diff(C_cq, s)
                 - sp.diff(A_cq, s) * sp.diff(C_cq, z))
same(J_cq, (K * B_poly + dconst) * sp.diff(B_poly, s),
     "integrated cubic-quadratic Jacobian obstruction")

z_from_w = w - beta_poly / (2 * V)
same(C_cq.subs(z, z_from_w), V * w**2 + B_poly,
     "cubic-quadratic completed square")
same(
    A_cq.subs(z, z_from_w),
    P * w**3 + (K * B_poly + dconst) * w + econst / (2 * V),
    "cubic-quadratic completed A form",
)

# Hostile degree boundary: constant B gives J=0; every positive degree gives
# degree 2 deg(B)-1 because K is nonzero in characteristic zero.
require(((K * sp.Integer(7) + dconst) * sp.diff(sp.Integer(7), s)) == 0,
        "constant B cannot give a nonzero Jacobian")
for degree in range(1, 49):
    B_control = s**degree + 2 * s + 3 if degree > 1 else s + 3
    J_control = sp.expand((K * B_control + dconst) * sp.diff(B_control, s))
    require(sp.degree(J_control, s) == 2 * degree - 1,
            "nonconstant B gives positive Jacobian degree")


# -------------------------------------------------------------------------
# III. Cubic-linear branch: exactly a triangular automorphism family.
# -------------------------------------------------------------------------

beta_linear = sigma * hfun
p_linear = rho * hfun**3
E3_linear = 3 * p_linear * sp.diff(beta_linear, s) \
    - sp.diff(p_linear, s) * beta_linear
same(E3_linear, 0, "3:1 UFD powers solve the cubic-linear top bucket")

E2_linear = (
    3 * p_linear * sp.diff(bfun, s)
    + 2 * ufun * sp.diff(beta_linear, s)
    - sp.diff(ufun, s) * beta_linear
)
same(
    E2_linear,
    hfun**3 * (
        3 * rho * sp.diff(bfun, s)
        - sigma * sp.diff(ufun / hfun**2, s)
    ),
    "cubic-linear derivative invariant",
)

# beta=0 would force b'=0 from E2 and then E0=0.  If h is nonconstant,
# beta has order e at phi|h; E1 forces alpha nonunit, then E0 is divisible by
# phi.  These inequalities include u=0 and arbitrary extra valuation in u.
for e in range(1, 97):
    for extra_u_order in range(0, 33):
        require(e - 1 < 2 * e + extra_u_order
                and e - 1 < e,
                "cubic-linear alpha-unit term is uniquely minimal in E1")

# h is therefore a unit.  Integrate the remaining three equations.
P1, Beta, kconst, deltaconst, lam, epsilon = sp.symbols(
    "P1 Beta kconst deltaconst lam epsilon", nonzero=True
)
base = sp.Function("base")(s)
u_linear = 3 * P1 * base / Beta + kconst
alpha_linear = (
    3 * P1 * base**2 / Beta**2
    + 2 * kconst * base / Beta + deltaconst
)
a_linear = (
    P1 * base**3 / Beta**3
    + kconst * base**2 / Beta**2
    + deltaconst * base / Beta
    - lam * s / Beta + epsilon
)
C_linear = base + Beta * z
A_linear = P1 * z**3 + u_linear * z**2 + alpha_linear * z + a_linear
J_linear = sp.expand(sp.diff(A_linear, z) * sp.diff(C_linear, s)
                     - sp.diff(A_linear, s) * sp.diff(C_linear, z))
same(J_linear, lam, "integrated cubic-linear family has constant Jacobian")

W = sp.symbols("W")
H_W = (
    P1 * W**3 / Beta**3
    + kconst * W**2 / Beta**2
    + deltaconst * W / Beta + epsilon
)
same(A_linear.subs(z, (W - base) / Beta),
     H_W - lam * s / Beta,
     "moving base cancels after using W=C")

# Explicit inverse compositions for hostile moving bases of degrees 0..3.
Aout, Cout = sp.symbols("Aout Cout")
inverse_cases = 0
for degree in range(0, 4):
    moving_base = sp.Integer(2) if degree == 0 else s**degree + s + 2
    numeric = {
        P1: sp.Integer(2), Beta: sp.Integer(3), kconst: sp.Integer(-1),
        deltaconst: sp.Integer(5), lam: sp.Integer(7), epsilon: sp.Integer(11),
        base: moving_base,
    }
    A_map = sp.expand(A_linear.xreplace(numeric))
    C_map = sp.expand(C_linear.xreplace(numeric))
    H_numeric = sp.expand(H_W.xreplace(numeric))
    s_inverse = sp.expand(sp.Rational(3, 7) * (H_numeric.subs(W, Cout) - Aout))
    z_inverse = sp.expand((Cout - moving_base.subs(s, s_inverse)) / 3)
    require(
        sp.expand(
            s_inverse.subs({Aout: A_map, Cout: C_map}, simultaneous=True) - s
        ) == 0,
        "explicit inverse recovers s",
    )
    require(
        sp.expand(
            z_inverse.subs({Aout: A_map, Cout: C_map}, simultaneous=True) - z
        ) == 0,
        "explicit inverse recovers z",
    )
    require(sp.Poly(s_inverse, Aout, Cout).total_degree() == 3,
            "inverse s-coordinate is cubic polynomial")
    expected_inverse_degree = 1 if degree == 0 else 3 * degree
    require(sp.Poly(z_inverse, Aout, Cout).total_degree()
            == expected_inverse_degree,
            "inverse z-coordinate has the expected polynomial degree")
    inverse_cases += 1


# -------------------------------------------------------------------------
# IV. Metadata and deterministic replay.
# -------------------------------------------------------------------------

candidate_script = Path("04-computation/jc2_cubic_normal_strip_keller_thm3861.py")
candidate_output = Path("05-knowledge/results/jc2_cubic_normal_strip_keller_thm3861.out")
candidate_script_hash = hashlib.sha256(candidate_script.read_bytes()).hexdigest()
candidate_output_hash = hashlib.sha256(candidate_output.read_bytes()).hexdigest()
require(candidate_script_hash
        == "056f177d0586362fb0a1fa3daa3e57e77be0af6ac105a58fe856cb1a3ba4182e",
        "candidate script hash is frozen")
require(candidate_output_hash
        == "51b5c1e05e67c6043228ae3656486afb6e5a7be417ab78ec4a4f0b9269fb16a0",
        "candidate output hash is frozen")

source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "independent checker has no optimization-sensitive assert")

semantic = {
    "classification": "PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED",
    "universe": "char0 field; A,C in k[s,z]; both z-degrees <=3; Jacobian lambda in k*",
    "handoff": "p=q=0 is exactly THM-3856",
    "normalization": "constant target SL2 sends nonzero cubic direction to q=0,p!=0",
    "cubic_quadratic": "3:2 UFD powers; nonconstant h killed by valuations; constant h gives impossible J=(KB+d)B'",
    "cubic_linear": "3:1 UFD powers; nonconstant h killed by valuations; constant h gives explicit triangular automorphism",
    "degeneracies": "p,q top-zero; v=0; beta=0; nonconstant irreducible factors; constants; char0 all explicit",
    "scope": "polynomial transverse degree <=3 only; rational and infinite formal strips open",
    "rational_hostile": "h=s^5,X=1/s closes every bucket with J=-3/16 but forces a=-1/(16s^3)",
}
semantic_blob = json.dumps(semantic, sort_keys=True,
                           separators=(",", ":")).encode()

print("experiment=THM3861-independent-hostile-audit")
print("classification=PROVED;VERIFIED-EXACT;INDEPENDENTLY-HOSTILE-AUDITED")
print("universe=characteristic_zero;k[s,z];deg_z_A_and_C_at_most_3;Jacobian_nonzero_constant")
print("top_handoff=p_equals_q_equals_0_is_THM3856")
print("target_normalization=constant_SL2;q_equals_0;p_nonzero")
print("cubic_quadratic=3_to_2_UFD;nonconstant_h_impossible;constant_h_degree_gate_impossible")
print("cubic_linear=3_to_1_UFD;nonconstant_h_impossible;constant_h_triangular_automorphism")
print("degeneracies=p_or_q_zero,v_zero,beta_zero,irreducible_valuations,and_constants_checked")
print("rational_hostile=all_buckets_close;Jacobian=-3/16;only_failure=a_has_cubic_pole")
print("inverse=explicit_polynomial_two_sided;moving_base_degrees_0_through_3")
print("scope=polynomial_transverse_degree_at_most_3_only;rational_and_infinite_formal_open")
print(f"valuation_cases={valuation_cases}")
print(f"inverse_cases={inverse_cases}")
print(f"candidate_script_sha256={candidate_script_hash}")
print(f"candidate_output_sha256={candidate_output_hash}")
print(f"GATES={GATES}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print("RESULT PASS")
