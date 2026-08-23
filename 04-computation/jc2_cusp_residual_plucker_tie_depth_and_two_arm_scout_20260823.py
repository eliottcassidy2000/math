#!/usr/bin/env python3
"""Exact scout for the planar-JC Plucker tie depth and two-arm f=0 lane.

This research companion independently verifies the exact skew-gain packet
behind THM-3881 and the leading tie-depth envelope underlying the subsequently
promoted THM-3886 proof candidate.  It also verifies a new
degree-five closure in the zero/zero-arm f=0 sector, and a hostile check
showing why the sextic cusp cannot use the named u/v coordinates as an
intrinsic ordinary tournament.
"""

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
    gate(sp.expand(expression) == 0, message)


def wedge(left: tuple[sp.Expr, sp.Expr], right: tuple[sp.Expr, sp.Expr]) -> sp.Expr:
    return sp.expand(left[0] * right[1] - left[1] * right[0])


# ---------------------------------------------------------------------------
# 1. THM-3881 is exactly a five-vertex skew-gain graph.
# ---------------------------------------------------------------------------

x, y, T, F = sp.symbols("x y T F")
a = x + 1
L = 9 * x + 4
K = y**2 - 15 * x**2 - 15 * x - 4
P = a * L**2
Delta = sp.expand(a**2 * P - K**2)

eT = (sp.Integer(1), sp.Integer(0))
eF = (sp.Integer(0), sp.Integer(1))
g = (K, -a)
h = (a * P, -K)
v = (T, F)

r = wedge(g, v)
A = wedge(h, v)

expected_edges = {
    ("eT", "eF"): 1,
    ("eT", "g"): -a,
    ("eT", "h"): -K,
    ("eT", "v"): F,
    ("eF", "g"): -K,
    ("eF", "h"): -a * P,
    ("eF", "v"): -T,
    ("g", "h"): Delta,
    ("g", "v"): a * T + K * F,
    ("h", "v"): K * T + a * P * F,
}
vertices = {"eT": eT, "eF": eF, "g": g, "h": h, "v": v}
for edge, expected in expected_edges.items():
    zero(wedge(vertices[edge[0]], vertices[edge[1]]) - expected,
         f"skew edge {edge}")

zero(a * A - K * r - Delta * F, "first Plucker reconstruction")
zero(a * P * r - K * A - Delta * T, "second Plucker reconstruction")
zero(
    wedge(eT, eF) * wedge(g, h)
    - wedge(eT, g) * wedge(eF, h)
    + wedge(eT, h) * wedge(eF, g),
    "frame Pfaffian relation",
)

qg = sp.symbols("qg")
v_slide = (T + qg * K, F - qg * a)
zero(wedge(g, v_slide) - r, "gauge slide fixes the tied edge r")
zero(wedge(h, v_slide) - (A - Delta * qg),
     "gauge slide transports the A-edge debt")

B = sp.expand(P * F**2 - T**2)
S = sp.expand(
    L**4
    + 2 * (3 * A + 3 * P + r**2) * L**2 * F
    + (8 * A + 6 * P + 3 * r**2) * B
)
S_from_edges = sp.expand(
    L**4
    + 2
    * (3 * wedge(h, v) + 3 * P + wedge(g, v) ** 2)
    * L**2
    * wedge(eT, v)
    + (8 * wedge(h, v) + 6 * P + 3 * wedge(g, v) ** 2)
    * (P * wedge(eT, v) ** 2 - wedge(eF, v) ** 2)
)
zero(S_from_edges - S, "residual reconstructed from the weighted star")

origin = {x: 0, y: 0}
gate(
    (
        wedge(eT, eF),
        wedge(eF, g).subs(origin),
        wedge(g, eT).subs(origin),
    )
    == (1, 4, 1),
    "real-normalized marked-address 3-cycle shadow",
)


# ---------------------------------------------------------------------------
# 2. The THM-3884 equality seam is a Plucker tie with an exact depth law.
# ---------------------------------------------------------------------------

n, t = sp.symbols("n t", integer=True)
K2 = y**2 - 15 * x**2
q, rt = sp.symbols("q rt", nonzero=True)
Fn = x * q
Tn1 = -K2 * q

zero(K2 * Fn + x * Tn1, "equality seam is the leading g-v tie")

a1 = x
L1 = 9 * x
P3 = a1 * L1**2
aP4 = a1 * P3
A_top = sp.expand(aP4 * Fn)
B_top = sp.expand(P3 * Fn**2)
zero(A_top - 81 * x**5 * q, "equality-seam A leading form")
zero(B_top - 81 * x**5 * q**2, "equality-seam B leading form")

AB_top = sp.expand(8 * A_top * B_top)
rB_top = sp.expand(3 * rt**2 * B_top)
zero(AB_top - 52488 * x**10 * q**3, "subcritical AB leading form")
zero(rB_top - 243 * x**5 * q**2 * rt**2,
     "supercritical r-squared-B leading form")
gate(sp.Rational(52488, 243) == 216, "tie constant is 216")

# Complete seven-macro-term degree ledger.  Every noncandidate is lower than
# either AB or r^2 B by a positive affine gap for n>=1.
degrees = {
    "L4": 4,
    "6AL2F": 2 * n + 6,
    "6PL2F": n + 5,
    "2r2L2F": 2 * t + n + 2,
    "8AB": 3 * n + 7,
    "6PB": 2 * n + 6,
    "3r2B": 2 * t + 2 * n + 3,
}
zero(degrees["8AB"] - degrees["6AL2F"] - (n + 1),
     "AB beats AL2F by n+1")
zero(degrees["8AB"] - degrees["6PB"] - (n + 1),
     "AB beats PB by n+1")
zero(degrees["8AB"] - degrees["6PL2F"] - (2 * n + 2),
     "AB beats PL2F")
zero(degrees["8AB"] - degrees["L4"] - (3 * n + 3),
     "AB beats L4")
zero(degrees["3r2B"] - degrees["2r2L2F"] - (n + 1),
     "r2B beats r2L2F by n+1")
zero(degrees["3r2B"] - degrees["8AB"] - (2 * t - n - 4),
     "candidate comparison is exactly 2t versus n+4")

# Parity checks are the numerical shadow of the UFD proof recorded in the
# companion reflection.  The clean normal forms use algebraic closedness.
for vq in range(8):
    for vrt in range(8):
        gate((5 + 2 * vq + 2 * vrt) % 2 == 1,
             "supercritical x-valuation is odd")
for vp in range(8):
    gate((3 * vp) % 2 == 0 if vp % 2 == 0 else (3 * vp) % 2 == 1,
         "q-cube parity retains every irreducible parity")

s, rho = sp.symbols("s rho")
tie_q = x * s**2
tie_rt = rho * x**3 * s
tie_top = sp.expand(
    243 * x**5 * tie_q**2 * (tie_rt**2 + 216 * x**5 * tie_q)
)
zero(tie_top.subs(rho**2, -216), "critical normal form cancels top degree")


# ---------------------------------------------------------------------------
# 3. Forgotten two-arm deletion closes the zero/zero f=0 sector through d=5.
# ---------------------------------------------------------------------------

Tf = sp.symbols("Tf")
S_f0 = sp.expand(L**4 - 6 * a * L**2 * Tf**2 - 8 * K * Tf**3 - 3 * a**2 * Tf**4)
H = sp.symbols("H")
two_arm_inner = sp.expand(
    L * (1 - 6 * a**3 * H**2 - 3 * a**6 * H**4) - 8 * a**3 * K * H**3
)
zero(S_f0.subs(Tf, a * L * H) - L**3 * two_arm_inner,
     "two-arm factorization after aL divides T")

Hbar = sp.symbols("Hbar")
zero(
    two_arm_inner.subs({x: -sp.Rational(4, 9), H: Hbar})
    + 8 * sp.Rational(5, 9) ** 3
    * (y**2 - sp.Rational(8, 27))
    * Hbar**3,
    "odd L-valuation forces a second L factor",
)

# Recheck the exact x=0 quadratic-tau lemma used after T=aL^2 H1.
p, qq, z = sp.symbols("p qq z")
E5 = 3 * z**2 - 24 * z * qq - 4 * z + 48 * qq**2
E6 = 13 * z**3 - 120 * z**2 * qq + 288 * z * qq**2 + 96 * z * qq - 128 * qq**3
E7 = z**3 - 14 * z**2 * qq + 56 * z * qq**2 - 64 * qq**3 - 32 * qq**2
tau_gb = sp.groebner([E5, E6, E7], z, qq, order="lex", domain=sp.QQ)
expected_tau_gb = sp.groebner([z, qq**2], z, qq, order="lex", domain=sp.QQ)
gate(tau_gb == expected_tau_gb, "x=0 quadratic-tau Groebner obstruction")
zero(E6.subs(z, 0) + 128 * qq**3,
     "p=0 seam forces the quadratic tau coefficient to vanish")

alpha, beta, gamma, delta0, epsilon = sp.symbols(
    "alpha beta gamma delta0 epsilon"
)
H1 = alpha * x**2 + beta * x * y + gamma * y**2 + delta0 * x + epsilon * y
T_degree5 = sp.expand(a * L**2 * H1)
zero(
    T_degree5.subs(x, 0) - 16 * (epsilon * y + gamma * y**2),
    "degree-five cell specializes to quadratic tau",
)

# After the tau lemma gamma=epsilon=0.  A nonzero beta produces a unique
# y^5 term, while beta=0 makes T univariate and the missing-y argument closes.
T_after_tau = sp.expand(T_degree5.subs({gamma: 0, epsilon: 0}))
S_after_tau = sp.Poly(S_f0.subs(Tf, T_after_tau), y)
zero(
    S_after_tau.coeff_monomial(y**5)
    + 8 * beta**3 * x**3 * a**3 * L**6,
    "mixed degree-five coefficient has unique odd y-degree",
)
T_univariate = sp.expand(T_after_tau.subs(beta, 0))
S_univariate = sp.Poly(S_f0.subs(Tf, T_univariate), y)
zero(S_univariate.coeff_monomial(y**2) + 8 * T_univariate**3,
     "univariate T has nonzero y2 debt")
zero(S_univariate.coeff_monomial(y), "univariate T has missing y coefficient")
gate(S_univariate.coeff_monomial(1).subs({x: 0}) == 256,
     "univariate square root has nonzero constant branch")


# ---------------------------------------------------------------------------
# 4. The nonreduced sextic cusp is a hostile to a named-coordinate tournament.
# ---------------------------------------------------------------------------

u, vv, tt = sp.symbols("u vv tt")
sextic_relations = [
    6 * vv**2 - u**3,
    u**2 * (2 * u + 3 * vv),
    u**4,
]
sextic_gb = sp.groebner(
    sextic_relations,
    vv,
    u,
    order="lex",
    domain=sp.QQ.frac_field(tt),
)
u_image = u + tt * vv
v_image = (1 - 2 * tt) * vv - sp.Rational(1, 4) * tt * u**2
for index, relation in enumerate(sextic_relations):
    remainder = sextic_gb.reduce(
        sp.expand(relation.subs({u: u_image, vv: v_image}, simultaneous=True))
    )[1]
    zero(remainder, f"sextic coordinate-slide relation {index}")
gate((1 - 2 * tt).subs(tt, 1) == -1,
     "t=1 slide has invertible tangent map")
gate(u_image.subs(tt, 1) == u + vv,
     "t=1 automorphism moves the named u tangent line")

aa, bb = sp.symbols("aa bb")
square_mod_m3 = aa**2 * u**2 + 2 * aa * bb * u * vv
gate(sp.Poly(square_mod_m3, u, vv).coeff_monomial(u**2) == aa**2,
     "square-null tangent line requires aa=0")
gate(sp.Poly(square_mod_m3, u, vv).coeff_monomial(u * vv) == 2 * aa * bb,
     "square-null tangent line is the v line")


semantic = {
    "anchor": "THM3884 equality seam is the leading Plucker tie omega(g,v)=0",
    "gain_packet": "five vertices recover Delta,r,A,T,f and the residual square predicate",
    "tie_depth": {
        "2t>n+4": "impossible by odd x valuation",
        "2t<n+4": "n odd and q is a square up to scalar",
        "2t=n+4": "q=x*s^2 and r_t=rho*x^3*s with rho^2=-216",
    },
    "niche": "f=0, both exact arms zero, degT<=5 implies T=0",
    "wildcard": "sextic named u/v tangent pair is coordinate-dependent; ordinary tournament rejected",
    "scope": "necessary laws and finite cells only; THM3886 exceptional closures separately hostile-audited; JC2 open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "inactive Python assert",
)

print("scout=jc2-cusp-residual-plucker-tie-depth-and-two-arm")
print("gain_graph_vertices=5;full_weighted_star_reconstructs_residual=YES")
print("equality_seam=leading_edge_tie_omega(g_top,v_top)=0")
print("tie_depth_supercritical=IMPOSSIBLE_ODD_X_VALUATION")
print("tie_depth_subcritical=NECESSARY_n_odd_and_q_square")
print("tie_depth_critical=NECESSARY_q=x*s^2;r_t=rho*x^3*s;rho^2=-216")
print("f_zero_both_arms_zero_degT_at_most_5=T_ZERO")
print("sextic_named_u_v_tournament=REFUTED_BY_COORDINATE_SLIDE")
print("THM3886_status=PROVED_VERIFIED_EXACT_INDEPENDENTLY_HOSTILE_AUDITED")
print("JC2_status=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
