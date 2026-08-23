#!/usr/bin/env python3
"""Exact hostile probe for the fixed-nonzero-product marked-root stratum.

This is a synthesis sidecar connecting THM-3851's reciprocal-cubic torus to
the reserved THM-3859 marked-root graph lane and to the coefficient-content
sidecar used by the provisional THM-3853/3855 inverse-discriminant work.
It makes no claim about the full marked-root grammar.
"""

from __future__ import annotations

import ast
import hashlib
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(newline="\n")


GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"gate failed: {label}")


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(expression) == 0, label)


U, V, T, c = sp.symbols("U V T c")
f = T**3 - U * T**2 + V * T - c
D_c = U**2 * V**2 - 4 * V**3 - 4 * c * U**3 + 18 * c * U * V - 27 * c**2

# Independent discriminant calculation by a Sylvester determinant.
sylvester = sp.Matrix(
    [
        [1, -U, V, -c, 0],
        [0, 1, -U, V, -c],
        [3, -2 * U, V, 0, 0],
        [0, 3, -2 * U, V, 0],
        [0, 0, 3, -2 * U, V],
    ]
)
zero(-sylvester.det() - D_c, "fixed-product cubic discriminant")

# For c nonzero, the marked root is a unit and the incidence ring is the
# Laurent graph k[U,T,T^-1].
T_inverse = (T**2 - U * T + V) / c
zero(T * T_inverse - 1 - f / c, "marked root inverse")
V_laurent = c * T**-1 - T**2 + U * T
zero(f.subs(V, V_laurent), "Laurent marked-root graph")
zero(c * T_inverse - T**2 + U * T - V, "quotient and Laurent maps compose")
gate(
    sp.Poly(c - T**3 + U * T**2, T).eval(0) == c,
    "the Laurent pole disappears only at the excluded boundary c=0",
)

# The coefficient map's Jacobian is the monogenic different up to the unit T.
source_jacobian = sp.diff(V_laurent, T)
different = sp.diff(f, T).subs(V, V_laurent)
zero(source_jacobian + different / T, "different and source Jacobian")

# Repeated-root normalization and the triple-root cusp packet.
q = sp.symbols("q", nonzero=True)
U_q = 2 * q + c * q**-2
V_q = q**2 + 2 * c * q**-1
zero(f.subs({U: U_q, V: V_q}) - (T - q) ** 2 * (T - c * q**-2), "repeated-root factorization")
zero(D_c.subs({U: U_q, V: V_q}), "normalization lies on the discriminant")
zero(sp.diff(U_q, q) - 2 * (1 - c * q**-3), "U normalization derivative")
zero(sp.diff(V_q, q) - 2 * q * (1 - c * q**-3), "V normalization derivative")
derivative_gcd = sp.monic(
    sp.gcd(
        sp.together(sp.diff(U_q, q)).as_numer_denom()[0],
        sp.together(sp.diff(V_q, q)).as_numer_denom()[0],
    ),
    q,
)
gate(derivative_gcd == q**3 - c, "exact triple-root cusp addresses")

# Saturate the two coefficient-difference equations by c*q*r*(q-r).
# A unit ideal proves injectivity for c,q,r nonzero and q!=r.
r, inverse = sp.symbols("r inverse")
same_U = 2 * q**2 * r**2 - c * (q + r)
same_V = q * r * (q + r) - 2 * c
injectivity_groebner = sp.groebner(
    [same_U, same_V, inverse * c * q * r * (q - r) - 1],
    inverse,
    c,
    q,
    r,
    order="grevlex",
)
gate(
    len(injectivity_groebner.polys) == 1
    and injectivity_groebner.polys[0].as_expr() == 1,
    "fixed-product branch normalization injective",
)

# The projective normalization has exactly q=0 and q=infinity above the line
# at infinity, and their images are distinct smooth points when c!=0.
s, t, Z = sp.symbols("s t Z")
U_st = t * (2 * s**3 + c * t**3)
V_st = s * (s**3 + 2 * c * t**3)
Z_st = s**2 * t**2
D_h = U**2 * V**2 - 4 * V**3 * Z - 4 * c * U**3 * Z + 18 * c * U * V * Z**2 - 27 * c**2 * Z**4
zero(D_h.subs({U: U_st, V: V_st, Z: Z_st}), "projective normalization")
gate(
    (U_st.subs({s: 0, t: 1}), V_st.subs({s: 0, t: 1}), Z_st.subs({s: 0, t: 1}))
    == (c, 0, 0),
    "q=0 projective address",
)
gate(
    (U_st.subs({s: 1, t: 0}), V_st.subs({s: 1, t: 0}), Z_st.subs({s: 1, t: 0}))
    == (0, 1, 0),
    "q=infinity projective address",
)
gate(sp.diff(D_h, Z).subs({U: 1, V: 0, Z: 0}) == -4 * c, "first infinity point smooth for c nonzero")
gate(sp.diff(D_h, Z).subs({U: 0, V: 1, Z: 0}) == -4, "second infinity point smooth")

# Scaling by a cube root of c identifies every nonzero c with THM-3851's
# c=1 reciprocal model.
gamma, u, v, tau = sp.symbols("gamma u v tau")
scaled_D = D_c.subs({U: gamma * u, V: gamma**2 * v, c: gamma**3})
D_one = u**2 * v**2 - 4 * v**3 - 4 * u**3 + 18 * u * v - 27
zero(scaled_D - gamma**6 * D_one, "cube-root scaling of the discriminant")
zero(
    f.subs({T: gamma * tau, U: gamma * u, V: gamma**2 * v, c: gamma**3})
    - gamma**3 * (tau**3 - u * tau**2 + v * tau - 1),
    "cube-root scaling of the marked-root cubic",
)

# Ordered roots retain the same toric quotient dictionary at fixed product c.
x, y = sp.symbols("x y", nonzero=True)
z = c / (x * y)
U_xy = x + y + z
V_xy = x * y + c / x + c / y
W_xy = (x - y) * (y - z) * (z - x)
zero(f.subs({U: U_xy, V: V_xy}) - (T - x) * (T - y) * (T - z), "ordered-root factorization")
zero(D_c.subs({U: U_xy, V: V_xy}) - W_xy**2, "fixed-product Vandermonde square")
log_jacobian = sp.det(
    sp.Matrix(
        [
            [x * sp.diff(U_xy, x), y * sp.diff(U_xy, y)],
            [x * sp.diff(V_xy, x), y * sp.diff(V_xy, y)],
        ]
    )
)
zero(log_jacobian + W_xy, "fixed-product logarithmic Jacobian")
fixed_cycle = sp.groebner([x - y, x**2 * y - c], y, x, order="lex")
expected_cycle = sp.groebner([y - x, x**3 - c], y, x, order="lex")
gate(fixed_cycle == expected_cycle, "A3 fixed locus is the three triple-root points")

# Exact cross-control from the incoming THM-3853/3855 lane.  The C-oriented
# one-place target has polynomial normalization A1 with degrees (5,4), so it
# cannot be an affine-coordinate reparametrization of the fixed-product Gm
# branch.  This verifies only the elementary target geometry, not either
# provisional inverse-discriminant theorem.
AA, CC, parameter = sp.symbols("AA CC parameter")
delta_zero = AA * (CC + 5 * AA) * (4 * CC + 19 * AA) * (3 * CC - 17 * AA)
D_parameter = -parameter * (5 * parameter + 1) * (17 * parameter - 3) * (19 * parameter + 4)
CC_parameter = -D_parameter
AA_parameter = parameter * CC_parameter
one_place_target = delta_zero + CC**5
zero(
    one_place_target.subs({AA: AA_parameter, CC: CC_parameter}),
    "incoming one-place target parametrization",
)
gate(sp.degree(AA_parameter, parameter) == 5, "one-place target A degree five")
gate(sp.degree(CC_parameter, parameter) == 4, "one-place target C degree four")
gate(
    sp.gcd(D_parameter, sp.diff(D_parameter, parameter)) == 1,
    "one-place target has four distinct glued finite addresses",
)

# Optimization-safety check.
source_text = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
    "no inactive Python asserts",
)

semantic_lines = [
    "FIXED_PRODUCT_MARKED_ROOT_SYNTHESIS PASS",
    "CUBIC f=T^3-UT^2+VT-c with c nonzero",
    "DISCRIMINANT U^2V^2-4V^3-4cU^3+18cUV-27c^2",
    "SOURCE k[U,V,T]/(f)=k[U,T,T^-1];V=cT^-1-T^2+UT",
    "UNIT T is a nonconstant global unit;the pole at T=0 is unavoidable for c nonzero",
    "BRANCH q->(2q+cq^-2,q^2+2cq^-1) is injective Gm",
    "CUSPS normalization ramifies exactly at q^3=c",
    "PLACES projective normalization adds exactly q=0 and q=infinity at distinct smooth points",
    "SCALING c=gamma^3 reduces the packet to THM-3851's reciprocal model",
    "ORDERED_ROOTS product c torus has Vandermonde^2=discriminant and log Jacobian=-Vandermonde",
    "CROSS_CONTROL incoming one-place target has A1 normalization of degrees (5,4),not Gm",
    "CONNECTION this closes only the fixed-nonzero-product stratum of reserved THM-3859",
    "SIDECAR discriminant forgets index-form unit representation and source character units",
    "SCOPE promotable scoped corollary;not the full THM-3859 classification and no JC conclusion",
]
semantic_sha = hashlib.sha256("\n".join(semantic_lines).encode("utf-8")).hexdigest()
transcript = "\n".join(semantic_lines + [f"GATES {GATES}", f"SEMANTIC_SHA256 {semantic_sha}"]) + "\n"

if len(sys.argv) == 1:
    sys.stdout.write(transcript)
elif len(sys.argv) == 3 and sys.argv[1] == "--verify-frozen":
    frozen = Path(sys.argv[2]).read_bytes()
    gate(frozen == transcript.encode("utf-8"), "frozen transcript byte match")
    print("FROZEN_REPLAY PASS")
    print(f"FROZEN_SHA256 {hashlib.sha256(frozen).hexdigest()}")
else:
    raise SystemExit("usage: constant_product_marked_root_probe.py [--verify-frozen PATH]")
