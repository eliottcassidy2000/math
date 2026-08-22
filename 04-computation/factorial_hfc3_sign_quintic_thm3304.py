"""Independent exact certificate for the homogeneous S3-sign quintic cell."""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import sympy as sp

sys.path.insert(0, str(Path(__file__).parent))
import factorial_hfc3_symmetry_cells_support_thm3304 as sparse


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


A, B, x, y = sp.symbols("A B x y")
z = 1 - x - y
delta = (x - y) * (y - z) * (z - x)
e2 = x * y + y * z + z * x
g = sp.expand(delta * (A + B * e2))

# Route 1: sparse homogeneous factorial/Dirichlet evaluation.
forms, gcd_sparse, infinity = sparse.sign_degree_five_probe()
p2 = sp.Poly(165 * A**2 + 66 * A * B + 7 * B**2, A, B)
p4 = sp.Poly(71060 * A**4 + 53295 * A**3 * B
             + 15675 * A**2 * B**2 + 2123 * A * B**3 + 111 * B**4,
             A, B)
require(forms[2].as_expr() == 1728 * p2.as_expr().subs({A: 1, B: sp.Symbol("t")}),
        "sparse second-moment numerator")
require(forms[4].as_expr() == 37623398400 * p4.as_expr().subs({A: 1, B: sp.Symbol("t")}),
        "sparse fourth-moment numerator")

# Route 2: direct iterated integration on the affine triangle.
direct2 = sp.factor(sp.integrate(sp.integrate(sp.expand(g**2),
                                              (y, 0, 1 - x)), (x, 0, 1)))
direct4 = sp.factor(sp.integrate(sp.integrate(sp.expand(g**4),
                                              (y, 0, 1 - x)), (x, 0, 1)))
require(sp.cancel(direct2 - p2.as_expr() / 277200) == 0,
        "direct second moment")
require(sp.cancel(direct4 - p4.as_expr() / 29875045200) == 0,
        "direct fourth moment")

# Projective exclusion.  A!=0 is t=B/A; A=0 is the infinity chart.
t = sp.symbols("t")
q2 = sp.Poly(7 * t**2 + 66 * t + 165, t, domain=sp.QQ)
q4 = sp.Poly(111 * t**4 + 2123 * t**3 + 15675 * t**2
             + 53295 * t + 71060, t, domain=sp.QQ)
resultant = sp.resultant(q2.as_expr(), q4.as_expr(), t)
require(resultant == 846709600, "primitive resultant")
require(sp.gcd(q2, q4).degree() == 0, "affine charts are disjoint")
require(infinity[2] != 0 and infinity[4] != 0
        and p2.coeff_monomial(B**2) == 7
        and p4.coeff_monomial(B**4) == 111,
        "infinity chart")

# Machinery positive control and the sign-selection rule.
require(sp.gcd(q2, q2).degree() == 2, "gcd positive control")
require(sp.expand(g.xreplace({x: y, y: x}) + g) == 0,
        "transposition changes sign")
cyclic = sp.expand(g.xreplace({x: y, y: z}))
# Simultaneous substitution needs a temporary-symbol route.
u, v, w = sp.symbols("u v w")
Delta3 = (u - v) * (v - w) * (w - u)
E23 = u*v + v*w + w*u
G3 = sp.expand(Delta3 * (A*(u+v+w)**2 + B*E23))
require(sp.expand(G3.subs({u: v, v: w, w: u}, simultaneous=True) - G3) == 0,
        "three-cycle fixes the alternating form")

payload = "\n".join([
    "universe=homogeneous degree-5 S3-sign eigenspace on Delta_2",
    "normal_form=Vandermonde*(A*e1^2+B*e2)",
    "odd_moments=zero by transposition sign",
    "I2=(165*A^2+66*A*B+7*B^2)/277200",
    ("I4=(71060*A^4+53295*A^3*B+15675*A^2*B^2+"
     "2123*A*B^3+111*B^4)/29875045200"),
    "primitive_resultant=846709600",
    "infinity_chart=(7,111)",
    "consequence=no nonzero sign quintic kills both surviving moments m=2,4",
    "routes=sparse factorial-Dirichlet numerator + direct iterated triangle integral",
    "controls=gcd-self positive; transposition sign; three-cycle invariance",
]) + "\n"
print(payload, end="")
print("payload_sha256=" + hashlib.sha256(payload.encode("ascii")).hexdigest())
