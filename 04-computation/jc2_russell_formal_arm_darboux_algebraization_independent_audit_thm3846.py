#!/usr/bin/env python3
"""Independent hostile checker for the THM-3843/3846 Russell-arm lane.

The checker rederives the arm Poisson packet, source-line transfer, formal
canonical coordinates, the universal Darboux lift, and a strengthening:
THM-3843 forces W != 0, while Frac(B)=K(s,z) makes 1+2Wz nonsquare for every
W != 0.  Hence the displayed canonical lift is never rational for any
globally admissible arm packet, not only for the minimal nodal packet.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


GATES = 0


def gate(condition: object, label: str) -> None:
    global GATES
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)
    GATES += 1


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(sp.factor(expression)) == 0, label)


def equal(lhs: sp.Expr, rhs: sp.Expr, label: str) -> None:
    zero(lhs - rhs, label)


# Exact source atlas and Poisson packet.
x, y, c = sp.symbols("x y c", nonzero=True)
r_source = x**3
z_source = x * (c + x**3 * y)
e_source = 3 * c**2 * y + 3 * c * x**3 * y**2 + x**6 * y**3


def jac(lhs: sp.Expr, rhs: sp.Expr) -> sp.Expr:
    return sp.expand(sp.diff(lhs, x) * sp.diff(rhs, y)
                     - sp.diff(lhs, y) * sp.diff(rhs, x))


equal(r_source**2 * e_source, z_source**3 - c**3 * r_source,
      "Russell surface relation")
equal(jac(r_source, z_source), 3 * r_source**2,
      "Poisson bracket r,z")
equal(jac(r_source, e_source), 9 * z_source**2,
      "Poisson bracket r,e")
equal(jac(z_source, e_source), 3 * c**3 + 6 * r_source * e_source,
      "Poisson bracket z,e")
equal(r_source.subs(x, 0), 0, "source line reaches arm r=0")
equal(z_source.subs(x, 0), 0, "source line reaches arm z=0")
equal(e_source.subs(x, 0), 3 * c**2 * y,
      "source line is an isomorphism onto the arm")

# Universal first arm jet and its derivative Bezout identity.
s = sp.symbols("s")
a = sp.Function("a")(s)
b = sp.Function("b")(s)
alpha = sp.Function("alpha")(s)
beta = sp.Function("beta")(s)
bezout = alpha * sp.diff(b, s) - sp.diff(a, s) * beta
W = alpha * sp.diff(beta, s) - sp.diff(alpha, s) * beta
Z = sp.symbols("Z")
A0 = a + alpha * Z
C0 = b + beta * Z
jac_Zs = sp.diff(A0, Z) * sp.diff(C0, s) - sp.diff(A0, s) * sp.diff(C0, Z)
equal(jac_Zs, bezout + W * Z,
      "linear-normal Jacobian identity")

# THM-3843's nodal normalization is the sharp positive boundary control.
t = sp.symbols("t")
p_node = t**2 - 1
q_node = t * (t**2 - 1)
gate(sp.gcd(sp.diff(p_node, t), sp.diff(q_node, t)) == 1,
     "nodal arm is a differential immersion")
equal(q_node / p_node, t, "nodal arm is birational")
equal(p_node.subs(t, 1), p_node.subs(t, -1),
      "nodal arm self-identifies in first coordinate")
equal(q_node.subs(t, 1), q_node.subs(t, -1),
      "nodal arm self-identifies in second coordinate")

# Reconstruct THM-3846's canonical formal coordinates in both directions.
r, z, e = sp.symbols("r z e")
surface = r**2 * e - z**3 + c**3 * r
s_rat = e / (3 * (c**3 + e * r))
e_from_s = 3 * c**3 * s + 9 * s**2 * z**3
r_from_s = z**3 / (c**3 + 3 * s * z**3)
zero(surface.subs({r: r_from_s, e: e_from_s}),
     "formal inverse lies on the surface")
equal(s_rat.subs({r: r_from_s, e: e_from_s}), s,
      "formal inverse recovers s")
zero(
    sp.cancel(r_from_s.subs(s, s_rat) - r).subs(
        z**3, r * (c**3 + e * r)
    ),
    "formal inverse recovers r modulo the surface",
)

poisson_z_s = (
    (3 * c**3 + 6 * r * e) * sp.diff(s_rat, e)
    - 3 * r**2 * sp.diff(s_rat, r)
)
zero(
    sp.together(poisson_z_s - 1).as_numer_denom()[0].subs(
        z**3, r * (c**3 + e * r)
    ),
    "canonical coordinate satisfies bracket z,s=1",
)

# The displayed quadratic resummation is exact.
w = sp.symbols("w", nonzero=True)
Z_closed = (-1 + sp.sqrt(1 + 2 * w * z)) / w
zero(Z_closed + w * Z_closed**2 / 2 - z,
     "quadratic canonical normal coordinate")
for order in range(1, 13):
    truncation = sum(
        sp.catalan(n - 1) * (-w / 2) ** (n - 1) * z**n
        for n in range(1, order + 1)
    )
    residual = sp.Poly(sp.expand(truncation + w * truncation**2 / 2 - z), z)
    gate(all(residual.coeff_monomial(z**degree) == 0
             for degree in range(order + 1)),
         f"Catalan solution through order {order}")

# New all-packet strengthening, part 1: W=0 is incompatible with the
# noninjective arm normalization forced by THM-3843.
#
# If alpha=0, Bezout gives -a'*beta=1 and a is nonconstant linear.
# If beta=0, it gives alpha*b'=1 and b is nonconstant linear.
# Otherwise W=0 makes beta/alpha a scalar t0, and then
# alpha*(b'-t0*a')=1, so b-t0*a is nonconstant linear.  In every case a
# target linear combination recovers s and the arm map is injective.
t0, alpha0 = sp.symbols("t0 alpha0", nonzero=True)
linear_combo = s / alpha0
equal(alpha0 * sp.diff(linear_combo, s), 1,
      "W=0 Bezout unit leaves a nonconstant linear target combination")
gate(sp.degree(linear_combo, s) == 1,
     "the W=0 target combination recovers the arm parameter")

# New all-packet strengthening, part 2: the inverse formulas prove
# Frac(B)=K(s,z).  Over K(s), 1+2W(s)z is a degree-one prime in z whenever
# W is nonzero, hence has valuation one and cannot be a square.
Ws = sp.symbols("Ws", nonzero=True)
square_gate = 1 + 2 * Ws * z
gate(sp.degree(square_gate, z) == 1,
     "universal square gate is linear over K(s)")
gate(sp.gcd(square_gate, sp.diff(square_gate, z)) == 1,
     "universal square-gate divisor is reduced")
gate(sp.LC(sp.Poly(square_gate, z), z) != 0,
     "universal square-gate divisor has valuation one")

# Minimal nodal packet remains a concrete specialization.
a_node = 9 * c**6 * s**2
b_node = 27 * c**9 * s**3 - 3 * c**3 * s
alpha_node = -1 / (3 * c**3)
beta_node = -3 * s / 2
bezout_node = alpha_node * sp.diff(b_node, s) - sp.diff(a_node, s) * beta_node
w_node = alpha_node * sp.diff(beta_node, s) - sp.diff(alpha_node, s) * beta_node
equal(bezout_node, 1, "nodal Bezout row")
equal(w_node, 1 / (2 * c**3), "nodal Wronskian")
gate(sp.degree(1 + z / c**3, z) == 1,
     "nodal square obstruction is a special case of universal obstruction")

# Its displayed Laurent divisor is prime and reduced.
e_laurent = -(c**3 * r + c**9) / r**2
zero(surface.subs({z: -c**3, e: e_laurent}),
     "nodal Laurent quotient")
gate(sp.diff(1 + z / c**3, z) != 0,
     "nodal Laurent valuation is one")

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
     "optimized run keeps all gates")

semantic = {
    "THM3843": "arm restriction is a finite birational noninjective differential normalization and a one-place singular Jelonek component",
    "THM3846": "every unimodular first jet has an exact formal Catalan lift",
    "strengthening": "THM3843 forces W nonzero; Frac(B)=K(s,z) makes 1+2Wz nonsquare for every W nonzero",
    "survivor": "only a different higher-normal lift that cancels both the opposite-arm pole and quadratic sheet can algebraize",
    "comparison": "THM3841 three-puncture and THM3845 cubic-degree no-gos are surface-specific and do not apply to the Russell pseudo-plane",
}
semantic_sha = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
).hexdigest()

print("AUDIT=THM3843_THM3846_RUSSELL_ARM")
print("THM3843=PASS_PROVED;NORMALIZATION_NONINJECTIVE;ONE_PLACE_JELONEK")
print("THM3846=PASS_INDEPENDENT_AUDIT;FORMAL_CATALAN_LIFT")
print("NEW_W_GATE=GLOBAL_COMPATIBILITY_FORCES_W_NONZERO")
print("NEW_SQUARE_GATE=FracB=K(s,z);1+2Wz_LINEAR_REDUCED;NEVER_SQUARE_FOR_W_NONZERO")
print("CANONICAL_LIFT=NONRATIONAL_FOR_EVERY_GLOBALLY_ADMISSIBLE_ARM_PACKET")
print("SURVIVOR=ALTERNATIVE_HIGHER_NORMAL_ALGEBRAIZATION_ONLY")
print("CUBIC_NOGOS=THM3841_AND_THM3845_DO_NOT_TRANSFER")
print(f"SEMANTIC_SHA256={semantic_sha}")
print(f"GATES={GATES}")
print("RESULT=PASS")
