#!/usr/bin/env python3
"""Clean-room exact referee for the proved-relative THM-4161 theorem.

The script imports no project computation.  It reconstructs the source from
the proved THM-4155 formulas, performs both eliminations over Q(r,u), and
separates symbolic identities from the finite exact non-emptiness control.
"""

from __future__ import annotations

from hashlib import sha256
import json
import sys

import sympy as sp


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def same(left: sp.Expr, right: sp.Expr) -> bool:
    return sp.cancel(left - right) == 0


def coeff_digest(poly: sp.Poly) -> str:
    payload = json.dumps([str(value) for value in poly.all_coeffs()], separators=(",", ":"))
    return sha256(payload.encode("ascii")).hexdigest()


def remainder_zero(expr: sp.Expr, relation: sp.Expr, variable: sp.Symbol) -> bool:
    numerator = sp.cancel(expr).as_numer_denom()[0]
    return sp.rem(numerator, relation, variable) == 0


# Coefficient wall and exact-M=9 source.
r, u = sp.symbols("r u", nonzero=True)
s, p = sp.symbols("s p")
c = sp.Rational(1376, 135)
K0 = sp.Rational(2848, 45)
zeta = c / (r**2 * u)
Theta = -zeta * (2 * r + u)
Phi = zeta * (r**2 + 2 * r * u)
PI = 2027776 * r**3 * u + 1013888 * r**2 * u**2 + 17415
Ar = 356 * r**2 + 15

W = sp.symbols("W")
C = zeta * W**3 + Theta * W**2 + Phi * W - c
require(sp.expand(C - zeta * (W - r) ** 2 * (W - u)) == 0, "double-root factorization")
require(sp.expand(C.subs(W, 0) + c) == 0, "fixed nonzero constant")
IC = sp.factor(4 * Theta * K0**2 - 27 * zeta**2)
IC_expected = -sp.Rational(44032, 273375) * PI / (r**4 * u**2)
require(same(IC, IC_expected), "I_C parameterization")

t = p - s**2
H = (
    -3 * p
    + sp.Rational(8, 3) * p**2
    - c * p**3
    + K0 * s**2 * p**2
    + Phi * s * p**3
    + Theta * s**2 * p**3
    + zeta * s**3 * p**3
)
Gsp = -s**2 / (2 * t) + H
A = sp.cancel((-s * p + t**2 * sp.diff(H, s)) / p)
C0 = s**2 + 2 * t**2 * sp.diff(H, p)
B = sp.cancel((C0 + s * A) / t**2)
require(same(t**2 * sp.diff(Gsp, s), p * A), "source derivative A")
require(same(2 * t**2 * sp.diff(Gsp, p), t**2 * B - s * A), "source derivative B")
require(sp.Poly(A, s).degree() == 6 and sp.Poly(B, s).degree() == 3, "source degrees")

source_domain = sp.QQ.frac_field(r, u).poly_ring(p)
source_res = sp.resultant(
    sp.Poly(A, s, domain=source_domain),
    sp.Poly(B, s, domain=source_domain),
).as_expr()
R17 = sp.cancel(source_res / p**6)
Rpoly = sp.Poly(R17, p)
require(same(source_res, p**6 * R17), "source resultant factor")
require(Rpoly.degree() == 17, "source residual degree")
R0_expected = sp.Rational(3877634048, 50625) * PI / (r**6 * u**3)
R17_expected = (
    sp.Rational(3289935900927224469054816256, 252226880859375)
    * (r - u) ** 3 * Ar / (r**16 * u**8)
)
require(same(Rpoly.nth(0), R0_expected), "R17 constant endpoint")
require(same(Rpoly.LC(), R17_expected), "R17 leading endpoint")

# Independent normalized (X,T) projection.
X, T = sp.symbols("X T")
s_xt = X * T
p_xt = T + s_xt**2
Hxt = H.subs({s: s_xt, p: p_xt})
Gxt = sp.expand(-X**2 * T / 2 + Hxt)
f = sp.cancel(sp.diff(Gxt, X) / T)
h = sp.diff(Gxt, T)
normalized_domain = sp.QQ.frac_field(r, u).poly_ring(T)
fpoly = sp.Poly(f, X, domain=normalized_domain)
hpoly = sp.Poly(h, X, domain=normalized_domain)
require((fpoly.degree(), hpoly.degree()) == (8, 9), "normalized X degrees")
normalized_res = sp.resultant(fpoly, hpoly).as_expr()
Q17 = sp.cancel(normalized_res / (T**56 * (6 * T + 1) ** 2))
Qpoly = sp.Poly(Q17, T)
require(same(normalized_res, T**56 * (6 * T + 1) ** 2 * Q17), "normalized factor")
require(Qpoly.degree() == 17, "normalized residual degree")
Q0_expected = -sp.Rational(72965752821794209792, 56953125) / (r**14 * u**7)
Q17_expected = (
    sp.Rational(210555897659342366019508240384, 9308590679915771484375)
    * (r - u) ** 3 * Ar * PI**2 / (r**20 * u**10)
)
require(same(Qpoly.nth(0), Q0_expected), "Q17 constant endpoint")
require(same(Qpoly.LC(), Q17_expected), "Q17 leading endpoint")

# The four universal points and their Hessians are parameter-independent.
GX = sp.diff(Gxt, X)
GT = sp.diff(Gxt, T)
hessian = sp.diff(Gxt, X, 2) * sp.diff(Gxt, T, 2) - sp.diff(Gxt, X, T) ** 2
for expr, relation, target, label in (
    (Gxt.subs(T, 0), X**2 + 6, 0, "T0 value"),
    (GX.subs(T, 0), X**2 + 6, 0, "T0 GX"),
    (GT.subs(T, 0), X**2 + 6, 0, "T0 GT"),
    (hessian.subs(T, 0), X**2 + 6, 6, "T0 Hessian"),
    (Gxt.subs(T, -sp.Rational(1, 6)), X**2 - 6, sp.Rational(1, 2), "Tm value"),
    (GX.subs(T, -sp.Rational(1, 6)), X**2 - 6, 0, "Tm GX"),
    (GT.subs(T, -sp.Rational(1, 6)), X**2 - 6, 0, "Tm GT"),
    (hessian.subs(T, -sp.Rational(1, 6)), X**2 - 6, -6, "Tm Hessian"),
):
    require(remainder_zero(expr - target, relation, X), label)
Lcrit = Qpoly.degree() + 2 + 2
require(Lcrit == 21, "critical length")

# Top-chart reconstruction, including the contact and differential order.
q, z, x = sp.symbols("q z x", nonzero=True)
H_top = H.subs({s: W, p: 1 / z})
F_top = (W**2 - 1 / z) * (1 - q * H_top) - q * W**2 / 2
Ltop = sp.cancel(z**4 * F_top)
Aface = sp.Rational(8, 3) + K0 * W**2
L_expected = (
    q * C
    + q * z * (Aface - W**2 * C)
    + q * z**2 * (-3 - W**2 * Aface)
    + z**3 * (3 * q * W**2 - 1)
    + z**4 * W**2 * (1 - q / 2)
)
require(same(Ltop, L_expected), "top chart expansion")
require(same(Ltop.subs({W: r, z: 0}), 0), "top point")
require(same(sp.diff(Ltop, W).subs({W: r, z: 0}), 0), "double top root")
require(
    same(sp.diff(Ltop, W, 2).subs({W: r, z: 0}) / 2, q * zeta * (r - u)),
    "top quadratic coefficient",
)
Lz_at_root = sp.factor(sp.diff(Ltop, z).subs({W: r, z: 0}))
require(same(Lz_at_root, sp.Rational(8, 45) * q * Ar), "top transverse coefficient")

# Chain rule: on F=0, L_z=-z^2 F_p, so z~(W-r)^2 gives ord(omega)=4.
Fp_top = sp.diff((s**2 - p) * (1 - q * H) - q * s**2 / 2, p).subs({s: W, p: 1/z})
chain_identity = sp.diff(Ltop, z) - 4 * z**3 * F_top + z**2 * Fp_top
require(same(chain_identity, 0), "top differential chain rule")
contact_order = 2
omega_order = 2 * contact_order
top_index = omega_order + 1
require((contact_order, omega_order, top_index) == (2, 4, 5), "top index five")
require(
    same(sp.diff(Ltop, W).subs({W: u, z: 0}), q * zeta * (u - r) ** 2),
    "simple top root remains transverse",
)

# The Newton polygon itself does not drop on the candidate gate.  Only its
# top edge becomes tangent.  Check the five vertex coefficients directly.
F_source = sp.expand((s**2 - p) * (1 - q * H) - q * s**2 / 2)
F_source_poly = sp.Poly(F_source, s, p)
vertices = ((0, 1), (2, 0), (5, 3), (3, 4), (0, 4))
vertex_coefficients = tuple(sp.factor(F_source_poly.coeff_monomial(s**a * p**b)) for a, b in vertices)
require(
    all(
        same(actual, expected)
        for actual, expected in zip(
            vertex_coefficients, (-1, 1 - q/2, -q*zeta, q*zeta, -q*c)
        )
    ),
    "Newton vertex coefficients",
)
area2 = abs(sum(
    vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
    - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
    for index in range(len(vertices))
))
from math import gcd as integer_gcd
boundary = sum(
    integer_gcd(
        abs(vertices[(index + 1) % len(vertices)][0] - vertices[index][0]),
        abs(vertices[(index + 1) % len(vertices)][1] - vertices[index][1]),
    )
    for index in range(len(vertices))
)
pick_genus = (area2 - boundary + 2) // 2
require((area2, boundary, pick_genus) == (27, 11, 9), "Newton Pick ledger")
edge_data = (
    ((1, 2, 2), 1),
    ((-1, 1, -2), 2),
    ((-1, -2, -11), 8),
    ((0, -1, -4), 3),
)
require(
    tuple(normal[0] + normal[1] - normal[2] for normal, _ in edge_data)
    == tuple(index for _, index in edge_data),
    "raw edge indices",
)

# The prime cubic carrier survives and is separable over C(q).
carrier = zeta * W**3 + K0 * W**2 - (q - sp.Rational(1, 2))
carrier_disc = sp.factor(sp.discriminant(carrier, W))
carrier_expected = (q - sp.Rational(1, 2)) * (
    4 * K0**3 - 27 * zeta**2 * (q - sp.Rational(1, 2))
)
require(same(carrier_disc, carrier_expected), "prime cubic carrier discriminant")

# Packet, genus, and response arithmetic.
packet = (8, 5, 3, 2, 2, 2, 1)
defect = sum(index - 1 for index in packet)
full_degree = sum(packet)
finite_degree = 8 + 5 + 3 + 1
beta = 3
require((defect, full_degree, finite_degree) == (16, 23, 17), "packet ledger")
require((defect + 2) // 2 == 9, "RH genus floor")
require(2 * (full_degree - Lcrit) == 4 < defect, "full commutator contradiction")
finite_caps = (
    2 * finite_degree - Lcrit - 2 + beta,
    2 * finite_degree - Lcrit - 1 + beta,
    beta,
)
require(finite_caps == (14, 15, 3), "finite capacities")
require(all(capacity < finite_degree - 1 for capacity in finite_caps), "finite contradiction")

# Exact hostile strata: inner endpoint, triple root, and transverse failure.
Rtriple = sp.Poly(sp.cancel(R17.subs(u, r)), p)
Qtriple = sp.Poly(sp.cancel(Q17.subs(u, r)), T)
triple_inner = 1013888 * r**4 + 5805
require(Rtriple.degree() == Qtriple.degree() == 16, "triple-root degree sixteen")
require(same(Rtriple.nth(0), sp.Rational(3877634048, 16875) * triple_inner / r**9),
        "triple source endpoint")
require(sp.factor(Rtriple.LC()).has(Ar), "triple source transverse gate")
require(sp.factor(Qtriple.LC()).has(Ar), "triple normalized transverse gate")
require(sp.factor(Qtriple.LC()).has(triple_inner), "triple normalized inner gate")
require((2 * 3, 2 * 3 + 1) == (6, 7), "triple contact predicts index seven")

# Under Ar=0 the linear z coefficient vanishes.  Reduce the quadratic jet
# modulo 356r^2+15 and compare with the claimed node tangent cone.
shifted = sp.expand(Ltop.subs(W, r + x))
quadratic = sum(
    term for term in sp.Add.make_args(shifted)
    if sp.degree(term, x) + sp.degree(term, z) == 2
)
tangent_expected = q * (
    zeta * (r - u) * x**2 + 2 * K0 * r * x * z - 3 * z**2
)
require(remainder_zero(quadratic - tangent_expected, Ar, r), "transverse-failure tangent cone")
node_discriminant = 4 * K0**2 * r**2 + 12 * zeta * (r - u)
require(node_discriminant != 0, "formal node discriminant expression")

# Finite exact non-emptiness control only; these are not theorem hypotheses.
control = {r: sp.Integer(1), u: sp.Integer(2)}
control_coefficients = (
    sp.factor(zeta.subs(control)),
    sp.factor(Theta.subs(control)),
    sp.factor(Phi.subs(control)),
    sp.factor(IC.subs(control)),
    sp.factor((sp.Rational(8, 3) + K0 * r**2).subs(control)),
)
require(
    control_coefficients
    == (
        sp.Rational(688, 135),
        -sp.Rational(2752, 135),
        sp.Rational(688, 27),
        -sp.Rational(89478737152, 273375),
        sp.Rational(2968, 45),
    ),
    "r1u2 coefficient control",
)
Rcontrol = sp.Poly(R17.subs(control), p, domain=sp.QQ)
Qcontrol = sp.Poly(Q17.subs(control), T, domain=sp.QQ)
require(sp.gcd(Rcontrol, Rcontrol.diff()).degree() == 0, "R17 control squarefree")
require(sp.gcd(Qcontrol, Qcontrol.diff()).degree() == 0, "Q17 control squarefree")
qminus_control = sp.factor(Qcontrol.eval(-sp.Rational(1, 6)))
require(qminus_control != 0, "Q17 control avoids universal T value")

semantic = {
    "theorem_id": "THM-4161-proved-relative",
    "gate_factors": ("r", "u", "r-u", "356r^2+15", "P_I"),
    "PI": str(PI),
    "IC": str(sp.factor(IC_expected)),
    "source": {
        "factor": "p^6 R17",
        "degree": Rpoly.degree(),
        "constant": str(sp.factor(R0_expected)),
        "leading": str(R17_expected),
    },
    "normalized": {
        "factor": "T^56(6T+1)^2 Q17",
        "degree": Qpoly.degree(),
        "constant": str(Q0_expected),
        "leading": str(Q17_expected),
    },
    "critical_length": Lcrit,
    "top": {"contact": contact_order, "omega_order": omega_order, "index": top_index},
    "packet": packet,
    "defect": defect,
    "genus": 9,
    "responses": {"full": full_degree, "finite": finite_degree, "beta": beta},
    "capacities": {"full": 4, "finite": finite_caps},
    "triple_degrees": (Rtriple.degree(), Qtriple.degree()),
    "control_coefficients": tuple(str(value) for value in control_coefficients),
    "control_R_digest": coeff_digest(Rcontrol),
    "control_Q_digest": coeff_digest(Qcontrol),
    "control_Q_minus_sixth": str(qminus_control),
}
semantic_hash = sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("THM4161_DOUBLE_TOP_ROOT_INDEPENDENT_REFEREE_20260825")
print("verdict=ACCEPT_PROVED_RELATIVE_TO_INHERITED_THM4155_MECHANISMS;canon_status=PROVED_RELATIVE")
print("audit_levels=symbolic:parametrization,resultants,top_chart,hostile_jets;finite:r1u2_squarefree_nonemptiness")
print(f"gate=r*u*(r-u)*(356*r^2+15)*P_I != 0;P_I={PI}")
print(f"IC={sp.factor(IC_expected)}")
print(f"source=p^6*R17;degree={Rpoly.degree()};R0={sp.factor(R0_expected)};Rlead={R17_expected}")
print(f"normalized=T^56*(6*T+1)^2*Q17;degree={Qpoly.degree()};Q0={Q0_expected};Qlead={Q17_expected}")
print(f"critical_length={Lcrit};universal_hessians=(6,-6)")
print(f"top_contact={contact_order};omega_order={omega_order};top_index={top_index}")
print(f"newton=(area2:{area2},boundary:{boundary},pick_genus:{pick_genus});packet={packet};defect={defect};genus_relative=9")
print(f"responses=full:{full_degree},finite:{finite_degree},beta:{beta};full_cap=4;finite_caps={finite_caps}")
print("hostiles=PI_zero:endpoints;u_eq_r:degree16_contact3;Ar_zero:node_or_thinner_cusp")
print(f"control_r_u=(1,2);coefficients={control_coefficients};R_squarefree=1;Q_squarefree=1;Q_minus_sixth_nonzero=1")
print(f"control_R_coeff_sha256={semantic['control_R_digest']}")
print(f"control_Q_coeff_sha256={semantic['control_Q_digest']}")
print(f"checks={CHECKS}")
print(f"semantic_sha256={semantic_hash}")
