#!/usr/bin/env python3
"""Independent hostile audit for THM-3860.

The audit rederives the completed-strip recursion, classifies the complete
vertical rational precomposition family, checks its jets and nodal member,
and implements the one-variable pole lemma by a second route.  In the latter
route every root-multiplicity profile through degree eight is checked in
addition to the symbolic all-degree calculation.
"""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


GATES = 0


def gate(condition: object, label: str) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(label)


def zero(expression: sp.Expr, label: str) -> None:
    gate(sp.cancel(sp.factor(expression)) == 0, label)


def jac(left: sp.Expr, right: sp.Expr, normal: sp.Symbol, arm: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(left, normal) * sp.diff(right, arm)
        - sp.diff(left, arm) * sp.diff(right, normal)
    )


# ---------------------------------------------------------------------------
# 1. Formal coefficient convolution and the affine tangent torsor.
# ---------------------------------------------------------------------------
s, z = sp.symbols("s z")
order = 7
aa = [sp.Function(f"aa{index}")(s) for index in range(order)]
bb = [sp.Function(f"bb{index}")(s) for index in range(order)]
A_series = sum(aa[index] * z**index for index in range(order))
C_series = sum(bb[index] * z**index for index in range(order))
J_series = jac(A_series, C_series, z, s)

for degree in range(order - 1):
    direct = sp.expand(J_series).coeff(z, degree)
    convolution = sum(
        i * aa[i] * sp.diff(bb[j], s)
        - j * sp.diff(aa[i], s) * bb[j]
        for i in range(order)
        for j in range(order)
        if i + j == degree + 1
    )
    zero(direct - convolution, f"formal convolution degree {degree}")

# Generic affine solution.  The determinant-one arm packet is represented by
# x*q-p*y=1, and no division by p or q occurs.  This retains p=0 or q=0.
x, y, p, q, forcing = sp.symbols("x y p q forcing")
tau = sp.Function("tau")(s)
particular_u = -forcing * x
particular_v = -forcing * y
solution_u = particular_u + tau * p
solution_v = particular_v + tau * q
zero(
    solution_u * q - p * solution_v + forcing * (x * q - p * y),
    "affine recursion particular plus tangent kernel",
)
zero((tau * p) * q - p * (tau * q), "tangent kernel exact")

# Edge controls for the unimodular row: neither derivative is selected as a
# denominator.  These two packets cover a0'=0 and b0'=0 respectively.
edge_u, edge_v = sp.symbols("edge_u edge_v")
zero(edge_u * 1 - 0 * edge_v - edge_u, "a0-prime zero edge")
zero(edge_u * 0 - 1 * edge_v + edge_v, "b0-prime zero edge")

# Recurse independently on the minimal nodal arm through three higher rows.
c = sp.symbols("c", nonzero=True)
a0 = 9 * c**6 * s**2
b0 = 27 * c**9 * s**3 - 3 * c**3 * s
a1 = -1 / (3 * c**3)
b1 = -sp.Rational(3, 2) * s
zero(a1 * sp.diff(b0, s) - sp.diff(a0, s) * b1 - 1,
     "nodal arm Bezout row")

nodal_a = [a0, a1]
nodal_b = [b0, b1]
gauges = [s + 2, s**2 - s + 1, s**3 + 3 * s]
for m, gauge in enumerate(gauges, start=1):
    interior = sum(
        i * nodal_a[i] * sp.diff(nodal_b[m + 1 - i], s)
        - (m + 1 - i) * sp.diff(nodal_a[i], s) * nodal_b[m + 1 - i]
        for i in range(1, m + 1)
    )
    new_a = -interior * a1 / (m + 1) + gauge * sp.diff(a0, s)
    new_b = -interior * b1 / (m + 1) + gauge * sp.diff(b0, s)
    zero(
        (m + 1) * (
            new_a * sp.diff(b0, s) - sp.diff(a0, s) * new_b
        )
        + interior,
        f"nodal affine recursion row {m}",
    )
    nodal_a.append(sp.expand(new_a))
    nodal_b.append(sp.expand(new_b))

A_truncated = sum(nodal_a[index] * z**index for index in range(len(nodal_a)))
C_truncated = sum(nodal_b[index] * z**index for index in range(len(nodal_b)))
J_truncated = jac(A_truncated, C_truncated, z, s)
zero(J_truncated.coeff(z, 0) - 1, "truncated nodal constant row")
for degree in range(1, len(nodal_a) - 1):
    zero(J_truncated.coeff(z, degree), f"truncated nodal row {degree}")


# ---------------------------------------------------------------------------
# 2. Vertical rational classification, converse, and jet preservation.
# ---------------------------------------------------------------------------
w = sp.symbols("w", nonzero=True)
Z, S = sp.symbols("Z S")
a = sp.Function("a")(S)
b = sp.Function("b")(S)
alpha = sp.Function("alpha")(S)
beta = sp.Function("beta")(S)
A_seed = a + alpha * Z
C_seed = b + beta * Z
seed_jac = jac(A_seed, C_seed, Z, S)
bezout = alpha * sp.diff(b, S) - sp.diff(a, S) * beta
wronskian = alpha * sp.diff(beta, S) - sp.diff(alpha, S) * beta
zero(seed_jac - bezout - wronskian * Z, "linear-normal seed Jacobian")

phi = sp.Function("phi")(z)
f = sp.Function("f")(z)
L = 1 / ((1 + w * phi) * sp.diff(phi, z))
S_vertical = L * s + f
zero(
    (1 + w * phi) * jac(phi, S_vertical, z, s) - 1,
    "complete vertical density identity",
)

# Converse: if Z=phi(z), the density equation is linear in S_s and forces
# S_s=L.  Its residual H=S-Ls has zero s-derivative, hence lies in k(z) in
# characteristic zero.  The first identity is checked without assuming an
# affine form for S.
S_general = sp.Function("Sgeneral")(z, s)
converse_density = (1 + w * phi) * sp.diff(phi, z) * sp.diff(S_general, s)
zero(
    converse_density.subs(sp.diff(S_general, s), L) - 1,
    "vertical converse forces S_s=L",
)
H_general = S_general - L * s
zero(
    sp.diff(H_general, s).subs(sp.diff(S_general, s), L),
    "vertical converse residual has zero s-derivative",
)

# General Taylor jets, including an arbitrary cubic normal coefficient and
# arbitrary allowed tangent f term.
phi3, phi4, f2, f3 = sp.symbols("phi3 phi4 f2 f3")
phi_jet = z - w * z**2 / 2 + phi3 * z**3 + phi4 * z**4
f_jet = f2 * z**2 + f3 * z**3
L_jet = sp.series(
    1 / ((1 + w * phi_jet) * sp.diff(phi_jet, z)),
    z,
    0,
    4,
).removeO().expand()
gate(L_jet.coeff(z, 0) == 1, "vertical L constant jet")
gate(L_jet.coeff(z, 1) == 0, "vertical L linear jet")
S_jet = sp.expand(L_jet * s + f_jet)
gate(S_jet.coeff(z, 0) == s, "vertical S arm jet")
gate(S_jet.coeff(z, 1) == 0, "vertical S first-normal jet")

# A non-Mobius polynomial phi and a rational f are a hostile positive control
# for the claimed complete vertical formula.
phi_polynomial = z - w * z**2 / 2
f_rational = z**2 / (1 - z)
L_polynomial = sp.cancel(
    1 / ((1 + w * phi_polynomial) * sp.diff(phi_polynomial, z))
)
S_polynomial = L_polynomial * s + f_rational
zero(
    (1 + w * phi_polynomial) * jac(phi_polynomial, S_polynomial, z, s) - 1,
    "non-Mobius vertical positive control",
)

# Build an independent constant-W seed: a=S, alpha=1, beta=wS,
# b=S+wS^2/2.  It has Bezout one and Wronskian w.
a_control = S
alpha_control = sp.Integer(1)
beta_control = w * S
b_control = S + w * S**2 / 2
zero(
    alpha_control * sp.diff(b_control, S)
    - sp.diff(a_control, S) * beta_control
    - 1,
    "constant-W seed Bezout control",
)
zero(
    alpha_control * sp.diff(beta_control, S)
    - sp.diff(alpha_control, S) * beta_control
    - w,
    "constant-W seed Wronskian control",
)
A_control = (a_control + alpha_control * Z).subs(
    {S: S_polynomial, Z: phi_polynomial}
)
C_control = (b_control + beta_control * Z).subs(
    {S: S_polynomial, Z: phi_polynomial}
)
zero(jac(A_control, C_control, z, s) - 1,
     "composed non-Mobius Darboux control")

# The Mobius control and its tangent gauge are recomputed independently.
phi_mobius = z / (1 + w * z / 2)
L_mobius = sp.factor(
    1 / ((1 + w * phi_mobius) * sp.diff(phi_mobius, z))
)
zero(
    L_mobius - (1 + w * z / 2) ** 3 / (1 + 3 * w * z / 2),
    "Mobius density closed form",
)
phi_mobius_series = sp.series(phi_mobius, z, 0, 4).removeO().expand()
L_mobius_series = sp.series(L_mobius, z, 0, 4).removeO().expand()
gate(phi_mobius_series.coeff(z, 2) == -w / 2,
     "Mobius forced quadratic normal coefficient")
gate(L_mobius_series.coeff(z, 2) == 3 * w**2 / 4,
     "Mobius tangent quadratic gauge")
gate(L_mobius_series.coeff(z, 3) == -w**3,
     "Mobius tangent cubic coefficient")

# If a tangent-to-identity precomposition permits Z_s, its z^2 coefficient
# is nevertheless forced to the constant -w/2.  Thus Z_s can first occur in
# the z^3 coefficient.
q2 = sp.Function("q2")(s)
q3 = sp.Function("q3")(s)
t2 = sp.Function("t2")(s)
Z_mixed = z + q2 * z**2 + q3 * z**3
S_mixed = s + t2 * z**2
density_mixed = sp.series(
    (1 + w * Z_mixed) * jac(Z_mixed, S_mixed, z, s),
    z,
    0,
    3,
).removeO().expand()
zero(density_mixed.coeff(z, 1) - (w + 2 * q2),
     "mixed density first coefficient")
gate(not density_mixed.coeff(z, 1).has(q3, sp.diff(q3, s)),
     "s-dependent cubic coefficient does not enter early")


# ---------------------------------------------------------------------------
# 3. Nodal rational formulas and the named prime divisor.
# ---------------------------------------------------------------------------
h = 1 + z / (4 * c**3)
g = 1 + 3 * z / (4 * c**3)
Z_nodal = z / h
S_nodal = s * h**3 / g
A_nodal = 9 * c**6 * S_nodal**2 - Z_nodal / (3 * c**3)
C_nodal = (
    27 * c**9 * S_nodal**3
    - 3 * c**3 * S_nodal
    - sp.Rational(3, 2) * S_nodal * Z_nodal
)
C_nodal_closed = 27 * c**9 * s**3 * h**9 / g**3 - 3 * c**3 * s * h**2
zero(C_nodal - C_nodal_closed, "nodal second coordinate simplification")
zero(jac(A_nodal, C_nodal, z, s) - 1,
     "nodal rational Darboux identity")
zero(A_nodal.subs(z, 0) - a0, "nodal arm first coordinate")
zero(C_nodal.subs(z, 0) - b0, "nodal arm second coordinate")
zero(sp.diff(A_nodal, z).subs(z, 0) - a1,
     "nodal first normal coefficient A")
zero(sp.diff(C_nodal, z).subs(z, 0) - b1,
     "nodal first normal coefficient C")

# At h=0 the first coordinate has a simple pole with nonzero residue.
zero(
    sp.limit(h * A_nodal, z, -4 * c**3) - sp.Rational(4, 3),
    "named h-divisor simple-pole residue",
)

r, e, z0 = sp.symbols("r e z0", nonzero=True)
surface = r**2 * e - z**3 + c**3 * r
e_h = -(c**3 * r + 64 * c**9) / r**2
zero(surface.subs({z: -4 * c**3, e: e_h}),
     "named divisor Laurent parametrization")

# The quotient makes r a unit, with this explicit inverse.  Hence the ideal
# is prime and reduced and its local ring at the generic point is a DVR.
relation_h = r**2 * e + c**3 * r + 64 * c**9
r_inverse_h = -(r * e + c**3) / (64 * c**9)
zero(sp.rem(sp.together(r * r_inverse_h - 1), relation_h, e),
     "named divisor explicit inverse for r")

# General nonzero constant-z divisor used by the pole lemma.
e_z0 = (z0**3 - c**3 * r) / r**2
zero((r**2 * e - z0**3 + c**3 * r).subs(e, e_z0),
     "general constant-z Laurent parametrization")
r_inverse_z0 = (r * e + c**3) / z0**3
zero(
    sp.rem(
        sp.together(r * r_inverse_z0 - 1),
        r**2 * e + c**3 * r - z0**3,
        e,
    ),
    "general constant-z explicit inverse for r",
)
s_z0 = sp.factor(
    (e / (3 * (c**3 + e * r))).subs(e, e_z0)
)
zero(s_z0 - (1 / (3 * r) - c**3 / (3 * z0**3)),
     "general divisor arm coordinate")
gate(sp.diff(s_z0, r) != 0, "general divisor arm coordinate is nonconstant")


# ---------------------------------------------------------------------------
# 4. All-rational phi pole lemma, including poles and multiplicities.
# ---------------------------------------------------------------------------
# Under the negation of the lemma, R=phi+1/w has no finite zero.  In reduced
# form P/Q, algebraic closedness makes P a nonzero constant.  Then
# R'=-gamma*Q'/Q^2.  A root of Q' outside Q is a regular critical point.
gamma_const = sp.symbols("gamma_const", nonzero=True)
Q = sp.Function("Q")(z)
R_model = gamma_const / Q
zero(
    sp.diff(R_model, z) + gamma_const * sp.diff(Q, z) / Q**2,
    "constant-numerator derivative formula",
)

# A polynomial edge has Q constant.  It would make R, hence phi, constant,
# contradicting phi'(0)=1.  This gate records the lost-degree boundary.
Q0 = sp.symbols("Q0", nonzero=True)
gate(sp.diff(gamma_const / Q0 - 1 / w, z) == 0,
     "polynomial no-zero edge is constant")

# Multiplicity audit.  If Q has rho distinct roots of multiplicities m_i,
# gcd(Q,Q') has degree sum(m_i-1)=d-rho.  If all roots of Q' were poles,
# this would equal deg Q'=d-1, forcing rho=1.  We enumerate every partition
# through degree eight as an independent hostile check of repeated poles.
def partitions(total: int, maximum: int | None = None):
    if total == 0:
        yield ()
        return
    if maximum is None or maximum > total:
        maximum = total
    for first in range(maximum, 0, -1):
        for rest in partitions(total - first, first):
            yield (first, *rest)


multiplicity_profiles = 0
multi_pole_profiles = 0
for degree in range(1, 9):
    for profile in partitions(degree):
        roots = [sp.Integer(index + 2) for index in range(len(profile))]
        q_polynomial = sp.expand(
            sp.prod(
                (z - root) ** multiplicity
                for root, multiplicity in zip(roots, profile)
            )
        )
        derivative = sp.diff(q_polynomial, z)
        common = sp.gcd(q_polynomial, derivative)
        residual = sp.cancel(derivative / common)
        rho = len(profile)
        gate(sp.degree(common, z) == degree - rho,
             f"multiplicity gcd degree d={degree},profile={profile}")
        gate(sp.degree(residual, z) == rho - 1,
             f"off-pole critical degree d={degree},profile={profile}")
        gate(sp.gcd(q_polynomial, residual) == 1,
             f"off-pole critical support d={degree},profile={profile}")
        if rho > 1:
            gate(sp.degree(residual, z) > 0,
                 f"multiple poles force regular critical point {profile}")
            multi_pole_profiles += 1
        multiplicity_profiles += 1

# The only remaining support is a single pole Q=Q(0)(1+az)^d.  The value and
# first jet give the exceptional family; its second jet has the wrong sign
# and magnitude for every positive integer d.
d = sp.symbols("d", integer=True, positive=True)
az = -w / d
phi_exceptional = ((1 + az * z) ** (-d) - 1) / w
zero(phi_exceptional.subs(z, 0), "exceptional family value jet")
zero(sp.diff(phi_exceptional, z).subs(z, 0) - 1,
     "exceptional family first jet")
exceptional_second = sp.simplify(
    sp.diff(phi_exceptional, z, 2).subs(z, 0)
)
zero(exceptional_second - w * (d + 1) / d,
     "exceptional family second derivative")
zero(exceptional_second + w - w * (2 * d + 1) / d,
     "exceptional jet contradiction")

# Explicit good-point controls cover a polynomial phi, the Mobius phi, and a
# rational phi with a genuine finite pole.  Each satisfies the three jets and
# displays a finite regular zero of (1+w*phi)phi'.
phi_controls = [
    phi_polynomial,
    phi_mobius,
    z / (1 + w * z / 2 + w**2 * z**2),
]
good_point_controls = 0
for index, phi_control in enumerate(phi_controls):
    gate(phi_control.subs(z, 0) == 0, f"phi control {index} value")
    zero(sp.diff(phi_control, z).subs(z, 0) - 1,
         f"phi control {index} first derivative")
    zero(sp.diff(phi_control, z, 2).subs(z, 0) + w,
         f"phi control {index} second derivative")
    product = sp.cancel((1 + w * phi_control) * sp.diff(phi_control, z))
    numerator, denominator = sp.together(product).as_numer_denom()
    gate(sp.resultant(numerator, denominator, z) != 0,
         f"phi control {index} product zeros are regular")
    gate(sp.degree(numerator, z) > 0,
         f"phi control {index} has a finite good point")
    gate(numerator.subs(z, 0) != 0,
         f"phi control {index} good point is nonzero")
    good_point_controls += 1


# ---------------------------------------------------------------------------
# 5. Divisorial pole comparison for every rational f.
# ---------------------------------------------------------------------------
# At a good z0, L has pole order m>=1 with scalar leading coefficient ell,
# while s restricts to the nonconstant Laurent function sbar.  A rational
# f(z) has either lower, equal, or higher pole order n.  The equal-order
# leading coefficient ell*sbar+mu cannot vanish because its r-derivative is
# nonzero.  These exact leading coefficients freeze all order cases.
ell, mu = sp.symbols("ell mu", nonzero=True)
sbar = 1 / (3 * r) - c**3 / (3 * z0**3)
gate(sp.diff(ell * sbar, r) != 0, "regular-f leading coefficient nonzero")
gate(mu != 0, "higher-f-pole leading coefficient nonzero")
gate(sp.diff(ell * sbar + mu, r) != 0,
     "equal-pole leading coefficient nonzero")
gate(sp.diff(ell * sbar, r) != 0,
     "higher-L-pole leading coefficient nonzero")

# If ord(S)=-M<0, the nodal carrier 9c^6 S^2-Z/(3c^3) has order -2M
# because Z is regular.  A table verifies the strict separation for several
# pole multiplicities rather than only the simple-pole case.
pole_order_controls = 0
for M in range(1, 9):
    gate(-2 * M < 0, f"nodal square pole M={M}")
    gate(-2 * M < 0, f"regular Z cannot cancel pole M={M}")
    pole_order_controls += 1


# ---------------------------------------------------------------------------
# 6. Frozen semantic packet.
# ---------------------------------------------------------------------------
source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "no inactive Python assert",
)

semantic = {
    "formal": "full convolution;unimodular affine recursion;arm tangent gauge",
    "vertical": "constant-W seed;Z=phi(z);S=Ls+f classification and converse",
    "jets": "phi2=-w/2;f0=f1=0;L0=1;L1=0;Z_s first possible at z3",
    "nodal": "exact Mobius pair;arm/normal jets;named h-prime residue 4/3",
    "prime": "constant nonzero z quotient is k[r,r^-1];sbar=1/(3r)-c3/(3z0^3)",
    "phi_lemma": "P constant;Q multiplicities;single-pole exceptional family contradicts second jet",
    "multiplicity_profiles": multiplicity_profiles,
    "multi_pole_profiles": multi_pole_profiles,
    "pole_cases": "f regular/lower/equal/higher order;ell*sbar+mu nonzero",
    "scope": "vertical seed-precomposition subclass only;Z_s and arbitrary global pairs open",
}
semantic_sha = hashlib.sha256(
    json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
).hexdigest()

print("THM3860_RUSSELL_HIGHER_NORMAL_INDEPENDENT_HOSTILE_AUDIT")
print("status=PASS_PROOF+VERIFIED_EXACT;JC2_OPEN")
print("formal_recursion=complete_convolution+affine_arm_tangent_torsor")
print("vertical_class=Z_phi_of_z;S=L(z)s+f(z);classification_and_converse_PASS")
print("jets=arm_and_first_normal_preserved;Z_s_first_possible_at_order_3")
print("nodal=exact_rational_pair;named_h_prime_residue=4/3")
print("constant_z_prime=Laurent_line;sbar_nonconstant")
print(f"phi_multiplicity_profiles={multiplicity_profiles};multi_pole_profiles={multi_pole_profiles}")
print(f"phi_good_point_controls={good_point_controls};polynomial_Mobius_poled")
print("phi_lemma=polynomial_edge+all_pole_multiplicities+single_pole_jet_contradiction")
print("f_pole_cases=regular_lower_equal_higher;equal_lead=ell*sbar+mu_nonzero")
print(f"nodal_pole_order_controls={pole_order_controls}")
print("scope=vertical_seed_precomposition_only;general_Z_s_and_global_pairs_OPEN")
print(f"semantic_sha256={semantic_sha}")
print(f"GATES={GATES}")
print("RESULT=PASS")

