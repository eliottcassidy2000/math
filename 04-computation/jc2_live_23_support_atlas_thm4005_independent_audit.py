#!/usr/bin/env python3
"""No-import hostile audit of the live reduced 2:3 invariant atlas.

This script is intentionally independent of the primary THM-4005 verifier.
It starts from the canonical THM-3992 boundary and the raw THM-3997 normal
row formulas, imposes gamma=-a^3/2, and reconstructs the source-weight table
by grouping actual x^j t^n monomials.  It also audits the limited triangular
target-gauge transfer and independently derives the finite companion packet.
"""

from __future__ import annotations

from hashlib import sha256
import json
import sys

import sympy as sp


sys.stdout.reconfigure(newline="\n")
CHECKS = 0
EXPECTED_SEMANTIC_SHA256 = "432395a4b9ae2e17325e93294f04fcb156fd756d1e86bd63c6fe533860be7fb7"


def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(label)


def equal(left, right, label):
    require(simp(left - right) == 0, label)


x, t, u = sp.symbols("x t u")
a, A5 = sp.symbols("a A5", nonzero=True)
b, d = sp.symbols("b d")

# ----------------------------------------------------------------------
# I. Reconstruct the invariant source-normal rows from canonical formulas.
# ----------------------------------------------------------------------

gamma = -a**3 / 2
beta = b * gamma
delta = d * gamma
theta0 = beta / (2 * gamma)

X = gamma * x
A_boundary_raw = X**2 + a
C_boundary_raw = X**3 + sp.Rational(3, 2) * a * X

N_raw = -sp.Rational(2, 3) / (gamma * a) + 2 * gamma**2 * x * (
    theta0 + a * x / gamma**2
)
K_raw = -x / a + 3 * gamma * (gamma**2 * x**2 + a / 2) * (
    theta0 + a * x / gamma**2
)

# THM-3997 equations (23), (29), and (31), before division by characters.
m0 = beta**2 / 4 + sp.Rational(1, 3) / gamma**3 - sp.Rational(2, 9) / (
    a**3 * gamma**2
)
m1 = gamma * delta + sp.Rational(5, 4) * a * beta / gamma
m2 = 2 * (-3 * a**3 - 5 * gamma) / (5 * a * gamma**2)
m3 = 3 * gamma**2 * theta0

l0 = 3 * a * m1 / (4 * gamma) - theta0 * (3 * a**3 + 2 * gamma) / (
    2 * a * gamma
)
l1 = 3 * a * m2 / (4 * gamma) + (
    -9 * a**6
    + 36 * a**3 * theta0**2 * gamma**6
    - 6 * a**3 * gamma
    - 4 * gamma**2
) / (12 * a**3 * gamma**3)
l2 = 3 * a * theta0 * gamma + 3 * a * m3 / (4 * gamma) + 3 * gamma * m1 / 2
l3 = 3 * a**2 / (2 * gamma) + 3 * gamma * m2 / 2
l4 = 3 * gamma * m3 / 2

M_raw = m0 + m1 * x + m2 * x**2 + m3 * x**3
L_raw = l0 + l1 * x + l2 * x**2 + l3 * x**3 + l4 * x**4

A_rows_expected = (
    1 + A5 * x**2 / 4,
    sp.Rational(4, 3) / A5 + A5 * b * x / 4 + 2 * x**2,
    (9 * A5**3 * b**2 - 512) / (144 * A5**2)
    + (A5 * d + 5 * b) * x / 4
    - 4 * x**2 / (5 * A5)
    + sp.Rational(3, 8) * A5 * b * x**3,
)
C_rows_expected = (
    -sp.Rational(3, 4) * x - A5 * x**3 / 8,
    -sp.Rational(3, 8) * b
    - 4 * x / A5
    - sp.Rational(3, 16) * A5 * b * x**2
    - sp.Rational(3, 2) * x**3,
    -(3 * A5 * d + 7 * b) / (8 * A5)
    + (2816 - 45 * A5**3 * b**2) * x / (480 * A5**2)
    - sp.Rational(3, 16) * (A5 * d + 12 * b) * x**2
    - 12 * x**3 / (5 * A5)
    - sp.Rational(9, 32) * A5 * b * x**4,
)

A_rows_raw = (A_boundary_raw / a, N_raw / a, M_raw / a)
C_rows_raw = (C_boundary_raw / a**4, K_raw / a**4, L_raw / a**4)
for row, (raw, expected) in enumerate(zip(A_rows_raw, A_rows_expected)):
    equal(raw, expected.subs(A5, a**5), f"raw-to-invariant A row {row}")
for row, (raw, expected) in enumerate(zip(C_rows_raw, C_rows_expected)):
    equal(raw, expected.subs(A5, a**5), f"raw-to-invariant C row {row}")

# Residual mu_5 bookkeeping.  The live seam is stable because zeta^6=zeta.
zeta = sp.symbols("zeta", nonzero=True)
require(
    sp.rem(sp.expand((zeta**2 * a) ** 5 - a**5), zeta**5 - 1, zeta) == 0,
    "A5 invariant under residual mu5",
)
require(
    sp.rem(sp.expand(-(zeta**2 * a) ** 3 / 2 - zeta * gamma),
           zeta**5 - 1, zeta) == 0,
    "live seam stable under residual mu5",
)


def extract_table(rows):
    """Return {(t-order,x-degree): (weight, coefficient)}, dropping scalars."""
    table = {}
    for n, row in enumerate(rows):
        poly = sp.Poly(sp.expand(row), x)
        for (j,), coefficient in poly.terms():
            if n == 0 and j == 0:  # target scalar, not retained support
                continue
            table[(n, j)] = (j - 2 * n, simp(coefficient))
    return table


def group_components(table):
    """Represent each weight component as a polynomial in u=x^2*t."""
    components = {}
    for (n, _j), (weight, coefficient) in table.items():
        components[weight] = sp.expand(components.get(weight, 0) + coefficient * u**n)
    return {weight: simp(value) for weight, value in components.items()}


A_table = extract_table(A_rows_expected)
C_table = extract_table(C_rows_expected)
A_components = group_components(A_table)
C_components = group_components(C_table)
require((0, 0) not in A_table, "target scalar is deleted from retained A support")
require(A_table[(1, 0)][0] == A_table[(2, 2)][0] == -2,
        "two A monomials occupy one weight component")
require(C_table[(0, 1)][0] == C_table[(1, 3)][0] == 1,
        "two C monomials occupy one weight component")
require(len(A_table) > len(A_components) and len(C_table) > len(C_components),
        "monomial counting strictly overcounts retained weight support")

# The map (n,j) -> (weight=j-2n,n) is injective.  Future diagonals can add
# higher u powers to a component but cannot cancel a forced lower u coefficient.
tags = {(j - 2 * n, n) for n in range(8) for j in range(2 * n + 5)}
require(len(tags) == sum(2 * n + 5 for n in range(8)),
        "source monomial to weight/u-order tags are injective")

expected_A_entries = {
    (0, 2): (2, A5 / 4),
    (1, 0): (-2, sp.Rational(4, 3) / A5),
    (1, 1): (-1, A5 * b / 4),
    (1, 2): (0, sp.Integer(2)),
    (2, 0): (-4, (9 * A5**3 * b**2 - 512) / (144 * A5**2)),
    (2, 1): (-3, (A5 * d + 5 * b) / 4),
    (2, 2): (-2, -sp.Rational(4, 5) / A5),
    (2, 3): (-1, sp.Rational(3, 8) * A5 * b),
}
expected_C_entries = {
    (0, 1): (1, -sp.Rational(3, 4)),
    (0, 3): (3, -A5 / 8),
    (1, 0): (-2, -sp.Rational(3, 8) * b),
    (1, 1): (-1, -4 / A5),
    (1, 2): (0, -sp.Rational(3, 16) * A5 * b),
    (1, 3): (1, -sp.Rational(3, 2)),
    (2, 0): (-4, -(3 * A5 * d + 7 * b) / (8 * A5)),
    (2, 1): (-3, (2816 - 45 * A5**3 * b**2) / (480 * A5**2)),
    (2, 2): (-2, -sp.Rational(3, 16) * (A5 * d + 12 * b)),
    (2, 3): (-1, -sp.Rational(12, 5) / A5),
    (2, 4): (0, -sp.Rational(9, 32) * A5 * b),
}
require(set(A_table) == set(expected_A_entries), "complete A table keys")
require(set(C_table) == set(expected_C_entries), "complete C table keys")
for key, (weight, coefficient) in expected_A_entries.items():
    require(A_table[key][0] == weight, f"A weight at {key}")
    equal(A_table[key][1], coefficient, f"A coefficient at {key}")
for key, (weight, coefficient) in expected_C_entries.items():
    require(C_table[key][0] == weight, f"C weight at {key}")
    equal(C_table[key][1], coefficient, f"C coefficient at {key}")

# Base support and exact b,d strata.  These sets record forced weights only.
A_base = {2, 0, -2}
C_base = {3, 1, -1}
require(all(simp(A_components[w]) != 0 for w in A_base), "A base components nonzero")
require(all(simp(C_components[w]) != 0 for w in C_base), "C base components nonzero")

A_b_nonzero = A_base | {-1}
C_b_nonzero = C_base | {-2, 0}
equal(A_table[(1, 1)][1] / b, A5 / 4, "b!=0 forces A weight -1")
equal(C_table[(1, 0)][1] / b, -sp.Rational(3, 8), "b!=0 forces C weight -2")
equal(C_table[(1, 2)][1] / b, -sp.Rational(3, 16) * A5,
      "b!=0 forces C weight 0")

A_b0_d_nonzero = A_base | {-4, -3}
C_b0_d_nonzero = C_base | {-4, -3, -2}
equal(A_table[(2, 0)][1].subs(b, 0), -sp.Rational(32, 9) / A5**2,
      "b=0 forces A weight -4")
equal(A_table[(2, 1)][1].subs(b, 0) / d, A5 / 4,
      "b=0,d!=0 forces A weight -3")
equal(C_table[(2, 0)][1].subs(b, 0) / d, -sp.Rational(3, 8),
      "b=0,d!=0 forces C weight -4")
equal(C_table[(2, 1)][1].subs(b, 0), sp.Rational(88, 15) / A5**2,
      "b=0 forces C weight -3")
equal(C_table[(2, 2)][1].subs(b, 0) / d, -sp.Rational(3, 16) * A5,
      "b=0,d!=0 forces C weight -2")

A_b0_d0 = A_base | {-4}
C_b0_d0 = C_base | {-3}
require(all(weight % 2 for weight in C_b0_d0),
        "minimal C displayed weights are all odd")


def minimum_after_even_tail(forced):
    """THM-3987 supplies one even weight <=-4; compute worst collision."""
    collidable = {weight for weight in forced if weight <= -4 and weight % 2 == 0}
    return len(forced) if collidable else len(forced) + 1


strata = {
    "b_nonzero": (
        minimum_after_even_tail(A_b_nonzero),
        minimum_after_even_tail(C_b_nonzero),
    ),
    "b_zero_d_nonzero": (
        minimum_after_even_tail(A_b0_d_nonzero),
        minimum_after_even_tail(C_b0_d_nonzero),
    ),
    "b_zero_d_zero": (
        minimum_after_even_tail(A_b0_d0),
        minimum_after_even_tail(C_b0_d0),
    ),
}
require(strata == {
    "b_nonzero": (5, 6),
    "b_zero_d_nonzero": (5, 6),
    "b_zero_d_zero": (4, 5),
}, "exact support lower-bound strata")

# Sharp support-table hostile: at b=d=0 and a tail of weight -4 the invoice
# itself has size 4x5.  This is not an existence witness for a Darboux pair.
require(len(A_b0_d0 | {-4}) == 4, "formal A lower invoice can remain four")
require(len(C_b0_d0 | {-4}) == 5, "formal C lower invoice first reaches five")

# Special second-row cancellations catch any proof which illegally treats
# every displayed coefficient as generically nonzero in the b!=0 stratum.
b_A_cancel = 16 * sp.sqrt(2) / 3
equal(A_table[(2, 0)][1].subs({A5: 1, b: b_A_cancel}), 0,
      "hostile A weight -4 coefficient can cancel when b!=0")
equal(A_table[(2, 1)][1].subs({A5: 1, b: b_A_cancel, d: -5 * b_A_cancel}), 0,
      "hostile A weight -3 coefficient can cancel when b!=0")
b_C_cancel = sp.sqrt(sp.Rational(2816, 45))
equal(C_table[(2, 1)][1].subs({A5: 1, b: b_C_cancel}), 0,
      "hostile C weight -3 coefficient can cancel when b!=0")
equal(C_table[(2, 0)][1].subs({A5: 1, b: b_C_cancel, d: -7 * b_C_cancel / 3}), 0,
      "hostile C weight -4 coefficient can cancel when b!=0")

# ----------------------------------------------------------------------
# II. Exact gauge/shear transfer and its failure boundary.
# ----------------------------------------------------------------------

kappa, scale = sp.symbols("kappa scale", nonzero=True)
A_minimal_coeffs = {weight: sp.Integer(1) for weight in A_b0_d0}
C_minimal_coeffs = {weight: sp.Integer(1) for weight in C_b0_d0 | {-4}}

# Undoing the recorded triangular operation has C_old=scale*C_new+kappa*A.
C_old = {
    weight: simp(scale * C_minimal_coeffs.get(weight, 0)
                 + kappa * A_minimal_coeffs.get(weight, 0))
    for weight in set(A_minimal_coeffs) | set(C_minimal_coeffs)
}
for odd_weight in C_b0_d0:
    equal(C_old[odd_weight], scale, f"linear shear preserves forced C weight {odd_weight}")
require(len(C_b0_d0) == 4, "four odd C weights survive every recorded shear")

# Hence A-support 3 is impossible outright; if A-support is exactly four,
# b=d=0 and those four odd C pieces forbid C-support 3.
require(min(value[0] for value in strata.values()) == 4,
        "A has at least four weights in fixed gauge")
require(min(value[1] for value in strata.values()) == 5,
        "C has at least five weights in fixed gauge")

# Canceling the fifth (even-tail) weight of an *exact* normalized 4x5 packet
# introduces the other three A weights.  Thus that packet becomes 4x7, not
# 4x4.  This is a useful guard against an overaggressive cancellation story.
C_old_cancel = {
    weight: simp(C_minimal_coeffs.get(weight, 0) - A_minimal_coeffs.get(weight, 0))
    for weight in set(A_minimal_coeffs) | set(C_minimal_coeffs)
}
require(C_old_cancel[-4] == 0, "hostile linear shear can cancel the even C tail")
require(len({w for w, coefficient in C_old_cancel.items() if coefficient != 0}) == 7,
        "exact normalized 4x5 invoice shears to 4x7 in this control")

# The lower-bound atlas permits additional normalized C weights.  If C_new
# contains a copy of all four A components, undoing the shear can cancel that
# entire copy and leave only the four forced odd C pieces.  Therefore no
# pre-gauge 4x5 floor follows from the atlas; the robust transfer is precisely
# the exclusion of C-support three when A-support is four.
C_expanded_coeffs = dict(C_minimal_coeffs)
C_expanded_coeffs.update(A_minimal_coeffs)
C_expanded_old = {
    weight: simp(C_expanded_coeffs.get(weight, 0) - A_minimal_coeffs.get(weight, 0))
    for weight in set(A_minimal_coeffs) | set(C_expanded_coeffs)
}
require({w for w, coefficient in C_expanded_old.items() if coefficient != 0}
        == C_b0_d0, "expanded atlas-compatible packet can shear to pre-gauge 4x4")

# An arbitrary nonlinear triangular target automorphism preserves the source
# Jacobian but can create many support weights; it is outside the transfer.
A_host = x**2 + t
C_host = x**3 + x * t


def jac(first, second):
    return sp.expand(sp.diff(first, x) * sp.diff(second, t)
                     - sp.diff(first, t) * sp.diff(second, x))


equal(jac(A_host, C_host + A_host**2), jac(A_host, C_host),
      "nonlinear triangular target automorphism preserves Jacobian")
host_A_weights = {2, -2}
host_square_weights = {left + right for left in host_A_weights for right in host_A_weights}
require(host_square_weights == {4, 0, -4},
        "nonlinear target operation changes source support")

# ----------------------------------------------------------------------
# III. Independent finite companion derivation and report hostiles.
# ----------------------------------------------------------------------

z = 1 + x**2 * t
p = z * t
y = x * t * p
c20 = -sp.Rational(16, 3) / A5**2
c30 = b**2 / 4 + sp.Rational(2752, 135) / A5**3
R_known = c20 * p**2 + b * y + d * p * y + c30 * p**3
G_known = x**2 * t + 6 * p / A5 + R_known
Q_known = simp(G_known / t)
Q_rows = tuple(simp(sp.expand(Q_known).coeff(t, order)) for order in range(3))
Q_rows_expected = (
    x**2 + 6 / A5,
    b * x + 6 * x**2 / A5 - sp.Rational(16, 3) / A5**2,
    b**2 / 4 + b * x**3 + d * x
    - sp.Rational(32, 3) * x**2 / A5**2
    + sp.Rational(2752, 135) / A5**3,
)
for order in range(3):
    equal(Q_rows[order], Q_rows_expected[order], f"companion row {order}")

# Derive, rather than assume, the degree-two Weierstrass packet by Euclidean
# division in x at each t-order.
W0 = Q_rows[0]
U1 = sp.Poly(Q_rows[1], x).coeff_monomial(x**2)
W1 = simp(Q_rows[1] - U1 * W0)
require(sp.degree(W1, x) < 2, "Weierstrass order-one remainder degree")
U2, W2 = sp.div(sp.Poly(simp(Q_rows[2] - U1 * W1), x), sp.Poly(W0, x))
U2, W2 = simp(U2.as_expr()), simp(W2.as_expr())
require(sp.degree(W2, x) < 2, "Weierstrass order-two remainder degree")
U_packet = 1 + U1 * t + U2 * t**2
W_packet = W0 + W1 * t + W2 * t**2
packet_defect = sp.expand(Q_known - U_packet * W_packet)
for order in range(3):
    equal(packet_defect.coeff(t, order), 0, f"derived packet product row {order}")

B_packet = sp.expand(W_packet).coeff(x, 1)
C_packet = sp.expand(W_packet).coeff(x, 0)
disc = sp.expand(B_packet**2 - 4 * C_packet)
disc_rows = tuple(simp(disc.coeff(t, order)) for order in range(3))
disc_expected = (
    -24 / A5,
    sp.Rational(496, 3) / A5**2,
    -sp.Rational(179488, 135) / A5**3,
)
for order in range(3):
    equal(disc_rows[order], disc_expected[order], f"packet discriminant row {order}")

# Verify the two normalization germs directly modulo e^2+6/A5.
e = sp.symbols("e", nonzero=True)


def reduce_e_polynomial(expr):
    numerator = sp.together(expr * A5**12).as_numer_denom()[0]
    return simp(sp.rem(sp.Poly(numerator, e), sp.Poly(A5 * e**2 + 6, e)).as_expr())


r1 = -b / 2 - sp.Rational(31, 9) * e / A5
r2 = -d / 2 + 6 * b / A5 + sp.Rational(653, 30) * e / A5**2
branch = e + r1 * t + r2 * t**2
branch_value = sp.expand(Q_known.subs(x, branch))
for order in range(3):
    equal(reduce_e_polynomial(branch_value.coeff(t, order)), 0,
          f"normalization branch row {order}")
equal(r1 + r1.subs(e, -e), -b, "branch center row one")
equal(r2 + r2.subs(e, -e), -d + 12 * b / A5, "branch center row two")

# Endpoint and residual hostiles.
E_known = simp(1 - R_known.subs(p, 0).subs(y, sp.symbols("Y")))
# R_known has already had p,y replaced by source expressions, so rebuild the
# abstract residual for the endpoint projection.
Pabs, Yabs = sp.symbols("Pabs Yabs")
R_abstract = c20 * Pabs**2 + b * Yabs + d * Pabs * Yabs + c30 * Pabs**3
equal(1 - R_abstract.subs(Pabs, 0), 1 - b * Yabs,
      "known endpoint jet through the available layer")
require(simp(c20) != 0, "R=0 hostile contradicted by mandatory p2 coefficient")
R_disjoint = c20 * Pabs**2 + d * Pabs * Yabs + c30.subs(b, 0) * Pabs**3
equal(R_disjoint.subs(Pabs, 0), 0, "boundary-disjoint residual projection")
equal(R_disjoint,
      Pabs**2 * (c20 + c30.subs(b, 0) * Pabs) + Pabs * Yabs * d,
      "boundary-disjoint residual remains in (p2,py)")
Y2 = sp.symbols("Y2", nonzero=True)
require(sp.Poly(1 - Y2 * Yabs**2, Yabs).degree() == 2,
        "b=0 higher-y endpoint hostile is nonconstant")

# The first unknown residual layer really starts in the t^3 companion row.
c40, c21, c02 = sp.symbols("c40 c21 c02")
R_next = c40 * Pabs**4 + c21 * Pabs**2 * Yabs + c02 * Yabs**2
Q_next = simp(R_next.subs({Pabs: p, Yabs: y}) / t)
for order in range(3):
    equal(sp.expand(Q_next).coeff(t, order), 0, f"next residual absent row {order}")
equal(sp.expand(Q_next).coeff(t, 3), c40 + c21 * x + c02 * x**2,
      "first missing companion row")

# Fixed nodal defect is not invariant under an arbitrary target shear.
Avar, Cvar, Ivar = sp.symbols("Avar Cvar Ivar")
defect = Cvar**2 - Avar**3 + Ivar * Avar
defect_sheared = (Cvar + kappa * Avar) ** 2 - Avar**3 + Ivar * Avar
equal(defect_sheared - defect, 2 * kappa * Avar * Cvar + kappa**2 * Avar**2,
      "target shear changes fixed nodal defect")

D = sp.symbols("D", integer=True, positive=True)
require(all(4 * degree - 2 > 0 for degree in range(1, 100)),
        "THM3998 all-positive-degree gap hostile")

# Semantic record contains only exact, order-insensitive conclusions.
semantic_record = {
    "A_table": {
        f"{n},{j}": [weight, str(coefficient)]
        for (n, j), (weight, coefficient) in sorted(A_table.items())
    },
    "C_table": {
        f"{n},{j}": [weight, str(coefficient)]
        for (n, j), (weight, coefficient) in sorted(C_table.items())
    },
    "strata": {key: list(value) for key, value in strata.items()},
    "fixed_gauge_first_candidate": [4, 5],
    "pre_gauge_transferred_exclusions": [[3, 4], [4, 3]],
    "pre_gauge_floor_not_transferred": [4, 4],
    "Q_rows": [str(value) for value in Q_rows],
    "U": str(simp(U_packet)),
    "W": str(simp(W_packet)),
    "disc_rows": [str(value) for value in disc_rows],
    "next_row": "c40+c21*x+c02*x^2",
}
semantic_digest = sha256(
    json.dumps(semantic_record, sort_keys=True,
               separators=(",", ":")).encode("ascii")
).hexdigest()
if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
    require(semantic_digest == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

print("JC2_LIVE_34_INVARIANT_ATLAS_INDEPENDENT_AUDIT_20260824")
print("status=THM-4005_INDEPENDENT_AUDIT_PASS;VERIFIED_EXACT;JC2_OPEN")
print("scope=oriented_reduced_2:3_live_seam_in_fixed_THM3992_gauge;A5_nonzero;b,d_arbitrary")
print("weight_rule=x^j*t^n_has_weight_j-2n;target_scalars_deleted")
print("A_table=" + repr(tuple((key, *A_table[key]) for key in sorted(A_table))))
print("C_table=" + repr(tuple((key, *C_table[key]) for key in sorted(C_table))))
print("strata=" + repr(strata))
print("fixed_gauge_first_unrejected_support=4x5_only_on_b=d=0")
print("transfer=recorded_diagonal_scalings_translations_and_C_by_A_linear_shear_exclude_3x4_and_4x3")
print("transfer_firewall=4x5_floor_does_not_transfer;expanded_C_copy_of_A_can_cancel_to_pre_gauge_4x4")
print("hostiles=special_t2_cancellations;tail_collision;closed_support_monomial_overcount;nonlinear_target_change;R0;boundary_blind_p2")
print("companion=Q_rows_and_derived_Weierstrass_packet_and_two_normalization_germs_pass")
print("first_missing_row=t^3*(c40+c21*x+c02*x^2)")
print("quantifier_firewall=no_atlas_point_existence;no_higher_row_lift;no_global_factor_or_owner_census;no_arbitrary_target_invariance")
print("semantic_sha256=" + semantic_digest)
print(f"checks={CHECKS}")
print("RESULT=PASS;FIRST_FIXED_GAUGE_BOUNDARY=4x5;JC2=OPEN")
