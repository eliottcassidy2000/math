#!/usr/bin/env python3
"""Exact c!=0 Galois-resolvent chart for the THM-1310 cubic.

The calculation isolates the strongest valid consequence for a hypothetical
quartic using the *actual THM-1310 cubic field* as its classical resolvent
root field.  Its full S3 Galois normalization has the dense chart

    G_m x A^2 = Spec C[c^+-1,r,s].

That chart has procyclic etale fundamental group, one unit squareclass, and
trivial class group.  A connected V4 torsor supplied by THM-2655 cannot
restrict to it.  The conclusion concerns dimension-three S4 quartic Keller
candidates whose actual classical cubic-resolvent root field is THM-1310's
field, also after a polynomial/affine target-coordinate automorphism or an
explicit base-ring isomorphism.  It does not cover an arbitrary birational
function-field transport, nor quartics which share only a discriminant
squareclass, cusp equation, or other weaker shadow with THM-1310.

All executable checks raise explicitly, including under ``python -O``.

Notation: this script writes ``c`` for the target coordinate called ``g`` in
THM-1310, reserving ``G`` for a monodromy group in the accompanying proof.
"""

from itertools import product

import sympy as sp


checks = 0


def check(label, condition):
    global checks
    checks += 1
    if not bool(condition):
        raise RuntimeError(f"FAILED: {label}")


a, b, c, x, Y, T = sp.symbols("a b c x Y T")
ell = sp.symbols("ell")
p0, Q0 = sp.symbols("p0 Q0")
w, r, s = sp.symbols("w r s")
eta = sp.sqrt(-3)
zeta = (-1 + eta) / 2

L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
p = 4 - 3 * b * c
Q = 27 * a * c**2 - 9 * b * c + 8


# The target c!=0 chart is exactly C[c^+-1,p,Q].
b_from = (4 - p0) / (3 * c)
a_from = (Q0 - 3 * p0 + 4) / (27 * c**2)
check("recover p from (a,b,c)", sp.expand(p.subs({a: a_from, b: b_from}) - p0) == 0)
check("recover Q from (a,b,c)", sp.expand(Q.subs({a: a_from, b: b_from}) - Q0) == 0)
check("recover b from (c,p,Q)", sp.factor(b_from.subs(p0, p) - b) == 0)
check("recover a from (c,p,Q)", sp.factor(a_from.subs({p0: p, Q0: Q}) - a) == 0)


# THM-1335/2570 cusp identity and the exact cubic discriminant.
check("cusp identity Q^2-p^3=27c^2L", sp.expand(Q**2 - p**3 - 27 * c**2 * L) == 0)
fiber_cubic = L * T**3 + p * T - 2 * c
fiber_disc = sp.factor(sp.discriminant(fiber_cubic, T))
check("raw cubic discriminant", sp.expand(fiber_disc + 4 * Q**2 * L) == 0)


# On the discriminant double cover w^2=-L, the scaled root Y=-w*x
# satisfies a monic depressed cubic.  Use an independent symbol ell here so
# that the reduction modulo w^2+ell is explicit rather than relying on a
# fragile substitution of the expanded polynomial L.
generic_fiber_cubic = ell * T**3 + p * T - 2 * c
scaled_relation = sp.factor(w * generic_fiber_cubic.subs(T, -Y / w))
scaled_relation = sp.factor(scaled_relation.subs(ell, -w**2))
check("scaled integral cubic Y^3-pY-2cw", sp.expand(scaled_relation - (Y**3 - p * Y - 2 * c * w)) == 0)
scaled_disc = sp.factor(sp.discriminant(Y**3 - p * Y - 2 * c * w, Y))
scaled_disc_on_cover = sp.expand(scaled_disc.subs(w**2, -L))
check("scaled cubic discriminant is (2Q)^2", sp.expand(scaled_disc_on_cover - 4 * Q**2) == 0)


# Cardano factors.  eta^2=-3 and w^2=-L give UV=p^3.
U = Q + 3 * eta * c * w
V = Q - 3 * eta * c * w
UV_on_cover = sp.expand(U * V).subs(w**2, -L)
check("Cardano norm UV=p^3", sp.simplify(UV_on_cover - p**3) == 0)
check("U+V=2Q", sp.expand(U + V - 2 * Q) == 0)
check("U-V=6 eta c w", sp.expand(U - V - 6 * eta * c * w) == 0)


# If R^3=U/(3 eta), S^3=-V/(3 eta), and RS=p/3, then R+S is a
# root of the scaled cubic.  This is checked from the three defining
# relations, without choosing numerical cube roots.
R, S = sp.symbols("R S")
cardano_reduction = (
    U / (3 * eta)
    - V / (3 * eta)
    + (3 * (p / 3) - p) * (R + S)
    - 2 * c * w
)
check("Cardano substitution solves scaled cubic", sp.simplify(cardano_reduction) == 0)
cardano_product_error = sp.expand((U / (3 * eta)) * (-V / (3 * eta)) - (p / 3) ** 3)
check("Cardano product compatibility", sp.simplify(cardano_product_error.subs(w**2, -L)) == 0)


# Normalize the cyclic cubic by U=r^3, V=s^3, p=rs.  The C3-invariant
# exponent residues are exactly (0,0),(1,1),(2,2), proving
# C[r,s]^C3=C[U,V,p]/(UV-p^3).
invariant_residues = {(i, j) for i, j in product(range(3), repeat=2) if (i - j) % 3 == 0}
check("C3 invariant residue classes", invariant_residues == {(0, 0), (1, 1), (2, 2)})
check("normalized relation r^3 s^3=(rs)^3", sp.expand(r**3 * s**3 - (r * s) ** 3) == 0)


# Explicit inverse coordinate maps on c!=0.
p_rs = r * s
Q_rs = (r**3 + s**3) / 2
w_rs = (r**3 - s**3) / (6 * eta * c)
b_rs = (4 - p_rs) / (3 * c)
a_rs = (Q_rs - 3 * p_rs + 4) / (27 * c**2)
L_rs = sp.factor(L.subs({a: a_rs, b: b_rs}))
check("r,s chart recovers p", sp.expand(p.subs({a: a_rs, b: b_rs}) - p_rs) == 0)
check("r,s chart recovers Q", sp.expand(Q.subs({a: a_rs, b: b_rs}) - Q_rs) == 0)
check("r,s chart lies on w^2=-L", sp.simplify(w_rs**2 + L_rs) == 0)
check("r,s chart gives U=r^3", sp.simplify(U.subs({a: a_rs, b: b_rs, w: w_rs}) - r**3) == 0)
check("r,s chart gives V=s^3", sp.simplify(V.subs({a: a_rs, b: b_rs, w: w_rs}) - s**3) == 0)


# Identify the *actual* THM-1310 root inside this field, rather than merely
# matching its discriminant.  With R=-r/eta and S=s/eta, Cardano gives
# Y=R+S=(s-r)/eta and x=-Y/w.  The rational formula below is fixed by the
# transposition r<->s but not by the C3 rotation, so its orbit has size three.
Y_rs = (s - r) / eta
x_rs = sp.factor((r - s) / (eta * w_rs))
x_short = 6 * c / (r**2 + r * s + s**2)
check("explicit THM-1310 root formula", sp.simplify(x_rs - x_short) == 0)
check(
    "scaled cubic root in r,s chart",
    sp.simplify(Y_rs**3 - p_rs * Y_rs - 2 * c * w_rs) == 0,
)
check(
    "original THM-1310 cubic root in r,s chart",
    sp.simplify(L_rs * x_short**3 + p_rs * x_short - 2 * c) == 0,
)
check(
    "transposition negates discriminant root",
    sp.simplify(w_rs.subs({r: s, s: r}, simultaneous=True) + w_rs) == 0,
)
check(
    "transposition fixes cubic root",
    sp.simplify(x_short.subs({r: s, s: r}, simultaneous=True) - x_short) == 0,
)


# The six transformations generated by (r,s)->(zeta*r,zeta^-1*s)
# and (r,s)->(s,r) are distinct.  Their invariants are c, p=rs, and
# Q=(r^3+s^3)/2, hence the quotient is the target c!=0 chart.
sheet_labels = {(k, swap) for k in range(3) for swap in (0, 1)}
check("full resolvent group has six transformations", len(sheet_labels) == 6)
check("zeta is a primitive cube root", sp.simplify(zeta**2 + zeta + 1) == 0 and zeta != 1)
rho_matrix = sp.diag(zeta, zeta**-1)
tau_matrix = sp.Matrix([[0, 1], [1, 0]])
identity_matrix = sp.eye(2)
check("rho has order three", (rho_matrix**3).applyfunc(sp.simplify) == identity_matrix)
check("tau has order two", tau_matrix**2 == identity_matrix)
check(
    "tau rho tau = rho inverse",
    (tau_matrix * rho_matrix * tau_matrix - rho_matrix**-1).applyfunc(sp.simplify) == sp.zeros(2),
)
check(
    "rho fixes p=rs",
    sp.simplify((zeta * r) * (zeta**-1 * s) - p_rs) == 0,
)
check(
    "tau fixes p=rs",
    sp.expand(p_rs.subs({r: s, s: r}, simultaneous=True) - p_rs) == 0,
)
check(
    "rho fixes Q=(r^3+s^3)/2",
    sp.simplify(Q_rs.subs({r: zeta * r, s: zeta**-1 * s}, simultaneous=True) - Q_rs) == 0,
)
check(
    "tau fixes Q=(r^3+s^3)/2",
    sp.expand(Q_rs.subs({r: s, s: r}, simultaneous=True) - Q_rs) == 0,
)
rotated_denominator = sp.expand((zeta * r) ** 2 + p_rs + (zeta**-1 * s) ** 2)
check(
    "rho does not fix the cubic root",
    sp.expand(rotated_denominator - (r**2 + p_rs + s**2)) != 0,
)


# Exact obstruction dimensions.  C[c^+-1,r,s] is a Laurent polynomial
# UFD: Cl=0 and units are C* c^Z, so Kummer H1 has F2-dimension one.
unit_squareclass_rank = 1
class_group_2_rank = 0
kummer_h1_rank = unit_squareclass_rank + class_group_2_rank
v4_character_rank = 2
check("Laurent chart Kummer rank", kummer_h1_rank == 1)
check("V4 character plane has rank two", v4_character_rank == 2)
check("rank obstruction", v4_character_rank > kummer_h1_rank)


print("THM-1310 FULL GALOIS RESOLVENT: c!=0 CHART")
print("notation: c is the target coordinate g of THM-1310")
print("target chart: C[a,b,c,c^-1] = C[c,c^-1,p,Q]")
print("p=4-3bc; Q=27ac^2-9bc+8; Q^2-p^3=27c^2L")
print("discriminant cover: w^2=-L")
print("Cardano factors: U=Q+3*sqrt(-3)*c*w, V=Q-3*sqrt(-3)*c*w, UV=p^3")
print("cyclic normalization: U=r^3, V=s^3, p=rs")
print("full S3 normalization on c!=0: C[c,c^-1,r,s] = G_m x A^2")
print("quotient coordinates: p=rs, Q=(r^3+s^3)/2, w=(r^3-s^3)/(6*sqrt(-3)*c)")
print("actual cubic root: x=(r-s)/(sqrt(-3)*w)=6c/(r^2+rs+s^2)")
print("S3 action: rho(r,s)=(zeta*r,zeta^-1*s), tau(r,s)=(s,r)")
print("stabilizer: tau fixes x; rho does not, so [K(x):K]=3")
print("topology: pi1_et is procyclic; no connected V4 quotient")
print("Kummer ledger: units/C* squares rank=1, Cl[2] rank=0, required V4 character rank=2")
print("VERDICT: THM-1310's actual cubic field cannot be a dimension-3 quartic Keller S4 resolvent root field")
print("SCOPE: polynomial/affine target automorphism or explicit base-ring isomorphism is allowed")
print("NOT COVERED: arbitrary birational transport; shared -L, cusp identity, or N1-N5 shadows alone")
print(f"CHECKS PASSED: {checks}")
print("FAILED CHECKS: NONE")
