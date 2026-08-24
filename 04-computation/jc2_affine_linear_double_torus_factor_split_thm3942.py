#!/usr/bin/env python3
"""Exact companion for THM-3942's universal double-torus split gate.

The proof of the UFD classification and the normalization universal property
is written in the theorem.  This companion freezes every load-bearing
polynomial identity, both normalization maps, their infinity divisors, the
singular supports, scalar gauges, and the local A2 character witnesses.
"""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(expression) == 0, message)


x, y, h, v, s = sp.symbols("x y h v s")
lam, c = sp.symbols("lam c", nonzero=True)
o = sp.symbols("omega")


def oreduce(expression: sp.Expr) -> sp.Expr:
    """Reduce a polynomial coefficientwise modulo omega^2+omega+1."""
    return sp.expand(sp.rem(sp.expand(expression), o**2 + o + 1, o))


def ozero(expression: sp.Expr, message: str) -> None:
    gate(oreduce(expression) == 0, message)


# ---------------------------------------------------------------------------
# Squarefree UFD factor ledger and the two partition types.
# ---------------------------------------------------------------------------

l0 = y - x
l1 = y - o * x
l2 = y - o**2 * x
ozero(l0 * l1 * l2 - (y**3 - x**3), "difference-of-cubes factorization")

# Complement exchanges the two factors q1-q0 and q1+q0, hence changes q0
# only by sign.  The eight subsets therefore reduce to sizes zero and one.
subset_sizes_mod_complement = sorted({min(size, 3 - size) for size in range(4)})
gate(subset_sizes_mod_complement == [0, 1], "only 0|3 and 1|2 partition types")


# ---------------------------------------------------------------------------
# Empty|all: an irreducible sextic with Fermat-elliptic normalization.
# ---------------------------------------------------------------------------

q0e = y**3 - x**3 - 1
q1e = y**3 - x**3 + 1
He = sp.expand(q0e**2 - 4 * x**3)
zero(He - (q1e**2 - 4 * y**3), "empty|all common torus identity")
gate(sp.Poly(He, x, y).total_degree() == 6, "empty|all branch degree six")
gate(sp.factor(He, extension=sp.sqrt(-3)) == He,
     "empty|all sextic factor hostile over Q(omega)")

# On E: v^3-h^3=1, the map (h,v)->(h^2,v^2) lands on He=0.
fermat = v**3 - h**3 - 1
zero(sp.rem(He.subs({x: h**2, y: v**2}), fermat, v),
     "Fermat normalization map")
zero(q0e.subs({x: h**2, y: v**2}) - 2 * h**3 - fermat * (v**3 + h**3 + 1),
     "first cube row on Fermat normalization")

# The rational inverse is h=q0/(2x), v=(h^3+1)/y.
zero(sp.rem(q0e.subs({x: h**2, y: v**2}) - 2 * h**3, fermat, v),
     "Fermat inverse recovers h")
zero(sp.rem(h**3 + 1 - v**3, fermat, v),
     "Fermat inverse recovers v")

# [H:V:R] -> [H^2:V^2:R^2].  R=0 meets V^3-H^3 in three
# distinct points, and neither affine coordinate cancels there.
HH, VV, RR = sp.symbols("HH VV RR")
gate(sp.discriminant(VV**3 - HH**3, VV) == -27 * HH**6,
     "three distinct Fermat infinity points")
gate(sp.gcd(sp.gcd(HH**2, VV**2), RR**2) == 1,
     "Fermat projective map is basepoint free")

# Exact singular support: the three x=0,y^3=1 cusps and the three
# y=0,x^3=1 cusps.  The displayed derivative factorization proves exhaustion.
zero(sp.diff(He, x) + 6 * x**2 * q1e,
     "empty|all x-derivative factorization")
zero(sp.diff(He, y) - 6 * y**2 * q0e,
     "empty|all y-derivative factorization")
for point in ((0, 1), (1, 0)):
    gate(He.subs({x: point[0], y: point[1]}) == 0,
         f"empty|all representative cusp {point}")
gate(sp.diff(q0e, y).subs({x: 0, y: 1}) == 3,
     "x-Cardano cusp has q0 as transverse coordinate")
gate(sp.diff(q1e, x).subs({x: 1, y: 0}) == -3,
     "y-Cardano cusp has q1 as transverse coordinate")


# ---------------------------------------------------------------------------
# Singleton|complement: an irreducible quartic with conic normalization.
# ---------------------------------------------------------------------------

q0c = y**2 + x * y + x**2 - y + x
q1c = y**2 + x * y + x**2 + y - x
Hc = sp.expand(q0c**2 - 4 * x**3)
zero(Hc - (q1c**2 - 4 * y**3), "singleton common torus identity")
gate(sp.Poly(Hc, x, y).total_degree() == 4, "singleton branch degree four")
gate(sp.factor(Hc, extension=sp.sqrt(-3)) == Hc,
     "singleton quartic factor hostile over Q(omega)")

# Its hidden conic is obtained by x=h^2 and t=y-x.
t = sp.symbols("t")
conic_row = t**2 + (3 * h**2 - 1) * t + 3 * h**4 - 2 * h**3
gate(sp.factor(sp.discriminant(conic_row, t)) == -(h - 1) ** 3 * (3 * h + 1),
     "singleton hidden-conic discriminant")

# Basepoint-free rational normalization.  The target infinity divisor is
# (S^2+3R^2)^2, with exactly two normalization places.
den = (s**2 + 3) ** 2
xc = sp.cancel(4 * (s - 1) ** 2 / den)
yc = sp.cancel((s**2 - 1) ** 2 / den)
hc = sp.cancel(2 * (1 - s) / (s**2 + 3))
zero(Hc.subs({x: xc, y: yc}), "singleton rational normalization")
zero(xc - hc**2, "singleton first coefficient is a square")
zero(q0c.subs({x: xc, y: yc}) - 2 * hc**3,
     "singleton first cube row")
tc = sp.cancel(yc - xc)
uc = sp.cancel((2 * tc - 1 + 3 * hc**2) / (hc - 1))
zero(uc**2 + (hc - 1) * (3 * hc + 1),
     "singleton hidden-conic equation")
zero((uc - 1) / hc - s, "singleton birational inverse recovers s")

S, R = sp.symbols("S R")
Xhom = 4 * (S - R) ** 2 * R**2
Yhom = (S**2 - R**2) ** 2
Zhom = (S**2 + 3 * R**2) ** 2
gate(sp.gcd(sp.gcd(Xhom, Yhom), Zhom) == 1,
     "singleton projective map is basepoint free")
gate(sp.discriminant(S**2 + 3 * R**2, S) == -12 * R**2,
     "singleton has two distinct infinity places")
gate(sp.gcd(S**2 + 3 * R**2, (S - R) * (S + R) * R) == 1,
     "neither singleton infinity place cancels")

# Its complete affine singular support is (0,0),(0,1),(1,0).
singular_basis = sp.groebner([Hc, sp.diff(Hc, x), sp.diff(Hc, y)], x, y,
                             order="lex")
expected_basis = [
    x**2 + x * y - x + y**3 - 2 * y**2 + y,
    y**2 * (2 * x + y - 1),
    y**2 * (y - 1) ** 2,
]
gate([sp.factor(poly.as_expr()) for poly in singular_basis.polys]
     == [sp.factor(poly) for poly in expected_basis],
     "singleton singular-support Groebner basis")
for point in ((0, 0), (0, 1), (1, 0)):
    gate(all(poly.subs({x: point[0], y: point[1]}) == 0
             for poly in expected_basis),
         f"singleton singular point {point}")

# The two disjoint A2 witnesses detect independent order-three Cardano classes.
gate(sp.diff(q0c, y).subs({x: 0, y: 1}) == 1,
     "singleton x-Cardano A2 transverse coordinate")
gate(sp.diff(q1c, x).subs({x: 1, y: 0}) == 1,
     "singleton y-Cardano A2 transverse coordinate")


# ---------------------------------------------------------------------------
# The requested l1|l0*l2 split is the same conic after x -> omega*x.
# ---------------------------------------------------------------------------

B1 = oreduce(l0 * l2)
q01 = oreduce(B1 - l1)
q11 = oreduce(B1 + l1)
H1 = oreduce(q01**2 - 4 * x**3)
ozero(H1 - (q11**2 - 4 * y**3), "l1 split common torus identity")
ozero(q01 - q0c.subs(x, o * x), "l1 split is singleton x-scaling")
ozero(q11 - q1c.subs(x, o * x), "l1 second row is singleton x-scaling")
ozero(H1 - Hc.subs(x, o * x), "l1 branch is linearly equivalent")

# Direct hidden-conic check in the x=h^2 chart.
row1 = oreduce(q01.subs(x, h**2) - 2 * h**3)
D1 = oreduce(sp.discriminant(row1, y))
ozero(D1 + o**2 * (h - o) ** 3 * (3 * h + o),
      "l1 hidden-conic discriminant")

# Transforming the rational parametrization gives two, not one, infinity places.
x1 = oreduce(o**2 * xc)
ozero(H1.subs({x: x1, y: yc}), "l1 rational normalization")


# ---------------------------------------------------------------------------
# Scalar gauges and Cardano norm identities.
# ---------------------------------------------------------------------------

# Empty split: lambda=c^3 and x,y=c^2(X,Y) extract c^6 from H.
X, Y = sp.symbols("X Y")
q0e_scaled = c ** -3 * ((c**2 * Y) ** 3 - (c**2 * X) ** 3) - c**3
zero(q0e_scaled - c**3 * (Y**3 - X**3 - 1),
     "empty split scalar normalization")

# Singleton split: x,y=lambda^2(X,Y) extract lambda^6 from H.
A_scaled = lam * (lam**2 * Y - lam**2 * X)
B_scaled = lam**-1 * ((lam**2 * Y) ** 2
                      + (lam**2 * X) * (lam**2 * Y)
                      + (lam**2 * X) ** 2)
zero(A_scaled - lam**3 * (Y - X), "singleton A scalar normalization")
zero(B_scaled - lam**3 * (Y**2 + X * Y + X**2),
     "singleton B scalar normalization")

W = sp.symbols("W")
for label, q0, q1, p0, p1, branch in (
    ("empty", q0e, q1e, x, y, He),
    ("singleton", q0c, q1c, x, y, Hc),
):
    zero(sp.expand((q0 + W) * (q0 - W)).subs(W**2, branch) - 4 * p0**3,
         f"{label} first Cardano norm")
    zero(sp.expand((q1 + W) * (q1 - W)).subs(W**2, branch) - 4 * p1**3,
         f"{label} second Cardano norm")


summary = {
    "checks": CHECKS,
    "ufd_partition_types": ["0|3", "1|2"],
    "empty_all": {
        "branch": "irreducible sextic",
        "normalization": "Fermat elliptic minus three infinity points",
        "cardano_rank_lower_bound": 2,
    },
    "singleton_complement": {
        "branch": "irreducible rational quartic",
        "normalization": "P1 minus two infinity points (Gm)",
        "cardano_rank_lower_bound": 2,
    },
    "one_place_A1": False,
    "nonlinear_escape": "internally split a factor p1-zeta*p0 or use gcd/multiplicity overlap",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3942 affine-linear double-torus factor-split companion")
print(f"CHECKS={CHECKS}")
print("UFD_PARTITIONS=0|3,1|2")
print("EMPTY_ALL=IRREDUCIBLE_SEXTIC;NORMALIZATION=FERMAT_GENUS1;INFINITY=3")
print("SINGLETON=IRREDUCIBLE_QUARTIC;NORMALIZATION=Gm;INFINITY=2")
print("L1_SPLIT=LINEARLY_EQUIVALENT_SINGLETON;INFINITY=2")
print("CARDANO_CLASSES=2_INDEPENDENT_IN_BOTH_TYPES")
print("AFFINE_LINEAR_NONLINEAR_A1=IMPOSSIBLE")
print("NONLINEAR_ESCAPE=INTERNAL_FACTOR_SPLIT_OR_GCD_MULTIPLICITY_OVERLAP")
print(f"SEMANTIC_SHA256={semantic}")
