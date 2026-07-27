"""Exact symbolic certificate for THM-2570.

The Jelonek surface of the sporadic Keller map is globally normalized by
A^2_(c,lambda).  On c != 0 it is the cylinder over the cusp
QQ[theta^2,theta^3], theta=2-c*lambda.  This companion verifies the finite
birational parametrization, inverse ring formulas, conductor identities,
singular locus, c=0 boundary, and the normalized finite-survivor section.

All checks are exact over QQ and avoid ``assert`` so ``python3 -O`` is a
substantive independent replay.
"""

import sympy as sp


a, b, c, lam, th = sp.symbols("a b c lambda theta")
x, y, z = sp.symbols("x y z")

L = 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2
T = 4 - 3 * b * c
S = 27 * a * c**2 - 9 * b * c + 8
R = b - 9 * a * c


def require(condition, message):
    """Raise in ordinary and optimized Python if an exact check fails."""
    if not condition:
        raise RuntimeError(message)


def zero(expr, message):
    require(sp.cancel(expr) == 0, message)


print("=" * 78)
print("[1] Global finite birational parametrization by A^2_(c,lambda)")
print("=" * 78)

aN = lam**2 * (3 - c * lam) / 27
bN = lam * (4 - c * lam) / 3
zero(L.subs({a: aN, b: bN}, simultaneous=True), "normalization misses L")

# Eliminate lambda from the two graph equations.  The only target equation is L.
graph_a = 27 * a - lam**2 * (3 - c * lam)
graph_b = 3 * b - lam * (4 - c * lam)
graph_gb = sp.groebner([graph_a, graph_b], lam, a, b, c, order="lex")
eliminated = [
    sp.factor(poly.as_expr())
    for poly in graph_gb.polys
    if not poly.as_expr().has(lam)
]
require(eliminated == [L], "normalization graph has unexpected target kernel")

monic_relation = lam**2 - 3 * b * lam + 27 * a
fraction_relation = (3 * b * c - 4) * lam - (27 * a * c - 3 * b)
zero(monic_relation.subs({a: aN, b: bN}, simultaneous=True), "lambda not integral")
zero(fraction_relation.subs({a: aN, b: bN}, simultaneous=True), "fraction fields differ")

print("  a=lambda^2(3-c lambda)/27")
print("  b=lambda(4-c lambda)/3")
print("  elimination kernel = (L)  [PASS]")
print("  lambda^2-3b lambda+27a=0: normalization is finite.  [PASS]")
print("  (3bc-4)lambda=27ac-3b: fraction fields agree.  [PASS]")
print("  Hence QQ[c,lambda] is the integral closure of QQ[a,b,c]/(L).")
print()

print("=" * 78)
print("[2] On c!=0 the surface is exactly a cusp cylinder")
print("=" * 78)

thetaN = 2 - c * lam
zero(T.subs({a: aN, b: bN}, simultaneous=True) - thetaN**2, "T != theta^2")
zero(S.subs({a: aN, b: bN}, simultaneous=True) - thetaN**3, "S != theta^3")

# The normalization is an elementary subintegral extension globally, not just
# a set-bijection visible after localizing at c.  The recovery identity proves
# B=A[theta] even though c need not be invertible.
lambda_from_theta = (3 * R + 3 * b * (2 - th)) / 4
zero(
    lambda_from_theta.subs(
        {a: aN, b: bN, th: thetaN}, simultaneous=True
    ) - lam,
    "theta does not generate the global normalization",
)

b_theta = (4 - th**2) / (3 * c)
a_theta = (th**3 - 3 * th**2 + 4) / (27 * c**2)
zero(a_theta - (th - 2) ** 2 * (th + 1) / (27 * c**2), "a factor failed")
zero(4 - 3 * b_theta * c - th**2, "theta inverse for b failed")
zero(27 * a_theta * c**2 - 9 * b_theta * c + 8 - th**3,
     "theta inverse for a failed")

print("  theta=2-c lambda, T=theta^2, S=theta^3.  [PASS]")
print("  lambda=(3(b-9ac)+3b(2-theta))/4, so B=A[theta] globally.  [PASS]")
print("  Since theta^2,theta^3 lie in A, normalization is subintegral.")
print("  A[c^-1] = QQ[c,c^-1,theta^2,theta^3].")
print("  normalization = QQ[c,c^-1,theta], a literal G_m x cusp cylinder.")
print("  inverse formulas:")
print("    b=(4-theta^2)/(3c)")
print("    a=(theta^3-3theta^2+4)/(27c^2)")
print("     =(theta-2)^2(theta+1)/(27c^2).  [PASS]")
print()

print("=" * 78)
print("[3] Exact conductor and its support")
print("=" * 78)

TN = sp.factor(T.subs({a: aN, b: bN}, simultaneous=True))
SN = sp.factor(S.subs({a: aN, b: bN}, simultaneous=True))
RN = sp.factor(R.subs({a: aN, b: bN}, simultaneous=True))
zero(TN - thetaN**2, "normalized T failed")
zero(SN - thetaN**3, "normalized S failed")
zero(TN * lam - 3 * RN, "T lambda != 3R")
zero(S - 2 * T + 3 * c * R, "S/T/R relation failed")
zero(4 - T - 3 * b * c, "unit relation failed")
conductor_support_gb = sp.groebner([L, T, R], a, b, c, order="lex")
conductor_support_basis = [
    sp.factor(poly.as_expr()) for poly in conductor_support_gb.polys
]
require(
    conductor_support_basis == [12 * a - b**2, 3 * b * c - 4],
    "unexpected conductor-support quotient",
)

print("  normalized conductor = T*QQ[c,lambda]=(2-c lambda)^2.")
print("  T lambda=3(b-9ac), and S=2T-3c(b-9ac).  [PASS]")
print("  contracted conductor in A is")
print("    (T,b-9ac)=(T,S).")
print("  Its support is T=0, b-9ac=0, namely")
print("    E={(4/(27t^2),4/(3t),t): t!=0}.")
print("  exact support quotient basis: [12a-b^2,3bc-4],")
print("    so A/(T,b-9ac)=QQ[c,c^-1] and B/A is one copy of E.  [PASS]")
print("  On c!=0 this is the cusp conductor (theta^2,theta^3).")
print()

print("=" * 78)
print("[4] Singular locus equals the conductor support")
print("=" * 78)

dLa = sp.diff(L, a)
dLb = sp.diff(L, b)
dLc = sp.diff(L, c)
zero(dLa - 2 * S, "L_a != 2S")
jacobian_gb = sp.groebner([L, dLa, dLb, dLc], a, b, c, order="lex")
jacobian_basis = [sp.factor(poly.as_expr()) for poly in jacobian_gb.polys]
require(jacobian_basis == [16 * a - b**3 * c, (3 * b * c - 4) ** 2],
        "unexpected singular-locus basis")

print("  L_a=2S; on L=L_a=0, S^2-T^3=27c^2L forces T=0.")
print("  Then L_b=2(b-9ac), so every singular point lies on E.")
print("  Conversely E kills L and all three partial derivatives.  [PASS]")
print("  exact Jacobian-ideal Groebner basis:")
print("   ", jacobian_basis)
print("  Its reduced support is E; the squared T records the cusp thickness.")
print()

print("=" * 78)
print("[5] The c=0 boundary is smooth and normalization is an isomorphism")
print("=" * 78)

zero(L.subs(c, 0) - (16 * a - b**2), "c=0 boundary equation failed")
a_boundary = sp.factor(aN.subs(c, 0))
b_boundary = sp.factor(bN.subs(c, 0))
require(a_boundary == lam**2 / 9, "boundary a failed")
require(b_boundary == 4 * lam / 3, "boundary b failed")
require(sp.diff(L, a).subs(c, 0) == 16, "boundary should be smooth")
require(thetaN.subs(c, 0) == 2, "boundary theta should be 2")

print("  V(L,c)=V(16a-b^2), a smooth parabola (L_a=16).")
print("  normalization boundary: (a,b,c)=(lambda^2/9,4lambda/3,0),")
print("  inverse lambda=3b/4.  Thus the normalization is an isomorphism there.")
print("  theta=2 on the whole boundary; lambda is the first-order coordinate")
print("  discarded when Phi_c collapses that parabola to (S,T)=(8,4).")
print()

print("=" * 78)
print("[6] Exact normalized finite-survivor section")
print("=" * 78)

q_theta = th**2 * (th - 2) ** 2 / (9 * c**2)
zero(b_theta**2 - 12 * a_theta - q_theta, "q(theta) failed")

x_surv = 2 * c / th**2
y_surv = (2 - th) * (3 * th + 2) / (6 * c)
z_surv = -th**2 * (th - 2) ** 2 * (th**2 + 4 * th + 2) / (8 * c**2)

radical = 1 + x * y
F1 = sp.expand(radical**3 * z + y**2 * radical * (4 + 3 * x * y))
F2 = sp.expand(y + 3 * x * radical**2 * z + 3 * x * y**2 * (4 + 3 * x * y))
F3 = sp.expand(2 * x - 3 * x**2 * y - x**3 * z)
source_sub = {x: x_surv, y: y_surv, z: z_surv}
target_sub = {a: a_theta, b: b_theta}
for actual, expected, label in zip(
    (F1.subs(source_sub), F2.subs(source_sub), F3.subs(source_sub)),
    (a_theta, b_theta, c),
    ("F1", "F2", "F3"),
):
    zero(actual - expected, f"normalized survivor {label} failed")

boundary_source = {x: 0, y: 4 * lam / 3, z: -7 * lam**2}
boundary_target = (lam**2 / 9, 4 * lam / 3, 0)
for actual, expected, label in zip(
    (F1.subs(boundary_source), F2.subs(boundary_source), F3.subs(boundary_source)),
    boundary_target,
    ("boundary F1", "boundary F2", "boundary F3"),
):
    zero(actual - expected, f"normalized survivor {label} failed")

print("  q=b^2-12a=theta^2(theta-2)^2/(9c^2).  [PASS]")
print("  for c theta !=0, the unique finite survivor is")
print("    x=2c/theta^2")
print("    y=(2-theta)(3theta+2)/(6c)")
print("    z=-theta^2(theta-2)^2(theta^2+4theta+2)/(8c^2)")
print("  direct substitution into F gives (a(theta,c),b(theta,c),c).  [PASS]")
print("  at c=0 the polynomial boundary survivor is")
print("    (x,y,z)=(0,4lambda/3,-7lambda^2), also checked directly.  [PASS]")
print()

print("=" * 78)
print("[7] Three structural theta sections")
print("=" * 78)

def values_at(theta_value):
    target = tuple(sp.factor(expr.subs(th, theta_value)) for expr in (a_theta, b_theta, c))
    if theta_value == 0:
        return target, None
    source = tuple(sp.factor(expr.subs(th, theta_value)) for expr in
                   (x_surv, y_surv, z_surv))
    return target, source

target0, source0 = values_at(0)
target2, source2 = values_at(2)
targetm1, sourcem1 = values_at(-1)
require(target0 == (sp.Rational(4, 27) / c**2, sp.Rational(4, 3) / c, c),
        "theta=0 target failed")
require(source0 is None, "theta=0 must have no survivor formula")
require(target2 == (0, 0, c) and source2 == (c / 2, 0, 0), "theta=2 failed")
require(
    targetm1 == (0, 1 / c, c)
    and sourcem1 == (2 * c, -1 / (2 * c), sp.Rational(9, 8) / c**2),
    "theta=-1 failed",
)

print("  theta=0: E, the conductor/singular curve; F-fibre empty.")
print("  theta=2: target c-axis (0,0,c), survivor (c/2,0,0).")
print("  theta=-1: target (0,1/c,c), survivor")
print("            (2c,-1/(2c),9/(8c^2)).")
print("  The factor a=(theta-2)^2(theta+1)/(27c^2) makes the two a=0")
print("  components and their double/simple multiplicities explicit.")
print()

print("ALL EXACT THM-2570 NORMALIZATION CHECKS PASSED")
