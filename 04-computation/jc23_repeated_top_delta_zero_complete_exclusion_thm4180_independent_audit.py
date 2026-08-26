#!/usr/bin/env python3
"""Independent normalized-projection audit for THM-4180.

This file does not import the primary source-pair certificate.  It rebuilds
the source in (X,T), verifies the Morse-resultant bridge, and reads the
row-A D_A wall from the *top* of the T-resultant rather than the bottom of
the source p-resultant.  It also reconstructs the Delta-zero faces directly.
"""

import sympy as sp


CHECKS = 0


def need(condition, message):
    global CHECKS
    CHECKS += 1
    if not bool(condition):
        raise RuntimeError(message)


def valuation(poly, variable):
    terms = sp.Poly(poly, variable).terms()
    need(bool(terms), "zero polynomial valuation")
    return min(monomial[0] for monomial, coefficient in terms if coefficient)


def exact_quotient(numerator, denominator, variable, message):
    quotient = sp.cancel(numerator/denominator)
    need(sp.denom(quotient) == 1, message + " polynomial quotient")
    return sp.Poly(quotient, variable)


def selected_face(poly, x, y, u, v, level):
    result = 0
    for (i, j), coefficient in sp.Poly(poly, x, y).terms():
        if u*i + v*j == level:
            result += coefficient*x**i*y**j
    return sp.factor(result)


X, T = sp.symbols("X T")
theta, phi, eta = sp.symbols("theta phi eta")
u = sp.symbols("u")
k0 = sp.Rational(2848, 45)
epsilon = -sp.Rational(1376, 135)

P = T + X**2*T**2
Y = X*T*P
G = sp.expand(
    -X**2*T/2 - 3*P + sp.Rational(8, 3)*P**2 + epsilon*P**3
    + k0*Y**2 + phi*P**2*Y + theta*P*Y**2
    + eta*P**3*Y - eta*Y**3
)
f = sp.expand(sp.cancel(sp.diff(G, X)/T))
h = sp.expand(sp.diff(G, T))
Hess = sp.factor(sp.det(sp.hessian(G, (X, T))))
critical_jacobian = sp.factor(sp.det(sp.Matrix((
    (sp.diff(f, X), sp.diff(f, T)),
    (sp.diff(h, X), sp.diff(h, T)),
))))
need(sp.factor(T*critical_jacobian - Hess - f*sp.diff(G, X, T)) == 0,
     "normalized Morse-resultant bridge")
need(sp.factor(f.subs(T, 0) + X) == 0,
     "universal T-zero f")
need(sp.factor(h.subs(T, 0) + (X**2 + 6)/2) == 0,
     "universal T-zero h")

# Direct Delta-zero face reconstruction.  The five inequalities are the
# primitive inward half-spaces of the claimed polygon; every one is active.
s, p, Q = sp.symbols("s p Q")
Hs = sp.expand(
    -3*p + sp.Rational(8, 3)*p**2 + epsilon*p**3 + k0*s**2*p**2
    + phi*s*p**3 + theta*s**2*p**3 + eta*s*p**3*(p - s**2)
)
F = sp.expand((s**2 - p)*(1 - Q*Hs) - Q*s**2/2)
support = tuple(
    powers for powers, coefficient
    in sp.Poly(F.subs({theta: 0, phi: 0, eta: 1}), s, p).terms()
    if coefficient != 0
)
halfspaces = ((1, 2, 2), (-1, 1, -2), (-1, -2, -11),
              (1, -1, -4), (1, 0, 0))
for uu, vv, level in halfspaces:
    need(all(uu*i + vv*j >= level for i, j in support),
         "support escaped Delta-zero polygon")
    need(sum(uu*i + vv*j == level for i, j in support) >= 2,
         "inactive Delta-zero polygon edge")
need(sp.factor(
    selected_face(F, s, p, 1, -1, -4)
    - Q*p**4*(epsilon + eta*s*p)
) == 0, "independent diagonal replacement face")
need(sp.factor(
    selected_face(F, s, p, -1, -2, -11)
    - Q*eta*s*p**3*(p - s**2)**2
) == 0, "independent repeated face")
need(sp.factor(
    selected_face(F, s, p, -1, 1, -2)
    - s**2*((1 - Q/2) - k0*Q*(s*p)**2 + eta*Q*(s*p)**3)
) == 0, "independent cubic face")

btop = sp.factor(epsilon + k0)
need(btop == sp.Rational(7168, 135), "row-D is empty at Delta zero")

# Row A away from D_A, specialized only at phi=0 for a clean independent
# reconstruction.  The primary source-pair computation retains arbitrary phi.
DA = sp.factor(4*theta*k0**2 - 27*eta**2)
fA = sp.expand(f.subs(phi, 0))
hA = sp.expand(h.subs(phi, 0))
need((sp.degree(fA, X), sp.degree(hA, X)) == (7, 8),
     "row-A normalized degrees")
need(sp.factor(sp.Poly(fA, X).LC() - 8*theta*T**7) == 0
     and sp.factor(sp.Poly(hA, X).LC() - 8*theta*T**7) == 0,
     "row-A normalized infinity gate")
RA = sp.resultant(fA, hA, X)
need(valuation(RA, T) == 42, "row-A normalized T-artifact")
QA = exact_quotient(RA, T**42*(6*T + 1)**2, T, "row-A")
need(QA.degree() == 19, "row-A normalized residual degree")
need(sp.factor(QA.TC() + 12288*theta**6) == 0,
     "row-A normalized bottom endpoint")
need(sp.factor(QA.LC() + 1458*theta*eta**4*DA**2) == 0,
     "row-A normalized top endpoint")

# Row A on D_A=0.  This is the dual audit of the primary p-adic cascade.
wall = {theta: 3*u**2, eta: 2*k0*u/3}
fw = sp.expand(f.subs(wall))
hw = sp.expand(h.subs(wall))
need((sp.degree(fw, X), sp.degree(hw, X)) == (7, 8),
     "row-A wall normalized degrees")
need(sp.factor(sp.Poly(fw, X).LC() - 24*u**2*T**7) == 0
     and sp.factor(sp.Poly(hw, X).LC() - 24*u**2*T**7) == 0,
     "row-A wall normalized infinity gate")
Rw = sp.resultant(fw, hw, X)
need(valuation(Rw, T) == 42, "row-A wall normalized T-artifact")
Q18 = exact_quotient(Rw, T**42*(6*T + 1)**2, T, "row-A wall")
J0 = 8544*phi + 22784*u + 1215*u**3
S0 = 18225*u**4 - 1515136*u**2 - 129777664
need(Q18.degree() == 18, "row-A wall generic degree")
need(sp.factor(Q18.TC() + 8957952*u**12) == 0,
     "row-A wall nonzero T-origin")
expected_q18 = (-sp.Rational(1092873416397493970665472,
                             5605041796875)*u**6*J0**2)
need(sp.factor(Q18.LC() - expected_q18) == 0,
     "row-A wall J top endpoint")

phi_J = -(22784*u + 1215*u**3)/8544
QJ = sp.Poly(sp.expand(Q18.as_expr().subs(phi, phi_J)), T)
need(QJ.degree() == 17, "row-A normalized J degree")
expected_q17 = (-sp.Rational(134737936586375168, 69198046875)
                *u**6*S0**2)
need(sp.factor(QJ.nth(17) - expected_q17) == 0,
     "row-A normalized S top endpoint")
q16_numerator = sp.together(QJ.nth(16)).as_numer_denom()[0]
need(sp.gcd(sp.Poly(S0, u), sp.Poly(q16_numerator, u)).degree() == 0,
     "row-A normalized terminal top coefficient")

# Rows B and C are recomputed in the normalized projection rather than
# imported from the primary source pair.
fB = sp.expand(f.subs(theta, 0))
hB = sp.expand(h.subs(theta, 0))
need((sp.degree(fB, X), sp.degree(hB, X)) == (6, 7),
     "row-B normalized degrees")
need(sp.factor(sp.Poly(fB, X).LC() - 7*T**6*(phi + eta*T)) == 0
     and sp.factor(sp.Poly(hB, X).LC() - T**6*(7*phi + 8*eta*T)) == 0,
     "row-B normalized infinity gate")
RB = sp.resultant(fB, hB, X)
need(valuation(RB, T) == 30, "row-B normalized T-artifact")
QB = exact_quotient(RB, T**30*(6*T + 1)**2, T, "row-B")
need(QB.degree() == 17, "row-B normalized residual degree")
need(sp.factor(QB.TC() + sp.Rational(50421, 32)*phi**5) == 0,
     "row-B normalized bottom endpoint")
need(sp.factor(QB.LC() - sp.Rational(531441, 4)*eta**7) == 0,
     "row-B normalized top endpoint")

fC = sp.expand(fB.subs(phi, 0))
hC = sp.expand(hB.subs(phi, 0))
need((sp.degree(fC, X), sp.degree(hC, X)) == (6, 7),
     "row-C normalized degrees")
need(sp.factor(sp.Poly(fC, X).LC() - 7*eta*T**7) == 0
     and sp.factor(sp.Poly(hC, X).LC() - 8*eta*T**7) == 0,
     "row-C normalized infinity gate")
RC = sp.resultant(fC, hC, X)
need(valuation(RC, T) == 32, "row-C normalized T-artifact")
QC = exact_quotient(RC, T**32*(6*T + 1)**2, T, "row-C")
need(QC.degree() == 15, "row-C normalized residual degree")
need(sp.factor(
    QC.TC() - sp.Rational(37845999468607963136, 61509375)*eta
) == 0, "row-C normalized bottom endpoint")
need(sp.factor(QC.LC() - sp.Rational(531441, 4)*eta**7) == 0,
     "row-C normalized top endpoint")

# Exact degree-to-length and response consequences.
need((19 + 4, 17 + 4, 15 + 4) == (23, 21, 19),
     "open-row critical lengths")
need((18 + 4, 17 + 4, 16 + 4) == (22, 21, 20),
     "row-A wall critical lengths")
for full_n, finite_n, length, defect, finite_index in (
        (25, 19, 23, 18, 15),
        (23, 17, 21, 16, 13),
        (21, 15, 19, 14, 11)):
    need(2*(full_n - length) < defect, "independent full response")
    need(2*finite_n - length - 1 + 3 < finite_n - 1,
         "independent finite open-row response")
for length in (22, 21, 20):
    need(2*(25 - length) < 18, "independent row-A wall full response")
    need(2*(19 + 3 - length) + 3 < 15,
         "independent row-A wall finite response")

print("THM4180_INDEPENDENT_NORMALIZED_AUDIT")
print("geometry=Delta_zero_halfspaces_and_replacement_face_recomputed")
print("row_partition=A:theta!=0;B:theta=0,phi!=0;C:theta=phi=0;Btop=7168/135;D:empty")
print("bridge=T*detD(f,h)-detHess(G)-f*G_XT=0")
print("row_A_normalized=Q19;DA_wall_Q_degrees=18,17,16")
print("row_A_DA_terminal=gcd(S0,Q16_top)=1")
print("row_B_normalized=Q17;row_C_normalized=Q15")
print("critical_lengths=A:23_or_22,21,20;B:21;C:19")
print("responses=all_full_and_finite_bounds_strict")
print(f"checks={CHECKS}")
print("verdict=INDEPENDENT_DELTA_ZERO_ACCEPT")
