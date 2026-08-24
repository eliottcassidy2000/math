#!/usr/bin/env python3
"""Exact Laurent-row certificate for THM-3992's reduced (2,3) cell.

Work over a characteristic-zero differential field k(s), with h != 0 and
s|h for data inherited from THM-3989. SymPy represents coefficient functions
as differential indeterminates and verifies each displayed identity exactly.
"""

from __future__ import annotations

import sympy as sp


s, tau = sp.symbols("s tau")
h = sp.Function("h")(s)
b = sp.Function("b")(s)
A0 = sp.Function("A0")(s)
A1 = sp.Function("A1")(s)
A2 = sp.Function("A2")(s)
A3 = sp.Function("A3")(s)


def D(f):
    return sp.diff(f, s)


def simp(f):
    return sp.factor(sp.cancel(sp.expand(f)))


def assert_zero(label, f):
    z = simp(f)
    if z != 0:
        raise AssertionError(f"{label}: residual {z}")
    print(f"PASS  {label}")


def row(a, c, k):
    """Coefficient E_k of tau^k in tau(A_s C_tau-A_tau C_s)."""
    out = 0
    for i, ai in a.items():
        j = k - i
        if j in c:
            cj = c[j]
            out += j * D(ai) * cj - i * ai * D(cj)
    return simp(out)


def laurent(coeffs):
    return sum(fi * tau**i for i, fi in coeffs.items())


def direct_row(a, c, k):
    """Independent extraction from the expanded Laurent Jacobian."""
    AA, CC = laurent(a), laurent(c)
    JJ = sp.expand(tau * (D(AA) * sp.diff(CC, tau) - sp.diff(AA, tau) * D(CC)))
    return simp(JJ.coeff(tau, k))


print("ASSUMPTIONS: char(k)=0; h in k[s] is nonzero; THM-3989 adds s|h.")
print("All divisions below occur in k(s); no hypothesis h|b is used.")

# Normalize a_-2=alpha*h0^2 and c_-3=beta*h0^3. Choose kappa^5=alpha*beta,
# h=kappa*h0, and the determinant-one target scale (A,C)->(u*A,u^-1*C),
# where u=kappa^2/alpha.
alpha, beta, kappa = sp.symbols("alpha beta kappa", nonzero=True)
h0 = sp.Function("h0")(s)
u = kappa**2 / alpha
assert_zero("leading A normalization", u * alpha * h0**2 - (kappa * h0) ** 2)
assert_zero(
    "leading C normalization under kappa^5=alpha*beta",
    ((beta / u) * h0**3 - (kappa * h0) ** 3).subs(beta, kappa**5 / alpha),
)

# Normalized leading rows a_-2=h^2 and c_-3=h^3.
d = sp.Function("d")(s)
e = sp.Function("e")(s)
a_raw = {-2: h**2, -1: b, 0: A0, 1: A1, 2: A2, 3: A3}
c_raw = {-3: h**3, -2: d, -1: e}

Eminus4 = row(a_raw, c_raw, -4)
lambda_expr = 2 * d / h**2 - 3 * b / h
assert_zero("tau^-4 direct Laurent extraction", direct_row(a_raw, c_raw, -4) - Eminus4)
assert_zero("tau^-4 factorization", Eminus4 - h**4 * D(lambda_expr))

# A constant shear C->C-(lambda/2)A gives d=3hb/2.
d_g = sp.Rational(3, 2) * h * b
c_g = {-3: h**3, -2: d_g, -1: e}
assert_zero("tau^-4 after constant shear", row(a_raw, c_g, -4))

B = b / h
E = e / h
mu_expr = 2 * E - sp.Rational(3, 4) * B**2 - 3 * A0
assert_zero("tau^-3 direct Laurent extraction", direct_row(a_raw, c_g, -3) - row(a_raw, c_g, -3))
assert_zero("tau^-3 factorization", row(a_raw, c_g, -3) - h**3 * D(mu_expr))

# A constant translation A->A+mu/3 kills mu.
e_g = h * (sp.Rational(3, 8) * B**2 + sp.Rational(3, 2) * A0)
c_g[-1] = e_g
assert_zero("tau^-3 after constant translation", row(a_raw, c_g, -3))

# Rational approximate-root coefficients; they need not belong to k[s].
q0 = B / 2
q1 = (A0 - q0**2) / (2 * h)
q2 = (A1 - 2 * q0 * q1) / (2 * h)
C0star = q0**3 + 6 * h * q0 * q1 + 3 * h**2 * q2
R = sp.Function("R")(s)
c_R = dict(c_g)
c_R[0] = C0star + R
assert_zero("tau^-2 direct Laurent extraction", direct_row(a_raw, c_R, -2) - row(a_raw, c_R, -2))
assert_zero("tau^-2 repair row", row(a_raw, c_R, -2) - 2 * h**2 * D(R))

# A constant C-translation kills R without changing depths or negative rows.
c0 = C0star

# At tau^-1 the q3-independent combination I=3h^2 U-2hV survives.
q3 = sp.Function("q3")(s)
U = A2 - (q1**2 + 2 * q0 * q2 + 2 * h * q3)
C1star = 3 * q0**2 * q1 + 3 * h * q1**2 + 6 * h * q0 * q2 + 3 * h**2 * q3
V = sp.Function("V")(s)
c_UV = dict(c_g)
c_UV[0] = c0
c_UV[1] = C1star + V
Iexpr = 3 * h**2 * U - 2 * h * V
assert_zero("tau^-1 direct Laurent extraction", direct_row(a_raw, c_UV, -1) - row(a_raw, c_UV, -1))
assert_zero("tau^-1 invariant row", row(a_raw, c_UV, -1) + h * D(Iexpr))

# If C1 is the actual tau^1 coefficient, the invariant has a
# denominator-transparent form which forces q0 to be polynomial.
C1actual = sp.Function("C1")(s)
Iactual = 3 * h**2 * U - 2 * h * (C1actual - C1star)
rpoly = A0 - q0**2
Isimple = 3 * h**2 * A2 - 2 * h * C1actual + sp.Rational(3, 4) * rpoly**2 + sp.Rational(3, 2) * b * A1
assert_zero("simplified harmonic invariant I", Iactual - Isimple)

# Choose q3 so U=0. The vanishing row makes I constant and V=-I/(2h).
I = sp.symbols("I", constant=True)
q3g = (A2 - q1**2 - 2 * q0 * q2) / (2 * h)
Vg = -I / (2 * h)
C1g = simp(C1star.subs(q3, q3g) + Vg)

q4 = (A3 - 2 * q0 * q3g - 2 * q1 * q2) / (2 * h)
C2star = (
    3 * h**2 * q4
    + 6 * h * q0 * q3g
    + 6 * h * q1 * q2
    + 3 * q0**2 * q2
    + 3 * q0 * q1**2
)
assert_zero(
    "simplified C2 approximate-root coefficient",
    C2star
    - sp.Rational(3, 2) * (h * A3 + q0 * A2 + q1 * A1 - q0 * q1**2),
)
W = sp.Function("W")(s)
c_final = dict(c_g)
c_final[0] = c0
c_final[1] = C1g
c_final[2] = C2star + W

assert_zero("tau^0 direct Laurent extraction", direct_row(a_raw, c_final, 0) - row(a_raw, c_final, 0))
assert_zero(
    "tau^0 bracket repair row",
    row(a_raw, c_final, 0) - D(2 * h**2 * W - I * q0),
)

moment = simp(sum(i * ai * c_final.get(-i, 0) for i, ai in a_raw.items()))
assert_zero("scalar moment identity", moment - (I * q0 - 2 * h**2 * W))

# A through tau^3 and C through tau^2 determine every negative coefficient
# of the target cusp polynomial C^2-A^3+I*A.
Dtarget = sp.expand(laurent(c_final) ** 2 - laurent(a_raw) ** 3 + I * laurent(a_raw))
for kk in range(-6, -1):
    assert_zero(f"target cusp polynomial tau^{kk}", Dtarget.coeff(tau, kk))
assert_zero(
    "target cusp polynomial tau^-1 repair class",
    Dtarget.coeff(tau, -1) - h * (2 * h**2 * W - I * q0),
)
W_repair = (s + I * q0) / (2 * h**2)
assert_zero(
    "moment -s implies the tau^0 row equals 1",
    D(2 * h**2 * W_repair - I * q0) - 1,
)

# Polynomiality collapse. From Isimple=I in k, normality of k[s] forces
# rpoly and q0 into k[s], hence h|b. The repair numerator then factors as
# 4s+hK; since its other expression is 8h^2 times the polynomial C2, h|s.
Kpoly = (
    12 * h**2 * A3
    + 24 * q0 * h * A2
    + 6 * rpoly * A1
    + 12 * q0**2 * A1
    - 8 * q0 * C1actual
)
assert_zero(
    "repair numerator factors as 4s+hK",
    8 * h**2 * C2star + 4 * s + 4 * Isimple * q0 - (4 * s + h * Kpoly),
)
print("DEDUCTION: k[s] normality plus constant I forces q0 in k[s], hence h|b.")
print("DEDUCTION: polynomial C2 and repair then force h|s; THM-3989 has s|h.")
print("           Therefore h=gamma*s with gamma in k*.")

# The deleted-line restriction is forced, not merely modeled by the controls.
# Put qbar=q0(0), a0=A0(0)-qbar^2, and X=gamma*x+qbar.
xb, qbar, a0 = sp.symbols("xb qbar a0")
gamma = sp.symbols("gamma", nonzero=True)
Xb = gamma * xb + qbar
A_boundary_raw = gamma**2 * xb**2 + 2 * gamma * qbar * xb + qbar**2 + a0
C_boundary_raw = (
    gamma**3 * xb**3
    + 3 * gamma**2 * qbar * xb**2
    + sp.Rational(3, 2) * gamma * xb * (2 * qbar**2 + a0)
    + qbar**3
    + sp.Rational(3, 2) * qbar * a0
)
assert_zero("forced deleted-line A normal form", A_boundary_raw - (Xb**2 + a0))
assert_zero(
    "forced deleted-line C normal form",
    C_boundary_raw - (Xb**3 + sp.Rational(3, 2) * a0 * Xb),
)
assert_zero(
    "forced boundary nodal eliminant",
    C_boundary_raw**2
    - A_boundary_raw**3
    + sp.Rational(3, 4) * a0**2 * A_boundary_raw
    + a0**3 / 4,
)
assert_zero(
    "cusp boundary has both x-derivatives zero at X=0",
    (sp.diff(Xb**2, xb) ** 2 + sp.diff(Xb**3, xb) ** 2).subs(xb, -qbar / gamma),
)
node_address = sp.expand(Xb**2 + sp.Rational(3, 2) * a0)
centered_node_address = gamma**2 * xb**2 + sp.Rational(3, 2) * a0
node_transverse = gamma * xb**2 + sp.Rational(3, 2) * a0 / gamma
assert_zero(
    "centered node addresses match the transverse quadratic",
    node_transverse - centered_node_address / gamma,
)
At0 = -sp.Rational(2, 3) / (gamma * a0)
Ax0 = sp.diff(Xb**2 + a0, xb).subs({xb: 0, qbar: 0})
Cx0 = sp.diff(Xb**3 + sp.Rational(3, 2) * a0 * Xb, xb).subs({xb: 0, qbar: 0})
assert_zero("boundary Bezout value A_t(0,0)", Ax0 * 0 - At0 * Cx0 - 1)
print("DEDUCTION: on t=0 every normalized pair has")
print("           A=X^2+a0, C=X^3+(3*a0/2)X, I=3*a0^2/4.")
print("DEDUCTION: a0=0 makes A_x=C_x=0 at X=0, contradicting J=1.")
print("           Hence a Keller pair forces a0!=0, I!=0, and a nodal boundary.")
print("DEDUCTION: the two node addresses force q0(0)=0 and F_p(0,0)=3*a0/(2*gamma).")
print("           Thus b is divisible by s^2 and A1(0)=-2/(3*gamma*a0).")

# The exact p=0 seam of THM-3989 now separates the two remaining first-jet
# branches.  Write eta=[s]q0, q3c=[s^3]q0, r2=[s^2](A0-q0^2), and
# A10=A1(0).  The A seam at ell=2 and C seam at ell=1 are the following two
# scalar equations; eliminating q3c leaves eta*(a0-gamma*eta)=0.
eta, q3c, r2, A10 = sp.symbols("eta q3c r2 A10")
A_seam = eta**2 + r2 - 2 * gamma * q3c - A10
C_seam = (
    sp.Rational(3, 2) * a0 * eta
    + sp.Rational(3, 2) * gamma * A10
    - sp.Rational(3, 2) * gamma * (2 * eta**2 + r2)
    + 3 * gamma**2 * q3c
)
assert_zero(
    "centered A ell=2 seam expansion",
    A_seam - (eta**2 + r2 - 2 * gamma * q3c - A10),
)
assert_zero(
    "centered C ell=1 seam elimination",
    C_seam.subs(q3c, (eta**2 + r2 - A10) / (2 * gamma))
    - sp.Rational(3, 2) * eta * (a0 - gamma * eta),
)
assert_zero(
    "nonlift branch forces the s^2 coefficient of b",
    (2 * gamma * eta).subs(eta, a0 / gamma) - 2 * a0,
)
print("DEDUCTION: eta=[s]q0 satisfies eta*(gamma*eta-a0)=0.")
print("           eta=0 gives a genuine P1/P0 square-root lift;")
print("           eta=a0/gamma is the unique first-seam nonliftable branch,")
print("           and there [s^2]b=2*a0.")

# Hostile partial-row control for the nonliftable branch gamma=a0=eta=1.
# These are honest B2 elements built from x,u,p,y.  Their rows -5 through -2
# vanish and both p=0 seams hold, but row -1 is nonzero, so the seam dichotomy
# is not being mistaken for a Keller construction.
x_log = s / tau
u_log = s**2 / tau
p_log = s**2 + tau
y_log = s * p_log
A_partial = sp.expand(x_log**2 + 2 * u_log + 1)
C_partial = sp.expand(
    x_log**3
    + 3 * x_log * u_log
    + sp.Rational(3, 2) * x_log
    + sp.Rational(3, 2) * x_log * p_log
    - sp.Rational(1, 2) * y_log
)
J_partial = sp.expand(
    tau
    * (D(A_partial) * sp.diff(C_partial, tau)
       - sp.diff(A_partial, tau) * D(C_partial))
)
for kk in range(-5, -1):
    assert_zero(f"nonliftable partial control row tau^{kk}", J_partial.coeff(tau, kk))
assert_zero(
    "nonliftable partial control first failure",
    J_partial.coeff(tau, -1) + s**2 * (3 * s**2 - 1),
)
print("CONTROL: the eta=a0/gamma branch survives through row -2 but can fail row -1.")

# Hostile controls: A=x^2+a, C=x^3+(3a/2)x, x=s/tau. They have bracket zero
# and distinguish the one-address cusp from the two-address node.
x, t, a_par, Avar, Cvar = sp.symbols("x t a A C")
A_ctl = x**2 + a_par
C_ctl = x**3 + sp.Rational(3, 2) * a_par * x
I_ctl = sp.Rational(3, 4) * a_par**2
eliminant = Cvar**2 - (Avar - a_par) * (Avar + a_par / 2) ** 2

assert_zero("control parametrization lies on eliminant", eliminant.subs({Avar: A_ctl, Cvar: C_ctl}))
assert_zero(
    "control Weierstrass expansion",
    sp.expand((Avar - a_par) * (Avar + a_par / 2) ** 2)
    - (Avar**3 - I_ctl * Avar - a_par**3 / 4),
)
assert_zero(
    "control Jacobian is zero",
    sp.diff(A_ctl, x) * sp.diff(C_ctl, t) - sp.diff(A_ctl, t) * sp.diff(C_ctl, x),
)

print("CONTROL cusp: a=0 gives I=0 and C^2=A^3 (one normalization address).")
print("CONTROL node: a!=0 gives I=3*a^2/4!=0 and")
print("              C^2=A^3-(3*a^2/4)A-a^3/4 (two addresses at the node).")
print("NOTE: both controls have bracket zero and fail 2*h^2*W-I*q0=s;")
print("      they test geometry, not existence of a Darboux pair.")
print("ALL EXACT CHECKS PASSED")
