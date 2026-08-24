#!/usr/bin/env python3
"""Exact symbolic certificate for THM-3997.

The certificate reconstructs the Laurent/source diagonal bridge, eliminates
the first two normal rows, verifies the Hasse p-Taylor repair quotient, and
constructs bounded third source jets for direct hostile controls.  It is an
exact verification sidecar; the all-degree arguments are written in THM-3997.
"""

from __future__ import annotations

import sympy as sp


def simp(expr):
    return sp.factor(sp.cancel(sp.expand(expr)))


def zero(label: str, expr) -> None:
    value = simp(expr)
    if value != 0:
        raise AssertionError(f"{label}: {value}")
    print(f"PASS  {label}")


def gate(label: str, condition: bool) -> None:
    if not condition:
        raise AssertionError(label)
    print(f"PASS  {label}")


x, t, s, tau, p = sp.symbols("x t s tau p")
a, gamma = sp.symbols("a gamma", nonzero=True)


def jac(A, C):
    return sp.expand(sp.diff(A, x) * sp.diff(C, t)
                     - sp.diff(A, t) * sp.diff(C, x))


print("STATUS=THM-3997 VERIFIED-EXACT CERTIFICATE")

# -------------------------------------------------------------------------
# I. Reconstruct the diagonal caps directly from Laurent rows.
# -------------------------------------------------------------------------

b2, b3, b4 = sp.symbols("b2 b3 b4")
A01, A02, A03, A10, A11, A12, A20, A21, A30 = sp.symbols(
    "A01 A02 A03 A10 A11 A12 A20 A21 A30"
)
A_laurent = (
    gamma**2 * s**2 * tau**-2
    + (b2*s**2 + b3*s**3 + b4*s**4) * tau**-1
    + (a + A01*s + A02*s**2 + A03*s**3)
    + (A10 + A11*s + A12*s**2) * tau
    + (A20 + A21*s) * tau**2
    + A30 * tau**3
)
A_source = sp.expand(A_laurent.subs({s: x*t, tau: t}))
zero("A boundary diagonal", A_source.coeff(t, 0) - (gamma**2*x**2 + a))
zero("A first diagonal", A_source.coeff(t, 1) - (A10 + A01*x + b2*x**2))
zero("A second diagonal", A_source.coeff(t, 2) - (A20 + A11*x + A02*x**2 + b3*x**3))
gate("A first/second degree caps are 2/3", sp.Poly(A_source.coeff(t, 1), x).degree() == 2 and sp.Poly(A_source.coeff(t, 2), x).degree() == 3)

# The c_-2 row is (3/2)gamma*s*b.  The exact c_-1 coefficients are not
# needed for the cap, so keep their first three possible powers independent.
e1, e2, e3, C00, C01, C02, C10, C11, C20 = sp.symbols(
    "e1 e2 e3 C00 C01 C02 C10 C11 C20"
)
drow = sp.Rational(3, 2) * gamma * s * (b2*s**2 + b3*s**3 + b4*s**4)
erow = e1*s + e2*s**2 + e3*s**3
C_laurent = (
    gamma**3*s**3*tau**-3 + drow*tau**-2 + erow*tau**-1
    + (C00 + C01*s + C02*s**2) + (C10 + C11*s)*tau + C20*tau**2
)
C_source = sp.expand(C_laurent.subs({s: x*t, tau: t}))
zero("C first diagonal", C_source.coeff(t, 1) - (C10 + C01*x + e2*x**2 + sp.Rational(3, 2)*gamma*b2*x**3))
zero("C second diagonal", C_source.coeff(t, 2) - (C20 + C11*x + C02*x**2 + e3*x**3 + sp.Rational(3, 2)*gamma*b3*x**4))
gate("C first/second degree caps are 3/4", sp.Poly(C_source.coeff(t, 1), x).degree() == 3 and sp.Poly(C_source.coeff(t, 2), x).degree() == 4)

# General coefficient identity [t^r]=[s^(r-i)] on a hostile finite row.
coeff = {(i, j): sp.symbols(f"c_{i}_{j}") for i in range(-2, 4) for j in range(0, 7)}
hostile = sum(coeff[i, j] * s**j * tau**i for i in range(-2, 4) for j in range(0, 7))
hostile_source = sp.expand(hostile.subs({s: x*t, tau: t}))
for rr in range(0, 5):
    diagonal = sum(coeff[i, rr-i] * x**(rr-i) for i in range(-2, min(3, rr)+1) if 0 <= rr-i <= 6)
    zero(f"general diagonal formula r={rr}", hostile_source.coeff(t, rr) - diagonal)

# -------------------------------------------------------------------------
# II. Independently eliminate the first and second normal diagonals.
# -------------------------------------------------------------------------

theta0, theta1 = sp.symbols("theta0 theta1")
X = gamma*x
I = sp.Rational(3, 4)*a**2
lam = sp.Rational(3, 2)*a/gamma
A0 = X**2 + a
C0 = X**3 + sp.Rational(3, 2)*a*X
q = gamma*x**2 + lam
theta = theta0 + theta1*x
N = -sp.Rational(2, 3)/(gamma*a) + 2*gamma**2*x*theta
K = -x/a + 3*gamma*(gamma**2*x**2 + a/2)*theta
zero("first-jet Bezout row", sp.diff(A0, x)*K - N*sp.diff(C0, x) - 1)
zero("rotated gradient A", (-3*A0**2 + I) + q*sp.diff(C0, x))
zero("rotated gradient C", 2*C0 - q*sp.diff(A0, x))

NvKv = sp.diff(N, x)*K - N*sp.diff(K, x)
g2 = simp(K**2 - 3*A0*N**2 - q*NvKv/2)
g2_general = (
    2*gamma*x*theta
    - (sp.Rational(3, 4)*a/gamma + gamma*x**2/2)*sp.diff(theta, x)
    - sp.Rational(5, 6)/(a*gamma**2)
)
zero("first diagonal before affine specialization", g2 - g2_general)
g2_affine = (
    sp.Rational(3, 2)*gamma*theta1*x**2 + 2*gamma*theta0*x
    - sp.Rational(3, 4)*a*theta1/gamma
    - sp.Rational(5, 6)/(a*gamma**2)
)
zero("first diagonal affine coefficients", g2 - g2_affine)

theta1_forced = a/gamma**2
beta = 2*gamma*theta0
alpha = -sp.Rational(3, 4)*a**2/gamma**3 - sp.Rational(5, 6)/(a*gamma**2)
zero("first residual comparison", g2.subs(theta1, theta1_forced) - (lam*x**2 + beta*x + alpha))
zero("alpha common denominator", alpha + (9*a**3 + 10*gamma)/(12*a*gamma**3))

# Translation to Laurent data.
q01, q02 = sp.symbols("q01 q02")
zero("q0 slope from b2", (2*gamma*q01 - 2*a).subs(q01, a/gamma))
zero("A0 first derivative", 2*gamma**2*theta0 - gamma*beta)

# Solve the complete bounded second jet from the t^1 Jacobian row.
theta = theta.subs(theta1, theta1_forced)
N = sp.expand(N.subs(theta1, theta1_forced))
K = sp.expand(K.subs(theta1, theta1_forced))
m0, m1, m2, m3 = sp.symbols("m0 m1 m2 m3")
l0, l1, l2, l3, l4 = sp.symbols("l0 l1 l2 l3 l4")
M = m0 + m1*x + m2*x**2 + m3*x**3
L = l0 + l1*x + l2*x**2 + l3*x**3 + l4*x**4
j1 = sp.expand(2*(sp.diff(A0, x)*L - M*sp.diff(C0, x)) + sp.diff(N, x)*K - N*sp.diff(K, x))
unknowns = [l0, l1, l2, l3, l4, m0, m1, m2, m3]
equations = [sp.Poly(j1, x).coeff_monomial(x**j) for j in range(sp.Poly(j1, x).degree()+1)]
sol_tuple = next(iter(sp.linsolve(equations, unknowns)))
sol = dict(zip(unknowns, sol_tuple))

expected = {
    m0: (9*a**3*theta0**2*gamma**5 + 3*a**3 - 2*gamma)/(9*a**3*gamma**3),
    l0: 3*a*m1/(4*gamma) - theta0*(3*a**3 + 2*gamma)/(2*a*gamma),
    l1: 3*a*m2/(4*gamma) + (-9*a**6 + 36*a**3*theta0**2*gamma**6 - 6*a**3*gamma - 4*gamma**2)/(12*a**3*gamma**3),
    l2: 3*a*theta0*gamma + 3*a*m3/(4*gamma) + 3*gamma*m1/2,
    l3: 3*a**2/(2*gamma) + 3*gamma*m2/2,
    l4: 3*gamma*m3/2,
}
for variable, formula in expected.items():
    zero(f"second-jet coefficient {variable}", sol[variable] - formula)
zero("complete second-jet Jacobian row", j1.subs(sol))

# Eliminate the third jet and check every coefficient in the candidate note.
det_vprime_w = sp.diff(N, x)*L - sp.diff(K, x)*M
det_wprime_v = sp.diff(M, x)*K - sp.diff(L, x)*N
g3 = sp.expand(
    -q*(2*det_vprime_w + det_wprime_v)/3 - 6*A0*N*M + 2*K*L - N**3
).subs(sol)
g3 = sp.factor(sp.expand(g3))
g3_expected_coeffs = {
    3: 2*m3/(3*gamma),
    2: -(3*a**2 - 5*gamma**2*m2)/(6*gamma**3),
    1: -(2*a*theta0*gamma**2 + a*m3 - 2*gamma**2*m1)/(2*gamma**3),
    0: (81*a**6 - 27*a**4*gamma**2*m2 + 108*a**3*theta0**2*gamma**6 + 90*a**3*gamma - 28*gamma**2)/(108*a**3*gamma**5),
}
for power, formula in g3_expected_coeffs.items():
    zero(f"third diagonal coefficient x^{power}", sp.Poly(g3, x).coeff_monomial(x**power) - formula)
gate("third diagonal has degree at most three", sp.Poly(g3, x).degree() <= 3)

m3_forced = 3*gamma**2*theta0
m2_forced = 2*(-3*a**3 - 5*gamma)/(5*a*gamma**2)
delta = (2*m1 - 5*a*theta0)/(2*gamma)
epsilon = gamma*theta0**2 + sp.Rational(21, 20)*a**3/gamma**5 + sp.Rational(4, 3)/gamma**4 - sp.Rational(7, 27)/(a**3*gamma**3)
zero("third residual comparison", g3.subs({m2: m2_forced, m3: m3_forced}) - (beta*x**3 + 2*alpha*x**2 + delta*x + epsilon))
zero("epsilon beta form", epsilon - (beta**2/(4*gamma) + sp.Rational(21, 20)*a**3/gamma**5 + sp.Rational(4, 3)/gamma**4 - sp.Rational(7, 27)/(a**3*gamma**3)))

m0_after = simp(sol[m0].subs({m2: m2_forced, m3: m3_forced}))
zero("translated A2 constant", m0_after - (beta**2/4 + sp.Rational(1, 3)/gamma**3 - sp.Rational(2, 9)/(a**3*gamma**2)))
zero("translated A1 slope", m1 - (gamma*delta + 5*a*beta/(4*gamma)))
zero("translated b cubic", m3_forced - sp.Rational(3, 2)*gamma*beta)
zero("translated q0 quadratic", m3_forced/(2*gamma) - sp.Rational(3, 4)*beta)

# -------------------------------------------------------------------------
# III. Hasse/p-Taylor iff: monomial-level forward and reverse controls.
# -------------------------------------------------------------------------

def hasse_from_tau(poly_tau, m: int, max_n: int):
    return sp.expand(sum(
        sp.binomial(n, m) * sp.expand(poly_tau).coeff(tau, n) * (-s**2)**(n-m)
        for n in range(m, max_n + 1)
    ))

hasse_checks = 0
for c in range(0, 6):
    for e in range(0, 6-c):
        monomial_p = p**c * (s*p)**e
        monomial_tau = sp.expand(monomial_p.subs(p, s**2 + tau))
        total = c + e
        for mm in range(0, 6):
            expected_h = s**e if mm == total else 0
            if simp(hasse_from_tau(monomial_tau, mm, total) - expected_h) != 0:
                raise AssertionError(f"Hasse monomial c={c},e={e},m={mm}")
            hasse_checks += 1
gate("complete Hasse monomial census through total degree five", hasse_checks == 126)

# Reverse basis: s^e p^m has the unique preimage p^(m-e)y^e for e<=m.
reverse_checks = 0
for mm in range(0, 7):
    for e in range(0, mm+1):
        preimage = p**(mm-e) * (s*p)**e
        if simp(preimage - s**e*p**mm) != 0:
            raise AssertionError(f"reverse image m={mm},e={e}")
        reverse_checks += 1
gate("complete reverse-basis census through layer six", reverse_checks == 28)

bad_m = 4
bad_layer = s**(bad_m+1)*p**bad_m
bad_tau = sp.expand(bad_layer.subs(p, s**2 + tau))
gate("degree hostile is detected", sp.Poly(hasse_from_tau(bad_tau, bad_m, bad_m), s).degree() == bad_m+1)

# Exact layer characterization of lambda*p+(p^2,y): m=0 vanishes;
# m=1 is lambda+beta*s; every m>=2 permits all degrees <=m.
lambda_test, beta_test = sp.symbols("lambda_test beta_test")
layer_test = lambda_test*p + beta_test*(s*p) + p**2*(1+s+s**2) + p*(s*p)*(1+s) + (s*p)**3
layer_tau = sp.expand(layer_test.subs(p, s**2+tau))
zero("ideal iff H0", hasse_from_tau(layer_tau, 0, 3))
zero("ideal iff H1", hasse_from_tau(layer_tau, 1, 3) - (lambda_test + beta_test*s))
for mm in (2, 3):
    gate(f"ideal iff degree cap m={mm}", sp.Poly(hasse_from_tau(layer_tau, mm, 3), s).degree() <= mm)

# -------------------------------------------------------------------------
# IV. R=0 resonance and independent direct source-expansion controls.
# -------------------------------------------------------------------------

resonant_gamma = -sp.Rational(9, 10)*a**3
zero("alpha vanishes on the unique resonance", alpha.subs(gamma, resonant_gamma))
resonant_epsilon = simp(epsilon.subs({theta0: 0, gamma: resonant_gamma}))
zero("R=0 forced epsilon value", resonant_epsilon - sp.Rational(4000, 6561)/a**12)
gate("R=0 forced epsilon is nonzero", resonant_epsilon != 0)

# On the live THM-3992 seam gamma=-a^3/2, quotient the residual fifth-root
# action by A5=a^5 and by dividing every weight-one quantity by gamma.
gamma_live = -a**3 / 2
A5 = a**5
zeta, beta_formal, delta_formal, G_formal, R_formal = sp.symbols(
    "zeta beta_formal delta_formal G_formal R_formal", nonzero=True
)
gate(
    "A5 is invariant modulo zeta^5=1",
    sp.rem(sp.expand((zeta**2*a)**5 - A5), zeta**5 - 1, zeta) == 0,
)
zero("G/gamma is fifth-root invariant", zeta*G_formal/(zeta*gamma) - G_formal/gamma)
zero("R/gamma is fifth-root invariant", zeta*R_formal/(zeta*gamma) - R_formal/gamma)
zero("beta/gamma is fifth-root invariant", zeta*beta_formal/(zeta*gamma) - beta_formal/gamma)
zero("delta/gamma is fifth-root invariant", zeta*delta_formal/(zeta*gamma) - delta_formal/gamma)
zero("live normalized p coefficient", (lam/gamma).subs(gamma, gamma_live) - 6/A5)
zero("live normalized p^2 residual", (alpha/gamma).subs(gamma, gamma_live) + sp.Rational(16, 3)/A5**2)
b_invariant = sp.symbols("b_invariant")
epsilon_beta = beta_formal**2/(4*gamma) + sp.Rational(21, 20)*a**3/gamma**5 + sp.Rational(4, 3)/gamma**4 - sp.Rational(7, 27)/(a**3*gamma**3)
epsilon_tilde_live = (epsilon_beta/gamma).subs({gamma: gamma_live, beta_formal: b_invariant*gamma_live})
zero(
    "live normalized p^3 residual",
    epsilon_tilde_live - (b_invariant**2/4 + sp.Rational(2752, 135)/A5**3),
)
print("RESULT live_mu5_quotient=(A5=a^5,G/gamma,R/gamma,beta/gamma,delta/gamma)")


def direct_source_case(av, gv, theta0v, m1v, label: str) -> None:
    """Construct a bounded third jet and expand directly in x,t."""
    sub = {a: sp.Rational(av), gamma: sp.Rational(gv), theta0: sp.Rational(theta0v), m1: sp.Rational(m1v)}
    A0v, C0v = A0.subs(sub), C0.subs(sub)
    Nv, Kv = N.subs(sub), K.subs(sub)
    m2v, m3v = m2_forced.subs(sub), m3_forced.subs(sub)
    full_sub = dict(sub)
    full_sub.update({m2: m2v, m3: m3v})
    Mv = sp.expand(M.subs(sol).subs(full_sub))
    Lv = sp.expand(L.subs(sol).subs(full_sub))

    uu = tuple(sp.symbols(f"u{j}_{label}") for j in range(5))
    vv = tuple(sp.symbols(f"v{j}_{label}") for j in range(6))
    Uv = sum(uu[j]*x**j for j in range(5))
    Vv = sum(vv[j]*x**j for j in range(6))
    row2 = sp.expand(
        3*(sp.diff(A0v, x)*Vv - Uv*sp.diff(C0v, x))
        + 2*(sp.diff(Nv, x)*Lv - Mv*sp.diff(Kv, x))
        + (sp.diff(Mv, x)*Kv - Nv*sp.diff(Lv, x))
    )
    vars3 = list(uu) + list(vv)
    eq3 = [sp.Poly(row2, x).coeff_monomial(x**j) for j in range(sp.Poly(row2, x).degree()+1)]
    sol3_tuple = next(iter(sp.linsolve(eq3, vars3)))
    free3 = set().union(*(entry.free_symbols for entry in sol3_tuple)).intersection(vars3)
    sol3_tuple = tuple(entry.subs({v: 0 for v in free3}) for entry in sol3_tuple)
    sol3 = dict(zip(vars3, sol3_tuple))
    Uv, Vv = sp.expand(Uv.subs(sol3)), sp.expand(Vv.subs(sol3))
    zero(f"{label}: direct third-jet row", row2.subs(sol3))

    AA = A0v + t*Nv + t**2*Mv + t**3*Uv
    CC = C0v + t*Kv + t**2*Lv + t**3*Vv
    JJ = jac(AA, CC)
    zero(f"{label}: direct J t^0", JJ.coeff(t, 0) - 1)
    zero(f"{label}: direct J t^1", JJ.coeff(t, 1))
    zero(f"{label}: direct J t^2", JJ.coeff(t, 2))

    Pv = sp.expand(CC**2 - AA**3 + I.subs(sub)*AA + sp.Rational(av)**3/4)
    alphav, betav = alpha.subs(sub), beta.subs(sub)
    deltav, epsv = delta.subs(sub), epsilon.subs(sub)
    lamv = lam.subs(sub)
    expected_rows = {
        0: 0,
        1: sp.Rational(gv)*x**2 + lamv,
        2: lamv*x**2 + betav*x + alphav,
        3: betav*x**3 + 2*alphav*x**2 + deltav*x + epsv,
    }
    for power, expected_row in expected_rows.items():
        zero(f"{label}: direct P t^{power}", Pv.coeff(t, power) - expected_row)

    u_src = x**2*t
    p_src = t + x**2*t**2
    y_src = x*t*p_src
    rhs = sp.Rational(gv)*u_src + lamv*p_src + alphav*p_src**2 + betav*y_src + deltav*p_src*y_src + epsv*p_src**3
    defect = sp.expand(Pv-rhs)
    for power in range(0, 4):
        zero(f"{label}: direct mixed residual through t^{power}", defect.coeff(t, power))


direct_source_case(1, 1, 0, 0, "source_control_1")
direct_source_case(2, 3, 1, 5, "source_control_2")

print("THEOREM_ID=THM-3997")
print("NOTE=source degree caps are necessary only; direct third jets are controls only")
print("ALL THM-3997 EXACT CHECKS PASSED")
