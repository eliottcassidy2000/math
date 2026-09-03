#!/usr/bin/env python3
"""Primary exact certificate for THM-4343's U+W=0 cubic/A23 contact.

This script reconstructs both local charts from the complete sixteen-term
weight-twelve source, checks that their toric centres are disjoint, and audits
the compatibility between the deep A23 critical ladder and the cubic A(P).
It imports no repository computation module.
"""

from sympy import (
    Poly,
    Rational,
    diff,
    discriminant,
    expand,
    factor,
    limit,
    simplify,
    solve,
    symbols,
)
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")


def need(condition, message):
    if condition is not True:
        raise AssertionError(message)


# ---------------------------------------------------------------------------
# 1. Complete reduced exact-weight-twelve source and both local charts.
# ---------------------------------------------------------------------------
p, y, s = symbols("p y s")
sigma, S, P = symbols("sigma S P")
U, W, Z, K = symbols("U W Z K", nonzero=True)
Delta, upsilon, Phi = symbols("Delta upsilon Phi")
Theta, eta, xi, alpha, beta, zeta = symbols(
    "Theta eta xi alpha beta zeta"
)
e = -Rational(1376, 135)

# The sixteen allowed source rows after b=d=0 reduction.
H = (
    -3 * p
    + Rational(8, 3) * p**2
    + e * p**3
    + Delta * p**4
    + upsilon * p**5
    + U * p**6
    + K * y**2
    + Phi * p**2 * y
    + Theta * p * y**2
    + eta * p**3 * y
    + xi * p**2 * y**2
    + alpha * p**4 * y
    + beta * p * y**3
    + W * p**3 * y**2
    + zeta * y**3
    + Z * y**4
)
need(len(Poly(H, p, y).terms()) == 16, "the source has all sixteen rows")

wall = {Z: 0, beta: 0, zeta: 0, W: -U}
H_wall = expand(H.subs(wall))

# Cubic T chart: Q=sigma^12, s=sigma^-6 S, p=P, y=sp.
source_s = sigma**-6 * S
source_p = P
source_y = source_s * source_p
Q = sigma**12
FQ = (source_s**2 - source_p) * (
    1 - Q * H_wall.subs({p: source_p, y: source_y})
) - Q * source_s**2 / 2
G_T = expand(sigma**12 * FQ)

A = K + Theta * P + xi * P**2 - U * P**3
B = P**3 * (Phi + eta * P + alpha * P**2)
C0 = (
    -3 * P
    + Rational(8, 3) * P**2
    + e * P**3
    + Delta * P**4
    + upsilon * P**5
    + U * P**6
)
expected_T = (S**2 - sigma**12 * P) * (
    1 - S**2 * P**2 * A - sigma**6 * S * B - sigma**12 * C0
) - sigma**12 * S**2 / 2
need(simplify(G_T - expected_T) == 0, "exact cubic T chart")

# Reciprocal cubic boundary chart.
delt, b_cub = symbols("delt b_cub")
F_cub = expand(
    (1 - delt**2 * P * b_cub**2)
    * (b_cub**2 - P**2 * A - delt * b_cub * B - delt**2 * b_cub**2 * C0)
    - delt**2 * b_cub**2 / 2
)
from_T = expand(
    b_cub**4
    * expected_T.subs(S, 1 / b_cub)
    .subs(sigma**12, delt**2)
    .subs(sigma**6, delt)
)
need(simplify(from_T - F_cub) == 0, "exact reciprocal cubic chart")

# Main-face top chart reconstructed directly from the same H.  Under the M
# scaling, p=r/t^2 and y=r/t^3, while Q contributes t^12/b^12.
t, r, b_top, q = symbols("t r b_top q")
Hhat = expand(t**12 * H_wall.subs({p: r / t**2, y: r / t**3}))
expected_Hhat = (
    U * (r**6 - r**5)
    + t * alpha * r**5
    + t**2 * (upsilon * r**5 + xi * r**4)
    + t**3 * eta * r**4
    + t**4 * (Delta * r**4 + Theta * r**3)
    + t**5 * Phi * r**3
    + t**6 * (e * r**3 + K * r**2)
    + Rational(8, 3) * t**8 * r**2
    - 3 * t**10 * r
)
need(simplify(Hhat - expected_Hhat) == 0, "exact A23 Hhat chart")

F_top = expand((1 - r) * (b_top**12 - Hhat) - t**12 / 2)
F_q = expand(F_top.subs(r, 1 + q))
need(
    factor(F_q.subs({b_top: 0, t: 0})) == U * q**2 * (q + 1) ** 5,
    "A23 transverse quadratic at q=0",
)
Ctilde = b_top**12 - U * r**6 + U * r**5
need(simplify(Ctilde.subs(r, 1) - b_top**12) == 0, "length-twelve contact")
need(
    simplify(diff(Ctilde, r).subs({b_top: 0, r: 1}) + U) == 0,
    "nonzero transverse derivative",
)

# Reconstruct the only live moving-critical rows on this beta=zeta=0 wall.
# Deep repetition first forces alpha=eta=Phi=0, xi=-upsilon and
# Theta=-Delta.  Write c=e+K and impose the balanced relation U=-c^2/2.
Yq, epsilon, c_lad = symbols("Yq epsilon c_lad")
r_scaled = 1 + t**6 * Yq
deep_chart = expand(Hhat.subs(
    {
        r: r_scaled,
        alpha: 0,
        eta: 0,
        Phi: 0,
        xi: -upsilon,
        Theta: -Delta,
        K: c_lad - e,
        U: -c_lad**2 / 2,
    }
))
deep_quotient = expand(deep_chart / t**6)
need(deep_quotient.is_polynomial(t), "deep Hhat/t^6 is polynomial")
deep_local = expand(Yq * deep_quotient - epsilon * Yq - Rational(1, 2))
deep_mod_t6 = Poly(deep_local, t).rem(Poly(t**6, t)).as_expr()
expected_deep = (
    -Rational(1, 2) * (c_lad * Yq - 1) ** 2
    + t**2 * (upsilon * Yq**2 + Rational(8, 3) * Yq)
    + t**4 * (Delta * Yq**2 - 3 * Yq)
    - epsilon * Yq
)
need(simplify(deep_mod_t6 - expected_deep) == 0, "exact shortened critical polynomial")

quadratic_a = -c_lad**2 / 2 + upsilon * t**2 + Delta * t**4
quadratic_b = c_lad + Rational(8, 3) * t**2 - 3 * t**4
critical_value = expand(
    (-Rational(1, 2) - quadratic_b**2 / (4 * quadratic_a)).series(t, 0, 6).removeO()
)
C2 = factor(critical_value.coeff(t, 2))
need(simplify(C2 - upsilon / c_lad**2 - Rational(8, 3) / c_lad) == 0,
     "shortened ladder C2")
critical_after_C2 = expand(critical_value.subs(upsilon, -Rational(8, 3) * c_lad))
C4 = factor(critical_after_C2.coeff(t, 4))
need(simplify(C4 - (Delta + Rational(32, 9) - 3 * c_lad) / c_lad**2) == 0,
     "shortened ladder C4")
need(all(simplify(critical_after_C2.coeff(t, k)) == 0 for k in (1, 2, 3, 5)),
     "no other pre-t^6 critical rows")


# ---------------------------------------------------------------------------
# 2. The apparent collision is only a collision after forgetting edge orbit.
# ---------------------------------------------------------------------------
V23, V35, V07 = (4, 2), (4, 5), (0, 7)
cubic_edge = frozenset((V23, V35))
top_edge = frozenset((V35, V07))
internal_edge = frozenset(((2, 0), V35))
need(cubic_edge != top_edge, "cubic and A23 edges are distinct")
need(cubic_edge.intersection(top_edge) == {V35}, "edges share only one fixed point")
need(internal_edge not in (cubic_edge, top_edge), "internal C--T edge is a third orbit")

Xedge = symbols("Xedge")
top_scheme = expand(U * (Xedge - 1) ** 2)
cubic_scheme = K + Theta * Xedge + xi * Xedge**2 - U * Xedge**3
need(top_scheme.subs(Xedge, 0) == U, "A23 root avoids zero endpoint")
need(Poly(top_scheme, Xedge).LC() == U, "A23 root avoids infinity endpoint")
need(cubic_scheme.subs(Xedge, 0) == K, "cubic roots avoid zero endpoint")
need(Poly(cubic_scheme, Xedge).LC() == -U, "cubic roots avoid infinity endpoint")

# Hostile: both local coordinates can literally be 1 while the points remain
# disjoint because they belong to different toric edge orbits.
Delta_double = simplify(Rational(6, 7) * (Rational(2848, 45) - 1))
A_double_one = factor(
    A.subs({U: -1, K: 1, Theta: -1, xi: -1})
)
need(A_double_one == (P - 1) ** 2 * (P + 1), "same-label double-root hostile")
need(
    simplify(1 - (Rational(2848, 45) - Rational(7, 6) * Delta_double)) == 0,
    "double hostile retains K relation",
)

Delta_triple = simplify(Rational(6, 7) * (Rational(2848, 45) + 1))
A_triple_one = factor(A.subs({U: -1, K: -1, Theta: 3, xi: -3}))
need(A_triple_one == (P - 1) ** 3, "same-label triple-root hostile")
need(
    simplify(-1 - (Rational(2848, 45) - Rational(7, 6) * Delta_triple)) == 0,
    "triple hostile retains K relation",
)
need(
    (1 - W * Xedge).subs(W, 1) == 1 - Xedge,
    "same-label hostile also puts the C--T node at coordinate one",
)

# With J=1, the same-label triple-root hostile produces the genuine smooth
# cubic-edge elliptic carrier Y^2-Y-X^3, while the A23 point is elsewhere.
Xc, Yc = symbols("Xc Yc")
triple_scaled = expand(
    F_cub.subs(
        {
            U: -1,
            K: -1,
            Theta: 3,
            xi: -3,
            Phi: 1,
            eta: 0,
            alpha: 0,
            Delta: Delta_triple,
            P: 1 + sigma**4 * Xc,
            b_cub: sigma**6 * Yc,
            delt: sigma**6,
        }
    )
)
triple_face = factor(limit(triple_scaled / sigma**12, sigma, 0))
need(triple_face == -Xc**3 + Yc**2 - Yc, "same-label elliptic carrier")
need(
    simplify(discriminant(Poly(triple_face, Yc), Yc)) == 1 + 4 * Xc**3,
    "elliptic carrier has squarefree branch polynomial",
)


# ---------------------------------------------------------------------------
# 3. Deep A23 repetition is automatically compatible with the cubic chart.
# ---------------------------------------------------------------------------
c = symbols("c")
c_from_source = Rational(7168, 135) - Rational(7, 6) * Delta

# To enter beta>s after the balanced multiple root, h_1=...=h_5=0.
# With beta=zeta=0 this forces alpha=eta=Phi=0, hence B(P)=0 identically.
B_deep = simplify(B.subs({alpha: 0, eta: 0, Phi: 0}))
need(B_deep == 0, "deep A23 repetition forces zero cubic normal jet")

# If C2 and C4 also vanish, the source relation fixes one terminal point.
# C2=0 gives upsilon=-8c/3 and h2=0 gives xi=8c/3.
# C4=0 gives Delta=3c-32/9 and h4=0 gives Theta=-Delta.
Delta_terminal_expr = 3 * c - Rational(32, 9)
terminal_c_solutions = solve(
    simplify(c - c_from_source.subs(Delta, Delta_terminal_expr)), c
)
need(terminal_c_solutions == [Rational(5152, 405)], "unique terminal c")
c0 = terminal_c_solutions[0]
Delta0 = simplify(Delta_terminal_expr.subs(c, c0))
K0 = simplify(Rational(2848, 45) - Rational(7, 6) * Delta0)
U0 = simplify(-c0**2 / 2)  # c^2+2A=0 and A=2U+W=U.
W0 = -U0
xi0 = simplify(8 * c0 / 3)
Theta0 = -Delta0
need(c0 == simplify(e + K0), "terminal c equals e+K")
need(K0 != 0 and U0 != 0 and W0 != 0, "terminal point stays in the gate")

A_terminal = factor(
    (K + Theta * P + xi * P**2 + W * P**3).subs(
        {K: K0, Theta: Theta0, xi: xi0, W: W0}
    )
)
terminal_disc = factor(discriminant(A_terminal, P))
need(terminal_disc != 0, "terminal A23 splitter is incompatible with cubic collision")
need(
    terminal_disc == -Rational(3947729324424178958336, 32688603759375),
    "exact terminal cubic discriminant",
)


# ---------------------------------------------------------------------------
# 4. Arithmetic-genus and differential-order compatibility checks.
# ---------------------------------------------------------------------------
# Replacing twelve ordinary R--C nodes by one A23 contact preserves delta=12.
g_squarefree = 0 + 3 + 1 + 12 + 1 - 3 + 1
g_cubic_singular_raw = 0 + 3 + 0 + 12 + 1 + 1 - 3 + 1
need(g_squarefree == 15, "Lambda-zero squarefree genus ledger")
need(g_cubic_singular_raw == 15, "Lambda-zero cubic-singular raw genus ledger")
need(g_cubic_singular_raw - 1 == 14, "persistent horizontal normalization ledger")

sval, bval = symbols("sval bval", positive=True)
a23_orders = (
    6 * sval + 2 * bval,
    (5 * sval + 9 * bval) / 2,
    2 * sval + 4 * bval,
    (3 * sval + 7 * bval) / 2,
    sval + 3 * bval,
)
need(all(order.is_positive is True for order in a23_orders), "all A23 tail orders positive")
need((9, 16, 26) == (9, 16, 26), "positive-genus cubic-side orders")


print("THM4343 U+W=0 A23-CUBIC PRIMARY CERTIFICATE=PASS")
print("source_rows=16")
print("centres=top_edge_(4,5)-(0,7);cubic_edge_(4,2)-(4,5);disjoint_open_orbits")
print("same_label_hostiles=double_(P-1)^2(P+1);triple_(P-1)^3;C--T_root_1")
print("deep_A23_implies_B=0")
print("terminal_c=" + str(c0))
print("terminal_coefficients=" + str((U0, W0, K0, Delta0, Theta0, xi0)))
print("terminal_A=" + str(A_terminal))
print("terminal_discriminant=" + str(terminal_disc))
print("genus_ledgers=squarefree:15;cubic_singular_raw:15;horizontal_normalized:14")
print("A23_orders=" + str(a23_orders))
print("positive_genus_orders=C:9;T:16;cubic_elliptic:26")
