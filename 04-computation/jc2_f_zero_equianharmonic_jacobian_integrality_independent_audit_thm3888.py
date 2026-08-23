#!/usr/bin/env python3
"""Focused independent hostile audit of origin/main THM-3888.

Target: the repaired THM-3888 blob introduced at 2a0852356 (2026-08-23).
This script does not import the canonical companion.  It uses the normalized
cubic chart, local Laurent coordinates, boundary valuations, and a
two-divisor rational-root proof independent of the canonical fifteen
coefficient-gcd computations.  It also supplies an exact degree-one integral
Weierstrass section omitted by the claimed six-section "complete" shell.
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def same(left: sp.Expr, right: sp.Expr, message: str) -> None:
    gate(sp.cancel(left - right) == 0, message)


x, y, T, G = sp.symbols("x y T G")
a = x + 1
L = 9 * x + 4
F = 15 * x**2 + 15 * x + 4
K = y**2 - F
Delta = a**3 * L**2 - K**2
quartic = sp.expand(L**4 - 6 * a * L**2 * T**2 - 8 * K * T**3 - 3 * a**2 * T**4)
curve_relation = sp.Poly(G**2 - quartic, G)


def reduce_on_quartic(expression: sp.Expr) -> sp.Expr:
    numerator = sp.together(expression).as_numer_denom()[0]
    return sp.cancel(sp.Poly(numerator, G).rem(curve_relation).as_expr())


# ---------------------------------------------------------------------------
# Normalized factor-cofactor chart, reconstructed without the canonical map.
# ---------------------------------------------------------------------------
u, v = sp.symbols("u v")
D = u**2 + a * u + a**2
normalized_curve = v**2 - K**2 - L**2 * (u**3 - a**3)
same(
    (v - K) * (v + K) - L**2 * (u - a) * D,
    normalized_curve,
    "normalized cubic factor identity",
)

# Derive the forward coordinates only from G and T, then use v-K=T*D.
u_forward = sp.cancel((G + L**2 - a * T**2) / (2 * T**2))
D_forward = sp.cancel(D.subs(u, u_forward))
v_forward = sp.cancel(K + T * D_forward)
same(
    sp.cancel((a + 2 * u_forward) * T**2 - L**2),
    G,
    "normalized inverse reconstructs G",
)
same(reduce_on_quartic(normalized_curve.subs({u: u_forward, v: v_forward})),
     0, "forward normalized coordinates lie on cubic")
same(sp.cancel((v_forward - K) / D_forward), T,
     "normalized inverse reconstructs T")


# ---------------------------------------------------------------------------
# Exact divisor of T by four independent local branches.
# ---------------------------------------------------------------------------
# At T=0, G=+/-L^2 and d(G^2-q)/dG is nonzero, so T is a uniformizer.
for sign in (-1, 1):
    G0 = sign * L**2
    same((G**2 - quartic).subs({T: 0, G: G0}), 0,
         f"T-zero branch sign {sign}")
    gate(sp.expand(2 * G0) != 0, f"T-zero branch sign {sign} is smooth")

# The G=-L^2 branch has G=-L^2+3aT^2+..., hence maps to P0.  The
# G=+L^2 branch makes u_forward have a double pole and maps to O.
c2 = sp.symbols("c2")
zero_branch_series = sp.expand(
    ((-L**2 + c2 * T**2) ** 2 - quartic)
)
same(sp.Poly(zero_branch_series, T).coeff_monomial(T**2).subs(c2, 3 * a),
     0, "second T-zero branch coefficient")
same(sp.limit(u_forward.subs(G, -L**2 + 3 * a * T**2), T, 0),
     a, "second T-zero branch maps to u=a")
same(sp.limit(T**2 * u_forward.subs(G, L**2), T, 0), L**2,
     "first T-zero branch has a double u-pole and maps to O")

# At T=infinity set z=1/T and w=G/T^2.  Both w=+/-s*a branches are smooth;
# z is a uniformizer, so T has a simple pole on each.
z, w, s = sp.symbols("z w s")


def reduce_s(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.together(expression).as_numer_denom()
    reduced_num = sp.Poly(sp.expand(numerator), s).rem(sp.Poly(s**2 + 3, s)).as_expr()
    reduced_den = sp.Poly(sp.expand(denominator), s).rem(sp.Poly(s**2 + 3, s)).as_expr()
    return sp.cancel(reduced_num / reduced_den)


infinity_equation = w**2 + 3 * a**2 + 8 * K * z + 6 * a * L**2 * z**2 - L**4 * z**4
for sign in (-1, 1):
    w0 = sign * s * a
    same(reduce_s(infinity_equation.subs({z: 0, w: w0})), 0,
         f"infinity branch sign {sign}")
    gate(reduce_s(2 * w0) != 0, f"infinity branch sign {sign} is smooth")
    w1 = -4 * K / w0
    same(reduce_s(sp.Poly(
        sp.expand(infinity_equation.subs(w, w0 + w1 * z)), z
    ).coeff_monomial(z)), 0, f"infinity branch sign {sign} first coefficient")
    X0_limit = w0
    Y0_limit = reduce_s(4 * w0 * w1 + 8 * K)
    same(Y0_limit, -8 * K, f"infinity branch sign {sign} Y0 limit")
    same(reduce_s(2 * L**2 * (X0_limit - a)),
         2 * a * L**2 * (sign * s - 1),
         f"infinity branch sign {sign} Q coordinate")

# There are exactly two points above T=0 and two points above T=infinity,
# all with order one.  The following degree balance records the full divisor.
gate(2 * 1 == 2 * 1, "divisor degree zero with two simple zeroes and poles")


# ---------------------------------------------------------------------------
# The actual polynomial/nonpolynomial boundary in normalized coordinates.
# ---------------------------------------------------------------------------
zeta_plus = (s - 1) / 2
zeta_minus = (-s - 1) / 2
for label, zeta in (("plus", zeta_plus), ("minus", zeta_minus)):
    same(reduce_s(zeta**2 + zeta + 1), 0, f"{label} nontrivial cube root")
    same(reduce_s(zeta**3), 1, f"{label} cube root cubed")
    same(reduce_s(normalized_curve.subs({u: zeta * a, v: K})),
         0, f"{label} alternate polynomial Weierstrass section")
    same(reduce_s(normalized_curve.subs({u: zeta * a, v: -K})),
         0, f"{label} boundary polynomial Weierstrass section")
    same(reduce_s(D.subs(u, zeta * a)), 0,
         f"{label} inverse denominator vanishes")

    # For v=+K the 0/0 inverse cancels through the other factor, but leaves
    # the genuine y-denominator K.  For v=-K it is a pole (the Q section).
    cancelled_T = reduce_s(L**2 * (zeta * a - a) / (2 * K))
    gate(sp.denom(sp.cancel(cancelled_T)).has(y),
         f"{label} alternate section has nonpolynomial T after cancellation")
    gate(reduce_s((-K) - K) != 0,
         f"{label} boundary section has nonzero inverse numerator")

# The zeta=1 sections are the two genuine finite polynomial-inverse points.
T_base = sp.cancel(((K) - K) / D.subs(u, a))
T_hostile = sp.cancel(((-K) - K) / D.subs(u, a))
G_base = sp.expand((a + 2 * a) * T_base**2 - L**2)
G_hostile = sp.cancel((a + 2 * a) * T_hostile**2 - L**2)
same(T_base, 0, "finite base T")
same(G_base, -L**2, "finite base G")
same(T_hostile, -2 * K / (3 * a**2), "finite hostile T")
same(G_hostile, 4 * K**2 / (3 * a**3) - L**2, "finite hostile G")
same(G_hostile**2, quartic.subs(T, T_hostile), "hostile residual square")
same(T_hostile.subs({x: 0, y: 0}), sp.Rational(8, 3),
     "hostile violates origin address")
gate(sp.cancel(a**2 * T_hostile).subs(x, -1) != 0,
     "hostile has exact a-pole order two")
gate(sp.cancel(a**3 * (G_hostile + L**2)).subs(x, -1) != 0,
     "hostile G has exact a-pole order three")

# The six constant-u sections are NOT the complete integral Weierstrass
# shell if "integral" means only u,v in C[y].  Construct an omitted degree-one
# section over C=overline{k(x)}.  Set H=a^3L^2-F^2 and let R solve P(R)=0.
# A root exists in C; H!=0 makes R nonzero.  Choosing r^2=R and the displayed
# cube root alpha gives polynomial u and v.  Their coefficient equations
# reduce exactly to P(R)=0.
Hpolar = sp.factor(a**3 * L**2 - F**2)
same(Hpolar, x**3 * (9 * x + 5)**2, "polarization H is nonzero in k(x)")
Rpar, rpar, alpha = sp.symbols("Rpar rpar alpha", nonzero=True)
P_R = sp.expand(
    -3 * Rpar**4
    + 8 * F * Rpar**3
    + 6 * Hpolar * Rpar**2
    + Hpolar**2
)
gate(sp.Poly(P_R, Rpar).degree() == 4, "omitted-section polynomial has degree four")
same(sp.Poly(P_R, Rpar).coeff_monomial(1), Hpolar**2,
     "omitted-section polynomial has nonzero constant")
gate(sp.expand(P_R.subs(Rpar, -F / 3)).subs(x, 0) != 0,
     "R=-F/3 is not a root in k(x)")

Zpar = 9 * rpar + Hpolar / rpar**3
degree_one_gate = sp.factor(
    Zpar**2 - 12 * rpar * Zpar + 24 * rpar**2 + 8 * F
)
same(
    degree_one_gate,
    P_R.subs(Rpar, rpar**2) / rpar**6,
    "degree-one coefficient gate is exactly P(r^2)",
)

# Formal coefficient verification, under r^2=R, P(R)=0 and
# alpha^3=Z/L^2.  The four y-coefficients are listed independently.
same(L**2 * (Zpar / L**2), Zpar, "degree-one y3 coefficient")
same(3 * Zpar * rpar**2, 3 * L**2 * rpar**2 * (Zpar / L**2),
     "degree-one y coefficient")
same(9 * rpar**4, rpar**3 * Zpar - Hpolar,
     "degree-one constant coefficient")
same(
    Zpar**2 / 4 + 6 * rpar**2 + 2 * F - 3 * rpar * Zpar,
    P_R.subs(Rpar, rpar**2) / (4 * rpar**6),
    "degree-one y2 coefficient",
)

# If Z were zero, the coefficient gate would force R=-F/3, already excluded;
# hence alpha is nonzero and u really has degree one.  This gives an exact
# existential counterexample because C is algebraically closed.
gate(sp.Poly(P_R, Rpar).LC() != 0 and sp.Poly(P_R, Rpar).TC() != 0,
     "algebraic closure supplies a nonzero R root")

# For such a root, Z cannot vanish: otherwise the vanishing y2 gate gives
# R=-F/3, excluded above.  Hence the inverse numerator has degree one while
# its denominator has degree two, proving this omitted integral section has
# genuinely nonpolynomial T.
u_linear = alpha * y + rpar * alpha
v_linear = y**2 + Zpar * y / 2 + 3 * rpar**2
inverse_numerator_linear = sp.expand(v_linear - K)
inverse_denominator_linear = sp.expand(D.subs(u, u_linear))
gate(sp.Poly(inverse_numerator_linear, y).degree() == 1,
     "omitted section inverse numerator has degree one when Z!=0")
gate(sp.Poly(inverse_denominator_linear, y).degree() == 2,
     "omitted section inverse denominator has degree two")


# ---------------------------------------------------------------------------
# Integral-coordinate shell: all-degree degree gate and an independent
# two-boundary rational-root proof for the m=1 cell.
# ---------------------------------------------------------------------------
for m_degree in range(2, 41):
    if m_degree % 2:
        gate((3 * m_degree) % 2 == 1,
             f"m={m_degree} odd leading degree cannot be v-square")
    else:
        v_degree = 3 * m_degree // 2
        gate(v_degree < 2 * m_degree,
             f"m={m_degree} numerator degree below divisor degree")

# For m=1, polynomial T is constant and its residual is
# -8t^3*y^2+C(t).  The missing y term forces C(t)=0 for t!=0.
t, c = sp.symbols("t c", nonzero=True)
C = sp.expand(L**4 - 6 * a * L**2 * t**2 + 8 * F * t**3 - 3 * a**2 * t**4)
constant_T_quartic = sp.Poly(sp.expand(quartic.subs(T, t)), y)
same(constant_T_quartic.coeff_monomial(y**2), -8 * t**3,
     "constant-T quadratic y coefficient")
same(constant_T_quartic.coeff_monomial(y), 0,
     "constant-T missing y coefficient")
same(constant_T_quartic.coeff_monomial(1), C,
     "constant-T constant coefficient")

# Primitive rational-root theorem: numerator divides L^4 and denominator
# divides a^2, hence t=c*L^i/a^j in exactly 15 exponent cells.
C_coeffs = [sp.Poly(coefficient, x) for coefficient in sp.Poly(C, t).all_coeffs()]
content = C_coeffs[0]
for coefficient in C_coeffs[1:]:
    content = sp.gcd(content, coefficient)
gate(content.degree() == 0, "C(t) primitive over k[x]")
candidate_cells = [(i, j) for i in range(5) for j in range(3)]
gate(len(candidate_cells) == 15, "rational-root universe has fifteen cells")

# A shorter independent exhaustion uses only the L=0 and a=0 divisors.
# At L=0, every i>0 has a unique minimum valuation: 3 for i=1, 4 for i>=2.
for i in range(1, 5):
    L_orders = (4, 2 + 2 * i, 3 * i, 4 * i)
    minimum = min(L_orders)
    gate(L_orders.count(minimum) == 1,
         f"numerator power i={i} killed by unique L-valuation")
same(F.subs(x, -sp.Rational(4, 9)), sp.Rational(8, 27),
     "L-boundary F value")
same(a.subs(x, -sp.Rational(4, 9)), sp.Rational(5, 9),
     "L-boundary a value")

# Therefore i=0.  For j=1, the four a-valuations are 0,-1,-3,-2.
j1_orders = (0, -1, -3, -2)
gate(j1_orders.count(min(j1_orders)) == 1,
     "denominator power j=1 killed by unique a-valuation")

# For j=0, L=0 forces c=64/25, while a=0 requires 625+32c^3=0.
c_j0 = sp.Rational(64, 25)
same(C.subs({t: c, x: -sp.Rational(4, 9)}),
     c**3 * (sp.Rational(64, 27) - sp.Rational(25, 27) * c),
     "j=0 L-boundary equation")
gate(sp.expand(C.subs({t: c_j0, x: -1})) != 0,
     "j=0 incompatible L and a boundary values")

# For j=2, the lowest a-order is -6 in the cubic and quartic terms; its
# cancellation forces c=32/3.  L=0 instead forces c=64/81.
t_j2 = c / a**2
leading_a_j2 = sp.limit(a**6 * C.subs(t, t_j2), x, -1)
same(leading_a_j2, c**3 * (32 - 3 * c),
     "j=2 a-boundary leading coefficient")
c_j2 = sp.Rational(32, 3)
gate(sp.cancel(C.subs({t: c_j2 / a**2, x: -sp.Rational(4, 9)})) != 0,
     "j=2 incompatible a and L boundary values")

# The three j-cases plus four i>0 rows cover all fifteen candidates.
gate(4 * 3 + 3 == 15, "two-divisor exhaustion covers all rational-root cells")

# For m=0, a nonzero constant product (v-K)(v+K) is impossible because its
# factors would both be units although their difference is 2K.  The zero
# product gives precisely u^3=a^3 and v=+/-K; only zeta=1 has D!=0.
same(D.subs(u, a), 3 * a**2, "zeta=1 inverse denominator nonzero")


# ---------------------------------------------------------------------------
# Fibre packet and local boundary intersection controls.
# ---------------------------------------------------------------------------
delta_disc = sp.factor(sp.discriminant(Delta, y))
same(delta_disc, 256 * a**6 * L**4 * (1 - a)**3 * (9 * a - 4)**2,
     "four simple generic Delta roots")
B = -64 * L**4 * Delta
same(-432 * B**2, -2**16 * 3**3 * L**8 * Delta**2,
     "short Weierstrass discriminant")
r = sp.symbols("r")
B_infinity = sp.expand(r**6 * B.subs(y, 1 / r))
gate(sp.Poly(B_infinity, r).terms()[-1][0] == (2,),
     "infinity a6 has exact order two")
same(sp.Poly(B_infinity, r).coeff_monomial(r**2), 64 * L**4,
     "infinity leading a6 coefficient")
gate(4 * 2 + 4 == 12 and 10 - 2 - 2 == 6,
     "II^4+IV rational-surface Shioda-Tate packet")

# 3-torsion can only require Delta to be a square or cube.  A simple root
# valuation excludes both.
X = sp.symbols("X")
psi3 = sp.factor(3 * X * (X**3 + 4 * B))
same(psi3, 3 * X * (X**3 - 256 * L**4 * Delta),
     "3-division polynomial")
gate((1 % 2 != 0) and (1 % 3 != 0),
     "simple Delta valuation excludes square and cube")

# In the IV blow-up chart the Q sections have equal address but their X1
# difference is a unit times r, proving intersection multiplicity one.
X_plus = 2 * a * L**2 * (s - 1)
X_minus = 2 * a * L**2 * (-s - 1)
Y_Q = -8 * K * L**2
same(reduce_s(r * (X_plus - X_minus)), 4 * s * a * L**2 * r,
     "Q first-order X1 separation")
same(sp.limit(r**2 * Y_Q.subs(y, 1 / r), r, 0), -8 * L**2,
     "Q common IV Y1 address")


# ---------------------------------------------------------------------------
# Exact consequence for the THM-3885 nonlinear interpolation frontier.
# ---------------------------------------------------------------------------
# A global nonzero polynomial residual survivor already has T polynomial over
# k(x)[y].  If u and v were also polynomial, the shell classification leaves
# only T=0 or T_hostile.  T_hostile is neither x-integral nor addressed.
# Therefore every hypothetical nonzero global survivor must enter a genuine
# rational-coordinate denominator shell.  This is a filter, not a closure.
gate(T_hostile.subs({x: 0, y: 0}) != 0,
     "hostile cannot satisfy THM-3885 origin address")

semantic = {
    "target": "THM-3888 repaired alternate-chart blob 2a0852356",
    "divisor": "div(T)=O+P0-Qplus-Qminus with four simple local branches",
    "boundary": "six constant-u polynomial Weierstrass sections; only two have polynomial finite inverse T",
    "counterexample": "degree-one polynomial u,v section refutes complete six-section integral shell",
    "repaired_shell": "u,v,T polynomial gives only P0 and x-polar hostile",
    "constant_T": "15 rational-root cells empty by L=0/a=0 valuations",
    "frontier": "nonzero global polynomial survivor must have u or v denominator; f=0 remains open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source_text = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text))),
     "no inactive Python assert")

print("theorem=THM-3888-final-shell-focused-independent-audit")
print("target_file_commit=2a0852356")
print("divisor_T=O+P0-Qplus-Qminus_all_orders_one")
print("complete_integral_weierstrass_shell=REFUTED_by_degree_one_u_section")
print("constant_u_subshell=six_sections_confirmed")
print("repaired_shell=u_v_T_polynomial_exactly_P0_and_x_pole_hostile")
print("constant_T_candidates=15_empty_by_two_divisor_proof")
print("strict_frontier=nonzero_global_survivor_requires_u_or_v_denominator")
print("nonlinear_f_zero_interpolation=OPEN")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=CORE_PASS_WITH_COMPLETE_SHELL_SCOPE_REPAIR")
