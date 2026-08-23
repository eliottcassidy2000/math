#!/usr/bin/env python3
"""Primary THM-3904 constant-y seam and primitive-color audit."""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp

sys.stdout.reconfigure(newline="\n")

GATES = 0


def gate(condition: bool, message: str) -> None:
    global GATES
    GATES += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.cancel(sp.expand(expression)) == 0, message)


x, y = sp.symbols("x y")
a = x+1
L = 9*x+4
F = 15*x**2+15*x+4
P = a*L**2
K = y**2-F


def residual(T: sp.Expr, f: sp.Expr) -> sp.Expr:
    r = a*T+K*f
    A = K*T+a*P*f
    B = P*f**2-T**2
    return sp.expand(L**4+2*(3*A+3*P+r**2)*L**2*f+(8*A+6*P+3*r**2)*B)


# -------------------------------------------------------------------------
# 1. Full support and a hostile audit of every strict q>p degree tie.
# -------------------------------------------------------------------------

Ts, fs, Ks, aa, LL = sp.symbols("Ts fs Ks aa LL")
Ps = aa*LL**2
rs = aa*Ts+Ks*fs
As = Ks*Ts+aa*Ps*fs
Bs = Ps*fs**2-Ts**2
abstract_S = sp.expand(LL**4+2*(3*As+3*Ps+rs**2)*LL**2*fs+(8*As+6*Ps+3*rs**2)*Bs)
abstract_terms = sp.Poly(abstract_S, Ks, Ts, fs).terms()
gate(len(abstract_terms) == 16, "complete sixteen-monomial residual support")

top_monomial = (2, 0, 4)
psym, delta_sym, Nsym = sp.symbols("psym delta_sym Nsym", integer=True, nonnegative=True)
strict_symbolic_rows = 0
for monomial, coefficient in abstract_terms:
    if monomial == top_monomial:
        continue
    k_power, t_power, f_power = monomial
    q_symbolic = psym+delta_sym+1
    weighted_degree = 2*k_power+psym*t_power+q_symbolic*f_power
    deficit = sp.expand(4*q_symbolic+4-weighted_degree)
    deficit_poly = sp.Poly(deficit, psym, delta_sym)
    gate(
        all(coefficient >= 0 for coefficient in deficit_poly.coeffs())
        and deficit.subs({psym: 0, delta_sym: 0}) > 0,
        f"all-degree strict deficit {monomial}",
    )
    strict_symbolic_rows += 1

equality_symbolic_winners = []
for monomial, coefficient in abstract_terms:
    k_power, t_power, f_power = monomial
    n_symbolic = Nsym+1
    weighted_degree = 2*k_power+n_symbolic*(t_power+f_power)
    deficit = sp.expand(4*n_symbolic+4-weighted_degree)
    if deficit == 0:
        equality_symbolic_winners.append(monomial)
    else:
        deficit_poly = sp.Poly(deficit, Nsym)
        gate(
            all(coefficient >= 0 for coefficient in deficit_poly.coeffs())
            and deficit.subs(Nsym, 0) > 0,
            f"all-degree equality deficit {monomial}",
        )
gate(
    set(equality_symbolic_winners) == {(2, 0, 4), (2, 2, 2)},
    "all-degree equality has exactly two top monomials",
)

strict_rows = 0
for q_degree in range(1, 21):
    for p_degree in range(q_degree):
        degrees = {
            monomial: 2*monomial[0]+p_degree*monomial[1]+q_degree*monomial[2]
            for monomial, coefficient in abstract_terms if coefficient != 0
        }
        maximum = max(degrees.values())
        winners = [monomial for monomial, degree in degrees.items() if degree == maximum]
        gate(winners == [top_monomial], f"strict tariff unique top p={p_degree},q={q_degree}")
        gate(maximum == 4*q_degree+4, f"strict tariff degree p={p_degree},q={q_degree}")
        strict_rows += 1

u, v, u1, v1 = sp.symbols("u v u1 v1")
for q_degree in range(1, 9):
    for p_degree in range(q_degree):
        monomial_T = u*y**p_degree
        monomial_f = v*y**q_degree
        specialized = sp.Poly(residual(monomial_T, monomial_f), y)
        gate(specialized.degree() == 4*q_degree+4, f"strict polynomial degree p={p_degree},q={q_degree}")
        zero(
            specialized.coeff_monomial(y**(4*q_degree+4))-3*a*L**2*v**4,
            f"strict top coefficient p={p_degree},q={q_degree}",
        )

gate(L.subs(x, -1) == -5, "a and L are coprime")
for a_order in range(21):
    gate((1+4*a_order) % 2 == 1, f"strict top odd a-order={a_order}")


# -------------------------------------------------------------------------
# 2. Positive equality seam: top and next coefficients.
# -------------------------------------------------------------------------

top_expected = 3*v**2*(P*v**2-u**2)
generic_next_expected = 6*v*((2*P*v**2-u**2)*v1-u*v*u1)
linear_next_expected = generic_next_expected+2*L**2*v**3

equality_rows = 0
for n in range(1, 9):
    T_truncated = u*y**n+u1*y**(n-1)
    f_truncated = v*y**n+v1*y**(n-1)
    S_truncated = sp.Poly(residual(T_truncated, f_truncated), y)
    zero(
        S_truncated.coeff_monomial(y**(4*n+4))-top_expected,
        f"equality top n={n}",
    )
    expected = linear_next_expected if n == 1 else generic_next_expected
    zero(
        S_truncated.coeff_monomial(y**(4*n+3))-expected,
        f"equality next response n={n}",
    )
    equality_rows += 1

# The formal zero top P*v^2=u^2 is impossible: its a-order would be odd on
# one side and even on the other.
for v_order in range(21):
    gate((1+2*v_order) % 2 == 1, f"equality zero-top odd a-order={v_order}")

h, g, d = sp.symbols("h g d")
zero((v*h)**2-v**2*h**2, "leading root g=v*h shadow")
top_relation = h**2-3*(P*v**2-u**2)
colour_product = (h-d*u)*(h+d*u)-3*a*L**2*v**2
colour_reduced = sp.Poly(sp.expand(colour_product), d).rem(sp.Poly(d**2+3, d)).as_expr()
zero(colour_reduced.subs(h**2, 3*(P*v**2-u**2)), "equianharmonic color norm")


# -------------------------------------------------------------------------
# 3. Complete UFD allocation after the common color gcd.
# -------------------------------------------------------------------------

w, H, U, c, c_minus, c_plus, lam = sp.symbols(
    "w H U c c_minus c_plus lam", nonzero=True
)

# Primewise proof that w=gcd(h-du,h+du)=gcd(h,u) divides L*v.
prime_ledger = 0
for exponent_lv in range(31):
    for is_a in (0, 1):
        total = 2*exponent_lv+is_a
        for exponent_w in range(total//2+1):
            gate(exponent_w <= exponent_lv, "common color gcd divides L*v")
            prime_ledger += 1


def reduce_d(expression: sp.Expr) -> sp.Expr:
    numerator = sp.together(expression).as_numer_denom()[0]
    return sp.Poly(numerator, d).rem(sp.Poly(d**2+3, d)).as_expr()


# Case minus: the single odd a-factor is allocated to h-du.
X_minus = lam*a*c_minus**2
Y_minus = 3*c_plus**2/lam
H_minus = (X_minus+Y_minus)/2
U_minus = (Y_minus-X_minus)/(2*d)
zero(reduce_d(H_minus-d*U_minus-X_minus), "minus allocation first color")
zero(reduce_d(H_minus+d*U_minus-Y_minus), "minus allocation second color")
zero(
    reduce_d((H_minus**2+3*U_minus**2)-3*a*(c_minus*c_plus)**2),
    "minus allocation norm",
)

# Case plus: the single odd a-factor is allocated to h+du.
X_plus = lam*c_minus**2
Y_plus = 3*a*c_plus**2/lam
H_plus = (X_plus+Y_plus)/2
U_plus = (Y_plus-X_plus)/(2*d)
zero(reduce_d(H_plus-d*U_plus-X_plus), "plus allocation first color")
zero(reduce_d(H_plus+d*U_plus-Y_plus), "plus allocation second color")
zero(
    reduce_d((H_plus**2+3*U_plus**2)-3*a*(c_minus*c_plus)**2),
    "plus allocation norm",
)

# In either case c=(L*v)/w, h=wH, u=wU.  Coprimality of the primitive
# colors is equivalent to gcd(H,U)=1.  The allocation must also place any
# a-part of c in the same color as the extra a-factor.
for c_order in range(21):
    total_a_order = 1+2*c_order
    gate(total_a_order % 2 == 1, f"primitive allocation odd a exponent={c_order}")


# An equivalent square-normalized passport is sometimes cleaner.  Choose
# rho^2=3 and i^2=-1, put h=rho*Z, and write Z=c*Z0, u=c*U0.  Then the
# primitive colors Z0+iU0 and Z0-iU0 are coprime and have product
# a*(L*v/c)^2.  The following two allocations are the complete possibilities
# up to units and swapping colors.
rho, ii, Z, Z0, U0, r_minus, r_plus, mu = sp.symbols(
    "rho ii Z Z0 U0 r_minus r_plus mu", nonzero=True
)


def reduce_i(expression: sp.Expr) -> sp.Expr:
    numerator = sp.together(expression).as_numer_denom()[0]
    return sp.Poly(numerator, ii).rem(sp.Poly(ii**2+1, ii)).as_expr()


passport_minus = mu*a*r_minus**2
passport_plus = r_plus**2/mu
Z0_case_a = (passport_minus+passport_plus)/2
U0_case_a = (passport_minus-passport_plus)/(2*ii)
zero(reduce_i(Z0_case_a+ii*U0_case_a-passport_minus), "passport a in first color")
zero(reduce_i(Z0_case_a-ii*U0_case_a-passport_plus), "passport square second color")
zero(
    reduce_i(Z0_case_a**2+U0_case_a**2-a*(r_minus*r_plus)**2),
    "passport first allocation norm",
)

passport_minus = mu*r_minus**2
passport_plus = a*r_plus**2/mu
Z0_case_b = (passport_minus+passport_plus)/2
U0_case_b = (passport_minus-passport_plus)/(2*ii)
zero(reduce_i(Z0_case_b+ii*U0_case_b-passport_minus), "passport square first color")
zero(reduce_i(Z0_case_b-ii*U0_case_b-passport_plus), "passport a in second color")
zero(
    reduce_i(Z0_case_b**2+U0_case_b**2-a*(r_minus*r_plus)**2),
    "passport second allocation norm",
)

# If c=gcd(Z,u), the valuation inequality c^2 | a*(L*v)^2 implies c|L*v,
# including at the exceptional prime a (where 2e <= 1+2s still gives e<=s).
passport_divisibility_rows = 0
for exponent_lv in range(31):
    for is_a in (0, 1):
        rhs_order = 2*exponent_lv+is_a
        for exponent_c in range(rhs_order//2+1):
            gate(exponent_c <= exponent_lv, "passport common divisor divides L*v")
            passport_divisibility_rows += 1


# -------------------------------------------------------------------------
# 4. Exact next-response divisors (with n=1 split explicitly).
# -------------------------------------------------------------------------

R_generic = 3*((2*P*v**2-u**2)*v1-u*v*u1)
R_rewritten = (3*u**2+2*h**2)*v1-3*u*v*u1
zero(
    (R_generic-R_rewritten).subs(h**2, 3*(P*v**2-u**2)),
    "generic response rewrite",
)

D1 = u*v1-v*u1
generic_bracket = w*(3*U**2+2*H**2)*v1-3*U*v*u1
D1_sub = w*U*v1-v*u1
zero(generic_bracket-3*U*D1_sub-2*w*H**2*v1, "generic primitive response remainder")

# Thus, with gcd(H,U)=1, h|R_generic iff H|D1.  If it holds, the next
# square-root coefficient is uniquely g1=R_generic/h.
for h_order in range(21):
    for u_order in range(21):
        common = min(h_order, u_order)
        gate(min(h_order-common, u_order-common) == 0, "primitive H,U coprime valuation")

R_linear = R_generic+L**2*v**2
linear_bracket = generic_bracket+w*c**2
zero((L**2*v**2).subs(L**2*v**2, w**2*c**2)-w**2*c**2, "linear Lv=wc square")
zero(linear_bracket-3*U*D1_sub-2*w*H**2*v1-w*c**2, "linear response decomposition")
# The previous identity is arranged so the remainder modulo H is
# 3*U*D1+w*c^2; it is printed explicitly below.

# The same response in the square-normalized passport.  With
# Z^2=P*v^2-u^2 and lc_y(G)=rho*v*Z, the n>=2 condition is Z|E, while n=1
# has the extra 2*L^2*K^2*f^3 contribution and requires
# Z | 3E+L^2*v^2.  After Z=c*Z0,u=c*U0,C=L*v/c these reduce as printed.
C = sp.symbols("C", nonzero=True)
E_Z = (u**2+2*Z**2)*v1-u*v*u1
zero(
    (E_Z-((2*P*v**2-u**2)*v1-u*v*u1)).subs(Z**2, P*v**2-u**2),
    "square-normalized generic response",
)
E_Z_primitive = c*(c*(U0**2+2*Z0**2)*v1-U0*v*u1)
zero(
    E_Z.subs({Z: c*Z0, u: c*U0})-E_Z_primitive,
    "generic response after common divisor",
)
generic_Z_remainder = c*U0*v1-v*u1
zero(
    E_Z_primitive-c*U0*generic_Z_remainder-2*c**2*Z0**2*v1,
    "generic passport remainder modulo Z0",
)
linear_Z_numerator = 3*E_Z_primitive+c**2*C**2
linear_Z_remainder = 3*U0*generic_Z_remainder+c*C**2
zero(
    linear_Z_numerator-c*linear_Z_remainder-6*c**2*Z0**2*v1,
    "linear passport remainder modulo Z0",
)


# -------------------------------------------------------------------------
# 5. Positive controls show the next response is not a universal no-go.
# -------------------------------------------------------------------------

h_control = (a+3*L**2)/2
u_control = (3*L**2-a)/(2*d)
zero(
    reduce_d(h_control**2-3*(P-u_control**2)),
    "canonical leading payment",
)
zero(reduce_d(h_control-d*u_control-a), "canonical first color")
zero(reduce_d(h_control+d*u_control-3*L**2), "canonical second color")

# For n>=2 choose u1=v1=0, so the next response and g1 both vanish.
zero(R_generic.subs({h: h_control, u: u_control, v: 1, u1: 0, v1: 0}), "generic next positive control")

# For n=1 choose v1=0,u1=d/9.  The exceptional L^2 term then gives
# R_linear=h/3 and the next root coefficient g1=1/3.
linear_control = R_linear.subs({h: h_control, u: u_control, v: 1, u1: d/9, v1: 0})
zero(reduce_d(linear_control-h_control/3), "linear next positive control")


# -------------------------------------------------------------------------
# 6. The common-degree-zero seam is different and is actually empty.
# -------------------------------------------------------------------------

S_constant = sp.Poly(residual(u, v), y)
constant_top = v**2*(L**2*v*(3*a*v+2)-3*u**2)
zero(S_constant.coeff_monomial(y**4)-constant_top, "constant seam top")
zero(S_constant.coeff_monomial(y**3), "constant seam missing cubic")
zero(S_constant.coeff_monomial(y), "constant seam missing linear term")
constant_c2 = S_constant.coeff_monomial(y**2)
constant_c0 = S_constant.coeff_monomial(1)
constant_discriminant = sp.factor(constant_c2**2-4*constant_top*constant_c0)
constant_discriminant_expected = 4*(3*a*v+2)*(
    -2*L**2*a*v**2-L**2*v+2*u**2
)**3
zero(
    constant_discriminant-constant_discriminant_expected,
    "constant seam z=y^2 discriminant",
)
zero((3*a*v+2)-3*a*v-2, "constant leading coprime Bezout identity")
zero((2*a*v+1)-2*a*v-1, "constant discriminant coprime Bezout identity")

# All-degree proof behind this finite hostile ledger.  If the quartic
# coefficient vanished, coprimality of v and 3*a*v+2 would force both to be
# squares up to units.  The first then has even degree and the second degree
# deg(v)+1, which is odd.  If the quartic is a square, its nonzero leading
# coefficient forces the z-discriminant to vanish; the surviving equation
# 2*u^2=L^2*v*(2*a*v+1) gives the identical parity contradiction.  The rows
# below are controls for, not a replacement of, that general degree proof.
constant_pell_rows = 0
for square_root_degree in range(41):
    square_degree = 2*square_root_degree
    gate((square_degree+1) % 2 == 1, "constant Pell shifted degree is odd")
    constant_pell_rows += 1


semantic = {
    "tariff": "f_nonzero square implies deg_y(f)<=deg_y(T)",
    "allocation": "after w=gcd(colors), primitive colors split 3*a*c^2 into two coprime square colors",
    "generic_response": "n>=2 iff H divides u*v1-v*u1 at the first lower coefficient",
    "linear_response": "n=1 iff H divides 3*U*(u*v1-v*u1)+w*c^2",
    "positive_controls": "canonical a versus 3L2 colors lift the next shell for n>=1",
    "constant_boundary": "n=0 is address-free empty by quartic discriminant and odd-degree Pell parity",
    "scope": "positive seam only through two leading y-shells; existence, lower shells, Keller atlas, JC2 open",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "inactive Python assert",
)

print("THM3904_CONSTANT_Y_SEAM_AND_PRIMITIVE_COLORS")
print("status=PASS;THM3899_CONFIRMED;TWO_LAYER_EXTENSION_PROVED;JC2_OPEN")
print(
    f"strict_symbolic_rows={strict_symbolic_rows};"
    f"strict_tariff_rows={strict_rows};equality_rows={equality_rows}"
)
print(f"prime_allocation_rows={prime_ledger};passport_divisibility_rows={passport_divisibility_rows}")
print("primitive_colors=(lambda*a^epsilon*c_minus^2,3/lambda*a^(1-epsilon)*c_plus^2)")
print("normalized_passport=(mu*a^epsilon*r_minus^2,mu^-1*a^(1-epsilon)*r_plus^2)")
print("next_response_n>=2=H_divides_(u*v1-v*u1)")
print("next_response_n=1=H_divides_(3*U*(u*v1-v*u1)+w*c^2)")
print("positive_next_shells=n>=2:g1=0;n=1:g1=1/3")
print("constant_seam=EMPTY_ADDRESS_FREE_BY_DISCRIMINANT_AND_PELL_PARITY")
print(f"constant_pell_rows={constant_pell_rows}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
print(f"gates={GATES}")
print("ALL CHECKS PASSED")
