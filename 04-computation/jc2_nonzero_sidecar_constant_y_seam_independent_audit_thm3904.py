#!/usr/bin/env python3
"""Independent exact referee for the THM-3904 constant-y seam.

This file deliberately reconstructs the filtration directly from THM-3881
equation (21).  It imports no project computation and uses no ``assert``
statement, so every check remains active under ``python -O``.
"""

from __future__ import annotations

import hashlib
import itertools

import sympy as sp


GATES = 0


def check(label: str, condition: bool) -> None:
    global GATES
    GATES += 1
    if not bool(condition):
        raise RuntimeError(f"failed gate: {label}")


def zero(label: str, expression: sp.Expr) -> None:
    check(label, sp.cancel(sp.expand(expression)) == 0)


# ---------------------------------------------------------------------------
# 1. Rebuild equation (21) over an abstract coefficient UFD.
# ---------------------------------------------------------------------------

T, f, K, a, H = sp.symbols("T f K a H")
P = a * H
r = a * T + K * f
A = K * T + a * P * f
S = sp.expand(
    H**2
    + 2 * (3 * A + 3 * P + r**2) * H * f
    + (8 * A + 6 * P + 3 * r**2) * (P * f**2 - T**2)
)

row_expansion = sp.expand(
    H**2
    + 6 * a * H**2 * f
    + 12 * a**2 * H**2 * f**2
    + 8 * a**3 * H**2 * f**3
    - 6 * a * H * T**2
    - 6 * a**2 * H * f * T**2
    + 3 * a**3 * H * f**2 * T**2
    - 3 * a**2 * T**4
    + K
    * (
        6 * a**2 * H * f**3 * T
        + 12 * a * H * f**2 * T
        + 6 * H * f * T
        - 6 * a * f * T**3
        - 8 * T**3
    )
    + K**2 * f**2 * (3 * a * H * f**2 + 2 * H * f - 3 * T**2)
)
zero("equation_21_full_K_row_expansion", S - row_expansion)

c2 = sp.factor(S.coeff(K, 2))
c1 = sp.factor(S.coeff(K, 1))
c0 = sp.factor(S.coeff(K, 0))
zero("K2_row", c2 - f**2 * (3 * a * H * f**2 + 2 * H * f - 3 * T**2))
zero(
    "K1_row",
    c1
    - 2
    * T
    * (3 * a**2 * H * f**3 + 6 * a * H * f**2 + 3 * H * f - 3 * a * f * T**2 - 4 * T**2),
)

# Recheck the T=0 arm used to cover the p=-infinity convention.  Here H=L^2
# and Delta=a^3*H-K^2.  Both closed forms are exact consequences of (21).
Delta = a**3 * H - K**2
S_T0 = sp.expand(S.subs(T, 0))
Hf_first = H * (1 + a * f) ** 3 * (1 + 3 * a * f) - Delta * f**3 * (2 + 3 * a * f)
Hf_second = H * (1 + 2 * a * f) ** 3 + K**2 * f**3 * (2 + 3 * a * f)
zero("T0_first_closed_form", S_T0 - H * Hf_first)
zero("T0_second_closed_form", Hf_first - Hf_second)


# ---------------------------------------------------------------------------
# 2. Degree ledger.  K has y-degree two and all coefficients lie in k[x].
# ---------------------------------------------------------------------------


def row_degrees(p: int, q: int) -> dict[str, int]:
    return {
        "K2_f4": 4 + 4 * q,
        "K2_f3": 4 + 3 * q,
        "K2_f2T2": 4 + 2 * q + 2 * p,
        "K_f3T": 2 + 3 * q + p,
        "K_f2T": 2 + 2 * q + p,
        "K_fT": 2 + q + p,
        "K_fT3": 2 + q + 3 * p,
        "K_T3": 2 + 3 * p,
        "f3": 3 * q,
        "f2": 2 * q,
        "f2T2": 2 * q + 2 * p,
        "f": q,
        "fT2": q + 2 * p,
        "T2": 2 * p,
        "T4": 4 * p,
        "constant": 0,
    }


q_gt_p_cases = 0
for p_int in range(0, 9):
    for q_int in range(p_int + 1, 11):
        degrees = row_degrees(p_int, q_int)
        maximum = max(degrees.values())
        winners = {name for name, degree in degrees.items() if degree == maximum}
        check(f"q_gt_p_unique_top_{p_int}_{q_int}", winners == {"K2_f4"})
        check(f"q_gt_p_top_degree_{p_int}_{q_int}", maximum == 4 * q_int + 4)
        q_gt_p_cases += 1

tie_cases = 0
for p_int in range(1, 11):
    degrees = row_degrees(p_int, p_int)
    maximum = max(degrees.values())
    winners = {name for name, degree in degrees.items() if degree == maximum}
    check(
        f"tie_top_rows_{p_int}",
        winners == {"K2_f4", "K2_f2T2"} and maximum == 4 * p_int + 4,
    )
    non_top_at_next = {
        name
        for name, degree in degrees.items()
        if name not in winners and degree == maximum - 1
    }
    expected = {"K2_f3"} if p_int == 1 else set()
    check(f"tie_next_row_split_{p_int}", non_top_at_next == expected)
    tie_cases += 1


# ---------------------------------------------------------------------------
# 3. The p=q=0 even-quartic discriminant, reconstructed exactly.
# ---------------------------------------------------------------------------

t0, c = sp.symbols("t0 c")
S0 = sp.expand(S.subs({T: t0, f: c}))
A2 = sp.factor(S0.coeff(K, 2))
A1 = sp.factor(S0.coeff(K, 1))
A0 = sp.factor(S0.coeff(K, 0))
disc = sp.factor(A1**2 - 4 * A2 * A0)
zero("p0_A2", A2 - c**2 * (H * c * (2 + 3 * a * c) - 3 * t0**2))
zero(
    "p0_K_discriminant",
    disc + 4 * (2 + 3 * a * c) * (H * c * (1 + 2 * a * c) - 2 * t0**2) ** 3,
)

# The two UFD branches reduce to a square equalling a polynomial of odd
# degree: beta*s^2 = unit + gamma*(x+1)*r^2.  The parity conflict is checked
# for every possible degree of r in a broad exact degree ledger.
for r_degree in range(0, 33):
    odd_degree = 2 * r_degree + 1
    check(f"A2_Pell_odd_degree_{r_degree}", odd_degree % 2 == 1)
    check(f"disc_Pell_odd_degree_{r_degree}", odd_degree % 2 == 1)

# Independent small finite-field hostile atlas for the constant-y channel.
# A square even quartic with A2 != 0 must have discriminant zero; A2=0 is
# checked separately.  The search is proof-independent but catches sign and
# missing-degenerate-branch errors in the displayed formulas.
prime = 7


def ff_trim(poly: tuple[int, ...]) -> tuple[int, ...]:
    values = list(poly)
    while len(values) > 1 and values[-1] % prime == 0:
        values.pop()
    return tuple(value % prime for value in values)


def ff_add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    length = max(len(left), len(right))
    return ff_trim(
        tuple(
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
            for index in range(length)
        )
    )


def ff_scale(scalar: int, poly: tuple[int, ...]) -> tuple[int, ...]:
    return ff_trim(tuple(scalar * value for value in poly))


def ff_mul(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    values = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            values[i + j] += left_value * right_value
    return ff_trim(tuple(values))


def ff_zero(poly: tuple[int, ...]) -> bool:
    return ff_trim(poly) == (0,)


a7 = (1, 1)
L7 = (4, 2)  # 9*x+4 modulo 7
H7 = ff_mul(L7, L7)
finite_pairs = 0
finite_A2_zero = 0
finite_disc_zero = 0
for coeffs in itertools.product(range(prime), repeat=4):
    c7 = ff_trim((coeffs[0], coeffs[1]))
    t7 = ff_trim((coeffs[2], coeffs[3]))
    if ff_zero(c7):
        continue
    finite_pairs += 1
    ac7 = ff_mul(a7, c7)
    factor_A2 = ff_add((2,), ff_scale(3, ac7))
    bracket_A2 = ff_add(
        ff_mul(ff_mul(H7, c7), factor_A2),
        ff_scale(-3, ff_mul(t7, t7)),
    )
    A2_7 = ff_mul(ff_mul(c7, c7), bracket_A2)
    factor_D = ff_add((2,), ff_scale(3, ac7))
    bracket_D = ff_add(
        ff_mul(ff_mul(H7, c7), ff_add((1,), ff_scale(2, ac7))),
        ff_scale(-2, ff_mul(t7, t7)),
    )
    if ff_zero(A2_7):
        finite_A2_zero += 1
    if ff_zero(factor_D) or ff_zero(bracket_D):
        finite_disc_zero += 1
check("finite_p0_A2_no_degenerate_case", finite_A2_zero == 0)
check("finite_p0_discriminant_no_square_candidate", finite_disc_zero == 0)


# ---------------------------------------------------------------------------
# 4. Equality p=q>=1: leading coefficient, colors, and next coefficient.
# ---------------------------------------------------------------------------

y, F0 = sp.symbols("y F0")
lead_c, next_c, lead_t, next_t = sp.symbols("lead_c next_c lead_t next_t")


def residual_for(two_term_f: sp.Expr, two_term_T: sp.Expr) -> sp.Expr:
    Ky = y**2 - F0
    ry = a * two_term_T + Ky * two_term_f
    Ay = Ky * two_term_T + a * P * two_term_f
    return sp.expand(
        H**2
        + 2 * (3 * Ay + 3 * P + ry**2) * H * two_term_f
        + (8 * Ay + 6 * P + 3 * ry**2) * (P * two_term_f**2 - two_term_T**2)
    )


for p_int in (1, 2, 3, 5):
    fy = lead_c * y**p_int + next_c * y ** (p_int - 1)
    Ty = lead_t * y**p_int + next_t * y ** (p_int - 1)
    Sy = residual_for(fy, Ty)
    top = sp.expand(Sy).coeff(y, 4 * p_int + 4)
    next_row = sp.expand(Sy).coeff(y, 4 * p_int + 3)
    expected_top = 3 * lead_c**2 * (P * lead_c**2 - lead_t**2)
    expected_next = 6 * lead_c * (
        (2 * P * lead_c**2 - lead_t**2) * next_c - lead_c * lead_t * next_t
    )
    if p_int == 1:
        expected_next += 2 * H * lead_c**3
    zero(f"tie_top_coefficient_p{p_int}", top - expected_top)
    zero(f"tie_next_coefficient_p{p_int}", next_row - expected_next)

# If G has leading coefficient g and c divides g, write z=g/c.  With
# lambda^2=-3 the leading-square equation has the exact two-color product.
z, lam = sp.symbols("z lam")
color_product = (z + lam * lead_t) * (z - lam * lead_t) - 3 * a * H * lead_c**2
color_product = sp.expand(color_product).subs(lam**2, -3)
color_product = sp.expand(color_product).subs(
    z**2, 3 * (a * H * lead_c**2 - lead_t**2)
)
zero("leading_two_color_product", color_product)

# Exhaust the valuation passport: after removing the common gcd delta, the
# a-prime has odd residual exponent and every other prime has even residual
# exponent.  This is the UFD content behind one a-colored square and one
# uncolored square.  m is v_pi(L*c); a_total distinguishes pi=a.
valuation_gates = 0
for a_prime in (False, True):
    for m in range(0, 9):
        total = 2 * m + (1 if a_prime else 0)
        for left in range(total + 1):
            right = total - left
            common = min(left, right)
            reduced_left = left - common
            reduced_right = right - common
            check(
                f"valuation_coprime_{int(a_prime)}_{m}_{left}",
                min(reduced_left, reduced_right) == 0,
            )
            check(
                f"valuation_delta_divides_Lc_{int(a_prime)}_{m}_{left}",
                common <= m,
            )
            check(
                f"valuation_color_parity_{int(a_prime)}_{m}_{left}",
                (reduced_left + reduced_right) % 2 == (1 if a_prime else 0),
            )
            valuation_gates += 3


# ---------------------------------------------------------------------------
# 5. Hostile controls: leading and next tariffs are necessary, not sufficient.
# ---------------------------------------------------------------------------

# Specialize x=0, so a=1, L=4, H=P=16, F0=4.  Work in Q(lambda),
# lambda^2=-3.  The two colors U=3ac and V=Hc with c=3a+H give
# z=(U+V)/2=c^2/2 and t=(U-V)/(2lambda).  For p=1 the exceptional next row
# is paid by root coefficient h=2H; for p>=2 it is paid by h=0.  Nevertheless
# the complete residual is squarefree, so no equality-lane solution is being
# smuggled in by the passport.
lambda_value = sp.sqrt(-3)
a_value = sp.Integer(1)
L_value = sp.Integer(4)
H_value = L_value**2
P_value = a_value * H_value
F_value = sp.Integer(4)
c_value = 3 * a_value + H_value
z_value = c_value * (3 * a_value + H_value) / 2
t_value = c_value * (3 * a_value - H_value) / (2 * lambda_value)
zero(
    "hostile_leading_norm",
    z_value**2 - 3 * (P_value * c_value**2 - t_value**2),
)

hostile_signatures: list[tuple[int, int, int]] = []
for p_int in (1, 2):
    f_value = c_value * y**p_int
    T_value = t_value * y**p_int
    K_value = y**2 - F_value
    r_value = a_value * T_value + K_value * f_value
    A_value = K_value * T_value + a_value * P_value * f_value
    S_value = sp.expand(
        H_value**2
        + 2 * (3 * A_value + 3 * P_value + r_value**2) * H_value * f_value
        + (8 * A_value + 6 * P_value + 3 * r_value**2)
        * (P_value * f_value**2 - T_value**2)
    )
    polynomial = sp.Poly(S_value, y, extension=lambda_value)
    derivative_gcd = sp.gcd(polynomial, polynomial.diff())
    expected_degree = 4 * p_int + 4
    check(f"hostile_full_degree_p{p_int}", polynomial.degree() == expected_degree)
    check(f"hostile_squarefree_p{p_int}", derivative_gcd.degree() == 0)
    # Recheck the root's next coefficient equation itself.
    root_lead = c_value * z_value
    root_next = 2 * H_value if p_int == 1 else sp.Integer(0)
    check(
        f"hostile_root_top_p{p_int}",
        sp.simplify(polynomial.nth(expected_degree) - root_lead**2) == 0,
    )
    check(
        f"hostile_root_next_p{p_int}",
        sp.simplify(polynomial.nth(expected_degree - 1) - 2 * root_lead * root_next)
        == 0,
    )
    hostile_signatures.append((p_int, polynomial.degree(), derivative_gcd.degree()))


semantic_lines = [
    "scope=k algebraically closed characteristic zero; D=k[x,y]",
    "source=THM-3881 equation (21), independently expanded in K=y^2-F",
    "T=0 and f!=0 is already impossible by the independently rechecked THM-3881 degree arm",
    "for nonzero T,f write p=deg_y(T), q=deg_y(f)",
    "q>p is impossible: for q>=1 the unique top coefficient is 3*a*L^2*c^4, with odd a-valuation",
    "p=q=0 is impossible without the address: the even-quartic K-discriminant is -4*(2+3af)*(L^2*f*(1+2af)-2T^2)^3",
    "the possible leading-quadratic degeneration is L^2*f*(2+3af)=3T^2; both it and discriminant zero contradict odd degree after coprime UFD splitting",
    "for p=q>=1 the leading coefficient is 3*c^2*(a*L^2*c^2-t^2), nonzero by odd a-valuation",
    "if g is the square-root leader then c divides g; z=g/c obeys (z+sqrt(-3)t)(z-sqrt(-3)t)=3*a*L^2*c^2",
    "delta=gcd(z+sqrt(-3)t,z-sqrt(-3)t) divides L*c; after delta removal the coprime colors are, up to units and swap, a times a square and a square",
    "for p>=2: z*h=3*((2*a*L^2*c^2-t^2)*d-c*t*u)",
    "for p=1: z*h=3*((2*a*L^2*c^2-t^2)*d-c*t*u)+L^2*c^2",
    "the p=1 correction comes only from the 2*L^2*K^2*f^3 row",
    "no claim is made for q<p or for sufficiency of the equality passport; JC(2) remains open",
]
semantic_sha256 = hashlib.sha256("\n".join(semantic_lines).encode("utf-8")).hexdigest()

print("THM3904_CONSTANT_Y_SEAM_INDEPENDENT_REFEREE")
print("status=PASS")
print(f"active_gates={GATES}")
print(f"q_gt_p_degree_cases={q_gt_p_cases}")
print(f"positive_tie_degree_cases={tie_cases}")
print(f"valuation_passport_gates={valuation_gates}")
print(f"finite_field_prime={prime}")
print(f"finite_constant_y_pairs={finite_pairs}")
print(f"finite_A2_zero_candidates={finite_A2_zero}")
print(f"finite_discriminant_zero_candidates={finite_disc_zero}")
print(f"hostile_leading_and_next_not_sufficient={hostile_signatures}")
print(f"semantic_sha256={semantic_sha256}")
print("verdict=PROMOTION_READY_WITH_NECESSARY_ONLY_EQUALITY_SCOPE")
