#!/usr/bin/env python3
"""Exact degree-26 Faber and valuation controls for THM-2823."""

import ast
import math
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def has_asserts(path):
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    return any(isinstance(node, ast.Assert) for node in ast.walk(tree))


def zero(expr):
    return sp.cancel(sp.together(expr)) == 0


q, d, s, T = sp.symbols("q d s T")
c18, c14, c10, c6, c2 = sp.symbols("c18 c14 c10 c6 c2")
seeds = (c18, c14, c10, c6, c2)


def recurrence_coefficients(m):
    """Coefficients of (1+2dt^2+qt^3+(d^2-s)t^4)^(m/4)."""
    alpha = sp.Rational(m, 4)
    coefficients = [sp.Integer(1)]
    base = {2: 2 * d, 3: q, 4: d**2 - s}
    for n in range(1, m + 4):
        value = sum(
            ((alpha + 1) * k - n) * base[k] * coefficients[n - k]
            for k in base
            if k <= n
        ) / n
        coefficients.append(sp.expand(value))
    return coefficients


def multinomial_coefficient(m, n):
    """Independent direct coefficient extraction at one exponent n."""
    alpha = sp.Rational(m, 4)
    answer = 0
    for i in range(n // 2 + 1):
        for j in range(n // 3 + 1):
            remainder = n - 2 * i - 3 * j
            if remainder < 0 or remainder % 4:
                continue
            k = remainder // 4
            total = i + j + k
            multinomial = math.factorial(total) // (
                math.factorial(i) * math.factorial(j) * math.factorial(k)
            )
            answer += (
                sp.binomial(alpha, total)
                * multinomial
                * (2 * d) ** i
                * q**j
                * (d**2 - s) ** k
            )
    return sp.expand(answer)


def replace_even_q(expr):
    polynomial = sp.Poly(sp.expand(expr), q)
    answer = 0
    for (exponent,), coefficient in polynomial.terms():
        require(exponent % 2 == 0, "unexpected odd q power")
        answer += coefficient * T ** (exponent // 2)
    return sp.expand(answer)


rows = {}
for m in (26, 18, 14, 10, 6, 2):
    coefficients = recurrence_coefficients(m)
    for index in (m + 1, m + 2, m + 3):
        require(
            zero(coefficients[index] - multinomial_coefficient(m, index)),
            f"recurrence/direct mismatch m={m}, n={index}",
        )
    rows[m] = {
        "phi": replace_even_q(sp.cancel(4 * coefficients[m + 1] / q)),
        "psi": replace_even_q(4 * coefficients[m + 2]),
        "K": replace_even_q(
            sp.cancel((4 * coefficients[m + 3] + 2 * d * coefficients[m + 1]) / q)
        ),
    }

bank_coefficients = {
    26: sp.Integer(1),
    18: c18,
    14: c14,
    10: c10,
    6: c6,
    2: c2,
}
Phi = sp.expand(
    sum(bank_coefficients[m] * rows[m]["phi"] for m in bank_coefficients)
)
Psi = sp.expand(
    sum(bank_coefficients[m] * rows[m]["psi"] for m in bank_coefficients)
)
K = sp.expand(
    sum(bank_coefficients[m] * rows[m]["K"] for m in bank_coefficients)
)

H = sp.expand(
    (
        72 * c18 * T**2
        - 4032 * c18 * T * d * s
        + 6720 * c18 * s**3
        + 896 * c14 * T * d
        - 4480 * c14 * s**2
        + 2560 * c10 * s
        - 1024 * c6
        - 143 * T**3 * d
        + 5148 * T**2 * d**2 * s
        + 1287 * T**2 * s**2
        - 24024 * T * d * s**3
        + 12012 * s**5
    )
    / 4096
)
require(zero(K + d * Phi / 2 - T * H), "K decomposition")

# A legal constant translation P_c=P+c kills the E22 coefficient.
translation_c, a22 = sp.symbols("translation_c a22")
translated_e22_coefficient = a22 - sp.Rational(13, 2) * translation_c
require(
    zero(translated_e22_coefficient.subs(translation_c, 2 * a22 / 13)),
    "degree-26 target translation",
)


def initial_form(polynomial, weights):
    """Return the minimum weighted part; seed coefficients have weight zero."""
    variables = (T, d, s) + seeds
    poly = sp.Poly(sp.expand(polynomial), *variables)
    weighted_terms = []
    for exponents, coefficient in poly.terms():
        weight = sum(exponents[index] * weights[index] for index in range(3))
        monomial = coefficient
        for variable, exponent in zip(variables, exponents):
            monomial *= variable**exponent
        weighted_terms.append((weight, monomial))
    minimum = min(weight for weight, _ in weighted_terms)
    answer = sp.expand(
        sum(monomial for weight, monomial in weighted_terms if weight == minimum)
    )
    return minimum, answer


r, zeta = sp.symbols("r zeta")
f = r**3 - 21 * r**2 + 35 * r - 7
g = 7 * r**3 - 35 * r**2 + 21 * r - 1
generic_cases = ((0, 0), (1, 0), (2, 0), (1, 1))
for a, b in generic_cases:
    weights = (2 * a - 3, 2 * b - 3, a + b - 3)
    phi_weight, phi_initial = initial_form(Phi, weights)
    psi_weight, psi_initial = initial_form(Psi, weights)
    r_expr = T * d / s**2
    expected_phi = -sp.Rational(13728, 16384) * s**6 * f.subs(r, r_expr)
    expected_psi = sp.Rational(13728, 16384) * s**7 * g.subs(r, r_expr)
    require(zero(phi_initial - expected_phi), f"generic Phi initial ({a},{b})")
    require(zero(psi_initial - expected_psi), f"generic Psi initial ({a},{b})")
    require(phi_weight == 6 * weights[2], f"generic Phi weight ({a},{b})")
    require(psi_weight == 7 * weights[2], f"generic Psi weight ({a},{b})")

require(sp.resultant(f, g, r) == -(2**21), "generic cubic resultant")

# Exceptional (a,b)=(0,1): two independent leading ratios survive.
exceptional_weights = (-3, -1, -2)
phi_weight, phi_initial = initial_form(Phi, exceptional_weights)
psi_weight, psi_initial = initial_form(Psi, exceptional_weights)
h_weight, h_initial = initial_form(H, exceptional_weights)

P_exceptional = (
    143 * zeta**2
    - 13728 * r**3
    - 20592 * zeta * r
    + 288288 * r**2
    + 48048 * zeta
    - 480480 * r
    + 96096
)
G_exceptional = (
    1287 * zeta**2
    + 5148 * zeta * r**2
    - 96096 * r**3
    - 72072 * zeta * r
    + 480480 * r**2
    + 60060 * zeta
    - 288288 * r
    + 13728
)
H_exceptional = 143 * (
    (9 - r) * zeta + 36 * r**2 - 168 * r + 84
)
r_expr = T * d / s**2
zeta_expr = T**2 / s**3
require(
    zero(
        phi_initial
        - s**6
        * P_exceptional.subs({r: r_expr, zeta: zeta_expr})
        / 16384
    ),
    "exceptional Phi initial",
)
require(
    zero(
        psi_initial
        + s**7
        * G_exceptional.subs({r: r_expr, zeta: zeta_expr})
        / 16384
    ),
    "exceptional Psi initial",
)
require(
    zero(
        h_initial
        - s**5
        * H_exceptional.subs({r: r_expr, zeta: zeta_expr})
        / 4096
    ),
    "exceptional H initial",
)
require((phi_weight, psi_weight, h_weight) == (-12, -14, -10),
        "exceptional initial weights")

exceptional_groebner = sp.groebner(
    [P_exceptional, G_exceptional, H_exceptional],
    zeta,
    r,
    order="lex",
    domain=sp.QQ,
)
require(
    len(exceptional_groebner.polys) == 1
    and exceptional_groebner.polys[0].as_expr() == 1,
    "exceptional unit ideal",
)

p_quintic = (
    2 * r**5
    + 3 * r**4
    - 488 * r**3
    + 2842 * r**2
    - 6930 * r
    + 4011
)
q_quintic = (
    13 * r**5
    - 182 * r**4
    + 1974 * r**3
    - 8776 * r**2
    + 13173 * r
    - 5130
)
quintic_resultant = sp.resultant(p_quintic, q_quintic, r)
require(
    quintic_resultant == -(2**36) * 3**2 * 31**9 * 37,
    "exceptional quintic resultant",
)
require(
    H_exceptional.subs(r, 9) != 0,
    "r=9 exceptional denominator boundary",
)
zeta_solution = 12 * (3 * r**2 - 14 * r + 7) / (r - 9)
require(
    zero(
        P_exceptional.subs(zeta, zeta_solution)
        + 6864 * p_quintic / (r - 9) ** 2
    ),
    "exceptional P elimination",
)
require(
    zero(
        G_exceptional.subs(zeta, zeta_solution)
        - 6864 * q_quintic / (r - 9) ** 2
    ),
    "exceptional G elimination",
)

# Exact monomial-weight controls for the three regular valuation regions.
minimum, _ = initial_form(H, (3, -3, 0))
require(minimum >= 0, "a>=3 worst H bound")
minimum, _ = initial_form(H, (1, -1, 0))
require(minimum >= 0, "a=2,b>=1 worst H bound")
for label, weights in (
    ("a=0,b=2", (-3, 0, -1)),
    ("a=0,b>=3", (-3, 0, 0)),
    ("a=1,b>=2", (-1, 0, 0)),
):
    minimum, leading = initial_form(Phi, weights)
    require(
        zero(leading - sp.Rational(143, 16384) * T**4),
        f"{label}: unique T^4 pole",
    )

# Both sextic carriers have an explicit finite root with (v(V),v(M))=(3,1).
x, y = sp.symbols("x y")
V_power = sp.Rational(16, 9) * x**5 * (x**3 - 1) ** 3
M_power = x * (x**3 - 1) * (x**3 - sp.Rational(1, 2))
require(
    sp.Poly(V_power, x).exquo(sp.Poly((x**3 - 1) ** 3, x)).as_expr() != 0,
    "power V triple factor",
)
require(
    sp.rem(sp.Poly(M_power, x), sp.Poly(x**3 - 1, x)).is_zero,
    "power M simple factor",
)
require(
    not sp.rem(
        sp.Poly(M_power / (x**3 - 1), x),
        sp.Poly(x**3 - 1, x),
    ).is_zero,
    "power M factor is exactly simple",
)

T3 = 4 * y**3 - 3 * y
V_chebyshev = (
    sp.Rational(256, 9) * (y**2 - sp.Rational(1, 4)) ** 4 * (y**2 - 1) ** 3
)
M_chebyshev = (
    T3
    / 4
    * (y**2 - sp.Rational(1, 4))
    * (y**2 - 1)
)
require(
    sp.rem(sp.Poly(V_chebyshev, y), sp.Poly((y**2 - 1) ** 3, y)).is_zero,
    "Chebyshev V triple factor",
)
require(
    sp.rem(sp.Poly(M_chebyshev, y), sp.Poly(y**2 - 1, y)).is_zero,
    "Chebyshev M simple factor",
)
require(
    not sp.rem(
        sp.Poly(M_chebyshev / (y**2 - 1), y),
        sp.Poly(y**2 - 1, y),
    ).is_zero,
    "Chebyshev M factor is exactly simple",
)

require(not has_asserts(__file__), "Python assert nodes are forbidden")

print("THM-2823 DEGREE-26 TRIPLE-POLE AUDIT -- exact")
print("assert_nodes=0")
print("Faber rows 26,18,14,10,6,2: recurrence = direct multinomial")
print("target translation kills E22 with c=2*a22/13")
print("K=-(d/2)*(Phi/q)+T*H26; H26 terms =", len(sp.Poly(H, T, d, s, *seeds).terms()))
print("generic polar cases checked:", generic_cases)
print("generic cubic resultant =", sp.factor(sp.resultant(f, g, r)))
print("exceptional weights Phi/Psi/H =", (phi_weight, psi_weight, h_weight))
print("exceptional Groebner basis =", [p.as_expr() for p in exceptional_groebner.polys])
print("exceptional quintic resultant =", sp.factor(quintic_resultant))
print("regular regions: H bounds and unique 143*T^4/16384 pole checked")
print("power carrier: every root of x^3-1 has (v(V),v(M))=(3,1)")
print("Chebyshev carrier: y=+/-1 has (v(V),v(M))=(3,1)")
print("ALL EXACT CONTROLS PASS")
