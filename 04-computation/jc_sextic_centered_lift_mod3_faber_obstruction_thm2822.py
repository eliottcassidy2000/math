#!/usr/bin/env python3
"""Exact carrier and Faber-row controls for THM-2822."""

import ast
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


def carrier(name, variable, D, T_pole, E, A0, C):
    """Verify the complete THM-2817/2796 response packet."""
    F = sp.cancel(E**2 / D)
    V = sp.factor(4 * D * T_pole**2 / C**2)
    G = sp.cancel(C * E / (2 * D * T_pole))
    M = sp.factor(E * T_pole)
    response = sp.factor(V * G)
    scale = sp.cancel(4 / C**2)

    require(zero(E**2 - D - A0), f"{name}: square defect")
    require(zero(sp.diff(F, variable) - 2 * G), f"{name}: F'=2G")
    require(zero(F - V * G**2), f"{name}: F=VG^2")
    require(
        zero(2 * V * sp.diff(G, variable) + sp.diff(V, variable) * G - 2),
        f"{name}: response ODE",
    )
    require(zero(response - (2 / C) * M), f"{name}: response carrier")
    require(
        zero(V - scale * M**2 + scale * A0 * T_pole**2),
        f"{name}: degree-eight square defect",
    )
    require(zero(1 / (1 - F) + D / A0), f"{name}: composed denominator")

    # If R=E/sqrt(D) and U=(2/C)T_pole sqrt(D), then R'=1/U.
    radical_free_rprime = (
        sp.diff(E, variable) - E * sp.diff(D, variable) / (2 * D)
    )
    require(
        zero(radical_free_rprime - C / (2 * T_pole)),
        f"{name}: exact deck differential",
    )

    require(sp.degree(D, variable) == 6, f"{name}: deg D")
    require(sp.degree(T_pole, variable) == 4, f"{name}: deg T_pole")
    require(sp.degree(E, variable) == 3, f"{name}: deg E")
    require(sp.degree(V, variable) == 14, f"{name}: deg V")
    require(sp.degree(M, variable) == 7, f"{name}: deg M")
    require(sp.degree(response, variable) == 7, f"{name}: deg response")

    return {
        "F": F,
        "V": V,
        "G": G,
        "M": M,
        "response": response,
        "square_defect": sp.factor(V - scale * M**2),
        "denominator": sp.factor(1 / (1 - F)),
    }


x, y, z, t, q = sp.symbols("x y z t q")

power = carrier(
    "power",
    x,
    x**3 * (x**3 - 1),
    x * (x**3 - 1),
    x**3 - sp.Rational(1, 2),
    sp.Rational(1, 4),
    -sp.Rational(3, 2),
)

T3 = 4 * y**3 - 3 * y
chebyshev = carrier(
    "chebyshev",
    y,
    (y**2 - sp.Rational(1, 4)) ** 2 * (y**2 - 1),
    (y**2 - sp.Rational(1, 4)) * (y**2 - 1),
    T3 / 4,
    sp.Rational(1, 16),
    -sp.Rational(3, 8),
)

require(
    zero(power["F"] - (2 * x**3 - 1) ** 2 / ((2 * x**3 - 1) ** 2 - 1)),
    "power composition",
)
require(
    zero(chebyshev["F"] - T3**2 / (T3**2 - 1)),
    "Chebyshev composition",
)


def direct_faber_coefficient(j, index):
    alpha = sp.Rational(2 * j - 1, 2)
    m = 4 * j - 2
    series = sp.series((1 + q * t**3) ** alpha, t, 0, m + 4)
    return sp.expand(series.removeO()).coeff(t, index)


def expected_faber_coefficient(j, index):
    if index % 3:
        return sp.Integer(0)
    alpha = sp.Rational(2 * j - 1, 2)
    return sp.binomial(alpha, index // 3) * q ** (index // 3)


phi_support = []
psi_support = []
response_support = []
for j in range(1, 31):
    m = 4 * j - 2
    values = []
    for offset in (1, 2, 3):
        direct = direct_faber_coefficient(j, m + offset)
        expected = expected_faber_coefficient(j, m + offset)
        require(zero(direct - expected), f"Faber coefficient j={j}, offset={offset}")
        values.append(sp.factor(4 * direct))

    phi, psi, response_row = values
    require((phi != 0) == (j % 3 == 1), f"Phi support j={j}")
    require((psi != 0) == (j % 3 == 0), f"Psi support j={j}")
    require((response_row != 0) == (j % 3 == 2), f"R support j={j}")

    if phi != 0:
        exponent = sp.Poly(phi, q).degree()
        require(exponent == (4 * j - 1) // 3, f"Phi exponent j={j}")
        phi_support.append(j)
    if psi != 0:
        exponent = sp.Poly(psi, q).degree()
        require(exponent == 4 * j // 3, f"Psi exponent j={j}")
        psi_support.append(j)
    if response_row != 0:
        quotient = sp.cancel(response_row / q)
        exponent = sp.Poly(quotient, q).degree()
        require(exponent == (4 * j - 2) // 3, f"K exponent j={j}")
        require(exponent % 4 == 2, f"K mod-four exponent j={j}")
        response_support.append(j)

require(phi_support == list(range(1, 31, 3)), "Phi residue class")
require(psi_support == list(range(3, 31, 3)), "Psi residue class")
require(response_support == list(range(2, 31, 3)), "R residue class")

# Finite full-bank linear-algebra control: after Phi=0 and Psi constant,
# exactly the j=2 mod 3 rows survive, and K=R_Q/q is in q^2 Q[q^4].
for row_count in range(1, 19):
    kept_exponents = []
    for j in range(1, row_count + 1):
        if j % 3 == 2:
            m = 4 * j - 2
            response_row = 4 * direct_faber_coefficient(j, m + 3)
            exponent = sp.Poly(sp.cancel(response_row / q), q).degree()
            kept_exponents.append(exponent)
    require(all(exponent % 4 == 2 for exponent in kept_exponents),
            f"full-bank q^2 Q[q^4], R={row_count}")

sample_rows = {}
for j in (2, 5, 8):
    m = 4 * j - 2
    sample_rows[j] = sp.factor(
        4 * direct_faber_coefficient(j, m + 3) / q
    )


def centered_zero_section(variable, data):
    V = data["V"]
    A_src = data["response"]
    P = sp.expand(V**2 * z**4 + A_src * z)
    px0 = sp.expand(sp.diff(P, variable).subs(z, 0))
    pz0 = sp.expand(sp.diff(P, z).subs(z, 0))
    require(px0 == 0, "centered P_x zero section")
    require(zero(pz0 - A_src), "centered P_z zero section")
    require(sp.degree(A_src, variable) == 7, "nonconstant centered multiplier")


centered_zero_section(x, power)
centered_zero_section(y, chebyshev)

# The zero-section shortcut is genuinely centered.  In the general source
# H=Vz^2+Bz+C0, L=A_src z+E0, its two coefficients are as stated.
v0, v1, v2, b0, b1, c0, c1, e0, e1, a0, a1 = sp.symbols(
    "v0 v1 v2 b0 b1 c0 c1 e0 e1 a0 a1"
)
V_generic = v0 + v1 * x + v2 * x**2
B_generic = b0 + b1 * x
C0_generic = c0 + c1 * x
E0_generic = e0 + e1 * x
A_generic = a0 + a1 * x
H_generic = V_generic * z**2 + B_generic * z + C0_generic
P_generic = sp.expand(H_generic**2 + A_generic * z + E0_generic)
require(
    zero(
        sp.diff(P_generic, x).subs(z, 0)
        - (2 * C0_generic * sp.diff(C0_generic, x)
           + sp.diff(E0_generic, x))
    ),
    "general zero-section P_x",
)
require(
    zero(
        sp.diff(P_generic, z).subs(z, 0)
        - (2 * C0_generic * B_generic + A_generic)
    ),
    "general zero-section P_z",
)

# A nonzero lower term can make the necessary zero-section gcd equal to one.
boundary_px = sp.Integer(1)  # C0=0 and E0=x.
boundary_pz = power["response"]
require(sp.gcd(sp.Poly(boundary_px, x), sp.Poly(boundary_pz, x)).degree() == 0,
        "lower-term boundary gcd")

require(not has_asserts(__file__), "Python assert nodes are forbidden")

print("THM-2822 CENTERED SEXTIC LIFT AUDIT -- exact")
print("assert_nodes=0")
print("power degrees D/T/E/V/M/A_src = 6/4/3/14/7/7")
print("power A_src =", sp.factor(power["response"]))
print("power square defect =", power["square_defect"])
print("power 1/(1-F) =", power["denominator"])
print("Chebyshev degrees D/T/E/V/M/A_src = 6/4/3/14/7/7")
print("Chebyshev A_src =", sp.factor(chebyshev["response"]))
print("Chebyshev square defect =", chebyshev["square_defect"])
print("Chebyshev 1/(1-F) =", chebyshev["denominator"])
print("Faber direct/formula rows checked: j=1..30")
print("Phi support: j=1 mod 3; Psi support: j=0 mod 3")
print("R support: j=2 mod 3; R/q lies in q^2 Q[q^4]")
print("sample K_2 =", sample_rows[2])
print("sample K_5 =", sample_rows[5])
print("sample K_8 =", sample_rows[8])
print("centered zero section: P_x(x,0)=0, P_z(x,0)=A_src, deg A_src=7")
print("general boundary: gcd(2C0*C0'+E0', 2C0*B+A_src)=1 is only necessary")
print("lower-term hostile control: C0=0, E0=x makes that gcd equal to one")
print("ALL EXACT CONTROLS PASS")
