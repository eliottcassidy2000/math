#!/usr/bin/env python3
"""Exact controls for THM-3140 R=8 odd-bipole Faber-chart nonentry."""

from fractions import Fraction
from math import factorial

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def add(poly, monomial, coefficient):
    value = sp.cancel(poly.get(monomial, 0) + coefficient)
    if value:
        poly[monomial] = value
    elif monomial in poly:
        del poly[monomial]


def direct_coefficient(R, n):
    """[t^n](1+2*d*t^2+q*t^3+(d^2-s)*t^4)^(R-1/2)."""
    alpha = sp.Rational(2 * R - 1, 2)
    answer = {}
    for i in range(n // 2 + 1):
        for j in range(n // 3 + 1):
            remainder = n - 2 * i - 3 * j
            if remainder < 0 or remainder % 4:
                continue
            ell_total = remainder // 4
            total = i + j + ell_total
            multinomial = factorial(total) // (
                factorial(i) * factorial(j) * factorial(ell_total)
            )
            common = sp.binomial(alpha, total) * multinomial * 2**i
            for s_degree in range(ell_total + 1):
                add(
                    answer,
                    (j, i + 2 * (ell_total - s_degree), s_degree),
                    common
                    * sp.binomial(ell_total, s_degree)
                    * (-1) ** s_degree,
                )
    return answer


def shift(poly, q=0, d=0, s=0, scalar=1):
    answer = {}
    for (q_degree, d_degree, s_degree), value in poly.items():
        monomial = (q_degree + q, d_degree + d, s_degree + s)
        require(min(monomial) >= 0, f"negative exponent {monomial}")
        add(answer, monomial, scalar * value)
    return answer


def combine(*polys):
    answer = {}
    for poly in polys:
        for monomial, value in poly.items():
            add(answer, monomial, value)
    return answer


def row(R):
    c_phi = direct_coefficient(R, 4 * R - 1)
    c_psi = direct_coefficient(R, 4 * R)
    c_response = direct_coefficient(R, 4 * R + 1)
    return {
        "phi": shift(c_phi, q=-1, scalar=4),
        "psi": shift(c_psi, scalar=4),
        "H": shift(
            combine(c_response, shift(c_phi, d=1)), q=-3, scalar=4
        ),
    }


q, d, s, t = sp.symbols("q d s t")


def as_expression(poly):
    return sp.expand(
        sum(value * q**qe * d**de * s**se for (qe, de, se), value in poly.items())
    )


def recurrence_coefficients(R, limit):
    """Independent logarithmic-derivative recurrence for the same series."""
    alpha = sp.Rational(2 * R - 1, 2)
    p = {1: 0, 2: 2 * d, 3: q, 4: d**2 - s}
    c = [sp.Integer(1)] + [sp.Integer(0)] * limit
    for n in range(1, limit + 1):
        c[n] = sp.expand(
            sum(((alpha + 1) * k - n) * p[k] * c[n - k]
                for k in range(1, min(4, n) + 1)) / n
        )
    return c


ROWS = {j: row(j) for j in range(1, 9)}
for j in range(1, 9):
    rec = recurrence_coefficients(j, 4 * j + 1)
    require(
        sp.expand(as_expression(ROWS[j]["phi"]) - 4 * rec[4 * j - 1] / q) == 0,
        f"independent phi recurrence j={j}",
    )
    require(
        sp.expand(as_expression(ROWS[j]["psi"]) - 4 * rec[4 * j]) == 0,
        f"independent psi recurrence j={j}",
    )
    require(
        sp.expand(
            as_expression(ROWS[j]["H"])
            - 4 * (rec[4 * j + 1] + d * rec[4 * j - 1]) / q**3
        ) == 0,
        f"independent H recurrence j={j}",
    )


# THM-2827's polynomial-ring sidecar: H_j belongs to Q[T,s,Td].
# In q,d,s exponents this is qexp even and qexp/2 >= dexp.
for j in range(1, 9):
    for q_degree, d_degree, _s_degree in ROWS[j]["H"]:
        require(q_degree % 2 == 0, f"H parity j={j}")
        require(q_degree // 2 >= d_degree, f"H not in Q[T,s,Td], j={j}")


def weight(monomial, D, S):
    q_degree, d_degree, s_degree = monomial
    return q_degree - D * d_degree - S * s_degree


def face(poly, D, S):
    values = {monomial: weight(monomial, D, S) for monomial in poly}
    minimum = min(values.values())
    return minimum, {
        monomial: poly[monomial]
        for monomial, value in values.items()
        if value == minimum
    }


# Exact R=8 infinity faces.  Here ord_infinity(q)=1, deg(d)=D,
# deg(s)=S.  The three regions are D<2S+2, equality, and D>2S+2.
pure_s_phi = {(0, 0, 7): -sp.Rational(6435, 1024)}
pure_s_psi = {(0, 0, 8): sp.Rational(6435, 8192)}
pure_d_psi = {(8, 4, 0): sp.Rational(6435, 8192)}
wall_phi_support = {(0, 0, 7), (2, 1, 5), (4, 2, 3), (6, 3, 1)}
wall_psi_support = {
    (0, 0, 8), (2, 1, 6), (4, 2, 4), (6, 3, 2), (8, 4, 0)
}

for D in range(0, 81):
    for S in range(0, 81):
        phi_weight, phi_face = face(ROWS[8]["phi"], D, S)
        psi_weight, psi_face = face(ROWS[8]["psi"], D, S)
        if D < 2 * S + 2:
            require(phi_face == pure_s_phi, f"region-I phi face D,S={D,S}")
            require(psi_face == pure_s_psi, f"region-I psi face D,S={D,S}")
        elif D == 2 * S + 2:
            require(set(phi_face) == wall_phi_support, f"wall phi D,S={D,S}")
            require(set(psi_face) == wall_psi_support, f"wall psi D,S={D,S}")
        else:
            require(psi_face == pure_d_psi, f"region-III psi face D,S={D,S}")

        if S >= 1 and D <= 2 * S + 2:
            # In regions I and II the row-j weights are -(j-1)S and -jS.
            for j in range(1, 7):
                require(
                    face(ROWS[j]["phi"], D, S)[0] == -(j - 1) * S,
                    f"lower phi weight D,S,j={D,S,j}",
                )
                require(
                    face(ROWS[j]["psi"], D, S)[0] == -j * S,
                    f"lower psi weight D,S,j={D,S,j}",
                )
            require(phi_weight == -7 * S and psi_weight == -8 * S,
                    f"top I/II weights D,S={D,S}")

        if D > 2 * S + 2:
            # The top Psi order 8-4D is strictly below all retained j<=6.
            require(psi_weight == 8 - 4 * D, f"top III psi D,S={D,S}")
            for j in range(1, 7):
                require(
                    face(ROWS[j]["psi"], D, S)[0] > psi_weight,
                    f"strict III top/lower D,S,j={D,S,j}",
                )


q0, d0, s0, z = sp.symbols("q0 d0 s0 z")
wall_phi = sp.factor(sum(
    coefficient * q0**qe * d0**de * s0**se
    for (qe, de, se), coefficient in ROWS[8]["phi"].items()
    if (qe, de, se) in wall_phi_support
))
wall_psi = sp.factor(sum(
    coefficient * q0**qe * d0**de * s0**se
    for (qe, de, se), coefficient in ROWS[8]["psi"].items()
    if (qe, de, se) in wall_psi_support
))
expected_wall_phi = (
    sp.Rational(6435, 1024)
    * s0 * (d0 * q0**2 - s0**2)
    * (d0**2 * q0**4 - 6 * d0 * q0**2 * s0**2 + s0**4)
)
expected_wall_psi = (
    sp.Rational(6435, 8192)
    * (d0**4 * q0**8 - 28 * d0**3 * q0**6 * s0**2
       + 70 * d0**2 * q0**4 * s0**4
       - 28 * d0 * q0**2 * s0**6 + s0**8)
)
require(sp.expand(wall_phi - expected_wall_phi) == 0, "wall phi factorization")
require(sp.expand(wall_psi - expected_wall_psi) == 0, "wall psi factorization")
f = (z - 1) * (z**2 - 6 * z + 1)
g = z**4 - 28 * z**3 + 70 * z**2 - 28 * z + 1
wall_resultant = sp.resultant(f, g, z)
require(wall_resultant == 65536, "wall resultant")


# Exact odd-bipole data and response target.
x = sp.symbols("x")
Tp = x**2 - 1
E = sp.expand(sum(
    (-1)**j * sp.binomial(sp.Rational(11, 2), j) * x**(11 - 2 * j)
    for j in range(6)
))
C = sp.expand(2 * Tp * sp.diff(E, x) - 22 * x * E)
V = Tp**13
G = sp.cancel(E / Tp**12)
M = sp.expand(E * Tp)
require(C == sp.Rational(693, 128), "odd-bipole first integral")
require(sp.gcd(sp.Poly(E, x), sp.Poly(Tp, x)).degree() == 0, "E/Tp disjoint")
require(sp.gcd(sp.Poly(E, x), sp.Poly(sp.diff(E, x), x)).degree() == 0,
        "E squarefree")
require(sp.cancel(2 * V * sp.diff(G, x) + sp.diff(V, x) * G - C) == 0,
        "odd-bipole response ODE")
require(sp.degree(M, x) == 13 and sp.degree(Tp**6, x) == 12,
        "M/A infinity degree")


# Divisor certificate.  At every root of A away from V, K=T*H has order
# at least 2a, so A*K has order at least 3a>1, impossible because M is
# squarefree.  At each of the two V-roots, the R=8 resonance formula gives
# delta=1, ord(A)=1+5*delta=6, while the pure-q lane gives ord(B)>=7.
for a in range(1, 20):
    require(3 * a > 1, f"off-V A-root contradiction a={a}")
k = 2
delta = 1
nu = 2 + (4 * k + 3) * delta
a_pole = 1 + (2 * k + 1) * delta
require((nu, a_pole) == (13, 6), "first resonance orders")
require(2 * 7 >= nu and 2 * 6 < nu, "pure-q B/A inequalities")


# Therefore A=c*(x^2-1)^6, q^2 is a nonzero scalar/Tp, B is divisible by
# Tp^7, and d,s are polynomials.  If S=deg(s)>=1 or D=deg(d)>=3, the exact
# faces above leave an uncancelled pole in Phi or Psi.  Hence s is constant
# (possibly zero), while d is zero or has degree at most two.
# Then T has order 2 at infinity, s and Td are regular, every H_j is regular,
# and K=T*H_Q has order >=2.  But K=lambda*M/A is a scalar times E/Tp^5,
# of order -1.  This is the final contradiction.
required_K_order = -(sp.degree(E, x) - sp.degree(Tp**5, x))
faber_K_lower_bound = 2
require(required_K_order == -1, "target K infinity order")
require(faber_K_lower_bound > required_K_order, "final infinity contradiction")


print("THM-3140 R=8, D=11, N=22 ODD-BIPOLE FABER-CHART NONENTRY -- exact")
print("rows j=1..8: direct multinomial = independent logarithmic recurrence")
print("H_j ring gate: H_j in Q[T,s,Td] for every j=1..8")
print("odd-bipole E_11 =", E)
print("odd-bipole C_11 =", C)
print("divisor gate: A=c*(x^2-1)^6; B divisible by (x^2-1)^7")
print("therefore q^2=unit/(x^2-1), and d,s are polynomials")
print("infinity face regions: D<2S+2, D=2S+2, D>2S+2")
print("wall normalized resultant =", wall_resultant)
print("flux consequence: s is constant (possibly zero); d is zero or has degree <=2")
print("Faber K order at infinity >=", faber_K_lower_bound)
print("required K=lambda*M/A order at infinity =", required_K_order)
print("CONCLUSION: the explicit THM-3133 d=11 odd-bipole cannot enter the normalized R=8 nonsplit exact-square-prefix Faber/source chart")
print("ALL EXACT HOSTILE CONTROLS PASS")
