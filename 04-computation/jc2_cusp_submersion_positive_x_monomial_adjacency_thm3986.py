#!/usr/bin/env python3
"""Exact companion for THM-3986's height-two cusp adjacency wall."""

from __future__ import annotations

import hashlib
import json
from math import gcd

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    """Assertion-free exact gate, retained under ``python -O``."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(sp.cancel(expression)) == 0, message)


def jacobian(f: sp.Expr, g: sp.Expr, first: sp.Symbol,
             second: sp.Symbol) -> sp.Expr:
    return sp.expand(
        sp.diff(f, first) * sp.diff(g, second)
        - sp.diff(f, second) * sp.diff(g, first)
    )


x, t, u = sp.symbols("x t u")
alpha, gamma, lam = sp.symbols("alpha gamma lam", nonzero=True)
M, rr, ss = sp.symbols("M rr ss", integer=True)
rho = sp.expand(u * (u + 1))
sigma = sp.expand(u**2 * (u + 1))


# Panel I: freeze the height-two tower and every positive unperturbed row.
u_source = x**2 * t
p_source = sp.expand(t * (1 + u_source))
y_source = sp.expand(x * t**2 * (1 + u_source))
zero(jacobian(p_source, y_source, x, t) + t * p_source,
     "J(p,y)=-tp")
for m_value in range(2, 11):
    base = sp.expand(alpha * p_source + gamma * y_source**m_value)
    zero(jacobian(base, y_source, x, t) + alpha * t * p_source,
         f"m={m_value}: base Hamiltonian row")
    gate(sp.diff(base, t).subs(t, 0) == alpha,
         f"m={m_value}: t=0 source control")


# Panel II: the universal logarithmic kernel and compatibility exponents.
# In the x!=0 chart, P=alpha*p, V=gamma*y^m and
# W=lambda*x*p^r*y^s.  Put w=2r+3s-1.
w = 2 * rr + 3 * ss - 1
N = 3 * M - 2
L = 2 * rr + 3 * ss - 3
row_u_P = 2 * u + 1
row_u_V = M * (3 * u + 2)
row_u_W = (w + 1) * u + rr + 2 * ss
kernel = (
    sp.expand(3 * M * row_u_W - w * row_u_V),
    sp.expand(w * row_u_P - 2 * row_u_W),
    sp.expand(2 * row_u_V - 3 * M * row_u_P),
)
expected_kernel = (
    M * (3 * u + 2 - rr),
    -(2 * u + ss + 1),
    M,
)
for index, (actual, expected) in enumerate(zip(kernel, expected_kernel)):
    zero(actual - expected, f"universal kernel component {index}")

raw_exponents = (
    sp.expand((2 * M - 1) * L - N * (rr + 2 * ss - 1)),
    sp.expand((M - 1) * L - N * (rr + ss - 1)),
    sp.expand(L - N),
    sp.expand(-L),
)
expected_exponents = (
    M * rr + ss - 3 * M + 1,
    1 - M * rr - ss,
    2 * rr + 3 * ss - 3 * M - 1,
    -L,
)
for index, (actual, expected) in enumerate(
        zip(raw_exponents, expected_exponents)):
    zero(actual - expected, f"compatibility exponent {index}")
zero(sum(raw_exponents) + 2 * N,
     "infinity exponent is -2N/delta")


# Panel III: all L>0 invalid-address orders on a hostile exponent grid.
# The proof is symbolic; the grid independently checks divisibility and the
# complete list of regular forbidden addresses for many values of m.
grid_controls: dict[str, dict[str, object]] = {}
for m_value in range(2, 13):
    n_value = 3 * m_value - 2
    max_exp = 3 * m_value + 4
    regular_zero: set[tuple[int, int]] = set()
    regular_A: set[tuple[int, int]] = set()
    for r_value in range(0, max_exp + 1):
        for s_value in range(0, max_exp + 1):
            l_value = 2 * r_value + 3 * s_value - 3
            if l_value <= 0:
                continue
            delta = gcd(n_value, l_value)
            numerators = (
                m_value * r_value + s_value - 3 * m_value + 1,
                1 - m_value * r_value - s_value,
                2 * r_value + 3 * s_value - 3 * m_value - 1,
                -l_value,
            )
            gate(all(value % delta == 0 for value in numerators),
                 f"m={m_value},r={r_value},s={s_value}: integral powers")
            e0, e1, eA, eB = (value // delta for value in numerators)
            gate(e0 + e1 + eA + eB == -2 * n_value // delta,
                 f"m={m_value},r={r_value},s={s_value}: infinity exponent")

            # A=3u+2-r collides with u only for r=2.  B=2u+s+1
            # collides with u+1 only for s=1.  A and B never collide.
            order_zero = e0 + (eA if r_value == 2 else 0)
            order_minus_one = e1 + (eB if s_value == 1 else 0)
            if order_zero == 0:
                regular_zero.add((r_value, s_value))
            if eA == 0:
                regular_A.add((r_value, s_value))
            gate(order_minus_one != 0,
                 f"m={m_value},r={r_value},s={s_value}: -1 excluded")
            gate(eB < 0,
                 f"m={m_value},r={r_value},s={s_value}: B-root pole")
            gate(2 * r_value + 3 * s_value != 1,
                 f"m={m_value},r={r_value},s={s_value}: A/B separated")

    expected_zero = {
        (2, m_value - 1),
        (1, 2 * m_value - 1),
        (0, 3 * m_value - 1),
    }
    expected_A = {
        (r_value, s_value)
        for r_value in range(0, max_exp + 1)
        for s_value in range(0, max_exp + 1)
        if 2 * r_value + 3 * s_value == 3 * m_value + 1
    }
    gate(regular_zero == expected_zero,
         f"m={m_value}: exactly three regular u=0 rows")
    gate(regular_A == expected_A,
         f"m={m_value}: regular A-roots are exactly L=N")
    grid_controls[str(m_value)] = {
        "max_r_s": max_exp,
        "regular_zero": sorted(expected_zero),
        "regular_A": sorted(expected_A),
    }


# Panel IV: uniform repairs of all regular forbidden fibres.
# At u=0, after absorbing a nonzero scalar into kappa, the three equations
# have degrees 2, 3, 4.  The tuned root u=0 is simple for every M>=2.
repair_zero_2 = sp.expand((u + 1) * (2 * u + M) - M)
zero(repair_zero_2 - u * (2 * u + M + 2),
     "r=2,s=m-1 regular-zero repair")
gate(sp.diff(repair_zero_2, u).subs(u, 0) == M + 2,
     "r=2,s=m-1 tuned root is simple")

repair_zero_1 = sp.expand(
    (u + 1) * (2 * (u + M))**2 - 4 * M**2 * (3 * u + 1)
)
gate(sp.factor(repair_zero_1).subs(u, 0) == 0,
     "r=1,s=2m-1 tuned root")
gate(sp.factor(sp.diff(repair_zero_1, u).subs(u, 0)
               - 8 * M * (1 - M)) == 0,
     "r=1,s=2m-1 tuned root is simple")
gate(sp.degree(repair_zero_1, u) == 3,
     "r=1,s=2m-1 repair degree")

repair_zero_0 = sp.expand(
    (u + 1) * (2 * u + 3 * M)**3
    - sp.Rational(27, 4) * M**3 * (3 * u + 2)**2
)
gate(sp.factor(repair_zero_0).subs(u, 0) == 0,
     "r=0,s=3m-1 tuned root")
gate(sp.factor(sp.diff(repair_zero_0, u).subs(u, 0)
               - 54 * M**2 * (1 - M)) == 0,
     "r=0,s=3m-1 tuned root is simple")
gate(sp.degree(repair_zero_0, u) == 4,
     "r=0,s=3m-1 repair degree")

# If A=3u+2-r has exponent zero, then L=N.  Write
# r=3q+2, s=M-1-2q.  The forbidden root is u=q.  Its logarithmic
# derivative is -2/M, so a tuned fibre has many other roots.
q = sp.symbols("q", integer=True, positive=True)
s_A = M - 1 - 2 * q
log_derivative_A = sp.cancel(
    q / u - (q + 1) / (u + 1) - 2 / (2 * u + s_A + 1)
)
zero(log_derivative_A.subs(u, q) + sp.Rational(2, 1) / M,
     "regular A-root is simple")
for m_value in range(3, 13):
    for q_value in range(1, (m_value - 1) // 2 + 1):
        s_value = m_value - 1 - 2 * q_value
        tuned = sp.Rational(q_value**q_value,
                            (q_value + 1)**(q_value + 1) * m_value)
        fibre = sp.expand(
            u**q_value
            - tuned * (u + 1)**(q_value + 1)
            * (2 * u + s_value + 1)
        )
        gate(sp.degree(fibre, u) == q_value + 2,
             f"m={m_value},q={q_value}: A-repair fibre degree")
        gate(fibre.subs(u, q_value) == 0,
             f"m={m_value},q={q_value}: tuned A-root")
        gate(sp.diff(fibre, u).subs(u, q_value) != 0,
             f"m={m_value},q={q_value}: tuned A-root simple")
        for invalid in (0, -1, -sp.Rational(s_value + 1, 2)):
            gate(fibre.subs(u, invalid) != 0,
                 f"m={m_value},q={q_value}: other invalid {invalid}")


# Panel V: the three L<=0 rows are x, xp and xy.
# For x, gcd(N,3)=1 and the common-power gate has factor exponents
# (3m-1,-1,3m+1,-3).  Every forbidden point is a strict zero or pole.
x_gate_exponents = (3 * M - 1, -1, 3 * M + 1, -3)
zero(sum(x_gate_exponents) - 2 * N,
     "x-row infinity exponent is +2N")
for m_value in range(2, 21):
    n_value = 3 * m_value - 2
    gate(gcd(n_value, 3) == 1, f"m={m_value}: x-row coprime powers")
    gate(all(value != 0 for value in (
        3 * m_value - 1, -1, 3 * m_value + 1, -3
    )), f"m={m_value}: x-row excludes all finite invalid addresses")

# xp has an explicit critical point at both colors, x=-alpha/lambda.
p_laurent = rho / x**2
y_laurent = sigma / x**3
x_xp = -alpha / lam
for m_value in range(2, 9):
    xp_row = sp.expand(
        alpha * p_laurent + gamma * y_laurent**m_value
        + lam * x * p_laurent
    )
    for color in (0, -1):
        zero(sp.diff(xp_row, x).subs({x: x_xp, u: color}),
             f"m={m_value}: xp x-row at color {color}")
        zero(sp.diff(xp_row, u).subs({x: x_xp, u: color}),
             f"m={m_value}: xp u-row at color {color}")

# xy combines with p.  Its compatibility address is independent of m.
f_xy = sp.expand(rho * (alpha + lam * u))
xy_address_full = sp.factor(
    3 * sp.diff(f_xy, u) * sigma
    - 2 * f_xy * sp.diff(sigma, u)
)
xy_address = alpha - 2 * lam * u - 3 * lam * u**2
zero(xy_address_full + u**2 * (u + 1) * xy_address,
     "xy-row quadratic address")
gate(xy_address.subs(u, 0) == alpha, "xy-row avoids u=0")
gate(sp.factor(xy_address.subs(u, -1) - (alpha - lam)) == 0,
     "xy-row u=-1 value")
zero(xy_address.subs({lam: alpha, u: sp.Rational(1, 3)}),
     "xy-row tuned off-color repair")


# Panel VI: independent binomial-resultant controls for the common-x gate.
XN, XL = sp.symbols("XN XL", nonzero=True)
power_controls: dict[str, str] = {}
for n_value, l_value in ((4, 1), (4, 2), (4, 4),
                         (7, 3), (7, 7), (10, 6), (13, 9)):
    delta = gcd(n_value, l_value)
    resultant = sp.factor(sp.resultant(
        x**n_value - XN, x**l_value - XL, x
    ))
    relation = sp.expand(
        XN**(l_value // delta) - XL**(n_value // delta)
    )
    # The resultant is a nonzero sign times relation^delta.
    quotient = sp.factor(resultant / relation**delta)
    gate(quotient in (1, -1),
         f"N={n_value},L={l_value}: binomial compatibility resultant")
    power_controls[f"N{n_value}_L{l_value}"] = str(quotient)


summary = {
    "checks": CHECKS,
    "family": "alpha*p+gamma*y^m+lambda*x*p^r*y^s",
    "classification": (
        "base is a submersion; every nonzero one-monomial x-adjacency "
        "is critical for m>=2,r,s>=0"
    ),
    "kernel": "P:V:W=m(3u+2-r):-(2u+s+1):m",
    "positive_L": {
        "infinity_exponent": "-2(3m-2)/gcd(3m-2,L)",
        "regular_zero": ["(2,m-1)", "(1,2m-1)", "(0,3m-1)"],
        "regular_A": "L=3m-2; tuned A-root simple",
        "grid": grid_controls,
    },
    "low_L": {
        "x": "coprime torus gate; all forbidden addresses strict",
        "xp": "two explicit color criticals x=-alpha/lambda",
        "xy": "address alpha-2lambda*u-3lambda*u^2",
    },
    "power_controls": power_controls,
    "scope": "one positive-x monomial correction; rational mates not excluded",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3986 cusp-submersion positive-x monomial adjacency companion")
print(f"CHECKS={CHECKS}")
print("FAMILY=alpha*p+gamma*y^m+lambda*x*p^r*y^s;m>=2;r,s>=0")
print("BASE=AFFINE_SUBMERSION;EVERY_NONZERO_ADJACENCY=AFFINE_CRITICAL")
print("KERNEL=P:V:W=m(3u+2-r):-(2u+s+1):m")
print("LOW_L=x_TORUS_GATE;xp_TWO_COLOR;xy_QUADRATIC")
print("REGULAR_U0=THREE_FAMILIES_REPAIRED;REGULAR_A=L=N_SIMPLE_REPAIR")
print("SCOPE=ONE_POSITIVE_X_MONOMIAL;REGULAR_MATES_ONLY;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
