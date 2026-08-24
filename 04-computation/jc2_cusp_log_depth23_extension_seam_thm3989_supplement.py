#!/usr/bin/env python3
"""Independent exact probe of the first THM-3989 depth 2:3 cell.

This is a research companion, not a JC(2) proof.  It independently checks
the p=0 cusp seam, the depth-one extension class, the four lowest 2:3
convolution integrations, and the forced p-jet in the canonical simple-base
slice.
"""

from __future__ import annotations

import hashlib
import itertools
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    """Assertion-free exact gate, retained by ``python -O``."""
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


def zero(expression: sp.Expr, message: str) -> None:
    gate(sp.factor(sp.cancel(expression)) == 0, message)


s, tau = sp.symbols("s tau")
x = s / tau
u = s**2 / tau
p = s**2 + tau
y = s * (s**2 + tau)


# -------------------------------------------------------------------------
# Panel I.  The exact p=0 seam.
# -------------------------------------------------------------------------

zero(x.subs(tau, -s**2) + 1 / s, "p=0 color: x=-1/s")
zero(u.subs(tau, -s**2) + 1, "p=0 color: u=-1")
zero(p.subs(tau, -s**2), "p=0 color: p=0")
zero(y.subs(tau, -s**2), "p=0 color: y=0")


def generator(a: int, b: int, c: int, e: int) -> sp.Expr:
    return sp.expand(x**a * u**b * p**c * y**e)


def coeff(poly: sp.Expr, variable: sp.Symbol, exponent: int) -> sp.Expr:
    return sp.expand(poly).coeff(variable, exponent)


def parity_sign(exponent: int) -> int:
    """Return (-1)^exponent without introducing floats for negative powers."""
    return -1 if exponent % 2 else 1


seam_cells = 0
for total in range(0, 9):
    for a, b, c, e in itertools.product(range(total + 1), repeat=4):
        if a + b + c + e > total:
            continue
        expression = generator(a, b, c, e)
        specialized = sp.expand(expression.subs(tau, -s**2))
        # Direct specialization must be a Laurent polynomial with no
        # positive s powers.  It is (-1)^(a+b) s^(-a) when no p or y
        # factor is present, and zero otherwise.
        expected_specialization = (
            parity_sign(a + b) * s ** (-a) if c + e == 0 else sp.Integer(0)
        )
        zero(
            specialized - expected_specialization,
            f"p=0 specialization ({a},{b},{c},{e})",
        )

        min_tau = -(a + b)
        max_tau = c + e - (a + b)
        for ell in range(1, 10):
            seam_sum = sp.Integer(0)
            for i in range(min_tau, max_tau + 1):
                fi = coeff(expression, tau, i)
                requested = ell - 2 * i
                if requested >= 0:
                    seam_sum += parity_sign(i) * coeff(fi, s, requested)
            zero(seam_sum, f"seam ell={ell}, monomial ({a},{b},{c},{e})")
            seam_cells += 1


# The first seam, isolated on generic finite Laurent coefficients.
g01, g13, g25, g37 = sp.symbols("g01 g13 g25 g37")
first_seam_symbolic = coeff(g01 * s, s, 1)
first_seam_symbolic -= coeff(g13 * s**3, s, 3)
first_seam_symbolic += coeff(g25 * s**5, s, 5)
first_seam_symbolic -= coeff(g37 * s**7, s, 7)
zero(first_seam_symbolic - (g01 - g13 + g25 - g37),
     "first seam has the expected alternating coefficient shape")


# -------------------------------------------------------------------------
# Panel II.  The depth-one extension class.
# -------------------------------------------------------------------------

# Canonical witnesses of all monomial symbols s^n in s*k[s].  Their
# tau^0 residue modulo k[s^2,s^3] is nonzero only at n=3.
theta_controls: dict[int, int] = {}
for n in range(1, 31):
    if n == 1:
        witness = x
    elif n == 2:
        witness = u
    elif n % 2 == 1:
        # n=1+2c, c>=1
        c = (n - 1) // 2
        witness = x * p**c
    else:
        # n=1+3+2c, c>=0
        c = (n - 4) // 2
        witness = x * p**c * y
    zero(coeff(witness, tau, -1) - s**n, f"P1 symbol s^{n}")
    residue = coeff(witness, tau, 0)
    theta = coeff(residue, s, 1)
    gate(theta == (1 if n == 3 else 0), f"P1 theta on s^{n}")
    theta_controls[n] = int(theta)


# Every polynomial with vanishing s coefficient lies in k[s^2,s^3].
q_coefficients = sp.symbols("q0:13")
q_no_s = sum(q_coefficients[n] * s**n for n in range(13) if n != 1)
P, Y = sp.symbols("P Y")
preimage = q_coefficients[0]
for n in range(2, 13):
    if n % 2 == 0:
        preimage += q_coefficients[n] * P ** (n // 2)
    else:
        preimage += q_coefficients[n] * P ** ((n - 3) // 2) * Y
zero(preimage.subs({P: s**2, Y: s**3}) - q_no_s,
     "k[s^2,s^3] is exactly the zero-s-coefficient subspace")


# Finite exact checks of the square approximate-root criterion.  Given h,
# prescribe q with [s]q=[s^3]h, build a P1 lift coefficient-by-coefficient,
# and verify that its square has the requested two negative coefficients.
for degree in range(1, 12):
    h_coeffs = [sp.Integer((7 * j + 3 * degree) % 11 - 5)
                for j in range(degree + 1)]
    h_coeffs[0] = sp.Integer(0)
    if all(value == 0 for value in h_coeffs[1:]):
        h_coeffs[1] = sp.Integer(1)
    h_poly = sum(h_coeffs[j] * s**j for j in range(degree + 1))
    lift = sp.Integer(0)
    for n in range(1, degree + 1):
        if h_coeffs[n] == 0:
            continue
        if n == 1:
            witness = x
        elif n == 2:
            witness = u
        elif n % 2 == 1:
            witness = x * p ** ((n - 1) // 2)
        else:
            witness = x * p ** ((n - 4) // 2) * y
        lift += h_coeffs[n] * witness
    lift = sp.expand(lift)
    lift_q = coeff(lift, tau, 0)
    zero(coeff(lift, tau, -1) - h_poly, f"constructed P1 lift degree {degree}")
    gate(coeff(lift_q, s, 1) == coeff(h_poly, s, 3),
         f"theta criterion degree {degree}")
    square = sp.expand(lift**2)
    zero(coeff(square, tau, -2) - h_poly**2,
         f"square leading coefficient degree {degree}")
    zero(coeff(square, tau, -1) - 2 * h_poly * lift_q,
         f"square next coefficient degree {degree}")


# The minimal hostile: h=s and q=s violate theta.  The tempting square
# packet A=x^2+2u has precisely that nonliftable next coefficient.
A_hostile = sp.expand(x**2 + 2 * u)
zero(coeff(A_hostile, tau, -2) - s**2, "hostile A leading h^2")
zero(coeff(A_hostile, tau, -1) - 2 * s**2, "hostile A next coefficient")
gate(coeff(s, s, 1) != coeff(s, s, 3),
     "h=s,q=s fail the P1 extension class")


# -------------------------------------------------------------------------
# Panel III.  The four lowest normalized 2:3 coefficient integrations.
# -------------------------------------------------------------------------

z = sp.symbols("z")
D = lambda expression: sp.diff(expression, z)
h, q, r, w, v = [sp.Function(name)(z) for name in ("h", "q", "r", "w", "v")]
a2 = sp.Function("a2")(z)
kappa = sp.symbols("kappa")

a_m2 = h**2
a_m1 = 2 * h * q
a_0 = q**2 + 2 * h * r
a_1 = 2 * h * w + 2 * q * r
c_m3 = h**3
c_m2 = 3 * h**2 * q
c_m1 = 3 * h * q**2 + 3 * h**2 * r
c_0 = q**3 + 6 * h * q * r + 3 * h**2 * w
c_1_principal = (
    sp.Rational(3, 2) * h * a2
    + sp.Rational(3, 2) * h * r**2
    + 3 * h * q * w
    + 3 * q**2 * r
)


def convolution(weight: int, aa: dict[int, sp.Expr], cc: dict[int, sp.Expr]) -> sp.Expr:
    return sp.expand(sum(
        j * D(ai) * cj - i * ai * D(cj)
        for i, ai in aa.items()
        for j, cj in cc.items()
        if i + j == weight
    ))


aa = {-2: a_m2, -1: a_m1, 0: a_0, 1: a_1, 2: a2}
cc = {-3: c_m3, -2: c_m2, -1: c_m1, 0: c_0,
      1: c_1_principal}
for weight in range(-5, 0):
    zero(convolution(weight, aa, cc), f"formal common-root row {weight}")

# Adding delta to c_1 leaves row -1 equal to 2h(h delta)'.  Hence the
# first free mismatch is delta=kappa/h.
delta = sp.Function("delta")(z)
cc_delta = dict(cc)
cc_delta[1] = c_1_principal + delta
zero(
    convolution(-1, aa, cc_delta) - 2 * h * D(h * delta),
    "row -1 mismatch is 2h(h delta)'",
)
zero((h * (kappa / h)) - kappa, "delta=kappa/h is the first mismatch")


# Check the integrated invariants themselves on independent generic
# coefficient functions, rather than only their common-root parametrization.
bfun, dfun, efun, ffun, gfun, c0fun, c1fun = [
    sp.Function(name)(z)
    for name in ("b", "d", "e", "f", "g", "c0", "c1")
]
raw_aa = {-2: h**2, -1: bfun, 0: ffun, 1: gfun, 2: a2}
raw_cc = {-3: h**3, -2: dfun, -1: efun, 0: c0fun, 1: c1fun}
zero(
    convolution(-4, raw_aa, raw_cc)
    - h**4 * D(2 * dfun / h**2 - 3 * bfun / h),
    "row -4 integrated invariant",
)


# -------------------------------------------------------------------------
# Panel IV.  Canonical simple-base slice and the forced p-jet.
# -------------------------------------------------------------------------

P0, Y0 = sp.symbols("P0 Y0")
lam, beta, K = sp.symbols("lam beta K", nonzero=True)
f_symbols: dict[tuple[int, int], sp.Symbol] = {}
F = sp.Integer(0)
for apow in range(5):
    for bpow in range(5 - apow):
        symbol = sp.symbols(f"f{apow}{bpow}")
        f_symbols[apow, bpow] = symbol
        F += symbol * P0**apow * Y0**bpow

F_log = sp.expand(F.subs({P0: s**2 + tau, Y0: s * (s**2 + tau)}))
phi = coeff(F_log, tau, 0)
psi = coeff(F_log, tau, 1)
chi = coeff(F_log, tau, 2)

slice_a = {
    -2: s**2,
    -1: lam * s**2,
    0: phi,
    1: psi,
    2: chi,
}
slice_c = {
    -3: beta * s**3,
    -2: sp.Rational(3, 2) * beta * lam * s**3,
    -1: (
        sp.Rational(3, 2) * beta * s
        * (phi + lam**2 * s**2 / 4)
        + K * s
    ),
    0: (
        beta * (
            sp.Rational(3, 2) * s * psi
            + sp.Rational(3, 4) * lam * s * phi
            - sp.Rational(1, 16) * lam**3 * s**3
        )
        + K * lam * s / 2
    ),
}

# The displayed coefficients solve rows -5 through -2 identically.
Ds = lambda expression: sp.diff(expression, s)


def s_convolution(weight: int, aa0: dict[int, sp.Expr],
                  cc0: dict[int, sp.Expr]) -> sp.Expr:
    return sp.expand(sum(
        j * Ds(ai) * cj - i * ai * Ds(cj)
        for i, ai in aa0.items()
        for j, cj in cc0.items()
        if i + j == weight
    ))


for weight in range(-5, -1):
    zero(s_convolution(weight, slice_a, slice_c),
         f"simple-base integrated row {weight}")

# Solve row -1 coefficientwise for c_1 through the range needed by the
# ell=3 seam.  The ODE has zero constant coefficient and determines [s]c_1.
c1_coefficients = sp.symbols("c1_0:8")
c1_series = sum(c1_coefficients[j] * s**j for j in range(8))
slice_c_with_c1 = dict(slice_c)
slice_c_with_c1[1] = c1_series
row_m1 = sp.expand(s_convolution(-1, slice_a, slice_c_with_c1))
zero(row_m1.coeff(s, 0), "row -1 has no s^0 equation")

solutions: dict[sp.Symbol, sp.Expr] = {}
for degree in range(1, 5):
    equation = sp.expand(row_m1.subs(solutions)).coeff(s, degree)
    variable = c1_coefficients[degree - 1]
    solved = sp.solve(equation, variable, dict=False)
    gate(len(solved) == 1, f"row -1 uniquely solves c1 degree {degree-1}")
    solutions[variable] = sp.factor(solved[0])
c1_solved = sp.expand(c1_series.subs(solutions))
gate(c1_solved.coeff(s, 0) == 0, "row -1 forces c1(0)=0")


def seam(expression_coefficients: dict[int, sp.Expr], ell: int) -> sp.Expr:
    total = sp.Integer(0)
    for index, polynomial in expression_coefficients.items():
        requested = ell - 2 * index
        if requested >= 0:
            total += parity_sign(index) * coeff(polynomial, s, requested)
    return sp.factor(total)


slice_c_with_c1[1] = c1_solved
seam_1 = seam(slice_c_with_c1, 1)
seam_2 = seam(slice_c_with_c1, 2)
seam_3 = seam(slice_c_with_c1, 3)
f00 = f_symbols[0, 0]
f10 = f_symbols[1, 0]
zero(
    seam_1 - lam * (-3 * beta * lam + 6 * beta * f00 + 4 * K) / 8,
    "simple-base ell=1 seam",
)
zero(seam_2, "simple-base ell=2 seam is automatic")
zero(
    seam_3
    + (beta * lam**3 - 3 * beta * lam**2 * f00
       + 12 * beta * f00 * f10 - 2 * K * lam**2
       + 8 * K * f10) / 16,
    "simple-base ell=3 seam",
)

K_tuned = sp.Rational(3, 4) * beta * lam - sp.Rational(3, 2) * beta * f00
zero(
    sp.factor(seam_3.subs(K, K_tuned))
    - beta * lam * (lam**2 - 12 * f10) / 32,
    "ell=1 tuning reduces ell=3 to F_p=lambda^2/12",
)

# Retraction control.  At ell=4 the positive tau^2 coefficient contributes
# +[s^0]c_2.  Therefore ell=4 cannot force F_y or any other new A-jet.
c2_0 = sp.symbols("c2_0")
slice_c_with_c2 = dict(slice_c_with_c1)
slice_c_with_c2[2] = c2_0
seam_4_without_c2 = seam(slice_c_with_c1, 4)
seam_4_with_c2 = seam(slice_c_with_c2, 4)
zero(seam_4_with_c2 - seam_4_without_c2 - c2_0,
     "ell=4 acquires the free +c2(0) term")
zero(seam_4_with_c2.subs(c2_0, -seam_4_without_c2),
     "c2(0) absorbs the entire ell=4 seam")


summary = {
    "checks": CHECKS,
    "seam_cells": seam_cells,
    "p0_seam": "sum_i (-1)^i [s^(ell-2i)] f_i = 0 for ell>=1",
    "P1_extension": "theta(h)=[s^3]h times [s] in k[s]/k[s^2,s^3]",
    "square_root_lift": (
        "h|a_-1 and [s](a_-1/(2h))=[s^3]h"
    ),
    "rows": {
        "-4": "(2c_-2/h^2-3a_-1/h)'=0",
        "-3": "formal root coefficient r after A translation",
        "-2": "formal root coefficient w after C translation",
        "-1": "c_1=common-root coefficient+kappa/h",
    },
    "simple_base": (
        "A=x^2+lambda*u+F(p,y), lambda!=0 => "
        "F_p(0,0)=lambda^2/12 if all negative bracket rows vanish"
    ),
    "retraction": "ell=4 has free +c_2(0); no F_y obstruction",
    "scope": "filtered 2:3 obstruction and stopping reason; JC2 open",
}
semantic = hashlib.sha256(
    json.dumps(summary, sort_keys=True).encode("utf-8")
).hexdigest()

print("THM-3989 independent depth-2:3 seam probe")
print(f"CHECKS={CHECKS}")
print(f"SEAM_CELLS={seam_cells}")
print("P0_SEAM=SUM_I_(-1)^I*[S^(ELL-2I)]F_I=0")
print("P1_EXTENSION=THETA(H)=[S^3]H*[S]")
print("SQUARE_ROOT_LIFT=H_DIVIDES_A_-1_AND_[S]Q=[S^3]H")
print("LOW_ROWS=-4,-3,-2_FORMAL_COMMON_ROOT;-1_FIRST_KAPPA/H_MISMATCH")
print("SIMPLE_BASE_FORCED_JET=F_P(0,0)=LAMBDA^2/12")
print("RETRACTION=ELL4_HAS_FREE_C2(0);NO_F_Y_OBSTRUCTION")
print("SCOPE=FILTERED_2_TO_3_STOPPING_REASON;JC2_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
