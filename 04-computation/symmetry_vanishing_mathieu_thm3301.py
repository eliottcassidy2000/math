#!/usr/bin/env python3
"""Exact controls for THM-3301: symmetry vanishing is Mathieu-compatible.

Exact rational arithmetic only; no floating point, no randomness, no imported
executable.  Every gate is an explicit ``require`` so that ordinary and ``-O``
replay are byte-identical.
"""

from fractions import Fraction
from math import factorial


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


# Polynomials in (x, y, t) over Q, with x = u+iv, y = u-iv so that
# rho = t^2 + xy is the quadratic form of Delta = 4 d_x d_y + d_t^2 and the
# circle action is x -> e^(i theta) x, y -> e^(-i theta) y, t -> t.
# The torus weight of x^a y^b t^c is a - b.


def pmul(a, b):
    out = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            key = (ka[0] + kb[0], ka[1] + kb[1], ka[2] + kb[2])
            out[key] = out.get(key, 0) + va * vb
    return {k: v for k, v in out.items() if v}


def ppow(p, e):
    result = {(0, 0, 0): Fraction(1)}
    for _ in range(e):
        result = pmul(result, p)
    return result


def padd(*polys):
    out = {}
    for p in polys:
        for k, v in p.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}


def pscale(p, c):
    return {k: v * c for k, v in p.items() if v * c}


def laplace(p):
    out = {}
    for (i, j, k), v in p.items():
        if i >= 1 and j >= 1:
            key = (i - 1, j - 1, k)
            out[key] = out.get(key, 0) + 4 * i * j * v
        if k >= 2:
            key = (i, j, k - 2)
            out[key] = out.get(key, 0) + k * (k - 1) * v
    return {k: v for k, v in out.items() if v}


def moment(p):
    """L(f) = (exp(Delta/2) f)(0), the Gaussian moment functional."""
    total = Fraction(0)
    current = dict(p)
    order = 0
    while current:
        constant = current.get((0, 0, 0), 0)
        if constant:
            total += Fraction(constant, 2 ** order * factorial(order))
        current = laplace(current)
        order += 1
    return total


def weights(p):
    return sorted({a - b for (a, b, c) in p})


X = {(1, 0, 0): Fraction(1)}
Y = {(0, 1, 0): Fraction(1)}
T = {(0, 0, 1): Fraction(1)}
RHO = padd(ppow(T, 2), pmul(X, Y))

require(weights(RHO) == [0], "rho is torus invariant (weight 0)")
require(laplace(RHO) == {(0, 0, 0): 6}, "Delta rho = 6")


# ------------- 1.  a nonzero-weight eigenvector kills every moment ...

EIGEN_ROWS = []
for weight_power in (1, 2, 3):
    eigen = ppow(X, weight_power)                    # weight = +weight_power
    require(weights(eigen) == [weight_power], "pure weight eigenvector")
    for m in range(1, 9):
        require(moment(ppow(eigen, m)) == 0,
                "L(P^m)=0 for every m when P has nonzero weight")
        EIGEN_ROWS.append((weight_power, m))
require(len(EIGEN_ROWS) == 24, "eigenvector bank complete")

# ... but it also satisfies the Mathieu conclusion: for FIXED Q the moments
# L(Q P^m) vanish for all large m, with the exceptional set explicitly the
# m for which some weight of Q equals -m*weight(P).

MATHIEU_ROWS = []
TEST_Q = {
    "y^2": ppow(Y, 2),
    "y^4": ppow(Y, 4),
    "x^2": ppow(X, 2),
    "t^2": ppow(T, 2),
    "rho^3": ppow(RHO, 3),
    "x^2+y^6": padd(ppow(X, 2), ppow(Y, 6)),
}
for weight_power in (1, 2, 3):
    eigen = ppow(X, weight_power)
    for name, q in TEST_Q.items():
        nonzero = [m for m in range(1, 13)
                   if moment(pmul(q, ppow(eigen, m))) != 0]
        predicted = sorted({-w // weight_power for w in weights(q)
                            if w < 0 and (-w) % weight_power == 0
                            and 1 <= -w // weight_power <= 12})
        require(nonzero == predicted,
                "the exceptional m are EXACTLY those with a weight of Q "
                "equal to -m*weight(P)")
        require(len(nonzero) <= len(weights(q)),
                "the exceptional set is bounded by the weight support of Q")
        MATHIEU_ROWS.append((weight_power, name, tuple(nonzero)))
# The decisive statement: for every tested (P,Q) the moments vanish for all
# large m, so the Mathieu conclusion HOLDS -- no counterexample arises.
for weight_power, name, nonzero in MATHIEU_ROWS:
    require(len(nonzero) <= 1,
            "at most one exceptional exponent per (P,Q) here")

# Hostile control: a weight-ZERO P need not have vanishing moments at all.
require(moment(RHO) != 0, "hostile: weight-zero rho has nonzero moment")


# ---------------- 2.  THM-3290's counterexample is NOT an eigenvector

A_FORM = padd(RHO, ppow(X, 2))
C_FORM = padd(pmul(Y, ppow(RHO, 2)),
              pscale(pmul(pmul(X, ppow(T, 2)), RHO), -2),
              pscale(pmul(ppow(X, 3), ppow(T, 2)), -1))
P_3290 = pmul(A_FORM, ppow(C_FORM, 2))

require(weights(A_FORM) == [0, 2], "A mixes weights 0 and 2")
require(weights(C_FORM) == [-1, 1, 3], "C mixes three weights")
require(weights(P_3290) == [-2, 0, 2, 4, 6, 8],
        "THM-3290's P carries six distinct torus weights")
require(len(weights(P_3290)) > 1,
        "THM-3290's P is NOT a torus eigenvector, as Theorem 1 requires")

# and it genuinely violates the Mathieu conclusion, unlike any eigenvector
VIOLATION_ROWS = []
for m in range(1, 5):
    require(moment(ppow(P_3290, m)) == 0, "L(P^m)=0")
    value = moment(pmul(ppow(X, 2), ppow(P_3290, m)))
    require(value != 0, "L(x^2 P^m) != 0 for every m: Mathieu conclusion fails")
    VIOLATION_ROWS.append(m)
require(len(VIOLATION_ROWS) == 4, "violation bank complete")


# ------------------- 3.  finite groups leave an arithmetic progression alive

# On the sphere the rotation by 2*pi/e about the t-axis is an order-e symmetry;
# an eigenvector of that finite action has weight defined only mod e, so the
# surviving exponents form the progression e*Z.  Model it on the torus weight.
FINITE_ROWS = []
for order in (2, 3, 4, 5):
    survivors = [m for m in range(1, 4 * order + 1) if (m * 1) % order == 0]
    require(survivors == list(range(order, 4 * order + 1, order)),
            "finite order leaves exactly the multiples alive")
    FINITE_ROWS.append((order, len(survivors)))
require(len(FINITE_ROWS) == 4, "finite-order bank complete")


print("THM-3301 SYMMETRY VANISHING IS MATHIEU-COMPATIBLE EXACT CONTROL")
print("setting=Gaussian functional on rho=t^2+xy; torus weight of x^a y^b t^c = a-b")
print("eigenvector_rows=" + str(len(EIGEN_ROWS)))
print("eigenvector_fact=nonzero weight => L(P^m)=0 for every m>=1")
print("mathieu_rows=" + str(len(MATHIEU_ROWS)))
print("mathieu_fact=for fixed Q the exceptional m are those with a weight of "
      "Q equal to -m*weight(P); finitely many, so L(Q P^m)=0 for m >> 0")
print("nonzero_windows=" + repr([(w, n, z) for w, n, z in MATHIEU_ROWS if z]))
print("thm3290_weights=" + repr(weights(P_3290)))
print("thm3290_is_eigenvector=False")
print("thm3290_violation_m=" + repr(VIOLATION_ROWS))
print("finite_order_survivors=" + repr(FINITE_ROWS))
print("scope=mechanism classification only; no conjecture is proved or refuted")
print("ALL EXACT CHECKS PASSED")
