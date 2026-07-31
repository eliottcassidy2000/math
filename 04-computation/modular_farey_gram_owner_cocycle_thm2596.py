#!/usr/bin/env python3
"""Exact companion for THM-2596.

Standard-library only.  Every truth-bearing check uses ``require`` so that
``python`` and ``python -O`` execute the same verification.
"""

from fractions import Fraction
from itertools import product


def require(condition, message):
    if not condition:
        raise RuntimeError("CHECK FAILED: " + message)


def mmul(a, b):
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(len(b)))
              for j in range(len(b[0])))
        for i in range(len(a))
    )


def mvec(a, x):
    return tuple(sum(a[i][j] * x[j] for j in range(len(x)))
                 for i in range(len(a)))


def trans(a):
    return tuple(zip(*a))


def madd(a, b):
    return tuple(tuple(a[i][j] + b[i][j] for j in range(len(a[0])))
                 for i in range(len(a)))


def msub(a, b):
    return tuple(tuple(a[i][j] - b[i][j] for j in range(len(a[0])))
                 for i in range(len(a)))


def smul(c, a):
    return tuple(tuple(c * x for x in row) for row in a)


def mpow(a, n):
    out = tuple(tuple(int(i == j) for j in range(len(a)))
                for i in range(len(a)))
    for _ in range(n):
        out = mmul(out, a)
    return out


def det2(a):
    return a[0][0] * a[1][1] - a[0][1] * a[1][0]


def inv2(a):
    d = det2(a)
    require(d in (-1, 1), "inverse requested outside GL_2(Z)")
    return ((a[1][1] // d, -a[0][1] // d),
            (-a[1][0] // d, a[0][0] // d))


def det3(a):
    return (
        a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
    )


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def defect(w, d):
    return dot(d, d) - 91 * dot(w, d)


I2 = ((1, 0), (0, 1))
NEG_I2 = ((-1, 0), (0, -1))
S = ((0, -1), (1, 0))
C = ((0, -1), (1, 1))
L = ((1, 1), (0, 1))
R = ((1, 0), (1, 1))

require(mpow(S, 2) == NEG_I2, "S^2=-I")
require(mpow(C, 3) == NEG_I2, "C^3=-I")
require(mmul(S, C) == smul(-1, L), "L=-SC")
require(mmul(S, inv2(C)) == R, "R=S C^{-1}")
print("PSL generators: S^2=C^3=1 projectively; L=-SC and R=SC^-1: PASS")


# Exhaust the passive Gram--owner covariance on a nontrivial integer bank.
gs = (S, C, L, R)
audit_count = 0
for u in product(range(-3, 4), repeat=2):
    for v in product(range(-3, 4), repeat=2):
        U = ((u[0], v[0]), (u[1], v[1]))
        if abs(det2(U)) != 1:
            continue
        G = mmul(trans(U), U)
        for w in product(range(-2, 3), repeat=2):
            ell = mvec(trans(U), w)
            for g in gs:
                Up = mmul(U, g)
                Gp = mmul(mmul(trans(g), G), g)
                ellp = mvec(trans(g), ell)
                require(Gp == mmul(trans(Up), Up), "Gram update")
                require(ellp == mvec(trans(Up), w), "owner update")
                gi = inv2(g)
                for z in product(range(-2, 3), repeat=2):
                    zp = mvec(gi, z)
                    d = mvec(U, z)
                    require(mvec(Up, zp) == d, "passive coordinate equality")
                    lhs = dot(z, mvec(G, z)) - 91 * dot(ell, z)
                    rhs = dot(zp, mvec(Gp, zp)) - 91 * dot(ellp, zp)
                    require(lhs == rhs == defect(w, d), "defect covariance")
                    audit_count += 1
print(f"Gram-owner passive covariance exact cases: {audit_count}")


# Active modular action is not Euclidean and therefore is not an LRC symmetry.
w = (1, 0)
d = (91, -45)
dp = mvec(C, d)
wp = mvec(trans(inv2(C)), w)
require((dp, wp) == ((45, 46), (1, 1)), "active hostile coordinates")
require(dot(w, d) == dot(wp, dp) == 91, "active hostile pairing")
require(defect(w, d) == 2025 and defect(wp, dp) == -4140,
        "active hostile changes safe to bad")
print("active C3 hostile: F_(1,0)(91,-45)=2025, F_(1,1)(45,46)=-4140")


# Same endpoint defects and same unimodular/acute combinatorics, opposite child.
u = (1, 0)
v0 = (1, 1)
v1 = (90, 1)
require(defect(w, u) == -90, "endpoint u")
require(defect(w, v0) == defect(w, v1) == -89, "matched endpoints")
require(abs(u[0] * v0[1] - u[1] * v0[0]) == 1, "first unimodular")
require(abs(u[0] * v1[1] - u[1] * v1[0]) == 1, "second unimodular")
require(dot(u, v0) > 0 and dot(u, v1) > 0, "both acute")
m0 = tuple(u[i] + v0[i] for i in range(2))
m1 = tuple(u[i] + v1[i] for i in range(2))
require(defect(w, m0) == -177 and defect(w, m1) == 1,
        "matched endpoints have opposite mediant status")
print("endpoint/V4/partial-cube hostile: (-90,-89) -> child defects -177 and +1")


# Berggren/Pythagorean ternary cross-section.
A = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
B = ((1, 2, 2), (2, 1, 2), (2, 2, 3))
D = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))
Q = ((1, 0, 0), (0, 1, 0), (0, 0, -1))
I3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

for name, M, expected_det in (("A", A, 1), ("B", B, -1), ("C", D, 1)):
    require(det3(M) == expected_det, f"{name} determinant")
    require(mmul(mmul(trans(M), Q), M) == Q, f"{name} Lorentz form")
    require(tuple(tuple(x % 2 for x in row) for row in M) == I3,
            f"{name} is identity mod 2")
    require(mpow(M, 3) != I3, f"{name} is not order three")

NA = msub(A, I3)
ND = msub(D, I3)
zero3 = ((0, 0, 0), (0, 0, 0), (0, 0, 0))
require(mpow(NA, 2) != zero3 and mpow(NA, 3) == zero3,
        "A nontrivial unipotent")
require(mpow(ND, 2) != zero3 and mpow(ND, 3) == zero3,
        "C nontrivial unipotent")
require(mmul(madd(B, I3), madd(msub(mpow(B, 2), smul(6, B)), I3)) == zero3,
        "B characteristic factor")
require(sum(B[i][i] for i in range(3)) == 5, "B has non-torsion trace")

mobius = (
    (A, lambda x: Fraction(1, 2 - x), lambda y: 2 - Fraction(1, y),
     Fraction(1, 2), Fraction(1, 1)),
    (B, lambda x: Fraction(1, 2 + x), lambda y: Fraction(1, y) - 2,
     Fraction(1, 3), Fraction(1, 2)),
    (D, lambda x: x / (2 * x + 1), lambda y: y / (1 - 2 * y),
     Fraction(0, 1), Fraction(1, 3)),
)
branch_cases = 0
for n in range(2, 41):
    for m in range(1, n):
        x = Fraction(m, n)
        triple = (n * n - m * m, 2 * m * n, n * n + m * m)
        parameter_children = (
            (n, 2 * n - m),
            (n, 2 * n + m),
            (m, n + 2 * m),
        )
        for (M, forward, inverse, lo, hi), (mp, np) in zip(mobius, parameter_children):
            y = forward(x)
            require(lo < y < hi, "disjoint ternary cylinder")
            require(inverse(y) == x, "piecewise inverse")
            child = (np * np - mp * mp, 2 * mp * np, np * np + mp * mp)
            require(mvec(M, triple) == child, "Berggren/parameter intertwining")
            branch_cases += 1
print(f"Berggren ternary cross-section exact rational cases: {branch_cases}")

farey_left = lambda x: x / (1 + x)
farey_right = lambda x: Fraction(1, 2 - x)
reflection = lambda x: 1 - x
prefix_cases = 0
for n in range(2, 101):
    for m in range(1, n):
        x = Fraction(m, n)
        require(mobius[0][1](x) == farey_right(x), "A is the right Farey branch")
        require(mobius[1][1](x) == farey_left(farey_right(reflection(x))),
                "B is the reflected middle Farey prefix")
        require(mobius[2][1](x) == farey_left(farey_left(x)),
                "C is the double-left Farey prefix")
        prefix_cases += 1
print(f"binary/ternary PGL prefix-code exact rational cases: {prefix_cases}")
print("all Berggren matrices are I mod 2 and none is order 3: V4/S3 torsor reading fails")
print("ALL CHECKS PASSED")
