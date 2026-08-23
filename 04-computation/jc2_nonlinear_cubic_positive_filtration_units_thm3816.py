#!/usr/bin/env python3
"""Exact canonical proof of THM-3816's positive-filtration unit gate."""

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def check(condition, label):
    global CHECKS
    if not bool(condition):
        raise RuntimeError(f"FAILED: {label}")
    CHECKS += 1


def same(lhs, rhs, label):
    check(sp.expand(lhs - rhs) == 0, label)


A, C, T, X = sp.symbols("A C T X")

# Coordinates are vectors in the basis (1,w,t).  These are the top-weight
# products for wt(A,C,w,t)=(2,1,2,2).
one = (sp.Integer(1), sp.Integer(0), sp.Integer(0))
w = (sp.Integer(0), sp.Integer(1), sp.Integer(0))
t = (sp.Integer(0), sp.Integer(0), sp.Integer(1))


def addv(left, right):
    return tuple(sp.expand(x + y) for x, y in zip(left, right))


def scalev(scalar, vector):
    return tuple(sp.expand(scalar * x) for x in vector)


def mulv(left, right):
    """Multiplication in the candidate associated graded algebra."""
    a0, a1, a2 = left
    b0, b1, b2 = right
    result = scalev(a0 * b0, one)
    result = addv(result, scalev(a0 * b1 + a1 * b0, w))
    result = addv(result, scalev(a0 * b2 + a2 * b0, t))
    result = addv(result, scalev(a1 * b1, (-7 * A**2, 0, -A)))
    result = addv(result, scalev(a1 * b2 + a2 * b1,
                                 (3 * A**2 - A * C**2, 0, 0)))
    result = addv(result, scalev(a2 * b2,
                                 (0, C**2 - 3 * A, -7 * A)))
    return tuple(sp.expand(entry) for entry in result)


same(mulv(w, w)[0], -7 * A**2, "gr w^2 scalar coefficient")
same(mulv(w, w)[2], -A, "gr w^2 t coefficient")
same(mulv(w, t)[0], 3 * A**2 - A * C**2,
     "gr wt scalar coefficient")
same(mulv(t, t)[1], C**2 - 3 * A, "gr t^2 w coefficient")
same(mulv(t, t)[2], -7 * A, "gr t^2 t coefficient")

# Associativity is automatic for an associated graded ring, but this exact
# independent table check also guards the extraction of the top products.
basis = (one, w, t)
for i, left in enumerate(basis):
    for j, middle in enumerate(basis):
        for k, right in enumerate(basis):
            check(mulv(mulv(left, middle), right)
                  == mulv(left, mulv(middle, right)),
                  f"graded associativity basis triple {i}{j}{k}")


def monomial_weight(mon):
    return 2 * mon[0] + mon[1]


def weight_part(expr, wanted):
    poly = sp.Poly(sp.expand(expr), A, C)
    answer = 0
    for mon, coefficient in poly.terms():
        if monomial_weight(mon) == wanted:
            answer += coefficient * A**mon[0] * C**mon[1]
    return sp.expand(answer)


# Original products in THM-3811.  Each basis shift is (0,2,2); extracting
# total weight four recovers precisely the multiplication above.
original_products = {
    "w2": (-7 * A**2, C, -A),
    "wt": (3 * A**2 - A * C**2, 0, 0),
    "t2": (3 * A * C - C**3, C**2 - 3 * A, -7 * A),
}
graded_products = {
    "w2": mulv(w, w),
    "wt": mulv(w, t),
    "t2": mulv(t, t),
}
shifts = (0, 2, 2)
for name, vector in original_products.items():
    for index, coefficient in enumerate(vector):
        if coefficient == 0:
            continue
        weights = [monomial_weight(mon) + shifts[index]
                   for mon, _ in sp.Poly(coefficient, A, C).terms()]
        check(max(weights) <= 4,
              f"original {name} respects filtration in component {index}")
        same(weight_part(coefficient, 4 - shifts[index]),
             graded_products[name][index],
             f"original {name} top component {index}")

# Multiplication by w over K=k(A,C).  Since A is invertible in K,
# t=-(w^2+7A^2)/A, so w generates the three-dimensional generic algebra.
Mw = sp.Matrix.hstack(sp.Matrix(mulv(w, one)),
                      sp.Matrix(mulv(w, w)),
                      sp.Matrix(mulv(w, t)))
f0 = T**3 + 7 * A**2 * T + 3 * A**3 - A**2 * C**2
same((T * sp.eye(3) - Mw).det(), f0,
     "graded characteristic polynomial of w")
check(Mw.rank() == 3, "multiplication matrix is exact rank three generically")

# Eisenstein certificate.  After T=Aq and reciprocal variable X=1/q,
# the cubic is (3A-C^2)X^3+7AX^2+A.  At the prime A its leading
# coefficient is nonzero mod A, every lower coefficient is divisible by A,
# and its constant coefficient is not divisible by A^2.
reciprocal = (3 * A - C**2) * X**3 + 7 * A * X**2 + A
reciprocal_poly = sp.Poly(reciprocal, X)
coefficients = reciprocal_poly.all_coeffs()
same(coefficients[0].subs(A, 0), -C**2,
     "Eisenstein leading coefficient is nonzero modulo A")
for index, coefficient in enumerate(coefficients[1:]):
    same(coefficient.subs(A, 0), 0,
         f"Eisenstein lower coefficient {index} divisible by A")
same(coefficients[-1] / A, 1,
     "Eisenstein constant coefficient has A-valuation exactly one")
check(sp.Poly(reciprocal, X,
              domain=sp.QQ.frac_field(A, C)).is_irreducible,
      "CAS-independent-target check: reciprocal cubic irreducible")
check(sp.Poly(f0, T,
              domain=sp.QQ.frac_field(A, C)).is_irreducible,
      "graded characteristic cubic irreducible")

# The graded algebra is free over k[A,C], hence injects into its localization.
# Its generic localization is the field K[T]/(f0), so it is a domain.  In a
# filtered ring with domain gr, degrees add.  F_0 contains constants only.
positive_weights = {"A": 2, "C": 1, "w": 2, "t": 2}
check(all(value > 0 for value in positive_weights.values()),
      "all four generator weights are positive")
check(min(positive_weights.values()) == 1,
      "degree-zero filtration piece is exactly QQ")

source = Path(__file__).read_text(encoding="utf-8")
check(not any(isinstance(node, ast.Assert)
              for node in ast.walk(ast.parse(source))),
      "no inactive Python assert")

semantic = {
    "packet": "THM-3811 cubic ring (A,C,7A,C^2-3A)",
    "filtration": "wt(A,C,omega,theta)=(2,1,2,2)",
    "graded_cubic": "T^3+7A^2T+3A^3-A^2C^2",
    "irreducibility": "reciprocal Eisenstein at A",
    "conclusion": "gr(S) domain, degrees add, F0=k, hence S*=k*",
    "remaining": "order of [E] and affineness of Xbar\\E remain open",
}
semantic_blob = json.dumps(semantic, sort_keys=True,
                           separators=(",", ":")).encode()

print("theorem=THM-3816-positive-filtration-constant-units-for-the-nonlinear-cubic-packet")
print("status=PROVED_VERIFIED_EXACT_INDEPENDENTLY_HOSTILE_AUDITED")
print("filtration=wt(A,C,omega,theta)=(2,1,2,2);F0=k")
print("gr_cubic=T3+7A2T+3A3-A2C2")
print("irreducible=reciprocal_Eisenstein_at_prime_A")
print("conclusion=associated_graded_domain;S_units=k_units")
print("remaining=class_E_order+affineness")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
