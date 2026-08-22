"""Exact certificate for the homogeneous alternating-sextic HFC(3) cell."""

from __future__ import annotations

import hashlib
import math
import sys
from fractions import Fraction
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))
import factorial_hfc3_symmetry_cells_support_thm3304 as sparse


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


A, B, C = sp.symbols("A B C")
PARAMETERS = (A, B, C)
MOMENTS = (2, 4, 6)
DENOMINATORS = {
    2: 25225200,
    4: 133989577722000,
    6: 6909970509968180054400,
}


def primitive_form(poly):
    cleared = sp.Poly(poly, *PARAMETERS, domain=sp.QQ).clear_denoms()[1]
    _, primitive = sp.Poly(cleared, *PARAMETERS, domain=sp.ZZ).primitive()
    return -primitive if primitive.LC() < 0 else primitive


def homogeneous_route():
    e1 = {(1, 0, 0): 1, (0, 1, 0): 1, (0, 0, 1): 1}
    e2 = {(1, 1, 0): 1, (1, 0, 1): 1, (0, 1, 1): 1}
    e3 = {(1, 1, 1): 1}
    delta = sparse.vandermonde()
    basis = [
        sparse.mul(delta, sparse.power(e1, 3)),
        sparse.mul(delta, sparse.mul(e1, e2)),
        sparse.mul(delta, e3),
    ]
    forms = {}
    raw = {}
    for exponent in MOMENTS:
        raw[exponent] = sparse.moment_form(basis, exponent)
        forms[exponent] = primitive_form(
            sparse.form_to_sympy(raw[exponent], PARAMETERS).as_expr())
    return basis, raw, forms


def mul2(left, right):
    out = {}
    for (i, j), av in left.items():
        for (k, ell), bv in right.items():
            key = (i + k, j + ell)
            out[key] = out.get(key, 0) + av * bv
    return {key: value for key, value in out.items() if value}


def pow2(poly, exponent):
    result = {(0, 0): 1}
    for _ in range(exponent):
        result = mul2(result, poly)
    return result


def integrate_triangle(poly):
    total = Fraction(0)
    for (i, j), coefficient in poly.items():
        total += Fraction(coefficient * math.factorial(i) * math.factorial(j),
                          math.factorial(i + j + 2))
    return total


def affine_route():
    x, y = sp.symbols("x y")
    z = 1 - x - y
    delta = (x - y) * (y - z) * (z - x)
    e2 = x*y + x*z + y*z
    e3 = x*y*z
    expressions = (delta, delta*e2, delta*e3)
    basis = []
    for expression in expressions:
        poly = sp.Poly(sp.expand(expression), x, y, domain=sp.ZZ)
        basis.append({powers: int(coefficient)
                      for powers, coefficient in poly.terms()})
    forms = {}
    for exponent in MOMENTS:
        powers = [[pow2(item, e) for e in range(exponent + 1)]
                  for item in basis]
        expression = 0
        for counts in sparse.weak_compositions(exponent, 3):
            product = {(0, 0): 1}
            for index, count in enumerate(counts):
                product = mul2(product, powers[index][count])
            coefficient = (sparse.multinomial(counts)
                           * integrate_triangle(product))
            monomial = sp.Rational(coefficient.numerator,
                                   coefficient.denominator)
            for variable, count in zip(PARAMETERS, counts):
                monomial *= variable**count
            expression += monomial
        forms[exponent] = sp.Poly(expression, *PARAMETERS, domain=sp.QQ)
    return basis, forms


def primitive_univariate(poly, variable):
    cleared = sp.Poly(poly, variable, domain=sp.QQ).clear_denoms()[1]
    _, primitive = sp.Poly(cleared, variable, domain=sp.ZZ).primitive()
    return -primitive if primitive.LC() < 0 else primitive


basis_homogeneous, raw_forms, FORMS = homogeneous_route()
basis_affine, affine_forms = affine_route()

require([len(item) for item in basis_homogeneous] == [18, 12, 6],
        "homogeneous supports")
require([len(item) for item in basis_affine] == [8, 16, 12],
        "affine supports")
require({m: len(FORMS[m].terms()) for m in MOMENTS} == {2: 6, 4: 15, 6: 28},
        "complete moment forms")

for exponent in MOMENTS:
    expected = sp.Poly(FORMS[exponent].as_expr() / DENOMINATORS[exponent],
                       *PARAMETERS, domain=sp.QQ)
    require(affine_forms[exponent] == expected,
            f"affine/homogeneous agreement at {exponent}")

# Symmetry controls: the normal form is alternating under a transposition and
# fixed by a three-cycle.
u, v, w = sp.symbols("u v w")
e1 = u + v + w
e2 = u*v + u*w + v*w
e3 = u*v*w
vandermonde = (u - v)*(v - w)*(w - u)
normal_form = sp.expand(vandermonde*(A*e1**3 + B*e1*e2 + C*e3))
require(sp.expand(normal_form.subs({u: v, v: u}, simultaneous=True)
                  + normal_form) == 0,
        "transposition sign")
require(sp.expand(normal_form.subs({u: v, v: w, w: u}, simultaneous=True)
                  - normal_form) == 0,
        "three-cycle invariance")

# Affine projective chart A=1.  Eliminate B from (P2,P4) and (P2,P6).
p2 = FORMS[2].as_expr().subs(A, 1)
p4 = FORMS[4].as_expr().subs(A, 1)
p6 = FORMS[6].as_expr().subs(A, 1)
require([sp.Poly(p, B).degree() for p in (p2, p4, p6)] == [2, 4, 6],
        "constant B-degrees on A=1")
R24 = primitive_univariate(sp.resultant(p2, p4, B), C)
R26 = primitive_univariate(sp.resultant(p2, p6, B), C)
require((R24.degree(), R26.degree()) == (8, 12),
        "affine resultant degrees")
factorization24 = sp.factor_list(R24.as_expr())[1]
require(len(factorization24) == 1
        and factorization24[0][1] == 1
        and sp.Poly(factorization24[0][0], C).degree() == 8,
        "I2/I4 hostile has a genuine degree-eight survivor field")

# Compact characteristic-zero coprimality certificate.  Full degree survives
# modulo 17, and the modular gcd is one; therefore the Q-gcd is one.
R24_17 = sp.Poly(R24.as_expr(), C, modulus=17)
R26_17 = sp.Poly(R26.as_expr(), C, modulus=17)
require((R24_17.degree(), R26_17.degree()) == (8, 12),
        "mod-17 degree preservation")
require((int(R24_17.LC()) % 17, int(R26_17.LC()) % 17) == (12, 9),
        "mod-17 leading coefficients")
require(sp.gcd(R24_17, R26_17).degree() == 0,
        "mod-17 coprimality")
require(sp.gcd(R24_17, R24_17).degree() == 8,
        "gcd machinery positive control")

# Independent exact Groebner route on the same affine chart.
groebner = sp.groebner((p2, p4, p6), B, C,
                       order="grevlex", method="f5b")
require(len(groebner.polys) == 1 and groebner.polys[0].as_expr() == 1,
        "affine Groebner unit ideal")

# Infinity line A=0.  B=1 covers its affine part; B=0,C=1 is the last point.
Q2 = sp.Poly(FORMS[2].as_expr().subs({A: 0, B: 1}), C, domain=sp.ZZ)
Q4 = sp.Poly(FORMS[4].as_expr().subs({A: 0, B: 1}), C, domain=sp.ZZ)
Q6 = sp.Poly(FORMS[6].as_expr().subs({A: 0, B: 1}), C, domain=sp.ZZ)
infinity_resultant = sp.resultant(Q2.as_expr(), Q4.as_expr(), C)
require(infinity_resultant == 14827264400,
        "infinity line excluded already by I2/I4")
require(sp.gcd(sp.gcd(Q2, Q4), Q6).degree() == 0,
        "infinity triple gcd")
require(tuple(FORMS[m].coeff_monomial(C**m) for m in MOMENTS) == (1, 1, 12),
        "last projective point")

def coefficient_hash(poly):
    payload = ";".join(
        ",".join(map(str, powers)) + "=" + str(int(coefficient))
        for powers, coefficient in poly.terms())
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


FORM_HASHES = {m: coefficient_hash(FORMS[m]) for m in MOMENTS}
RESULTANT_HASHES = {24: coefficient_hash(sp.Poly(R24.as_expr(), C)),
                    26: coefficient_hash(sp.Poly(R26.as_expr(), C))}

payload = "\n".join([
    "universe=homogeneous degree-6 S3-sign eigenspace on Delta_2",
    "normal_form=Vandermonde*(A*e1^3+B*e1*e2+C*e3)",
    "odd_moments=zero by transposition sign",
    "I2=(" + str(FORMS[2].as_expr()) + ")/" + str(DENOMINATORS[2]),
    "I4=(" + str(FORMS[4].as_expr()) + ")/" + str(DENOMINATORS[4]),
    "I6=(" + str(FORMS[6].as_expr()) + ")/" + str(DENOMINATORS[6]),
    "form_term_counts=(6,15,28)",
    "form_hashes=" + repr(FORM_HASHES),
    "A1_resultant_degrees=(8,12)",
    "A1_resultant_hashes=" + repr(RESULTANT_HASHES),
    "A1_mod17=(degrees 8,12; leading 12,9; gcd 1)",
    "A1_groebner=[1]",
    "A0B1_resultant_I2_I4=14827264400",
    "A0B0C1=(1,1,12)",
    "hostile=I2 and I4 alone have a genuine irreducible degree-8 survivor; I6 is load-bearing",
    "consequence=no nonzero alternating sextic kills I2,I4,I6",
    "routes=homogeneous factorial-Dirichlet + affine sparse triangle integration + resultant/mod17 + Groebner",
    "controls=transposition sign; C3 invariance; self-gcd degree 8; all projective charts",
]) + "\n"
print(payload, end="")
print("payload_sha256=" + hashlib.sha256(payload.encode("ascii")).hexdigest())
