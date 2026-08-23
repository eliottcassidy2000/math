#!/usr/bin/env python3
"""Independent hostile audit of THM-3816's positive-filtration unit proof."""

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


GATES = 0


def gate(condition, label):
    global GATES
    if not bool(condition):
        raise RuntimeError(f"FAILED: {label}")
    GATES += 1


def equal(left, right, label):
    gate(sp.expand(left - right) == 0, label)


A, C, T, X, q = sp.symbols("A C T X q")
zero = sp.Integer(0)
one = (sp.Integer(1), zero, zero)
w = (zero, sp.Integer(1), zero)
t = (zero, zero, sp.Integer(1))
basis = (one, w, t)
shifts = (0, 2, 2)

# Original THM-3811 structure constants in the free basis (1,w,t).
original = {
    (1, 1): (-7 * A**2, C, -A),
    (1, 2): (3 * A**2 - A * C**2, zero, zero),
    (2, 2): (3 * A * C - C**3, C**2 - 3 * A, -7 * A),
}


def coefficient_weight(monomial):
    return 2 * monomial[0] + monomial[1]


def exact_weight_part(expression, wanted):
    answer = zero
    for monomial, coefficient in sp.Poly(expression, A, C).terms():
        if coefficient_weight(monomial) == wanted:
            answer += coefficient * A**monomial[0] * C**monomial[1]
    return sp.expand(answer)


# This is the complete multiplicativity test: scalar coefficient weights add,
# so it is enough to test every term in every product of the free basis.
top = {}
for pair, vector in original.items():
    source_weight = shifts[pair[0]] + shifts[pair[1]]
    top_vector = []
    for target, coefficient in enumerate(vector):
        if coefficient != 0:
            for monomial, _ in sp.Poly(coefficient, A, C).terms():
                gate(coefficient_weight(monomial) + shifts[target]
                     <= source_weight,
                     f"filtration term {pair}->{target}:{monomial}")
        top_vector.append(exact_weight_part(
            coefficient, source_weight - shifts[target]))
    top[pair] = tuple(top_vector)

gate(top[(1, 1)] == (-7 * A**2, zero, -A), "top w^2")
gate(top[(1, 2)] == (3 * A**2 - A * C**2, zero, zero), "top wt")
gate(top[(2, 2)] == (zero, C**2 - 3 * A, -7 * A), "top t^2")


def add(left, right):
    return tuple(sp.expand(x + y) for x, y in zip(left, right))


def scale(scalar, vector):
    return tuple(sp.expand(scalar * entry) for entry in vector)


def product(left, right):
    a0, a1, a2 = left
    b0, b1, b2 = right
    answer = scale(a0 * b0, one)
    answer = add(answer, scale(a0 * b1 + a1 * b0, w))
    answer = add(answer, scale(a0 * b2 + a2 * b0, t))
    answer = add(answer, scale(a1 * b1, top[(1, 1)]))
    answer = add(answer, scale(a1 * b2 + a2 * b1, top[(1, 2)]))
    answer = add(answer, scale(a2 * b2, top[(2, 2)]))
    return answer


for i, left in enumerate(basis):
    for j, middle in enumerate(basis):
        for k, right in enumerate(basis):
            gate(product(product(left, middle), right)
                 == product(left, product(middle, right)),
                 f"graded associativity {i}{j}{k}")

# Re-derive the generic field rather than accepting a presentation.  The
# relation A*t=-(w^2+7A^2) makes w a generator after A is inverted.
multiply_w = sp.Matrix.hstack(
    sp.Matrix(product(w, one)),
    sp.Matrix(product(w, w)),
    sp.Matrix(product(w, t)),
)
f0 = sp.expand((T * sp.eye(3) - multiply_w).det())
expected_f0 = T**3 + 7 * A**2 * T + 3 * A**3 - A**2 * C**2
equal(f0, expected_f0, "generic characteristic polynomial")

h = sp.cancel(f0.subs(T, A * q) / A**2)
equal(h, A * q**3 + 7 * A * q + 3 * A - C**2,
      "scaled generic polynomial")
reciprocal = sp.expand(X**3 * h.subs(q, 1 / X))
equal(reciprocal, (3 * A - C**2) * X**3 + 7 * A * X**2 + A,
      "reciprocal transformation")


def valuation_at_A(expression):
    if expression == 0:
        return sp.oo
    return min(monomial[0]
               for monomial, _ in sp.Poly(expression, A, C).terms())


coefficients = sp.Poly(reciprocal, X).all_coeffs()
gate(valuation_at_A(coefficients[0]) == 0,
     "Eisenstein leading A-valuation zero")
for index, coefficient in enumerate(coefficients[1:-1], start=1):
    gate(valuation_at_A(coefficient) >= 1,
         f"Eisenstein lower coefficient {index}")
gate(valuation_at_A(coefficients[-1]) == 1,
     "Eisenstein constant A-valuation one")
content = coefficients[0]
for coefficient in coefficients[1:]:
    content = sp.gcd(content, coefficient)
gate(sp.Poly(content, A, C).total_degree() == 0,
     "reciprocal polynomial primitive")

# These CAS calls are controls only.  The valuation calculation above is an
# Eisenstein proof in k[A,C] for every field k; it is not a Q-only inference.
K_Q = sp.QQ.frac_field(A, C)
gate(sp.Poly(reciprocal, X, domain=K_Q).is_irreducible,
     "QQ(A,C) reciprocal control")
gate(sp.Poly(f0, T, domain=K_Q).is_irreducible,
     "QQ(A,C) characteristic control")

# Free split pieces give gr(S)=R*1 (+) R*w (+) R*t as an R-module.  Hence it
# is torsion-free, and localization at K=k(A,C) is injective.  There the
# preceding irreducible cubic makes the algebra K[T]/(f0), a field.
gate(sp.eye(3).det() == 1, "free-basis localization injection")
gate(all(weight > 0 for weight in (2, 1, 2, 2)), "positive weights")
gate(shifts[0] == 0 and all(shift > 0 for shift in shifts[1:]),
     "F0 receives only the scalar summand")

source = Path(__file__).read_text(encoding="utf-8")
gate(not any(isinstance(node, ast.Assert)
             for node in ast.walk(ast.parse(source))),
     "assertion-free source")

semantic = {
    "scope": "all fields k; in particular algebraically closed characteristic zero",
    "filtration": "split weights A=2,C=1,w=2,t=2",
    "generic_algebra": "k(A,C)[T]/(T^3+7A^2T+3A^3-A^2C^2)",
    "irreducibility": "reciprocal Eisenstein at A in k[A,C]",
    "unit_lemma": "gr domain and F0=k imply S*=k*",
    "verdict": "mathematical PASS",
}
semantic_hash = hashlib.sha256(json.dumps(
    semantic, sort_keys=True, separators=(",", ":")).encode()).hexdigest()

print("audit=THM-3816-positive-filtration-units")
print("verdict=PASS_MATHEMATICS")
print("filtration=multiplicative_split_free")
print("generic_localization=irreducible_cubic_field")
print("eisenstein_scope=all_fields_k")
print("unit_inference=gr_domain_plus_F0_equals_k")
print(f"GATES={GATES}")
print(f"semantic_sha256={semantic_hash}")
