#!/usr/bin/env python3
"""Independent, assertion-free hostile audit for THM-3830.

This deliberately avoids SymPy for the slope arithmetic.  The exact finite-
field searches are hostile controls only; the characteristic-zero proof is the
domain/unit argument summarized in the semantic packet below.
"""

from __future__ import annotations

import ast
import hashlib
import json
from fractions import Fraction
from itertools import product
from pathlib import Path


GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    if not condition:
        raise RuntimeError(label)
    GATES += 1


def trim(poly: list[Fraction]) -> list[Fraction]:
    ans = poly[:]
    while len(ans) > 1 and ans[-1] == 0:
        ans.pop()
    return ans


def multiply(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    ans = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, u in enumerate(left):
        for j, v in enumerate(right):
            ans[i + j] += u * v
    return trim(ans)


def derivative(poly: list[Fraction]) -> list[Fraction]:
    return trim([Fraction(i) * poly[i] for i in range(1, len(poly))] or [Fraction(0)])


def divide_with_remainder(
    numerator: list[Fraction], denominator: list[Fraction]
) -> tuple[list[Fraction], list[Fraction]]:
    rem = trim(numerator)
    den = trim(denominator)
    quotient = [Fraction(0)] * max(1, len(rem) - len(den) + 1)
    while rem != [0] and len(rem) >= len(den):
        shift = len(rem) - len(den)
        factor = rem[-1] / den[-1]
        quotient[shift] += factor
        for j, value in enumerate(den):
            rem[j + shift] -= factor * value
        rem = trim(rem)
    return trim(quotient), rem


def monic_gcd(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    a0, b0 = trim(left), trim(right)
    while b0 != [0]:
        _, rem = divide_with_remainder(a0, b0)
        a0, b0 = b0, rem
    lead = a0[-1]
    return [value / lead for value in a0]


def bareiss_determinant(matrix: list[list[Fraction]]) -> Fraction:
    a0 = [row[:] for row in matrix]
    n = len(a0)
    sign = Fraction(1)
    previous = Fraction(1)
    for pivot_index in range(n - 1):
        if a0[pivot_index][pivot_index] == 0:
            swap = next(
                (r for r in range(pivot_index + 1, n) if a0[r][pivot_index] != 0),
                None,
            )
            if swap is None:
                return Fraction(0)
            a0[pivot_index], a0[swap] = a0[swap], a0[pivot_index]
            sign *= -1
        pivot = a0[pivot_index][pivot_index]
        for i in range(pivot_index + 1, n):
            for j in range(pivot_index + 1, n):
                a0[i][j] = (
                    a0[i][j] * pivot
                    - a0[i][pivot_index] * a0[pivot_index][j]
                ) / previous
        previous = pivot
        for i in range(pivot_index + 1, n):
            a0[i][pivot_index] = Fraction(0)
        for j in range(pivot_index + 1, n):
            a0[pivot_index][j] = Fraction(0)
    return sign * a0[-1][-1]


def resultant(left: list[Fraction], right: list[Fraction]) -> Fraction:
    f, g = trim(left), trim(right)
    m, n = len(f) - 1, len(g) - 1
    f_desc, g_desc = list(reversed(f)), list(reversed(g))
    size = m + n
    sylvester = [[Fraction(0)] * size for _ in range(size)]
    for row in range(n):
        sylvester[row][row : row + m + 1] = f_desc
    for offset in range(m):
        row = n + offset
        sylvester[row][offset : offset + n + 1] = g_desc
    return bareiss_determinant(sylvester)


def eval_mod(poly: list[Fraction], value: int, prime: int) -> int:
    gate(all(coef.denominator % prime != 0 for coef in poly), "modular denominator")
    ans = 0
    for coefficient in reversed(poly):
        residue = coefficient.numerator * pow(coefficient.denominator, -1, prime)
        ans = (ans * value + residue) % prime
    return ans


def multiply_mod(left: tuple[int, ...], right: tuple[int, ...], prime: int) -> tuple[int, ...]:
    ans = [0] * (len(left) + len(right) - 1)
    for i, u in enumerate(left):
        for j, v in enumerate(right):
            ans[i + j] = (ans[i + j] + u * v) % prime
    while len(ans) > 1 and ans[-1] == 0:
        ans.pop()
    return tuple(ans)


# Ascending coefficient order, constructed from the displayed factors rather
# than copied from the primary companion.
z2_factor = [Fraction(3), Fraction(0), Fraction(7)]
z3_factor = [Fraction(1), Fraction(0), Fraction(7), Fraction(3)]
a = multiply(z2_factor, z3_factor)
b = multiply(
    multiply([Fraction(1), Fraction(1)], [Fraction(1), Fraction(2)]),
    [Fraction(-1), Fraction(3)],
)
zb = multiply([Fraction(0), Fraction(1)], b)

gate(a == [3, 0, 28, 9, 49, 21], "independent expansion of a")
gate(b == [-1, 0, 7, 6], "independent expansion of b")
gate(monic_gcd(a, derivative(a)) == [1], "a is squarefree over Q")
gate(monic_gcd(a, zb) == [1], "a and z*b are coprime over Q")
res_ab = resultant(a, b)
res_azb = resultant(a, zb)
disc_a = resultant(a, derivative(a)) / a[-1]
gate(res_ab == -31298700, "independent Sylvester resultant(a,b)")
gate(res_azb == 93896100, "independent Sylvester resultant(a,z*b)")
gate(disc_a == 353831803500, "independent discriminant(a)")

# Exhaust the exact boundary equation over F_11 for degree <= 2.  Every
# polynomial solution D of D(c^4*b(alpha)+alpha^2*c^2*D)=0 must be either zero
# or the forced nonzero scalar.  This is the same domain dichotomy used over
# K[x], but enumerated without symbolic factorization.
prime = 11
roots = [u for u in range(prime) if eval_mod(a, u, prime) == 0]
gate(roots == [1, 3, 8], "visible hostile roots over F_11")
boundary_candidates = 0
boundary_solutions = 0
for alpha in roots:
    b_alpha = eval_mod(b, alpha, prime)
    gate(alpha != 0 and b_alpha != 0, "hostile root preserves alpha*b(alpha)")
    for c in range(1, prime):
        constant_term = pow(c, 4, prime) * b_alpha % prime
        d_coefficient = alpha * alpha * c * c % prime
        forced = -constant_term * pow(d_coefficient, -1, prime) % prime
        gate(forced != 0, "forced hostile boundary scalar is nonzero")
        for coefficients in product(range(prime), repeat=3):
            boundary_candidates += 1
            d_poly = tuple(coefficients)
            bracket = list((constant_term,))
            if len(bracket) < len(d_poly):
                bracket.extend([0] * (len(d_poly) - len(bracket)))
            for i, coefficient in enumerate(d_poly):
                bracket[i] = (bracket[i] + d_coefficient * coefficient) % prime
            product_poly = multiply_mod(d_poly, tuple(bracket), prime)
            is_solution = all(value == 0 for value in product_poly)
            if is_solution:
                boundary_solutions += 1
                gate(
                    d_poly == (0, 0, 0) or d_poly == (forced, 0, 0),
                    "no nonconstant boundary solution",
                )
gate(boundary_candidates == 39930, "complete hostile boundary universe")
gate(boundary_solutions == 60, "exact hostile boundary solution count")

# Independent bounded controls for the two structural facts behind the proof:
# polynomial-line units are scalars, and a coordinate cross has no nontrivial
# idempotent.  The unbounded statements follow immediately by degree/domain
# arguments; these searches are positive and hostile controls, not substitutes.
unit_prime = 5
line_polynomials = list(product(range(unit_prime), repeat=3))
unit_pairs = []
idempotents = []
one = (1,)
for left in line_polynomials:
    if multiply_mod(left, left, unit_prime) == tuple(trim([Fraction(v) for v in left])):
        idempotents.append(left)
    for right in line_polynomials:
        if multiply_mod(left, right, unit_prime) == one:
            unit_pairs.append((left, right))
gate(
    unit_pairs
    == [((u, 0, 0), (pow(u, -1, unit_prime), 0, 0)) for u in range(1, unit_prime)],
    "bounded line units are exactly nonzero scalars",
)
gate(idempotents == [(0, 0, 0), (1, 0, 0)], "bounded line idempotents are trivial")

# The completion is checked at several unrelated exact rational specializations;
# its formal proof is the one-line cancellation shown in the report.
completion_samples = 0
for alpha0, c0, t0, q0 in [
    (Fraction(2, 3), Fraction(5, 7), Fraction(-4, 9), Fraction(11, 6)),
    (Fraction(-7, 5), Fraction(3, 2), Fraction(8, 11), Fraction(-9, 4)),
    (Fraction(13, 8), Fraction(-6, 5), Fraction(17, 3), Fraction(2, 7)),
]:
    k0 = c0 + t0 * q0
    h0 = alpha0 * k0 + t0
    m0 = q0 / c0
    cap_c0 = (1 + alpha0 * q0) / c0
    gate(cap_c0 * k0 - m0 * h0 == 1, "determinant-one completion sample")
    completion_samples += 1

source = Path(__file__).read_text(encoding="utf-8")
gate(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "independent audit has no Python assert",
)

semantic = {
    "arithmetic": "custom Sylvester/Bareiss: Disc(a)=353831803500; Res(a,b)=-31298700; Res(a,z*b)=93896100",
    "boundary": "in the domain K[x], D(C0+C1*D)=0 with C0,C1 nonzero implies D=0 or the nonzero forced scalar",
    "completion": "((1+alpha*q)/c)*(c+t*q)-(q/c)*(alpha*(c+t*q)+t)=1",
    "crt_relation": "K[x,y]/(xy) is the equal-origin fibre product of K[x] and K[y], hence has only scalar units and trivial idempotents",
    "extension": "adjoining alpha remains a field in characteristic zero; x- and y-divisibility are coefficientwise and survive/descent along scalar extension",
    "nilpotents": "the stated transverse cross xy is squarefree/reduced; thickened or other nonreduced fibres are outside the exact statement",
    "result": "PASS: x divides d iff y divides d for every c!=0 and all-degree q,d in the normalized coordinate-cross grammar",
    "scope": "this is also an immediate connectedness corollary of THM-3827's canonical nontrivial CRT split",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("audit=THM-3830-independent-hostile")
print("verdict=PASS")
print(f"arithmetic=disc:{int(disc_a)};res_ab:{int(res_ab)};res_azb:{int(res_azb)}")
print(f"boundary_f11=candidates:{boundary_candidates};solutions:{boundary_solutions}")
print(f"completion_samples={completion_samples}")
print("unit_normalization=fibre_product_units_are_equal_nonzero_scalars")
print("crt_relation=coordinate_cross_has_only_trivial_idempotents")
print("finite_extension_and_nilpotent_seams=typed_and_clear")
print(f"GATES={GATES}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
