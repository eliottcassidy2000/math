#!/usr/bin/env python3
"""Exact finite-sheet gate for the fixed Keller polynomial R_5.

The proved fixed-map normalizations are

    H   = 2^6 L N(L),
    J   = 2^35 L^7 N(H),
    G   = L^43 N(J),
    R_5 = L^271 N(G).

At the canonical finite inverse point q=(2,5/6,-7/8) above the old
boundary L=0, this probe evaluates R_5(q) by a tower of four cubic finite
etale algebras.  The primary route descends all the way to the five-term
polynomial L.  A disjoint route stops one level earlier and evaluates the
frozen 361-term polynomial H directly.  Thus no global expansion of J, G,
or R_5 is constructed.

Every division is guarded by a regular-representation unit test, every
cubic discriminant is a unit, and every universal inverse point is checked
by direct substitution into F.  Nonzero good reductions prove the rational
value R_5(q) is nonzero.  Consequently the finite old-L sheet is a generic
unit for R_5.

Scope: one fixed map, one named finite sheet, and the next old-L valuation.
This is not a renewal-face, irreducibility, fifth-image, degree-243, or
general Jacobian-conjecture computation.  ``require`` remains active under
``python -O``.
"""

from __future__ import annotations

import hashlib
import pickle
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def determinant_mod(matrix: list[list[int]], prime: int) -> int:
    """Determinant over F_p by fraction-free row elimination."""

    work = [[entry % prime for entry in row] for row in matrix]
    size = len(work)
    det = 1
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            det = -det
        pivot_value = work[column][column] % prime
        det = det * pivot_value % prime
        pivot_inverse = pow(pivot_value, prime - 2, prime)
        for row in range(column + 1, size):
            if work[row][column] == 0:
                continue
            factor = work[row][column] * pivot_inverse % prime
            for index in range(column, size):
                work[row][index] = (work[row][index] - factor * work[column][index]) % prime
    return det % prime


def solve_mod(matrix: list[list[int]], target: list[int], prime: int) -> list[int]:
    """Solve an invertible square system over F_p."""

    size = len(matrix)
    work = [
        [entry % prime for entry in matrix[row]] + [target[row] % prime]
        for row in range(size)
    ]
    for column in range(size):
        pivot = next((row for row in range(column, size) if work[row][column]), None)
        require(pivot is not None, "attempted to invert a zero divisor")
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
        pivot_inverse = pow(work[column][column], prime - 2, prime)
        work[column] = [entry * pivot_inverse % prime for entry in work[column]]
        for row in range(size):
            if row == column or work[row][column] == 0:
                continue
            factor = work[row][column]
            work[row] = [
                (work[row][index] - factor * work[column][index]) % prime
                for index in range(size + 1)
            ]
    return [work[row][-1] for row in range(size)]


@dataclass(frozen=True)
class PrimeField:
    prime: int

    @property
    def dimension(self) -> int:
        return 1

    def zero(self) -> int:
        return 0

    def one(self) -> int:
        return 1

    def const(self, value: int) -> int:
        return value % self.prime

    def add(self, left: int, right: int) -> int:
        return (left + right) % self.prime

    def neg(self, value: int) -> int:
        return (-value) % self.prime

    def sub(self, left: int, right: int) -> int:
        return (left - right) % self.prime

    def mul(self, left: int, right: int) -> int:
        return left * right % self.prime

    def scale(self, value: int, scalar: int) -> int:
        return value * scalar % self.prime

    def pow(self, value: int, exponent: int) -> int:
        return pow(value % self.prime, exponent, self.prime)

    def inv(self, value: int) -> int:
        require(value % self.prime != 0, "attempted to invert zero in F_p")
        return pow(value % self.prime, self.prime - 2, self.prime)

    def div(self, left: int, right: int) -> int:
        return self.mul(left, self.inv(right))

    def is_zero(self, value: int) -> bool:
        return value % self.prime == 0

    def embed(self, value: int) -> int:
        return value % self.prime

    def flatten(self, value: int) -> list[int]:
        return [value % self.prime]

    def unflatten(self, values: list[int]) -> int:
        require(len(values) == 1, "prime-field unflatten length changed")
        return values[0] % self.prime

    def basis(self) -> list[int]:
        return [1]


class CubicAlgebra:
    """A[W]/(W^3+pW+q), represented in the basis 1,W,W^2."""

    def __init__(self, base, cubic_p, cubic_q):
        self.base = base
        self.cubic_p = cubic_p
        self.cubic_q = cubic_q
        self.prime = base.prime
        self.dimension = 3 * base.dimension

    def zero(self):
        return (self.base.zero(), self.base.zero(), self.base.zero())

    def one(self):
        return (self.base.one(), self.base.zero(), self.base.zero())

    def const(self, value: int):
        return (self.base.const(value), self.base.zero(), self.base.zero())

    def embed(self, value):
        return (value, self.base.zero(), self.base.zero())

    def generator(self):
        return (self.base.zero(), self.base.one(), self.base.zero())

    def add(self, left, right):
        return tuple(self.base.add(left[index], right[index]) for index in range(3))

    def neg(self, value):
        return tuple(self.base.neg(entry) for entry in value)

    def sub(self, left, right):
        return tuple(self.base.sub(left[index], right[index]) for index in range(3))

    def scale(self, value, scalar: int):
        return tuple(self.base.scale(entry, scalar) for entry in value)

    def mul(self, left, right):
        raw = [self.base.zero() for _ in range(5)]
        for left_degree in range(3):
            for right_degree in range(3):
                degree = left_degree + right_degree
                raw[degree] = self.base.add(
                    raw[degree], self.base.mul(left[left_degree], right[right_degree])
                )
        for degree in (4, 3):
            coefficient = raw[degree]
            raw[degree] = self.base.zero()
            raw[degree - 2] = self.base.sub(
                raw[degree - 2], self.base.mul(coefficient, self.cubic_p)
            )
            raw[degree - 3] = self.base.sub(
                raw[degree - 3], self.base.mul(coefficient, self.cubic_q)
            )
        return tuple(raw[:3])

    def pow(self, value, exponent: int):
        require(exponent >= 0, "negative exponent passed to ring power")
        answer = self.one()
        base_value = value
        remaining = exponent
        while remaining:
            if remaining & 1:
                answer = self.mul(answer, base_value)
            base_value = self.mul(base_value, base_value)
            remaining >>= 1
        return answer

    def is_zero(self, value) -> bool:
        return all(self.base.is_zero(entry) for entry in value)

    def flatten(self, value) -> list[int]:
        answer: list[int] = []
        for coefficient in value:
            answer.extend(self.base.flatten(coefficient))
        return answer

    def unflatten(self, values: list[int]):
        require(len(values) == self.dimension, "cubic unflatten length changed")
        block = self.base.dimension
        return tuple(
            self.base.unflatten(values[index * block : (index + 1) * block])
            for index in range(3)
        )

    def basis(self):
        answer = []
        for cubic_degree in range(3):
            for base_vector in self.base.basis():
                coefficients = [self.base.zero(), self.base.zero(), self.base.zero()]
                coefficients[cubic_degree] = base_vector
                answer.append(tuple(coefficients))
        return answer

    def multiplication_matrix(self, value) -> list[list[int]]:
        columns = [self.flatten(self.mul(value, vector)) for vector in self.basis()]
        return [
            [columns[column][row] for column in range(self.dimension)]
            for row in range(self.dimension)
        ]

    def inv(self, value):
        matrix = self.multiplication_matrix(value)
        target = self.flatten(self.one())
        solution = solve_mod(matrix, target, self.prime)
        inverse = self.unflatten(solution)
        require(self.mul(value, inverse) == self.one(), "regular inverse check failed")
        return inverse

    def div(self, left, right):
        return self.mul(left, self.inv(right))

    def norm_to_base(self, value):
        one = self.base.one()
        zero = self.base.zero()
        columns = [
            self.mul(value, (one, zero, zero)),
            self.mul(value, (zero, one, zero)),
            self.mul(value, (zero, zero, one)),
        ]
        matrix = [[columns[column][row] for column in range(3)] for row in range(3)]
        positive = self.base.add(
            self.base.mul(matrix[0][0], self.base.mul(matrix[1][1], matrix[2][2])),
            self.base.add(
                self.base.mul(matrix[0][1], self.base.mul(matrix[1][2], matrix[2][0])),
                self.base.mul(matrix[0][2], self.base.mul(matrix[1][0], matrix[2][1])),
            ),
        )
        negative = self.base.add(
            self.base.mul(matrix[0][2], self.base.mul(matrix[1][1], matrix[2][0])),
            self.base.add(
                self.base.mul(matrix[0][1], self.base.mul(matrix[1][0], matrix[2][2])),
                self.base.mul(matrix[0][0], self.base.mul(matrix[1][2], matrix[2][1])),
            ),
        )
        return self.base.sub(positive, negative)


def ring_sum(ring, *values):
    answer = ring.zero()
    for value in values:
        answer = ring.add(answer, value)
    return answer


def ring_product(ring, *values):
    answer = ring.one()
    for value in values:
        answer = ring.mul(answer, value)
    return answer


def invariants(ring, point):
    a, b, c = point
    a2 = ring.pow(a, 2)
    b2 = ring.pow(b, 2)
    b3 = ring.mul(b2, b)
    c2 = ring.pow(c, 2)
    abc = ring_product(ring, a, b, c)
    L = ring_sum(
        ring,
        ring.scale(ring_product(ring, a2, c2), 27),
        ring.scale(abc, -18),
        ring.scale(a, 16),
        ring_product(ring, b3, c),
        ring.scale(b2, -1),
    )
    T = ring_sum(ring, ring.const(4), ring.scale(ring.mul(b, c), -3))
    S = ring_sum(
        ring,
        ring.scale(ring_product(ring, a, c2), 27),
        ring.scale(ring.mul(b, c), -9),
        ring.const(8),
    )
    K = ring_sum(ring, ring.scale(ring.mul(a, c), 9), ring.neg(b))
    M = ring_sum(
        ring,
        ring.scale(ring_product(ring, a, c2), 27),
        ring.scale(ring.mul(b, c), -9),
        ring.const(26),
    )
    Y0 = ring_sum(
        ring,
        ring.scale(ring_product(ring, a, b, c2), 81),
        ring.scale(ring_product(ring, a, c), -72),
        ring.scale(ring_product(ring, b2, c), -15),
        ring.scale(b, 16),
    )
    A1 = ring_sum(
        ring,
        ring.scale(ring_product(ring, a, b, c2), 27),
        ring.scale(ring_product(ring, a, c), 54),
        ring.scale(ring_product(ring, b2, c), -9),
        ring.scale(b, 2),
    )
    A2 = ring_sum(
        ring,
        ring.scale(ring_product(ring, a, b2, c2), 27),
        ring.scale(abc, 18),
        ring.scale(a, -48),
        ring.scale(ring_product(ring, b3, c), -9),
        ring.scale(b2, 10),
    )
    Z0 = ring_sum(
        ring,
        ring.scale(ring.mul(A2, T), -9),
        ring.scale(ring.mul(M, L), -4),
    )
    return {"L": L, "T": T, "S": S, "K": K, "Y0": Y0, "A1": A1, "A2": A2, "Z0": Z0}


def fixed_map(ring, point):
    x, y, z = point
    xy = ring.mul(x, y)
    u = ring.add(ring.one(), xy)
    four_plus = ring_sum(ring, ring.const(4), ring.scale(xy, 3))
    first = ring_sum(
        ring,
        ring.mul(ring.pow(u, 3), z),
        ring_product(ring, ring.pow(y, 2), u, four_plus),
    )
    second = ring_sum(
        ring,
        y,
        ring.scale(ring_product(ring, x, ring.pow(u, 2), z), 3),
        ring.scale(ring_product(ring, x, ring.pow(y, 2), four_plus), 3),
    )
    third = ring_sum(
        ring,
        ring.scale(x, 2),
        ring.scale(ring_product(ring, ring.pow(x, 2), y), -3),
        ring.neg(ring_product(ring, ring.pow(x, 3), z)),
    )
    return first, second, third


def absolute_norm(ring, value) -> int:
    if isinstance(ring, PrimeField):
        return value % ring.prime
    return absolute_norm(ring.base, ring.norm_to_base(value))


def flat_norm(ring, value) -> int:
    if isinstance(ring, PrimeField):
        return value % ring.prime
    return determinant_mod(ring.multiplication_matrix(value), ring.prime)


def coordinate_digest(ring, point) -> str:
    payload = ";".join(
        ",".join(str(entry) for entry in ring.flatten(coordinate)) for coordinate in point
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def inverse_step(ring, target, label: str, gate_ledger: list[tuple]):
    inv = invariants(ring, target)
    inverse_L = ring.inv(inv["L"])
    inverse_S = ring.inv(inv["S"])
    cubic_p = ring.mul(inv["T"], inverse_L)
    cubic_q = ring.scale(ring.mul(target[2], inverse_L), -2)
    discriminant = ring_sum(
        ring,
        ring.scale(ring.pow(cubic_p, 3), -4),
        ring.scale(ring.pow(cubic_q, 2), -27),
    )
    ring.inv(discriminant)

    extension = CubicAlgebra(ring, cubic_p, cubic_q)
    w = extension.generator()
    inverse_2S = ring.scale(inverse_S, pow(2, ring.prime - 2, ring.prime))
    inverse_8S = ring.scale(inverse_S, pow(8, ring.prime - 2, ring.prime))
    y = (
        ring.mul(inv["Y0"], inverse_2S),
        ring.mul(ring.scale(inv["L"], 6), inverse_2S),
        ring.mul(ring.scale(ring.mul(inv["K"], inv["L"]), -3), inverse_2S),
    )
    z = (
        ring.mul(inv["Z0"], inverse_8S),
        ring.mul(ring.scale(ring.mul(inv["L"], inv["A1"]), 6), inverse_8S),
        ring.mul(ring.scale(ring.mul(inv["L"], inv["A2"]), -9), inverse_8S),
    )
    source = (w, y, z)
    embedded_target = tuple(extension.embed(coordinate) for coordinate in target)
    require(fixed_map(extension, source) == embedded_target, f"{label}: inverse graph changed")

    gate_ledger.append(
        (
            label,
            ring.dimension,
            absolute_norm(ring, inv["L"]),
            absolute_norm(ring, inv["S"]),
            absolute_norm(ring, discriminant),
            coordinate_digest(extension, source),
        )
    )
    return extension, source, inv["L"]


ROOT = Path(__file__).resolve().parents[1]
H_PATH = ROOT / "05-knowledge/results/keller_L2_polynomial_opus_20260728.pkl"
H_RAW_SHA256 = "5a9459b3149e500c1b00b67bd804aa7e607de06bf4610c7cdf5fa26d41d74ce9"
H_LEDGER_SHA256 = "b146c11f33e895c08f72303d282e2668d955e0a58d9268a1b445d4d5202016c2"

raw_h = H_PATH.read_bytes()
require(hashlib.sha256(raw_h).hexdigest() == H_RAW_SHA256, "frozen H pickle changed")
H = pickle.loads(raw_h)
H_poly = sp.Poly(H)
H_terms = [(tuple(map(int, monomial)), int(coefficient)) for monomial, coefficient in H_poly.terms()]
H_ledger = "\n".join(f"{monomial}:{coefficient}" for monomial, coefficient in H_poly.terms())
require(hashlib.sha256(H_ledger.encode("ascii")).hexdigest() == H_LEDGER_SHA256, "H ledger changed")
require(len(H_terms) == 361 and H_poly.degree_list() == (14, 21, 12), "H support changed")


def evaluate_H(ring, point):
    maxima = H_poly.degree_list()
    power_tables = []
    for coordinate, maximum in zip(point, maxima):
        powers = [ring.one()]
        for _ in range(maximum):
            powers.append(ring.mul(powers[-1], coordinate))
        power_tables.append(powers)
    answer = ring.zero()
    for (i, j, k), coefficient in H_terms:
        term = ring_product(ring, power_tables[0][i], power_tables[1][j], power_tables[2][k])
        answer = ring.add(answer, ring.scale(term, coefficient))
    return answer


TRANSITIONS = {
    1: (2**6, 1),
    2: (2**35, 7),
    3: (1, 43),
    4: (1, 271),
}


def evaluate_from_L(
    ring,
    point,
    level: int,
    ledger: list[tuple],
    constants=TRANSITIONS,
    bottom_capture: list[int] | None = None,
):
    """Evaluate P_level, with P_0=L, through the complete norm recursion."""

    if level == 0:
        value = invariants(ring, point)["L"]
        if bottom_capture is not None:
            bottom_capture.append(absolute_norm(ring, value))
        return value
    extension, source, l_value = inverse_step(ring, point, f"L-route-{level}", ledger)
    lower_value = evaluate_from_L(
        extension, source, level - 1, ledger, constants, bottom_capture
    )
    coefficient, exponent = constants[level]
    return ring_product(
        ring,
        ring.const(coefficient),
        ring.pow(l_value, exponent),
        extension.norm_to_base(lower_value),
    )


def evaluate_from_H(ring, point, level: int, ledger: list[tuple], base_capture: list[tuple]):
    """Evaluate P_level, with explicit P_1=H, through the upper recursion."""

    require(level >= 1, "H route starts at level one")
    if level == 1:
        value = evaluate_H(ring, point)
        base_capture.append((ring, value))
        return value
    extension, source, l_value = inverse_step(ring, point, f"H-route-{level}", ledger)
    lower_value = evaluate_from_H(extension, source, level - 1, ledger, base_capture)
    coefficient, exponent = TRANSITIONS[level]
    return ring_product(
        ring,
        ring.const(coefficient),
        ring.pow(l_value, exponent),
        extension.norm_to_base(lower_value),
    )


FINITE_Q = (Fraction(2), Fraction(5, 6), Fraction(-7, 8))
OLD_BOUNDARY_TARGET = (Fraction(2, 27), Fraction(1), Fraction(1))


def exact_fixed_map(point):
    x, y, z = point
    u = 1 + x * y
    return (
        u**3 * z + y**2 * u * (4 + 3 * x * y),
        y + 3 * x * u**2 * z + 3 * x * y**2 * (4 + 3 * x * y),
        2 * x - 3 * x**2 * y - x**3 * z,
    )


def exact_L(point):
    a, b, c = point
    return 27 * a**2 * c**2 - 18 * a * b * c + 16 * a + b**3 * c - b**2


require(exact_fixed_map(FINITE_Q) == OLD_BOUNDARY_TARGET, "exact finite-sheet point changed")
require(exact_L(OLD_BOUNDARY_TARGET) == 0, "exact target left the old L boundary")
require(exact_L(FINITE_Q) == Fraction(241465, 1728), "exact L(q) ledger changed")


def fraction_mod(value: Fraction, prime: int) -> int:
    require(value.denominator % prime != 0, f"p={prime}: rational denominator vanished")
    return value.numerator * pow(value.denominator, prime - 2, prime) % prime


def evaluate_prime(prime: int):
    field = PrimeField(prime)
    q = tuple(fraction_mod(value, prime) for value in FINITE_Q)
    boundary_target = tuple(fraction_mod(value, prime) for value in OLD_BOUNDARY_TARGET)
    require(fixed_map(field, q) == boundary_target, f"p={prime}: canonical finite sheet changed")
    require(field.is_zero(invariants(field, boundary_target)["L"]), f"p={prime}: target left L=0")
    require(not field.is_zero(invariants(field, q)["L"]), f"p={prime}: finite point hit L=0")

    l_ledger: list[tuple] = []
    h_ledger: list[tuple] = []
    h_capture: list[tuple] = []
    l_bottom: list[int] = []
    value_from_l = evaluate_from_L(field, q, 4, l_ledger, bottom_capture=l_bottom)
    value_from_h = evaluate_from_H(field, q, 4, h_ledger, h_capture)
    require(value_from_l == value_from_h, f"p={prime}: L and explicit-H routes disagree")
    require(value_from_l != 0, f"p={prime}: R5 finite-sheet reduction vanished")
    require(len(l_ledger) == 4 and len(h_ledger) == 3, f"p={prime}: tower depth changed")
    require([row[1] for row in l_ledger] == [1, 3, 9, 27], f"p={prime}: L dimensions changed")
    require([row[1] for row in h_ledger] == [1, 3, 9], f"p={prime}: H dimensions changed")
    require(len(l_bottom) == 1, f"p={prime}: bottom L norm capture changed")

    # Multiplicativity of the cubic norm unrolls the four definitions to
    # R5=2^477 L^271 N(L)^43 N^2(L)^7 N^3(L) N^4(L).
    l_norm_orbit = tuple(row[2] for row in l_ledger) + (l_bottom[0],)
    unrolled_value = pow(2, 477, prime)
    for orbit_value, exponent in zip(l_norm_orbit, (271, 43, 7, 1, 1)):
        unrolled_value = unrolled_value * pow(orbit_value, exponent, prime) % prime
    require(unrolled_value == value_from_l, f"p={prime}: unrolled norm-orbit identity failed")

    h_ring, h_value = h_capture[0]
    h_norm_recursive = absolute_norm(h_ring, h_value)
    h_norm_flat = flat_norm(h_ring, h_value)
    require(h_norm_recursive == h_norm_flat, f"p={prime}: transitive/flat H norm mismatch")
    require(h_norm_flat != 0, f"p={prime}: explicit H became a zero divisor")

    hostile_constants = dict(TRANSITIONS)
    hostile_constants[1] = (1, 1)
    hostile_ledger: list[tuple] = []
    wrong_h_normalization = evaluate_from_L(field, q, 4, hostile_ledger, hostile_constants)
    expected_wrong = value_from_l * pow(pow(64, 27, prime), prime - 2, prime) % prime
    require(wrong_h_normalization == expected_wrong, f"p={prime}: scalar norm exponent changed")
    require(wrong_h_normalization != value_from_l, f"p={prime}: normalization hostile did not fire")

    return {
        "prime": prime,
        "value": value_from_l,
        "wrong_h_normalization": wrong_h_normalization,
        "h_flat_norm": h_norm_flat,
        "l_norm_orbit": l_norm_orbit,
        "l_ledger": l_ledger,
        "h_ledger": h_ledger,
    }


records = [evaluate_prime(prime) for prime in (101, 103, 107)]

R5_FACE_PAIR = (1699, 615)
R6_FACE_PAIR = (
    7 * R5_FACE_PAIR[0] - 2 * R5_FACE_PAIR[1],
    3 * R5_FACE_PAIR[0] - 2 * R5_FACE_PAIR[1],
)
require(R6_FACE_PAIR == (10663, 3867), "next exposed face pair changed")

semantic_lines = []
for record in records:
    semantic_lines.append(
        f"p={record['prime']};R5={record['value']};wrong64={record['wrong_h_normalization']};"
        f"Hflat={record['h_flat_norm']};Lorbit={record['l_norm_orbit']}"
    )
    for row in record["l_ledger"]:
        semantic_lines.append("L:" + ":".join(map(str, row)))
    for row in record["h_ledger"]:
        semantic_lines.append("H:" + ":".join(map(str, row)))
semantic_sha256 = hashlib.sha256("\n".join(semantic_lines).encode("ascii")).hexdigest()

print("== fixed Keller R5 finite-sheet recursive norm gate ==")
print(f"H raw sha256={H_RAW_SHA256}; terms={len(H_terms)}; degrees={H_poly.degree_list()}")
print("point q=(2,5/6,-7/8); F(q)=(2/27,1,1); L(F(q))=0; L(q)=241465/1728")
print("primary route: four cubic algebras, descending to L")
print("independent route: three cubic algebras, evaluating frozen H directly")
for record in records:
    print(
        f"p={record['prime']}: R5(q)={record['value']}; "
        f"explicit-H flat norm={record['h_flat_norm']}; "
        f"wrong-H-normalization control={record['wrong_h_normalization']}"
    )
    print(f"  (L,N(L),N^2(L),N^3(L),N^4(L))={record['l_norm_orbit']}")
    for route_name in ("l_ledger", "h_ledger"):
        compact = [
            (row[1], row[2], row[3], row[4], row[5][:12]) for row in record[route_name]
        ]
        print(f"  {route_name}: (base_dim,Norm(L),Norm(S),Norm(disc),source_sha12)={compact}")
print(f"semantic sha256={semantic_sha256}")
print("finite-sheet verdict: R5(2,5/6,-7/8) != 0 over Q")
print("consequence: v_L(N(R5))=-1699 and L^1699*N(R5) is polynomial and coprime to L")
print("next exposed face: in_max-lambda(L^1699*N(R5))=C*x^10663*(3*x*z-2*y)^3867, C!=0")
print("scope: fixed map/old-L gate only; renewal, fifth image, degree 243, and general JC remain open")
print("all exact checks passed")
