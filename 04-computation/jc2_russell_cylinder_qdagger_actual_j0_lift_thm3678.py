"""Exact minimal-packet J0=1 target-ring lift for Q_dagger (THM-3678)."""

import hashlib
import math

import sympy as sp
from flint import fmpq, fmpq_mat


x, q = sp.symbols("x q")
Q = (
    22868 * x**8
    - 89583 * x**6
    + 2916 * x**5
    + 123684 * x**4
    - 5832 * x**3
    - 63530 * x**2
    + 2916 * x
    - 2187
) / sp.Integer(2916)
D = 1 + x**2 * q
b = sp.expand((D - 1) * (D + 2) ** 2)
c = sp.expand(x * D * (D + 2))
e = sp.expand(q * (D + 3))

B = sp.Poly(b.subs(q, Q), x, domain=sp.QQ)
C = sp.Poly(c.subs(q, Q), x, domain=sp.QQ)
E = sp.Poly(e.subs(q, Q), x, domain=sp.QQ)
dB = sp.Poly(sp.diff(b, q).subs(q, Q), x, domain=sp.QQ)
dC = sp.Poly(sp.diff(c, q).subs(q, Q), x, domain=sp.QQ)
dE = sp.Poly(sp.diff(e, q).subs(q, Q), x, domain=sp.QQ)
GENERATORS = (B, C, E)
DELTA_GENERATORS = (dB, dC, dE)
WEIGHTS = tuple(generator.degree() for generator in GENERATORS)

EXPECTED = {
    "F_target": "587dc21dbba0b2bd48356bac7d8ab9bfbd7e44923aee8e35fe2237b00068cf72",
    "G_target": "2aed39df3314c9e74040c70afe9fc63753bcec0f00e0bd5276a2293fe2cdca6c",
    "F": "1180c242b841d42de48806f962deb9fe2caa182b6435929ec4612bc4859c36cb",
    "G": "942601ed5f8b7ceabdd15221affb3cb5b0ef98acc331859db0fd2527eeb5bcc3",
    "deltaF": "fca72977cd805fb4cb8253115d74c61e34b45b9c924711bec2baae916ce78aba",
    "deltaG": "0e3b6b7ae43928d7773d21839432c2faed7861fe5772270b20977a58fa683d10",
}

gates = 0


def require(condition, label):
    global gates
    if not condition:
        raise RuntimeError(label)
    gates += 1


def as_fmpq(value):
    value = sp.Rational(value)
    return fmpq(int(value.p), int(value.q))


def as_rational(value):
    return sp.Rational(int(value.numerator), int(value.denominator))


def packet(cutoff):
    """Raw b^i c^j e^k target packet and its restriction/normal derivative."""
    powers = []
    for generator, weight in zip(GENERATORS, WEIGHTS):
        row = [sp.Poly(1, x, domain=sp.QQ)]
        for _ in range(cutoff // weight):
            row.append(row[-1] * generator)
        powers.append(row)

    answer = []
    for i in range(cutoff // WEIGHTS[0] + 1):
        for j in range(cutoff // WEIGHTS[1] + 1):
            for k in range(cutoff // WEIGHTS[2] + 1):
                metadata = (i, j, k)
                if sum(exponent * weight for exponent, weight in zip(metadata, WEIGHTS)) > cutoff:
                    continue
                restriction = powers[0][i] * powers[1][j] * powers[2][k]
                delta = sp.Poly(0, x, domain=sp.QQ)
                for coordinate, exponent in enumerate(metadata):
                    if not exponent:
                        continue
                    term = sp.Poly(exponent, x, domain=sp.QQ)
                    for index, coordinate_exponent in enumerate(metadata):
                        term *= powers[index][coordinate_exponent - int(index == coordinate)]
                    delta += term * DELTA_GENERATORS[coordinate]
                answer.append((metadata, restriction, delta))
    return answer


def matrix_at(cutoff):
    monomials = packet(cutoff)
    restrictions = [item[1] for item in monomials]
    columns = [C.diff() * polynomial for polynomial in restrictions]
    columns += [-E.diff() * polynomial for polynomial in restrictions]
    row_count = max(column.degree() for column in columns) + 1
    operator_data = []
    augmented_data = []
    for row in range(row_count):
        row_data = [as_fmpq(column.nth(row)) for column in columns]
        operator_data.extend(row_data)
        augmented_data.extend(row_data)
        augmented_data.append(fmpq(int(row == 0)))
    operator = fmpq_mat(row_count, len(columns), operator_data)
    augmented = fmpq_mat(row_count, len(columns) + 1, augmented_data)
    return monomials, columns, operator, augmented


def polynomial_hash(polynomial):
    polynomial = sp.Poly(polynomial, x, domain=sp.QQ)
    payload = ";".join(
        f"{degree}:{polynomial.nth(degree)}" for degree in range(polynomial.degree() + 1)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def target_hash(values, monomials):
    payload = ";".join(
        f"{i},{j},{k}:{value}"
        for value, ((i, j, k), _restriction, _delta) in zip(values, monomials)
        if value
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


# Compiler, retained-fold, and grading controls.
require(C * C * E == B * (B + 4), "target-cylinder relation failed")
require(tuple(sp.Rational(Q.subs(x, point)) for point in (-1, 0, 1)) == (-3, -sp.Rational(3, 4), -3), "retained values")
require(tuple(sp.Rational(sp.diff(Q, x).subs(x, point)) for point in (-1, 0, 1)) == (-sp.Rational(9, 2), 1, sp.Rational(9, 2)), "retained slopes")
require(WEIGHTS == (30, 21, 18), "restriction weights")
require(math.gcd(*WEIGHTS) == 3, "packet grading gcd")
require(C.diff().gcd(E.diff()).degree() == 0, "Bezout necessary condition")

# Exact sharp lower packet: augmented rank exceeds operator rank.
lower_monomials, _lower_columns, lower_operator, lower_augmented = matrix_at(192)
lower_operator_rank = lower_operator.rank()
_lower_reduced, lower_augmented_rank = lower_augmented.rref()
require(len(lower_monomials) == 173, "cutoff-192 packet size")
require((lower_operator.nrows(), lower_operator.ncols()) == (213, 346), "cutoff-192 shape")
require((lower_operator_rank, lower_augmented_rank) == (208, 209), "cutoff-192 exact rank obstruction")

# First larger packet and deterministic exact RREF certificate.
monomials, columns, operator, augmented = matrix_at(195)
operator_rank = operator.rank()
reduced, augmented_rank = augmented.rref()
require(len(monomials) == 179, "cutoff-195 packet size")
require((operator.nrows(), operator.ncols()) == (216, 358), "cutoff-195 shape")
require((operator_rank, augmented_rank) == (214, 214), "cutoff-195 exact feasibility")

solution = [fmpq(0)] * len(columns)
for row in range(augmented_rank):
    pivot = next((column for column in range(len(columns) + 1) if reduced[row, column]), None)
    require(pivot is not None and pivot < len(columns), f"invalid RREF pivot row {row}")
    solution[pivot] = reduced[row, len(columns)]

G_values = solution[: len(monomials)]
F_values = solution[len(monomials) :]
F = sp.Poly(0, x, domain=sp.QQ)
G = sp.Poly(0, x, domain=sp.QQ)
delta_F = sp.Poly(0, x, domain=sp.QQ)
delta_G = sp.Poly(0, x, domain=sp.QQ)
for value, (_metadata, restriction, delta) in zip(F_values, monomials):
    if value:
        coefficient = as_rational(value)
        F += coefficient * restriction
        delta_F += coefficient * delta
for value, (_metadata, restriction, delta) in zip(G_values, monomials):
    if value:
        coefficient = as_rational(value)
        G += coefficient * restriction
        delta_G += coefficient * delta

require(C.diff() * G - F * E.diff() == sp.Poly(1, x, domain=sp.QQ), "exact J0 residual")
require(sum(bool(value) for value in F_values) == 107, "F target support")
require(sum(bool(value) for value in G_values) == 104, "G target support")
require((F.degree(), G.degree(), delta_F.degree(), delta_G.degree()) == (195, 192, 187, 184), "restriction degrees")
require(target_hash(F_values, monomials) == EXPECTED["F_target"], "F target hash")
require(target_hash(G_values, monomials) == EXPECTED["G_target"], "G target hash")
require(polynomial_hash(F) == EXPECTED["F"], "F restriction hash")
require(polynomial_hash(G) == EXPECTED["G"], "G restriction hash")
require(polynomial_hash(delta_F) == EXPECTED["deltaF"], "delta F hash")
require(polynomial_hash(delta_G) == EXPECTED["deltaG"], "delta G hash")
require(tuple(F.eval(point) for point in (-1, 0, 1)) == (0, 0, 0), "retained F")
require(tuple(G.eval(point) for point in (-1, 0, 1)) == (sp.Rational(1, 3),) * 3, "retained G")

print(f"weights={WEIGHTS};gcd={math.gcd(*WEIGHTS)};derivative_gcd_degree={C.diff().gcd(E.diff()).degree()}")
print(
    f"cutoff=192;monomials={len(lower_monomials)};shape={lower_operator.nrows()}x{lower_operator.ncols()};"
    f"operator_rank={lower_operator_rank};augmented_rank={lower_augmented_rank};feasible=False"
)
print(
    f"cutoff=195;monomials={len(monomials)};shape={operator.nrows()}x{operator.ncols()};"
    f"operator_rank={operator_rank};augmented_rank={augmented_rank};feasible=True"
)
print(f"F_support=107;G_support=104;F_target_hash={EXPECTED['F_target']};G_target_hash={EXPECTED['G_target']}")
print(f"F_degree=195;F_hash={EXPECTED['F']};G_degree=192;G_hash={EXPECTED['G']}")
print(f"deltaF_degree=187;deltaF_hash={EXPECTED['deltaF']}")
print(f"deltaG_degree=184;deltaG_hash={EXPECTED['deltaG']}")
print("retained_F=(0, 0, 0);retained_G=(1/3, 1/3, 1/3)")
print(f"RESULT=PASS;gates={gates}")
