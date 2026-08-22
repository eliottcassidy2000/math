"""Exact quartic-field J0 lift at the THM-3683 exceptional folds.

The finite-field split fibre (p,r)=(137,44) chooses a square only.  Every
entry of that square is rebuilt in K=Q[a]/(F(a)), expanded in the power basis
(1,a,a^2,a^3), solved over Q, and substituted into the full K[x] identity.
"""

import hashlib
import importlib.util
from pathlib import Path

import sympy as sp
from flint import fmpq, fmpq_mat
from sympy.polys.rings import ring


PROBE_PATH = Path("04-computation/jc2_russell_cylinder_exceptional_quartic_modular_lift_thm3687.py")
spec = importlib.util.spec_from_file_location("exceptional_modular_probe", PROBE_PATH)
probe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(probe)


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


alpha_symbol = sp.symbols("alpha")
field = sp.QQ.alg_field_from_poly(sp.Poly(probe.F_QUARTIC.subs(probe.r, alpha_symbol), alpha_symbol), "alpha")
alpha = field.from_sympy(field.ext.as_expr())
polynomial_ring, X = ring("X", field)


def field_coefficient(expression):
    return field.from_sympy(sp.cancel(expression).subs(probe.r, field.ext.as_expr()))


Q = polynomial_ring.zero
for (degree,), coefficient in probe.Q_R.terms():
    Q += field_coefficient(coefficient.as_expr()) * X**degree

one_poly = polynomial_ring.one
D = one_poly + X**2 * Q
B = (D - one_poly) * (D + 2) ** 2
C = X * D * (D + 2)
E = Q * (D + 3)
dB = 3 * X**2 * D * (D + 2)
dC = 2 * X**3 * (D + 1)
dE = 2 * (D + 1)

require(C * C * E == B * (B + 4), "quartic-field compiler relation")
require((B.degree(), C.degree(), E.degree()) == (30, 21, 18), "quartic-field restriction degrees")
require(
    tuple(Q.evaluate(X, field.convert(point)) for point in (-1, 0, 1))
    == (field.convert(-3), field.convert(sp.Rational(-3, 4)), field.convert(-3)),
    "quartic-field retained values",
)
require(
    tuple(Q.diff(X).evaluate(X, field.convert(point)) for point in (-1, 0, 1))
    == (field.convert(sp.Rational(-9, 2)), field.one, field.convert(sp.Rational(9, 2))),
    "quartic-field retained slopes",
)
require(
    5 * Q.diff(X).diff(X).evaluate(X, field.convert(-1))
    + 13 * Q.diff(X).diff(X).evaluate(X, field.convert(1))
    + field.convert(243)
    == field.zero,
    "quartic-field zero second debt",
)
p_alpha = field_coefficient(probe.P_OF_R)
require(
    729 * p_alpha - 42120 * alpha**2 + 15192 * alpha + field.convert(5717)
    == field.zero,
    "quartic-field zero fourth debt",
)
require(field_coefficient(probe.F_QUARTIC) == field.zero, "quartic-field zero sixth debt")


def packet(cutoff):
    generators = (B, C, E)
    deltas = (dB, dC, dE)
    weights = (30, 21, 18)
    powers = []
    for generator, weight in zip(generators, weights):
        row = [one_poly]
        for _ in range(cutoff // weight):
            row.append(row[-1] * generator)
        powers.append(row)
    answer = []
    for i in range(cutoff // weights[0] + 1):
        for j in range(cutoff // weights[1] + 1):
            for k in range(cutoff // weights[2] + 1):
                metadata = (i, j, k)
                if sum(a * b for a, b in zip(metadata, weights)) > cutoff:
                    continue
                restriction = powers[0][i] * powers[1][j] * powers[2][k]
                delta = polynomial_ring.zero
                for coordinate, exponent in enumerate(metadata):
                    if not exponent:
                        continue
                    term = polynomial_ring(exponent)
                    for index, coordinate_exponent in enumerate(metadata):
                        term *= powers[index][coordinate_exponent - int(index == coordinate)]
                    delta += term * deltas[coordinate]
                answer.append((metadata, restriction, delta))
    return answer


monomials = packet(198)
restrictions = [item[1] for item in monomials]
columns = [C.diff(X) * polynomial for polynomial in restrictions]
columns += [-E.diff(X) * polynomial for polynomial in restrictions]
require(len(monomials) == 187 and len(columns) == 374, "quartic-field packet shape")
require(max(column.degree() for column in columns) == 218, "quartic-field row degree")

# Select the same square in one split good fibre, but do not reuse any lifted
# coefficients from that fibre.
selector = probe.solve_j0(137, 44, 198)
pivot_columns = selector["pivot_columns"]
pivot_rows = selector["pivot_rows"]
require(len(pivot_columns) == len(pivot_rows) == 218, "selector square size")

# Two completely split good-prime fibres provide hostile finite-field controls.
# They are not used to infer characteristic-zero feasibility or minimality.
for split_prime in (137, 163):
    split_roots = probe.roots_mod_prime(split_prime)
    require(len(split_roots) == 4, f"split roots {split_prime}")
    for split_root in split_roots:
        lower = probe.solve_j0(split_prime, split_root, 195)
        upper = probe.solve_j0(split_prime, split_root, 198)
        require(
            (lower["rank"], lower["augmented"]) == (214, 215),
            f"modular lower obstruction {split_prime}:{split_root}",
        )
        require(
            (upper["rank"], upper["augmented"]) == (218, 218),
            f"modular upper feasibility {split_prime}:{split_root}",
        )


def to_fmpq(value):
    value = sp.Rational(value)
    return fmpq(int(value.p), int(value.q))


def from_fmpq(value):
    return sp.Rational(int(value.numerator), int(value.denominator))


def field_coordinates(value):
    """Return low-to-high rational coordinates in 1,alpha,alpha^2,alpha^3."""
    entries = list(value.to_list())
    entries = [sp.Rational(0)] * (4 - len(entries)) + [sp.Rational(entry) for entry in entries]
    return list(reversed(entries))


alpha_powers = [field.one, alpha, alpha**2, alpha**3]


def multiplication_block(value):
    """Four rows of the power-basis matrix for multiplication by value."""
    columns_as_coordinates = [field_coordinates(value * power) for power in alpha_powers]
    return [
        [columns_as_coordinates[column][row] for column in range(4)]
        for row in range(4)
    ]


# Expand the 218-square over K to an 872-square over Q.
dimension = 4 * len(pivot_rows)
matrix_data = []
rhs_data = []
zero_block = [[sp.Rational(0)] * 4 for _ in range(4)]
blocks = []
for row in pivot_rows:
    block_row = []
    for column in pivot_columns:
        coefficient = columns[column].get((row,), field.zero)
        block_row.append(multiplication_block(coefficient) if coefficient else zero_block)
    blocks.append(block_row)

for block_row_index, source_row in enumerate(pivot_rows):
    for coordinate_row in range(4):
        for block in blocks[block_row_index]:
            matrix_data.extend(to_fmpq(value) for value in block[coordinate_row])
        rhs_data.append(fmpq(int(source_row == 0 and coordinate_row == 0)))

matrix = fmpq_mat(dimension, dimension, matrix_data)
rhs = fmpq_mat(dimension, 1, rhs_data)
solution_rational = matrix.solve(rhs)


def field_from_solution(block_index):
    answer = field.zero
    for coordinate, power in enumerate(alpha_powers):
        answer += field.convert(from_fmpq(solution_rational[4 * block_index + coordinate, 0])) * power
    return answer


solution = [field.zero] * len(columns)
for block_index, column in enumerate(pivot_columns):
    solution[column] = field_from_solution(block_index)

size = len(monomials)
G_values = solution[:size]
F_values = solution[size:]
F1 = polynomial_ring.zero
G1 = polynomial_ring.zero
delta_F1 = polynomial_ring.zero
delta_G1 = polynomial_ring.zero
for value, (_metadata, restriction, delta) in zip(F_values, monomials):
    F1 += value * restriction
    delta_F1 += value * delta
for value, (_metadata, restriction, delta) in zip(G_values, monomials):
    G1 += value * restriction
    delta_G1 += value * delta

residual = C.diff(X) * G1 - F1 * E.diff(X) - one_poly
require(not residual, "full quartic-field J0 residual")


def serialize_field(value):
    return ",".join(str(item) for item in field_coordinates(value))


def target_hash(values):
    payload = ";".join(
        f"{i},{j},{k}:{serialize_field(value)}"
        for value, ((i, j, k), _restriction, _delta) in zip(values, monomials)
        if value
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def polynomial_hash(value):
    payload = ";".join(
        f"{degree}:{serialize_field(value.get((degree,), field.zero))}"
        for degree in range(value.degree() + 1)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


print("field=Q[alpha]/(72783360*alpha^4-77822208*alpha^3-28419741*alpha^2+7849770*alpha-1276420)")
print(f"cutoff=198;monomials={len(monomials)};K_square=218;Q_square={dimension};full_rows=219")
print("modular_controls=p137,p163;four_roots_each;cutoff195=214/215;cutoff198=218/218;scope=FINITE_FIELD_ONLY")
print(
    f"F1_support={sum(bool(value) for value in F_values)};G1_support={sum(bool(value) for value in G_values)};"
    f"F1_degree={F1.degree()};G1_degree={G1.degree()}"
)
print(f"F1_target_hash={target_hash(F_values)};G1_target_hash={target_hash(G_values)}")
print(f"F1_hash={polynomial_hash(F1)};G1_hash={polynomial_hash(G1)}")
print(f"deltaF1_degree={delta_F1.degree()};deltaG1_degree={delta_G1.degree()}")
print("J0=1;full_quartic_field_residual=0")
print("RESULT=PASS")
