#!/usr/bin/env python3
"""Finite-exact one-variable gate for the fourth Keller norm divisor.

Let J be THM-3495's primitive 66,146-term polynomial and N the cubic
function-field norm of the fixed sporadic Keller map.  This script computes

    G(A,1,1) = L(A,1,1)^43 * N(J)(A,1,1)  (mod 1009)

without expanding N(J) globally.  A sharp valuation-at-infinity bound gives
degree at most 542, so 543 regular target values determine the specialization
without finite-field function aliasing.  Extra values are held out.

The load-bearing evaluation is duplicated in two ways:

* J is aggregated using two different exponent-pair groupings; and
* the cubic norm is computed both by a closed formula and by a direct 3x3
  multiplication determinant.

The resulting slice is tested for squarefreeness and coprimality with the
specialized old factors L, H, and J.  This is a finite-exact specialization
gate.  The geometric image/multiplicity consequence belongs in the paired
proof document; this script itself does not infer a Jacobian-conjecture result.
"""

from __future__ import annotations

import contextlib
import gc
import hashlib
import io
import runpy
from pathlib import Path

import numpy as np
from flint import fmpz_mod_poly_ctx


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
GLOBAL_PROBE = ROOT / "04-computation/keller_level_three_global_norm_probe_20260816.py"

PRIME = 1009
B_VALUE = 1
C_VALUE = 1
DEGREE_BOUND = 542
INTERPOLATION_COUNT = DEGREE_BOUND + 1
HOLDOUT_COUNT = 12
EXPECTED_J_LEDGER_SHA256 = "9aca78e67d33351b2f2fb4dbe8ab5bdff06373fdbd8ef9ec73d29b15bffedefe"

EXPECTED_G_SLICE_SHA256 = "47fba77866ee50d00fcae28b834e8a0b0c18a4cf52c2cd1b9c05155410c91d00"
EXPECTED_FACTOR_DEGREES = [(1, 1), (2, 1), (2, 1), (4, 1), (12, 1), (21, 1), (500, 1)]


print("== F_1009 slice gate for G=L^43*Norm(J) ==", flush=True)
print("stage: reconstructing THM-3495 J from its frozen exact route", flush=True)
captured = io.StringIO()
with contextlib.redirect_stdout(captured):
    namespace = runpy.run_path(str(GLOBAL_PROBE))
require(
    captured.getvalue().rstrip().endswith("all exact checks passed"),
    "the imported global J reconstruction did not finish cleanly",
)
J = namespace["J"]
require(
    namespace["coefficient_hash"](J) == EXPECTED_J_LEDGER_SHA256,
    "THM-3495 J coefficient ledger changed",
)
J_dictionary = {
    (int(i), int(j), int(k)): int(coefficient)
    for (i, j, k), coefficient in J.to_dict().items()
}
H_terms = [
    ((int(i), int(j), int(k)), int(coefficient))
    for (i, j, k), coefficient in namespace["H_terms"]
]


def evaluate_integer_terms(terms, point: tuple[int, int, int]) -> int:
    maxima_local = tuple(max(monomial[axis] for monomial, _coefficient in terms) for axis in range(3))
    powers = []
    for value, maximum in zip(point, maxima_local):
        row = [1]
        for _ in range(maximum):
            row.append(row[-1] * value)
        powers.append(row)
    return sum(
        coefficient * powers[0][i] * powers[1][j] * powers[2][k]
        for (i, j, k), coefficient in terms
    )


def integer_fmap(point: tuple[int, int, int]) -> tuple[int, int, int]:
    x, y, z = point
    unit = 1 + x * y
    four_plus = 4 + 3 * x * y
    return (
        unit**3 * z + y**2 * unit * four_plus,
        y + 3 * x * unit**2 * z + 3 * x * y**2 * four_plus,
        2 * x - 3 * x**2 * y - x**3 * z,
    )


def integer_l(point: tuple[int, int, int]) -> int:
    a_value, b_value, c_value = point
    return (
        27 * a_value**2 * c_value**2
        - 18 * a_value * b_value * c_value
        + 16 * a_value
        + b_value**3 * c_value
        - b_value**2
    )


# A direct source-to-image witness makes the geometric consequence explicit.
# It is independent of the norm-support implication used in the proof.
q_witness = (3, -1, 0)
p_witness = integer_fmap(q_witness)
r_witness = integer_fmap(p_witness)
require(p_witness == (10, -46, 33), "first image witness changed")
require(
    r_witness == (-1854753363, 121225664, -19180),
    "second image witness changed",
)
require(evaluate_integer_terms(H_terms, q_witness) == 0, "H witness no longer vanishes")
require(
    evaluate_integer_terms(list(J_dictionary.items()), p_witness) == 0,
    "J image witness no longer vanishes",
)
require(
    (integer_l(p_witness), integer_l(r_witness)) == (-504, -69753247104),
    "image witness left the two required finite loci",
)
del J, namespace
gc.collect()
print("stage: J ledger reconstructed and released", flush=True)
print(
    "image witness: H(3,-1,0)=0; J(10,-46,33)=0; "
    "L(10,-46,33)/L(F(10,-46,33))=-504/-69753247104",
    flush=True,
)


def mod_inverse(value: int) -> int:
    value %= PRIME
    require(value != 0, "attempted inversion of zero")
    return pow(value, PRIME - 2, PRIME)


def scalar(value: int) -> np.ndarray:
    return np.array((value % PRIME, 0, 0), dtype=np.int64)


def algebra_mul(left: np.ndarray, right: np.ndarray, p_cubic: int, q_cubic: int) -> np.ndarray:
    """Multiply in F_p[w]/(w^3-p_cubic*w-q_cubic), vectorized by rows."""

    d0 = left[..., 0] * right[..., 0]
    d1 = left[..., 0] * right[..., 1] + left[..., 1] * right[..., 0]
    d2 = (
        left[..., 0] * right[..., 2]
        + left[..., 1] * right[..., 1]
        + left[..., 2] * right[..., 0]
    )
    d3 = left[..., 1] * right[..., 2] + left[..., 2] * right[..., 1]
    d4 = left[..., 2] * right[..., 2]
    return np.stack(
        (
            d0 + q_cubic * d3,
            d1 + p_cubic * d3 + q_cubic * d4,
            d2 + p_cubic * d4,
        ),
        axis=-1,
    ) % PRIME


def algebra_power(value: np.ndarray, exponent: int, p_cubic: int, q_cubic: int) -> np.ndarray:
    result = scalar(1)
    factor = value
    while exponent:
        if exponent & 1:
            result = algebra_mul(result, factor, p_cubic, q_cubic)
        exponent //= 2
        if exponent:
            factor = algebra_mul(factor, factor, p_cubic, q_cubic)
    return result


def power_table(value: np.ndarray, maximum: int, p_cubic: int, q_cubic: int) -> np.ndarray:
    rows = [scalar(1)]
    for _ in range(maximum):
        rows.append(algebra_mul(rows[-1], value, p_cubic, q_cubic))
    return np.asarray(rows, dtype=np.int64)


def algebra_fmap(x: np.ndarray, y: np.ndarray, z: np.ndarray, p_cubic: int, q_cubic: int):
    one = scalar(1)
    xy = algebra_mul(x, y, p_cubic, q_cubic)
    unit = (one + xy) % PRIME
    four_plus = (scalar(4) + 3 * xy) % PRIME
    first = (
        algebra_mul(algebra_power(unit, 3, p_cubic, q_cubic), z, p_cubic, q_cubic)
        + algebra_mul(
            algebra_mul(algebra_power(y, 2, p_cubic, q_cubic), unit, p_cubic, q_cubic),
            four_plus,
            p_cubic,
            q_cubic,
        )
    ) % PRIME
    second = (
        y
        + 3
        * algebra_mul(
            algebra_mul(x, algebra_power(unit, 2, p_cubic, q_cubic), p_cubic, q_cubic),
            z,
            p_cubic,
            q_cubic,
        )
        + 3
        * algebra_mul(
            algebra_mul(x, algebra_power(y, 2, p_cubic, q_cubic), p_cubic, q_cubic),
            four_plus,
            p_cubic,
            q_cubic,
        )
    ) % PRIME
    third = (
        2 * x
        - 3 * algebra_mul(algebra_power(x, 2, p_cubic, q_cubic), y, p_cubic, q_cubic)
        - algebra_mul(algebra_power(x, 3, p_cubic, q_cubic), z, p_cubic, q_cubic)
    ) % PRIME
    return first, second, third


def determinant_3x3(matrix: np.ndarray) -> int:
    a00, a01, a02 = map(int, matrix[0])
    a10, a11, a12 = map(int, matrix[1])
    a20, a21, a22 = map(int, matrix[2])
    return (
        a00 * (a11 * a22 - a12 * a21)
        - a01 * (a10 * a22 - a12 * a20)
        + a02 * (a10 * a21 - a11 * a20)
    ) % PRIME


def norm_closed(value: np.ndarray, p_cubic: int, q_cubic: int) -> int:
    # The standard closed formula below is written for
    # w^3+p_polynomial*w+q_polynomial=0, whereas algebra_mul stores
    # w^3=p_cubic*w+q_cubic.
    p_polynomial = -p_cubic % PRIME
    q_polynomial = -q_cubic % PRIME
    h0, h1, h2 = map(int, value)
    return (
        p_polynomial**2 * h2**2 * h0
        - p_polynomial * q_polynomial * h2**2 * h1
        - 2 * p_polynomial * h2 * h0**2
        + p_polynomial * h0 * h1**2
        + q_polynomial**2 * h2**3
        + 3 * q_polynomial * h2 * h0 * h1
        - q_polynomial * h1**3
        + h0**3
    ) % PRIME


def norm_matrix(value: np.ndarray, p_cubic: int, q_cubic: int) -> int:
    basis = np.asarray((scalar(1), (0, 1, 0), (0, 0, 1)), dtype=np.int64)
    columns = algebra_mul(np.broadcast_to(value, basis.shape), basis, p_cubic, q_cubic)
    return determinant_3x3(columns.T)


exponents = np.asarray(list(J_dictionary), dtype=np.int64)
coefficients = np.asarray([J_dictionary[tuple(row)] % PRIME for row in exponents], dtype=np.int64)
maxima = tuple(int(exponents[:, axis].max()) for axis in range(3))
require(maxima == (86, 129, 76), "J multidegree changed")

# At A=infinity on b=c=1, v(w)=2/3, v(q_y)>=-2/3, and v(q_z)>=-2.
# Therefore a J monomial has pole order at most 2(-i+j+3k)/3.
infinity_weight = int(np.max(-exponents[:, 0] + exponents[:, 1] + 3 * exponents[:, 2]))
require(infinity_weight == 228, "J infinity weight changed")
require(2 * infinity_weight % 3 == 0, "slice pole bound unexpectedly fractional")
per_sheet_pole = 2 * infinity_weight // 3
computed_degree_bound = 3 * per_sheet_pole + 2 * 43
require((per_sheet_pole, computed_degree_bound) == (152, DEGREE_BOUND), "degree bound changed")
require(DEGREE_BOUND < PRIME, "finite-field aliasing guard failed")


class EvaluationPlan:
    """Aggregate coefficients along one exponent, grouped by the other two."""

    def __init__(self, remaining_axis: int):
        fixed_axes = tuple(axis for axis in range(3) if axis != remaining_axis)
        order = np.lexsort(
            (
                exponents[:, remaining_axis],
                exponents[:, fixed_axes[1]],
                exponents[:, fixed_axes[0]],
            )
        )
        ordered_fixed = exponents[order][:, fixed_axes]
        change = np.ones(len(order), dtype=bool)
        change[1:] = np.any(ordered_fixed[1:] != ordered_fixed[:-1], axis=1)
        self.starts = np.flatnonzero(change)
        self.fixed_axes = fixed_axes
        self.remaining_axis = remaining_axis
        self.remaining_exponents = exponents[order, remaining_axis]
        self.fixed_exponents = ordered_fixed[self.starts]
        self.coefficients = coefficients[order]
        self.group_count = len(self.starts)

    def evaluate(self, powers: tuple[np.ndarray, np.ndarray, np.ndarray], p_cubic: int, q_cubic: int):
        term_rows = (
            powers[self.remaining_axis][self.remaining_exponents]
            * self.coefficients[:, None]
        ) % PRIME
        grouped = np.add.reduceat(term_rows, self.starts, axis=0) % PRIME
        first = powers[self.fixed_axes[0]][self.fixed_exponents[:, 0]]
        second = powers[self.fixed_axes[1]][self.fixed_exponents[:, 1]]
        grouped = algebra_mul(grouped, first, p_cubic, q_cubic)
        grouped = algebra_mul(grouped, second, p_cubic, q_cubic)
        return np.sum(grouped, axis=0, dtype=np.int64) % PRIME


plans = sorted((EvaluationPlan(axis) for axis in range(3)), key=lambda plan: plan.group_count)
primary_plan, secondary_plan = plans[:2]


def target_data(A: int):
    """Return the inverse cubic and reduced inverse coordinates at (A,1,1)."""

    A %= PRIME
    L = A * (27 * A - 2) % PRIME
    S = (27 * A - 1) % PRIME
    require(L != 0 and S != 0, "singular target entered evaluator")
    inverse_L = mod_inverse(L)
    p_cubic = -inverse_L % PRIME  # T=1: w^3=(-1/L)w+(2/L)
    q_cubic = 2 * inverse_L % PRIME
    K = (9 * A - 1) % PRIME
    M = (27 * A + 17) % PRIME
    Y0 = (9 * A + 1) % PRIME
    A1 = (81 * A - 7) % PRIME
    A2 = (1 - 3 * A) % PRIME
    Z0 = (-9 * A2 - 4 * M * L) % PRIME
    inverse_2S = mod_inverse(2 * S)
    inverse_8S = mod_inverse(8 * S)
    x = np.array((0, 1, 0), dtype=np.int64)
    y = np.array((Y0, 6 * L, -3 * K * L), dtype=np.int64) * inverse_2S % PRIME
    z = np.array((Z0, 6 * L * A1, -9 * L * A2), dtype=np.int64) * inverse_8S % PRIME
    return L, p_cubic, q_cubic, x, y, z


def is_regular_target(A: int) -> bool:
    A %= PRIME
    L = A * (27 * A - 2) % PRIME
    S = (27 * A - 1) % PRIME
    # Discriminant of L*w^3+w-2 is -4*L*(1+27*L).
    return L != 0 and S != 0 and (1 + 27 * L) % PRIME != 0


def evaluate_G(A: int, crosscheck: bool = False) -> int:
    L, p_cubic, q_cubic, x, y, z = target_data(A)
    image = algebra_fmap(x, y, z, p_cubic, q_cubic)
    expected = (scalar(A), scalar(B_VALUE), scalar(C_VALUE))
    require(all(np.array_equal(row, target) for row, target in zip(image, expected)), "inverse graph failed")
    powers = (
        power_table(x, maxima[0], p_cubic, q_cubic),
        power_table(y, maxima[1], p_cubic, q_cubic),
        power_table(z, maxima[2], p_cubic, q_cubic),
    )
    value = primary_plan.evaluate(powers, p_cubic, q_cubic)
    if crosscheck:
        second = secondary_plan.evaluate(powers, p_cubic, q_cubic)
        require(np.array_equal(value, second), "the two J aggregation routes disagree")
    closed = norm_closed(value, p_cubic, q_cubic)
    matrix = norm_matrix(value, p_cubic, q_cubic)
    require(closed == matrix, "closed and matrix cubic norms disagree")
    return pow(L, 43, PRIME) * closed % PRIME


def trim(polynomial: list[int]) -> list[int]:
    result = [value % PRIME for value in polynomial]
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def poly_add_scaled(target: list[int], source: list[int], scale: int) -> None:
    for index, coefficient in enumerate(source):
        target[index] = (target[index] + scale * coefficient) % PRIME


def poly_multiply_linear(polynomial: list[int], root: int) -> list[int]:
    result = [0] * (len(polynomial) + 1)
    for index, coefficient in enumerate(polynomial):
        result[index] = (result[index] - root * coefficient) % PRIME
        result[index + 1] = (result[index + 1] + coefficient) % PRIME
    return trim(result)


def poly_evaluate(polynomial: list[int], value: int) -> int:
    result = 0
    for coefficient in reversed(polynomial):
        result = (result * value + coefficient) % PRIME
    return result


def interpolate(nodes: list[int], values: list[int]) -> list[int]:
    require(len(nodes) == len(values), "interpolation input mismatch")
    divided = [value % PRIME for value in values]
    for order in range(1, len(nodes)):
        for index in range(len(nodes) - 1, order - 1, -1):
            denominator = nodes[index] - nodes[index - order]
            divided[index] = (
                (divided[index] - divided[index - 1]) * mod_inverse(denominator)
            ) % PRIME
    polynomial = [0] * len(nodes)
    basis = [1]
    for order, coefficient in enumerate(divided):
        poly_add_scaled(polynomial, basis, coefficient)
        if order + 1 < len(nodes):
            basis = poly_multiply_linear(basis, nodes[order])
    return trim(polynomial)


regular_targets = [value for value in range(PRIME) if is_regular_target(value)]
require(
    len(regular_targets) >= INTERPOLATION_COUNT + HOLDOUT_COUNT,
    "not enough regular residues for interpolation and holdouts",
)
nodes = regular_targets[:INTERPOLATION_COUNT]
holdouts = regular_targets[INTERPOLATION_COUNT : INTERPOLATION_COUNT + HOLDOUT_COUNT]

print(
    "stage: evaluating 543 interpolation nodes "
    f"(primary remaining axis={primary_plan.remaining_axis}, groups={primary_plan.group_count})",
    flush=True,
)
values = [evaluate_G(value, crosscheck=index < 3) for index, value in enumerate(nodes)]
G_coefficients = interpolate(nodes, values)
require(len(G_coefficients) - 1 <= DEGREE_BOUND, "interpolated G slice exceeded its proof bound")
for value, expected in zip(nodes, values):
    require(poly_evaluate(G_coefficients, value) == expected, "interpolation node mismatch")

print("stage: checking 12 held-out determinants and the second aggregation route", flush=True)
for value in holdouts:
    require(
        poly_evaluate(G_coefficients, value) == evaluate_G(value, crosscheck=True),
        "held-out direct norm disagrees with the interpolation",
    )

# The exact top degree must be visible: one fewer interpolation node is not
# enough, and its result fails immediately at the omitted node.
short_polynomial = interpolate(nodes[:-1], values[:-1])
require(
    poly_evaluate(short_polynomial, nodes[-1]) != values[-1],
    "degree-542 hostile unexpectedly passed with only 542 nodes",
)

field = fmpz_mod_poly_ctx(PRIME)
G_mod = field(G_coefficients)
X_mod = field.gen()


def specialized_polynomial(terms, variable_axis: int = 0):
    degree = max(monomial[variable_axis] for monomial, _coefficient in terms)
    coefficients_out = [0] * (degree + 1)
    for monomial, coefficient in terms:
        # b=c=1, so only the a exponent remains.
        coefficients_out[monomial[variable_axis]] = (
            coefficients_out[monomial[variable_axis]] + coefficient
        ) % PRIME
    return field(trim(coefficients_out))


J_terms = list(J_dictionary.items())
H_mod = specialized_polynomial(H_terms)
J_mod = specialized_polynomial(J_terms)
L_mod = X_mod * (27 * X_mod - 2)

require((H_mod.degree(), J_mod.degree()) == (14, 86), "old slice degrees changed")
require(G_mod.degree() == DEGREE_BOUND, "G slice did not attain the degree-542 bound")
require(G_mod.gcd(G_mod.derivative()).degree() == 0, "G slice is not squarefree")
require(G_mod.gcd(L_mod).degree() == 0, "G slice meets old factor L")
require(G_mod.gcd(H_mod).degree() == 0, "G slice meets old factor H")
require(G_mod.gcd(J_mod).degree() == 0, "G slice meets old factor J")

factor_unit, factor_rows = G_mod.factor()
factor_degrees = sorted((factor.degree(), exponent) for factor, exponent in factor_rows)
require(factor_degrees == EXPECTED_FACTOR_DEGREES, "G slice factor-degree ledger changed")

# Hostile controls prove that the derivative and old-factor gcd tests really
# detect the two failure modes used in the geometric argument.
repeated_hostile = G_mod * G_mod
require(
    repeated_hostile.gcd(repeated_hostile.derivative()).degree() == G_mod.degree(),
    "proper-power hostile escaped the derivative gcd",
)
old_factor_hostile = G_mod * H_mod
require(
    old_factor_hostile.gcd(H_mod).degree() == H_mod.degree(),
    "old-component hostile escaped the gcd test",
)

ledger = "\n".join(f"{index}:{coefficient}" for index, coefficient in enumerate(G_coefficients))
digest = hashlib.sha256(ledger.encode("ascii")).hexdigest()
require(digest == EXPECTED_G_SLICE_SHA256, "G slice coefficient ledger changed")

print(f"J infinity weight={infinity_weight}; per-sheet pole bound={per_sheet_pole}")
print(f"G(A,1,1) degree bound={DEGREE_BOUND}; attained degree={G_mod.degree()}")
print(f"regular residue count={len(regular_targets)}; interpolation/holdout={len(nodes)}/{len(holdouts)}")
print(
    "J aggregation plans: "
    + ", ".join(
        f"remaining-axis-{plan.remaining_axis}:{plan.group_count}-groups" for plan in plans
    )
)
print("closed cubic norm = direct 3x3 multiplication determinant: PASS")
print("twelve held-out direct norms and second J aggregation route: PASS")
print("degree-541 shortened-interpolation hostile: FAILS as required")
print(f"G slice squarefree; factor degree/exponent rows={factor_degrees}")
print(
    "gcd degrees with specialized (L,H,J)="
    f"({G_mod.gcd(L_mod).degree()},{G_mod.gcd(H_mod).degree()},{G_mod.gcd(J_mod).degree()})"
)
print("proper-power and injected-old-H hostiles: detected")
print(f"ascending G-slice coefficient-ledger sha256={digest}")
print(
    "scope: finite-exact mod-1009 slice/multiplicity/coprimality gate for one fixed map; "
    "no global Norm(J) expansion and no general JC conclusion"
)
print("all exact checks passed")
