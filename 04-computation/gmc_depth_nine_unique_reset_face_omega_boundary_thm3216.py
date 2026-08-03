#!/usr/bin/env python3
"""Exact depth-nine unique-reset face and omega boundary for THM-3216."""

import contextlib
import hashlib
import importlib.util
import io
from collections import Counter
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations_with_replacement, product
from math import gcd, lcm
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BASE_SCRIPT = HERE / "gmc_depth_eight_complete_quotient_reset_thm3209.py"
BASE_OUTPUT = ROOT / "05-knowledge/results/gmc_depth_eight_complete_quotient_reset_thm3209.out"
DEPENDENCIES = (
    (BASE_SCRIPT,
     "c04b80a5e4a9cbfc4d7a59fac338147f29582a36751853e0da98752a83288fe8"),
    (BASE_OUTPUT,
     "cb8c20f1a8fffd1b594cf09283ee87269b71638cb8df089db8c74c147717345f"),
)


def lf_hash(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for path, expected in DEPENDENCIES:
    require(lf_hash(path) == expected, ("dependency hash drift", str(path)))


spec = importlib.util.spec_from_file_location("thm3209_exact", BASE_SCRIPT)
base = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(base)


def ones(count):
    return (1,) * count


# Each row is (degree, upset size, minimal generators, normalized LP weight).
# The exact multiplier used below is weight / M_degree.  A common positive
# scalar is irrelevant to the supporting inequality.
CERTIFICATE = (
    (5, 1, ((5,),), Fraction(5, 72755559)),
    (14, 130, ((3, 2) + ones(9), (2, 2, 2) + ones(8)),
     Fraction(1577093, 47102691)),
    (8, 21, ((2,) + ones(6),), Fraction(53234, 56508147)),
    (10, 40, ((3,) + ones(7), (2, 2) + ones(6)),
     Fraction(151685, 66318847)),
    (12, 76, ((2,) + ones(10),), Fraction(30203115, 76749488)),
    (12, 74, ((2, 2) + ones(8),), Fraction(57683, 3902927)),
    (14, 128, ((4,) + ones(10), (3, 2) + ones(9),
               (2, 2, 2, 2, 2, 2, 1, 1)),
     Fraction(283819, 63667394)),
    (14, 121, ((7,) + ones(7), (3, 3, 2) + ones(6),
               (2, 2, 2, 2) + ones(6)),
     Fraction(2861, 13045678)),
    (14, 132, ((2, 2) + ones(10),), Fraction(10395173, 58298789)),
    (14, 129, ((5,) + ones(9), (3, 3) + ones(8),
               (2, 2, 2) + ones(8)),
     Fraction(1951835, 98897507)),
    (14, 127, ((5,) + ones(9), (4, 2) + ones(8),
               (3, 3) + ones(8), (2, 2, 2, 2) + ones(6)),
     Fraction(3979, 1773976)),
    (10, 41, ((2,) + ones(8),), Fraction(4381327, 51259840)),
    (11, 55, ((2,) + ones(9),), Fraction(6670053, 32445565)),
    (9, 29, ((2,) + ones(7),), Fraction(1319165, 92555776)),
    (13, 97, ((3,) + ones(10), (2, 2, 2, 2) + ones(5)),
     Fraction(54743, 8093902)),
    (11, 54, ((3,) + ones(8), (2, 2) + ones(7)),
     Fraction(198899, 26117585)),
    (14, 113, ((7,) + ones(7), (6, 2) + ones(6),
               (5, 3) + ones(6), (5, 2, 2) + ones(5),
               (4, 4) + ones(6), (4, 3, 2) + ones(5),
               (3, 3, 3) + ones(5), (3, 3, 2, 2, 2, 1, 1),
               (2, 2, 2, 2, 2, 2, 1, 1)),
     Fraction(2057, 73733672)),
    (13, 98, ((3,) + ones(10), (2, 2, 2) + ones(7)),
     Fraction(1313225, 43330928)),
)


# These are the exact maximum absolute raw response coordinates on the full
# depth-at-most-nine bank, before the discovery oracle's harmless 10^8 scale.
NORMALIZERS = {
    5: 1300713120,
    6: 591017995680,
    7: 67528413334656,
    8: 9558239814808320,
    9: 1141773245941342464,
    10: 98177407096144199040,
    11: 6402822693369065077104,
    12: 623312546828577253020720,
    13: 39394150378793250375693600,
    14: 2528571670236939601303479360,
}


POLES = tuple(base.POLES)
RESET = tuple(base.RESET)
MULTIPLICITY = Counter(POLES)
VALUES = tuple(sorted(MULTIPLICITY))
LAYERS = tuple(
    tuple(
        state
        for state in combinations_with_replacement(VALUES, depth)
        if all(Counter(state)[value] <= MULTIPLICITY[value]
               for value in set(state))
    )
    for depth in range(1, 10)
)
require(tuple(map(len, LAYERS)) == (8, 33, 93, 200, 348, 507, 631, 678, 631),
        "through-depth-nine state census drift")
STATES = tuple(state for layer in LAYERS for state in layer)
require(len(STATES) == 3129 and len(set(STATES)) == len(STATES),
        "through-depth-nine universe drift")


UPSETS = []
for degree, expected_size, generators, _ in CERTIFICATE:
    bank = base.partitions(degree)
    upset = frozenset(
        shape for shape in bank
        if any(base.coarsens(generator, shape) for generator in generators)
    )
    require(len(upset) == expected_size, ("upset size drift", degree, generators))
    require((degree,) in upset and ones(degree) not in upset,
            ("upset is not nonempty proper", degree, generators))
    require(all(not base.coarsens(fine, coarse) or coarse in upset
                for fine in upset for coarse in bank),
            ("upset lost upward closure", degree, generators))
    UPSETS.append(upset)
UPSETS = tuple(UPSETS)


# Reconstruct every exact response once.  Importing THM-3209 has already
# populated the expensive depth-eight cache, but this check does not rely on
# a serialized discovery cache.
VECTORS = {
    state: base.coefficient_vectors(1, base.BANK, 1, 3, state)
    for state in STATES
}
for degree, normalizer in NORMALIZERS.items():
    observed = max(
        abs(VECTORS[state][degree][shape].numerator)
        for state in STATES for shape in base.partitions(degree)
    )
    require(observed == normalizer, ("degree normalizer drift", degree))
    require(all(VECTORS[state][degree][shape].denominator == 1
                for state in STATES for shape in base.partitions(degree)),
            ("nonintegral response coordinate", degree))


ROWS = tuple(
    tuple(
        sum(VECTORS[state][degree][shape] for shape in upset).numerator
        for state in STATES
    )
    for (degree, _, _, _), upset in zip(CERTIFICATE, UPSETS)
)
RATIONAL_MULTIPLIERS = tuple(
    weight / NORMALIZERS[degree]
    for degree, _, _, weight in CERTIFICATE
)
require(all(weight > 0 for weight in RATIONAL_MULTIPLIERS),
        "dual multiplier lost nonnegativity")


COMMON_DENOMINATOR = reduce(
    lcm, (weight.denominator for weight in RATIONAL_MULTIPLIERS), 1)
PRIMITIVE_MULTIPLIERS = tuple(
    weight.numerator * (COMMON_DENOMINATOR // weight.denominator)
    for weight in RATIONAL_MULTIPLIERS
)
DIVISOR = reduce(gcd, PRIMITIVE_MULTIPLIERS)
PRIMITIVE_MULTIPLIERS = tuple(value // DIVISOR
                              for value in PRIMITIVE_MULTIPLIERS)
require(reduce(gcd, PRIMITIVE_MULTIPLIERS) == 1,
        "cleared dual is not primitive")


COORDINATES = tuple(
    sum(multiplier * row[index]
        for multiplier, row in zip(PRIMITIVE_MULTIPLIERS, ROWS))
    for index in range(len(STATES))
)
ZERO_STATES = tuple(state for state, value in zip(STATES, COORDINATES)
                    if value == 0)
require(ZERO_STATES == (RESET,), "supporting face lost unique reset")
require(all(value < 0 for state, value in zip(STATES, COORDINATES)
            if state != RESET), "supporting cocircuit lost strict negativity")


SIGN_CENSUS_BY_DEPTH = tuple(
    (sum(COORDINATES[sum(map(len, LAYERS[:depth - 1])) + index] < 0
         for index in range(len(LAYERS[depth - 1]))),
     sum(COORDINATES[sum(map(len, LAYERS[:depth - 1])) + index] == 0
         for index in range(len(LAYERS[depth - 1]))))
    for depth in range(1, 10)
)
require(SIGN_CENSUS_BY_DEPTH
        == ((8, 0), (33, 0), (93, 0), (200, 0), (348, 0),
            (507, 0), (631, 0), (677, 1), (631, 0)),
        "depthwise sign census drift")
CLOSEST_INDEX = max(
    (index for index, value in enumerate(COORDINATES) if value < 0),
    key=lambda index: COORDINATES[index],
)
FURTHEST_INDEX = min(
    (index for index, value in enumerate(COORDINATES) if value < 0),
    key=lambda index: COORDINATES[index],
)


# The physical multiset complement is an exact depth-nine/depth-seven
# involution.  In signed-alphabet notation it negates both endpoints.
REMAINDER = tuple(sorted((Counter(POLES) - Counter(RESET)).elements()))
require(REMAINDER == (1, 1, 1, 2, 2, 2, 4, 5),
        "reset complement drift")


def signed_difference(left, right):
    left, right = Counter(left), Counter(right)
    return tuple(left[value] - right[value] for value in VALUES)


for state in LAYERS[8]:
    complement = tuple(sorted((Counter(POLES) - Counter(state)).elements()))
    require(len(complement) == 7, "depth complement drift")
    lhs = signed_difference(RESET, state)
    rhs = tuple(-value for value in signed_difference(REMAINDER, complement))
    require(lhs == rhs, ("quotient complement identity drift", state))


# Exact monomial-basis omega probe.  If p_mu=sum A_(mu,lambda)m_lambda,
# omega acts on p_mu by (-1)^(n-length(mu)); conjugating this diagonal action
# by A gives its integral monomial matrix.
def p_to_m_coefficient(mu, lam):
    count = 0
    for assignment in product(range(len(lam)), repeat=len(mu)):
        sums = [0] * len(lam)
        for part, slot in zip(mu, assignment):
            sums[slot] += part
        if tuple(sums) == tuple(lam):
            count += 1
    return count


def invert(matrix):
    size = len(matrix)
    work = [
        [Fraction(value) for value in row]
        + [Fraction(int(i == j)) for j in range(size)]
        for i, row in enumerate(matrix)
    ]
    for column in range(size):
        pivot = next(row for row in range(column, size)
                     if work[row][column])
        work[column], work[pivot] = work[pivot], work[column]
        scale = work[column][column]
        work[column] = [value / scale for value in work[column]]
        for row in range(size):
            if row == column or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [left - scale * right
                         for left, right in zip(work[row], work[column])]
    return tuple(tuple(row[size:]) for row in work)


@lru_cache(maxsize=None)
def omega_matrix(degree):
    bank = base.partitions(degree)
    transition = tuple(
        tuple(p_to_m_coefficient(mu, lam) for lam in bank)
        for mu in bank
    )
    inverse = invert(transition)
    answer = tuple(
        tuple(sum(inverse[i][k] * (-1) ** (degree - len(bank[k]))
                  * transition[k][j] for k in range(len(bank)))
              for j in range(len(bank)))
        for i in range(len(bank))
    )
    require(all(value.denominator == 1 for row in answer for value in row),
            ("omega matrix lost integrality", degree))
    require(all(sum(answer[i][k] * answer[k][j]
                    for k in range(len(bank))) == int(i == j)
                for i in range(len(bank)) for j in range(len(bank))),
            ("omega matrix lost involutivity", degree))
    return bank, answer


def omega_principal_image(degree, generator):
    bank, omega = omega_matrix(degree)
    indicator = tuple(int(base.coarsens(generator, shape)) for shape in bank)
    image = tuple(sum(indicator[i] * omega[i][j]
                      for i in range(len(bank)))
                  for j in range(len(bank)))
    return tuple((shape, int(value)) for shape, value in zip(bank, image)
                 if value)


OMEGA_DEGREE_FOUR = omega_principal_image(4, (2, 2))
OMEGA_DEGREE_SIX = omega_principal_image(6, (3, 1, 1, 1))
require(OMEGA_DEGREE_FOUR == (((2, 2), 1),),
        "degree-four omega cone hostile drift")
require(OMEGA_DEGREE_SIX == (
    ((6,), -1), ((4, 1, 1), 1), ((3, 3), -1),
    ((3, 1, 1, 1), 1)), "degree-six omega sign hostile drift")


primitive_text = ",".join(map(str, PRIMITIVE_MULTIPLIERS)).encode("ascii")
print("THM-3216 depth-nine unique reset face exact control")
print("dependency_hash_checks=" + repr(len(DEPENDENCIES)))
print("poles_reset_remainder=" + repr((POLES, RESET, REMAINDER)))
print("state_counts_depth1_through9_total="
      + repr((tuple(map(len, LAYERS)), len(STATES))))
print("certificate_rows=" + repr(tuple(
    (degree, size, generators, weight)
    for degree, size, generators, weight in CERTIFICATE
)))
print("degree_normalizers=" + repr(NORMALIZERS))
print("primitive_multiplier_count_gcd="
      + repr((len(PRIMITIVE_MULTIPLIERS),
              reduce(gcd, PRIMITIVE_MULTIPLIERS))))
print("primitive_multipliers=" + repr(PRIMITIVE_MULTIPLIERS))
print("primitive_multipliers_sha256="
      + hashlib.sha256(primitive_text).hexdigest())
print("sign_census_by_depth_negative_zero=" + repr(SIGN_CENSUS_BY_DEPTH))
print("unique_zero_state=" + repr(ZERO_STATES))
print("closest_negative_state_coordinate="
      + repr((STATES[CLOSEST_INDEX], COORDINATES[CLOSEST_INDEX])))
print("furthest_negative_state_coordinate="
      + repr((STATES[FURTHEST_INDEX], COORDINATES[FURTHEST_INDEX])))
print("depth9_depth7_complement_checks=" + repr(len(LAYERS[8])))
print("omega_degree4_principal_image=" + repr(OMEGA_DEGREE_FOUR))
print("omega_degree6_principal_image=" + repr(OMEGA_DEGREE_SIX))
print("scope=virtual_selector_face_and_exterior_complement_no_physical_carrier")
print("all_exact_checks=PASS")
