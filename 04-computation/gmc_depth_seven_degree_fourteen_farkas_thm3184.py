#!/usr/bin/env python3
"""Exact even-skeleton Farkas certificate for THM-3184.

For support ``(1,3)``, bank ``I2``, this companion reconstructs all 1,820
physical prefix states through depth seven.  Ten explicitly generated
partition-coarsening upsets in degrees 8, 10, 12, and 14 have a positive
rational combination which is strictly negative on every state.  Therefore
no probability law can make every nontrivial upset response nonnegative
through degree 14.

The floating cutting-plane/MILP searches which discovered and compressed the
certificate are not used below.  Every promoted check uses integers and
``Fraction`` only.
"""

import ast
import hashlib
from collections import Counter
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations_with_replacement
from math import gcd, lcm
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_pole_prefix_hasse_current_scout.py"
UPSTREAM_BANK = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
UPSTREAM_LF_SHA256 = (
    "151edb9b8cee4807d3f08dc17af32e45420021ba30dfd116c04d9fcaf8bbd5b7")
UPSTREAM_BANK_LF_SHA256 = (
    "15b94691d53afbcdc6aefda89fc7cd5497534ca70fb780a686892dabb5646d6f")


def lf_bytes(path):
    return path.read_bytes().replace(b"\r\n", b"\n")


def load_upstream_prefix(maximum_degree):
    actual = hashlib.sha256(lf_bytes(UPSTREAM)).hexdigest()
    require(actual == UPSTREAM_LF_SHA256,
            ("pole-prefix helper hash drift", actual))
    bank_actual = hashlib.sha256(lf_bytes(UPSTREAM_BANK)).hexdigest()
    require(bank_actual == UPSTREAM_BANK_LF_SHA256,
            ("signed-bank helper hash drift", bank_actual))
    tree = ast.parse(UPSTREAM.read_text(encoding="utf-8"))
    prefix = []
    for node in tree.body:
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name)
                        and target.id == "MAXIMUM_DEGREE"
                        for target in node.targets)):
            node.value = ast.Constant(maximum_degree)
        prefix.append(node)
        if (isinstance(node, ast.Assign)
                and any(isinstance(target, ast.Name)
                        and target.id == "UNIVERSE"
                        for target in node.targets)):
            break
    module = ast.fix_missing_locations(
        ast.Module(body=prefix, type_ignores=[]))
    namespace = {"__file__": str(UPSTREAM)}
    exec(compile(module, str(UPSTREAM), "exec"), namespace)
    return namespace


UP = load_upstream_prefix(14)
UP["residual_roots"] = lru_cache(maxsize=None)(UP["residual_roots"])
UP["all_monomial_values"] = lru_cache(maxsize=None)(
    UP["all_monomial_values"])
upstream_partitions = UP["partitions"]
upstream_coarsens = UP["coarsens"]
coefficient_vectors = UP["coefficient_vectors"]
reduced_poles = UP["reduced_poles"]
BANK = UP["BANKS"][1]


# ---------------------------------------------------------------------------
# 1. Independent partition/refinement implementation.
# ---------------------------------------------------------------------------

@lru_cache(maxsize=None)
def independent_partitions(total, maximum=None):
    if total == 0:
        return ((),)
    if maximum is None or maximum > total:
        maximum = total
    answer = []
    for first in range(maximum, 0, -1):
        for tail in independent_partitions(total - first, first):
            answer.append((first,) + tail)
    return tuple(answer)


@lru_cache(maxsize=None)
def independent_coarsens(fine, coarse):
    """Whether ``fine`` packs exactly into bins ``coarse``."""

    if sum(fine) != sum(coarse) or len(fine) < len(coarse):
        return False
    pieces = tuple(sorted(fine, reverse=True))
    bins = tuple(sorted(coarse, reverse=True))

    @lru_cache(maxsize=None)
    def place(index, capacities):
        if index == len(pieces):
            return not any(capacities)
        piece = pieces[index]
        tried = set()
        for slot, capacity in enumerate(capacities):
            if capacity < piece or capacity in tried:
                continue
            tried.add(capacity)
            updated = list(capacities)
            updated[slot] -= piece
            updated.sort(reverse=True)
            if place(index + 1, tuple(updated)):
                return True
        return False

    return place(0, bins)


PARTITIONS = {
    degree: independent_partitions(degree)
    for degree in range(5, 15)
}
PARTITION_CHECKS = 0
COARSENING_CHECKS = 0
for degree, shapes in PARTITIONS.items():
    require(set(shapes) == set(upstream_partitions(degree)),
            ("independent partition bank mismatch", degree))
    PARTITION_CHECKS += len(shapes)
    for fine in shapes:
        for coarse in shapes:
            require(independent_coarsens(fine, coarse)
                    == upstream_coarsens(fine, coarse),
                    ("independent coarsening mismatch", degree, fine, coarse))
            COARSENING_CHECKS += 1


# ---------------------------------------------------------------------------
# 2. Complete depth-seven bank and the readable normalized dual.
# ---------------------------------------------------------------------------

POLES, _ = reduced_poles(1, BANK, 1, 3)
MULTIPLICITY = Counter(POLES)
VALUES = tuple(sorted(MULTIPLICITY))
BY_DEPTH = tuple(
    tuple(
        state
        for state in combinations_with_replacement(VALUES, depth)
        if all(Counter(state)[value] <= MULTIPLICITY[value]
               for value in set(state))
    )
    for depth in range(1, 8)
)
STATES = tuple(state for layer in BY_DEPTH for state in layer)
require(POLES == (8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1),
        "support-(1,3) pole multiset drift")
require(tuple(map(len, BY_DEPTH)) == (8, 33, 93, 200, 348, 507, 631)
        and len(STATES) == 1820,
        "depth-seven physical-state census drift")


def ones(count):
    return (1,) * count


# (degree, expected upset size, minimal-generator antichain, q_i, M_i).
# The Farkas functional is sum_i q_i r_i/M_i.
CERTIFICATE = (
    (14, 130, ((3, 2) + ones(9), (2, 2, 2) + ones(8)),
     Fraction(429, 1354), 92609810408824812936364032),
    (8, 21, ((2,) + ones(6),),
     Fraction(101, 28748), 54841362155328),
    (10, 40, ((3,) + ones(7), (2, 2) + ones(6)),
     Fraction(507, 24251), 1093237378467812904),
    (12, 76, ((2,) + ones(10),),
     Fraction(203, 3189), 732937666219205291328),
    (12, 74, ((2, 2) + ones(8),),
     Fraction(915, 24238), 1773893816717680035480),
    (14, 128, ((4,) + ones(10), (3, 2) + ones(9),
               (2, 2, 2, 2, 2, 2, 1, 1)),
     Fraction(117, 26686), 11449587449741073916566528),
    (14, 121, ((7,) + ones(7), (3, 3, 2) + ones(6),
               (2, 2, 2, 2) + ones(6)),
     Fraction(30, 29581), 158608091980421055669473280),
    (14, 132, ((2, 2) + ones(10),),
     Fraction(913, 5536), 6483094663875178504128552),
    (14, 129, ((5,) + ones(9), (3, 3) + ones(8),
               (2, 2, 2) + ones(8)),
     Fraction(5792, 19469), 68383863660622477237420032),
    (14, 127, ((5,) + ones(9), (4, 2) + ones(8),
               (3, 3) + ones(8), (2, 2, 2, 2) + ones(6)),
     Fraction(1999, 22331), 136981330444844093510221824),
)

UPSETS = []
MINIMAL_GENERATOR_CHECKS = 0
for degree, expected_size, generators, _, _ in CERTIFICATE:
    shapes = PARTITIONS[degree]
    require(all(generator in shapes for generator in generators),
            ("illegal upset generator", degree, generators))
    upset = frozenset(
        shape for shape in shapes
        if any(independent_coarsens(generator, shape)
               for generator in generators)
    )
    require(len(upset) == expected_size,
            ("upset size drift", degree, len(upset), expected_size))
    require((degree,) in upset and ones(degree) not in upset,
            ("upset lost nontrivial anchors", degree, generators))
    require(all(not independent_coarsens(fine, coarse) or coarse in upset
                for fine in upset for coarse in shapes),
            ("generated set is not an upset", degree, generators))
    minimal = tuple(
        shape for shape in shapes
        if shape in upset and not any(
            other != shape and other in upset
            and independent_coarsens(other, shape)
            for other in shapes
        )
    )
    require(set(minimal) == set(generators)
            and len(minimal) == len(generators),
            ("displayed generators are not the full minimal antichain",
             degree, generators, minimal))
    MINIMAL_GENERATOR_CHECKS += 1
    UPSETS.append(upset)


# ---------------------------------------------------------------------------
# 3. Rebuild every response coordinate and verify the Farkas inequality.
# ---------------------------------------------------------------------------

ROWS = [[] for _ in CERTIFICATE]
BALANCE_CHECKS = 0
for state in STATES:
    vectors = coefficient_vectors(1, BANK, 1, 3, state)
    for degree in range(5, 15):
        require(sum(vectors[degree].values()) == 0,
                ("state response lost zero mass", state, degree))
        BALANCE_CHECKS += 1
    for index, ((degree, _, _, _, _), upset) in enumerate(
            zip(CERTIFICATE, UPSETS)):
        value = sum(vectors[degree][shape] for shape in upset)
        require(value.denominator == 1,
                ("response row lost integrality", index, state))
        ROWS[index].append(value.numerator)
ROWS = tuple(tuple(row) for row in ROWS)

for index, ((_, _, _, _, expected_scale), row) in enumerate(
        zip(CERTIFICATE, ROWS)):
    actual_scale = max(abs(value) for value in row)
    require(actual_scale == expected_scale,
            ("row sup-norm drift", index, actual_scale, expected_scale))

NORMALIZED = tuple(
    sum(q * Fraction(row[state_index], scale)
        for (_, _, _, q, scale), row in zip(CERTIFICATE, ROWS))
    for state_index in range(len(STATES))
)
NORMALIZED_RANGE = (min(NORMALIZED), max(NORMALIZED))
require(NORMALIZED_RANGE[1] <= -Fraction(1, 10**11),
        ("normalized Farkas margin lost", NORMALIZED_RANGE[1]))

# Clear the readable rational/scaled form to a primitive nonnegative integer
# dual and independently replay the strict coordinate inequality.
RAW_MULTIPLIERS = tuple(
    q / scale for _, _, _, q, scale in CERTIFICATE
)
COMMON_DENOMINATOR = reduce(
    lcm, (value.denominator for value in RAW_MULTIPLIERS), 1)
INTEGER_MULTIPLIERS = tuple(
    value.numerator * (COMMON_DENOMINATOR // value.denominator)
    for value in RAW_MULTIPLIERS
)
MULTIPLIER_GCD = reduce(gcd, INTEGER_MULTIPLIERS, 0)
INTEGER_MULTIPLIERS = tuple(
    value // MULTIPLIER_GCD for value in INTEGER_MULTIPLIERS
)
require(all(value > 0 for value in INTEGER_MULTIPLIERS)
        and reduce(gcd, INTEGER_MULTIPLIERS, 0) == 1,
        "integer Farkas multipliers lost positivity or primitivity")

COMBINED = tuple(
    sum(multiplier * row[state_index]
        for multiplier, row in zip(INTEGER_MULTIPLIERS, ROWS))
    for state_index in range(len(STATES))
)
COMBINED_RANGE = (min(COMBINED), max(COMBINED))
require(COMBINED_RANGE[1] < 0,
        ("primitive integer Farkas row is not strictly negative",
         COMBINED_RANGE))


# ---------------------------------------------------------------------------
# 4. Anatomy of the first degree-fourteen hard face.
# ---------------------------------------------------------------------------

HARD_COMPLEMENT = frozenset(
    shape for shape in PARTITIONS[14] if shape not in UPSETS[0]
)
EXPECTED_HARD_COMPLEMENT = frozenset((
    ones(14), (2,) + ones(12), (3,) + ones(11),
    (2, 2) + ones(10), (4,) + ones(10)))
require(HARD_COMPLEMENT == EXPECTED_HARD_COMPLEMENT,
    ("hard-face complement drift", HARD_COMPLEMENT))
HARD_POSITIVE = tuple(
    (state, value) for state, value in zip(STATES, ROWS[0]) if value > 0
)
HARD_NEGATIVE_COUNT = sum(value < 0 for value in ROWS[0])
HARD_ZERO_COUNT = sum(value == 0 for value in ROWS[0])
HARD_POSITIVE_BY_DEPTH = tuple(
    sum(value > 0 for state, value in HARD_POSITIVE if len(state) == depth)
    for depth in range(1, 8)
)
require(len(HARD_POSITIVE) == 30 and HARD_NEGATIVE_COUNT == 1790
        and HARD_ZERO_COUNT == 0
        and HARD_POSITIVE_BY_DEPTH == (1, 2, 4, 3, 3, 5, 12),
        "hard-face sign census drift")
HARD_POSITIVE_SHA256 = hashlib.sha256(
    repr(HARD_POSITIVE).encode("ascii")).hexdigest()


print("THM-3184 depth-seven degree-fourteen even-skeleton Farkas death")
print("pole_prefix_dependency_sha256=" + UPSTREAM_LF_SHA256)
print("signed_bank_dependency_sha256=" + UPSTREAM_BANK_LF_SHA256)
print("poles=" + repr(POLES))
print("state_counts_depth1_depth7_total="
      + repr((tuple(map(len, BY_DEPTH)), len(STATES))))
print("independent_partition_checks=" + repr(PARTITION_CHECKS))
print("independent_coarsening_checks=" + repr(COARSENING_CHECKS))
print("state_degree_zero_mass_checks=" + repr(BALANCE_CHECKS))
print("minimal_generator_antichain_checks="
      + repr(MINIMAL_GENERATOR_CHECKS))
print("certificate_degrees_sizes_generators=" + repr(tuple(
    (degree, size, generators)
    for degree, size, generators, _, _ in CERTIFICATE)))
print("normalized_multipliers=" + repr(tuple(
    q for _, _, _, q, _ in CERTIFICATE)))
print("row_sup_norms=" + repr(tuple(
    scale for _, _, _, _, scale in CERTIFICATE)))
print("normalized_range=" + repr(NORMALIZED_RANGE))
print("normalized_uniform_bound_le=-1/100000000000")
print("primitive_integer_multipliers=" + repr(INTEGER_MULTIPLIERS))
print("primitive_combined_range=" + repr(COMBINED_RANGE))
print("hard_face_complement=" + repr(tuple(sorted(HARD_COMPLEMENT))))
print("hard_face_sign_census="
      + repr((len(HARD_POSITIVE), HARD_NEGATIVE_COUNT, HARD_ZERO_COUNT,
              HARD_POSITIVE_BY_DEPTH)))
print("hard_face_positive_sha256=" + HARD_POSITIVE_SHA256)
print("barcode=(C13_depth7_nonempty,C14_depth7_empty)")
print("scope=fixed_support_bank_averaged_virtual_prefix_currents")
print("all_exact_checks=PASS")
