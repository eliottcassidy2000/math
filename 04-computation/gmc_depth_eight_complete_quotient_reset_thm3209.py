#!/usr/bin/env python3
"""Exact quotient-reset and tangent-star controls for THM-3209."""

import ast
import hashlib
from collections import Counter
from fractions import Fraction
from functools import lru_cache, reduce
from itertools import combinations_with_replacement
from math import factorial
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
UPSTREAM = HERE / "gmc_pole_prefix_hasse_current_scout.py"
UPSTREAM_BANK = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
THM3184_SCRIPT = HERE / "gmc_depth_seven_degree_fourteen_farkas_thm3184.py"
THM3184_OUTPUT = ROOT / "05-knowledge/results/gmc_depth_seven_degree_fourteen_farkas_thm3184.out"
DEPENDENCIES = (
    (UPSTREAM,
     "151edb9b8cee4807d3f08dc17af32e45420021ba30dfd116c04d9fcaf8bbd5b7"),
    (UPSTREAM_BANK,
     "15b94691d53afbcdc6aefda89fc7cd5497534ca70fb780a686892dabb5646d6f"),
    (THM3184_SCRIPT,
     "9d560fc04fd7772a0b92b5973cd85b7ed61a2ec76328badfe5420ae2fcb7d129"),
    (THM3184_OUTPUT,
     "d599057d44b12cd222a0d6ecf3755cb5d2ae59d1faa1f7bd6731ad295ab667b3"),
)


def lf_hash(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for path, expected in DEPENDENCIES:
    require(lf_hash(path) == expected, ("dependency hash drift", str(path)))


def load_upstream_prefix(maximum_degree):
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
    module = ast.fix_missing_locations(ast.Module(body=prefix, type_ignores=[]))
    namespace = {"__file__": str(UPSTREAM)}
    exec(compile(module, str(UPSTREAM), "exec"), namespace)
    return namespace


UP = load_upstream_prefix(14)
UP["residual_roots"] = lru_cache(maxsize=None)(UP["residual_roots"])
UP["all_monomial_values"] = lru_cache(maxsize=None)(UP["all_monomial_values"])
all_monomial_values = UP["all_monomial_values"]
coefficient_vectors = UP["coefficient_vectors"]
coarsest_components = UP["coarsest_components"]
dominant_row = UP["dominant_row"]
reduced_poles = UP["reduced_poles"]
residual_roots = UP["residual_roots"]
BANK = UP["BANKS"][1]

# A separate degree-five namespace makes the complete through-depth-eight
# uniqueness check cheap instead of rebuilding all ten degrees on the 1,820
# already-covered shallower states.
UP5 = load_upstream_prefix(5)
UP5["residual_roots"] = lru_cache(maxsize=None)(UP5["residual_roots"])
UP5["all_monomial_values"] = lru_cache(maxsize=None)(UP5["all_monomial_values"])
coefficient_vectors_degree_five = UP5["coefficient_vectors"]
BANK5 = UP5["BANKS"][1]


@lru_cache(maxsize=None)
def partitions(total, maximum=None):
    if total == 0:
        return ((),)
    if maximum is None or maximum > total:
        maximum = total
    answer = []
    for first in range(maximum, 0, -1):
        for tail in partitions(total - first, first):
            answer.append((first,) + tail)
    return tuple(answer)


@lru_cache(maxsize=None)
def coarsens(fine, coarse):
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
            changed = list(capacities)
            changed[slot] -= piece
            changed.sort(reverse=True)
            if place(index + 1, tuple(changed)):
                return True
        return False

    return place(0, bins)


def ones(count):
    return (1,) * count


# THM-3184's ten rows, with the readable normalized weights q_i/M_i.
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
for degree, expected_size, generators, _, _ in CERTIFICATE:
    bank = partitions(degree)
    upset = frozenset(
        shape for shape in bank
        if any(coarsens(generator, shape) for generator in generators)
    )
    require(len(upset) == expected_size, ("upset size drift", degree))
    UPSETS.append(upset)
UPSETS = tuple(UPSETS)


POLES, _ = reduced_poles(1, BANK, 1, 3)
QUOTIENT_ROOTS = residual_roots(1, dominant_row(1), 1, 3)
RESET = (1, 3, 3, 4, 5, 6, 7, 8)
require(QUOTIENT_ROOTS == RESET, "dominant quotient roots lost reset state")
require(POLES == (8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2,
                  1, 1, 1, 1), "reduced pole multiset drift")


MULTIPLICITY = Counter(POLES)
VALUES = tuple(sorted(MULTIPLICITY))
DEPTH_EIGHT = tuple(
    state
    for state in combinations_with_replacement(VALUES, 8)
    if all(Counter(state)[value] <= MULTIPLICITY[value]
           for value in set(state))
)
require(len(DEPTH_EIGHT) == 678 and RESET in DEPTH_EIGHT,
        "depth-eight state bank drift")
SHALLOW_BY_DEPTH = tuple(
    tuple(
        state
        for state in combinations_with_replacement(VALUES, depth)
        if all(Counter(state)[value] <= MULTIPLICITY[value]
               for value in set(state))
    )
    for depth in range(1, 8)
)
require(tuple(map(len, SHALLOW_BY_DEPTH)) == (8, 33, 93, 200, 348, 507, 631),
        "shallower state bank drift")
SHALLOW_DEGREE_FIVE_ZERO_STATES = tuple(
    state
    for layer in SHALLOW_BY_DEPTH
    for state in layer
    if all(value == 0 for value in coefficient_vectors_degree_five(
        1, BANK5, 1, 3, state)[5].values())
)
require(not SHALLOW_DEGREE_FIVE_ZERO_STATES,
        "a shallower state acquired zero complete degree-five response")


# The complete quotient virtual alphabet is empty at RESET.  The finite check
# through degree fourteen mirrors the all-degree lambda-ring proof.
RESET_QUOTIENT = all_monomial_values(QUOTIENT_ROOTS, RESET)
for degree in range(1, 15):
    require(all(RESET_QUOTIENT[shape] == 0
                for shape in partitions(degree)),
            ("empty quotient alphabet acquired a positive-degree monomial",
             degree))
RESET_VECTORS = coefficient_vectors(1, BANK, 1, 3, RESET)
require(all(all(value == 0 for value in RESET_VECTORS[degree].values())
            for degree in range(5, 15)),
        "reset response failed to vanish through degree fourteen")


# Extend THM-3184's cocircuit to the new depth-eight layer.  It has one and
# only one zero there, and the complete degree-five vector has the same unique
# zero.  THM-3184 already makes the cocircuit strict on every shallower state.
COCIRCUIT_SIGNS = Counter()
COCIRCUIT_ZERO_STATES = []
DEGREE_FIVE_ZERO_STATES = []
DEGREE_FIVE_COARSEST = {}
for state in DEPTH_EIGHT:
    vectors = coefficient_vectors(1, BANK, 1, 3, state)
    degree_five = vectors[5]
    DEGREE_FIVE_COARSEST[state] = degree_five[(5,)]
    if all(value == 0 for value in degree_five.values()):
        DEGREE_FIVE_ZERO_STATES.append(state)
    coordinate = sum(
        (q / scale) * sum(vectors[degree][shape] for shape in upset)
        for (degree, _, _, q, scale), upset in zip(CERTIFICATE, UPSETS)
    )
    sign = (coordinate > 0) - (coordinate < 0)
    COCIRCUIT_SIGNS[sign] += 1
    if not coordinate:
        COCIRCUIT_ZERO_STATES.append(state)

require((COCIRCUIT_SIGNS[1], COCIRCUIT_SIGNS[-1], COCIRCUIT_SIGNS[0])
        == (15, 662, 1), "depth-eight cocircuit sign census drift")
require(tuple(COCIRCUIT_ZERO_STATES) == (RESET,),
        "THM-3184 cocircuit acquired another depth-eight zero")
require(tuple(DEGREE_FIVE_ZERO_STATES) == (RESET,),
        "complete degree-five response acquired another depth-eight zero")
require(SHALLOW_DEGREE_FIVE_ZERO_STATES
        + tuple(DEGREE_FIVE_ZERO_STATES) == (RESET,),
        "complete degree-five zero lost through-depth-eight uniqueness")


# A legal one-exchange hostile makes the reset visibly isolated already in the
# coarsest degree-five coordinate.
reset_counter = Counter(RESET)
ONE_EXCHANGE = set()
for removed in reset_counter:
    for added in MULTIPLICITY:
        candidate = reset_counter.copy()
        candidate[removed] -= 1
        candidate[added] += 1
        if (candidate != reset_counter
                and all(candidate[value] <= MULTIPLICITY[value]
                        for value in candidate)):
            ONE_EXCHANGE.add(tuple(sorted(candidate.elements())))
require(len(ONE_EXCHANGE) == 25, "one-exchange neighbour census drift")
require(all(DEGREE_FIVE_COARSEST[state] != 0 for state in ONE_EXCHANGE),
        "one-exchange neighbour remained reset in degree five")
SHARP_NEIGHBOUR = (2, 3, 3, 4, 5, 6, 7, 8)
SHARP_NEIGHBOUR_RESPONSE = DEGREE_FIVE_COARSEST[SHARP_NEIGHBOUR]
require(SHARP_NEIGHBOUR_RESPONSE == -44640
        and min(abs(DEGREE_FIVE_COARSEST[state]) for state in ONE_EXCHANGE)
        == 44640, "sharp one-exchange hostile drift")
neighbour_components = coarsest_components(
    1, BANK, 1, 3, SHARP_NEIGHBOUR, 5)
require(neighbour_components == (1440, -1, 0, -31),
        "sharp neighbour quotient calculation drift")


# The first tangent star consists of adding one still-available pole.  Its
# quotient alphabet is the negative singleton -{r}; verify the universal
# partition vector and the degree-five principal-upset obstruction.
REMAINDER = MULTIPLICITY - Counter(RESET)
LEGAL_TANGENT_R = tuple(sorted(REMAINDER))
require(LEGAL_TANGENT_R == (1, 2, 4, 5),
        "legal tangent pole types drift")
TANGENT_PHI_FIVE = []
TANGENT_COARSEST_FIVE = []
NEGATIVE_SINGLETON_CHECKS = 0
for root in LEGAL_TANGENT_R:
    state = tuple(sorted(RESET + (root,)))
    quotient = all_monomial_values(QUOTIENT_ROOTS, state)
    vectors = coefficient_vectors(1, BANK, 1, 3, state)
    for degree in range(2, 15):
        require(sum(quotient[shape] for shape in partitions(degree)) == 0,
                ("negative-singleton complete function did not vanish",
                 root, degree))
        phi_h = None
        if degree >= 5:
            phi_h = coarsest_components(
                1, BANK, 1, 3, state, degree)[0]
        for shape in partitions(degree):
            counts = Counter(shape)
            denominator = reduce(
                lambda left, right: left * right,
                (factorial(value) for value in counts.values()), 1)
            expected_monomial = Fraction(
                (-1) ** len(shape) * factorial(len(shape)), denominator
            ) * root**degree
            require(quotient[shape] == expected_monomial,
                    ("negative-singleton monomial formula drift",
                     root, degree, shape))
            if degree >= 5:
                require(vectors[degree][shape]
                        == phi_h * expected_monomial,
                        ("negative-singleton response formula drift",
                         root, degree, shape))
            NEGATIVE_SINGLETON_CHECKS += 1
    phi_five = coarsest_components(1, BANK, 1, 3, state, 5)[0]
    coarsest = vectors[5][(5,)]
    TANGENT_PHI_FIVE.append((root, phi_five))
    TANGENT_COARSEST_FIVE.append((root, coarsest))

require(tuple(TANGENT_PHI_FIVE)
        == ((1, 1440), (2, 1440), (4, 1440), (5, 1440)),
        "tangent Phi(h5) values drift")
require(tuple(TANGENT_COARSEST_FIVE) == (
    (1, -1440), (2, -46080), (4, -1474560), (5, -4500000)),
    "tangent principal-upset responses drift")


print("THM-3209 depth-eight complete quotient reset exact control")
print("dependency_hash_checks=" + repr(len(DEPENDENCIES)))
print("poles=" + repr(POLES))
print("dominant_quotient_roots=" + repr(QUOTIENT_ROOTS))
print("depth_eight_state_count=" + repr(len(DEPTH_EIGHT)))
print("state_counts_depth1_depth8_total="
      + repr((tuple(map(len, SHALLOW_BY_DEPTH)) + (len(DEPTH_EIGHT),),
              sum(map(len, SHALLOW_BY_DEPTH)) + len(DEPTH_EIGHT))))
print("all_degree_reset_state=" + repr(RESET))
print("empty_quotient_monomial_checks_through14="
      + repr(sum(len(partitions(degree)) for degree in range(1, 15))))
print("reset_response_checks_degrees5_through14=PASS")
print("depth_eight_cocircuit_signs_positive_negative_zero="
      + repr((COCIRCUIT_SIGNS[1], COCIRCUIT_SIGNS[-1],
              COCIRCUIT_SIGNS[0])))
print("depth_eight_cocircuit_unique_zero="
      + repr(tuple(COCIRCUIT_ZERO_STATES)))
print("depth_eight_degree_five_unique_zero="
      + repr(tuple(DEGREE_FIVE_ZERO_STATES)))
print("through_depth_eight_degree_five_unique_zero="
      + repr(SHALLOW_DEGREE_FIVE_ZERO_STATES
             + tuple(DEGREE_FIVE_ZERO_STATES)))
print("one_exchange_neighbour_count=" + repr(len(ONE_EXCHANGE)))
print("sharp_one_exchange_state_response_components="
      + repr((SHARP_NEIGHBOUR, SHARP_NEIGHBOUR_RESPONSE,
              neighbour_components)))
print("legal_negative_singleton_tangent_roots=" + repr(LEGAL_TANGENT_R))
print("negative_singleton_formula_checks="
      + repr(NEGATIVE_SINGLETON_CHECKS))
print("tangent_phi_h5=" + repr(tuple(TANGENT_PHI_FIVE)))
print("tangent_coarsest_degree_five_responses="
      + repr(tuple(TANGENT_COARSEST_FIVE)))
print("tangent_star_degree_five_principal_upset=STRICTLY_NEGATIVE")
print("scope=reset_and_first_tangent_star_only_no_depth9_global_claim")
print("all_exact_checks=PASS")
