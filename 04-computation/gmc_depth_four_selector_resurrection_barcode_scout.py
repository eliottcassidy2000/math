#!/usr/bin/env python3
"""Exact depth-four selector resurrection/death barcode.

For support ``(1,3)``, bank ``I2``, this companion proves both sides of the
sharp depth-four cell:

* a four-state probability law is strictly positive on every nontrivial
  coarsening upset in degrees 5 through 11;
* a positive combination of five upset rows is strictly negative on every
  one of the 334 physical prefix states of depth at most four once degree 12
  is included.

Only integer and ``Fraction`` arithmetic is used.  Floating LPs were useful
for discovery but are not part of the promoted evidence.
"""

import ast
from collections import Counter
from fractions import Fraction
from functools import reduce
from itertools import combinations_with_replacement
from math import gcd
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_pole_prefix_hasse_current_scout.py"


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
    module = ast.fix_missing_locations(
        ast.Module(body=prefix, type_ignores=[]))
    namespace = {"__file__": str(UPSTREAM)}
    exec(compile(module, str(UPSTREAM), "exec"), namespace)
    return namespace


UP = load_upstream_prefix(12)
partitions = UP["partitions"]
coarsens = UP["coarsens"]
coefficient_vectors = UP["coefficient_vectors"]
reduced_poles = UP["reduced_poles"]
transport_diagnostic = UP["transport_diagnostic"]
BANK = UP["BANKS"][1]


def all_upsets(degree):
    shapes = partitions(degree)
    antichains = []

    def extend(start, chosen):
        antichains.append(chosen)
        for index in range(start, len(shapes)):
            candidate = shapes[index]
            if all(not coarsens(candidate, old)
                   and not coarsens(old, candidate) for old in chosen):
                extend(index + 1, chosen + (candidate,))

    extend(0, ())
    answer = {
        frozenset(
            high for low in antichain for high in shapes
            if coarsens(low, high)
        )
        for antichain in antichains
    }
    require(len(answer) == len(antichains),
            "antichain/upset correspondence lost injectivity")
    return tuple(sorted(answer, key=lambda upset: (
        len(upset), tuple(sorted(upset)))))


UPSETS = {degree: all_upsets(degree) for degree in range(5, 12)}
require(tuple(len(UPSETS[degree]) for degree in range(5, 12))
        == (10, 27, 47, 168, 573, 3588, 19542),
        "degree-5-through-11 upset census drift")


POLES, _ = reduced_poles(1, BANK, 1, 3)
MULTIPLICITY = Counter(POLES)
VALUES = tuple(sorted(MULTIPLICITY))
BY_DEPTH = []
for depth in range(1, 5):
    BY_DEPTH.append(tuple(
        state
        for state in combinations_with_replacement(VALUES, depth)
        if all(Counter(state)[value] <= MULTIPLICITY[value]
               for value in set(state))
    ))
BY_DEPTH = tuple(BY_DEPTH)
STATES = tuple(state for layer in BY_DEPTH for state in layer)
require(POLES == (8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1),
        "support-(1,3) pole multiset drift")
require(tuple(map(len, BY_DEPTH)) == (8, 33, 93, 200)
        and len(STATES) == 334,
        "depth-four physical-state census drift")

VECTORS = {
    state: coefficient_vectors(1, BANK, 1, 3, state)
    for state in STATES
}


def current(law, degree):
    return {
        shape: sum(weight * VECTORS[state][degree][shape]
                   for state, weight in law.items())
        for shape in partitions(degree)
    }


def complement_upset(degree, excluded):
    return frozenset(shape for shape in partitions(degree)
                     if shape not in excluded)


def response_row(degree, upset):
    row = tuple(
        sum(VECTORS[state][degree][shape] for shape in upset)
        for state in STATES
    )
    require(all(value.denominator == 1 for value in row),
            "response row lost integrality")
    return tuple(value.numerator for value in row)


# ---------------------------------------------------------------------------
# 1. Strict resurrection through degree eleven.
# ---------------------------------------------------------------------------

LAW = {
    (1,): Fraction(2291, 5000),
    (7,): Fraction(36, 5000),
    (1, 1, 1, 1): Fraction(1999, 5000),
    (1, 1, 1, 2): Fraction(674, 5000),
}
require(sum(LAW.values()) == 1 and all(value > 0 for value in LAW.values()),
        "depth-four law left the probability simplex")

EXPECTED_COUNTS = (8, 25, 45, 166, 571, 3586, 19540)
EXPECTED_MINIMA = (
    Fraction(11134658556, 125),
    Fraction(3281855316351, 125),
    Fraction(268382847153552, 125),
    Fraction(11849576869236, 625),
    Fraction(23788552125790506, 625),
    Fraction(42037071115846938, 625),
    Fraction(2271700518167474796, 625),
)
EXPECTED_FLOW_DEMAND = (
    Fraction(85514905026, 125),
    Fraction(388192995599559, 625),
    Fraction(20730977643207087, 125),
    Fraction(17166140949873362619, 625),
    Fraction(2238552161748534218751, 625),
    Fraction(49931929549172056006632, 125),
    Fraction(24924117539388030078404892, 625),
)

NONZERO_COUNTS = []
MINIMA = []
MINIMUM_TAGS = []
FLOW_DEMAND = []
for degree in range(5, 12):
    mixed = current(LAW, degree)
    diagnostic = transport_diagnostic(mixed)
    require(diagnostic[0] and diagnostic[1] == diagnostic[2],
            f"depth-four law failed exact transport at degree {degree}")
    FLOW_DEMAND.append(diagnostic[1])
    all_masses = tuple(
        (sum(mixed[shape] for shape in upset), upset)
        for upset in UPSETS[degree]
    )
    zero_upsets = {upset for mass, upset in all_masses if mass == 0}
    require(zero_upsets == {
        frozenset(), frozenset(partitions(degree))},
        f"degree-{degree} acquired a nontrivial zero upset")
    masses = tuple((mass, upset) for mass, upset in all_masses if mass)
    require(all(mass > 0 for mass, _ in masses),
            f"depth-four law lost strict upset positivity at degree {degree}")
    NONZERO_COUNTS.append(len(masses))
    minimum = min(mass for mass, _ in masses)
    witnesses = tuple(upset for mass, upset in masses if mass == minimum)
    require(len(witnesses) == 1,
            f"degree-{degree} minimum upset lost uniqueness")
    expected_upset = (frozenset({(degree,)}) if degree <= 6 else
                      complement_upset(degree, {(1,) * degree}))
    require(witnesses == (expected_upset,),
            f"degree-{degree} minimum upset drift")
    MINIMA.append(minimum)
    MINIMUM_TAGS.append(
        "top" if degree <= 6 else "singleton-complement")

require(tuple(NONZERO_COUNTS) == EXPECTED_COUNTS,
        "strict nontrivial-upset census drift")
require(tuple(MINIMA) == EXPECTED_MINIMA,
        "per-degree strict margin drift")
require(tuple(FLOW_DEMAND) == EXPECTED_FLOW_DEMAND,
        "per-degree exact flow/demand drift")
require(sum(NONZERO_COUNTS) == 23941,
        "degree-5-through-11 strict-row total drift")


# The same law fails immediately at degree twelve.  This does not by itself
# exclude other depth-four laws; the five-row certificate below does.
CURRENT_12 = current(LAW, 12)
DIAGNOSTIC_12 = transport_diagnostic(CURRENT_12)
require(CURRENT_12[(1,) * 12]
        == Fraction(872326118251778482488, 625),
        "degree-twelve singleton coefficient drift")
require(not DIAGNOSTIC_12[0]
        and DIAGNOSTIC_12[1]
        == Fraction(2299434408100109034608182728, 625)
        and DIAGNOSTIC_12[2]
        == Fraction(2299561656443457484826636016, 625)
        and DIAGNOSTIC_12[3]
        == ((9, 1, 1, 1),
            Fraction(127248343348450218453288, 625)),
        "degree-twelve failure diagnostic drift")


# ---------------------------------------------------------------------------
# 2. Exact death of the entire depth-at-most-four bank at degree twelve.
# ---------------------------------------------------------------------------

F8 = complement_upset(8, {(1,) * 8})
F10 = complement_upset(10, {(1,) * 10})
F11_THREE = complement_upset(11, {
    (1,) * 11,
    (2,) + (1,) * 9,
    (2, 2) + (1,) * 7,
})
F11_ONE = complement_upset(11, {(1,) * 11})
F12_ONE = complement_upset(12, {(1,) * 12})

ROWS = (
    response_row(8, F8),
    response_row(10, F10),
    response_row(11, F11_THREE),
    response_row(11, F11_ONE),
    response_row(12, F12_ONE),
)
COEFFICIENTS = (
    2154825177301047232015910202093994919810226072626362713968741,
    23533987357602417932067004154062702538970243225880783736334,
    4845881049964897990450105824662373913248601773502462674,
    517214042158781933902782482597749952061536320507147656820,
    5015113163945262699318503892659649582247080863625246273,
)
require(all(value > 0 for value in COEFFICIENTS)
        and reduce(gcd, COEFFICIENTS) == 1,
        "degree-twelve Farkas coefficients lost positivity/primitivity")
COMBINED = tuple(
    sum(coefficient * row[index]
        for coefficient, row in zip(COEFFICIENTS, ROWS))
    for index in range(len(STATES))
)
COMBINED_RANGE = (min(COMBINED), max(COMBINED))
MIN_STATES = tuple(
    state for state, value in zip(STATES, COMBINED)
    if value == COMBINED_RANGE[0]
)
MAX_STATES = tuple(
    state for state, value in zip(STATES, COMBINED)
    if value == COMBINED_RANGE[1]
)
require(all(value < 0 for value in COMBINED),
        "degree-twelve depth-four separator lost strict negativity")
require(COMBINED_RANGE == (
    -1026703906405262770537884068784375228334079877724481598818095786550552069856,
    -1281646212982665289100583072758234749313467771458224387568535766509728384,
), "degree-twelve depth-four separator range drift")
require(MIN_STATES == ((4, 5, 5, 7),)
        and MAX_STATES == (
            (1,), (3,), (1, 1, 1, 1), (1, 1, 1, 2), (2, 2, 2, 3)),
        "degree-twelve separator equality/extreme-state set drift")


print("THM-3155 depth-four selector resurrection/death barcode")
print("poles=" + repr(POLES))
print("state_counts_depth1_depth2_depth3_depth4_total="
      + repr(tuple(map(len, BY_DEPTH)) + (len(STATES),)))
print("law=" + repr(tuple(sorted(LAW.items()))))
print("strict_nontrivial_upset_counts_N5_N11=" + repr(tuple(NONZERO_COUNTS)))
print("strict_nontrivial_upset_total=" + repr(sum(NONZERO_COUNTS)))
print("minimum_masses_N5_N11=" + repr(tuple(MINIMA)))
print("minimum_upset_tags_N5_N11=" + repr(tuple(MINIMUM_TAGS)))
print("flow_equals_demand_N5_N11=" + repr(tuple(FLOW_DEMAND)))
print("law_N12_singleton=" + repr(CURRENT_12[(1,) * 12]))
print("law_N12_flow_demand_witness="
      + repr((DIAGNOSTIC_12[1], DIAGNOSTIC_12[2], DIAGNOSTIC_12[3])))
print("depth4_N12_Farkas_coefficients=" + repr(COEFFICIENTS))
print("depth4_N12_Farkas_range=" + repr(COMBINED_RANGE))
print("depth4_N12_Farkas_max_states=" + repr(MAX_STATES))
print("depth4_N12_Farkas_min_states=" + repr(MIN_STATES))
print("barcode=(C11_depth3_empty,C11_depth4_nonempty,C12_depth4_empty)")
print("scope=fixed_support_bank_averaged_virtual_prefix_currents")
print("all_exact_checks=PASS")
