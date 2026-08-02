#!/usr/bin/env python3
"""Exact depth-at-most-three pole-selector persistence wall.

This companion extends THM-3144 from the 41 physical one-/two-pole prefix
states to every multiplicity-valid prefix of length at most three.  It checks
two fixed support/bank fibres:

* at ``(a,b)=(1,3)``, bank ``I2``, the same four upset facets as THM-3144
  admit a new positive integer combination strictly negative on all 134
  states;
* at ``(a,b)=(1,2)``, bank ``I2``, two facets already separate all 42
  states.

The calculation uses integers and ``Fraction`` only.  No floating-point LP is
part of the promoted evidence.
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


UP = load_upstream_prefix(11)
partitions = UP["partitions"]
coefficient_vectors = UP["coefficient_vectors"]
reduced_poles = UP["reduced_poles"]
BANK = UP["BANKS"][1]


def physical_states(a_value, b_value, maximum_depth=3):
    poles, _ = reduced_poles(1, BANK, a_value, b_value)
    multiplicity = Counter(poles)
    values = tuple(sorted(multiplicity))
    by_depth = []
    for depth in range(1, maximum_depth + 1):
        states = tuple(
            state
            for state in combinations_with_replacement(values, depth)
            if all(Counter(state)[value] <= multiplicity[value]
                   for value in set(state))
        )
        by_depth.append(states)
    return poles, tuple(by_depth), tuple(
        state for states in by_depth for state in states)


def complement_upset(degree, excluded):
    return frozenset(shape for shape in partitions(degree)
                     if shape not in excluded)


def response_row(vectors, states, degree, upset):
    row = tuple(
        sum(vectors[state][degree][shape] for shape in upset)
        for state in states
    )
    require(all(isinstance(value, Fraction) and value.denominator == 1
                for value in row), "response row lost integrality")
    return tuple(value.numerator for value in row)


def combine_rows(coefficients, rows):
    require(len(coefficients) == len(rows)
            and all(value > 0 for value in coefficients),
            "Farkas coefficients are not strictly positive")
    require(reduce(gcd, coefficients) == 1,
            "Farkas coefficients are not primitive")
    return tuple(
        sum(coefficient * row[index]
            for coefficient, row in zip(coefficients, rows))
        for index in range(len(rows[0]))
    )


# ---------------------------------------------------------------------------
# 1. Support (1,3), bank I2: all 134 depth-at-most-three states.
# ---------------------------------------------------------------------------

POLES_13, BY_DEPTH_13, STATES_13 = physical_states(1, 3)
require(POLES_13 == (8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1),
        "support-(1,3) pole multiset drift")
require(tuple(map(len, BY_DEPTH_13)) == (8, 33, 93)
        and len(STATES_13) == 134,
        "support-(1,3) depth-three state census drift")

VECTORS_13 = {
    state: coefficient_vectors(1, BANK, 1, 3, state)
    for state in STATES_13
}

F8 = complement_upset(8, {(1,) * 8})
F10 = complement_upset(10, {(1,) * 10})
F11_THREE = complement_upset(11, {
    (1,) * 11,
    (2,) + (1,) * 9,
    (2, 2) + (1,) * 7,
})
F11_ONE = complement_upset(11, {(1,) * 11})

ROWS_13 = (
    response_row(VECTORS_13, STATES_13, 8, F8),
    response_row(VECTORS_13, STATES_13, 10, F10),
    response_row(VECTORS_13, STATES_13, 11, F11_THREE),
    response_row(VECTORS_13, STATES_13, 11, F11_ONE),
)

# The THM-3144 depth-{1,2} certificate has one and only one depth-three
# escape.  This is the first failed extension, not the promoted separator.
OLD_COEFFICIENTS = (461295, 3948, 1, 22)
OLD_COMBINED = tuple(
    sum(coefficient * row[index]
        for coefficient, row in zip(OLD_COEFFICIENTS, ROWS_13))
    for index in range(len(STATES_13))
)
OLD_NONNEGATIVE = tuple(
    (state, value) for state, value in zip(STATES_13, OLD_COMBINED)
    if value >= 0
)
require(OLD_NONNEGATIVE == (((1, 1, 2), 744310486097342352),),
        "old depth-two wall-crossing state drift")

COEFFICIENTS_13 = (
    3627302205077683560821210804482107686117129010,
    29864988664447894013528399014705164256913403,
    7514363586001882980386000831487329663540,
    314610406969924544691100890629234040301424,
)
COMBINED_13 = combine_rows(COEFFICIENTS_13, ROWS_13)
RANGE_13 = (min(COMBINED_13), max(COMBINED_13))
MAX_STATES_13 = tuple(
    state for state, value in zip(STATES_13, COMBINED_13)
    if value == RANGE_13[1]
)
MIN_STATES_13 = tuple(
    state for state, value in zip(STATES_13, COMBINED_13)
    if value == RANGE_13[0]
)
require(all(value < 0 for value in COMBINED_13),
        "depth-three support-(1,3) separator lost strict negativity")
require(RANGE_13 == (
    -674800650195958551767268010800617689469061950339284516178880,
    -1166613115072139157133415139738043609298077113840435285120,
), "depth-three support-(1,3) separator range drift")
require(MAX_STATES_13 == ((1,), (7,), (1, 1, 2), (1, 2, 2))
        and MIN_STATES_13 == ((5, 5, 6),),
        "support-(1,3) equality/extreme-state set drift")


# ---------------------------------------------------------------------------
# 2. Support (1,2), bank I2: the two-facet cross-support wall.
# ---------------------------------------------------------------------------

POLES_12, BY_DEPTH_12, STATES_12 = physical_states(1, 2)
require(POLES_12 == (5, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1),
        "support-(1,2) pole multiset drift")
require(tuple(map(len, BY_DEPTH_12)) == (5, 13, 24)
        and len(STATES_12) == 42,
        "support-(1,2) depth-three state census drift")

VECTORS_12 = {
    state: coefficient_vectors(1, BANK, 1, 2, state)
    for state in STATES_12
}
F6 = complement_upset(6, {(1,) * 6})
F9_THREE = complement_upset(9, {
    (1,) * 9,
    (2,) + (1,) * 7,
    (3,) + (1,) * 6,
})
ROWS_12 = (
    response_row(VECTORS_12, STATES_12, 6, F6),
    response_row(VECTORS_12, STATES_12, 9, F9_THREE),
)
COEFFICIENTS_12 = (14903617223, 641337)
COMBINED_12 = combine_rows(COEFFICIENTS_12, ROWS_12)
RANGE_12 = (min(COMBINED_12), max(COMBINED_12))
MAX_STATES_12 = tuple(
    state for state, value in zip(STATES_12, COMBINED_12)
    if value == RANGE_12[1]
)
MIN_STATES_12 = tuple(
    state for state, value in zip(STATES_12, COMBINED_12)
    if value == RANGE_12[0]
)
require(all(value < 0 for value in COMBINED_12),
        "depth-three support-(1,2) separator lost strict negativity")
require(RANGE_12 == (-108248637454886456, -2009645062627112),
        "depth-three support-(1,2) separator range drift")
require(MAX_STATES_12 == ((1,), (3, 4, 5))
        and MIN_STATES_12 == ((2, 3, 4),),
        "support-(1,2) equality/extreme-state set drift")


print("THM-3149 depth-three selector persistence wall")
print("support_1_3_poles=" + repr(POLES_13))
print("support_1_3_state_counts_depth1_depth2_depth3_total="
      + repr(tuple(map(len, BY_DEPTH_13)) + (len(STATES_13),)))
print("old_depth2_certificate_coefficients=" + repr(OLD_COEFFICIENTS))
print("old_certificate_nonnegative_depth3_states=" + repr(OLD_NONNEGATIVE))
print("support_1_3_depth3_certificate_coefficients="
      + repr(COEFFICIENTS_13))
print("support_1_3_combined_range=" + repr(RANGE_13))
print("support_1_3_max_states=" + repr(MAX_STATES_13))
print("support_1_3_min_states=" + repr(MIN_STATES_13))
print("support_1_3_persistent_selector_N5_N11_exists=0")
print("support_1_2_poles=" + repr(POLES_12))
print("support_1_2_state_counts_depth1_depth2_depth3_total="
      + repr(tuple(map(len, BY_DEPTH_12)) + (len(STATES_12),)))
print("support_1_2_depth3_certificate_coefficients="
      + repr(COEFFICIENTS_12))
print("support_1_2_combined_range=" + repr(RANGE_12))
print("support_1_2_max_states=" + repr(MAX_STATES_12))
print("support_1_2_min_states=" + repr(MIN_STATES_12))
print("support_1_2_persistent_selector_N5_N11_exists=0")
print("scope=fixed_bank_averaged_virtual_prefix_currents;not_stopping_process")
print("all_exact_checks=PASS")
