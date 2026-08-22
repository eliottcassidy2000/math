#!/usr/bin/env python3
"""Exact depth-five selector resurrection/death barcode.

For support ``(1,3)``, bank ``I2``, this companion proves both sides of the
sharp depth-five cell:

* a seven-state probability law is strictly positive on every nontrivial
  coarsening upset in degrees 5 through 12;
* a positive combination of nine upset rows is strictly negative on every
  one of the 682 physical prefix states of depth at most five once degree 13
  is included.

The 403,539 strict inequalities are scanned as antichains using integer
bitmasks, so the large degree-twelve upset bank is never materialized.  Only
integer and ``Fraction`` arithmetic is used.  Floating LPs were useful for
discovery but are not part of the promoted evidence.
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


UP = load_upstream_prefix(13)
partitions = UP["partitions"]
coarsens = UP["coarsens"]
coefficient_vectors = UP["coefficient_vectors"]
reduced_poles = UP["reduced_poles"]
transport_diagnostic = UP["transport_diagnostic"]
BANK = UP["BANKS"][1]


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
    for depth in range(1, 6)
)
STATES = tuple(state for layer in BY_DEPTH for state in layer)
require(POLES == (8, 7, 6, 5, 5, 4, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1),
        "support-(1,3) pole multiset drift")
require(tuple(map(len, BY_DEPTH)) == (8, 33, 93, 200, 348)
        and len(STATES) == 682,
        "depth-five physical-state census drift")

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


def generated_upset(degree, generators):
    return frozenset(
        shape for shape in partitions(degree)
        if any(coarsens(generator, shape) for generator in generators)
    )


def response_row(degree, upset):
    row = tuple(
        sum(VECTORS[state][degree][shape] for shape in upset)
        for state in STATES
    )
    require(all(value.denominator == 1 for value in row),
            "response row lost integrality")
    return tuple(value.numerator for value in row)


def scan_all_upsets(degree, mixed, denominator):
    """Stream every antichain/upset and return exact strict statistics."""

    shapes = partitions(degree)
    size = len(shapes)
    index = {shape: position for position, shape in enumerate(shapes)}
    weights = []
    for shape in shapes:
        scaled = mixed[shape] * denominator
        require(scaled.denominator == 1,
                f"degree-{degree} law denominator drift")
        weights.append(scaled.numerator)

    upper_masks = []
    comparable_masks = []
    for left in shapes:
        upper = 0
        comparable = 0
        for position, right in enumerate(shapes):
            if coarsens(left, right):
                upper |= 1 << position
            if coarsens(left, right) or coarsens(right, left):
                comparable |= 1 << position
        upper_masks.append(upper)
        comparable_masks.append(comparable)

    full_mask = (1 << size) - 1
    top_mask = 1 << index[(degree,)]
    singleton_complement = full_mask ^ (1 << index[(1,) * degree])
    count = 0
    nontrivial = 0
    minimum = None
    witnesses = []

    def bit_sum(mask):
        answer = 0
        while mask:
            bit = mask & -mask
            answer += weights[bit.bit_length() - 1]
            mask ^= bit
        return answer

    def visit(start, chosen_mask, upset_mask, mass):
        nonlocal count, nontrivial, minimum, witnesses
        count += 1
        if upset_mask in (0, full_mask):
            require(mass == 0,
                    f"degree-{degree} trivial upset lost zero mass")
        else:
            nontrivial += 1
            require(mass > 0,
                    f"degree-{degree} acquired a nonpositive upset")
            if minimum is None or mass < minimum:
                minimum = mass
                witnesses = [upset_mask]
            elif mass == minimum:
                witnesses.append(upset_mask)
        for position in range(start, size):
            if chosen_mask & comparable_masks[position]:
                continue
            new_bits = upper_masks[position] & ~upset_mask
            visit(position + 1,
                  chosen_mask | (1 << position),
                  upset_mask | upper_masks[position],
                  mass + bit_sum(new_bits))

    visit(0, 0, 0, 0)
    expected_witness = top_mask if degree <= 6 else singleton_complement
    require(witnesses == [expected_witness],
            f"degree-{degree} minimum upset lost uniqueness/type")
    return count, nontrivial, minimum


# ---------------------------------------------------------------------------
# 1. Strict resurrection through degree twelve.
# ---------------------------------------------------------------------------

DENOMINATOR = 10**6
LAW_NUMERATORS = {
    (1,): 97856,
    (1, 1, 1): 56951,
    (1, 1, 1, 1): 140643,
    (1, 1, 2, 2, 2): 194498,
    (1, 2, 2, 3, 3): 7398,
    (2, 2, 2, 3, 3): 118572,
    (5, 5, 6, 7, 8): 384082,
}
require(sum(LAW_NUMERATORS.values()) == DENOMINATOR
        and reduce(gcd, LAW_NUMERATORS.values()) == 1
        and all(value > 0 for value in LAW_NUMERATORS.values()),
        "depth-five law left the primitive probability simplex")
LAW = {
    state: Fraction(numerator, DENOMINATOR)
    for state, numerator in LAW_NUMERATORS.items()
}

EXPECTED_TOTAL_COUNTS = (10, 27, 47, 168, 573, 3588, 19542, 379600)
EXPECTED_NONTRIVIAL_COUNTS = (8, 25, 45, 166, 571, 3586, 19540, 379598)
EXPECTED_MINIMUM_NUMERATORS = (
    53951073740640,
    14713182072397560,
    447905708691282144,
    175399501380224640,
    2519721209269962000,
    4663094943523464240,
    303959323177696749696,
    29594579030356898846112,
)

TOTAL_COUNTS = []
NONTRIVIAL_COUNTS = []
MINIMUM_NUMERATORS = []
for degree in range(5, 13):
    mixed = current(LAW, degree)
    diagnostic = transport_diagnostic(mixed)
    require(diagnostic[0] and diagnostic[1] == diagnostic[2],
            f"depth-five law failed exact transport at degree {degree}")
    total, nontrivial, minimum = scan_all_upsets(
        degree, mixed, DENOMINATOR)
    TOTAL_COUNTS.append(total)
    NONTRIVIAL_COUNTS.append(nontrivial)
    MINIMUM_NUMERATORS.append(minimum)

require(tuple(TOTAL_COUNTS) == EXPECTED_TOTAL_COUNTS,
        "degree-5-through-12 total-upset census drift")
require(tuple(NONTRIVIAL_COUNTS) == EXPECTED_NONTRIVIAL_COUNTS
        and sum(NONTRIVIAL_COUNTS) == 403539,
        "degree-5-through-12 strict-upset census drift")
require(tuple(MINIMUM_NUMERATORS) == EXPECTED_MINIMUM_NUMERATORS,
        "degree-5-through-12 strict-margin drift")


# This particular law fails at degree thirteen.  The failure is not used to
# exclude other laws; the independent nine-row certificate below does that.
CURRENT_13 = current(LAW, 13)
DIAGNOSTIC_13 = transport_diagnostic(CURRENT_13)
N13_SINGLETON = Fraction(299630262562066957334637, 62500)
N13_DEFICIT = Fraction(471610761586630417771013553, 125000)
require(CURRENT_13[(1,) * 13] == N13_SINGLETON,
        "degree-thirteen singleton coefficient drift")
require(not DIAGNOSTIC_13[0]
        and DIAGNOSTIC_13[2] - DIAGNOSTIC_13[1] == N13_DEFICIT
        and DIAGNOSTIC_13[3] == ((10, 1, 1, 1), N13_DEFICIT),
        "degree-thirteen law-failure diagnostic drift")


# ---------------------------------------------------------------------------
# 2. Exact death of the entire depth-at-most-five bank at degree thirteen.
# ---------------------------------------------------------------------------

R0 = complement_upset(8, {(1,) * 8})
R1 = complement_upset(10, {(1,) * 10})
R2 = complement_upset(11, {(1,) * 11})
R3 = complement_upset(12, {(1,) * 12})
R4 = complement_upset(13, {(1,) * 13})
R5 = generated_upset(11, ((2, 2) + (1,) * 7,))
R6 = generated_upset(13, (
    (4, 2) + (1,) * 7,
    (3, 3) + (1,) * 7,
    (3, 2, 2) + (1,) * 6,
    (2, 2, 2, 2) + (1,) * 5,
))
R7 = generated_upset(13, (
    (4, 2) + (1,) * 7,
    (3, 3, 3) + (1,) * 4,
    (3, 2, 2, 2) + (1,) * 4,
    (2, 2, 2, 2, 2) + (1,) * 3,
))
R8 = generated_upset(13, (
    (4, 2) + (1,) * 7,
    (3, 3) + (1,) * 7,
    (3, 2, 2) + (1,) * 6,
    (2, 2, 2, 2, 2) + (1,) * 3,
))

ROW_SPECS = (
    (8, R0), (10, R1), (11, R2), (12, R3), (13, R4),
    (11, R5), (13, R6), (13, R7), (13, R8),
)
ROWS = tuple(response_row(degree, upset) for degree, upset in ROW_SPECS)
COEFFICIENTS = (
    79966346203432495238210892467836051822404864748631537035246352913301789115999332193913950905464025191074830103726371948258801337873186194518,
    4385616886032821959435837177631104375265963397917293088514990249982327757237457317872391016149793462027134777492436937280742150390259271553,
    323165012303885417971182559526221438425147324616172737787579961666908713240511160571617207173787027117769558974661863206105407364575577287,
    7744440757659416200591308784108711415468181697333643558792830072051161414809397905536549772356510263760667669854566940871988150246747864,
    64198555842605394277362557914437771970209622650006876621643073243158046817590826310245420481922857881786832819412787764831035747707567,
    1238919946181093630207021295990664443370069758642989075525265110060651909301010930905141007031425295949856900214021045838758589520118799,
    57173834483208673247111573734889531022707614332675275902916730146155429982624862452119240877347978257760901011795049436871545379049,
    3220669545142804044646293595442679616378197077725758366899484010296430542402560592047500094580557804747642540554956255672701117951,
    43322631946675983735399838865716397620757342948083795752004906378365718688259231686082777317064579418881565313974450622142635728100,
)
require(all(value > 0 for value in COEFFICIENTS)
        and reduce(gcd, COEFFICIENTS) == 1,
        "degree-thirteen Farkas coefficients lost positivity/primitivity")
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
        "degree-thirteen depth-five separator lost strict negativity")
require(COMBINED_RANGE == (
    -131646038360152726651793725918294945734506800899964245935473472676004061252003844475399955296629305866663476692762840778970142279919601748623581179536276608,
    -7905909118176927518212656514112818602784521715342679413460937962315990351772853589181871947048820698679626976418226870569912439470329653112885762708992,
), "degree-thirteen depth-five separator range drift")
require(MIN_STATES == ((8,),),
        "degree-thirteen separator minimum-state drift")
require(MAX_STATES == (
    (1,), (1, 1, 2), (1, 1, 1, 1), (1, 1, 2, 2, 2),
    (1, 2, 2, 3, 3), (2, 2, 2, 3, 3), (4, 5, 5, 6, 7),
    (4, 5, 6, 7, 8), (5, 5, 6, 7, 8),
), "degree-thirteen separator equality-state set drift")


print("THM-3158 depth-five selector resurrection/death barcode")
print("poles=" + repr(POLES))
print("state_counts_depth1_depth2_depth3_depth4_depth5_total="
      + repr(tuple(map(len, BY_DEPTH)) + (len(STATES),)))
print("law_numerators_over_1e6=" + repr(tuple(sorted(LAW_NUMERATORS.items()))))
print("total_upset_counts_N5_N12=" + repr(tuple(TOTAL_COUNTS)))
print("strict_nontrivial_upset_counts_N5_N12="
      + repr(tuple(NONTRIVIAL_COUNTS)))
print("strict_nontrivial_upset_total=" + repr(sum(NONTRIVIAL_COUNTS)))
print("minimum_numerators_over_1e6_N5_N12="
      + repr(tuple(MINIMUM_NUMERATORS)))
print("minimum_upset_tags_N5_N12="
      + repr(("top", "top") + ("singleton-complement",) * 6))
print("law_N13_singleton=" + repr(CURRENT_13[(1,) * 13]))
print("law_N13_flow_demand_deficit_witness="
      + repr((DIAGNOSTIC_13[1], DIAGNOSTIC_13[2],
              DIAGNOSTIC_13[2] - DIAGNOSTIC_13[1], DIAGNOSTIC_13[3])))
print("depth5_N13_Farkas_coefficients=" + repr(COEFFICIENTS))
print("depth5_N13_Farkas_range=" + repr(COMBINED_RANGE))
print("depth5_N13_Farkas_max_states=" + repr(MAX_STATES))
print("depth5_N13_Farkas_min_states=" + repr(MIN_STATES))
print("barcode=(C12_depth4_empty,C12_depth5_nonempty,C13_depth5_empty)")
print("scope=fixed_support_bank_averaged_virtual_prefix_currents")
print("all_exact_checks=PASS")
