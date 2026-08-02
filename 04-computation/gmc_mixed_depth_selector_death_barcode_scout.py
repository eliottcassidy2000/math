#!/usr/bin/env python3
"""Exact mixed-depth pole-selector persistence/death-barcode scout.

At support (1,3), bank I2, fixed depth-one and fixed depth-two probability
selectors are both excluded at degree ten.  This script allows one law on the
union of the physical one- and two-pole prefix states.  It proves:

* a three-state mixed-depth law is strictly Hasse-positive for N=5,...,10;
* that law fails at N=11;
* the N=11 slice itself is nonempty, via an exact four-state law; but
* no single depth-{1,2} law persists through every N=5,...,11, by a four-row
  exact Farkas certificate.

All promoted calculations use ``Fraction``.  Floating-point LP is not used.
"""

import ast
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
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
VALUES = tuple(range(1, 9))
DEPTH_ONE = tuple((value,) for value in VALUES)
DEPTH_TWO = tuple(
    pair for pair in combinations_with_replacement(VALUES, 2)
    if pair[0] != pair[1] or MULTIPLICITY[pair[0]] >= 2
)
STATES = DEPTH_ONE + DEPTH_TWO
require(len(DEPTH_ONE) == 8 and len(DEPTH_TWO) == 33
        and len(STATES) == 41, "mixed-depth physical-state census drift")

VECTORS = {
    state: coefficient_vectors(1, BANK, 1, 3, state)
    for state in STATES
}


def response_row(degree, upset):
    return tuple(
        sum(VECTORS[state][degree][shape] for shape in upset)
        for state in STATES
    )


def current(law, degree):
    return {
        shape: sum(weight * VECTORS[state][degree][shape]
                   for state, weight in law.items())
        for shape in partitions(degree)
    }


ROWS = {}
NONZERO_COUNTS = []
SIGN_INDEFINITE_COUNTS = []
for degree in range(5, 12):
    degree_rows = tuple(
        (upset, response_row(degree, upset))
        for upset in UPSETS[degree]
    )
    ROWS[degree] = degree_rows
    nonzero = tuple((upset, row) for upset, row in degree_rows if any(row))
    NONZERO_COUNTS.append(len(nonzero))
    SIGN_INDEFINITE_COUNTS.append(sum(min(row) < 0 for _, row in nonzero))

require(tuple(NONZERO_COUNTS)
        == (8, 25, 45, 166, 571, 3586, 19540),
        "nonzero mixed-depth response-row census drift")
require(tuple(SIGN_INDEFINITE_COUNTS)
        == (0, 1, 1, 7, 55, 871, 8642),
        "sign-indefinite mixed-depth row census drift")


# ---------------------------------------------------------------------------
# 1. A strictly positive stopping-depth law through degree ten.
# ---------------------------------------------------------------------------

LAW_10 = {
    (1,): Fraction(2, 9),
    (2,): Fraction(2, 9),
    (2, 2): Fraction(5, 9),
}
require(sum(LAW_10.values()) == 1, "degree-ten law lost total mass one")

LAW_10_FLOWS = {}
positive_masses = []
for degree in range(5, 11):
    mixed = current(LAW_10, degree)
    diagnostic = transport_diagnostic(mixed)
    require(diagnostic[0], "mixed-depth degree-ten law failed transport")
    LAW_10_FLOWS[degree] = diagnostic
    for upset, row in ROWS[degree]:
        if not any(row):
            continue
        mass = sum(LAW_10[state] * row[STATES.index(state)]
                   for state in LAW_10)
        require(mass > 0,
                "mixed-depth degree-ten law lost strict upset positivity")
        positive_masses.append((mass, degree, upset))

require(len(positive_masses) == 4401,
        "strict degree-5-through-10 row census drift")
minimum_mass, minimum_degree, minimum_upset = min(
    positive_masses, key=lambda item: item[0])
require((9 * minimum_mass, minimum_degree, minimum_upset)
        == (802795680, 5, frozenset({(5,)})),
        "sharp positive mixed-depth margin drift")

LAW_10_N11 = transport_diagnostic(current(LAW_10, 11))
require(not LAW_10_N11[0]
        and LAW_10_N11[2] - LAW_10_N11[1]
        == 218510378218973888
        and LAW_10_N11[3]
        == ((8, 1, 1, 1), Fraction(218510378218973888)),
        "degree-eleven death of the degree-ten law drifted")


# ---------------------------------------------------------------------------
# 2. The N=11 slice is nevertheless nonempty.
# ---------------------------------------------------------------------------

LAW_11 = {
    (1,): Fraction(
        36522454716418213109737517674976206575815,
        122266832198150801495362183520814156562098),
    (1, 2): Fraction(
        55718439327751924214134941704508925909799,
        366800496594452404486086550562442469686294),
    (2, 2): Fraction(
        66928634320223529182681309111751364743415,
        122266832198150801495362183520814156562098),
    (7, 8): Fraction(
        728790156775253394695128497750829818805,
        366800496594452404486086550562442469686294),
}
require(sum(LAW_11.values()) == 1 and all(value > 0 for value in LAW_11.values()),
        "degree-eleven fibre law left the probability simplex")

LAW_11_CURRENT = current(LAW_11, 11)
LAW_11_DIAGNOSTIC = transport_diagnostic(LAW_11_CURRENT)
require(LAW_11_DIAGNOSTIC[0],
        "degree-eleven fibre law failed exact Hasse transport")
law_11_masses = []
for upset, row in ROWS[11]:
    if not any(row):
        continue
    mass = sum(LAW_11[state] * row[STATES.index(state)]
               for state in LAW_11)
    require(mass >= 0, "degree-eleven fibre law violated an exact upset")
    law_11_masses.append(mass)
require(sum(value == 0 for value in law_11_masses) == 3,
        "degree-eleven fibre vertex active-facet census drift")

LAW_11_ALL_DEGREES = tuple(
    transport_diagnostic(current(LAW_11, degree))[0]
    for degree in range(5, 12)
)


# ---------------------------------------------------------------------------
# 3. No one law persists through all seven degrees.
# ---------------------------------------------------------------------------

def complement_upset(degree, excluded):
    return frozenset(shape for shape in partitions(degree)
                     if shape not in excluded)


F8 = complement_upset(8, {(1,) * 8})
F10 = complement_upset(10, {(1,) * 10})
F11_THREE = complement_upset(11, {
    (1,) * 11,
    (2,) + (1,) * 9,
    (2, 2) + (1,) * 7,
})
F11_ONE = complement_upset(11, {(1,) * 11})

FARKAS_ROWS = (
    response_row(8, F8),
    response_row(10, F10),
    response_row(11, F11_THREE),
    response_row(11, F11_ONE),
)
FARKAS_COEFFICIENTS = (461295, 3948, 1, 22)
FARKAS_COMBINED = tuple(
    sum(coefficient * row[index]
        for coefficient, row in zip(FARKAS_COEFFICIENTS, FARKAS_ROWS))
    for index in range(len(STATES))
)
require(all(value < 0 for value in FARKAS_COMBINED),
        "cumulative mixed-depth Farkas row lost strict negativity")
require((min(FARKAS_COMBINED), max(FARKAS_COMBINED))
        == (-213093302419995239808, -340278874544980224),
        "cumulative mixed-depth Farkas range drift")


print("exact mixed-depth selector death-barcode scout")
print("state_counts_depth1_depth2_total="
      + repr((len(DEPTH_ONE), len(DEPTH_TWO), len(STATES))))
print("upset_counts_N5_N11="
      + repr(tuple(len(UPSETS[n]) for n in range(5, 12))))
print("nonzero_row_counts_N5_N11=" + repr(tuple(NONZERO_COUNTS)))
print("sign_indefinite_counts_N5_N11="
      + repr(tuple(SIGN_INDEFINITE_COUNTS)))
print("law_N10=" + repr(tuple(sorted(LAW_10.items()))))
print("law_N10_strict_rows=4401")
print("law_N10_minimum_9xmass_degree_upset="
      + repr((9 * minimum_mass, minimum_degree,
              tuple(sorted(minimum_upset)))))
print("law_N10_flow_demand=" + repr(tuple(
    (degree, LAW_10_FLOWS[degree][1], LAW_10_FLOWS[degree][2])
    for degree in range(5, 11))))
print("law_N10_at_N11_flow_demand_witness=" + repr((
    LAW_10_N11[1], LAW_10_N11[2], LAW_10_N11[3])))
print("law_N11=" + repr(tuple(sorted(LAW_11.items()))))
print("law_N11_zero_active_rows=3")
print("law_N11_flow_demand="
      + repr((LAW_11_DIAGNOSTIC[1], LAW_11_DIAGNOSTIC[2])))
print("law_N11_passes_by_degree_N5_N11=" + repr(LAW_11_ALL_DEGREES))
print("cumulative_Farkas_degrees_excluded=" + repr((
    (8, ((1,) * 8,)),
    (10, ((1,) * 10,)),
    (11, ((1,) * 11, (2,) + (1,) * 9, (2, 2) + (1,) * 7)),
    (11, ((1,) * 11,)),
)))
print("cumulative_Farkas_coefficients=" + repr(FARKAS_COEFFICIENTS))
print("cumulative_Farkas_combined=" + repr(FARKAS_COMBINED))
print("cumulative_Farkas_range="
      + repr((min(FARKAS_COMBINED), max(FARKAS_COMBINED))))
print("persistent_selector_N5_N11_exists=0")
print("N11_slice_selector_exists=1")
