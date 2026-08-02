#!/usr/bin/env python3
"""Exact stochastic pole-selector and barycentric-shadow scout.

The deterministic pole-prefix scout found that every first physical pole at
support (1,3), bank I2, eventually leaves the Hasse cone.  This companion asks
whether a probability law on those eight pole deletions can repair the current.

All calculations are over ``Fraction``.  The only imported research code is
the definition prefix of ``gmc_pole_prefix_hasse_current_scout.py``; its long
top-level census is deliberately not executed.
"""

import ast
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import gcd, lcm
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_pole_prefix_hasse_current_scout.py"


def load_upstream_prefix(maximum_degree=None):
    tree = ast.parse(UPSTREAM.read_text(encoding="utf-8"))
    prefix = []
    for node in tree.body:
        if (maximum_degree is not None and isinstance(node, ast.Assign)
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
    namespace = {"__file__": str(UPSTREAM)}
    module = ast.fix_missing_locations(
        ast.Module(body=prefix, type_ignores=[]))
    exec(compile(module,
                 str(UPSTREAM), "exec"), namespace)
    return namespace


UP = load_upstream_prefix()
BANKS = UP["BANKS"]
UNIVERSE = UP["UNIVERSE"]
ALL_SHAPES = UP["ALL_SHAPES"]
partitions = UP["partitions"]
coarsens = UP["coarsens"]
coefficient_vectors = UP["coefficient_vectors"]
reduced_poles = UP["reduced_poles"]
residual_roots = UP["residual_roots"]
dominant_row = UP["dominant_row"]
transport_diagnostic = UP["transport_diagnostic"]


def all_upsets(degree):
    """Enumerate every upset through its unique minimal antichain."""

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
    upsets = {
        frozenset(
            high for low in antichain for high in shapes
            if coarsens(low, high)
        )
        for antichain in antichains
    }
    require(len(upsets) == len(antichains),
            "antichain/upset correspondence lost injectivity")
    return tuple(sorted(upsets, key=lambda upset: (
        len(upset), tuple(sorted(upset)))))


UPSETS = {degree: all_upsets(degree) for degree in range(5, 10)}
require(tuple(len(UPSETS[degree]) for degree in range(5, 10))
        == (10, 27, 47, 168, 573), "upset census drift")


def primitive_row(row):
    denominator = 1
    for value in row:
        denominator = lcm(denominator, value.denominator)
    integers = tuple(int(value * denominator) for value in row)
    divisor = 0
    for value in integers:
        divisor = gcd(divisor, abs(value))
    require(divisor, "zero row has no primitive normalization")
    return tuple(value // divisor for value in integers)


def upset_rows(vectors, values):
    rows = []
    for degree in range(5, 10):
        for upset in UPSETS[degree]:
            row = tuple(sum(vectors[value][degree][shape]
                            for shape in upset) for value in values)
            if any(row):
                rows.append((degree, upset, primitive_row(row), row))
    return rows


def mix_vector(vectors, weights, degree):
    return {
        shape: sum(weight * vectors[value][degree][shape]
                   for value, weight in weights.items())
        for shape in partitions(degree)
    }


def solve_square(matrix, right):
    size = len(matrix)
    work = [[Fraction(value) for value in row] for row in matrix]
    answer = [Fraction(value) for value in right]
    rank = 0
    for column in range(size):
        pivot = next((row for row in range(rank, size)
                      if work[row][column]), None)
        if pivot is None:
            return None
        work[rank], work[pivot] = work[pivot], work[rank]
        answer[rank], answer[pivot] = answer[pivot], answer[rank]
        scale = work[rank][column]
        work[rank] = [value / scale for value in work[rank]]
        answer[rank] /= scale
        for row in range(size):
            if row == rank or not work[row][column]:
                continue
            scale = work[row][column]
            work[row] = [left - scale * right_value
                         for left, right_value
                         in zip(work[row], work[rank])]
            answer[row] -= scale * answer[rank]
        rank += 1
    return tuple(answer)


def exact_vertices(first_facet, second_facet):
    """Vertices of simplex intersected with two homogeneous halfspaces."""

    inequalities = []
    for index in range(8):
        row = [Fraction(0)] * 8
        row[index] = Fraction(1)
        inequalities.append(tuple(row))
    inequalities.extend((tuple(map(Fraction, first_facet)),
                         tuple(map(Fraction, second_facet))))

    vertices = set()
    for active in combinations(range(10), 7):
        matrix = [[Fraction(1)] * 8]
        matrix.extend(inequalities[index] for index in active)
        solution = solve_square(matrix, [Fraction(1)] + [Fraction(0)] * 7)
        if solution is None:
            continue
        if all(sum(left * right for left, right in zip(row, solution)) >= 0
               for row in inequalities):
            vertices.add(solution)
    return tuple(sorted(vertices))


def rational_monomial_values(roots, removed):
    """Fraction-valued version of the upstream virtual-alphabet recurrence."""

    power_sums = tuple(
        sum(Fraction(root) ** degree for root in roots)
        - Fraction(removed) ** degree
        for degree in range(1, 10)
    )
    values = {(): Fraction(1)}
    for shape in ALL_SHAPES:
        exponent = shape[-1]
        remainder = shape[:-1]
        value = power_sums[exponent - 1] * values[remainder]
        for old_exponent in set(remainder):
            merged = list(remainder)
            merged.remove(old_exponent)
            merged.append(old_exponent + exponent)
            merged = tuple(sorted(merged, reverse=True))
            value -= Counter(merged)[old_exponent + exponent] * values[merged]
        values[shape] = value / Counter(shape)[exponent]
    return values


def rational_removed_vectors(invariant, bank, a_value, b_value, removed):
    atoms = tuple(
        (coefficient, rational_monomial_values(
            residual_roots(invariant, row, a_value, b_value), removed))
        for coefficient, row in bank
    )
    quotient = rational_monomial_values(
        residual_roots(invariant, dominant_row(invariant), a_value, b_value),
        removed,
    )
    answer = {}
    for degree in range(5, 10):
        shapes = partitions(degree)
        phi_monomial = {
            shape: sum(coefficient * values[shape]
                       for coefficient, values in atoms)
            for shape in shapes
        }
        phi_h = sum(phi_monomial.values())
        quotient_h = sum(quotient[shape] for shape in shapes)
        vector = {
            shape: phi_h * quotient[shape]
            - phi_monomial[shape] * quotient_h
            for shape in shapes
        }
        require(sum(vector.values()) == 0,
                "rational synthetic current lost zero mass")
        answer[degree] = vector
    return answer


# ---------------------------------------------------------------------------
# 1.  The exact fibre polytope at support (1,3), bank I2.
# ---------------------------------------------------------------------------

VALUES = tuple(range(1, 9))
POLES, _ = reduced_poles(1, BANKS[1], 1, 3)
require(set(POLES) == set(VALUES), "support-(1,3) pole alphabet drift")
VECTORS = {
    value: coefficient_vectors(1, BANKS[1], 1, 3, (value,))
    for value in VALUES
}
ROWS = upset_rows(VECTORS, VALUES)

nonzero_by_degree = tuple(
    sum(degree == target for degree, _, _, _ in ROWS)
    for target in range(5, 10)
)
require(nonzero_by_degree == (8, 25, 45, 166, 571),
        "nonzero upset-row census drift")
require(len(ROWS) == 815, "total nonzero upset-row census drift")
require(len({row for _, _, row, _ in ROWS}) == 815,
        "primitive upset rows unexpectedly collided")

potential = tuple(
    (degree, upset, row) for degree, upset, row, _ in ROWS
    if min(row) < 0
)
require(len(potential) == 3, "potential facet census drift")

A = (83547971350, 15688221032, -569502791, -2926894294,
     -14620631052, -44124052840, -85206825873, -116934942806)
B = (77268972364224, 14817264695349, -301139382544,
     -1779540733230, -11604925470240, -38159896577461,
     -74842829308296, -98601655955064)
C = (-133125723152, -49818680888, 2688744487, 18251699828,
     112450828200, 398996381528, 867825875007, 1262919484552)
require({row for _, _, row in potential} == {A, B, C},
        "three primitive exceptional rows drifted")

# B is implied by A and coordinatewise nonnegativity.
ALPHA = Fraction(8315869923144, 9467425097)
B_RESIDUAL = tuple(Fraction(left) - ALPHA * right
                   for left, right in zip(B, A))
require(all(value >= 0 for value in B_RESIDUAL),
        "middle exceptional row lost its exact redundancy certificate")

VERTICES = exact_vertices(A, C)
require(len(VERTICES) == 24, "selector-polytope vertex census drift")
require(all(sum(value != 0 for value in vertex) == 2
            for vertex in VERTICES), "selector vertex lost support two")
require(not any(A[index] >= 0 and C[index] >= 0 for index in range(8)),
        "a deterministic pole unexpectedly entered the feasible polytope")
for vertex in VERTICES:
    require(all(sum(Fraction(row[index]) * vertex[index]
                    for index in range(8)) >= 0
                for _, _, row, _ in ROWS),
            "exact vertex violated an original upset inequality")

group_left = {0, 1}
group_right = set(range(2, 8))
support_pairs = Counter(
    tuple(index for index, value in enumerate(vertex) if value)
    for vertex in VERTICES
)
require(set(support_pairs) == {
    (left, right) for left in group_left for right in group_right
}, "vertex support graph drift")
require(set(support_pairs.values()) == {2},
        "each cross edge should contribute exactly two endpoints")


# The advertised strict interior law on the (r=1,r=8) edge.
LAW = {1: Fraction(3, 5), 8: Fraction(2, 5)}
LAW_VECTOR = tuple(LAW.get(value, Fraction(0)) for value in VALUES)
require(sum(Fraction(left) * right for left, right in zip(A, LAW_VECTOR)) > 0,
        "3/5 law left the A halfspace")
require(sum(Fraction(left) * right for left, right in zip(C, LAW_VECTOR)) > 0,
        "3/5 law left the C halfspace")
LAW_FLOW = {}
for degree in range(5, 10):
    diagnostic = transport_diagnostic(mix_vector(VECTORS, LAW, degree))
    require(diagnostic[0], "3/5 law failed exact Hasse transport")
    LAW_FLOW[degree] = diagnostic

EDGE_LOWER = Fraction(3077235337, 5275866162)
EDGE_UPPER = Fraction(157864935569, 174505650963)
require(EDGE_LOWER < Fraction(3, 5) < EDGE_UPPER,
        "advertised edge law left its exact interval")
require(sum(Fraction(left) * right for left, right in zip(
    A, (EDGE_LOWER,) + (Fraction(0),) * 6
       + (1 - EDGE_LOWER,))) == 0,
        "lower edge endpoint lost its A equality")
require(sum(Fraction(left) * right for left, right in zip(
    C, (EDGE_UPPER,) + (Fraction(0),) * 6
       + (1 - EDGE_UPPER,))) == 0,
        "upper edge endpoint lost its C equality")


# ---------------------------------------------------------------------------
# 2.  Same scalar mean, different Young current.
# ---------------------------------------------------------------------------

MEAN = sum(weight * value for value, weight in LAW.items())
SECOND_MOMENT = sum(weight * value**2 for value, weight in LAW.items())
VARIANCE = SECOND_MOMENT - MEAN**2
require(MEAN == Fraction(19, 5) and VARIANCE == Fraction(294, 25),
        "pole-law moments drifted")

SYNTHETIC = rational_removed_vectors(1, BANKS[1], 1, 3, MEAN)
SYNTHETIC_FLOW = {
    degree: transport_diagnostic(SYNTHETIC[degree])
    for degree in range(5, 10)
}
require(tuple(SYNTHETIC_FLOW[degree][0] for degree in range(5, 10))
        == (True, True, True, False, False),
        "synthetic barycentric-pole hostile drifted")
require(any(mix_vector(VECTORS, LAW, 5)[shape] != SYNTHETIC[5][shape]
            for shape in partitions(5)),
        "measure-valued and deterministic-mean currents collapsed")


# ---------------------------------------------------------------------------
# 3.  Chamber portability at fixed support, and a support wall.
# ---------------------------------------------------------------------------

COMMON_LAW = {1: Fraction(9, 500), 3: Fraction(487, 500),
              8: Fraction(4, 500)}
COMMON_RESULTS = {}
for invariant, bank in enumerate(BANKS):
    vectors = {
        value: coefficient_vectors(invariant, bank, 1, 3, (value,))
        for value in VALUES
    }
    for degree in range(5, 10):
        mixed = mix_vector(vectors, COMMON_LAW, degree)
        require(all(sum(mixed[shape] for shape in upset) >= 0
                    for upset in UPSETS[degree]),
                "common chamber law violated an exact upset")
        diagnostic = transport_diagnostic(mixed)
        require(diagnostic[0], "common chamber law failed exact transport")
        COMMON_RESULTS[invariant, degree] = diagnostic


# At support (1,2), bank I2, no probability law on the five physical pole
# values can satisfy even the following two Hasse inequalities simultaneously.
SMALL_VALUES = tuple(range(1, 6))
SMALL_POLES, _ = reduced_poles(1, BANKS[1], 1, 2)
require(set(SMALL_POLES) == set(SMALL_VALUES),
        "support-(1,2) pole alphabet drift")
SMALL_VECTORS = {
    value: coefficient_vectors(1, BANKS[1], 1, 2, (value,))
    for value in SMALL_VALUES
}
SMALL_ROWS = upset_rows(SMALL_VECTORS, SMALL_VALUES)
R6 = (2469992, -986920, -2955376, -3435376, -2426920)
R9 = (-60532076544, 1282401120, 16315005312,
      -34474521120, -48239160000)
require(any(degree == 6 and raw == tuple(map(Fraction, R6))
            for degree, _, _, raw in SMALL_ROWS),
        "support-wall N=6 inequality disappeared")
require(any(degree == 9 and raw == tuple(map(Fraction, R9))
            for degree, _, _, raw in SMALL_ROWS),
        "support-wall N=9 inequality disappeared")
FARKAS = tuple(10000 * left + right for left, right in zip(R6, R9))
require(all(value < 0 for value in FARKAS),
        "two-row exact Farkas certificate lost strict negativity")


# A seductive min/max rule is almost, but not quite, portable on the full
# THM-3120 support bank.  Use upset inequalities for the finite exact census
# and replay max flow on every failing case.
MINMAX_FAILURES = []
for a_value, b_value in UNIVERSE:
    for invariant, bank in enumerate(BANKS):
        poles, _ = reduced_poles(invariant, bank, a_value, b_value)
        low, high = min(poles), max(poles)
        vectors = {
            low: coefficient_vectors(
                invariant, bank, a_value, b_value, (low,)),
            high: coefficient_vectors(
                invariant, bank, a_value, b_value, (high,)),
        }
        weights = {low: Fraction(3, 5), high: Fraction(2, 5)}
        for degree in range(5, 10):
            mixed = mix_vector(vectors, weights, degree)
            upset_ok = all(sum(mixed[shape] for shape in upset) >= 0
                           for upset in UPSETS[degree])
            if not upset_ok:
                diagnostic = transport_diagnostic(mixed)
                require(not diagnostic[0],
                        "upset failure unexpectedly passed max flow")
                MINMAX_FAILURES.append(
                    (a_value, b_value, invariant, low, high, degree))

require(len(MINMAX_FAILURES) == 9,
        "full-bank min/max failure census drift")


# ---------------------------------------------------------------------------
# 4.  Degree ten is a sharp universal wall for this physical selector space.
# ---------------------------------------------------------------------------

HIGH = load_upstream_prefix(maximum_degree=10)
HIGH_VECTORS = {
    value: HIGH["coefficient_vectors"](1, HIGH["BANKS"][1],
                                       1, 3, (value,))[10]
    for value in VALUES
}
N10_SHAPES = HIGH["partitions"](10)
N10_UPSETS = all_upsets(10)
require(len(N10_SHAPES) == 42 and len(N10_UPSETS) == 3588,
        "degree-ten partition/upset census drift")
N10_ROWS = tuple(
    tuple(sum(HIGH_VECTORS[value][shape] for shape in upset)
          for value in VALUES)
    for upset in N10_UPSETS
)
require(sum(min(row) < 0 for row in N10_ROWS) == 6,
        "degree-ten sign-indefinite row census drift")

U2 = frozenset(shape for shape in N10_SHAPES if len(shape) <= 8)
U1 = frozenset(shape for shape in N10_SHAPES if len(shape) <= 9)
R10 = (-105427024770557952, -41110096821098400, 1323165981189312,
       10490179768896384, 84003185593589760, 328063977223727328,
       734662852146861120, 1034128690334272512)
S10 = (2135515427122176, 1594234834827264, -128048320760256,
       -1149608741796864, -8750331832665600, -36573564507761664,
       -89737341671630016, -137769965901950976)
require(tuple(sum(HIGH_VECTORS[value][shape] for shape in U2)
              for value in VALUES) == tuple(map(Fraction, R10)),
        "degree-ten U2 row drift")
require(tuple(sum(HIGH_VECTORS[value][shape] for shape in U1)
              for value in VALUES) == tuple(map(Fraction, S10)),
        "degree-ten U1 row drift")
N10_FARKAS = tuple(left + 11 * right for left, right in zip(R10, S10))
require(all(value < 0 for value in N10_FARKAS),
        "degree-ten selector Farkas certificate lost strict negativity")

# Adding a lazy/no-deletion state cannot repair the degree-ten wall: its two
# rank-tail coordinates are both zero, so the Farkas row forces all physical
# one-pole mass to vanish.
EMPTY_VECTOR = HIGH["coefficient_vectors"](
    1, HIGH["BANKS"][1], 1, 3, ())[10]
require(sum(EMPTY_VECTOR[shape] for shape in U2) == 0
        and sum(EMPTY_VECTOR[shape] for shape in U1) == 0,
        "empty selector acquired a degree-ten Farkas coordinate")

# Nor does jumping directly to a physical two-pole prefix help.  Respect the
# multiplicities in the reduced pole multiset; this leaves 33 unordered
# submultisets.  One upset is strictly negative on every state.
pole_multiplicity = Counter(POLES)
TWO_POLE_STATES = tuple(
    pair for pair in combinations_with_replacement(VALUES, 2)
    if pair[0] != pair[1] or pole_multiplicity[pair[0]] >= 2
)
require(len(TWO_POLE_STATES) == 33,
        "physical two-pole state census drift")
TWO_POLE_VECTORS = {
    pair: HIGH["coefficient_vectors"](
        1, HIGH["BANKS"][1], 1, 3, pair)[10]
    for pair in TWO_POLE_STATES
}
U3_EXCLUDED = {(1,) * 10, (2,) + (1,) * 8, (2, 2) + (1,) * 6}
U3 = frozenset(shape for shape in N10_SHAPES
               if shape not in U3_EXCLUDED)
TWO_POLE_ROW = tuple(
    (pair, sum(TWO_POLE_VECTORS[pair][shape] for shape in U3))
    for pair in TWO_POLE_STATES
)
require(all(value < 0 for _, value in TWO_POLE_ROW),
        "degree-ten two-pole separating row lost strict negativity")
require((min(value for _, value in TWO_POLE_ROW),
         max(value for _, value in TWO_POLE_ROW))
        == (-2197919641631883264, -64486707449419008),
        "degree-ten two-pole separating range drift")

N10_LAW_VECTOR = {
    shape: sum(weight * HIGH_VECTORS[value][shape]
               for value, weight in LAW.items())
    for shape in N10_SHAPES
}
N10_LAW_DIAGNOSTIC = HIGH["transport_diagnostic"](N10_LAW_VECTOR)
require(not N10_LAW_DIAGNOSTIC[0]
        and N10_LAW_DIAGNOSTIC[2] - N10_LAW_DIAGNOSTIC[1]
        == Fraction(269133385522535424, 5)
        and N10_LAW_DIAGNOSTIC[3]
        == ((7, 1, 1, 1), Fraction(269133385522535424, 5)),
        "degree-ten 3/5-law hostile drift")


print("exact stochastic pole-selector polytope scout")
print("upset_counts=" + repr(tuple(len(UPSETS[n]) for n in range(5, 10))))
print("nonzero_rows_by_degree=" + repr(nonzero_by_degree))
print("unique_primitive_rows=815")
print("coordinatewise_nonnegative_rows=812")
print("exceptional_A=" + repr(A))
print("exceptional_B=" + repr(B))
print("exceptional_C=" + repr(C))
print("B_redundancy_alpha=" + repr(ALPHA))
print("B_minus_alpha_A_nonnegative=" + repr(B_RESIDUAL))
print("feasible_polytope=Delta_7 intersect {A.lambda>=0,C.lambda>=0}")
print("vertex_count=24 minimal_support=2")
for number, vertex in enumerate(VERTICES, 1):
    support = tuple((index + 1, value) for index, value in enumerate(vertex)
                    if value)
    facet = "A" if sum(Fraction(left) * right
                       for left, right in zip(A, vertex)) == 0 else "C"
    print(f"vertex_{number:02d}_{facet}=" + repr(support))
print("r1_r8_feasible_interval=" + repr((EDGE_LOWER, EDGE_UPPER)))
print("law_3_over_5_flows=" + repr(tuple(
    (degree, LAW_FLOW[degree][1], LAW_FLOW[degree][2])
    for degree in range(5, 10))))
print("law_mean_second_moment_variance="
      + repr((MEAN, SECOND_MOMENT, VARIANCE)))
print("power_sum_degree2_mixture_minus_synthetic=" + repr(-VARIANCE))
print("synthetic_mean_flow_passes=" + repr(tuple(
    SYNTHETIC_FLOW[degree][0] for degree in range(5, 10))))
print("synthetic_N8_flow_demand_witness=" + repr((
    SYNTHETIC_FLOW[8][1], SYNTHETIC_FLOW[8][2],
    SYNTHETIC_FLOW[8][3])))
print("common_I1_I2_law=" + repr(tuple(sorted(COMMON_LAW.items()))))
print("common_I1_I2_mean=" + repr(sum(
    value * weight for value, weight in COMMON_LAW.items())))
print("common_I1_I2_cases=10 all_pass=1")
print("support_1_2_I2_Farkas_R6=" + repr(R6))
print("support_1_2_I2_Farkas_R9=" + repr(R9))
print("10000_R6_plus_R9=" + repr(FARKAS))
print("support_1_2_I2_probability_selector_exists=0")
print("full_bank_minmax_cases=1150 passes=1141 failures=9")
print("full_bank_minmax_failures=" + repr(tuple(MINMAX_FAILURES)))
print("N10_partition_upset_sign_indefinite_counts="
      + repr((len(N10_SHAPES), len(N10_UPSETS), 6)))
print("N10_U2_row_R=" + repr(R10))
print("N10_U1_row_S=" + repr(S10))
print("N10_R_plus_11S=" + repr(N10_FARKAS))
print("N10_probability_selector_exists=0")
print("N10_empty_state_R_S=" + repr((0, 0)))
print("N10_lazy_selector_polytope=delta_empty")
print("N10_two_pole_state_count=33")
print("N10_two_pole_U3_row=" + repr(TWO_POLE_ROW))
print("N10_two_pole_row_range=" + repr((
    min(value for _, value in TWO_POLE_ROW),
    max(value for _, value in TWO_POLE_ROW))))
print("N10_two_pole_probability_selector_exists=0")
print("law_3_over_5_N10_flow_demand_witness=" + repr((
    N10_LAW_DIAGNOSTIC[1], N10_LAW_DIAGNOSTIC[2],
    N10_LAW_DIAGNOSTIC[3])))
