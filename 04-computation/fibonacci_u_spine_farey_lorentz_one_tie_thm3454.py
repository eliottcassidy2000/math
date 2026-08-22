#!/usr/bin/env python3
"""Exact deterministic companion for THM-3454.

The script checks the fixed-cusp Farey fan / Berggren U-spine / Lorentz
cross-determinant identity, the six edge-separation levels of four consecutive
Fibonacci spine indices (whose rooted depths are one lower), their sharp
small-index boundaries, the adjacent-versus-
opposite S4 incidence obstruction, and finite controls for the recurrence and
Pell converse.  It also freezes two information-loss hostiles: Farey graph
distance is not determinant separation, and a scalar harmonic sum is not a
faithful Boolean-subset carrier.

All truth-bearing arithmetic is integral or rational.  Explicit exceptions,
not ``assert``, enforce every gate, so normal and optimized runs agree.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path


EXPECTED_SEMANTIC_SHA256 = "ff0262595d3dea442d874294e5166c95c47c5b7514375d8b13ab604cdbd6b4f9"
FAN_BOUND = 40
PELL_BOUND = 500
QUADRUPLE_BOUND = 18

U = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
IDENTITY3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
ROOT = (3, 4, 5)
EDGE_LABELS = ("01", "02", "03", "12", "13", "23")


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def matrix_multiply(left: tuple[tuple[int, ...], ...], right: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(sum(left[i][k] * right[k][j] for k in range(len(right))) for j in range(len(right[0])))
        for i in range(len(left))
    )


def matrix_power(matrix: tuple[tuple[int, ...], ...], exponent: int) -> tuple[tuple[int, ...], ...]:
    require(exponent >= 0, ("negative matrix exponent", exponent))
    result = IDENTITY3
    base = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = matrix_multiply(result, base)
        base = matrix_multiply(base, base)
        remaining >>= 1
    return result


def matrix_vector(matrix: tuple[tuple[int, ...], ...], vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(row[j] * vector[j] for j in range(len(vector))) for row in matrix)


def bareiss_determinant(matrix: tuple[tuple[int, ...], ...]) -> int:
    work = [list(row) for row in matrix]
    size = len(work)
    require(all(len(row) == size for row in work), "nonsquare determinant")
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next((row for row in range(pivot_index + 1, size) if work[row][pivot_index] != 0), None)
            require(swap is not None, "singular Bareiss pivot")
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for i in range(pivot_index + 1, size):
            for j in range(pivot_index + 1, size):
                numerator = work[i][j] * pivot - work[i][pivot_index] * work[pivot_index][j]
                require(numerator % previous == 0, ("Bareiss divisibility", pivot_index, i, j))
                work[i][j] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def det2(left: tuple[int, int], right: tuple[int, int]) -> int:
    return left[0] * right[1] - left[1] * right[0]


def spinor(t: int) -> tuple[int, int]:
    require(t >= 1, ("nonpositive fan index", t))
    return (t + 1, t)


def farey_vector(t: int) -> tuple[int, int]:
    require(t >= 1, ("nonpositive fan index", t))
    return (t, t + 1)


def phi(parameter: tuple[int, int]) -> tuple[int, int, int]:
    m, n = parameter
    return (m * m - n * n, 2 * m * n, m * m + n * n)


def pythagorean_point(t: int) -> tuple[int, int, int]:
    return phi(spinor(t))


def lorentz(left: tuple[int, int, int], right: tuple[int, int, int]) -> int:
    return left[2] * right[2] - left[0] * right[0] - left[1] * right[1]


def q_label(t: int) -> int:
    point = pythagorean_point(t)
    return point[0] * point[0] + 2


def fibonacci_through(index: int) -> tuple[int, ...]:
    values = [0, 1]
    while len(values) <= index:
        values.append(values[-1] + values[-2])
    return tuple(values)


FIB = fibonacci_through(40)


def edge_costs(values: tuple[int, int, int, int]) -> dict[str, int]:
    return {
        f"{i}{j}": values[j] - values[i]
        for i, j in combinations(range(4), 2)
    }


def weak_levels(costs: dict[str, int]) -> tuple[tuple[str, ...], ...]:
    levels: dict[int, list[str]] = {}
    for label, value in costs.items():
        levels.setdefault(value, []).append(label)
    return tuple(tuple(sorted(levels[value])) for value in sorted(levels))


def strict_arc_count(costs: dict[str, int]) -> int:
    return sum(costs[a] != costs[b] for a, b in combinations(EDGE_LABELS, 2))


def linear_extension_count(costs: dict[str, int]) -> int:
    count = 0
    for order in permutations(EDGE_LABELS):
        position = {label: index for index, label in enumerate(order)}
        if all(
            position[left] < position[right]
            for left in EDGE_LABELS
            for right in EDGE_LABELS
            if costs[left] < costs[right]
        ):
            count += 1
    return count


def edge(label: str) -> tuple[int, int]:
    return tuple(sorted((int(label[0]), int(label[1]))))


def permute_edge(value: tuple[int, int], permutation: tuple[int, ...]) -> tuple[int, int]:
    return tuple(sorted((permutation[value[0]], permutation[value[1]])))


def pair_orbit(seed: frozenset[tuple[int, int]]) -> frozenset[frozenset[tuple[int, int]]]:
    return frozenset(
        frozenset(permute_edge(value, permutation) for value in seed)
        for permutation in permutations(range(4))
    )


def recurrence_conditions(values: tuple[int, int, int, int]) -> tuple[bool, bool]:
    x0, x1, x2, x3 = values
    return (x2 - x1 == x0, x2 - x0 == x3 - x2)


def is_recurrence_window(values: tuple[int, int, int, int]) -> bool:
    x0, x1, x2, x3 = values
    return x2 == x0 + x1 and x3 == x1 + x2


def pell_value(x0: int, x1: int) -> int:
    return x1 * x1 - x0 * x1 - x0 * x0


# Universal fan rows, checked over a nontrivial exact control interval.
fan_rows: list[tuple[int, tuple[int, int, int], int]] = []
for t in range(1, FAN_BOUND + 1):
    expected_point = (2 * t + 1, 2 * t * (t + 1), 2 * t * (t + 1) + 1)
    point = pythagorean_point(t)
    require(point == expected_point, ("point formula", t, point, expected_point))
    require(matrix_vector(matrix_power(U, t - 1), ROOT) == point, ("U word", t))
    require(point[0] * point[0] + point[1] * point[1] == point[2] * point[2], ("null", t))
    require(det2(farey_vector(t), (1, 1)) == -1, ("fixed cusp", t))
    require(q_label(t) == 4 * t * (t + 1) + 3 == (2 * t + 1) ** 2 + 2, ("Q", t))
    fan_rows.append((t, point, q_label(t)))

pair_checks = 0
for s, t in combinations(range(1, FAN_BOUND + 1), 2):
    separation = t - s
    require(abs(det2(farey_vector(s), farey_vector(t))) == separation, ("Farey determinant", s, t))
    require(abs(det2(spinor(s), spinor(t))) == separation, ("spinor determinant", s, t))
    require(lorentz(pythagorean_point(s), pythagorean_point(t)) == 2 * separation * separation, ("Lorentz", s, t))
    # On one rooted U-ray, depth(P_t)=t-1, so its tree path length is t-s.
    require((t - 1) - (s - 1) == separation, ("tree depth", s, t))
    pair_checks += 1

for t in range(1, FAN_BOUND - 1):
    require(q_label(t + 2) - 2 * q_label(t + 1) + q_label(t) == 8, ("Q second difference", t))


# Fibonacci sampling turns the quadratic branch labels into a six-mode linear
# recurrence.  The theorem proves minimality from the six distinct Binet roots.
q_fibonacci_rows = tuple(4 * value * (value + 1) + 3 for value in FIB)
q_recurrence_coefficients = (4, -2, -6, 4, 2, -1)
for n in range(6, len(q_fibonacci_rows)):
    expected = sum(q_recurrence_coefficients[j] * q_fibonacci_rows[n - 1 - j] for j in range(6))
    require(q_fibonacci_rows[n] == expected, ("Q Fibonacci recurrence", n, q_fibonacci_rows[n], expected))
q_hankel = tuple(tuple(q_fibonacci_rows[i + j] for j in range(6)) for i in range(6))
q_hankel_determinant = bareiss_determinant(q_hankel)
require(q_hankel_determinant == 393216, ("Q Fibonacci Hankel determinant", q_hankel_determinant))


# Sharp Fibonacci edge-level boundaries.
expected_stable_levels = (("01",), ("12",), ("02", "23"), ("13",), ("03",))
fibonacci_rows: list[tuple[int, tuple[int, ...], tuple[tuple[str, ...], ...], int]] = []
for k in range(2, 21):
    values = tuple(FIB[k - 1 + i] for i in range(4))
    costs = edge_costs(values)
    levels = weak_levels(costs)
    fibonacci_rows.append((k, values, levels, strict_arc_count(costs)))
    if k == 2:
        require(levels == (("01",), ("02", "12", "23"), ("03", "13")), ("k=2 boundary", levels))
        require(len(set(values)) == 3 and strict_arc_count(costs) == 11, ("k=2 counts", values, costs))
        require(linear_extension_count(costs) == 12, ("k=2 refinements", costs))
    elif k == 3:
        require(levels == (("01", "12"), ("02", "23"), ("13",), ("03",)), ("k=3 boundary", levels))
        require(len(set(values)) == 4 and strict_arc_count(costs) == 13, ("k=3 counts", values, costs))
        require(linear_extension_count(costs) == 4, ("k=3 refinements", costs))
    else:
        require(levels == expected_stable_levels, ("stable levels", k, levels))
        require(strict_arc_count(costs) == 14, ("one missing comparison", k, costs))
        require(linear_extension_count(costs) == 2, ("two T6 refinements", k, costs))
        require(costs == {
            "01": FIB[k - 2],
            "02": FIB[k],
            "03": 2 * FIB[k],
            "12": FIB[k - 1],
            "13": FIB[k + 1],
            "23": FIB[k],
        }, ("six formulas", k, costs))


# S4 preserves whether two K4 edges intersect; the new tie is adjacent while
# THM-3339's Cassini pair is opposite.
adjacent_pair = frozenset((edge("02"), edge("23")))
opposite_pair = frozenset((edge("03"), edge("12")))
adjacent_orbit = pair_orbit(adjacent_pair)
opposite_orbit = pair_orbit(opposite_pair)
require(len(adjacent_orbit) == 12, ("adjacent orbit", len(adjacent_orbit)))
require(len(opposite_orbit) == 3, ("opposite orbit", len(opposite_orbit)))
require(adjacent_orbit.isdisjoint(opposite_orbit), "S4 incidence orbits overlap")
relabel_count = sum(
    frozenset(permute_edge(value, permutation) for value in adjacent_pair) == opposite_pair
    for permutation in permutations(range(4))
)
require(relabel_count == 0, ("forbidden relabelings", relabel_count))


# The two edge conditions are exactly the two-step additive recurrence.
quadruple_checks = 0
for values in combinations(range(1, QUADRUPLE_BOUND + 1), 4):
    conditions = recurrence_conditions(values)
    require((conditions == (True, True)) == is_recurrence_window(values), ("converse", values, conditions))
    quadruple_checks += 1

recurrence_order_checks = 0
for x0 in range(1, 50):
    for x1 in range(x0 + 1, 60):
        values = (x0, x1, x0 + x1, x0 + 2 * x1)
        has_target_order = weak_levels(edge_costs(values)) == expected_stable_levels
        require(has_target_order == (x1 < 2 * x0), ("recurrence order boundary", values))
        recurrence_order_checks += 1

long_window_checks = 0
for values in combinations(range(1, 16), 5):
    gaps = tuple(values[i + 1] - values[i] for i in range(4))
    vertex_recurrence = all(values[i + 1] == values[i] + values[i - 1] for i in range(1, 4))
    gap_transplant = gaps[1] == values[0] and all(gaps[i] == gaps[i - 1] + gaps[i - 2] for i in range(2, 4))
    require(vertex_recurrence == gap_transplant, ("long-window gap transplant", values, gaps))
    long_window_checks += 1


# Finite hostile control for the Pell/Cassini classification inherited from
# THM-3339.  The theorem file carries the universal descent proof.
fib_pairs = {(FIB[i], FIB[i + 1]) for i in range(2, len(FIB) - 1)}
pell_solutions: list[tuple[int, int, int]] = []
for x0 in range(1, PELL_BOUND + 1):
    for x1 in range(x0 + 1, PELL_BOUND + 1):
        value = pell_value(x0, x1)
        if abs(value) == 1:
            require((x0, x1) in fib_pairs, ("Pell false positive", x0, x1, value))
            pell_solutions.append((x0, x1, value))
for x0, x1 in fib_pairs:
    if x1 <= PELL_BOUND:
        require(abs(pell_value(x0, x1)) == 1, ("Fibonacci false negative", x0, x1))


# Hostiles for overclaims.
first_condition_only = (1, 2, 3, 4)
require(recurrence_conditions(first_condition_only) == (True, False), "first-condition hostile")

tie_condition_only = (2, 3, 4, 6)
require(recurrence_conditions(tie_condition_only) == (False, True), "tie-condition hostile")

nonpell_recurrence = (1, 3, 4, 7)
require(is_recurrence_window(nonpell_recurrence), "non-Pell hostile recurrence")
require(abs(pell_value(*nonpell_recurrence[:2])) != 1, "non-Pell hostile passed Pell gate")

ordered_nonpell_recurrence = (3, 4, 7, 11)
require(is_recurrence_window(ordered_nonpell_recurrence), "ordered non-Pell hostile recurrence")
require(weak_levels(edge_costs(ordered_nonpell_recurrence)) == expected_stable_levels, "ordered non-Pell hostile order")
require(abs(pell_value(*ordered_nonpell_recurrence[:2])) != 1, "ordered non-Pell hostile passed Pell gate")

true_window = (2, 3, 5, 8)
translated_false_window = (1, 2, 4, 7)
require(edge_costs(true_window) == edge_costs(translated_false_window), "translation hostile costs")
require(abs(pell_value(*true_window[:2])) == abs(pell_value(*translated_false_window[:2])) == 1, "translation hostile Pell")
require(is_recurrence_window(true_window) and not is_recurrence_window(translated_false_window), "translation hostile recurrence")

pell_order_nonrecurrence = (2, 3, 6, 10)
require(abs(pell_value(*pell_order_nonrecurrence[:2])) == 1, "Pell-order hostile Pell")
require(weak_levels(edge_costs(pell_order_nonrecurrence)) == expected_stable_levels, "Pell-order hostile order")
require(not is_recurrence_window(pell_order_nonrecurrence), "Pell-order hostile recurrence")

# rho_1 and rho_4 are not Farey adjacent (determinant 3), but each is adjacent
# to 1/1, so their Farey graph distance is exactly two, not their U distance 3.
farey_graph_hostile = (1, 4)
require(abs(det2(farey_vector(1), farey_vector(4))) == 3, "Farey hostile determinant")
require(abs(det2(farey_vector(1), (1, 1))) == abs(det2(farey_vector(4), (1, 1))) == 1, "Farey hostile cusp")

# A subset of N is faithfully carried by its labelled harmonic terms, but the
# scalar valuation is not faithful.
harmonic_left = sum((Fraction(1, n) for n in {2}), Fraction(0))
harmonic_right = sum((Fraction(1, n) for n in {3, 6}), Fraction(0))
require(harmonic_left == harmonic_right == Fraction(1, 2), "harmonic collision")

# Every Fibonacci-indexed harmonic subseries converges: F_(n+2)>=2F_n gives
# a geometric two-step majorant.  These finite rows freeze the exact premise.
for n in range(2, len(FIB) - 2):
    require(FIB[n + 2] >= 2 * FIB[n], ("two-step growth", n))
reciprocal_fibonacci_partial = sum((Fraction(1, FIB[n]) for n in range(2, 31)), Fraction(0))
require(reciprocal_fibonacci_partial < 4, reciprocal_fibonacci_partial)


repo_root = Path(__file__).resolve().parents[1]
dependency_paths = (
    Path("01-canon/theorems/THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone.md"),
    Path("01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md"),
    Path("01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md"),
)
dependency_hashes = tuple((path.as_posix(), lf_sha256(repo_root / path)) for path in dependency_paths)

semantic_payload = {
    "bounds": (FAN_BOUND, PELL_BOUND, QUADRUPLE_BOUND),
    "fan_rows": fan_rows,
    "pair_checks": pair_checks,
    "q_fibonacci_rows": q_fibonacci_rows,
    "q_recurrence_coefficients": q_recurrence_coefficients,
    "q_hankel_determinant": q_hankel_determinant,
    "fibonacci_rows": fibonacci_rows,
    "incidence_orbits": (len(adjacent_orbit), len(opposite_orbit), relabel_count),
    "quadruple_checks": quadruple_checks,
    "recurrence_order_checks": recurrence_order_checks,
    "long_window_checks": long_window_checks,
    "pell_solutions": pell_solutions,
    "hostiles": {
        "first_condition_only": first_condition_only,
        "tie_condition_only": tie_condition_only,
        "nonpell_recurrence": nonpell_recurrence,
        "ordered_nonpell_recurrence": ordered_nonpell_recurrence,
        "translation_pair": (true_window, translated_false_window),
        "pell_order_nonrecurrence": pell_order_nonrecurrence,
        "farey_graph": farey_graph_hostile,
        "harmonic": (str(harmonic_left), str(harmonic_right)),
    },
    "reciprocal_fibonacci_partial": str(reciprocal_fibonacci_partial),
    "dependency_hashes": dependency_hashes,
}
semantic_sha256 = sha256(
    json.dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
).hexdigest()
if EXPECTED_SEMANTIC_SHA256:
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, (semantic_sha256, EXPECTED_SEMANTIC_SHA256))


print("THM-3454 exact companion")
print("status=FINITE-EXACT SUPPORT FOR AN ELEMENTARY UNIVERSAL PROOF")
print(f"universe=fan_t_1_to_{FAN_BOUND};all_{pair_checks}_pairs;increasing_quadruples_1_to_{QUADRUPLE_BOUND};pell_pairs_to_{PELL_BOUND}")
print("fan_map=t -> rho_t=t/(t+1) -> P_t=U^(t-1)(3,4,5) -> q_t=4t(t+1)+3=Q_(t-1)")
print("metric_identity=d_U(P_s,P_t)=abs(det(rho_s,rho_t))=sqrt(<P_s,P_t>_L/2)=abs(t-s)")
print("q_branch_second_difference=8")
print("q_Fibonacci_initial=3,11,11,27,51,123,291,731,1851,4763,12323,32043")
print("q_Fibonacci_boundary=R_0_is_polynomial_extension_at_degenerate_depth_-1;geometric_sequence_starts_n=1")
print("q_Fibonacci_minimal_recurrence=R_n=4R_(n-1)-2R_(n-2)-6R_(n-3)+4R_(n-4)+2R_(n-5)-R_(n-6)")
print("q_Fibonacci_characteristic=x^6-4x^5+2x^4+6x^3-4x^2-2x+1")
print(f"q_Fibonacci_Hankel6_determinant={q_hankel_determinant}")
print("farey_graph_warning=determinant_separation_is_not_Farey_graph_distance")
print("k2_levels=01 < {02,12,23} < {03,13};node_duplicate;strict_arcs=11")
print("k3_levels={01,12} < {02,23} < 13 < 03;strict_arcs=13")
print("k_ge_4_levels=01 < 12 < {02,23} < 13 < 03;strict_arcs=14")
print("k_ge_4_costs=(F_(k-2),F_(k-1),F_k,F_k,F_(k+1),2F_k)")
print("strict_comparison=transitive_orientation_of_K6_minus_edge_{02,23}")
print("weak_comparison=transitive_semicomplete_digraph_with_one_bidirected_equality")
print("T6_refinements=exactly_two;orient_the_tied_pair_either_way")
print(f"S4_pair_orbits=adjacent:{len(adjacent_orbit)};opposite:{len(opposite_orbit)};tie_to_Cassini_relabelings:{relabel_count}/24")
print("recurrence_iff=(d12=x0 and d02=d23)")
print("long_window_iff=g1=x0_and_the_consecutive_gap_sequence_is_Fibonacci_recurrent")
print("Fibonacci_iff=recurrence_and_abs(x1^2-x0*x1-x0^2)=1")
print("stable_recurrence_order_iff=x1<2*x0;boundary_x1=2*x0_adds_tie_{01,12}")
print(f"hostile_first_condition_only={first_condition_only}")
print(f"hostile_tie_condition_only={tie_condition_only}")
print(f"hostile_recurrence_nonPell={nonpell_recurrence}")
print(f"hostile_stable_order_recurrence_nonPell={ordered_nonpell_recurrence}")
print(f"hostile_same_costs_and_absPell_true_false={true_window}|{translated_false_window}")
print(f"hostile_absPell_and_stable_order_nonrecurrence={pell_order_nonrecurrence}")
print("hostile_Farey_graph=rho_1,rho_4 have determinant/U-distance 3 but Farey graph distance 2")
print("Boolean_carrier=A subset N maps injectively to labelled U-spine and harmonic-term subsets")
print("harmonic_scalar_hostile={2} and {3,6} both sum to 1/2")
print("Fibonacci_harmonic_boundary=every Fibonacci-indexed subseries converges by two-step geometric majorant")
print(f"quadruple_checks={quadruple_checks}")
print(f"recurrence_order_checks={recurrence_order_checks}")
print(f"long_window_checks={long_window_checks}")
print(f"pell_solutions_to_{PELL_BOUND}={pell_solutions}")
print(f"reciprocal_Fibonacci_partial_n2_to_n30={reciprocal_fibonacci_partial}")
for path, value in dependency_hashes:
    print(f"dependency_sha256[{path}]={value}")
print(f"semantic_sha256={semantic_sha256}")
print("scope=no_Farey_graph_isometry;no_product_T6_identification;no_full_tree_equivalence;no_LRC_or_JC_consequence")
