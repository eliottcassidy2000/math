#!/usr/bin/env python3
r"""All-Q bridge from one addressed repair cell to a carrier-averaged tree.

This exact scratch theorem compares local and carrier-averaged overlap objects
for the persistent
reflected ray

    E=(1,2,4,9,11,12),  L=5544,
    q=(2Q,Q,Q,Q,Q,2Q),  z_i=q_i L-e_i.                    (1)

The local pair referee gives a pointwise lower bound ``w_ij`` for every
body-safe cell overlap.  Averaging preserves it, so ``omega_ij>=w_ij`` for

    omega_ij = |J|^{-1} sum_(j in J) mu(A_i(j) intersect A_j(j)).

The map loses the address of a repair cell but retains the Hunter forest
predicate.  The body carrier is the sidecar which recovers low-level overlap.
On (1), every local translation-free ``w_ij`` is zero, yet the actual carrier
weights close for every ``Q>=1``.

The relevant combinatorial object is a spanning tree, not merely a matching.
For any forest ``T`` and active clause set ``S``,

    |E(T[S])| <= |S|-1              when S is nonempty,

so pointwise

    1_(union A_i) <= sum_i 1_A_i-sum_(ij in T)1_(A_i intersect A_j).  (F)

Every edge belongs to 432 of the 1296 labelled spanning trees of K6.  Thus
the average tree contains one third of the total pair overlap.  This both
strictly strengthens matching subtraction and removes the old level-
multiplicity obstruction: the unequal-level graph is complete multipartite
and connected for every nonconstant level word.

The same replacement sharpens the body-independent local tail.  If
``m=min q_i`` and ``D=max q_i-m>=1``, choose a minimum vertex ``a`` and a
higher vertex ``b``; join every other minimum to ``b`` and every other higher
vertex to ``a``.  All five double-star edges have low endpoint ``m`` and gap
at most ``D``.  Choosing the two centres to minimize the signed label
correction and exhausting all 3003 bodies gives the sharp uniform label tax
``11/42``.  The pair theorem, ``G>=35/2976``, and ``epsilon<=1/(39m)``
therefore close whenever

    m > (13392/35)D + 65224/3185.                       (C0)

The joint certificate is sharper.  For every minimum/higher label partition,
give an edge ``a-b`` its exact all-gap floor
``G(1+(a-b)/L)`` and choose a maximum-floor bipartite spanning tree, breaking
ties by minimum signed label correction.  An exact 186186-row census proves
that every resulting slope is below 313 and every ``D=1`` threshold is below
297.  Since each affine threshold minus ``313D-16`` decreases for ``D>=1``,

    m >= 313D-16

closes.  Failure of this sufficient local tree certificate therefore forces
``m<=313D-17``.  This is again not a physical survivor claim.

Retaining the five individual gaps sharpens the large-spread branch further.
For ``a=11/168`` the exact endpoint minimum is

    g_*(d)=G(d-a)=(2304d-935)/[672(168d-11)].

Both ``g_*(d)`` and ``d/g_*(d)`` are increasing.  On the label-optimized
double-star, ``C<=11/42``, ``sum G>=5g_*(1)``, and
``sum d/sum G<=D/g_*(D)``.  Hence it closes above

    R(D)=3024D(168D-11)/(2304D-935)+330328/17797.

For ``D>=2``, ``R(D)-221D`` is strictly decreasing and ``R(2)<554``.
Combining this with the joint sharp ``D=1`` row gives the final clean cone

    D=1: m>=297;              D>=2: m>=221D+112.

Thus certificate failure forces ``m<=296`` at ``D=1`` and
``m<=221D+111`` for ``D>=2``.

The finite reduction uses the addressed endpoint-owner theorem.  At its good
cell ``j_*=2L/7``, the six-clause union has mass ``U_Q<1/2`` for every Q, while
the singleton sum is ``6/7+epsilon_Q``.  If ``k(x)`` clauses cover x, then

    binom(k(x),2) >= k(x)-1_(k(x)>0).

The bivariate endpoint theorem gives more than a large all-pair total.  Its
three component words admit the fixed tree

    T_*={1-9,1-11,2-9,4-11,11-12}.

On every open endpoint segment, ``T_*`` induces a tree on the active labels.
Consequently its five pair overlaps telescope exactly to singleton mass minus
union mass, hence are greater than ``5/14`` at ``j_*``.  All other cell
overlaps are nonnegative, so its carrier-averaged weight is greater than

    5/(14|J|),   |J|=2260.                                (2)

Moreover ``epsilon_Q<=epsilon_1/Q`` and

    5 < (14|J|/5) epsilon_1 < 6.

Therefore (2) beats epsilon for all ``Q>=6``.  On the finite prefix the script
computes the exact average six-clause union.  Cellwise pair overcount and K6
tree averaging give the certified tree lower bound

    (6/7+epsilon_Q-average_union_Q)/3,

which exceeds epsilon throughout ``1<=Q<=5``.  This avoids constructing huge
pair common rulers while proving the stronger existential tree conclusion.
Together these arguments prove the carrier-averaged spanning-tree certificate on
the whole hostile ray.  It does not close arbitrary reflected
``k=1`` packets: the all-Q step depends on this ray's addressed good word, and
failure of either tree certificate is never a physical survivor theorem.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, product
from math import lcm
from pathlib import Path


HERE = Path(__file__).resolve()
BASE_RELATIVE = Path("04-computation") / "lrc14_j7_critical_scalar_wall_independent_thm2941.py"
ROOT = next(parent for parent in HERE.parents if (parent / BASE_RELATIVE).is_file())
CARRIER_SOURCE = ROOT / BASE_RELATIVE
LOCAL_PAIR_SOURCE = (
    ROOT
    / "04-computation"
    / "lrc14_j7_reflected_nearest_level_matching_tail_referee_thm2941.py"
)
DILATION_SOURCE = HERE.with_name(
    "lrc14_j7_reflected_binary_dilation_resonance_referee_thm2941.py"
)
OUTPUT = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_j7_reflected_hostile_carrier_tree_allq_bridge_thm2941.out"
)
EXPECTED_CARRIER_SOURCE_SHA256 = "5d25a955fe184d6c1a3d8b632b4bbf901dc996ee46ad67c5748836fcc7134404"
EXPECTED_LOCAL_PAIR_SOURCE_SHA256 = "173b0edc01159be5ae7ec8f2c6a0d7d36bae347c67c8d1592c1f0976af6c1fb5"
EXPECTED_DILATION_SOURCE_SHA256 = "207131cb3e8d35902c626415ffa2edd3bf51e50e1665291cc45b6fddd5bf44c3"
EXPECTED_SEMANTIC_SHA256 = "08aa4a394b7eecdb567bd3896aee70dd8193da439ebfe7d0607b8069a3f31b01"
FINITE_MAX_Q = 5
TENT_FLOOR = F(35, 2976)
TREE_TARIFF = F(1283, 174096)
MAX_LABEL_DELTA = F(11, 168)
EXPECTED_SYMBOLIC_ENDPOINT_DIGEST = "f77f43a29782f354a37f07056f540f67a59babf658ba51abf24f7a591cee4e15"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CARRIER = load("hostile_carrier_engine", CARRIER_SOURCE)
DILATION = load("hostile_carrier_dilation", DILATION_SOURCE)
E = DILATION.E
TYPES = DILATION.TYPES
L = DILATION.L
GOOD_J = DILATION.GOOD_J


def canonical_edge(i: int, j: int) -> tuple[int, int]:
    return (i, j) if i < j else (j, i)


def prufer_tree(word: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    """Decode a length-four Pruefer word as a canonical labelled K6 tree."""
    degree = [1] * 6
    for vertex in word:
        degree[vertex] += 1
    edges = []
    for vertex in word:
        leaf = min(index for index, value in enumerate(degree) if value == 1)
        edges.append(canonical_edge(leaf, vertex))
        degree[leaf] -= 1
        degree[vertex] -= 1
    last = tuple(index for index, value in enumerate(degree) if value == 1)
    require(len(last) == 2, ("Pruefer decoder failed", word, degree))
    edges.append(canonical_edge(*last))
    return tuple(sorted(edges))


TREES = tuple(
    sorted(
        {
            prufer_tree(word)
            for word in product(range(6), repeat=4)
        }
    )
)

# Labels are indexed by E=(1,2,4,9,11,12).  The event orders below are exactly
# the three good components proved by the dilation referee's bivariate
# endpoint chains.  A True event opens an arc and a False event closes it.
ADDRESS_TREE = (
    (0, 3),  # 1--9
    (0, 4),  # 1--11
    (1, 3),  # 2--9
    (2, 4),  # 4--11
    (4, 5),  # 11--12
)
GOOD_COMPONENT_EVENTS = (
    (
        ("l4", 2, True),
        ("l11", 4, True),
        ("l1a", 0, True),
        ("r1a", 0, False),
        ("l12a", 5, True),
        ("r4", 2, False),
        ("r11", 4, False),
        ("r12a", 5, False),
    ),
    (
        ("l2", 1, True),
        ("l9", 3, True),
        ("l1b", 0, True),
        ("r2", 1, False),
        ("r9", 3, False),
        ("r1b", 0, False),
    ),
    (("l12c", 5, True), ("r12c", 5, False)),
)


def endpoint_active_sets():
    rows = []
    for component in GOOD_COMPONENT_EVENTS:
        active: set[int] = set()
        for _, vertex, opening in component:
            if opening:
                require(vertex not in active, ("duplicate opening", component, vertex))
                active.add(vertex)
            else:
                require(vertex in active, ("closing inactive label", component, vertex))
                active.remove(vertex)
            if active:
                rows.append(tuple(sorted(active)))
        require(not active, ("component did not close", component, active))
    return tuple(rows)


ADDRESS_ACTIVE_SETS = endpoint_active_sets()


def unequal_level_graph_connected(levels: tuple[int, ...]) -> bool:
    """Connectivity of the complete multipartite graph q_i != q_j."""
    seen = {0}
    frontier = [0]
    while frontier:
        vertex = frontier.pop()
        for other in range(6):
            if other not in seen and levels[vertex] != levels[other]:
                seen.add(other)
                frontier.append(other)
    return len(seen) == 6


def minimum_to_higher_double_star(levels: tuple[int, ...]):
    """An unequal tree whose every edge joins a minimum to a higher level."""
    minimum = min(levels)
    require(any(level > minimum for level in levels), ("constant level word", levels))
    low_center = next(index for index, level in enumerate(levels) if level == minimum)
    high_center = next(index for index, level in enumerate(levels) if level > minimum)
    edges = {canonical_edge(low_center, high_center)}
    for vertex, level in enumerate(levels):
        if vertex in (low_center, high_center):
            continue
        parent = high_center if level == minimum else low_center
        edges.add(canonical_edge(vertex, parent))
    tree = tuple(sorted(edges))
    require(len(tree) == 5 and tree in TREES, ("double-star is not a tree", levels, tree))
    require(
        all(min(levels[i], levels[j]) == minimum < max(levels[i], levels[j]) for i, j in tree),
        ("double-star edge type failed", levels, tree),
    )
    return tree


def tree_tariff_value(levels: tuple[int, ...], tree) -> F:
    total = F(0)
    for i, j in tree:
        low, high = sorted((levels[i], levels[j]))
        require(low < high, ("equal edge in tree tariff", levels, tree, (i, j)))
        total += F(high - low, low) + F(1, 14 * low)
    return total


def optimized_label_double_star(body: tuple[int, ...], mask: int):
    """Minimize the signed label tax over the two double-star centres."""
    low = tuple(index for index in range(6) if mask & (1 << index))
    high = tuple(index for index in range(6) if not (mask & (1 << index)))
    require(low and high, ("trivial label partition", body, mask))
    low_center = min(low, key=lambda index: body[index])
    high_center = max(high, key=lambda index: body[index])
    edges = {canonical_edge(low_center, high_center)}
    edges.update(canonical_edge(index, high_center) for index in low if index != low_center)
    edges.update(canonical_edge(low_center, index) for index in high if index != high_center)
    tree = tuple(sorted(edges))
    require(len(tree) == 5 and tree in TREES, ("label double-star failed", body, mask, tree))
    correction = sum(
        (body[i] - body[j] if i in low else body[j] - body[i])
        for i, j in tree
    )
    formula = (
        sum(body[index] for index in low)
        - sum(body[index] for index in high)
        + (len(high) - 1) * body[low_center]
        - (len(low) - 1) * body[high_center]
    )
    require(correction == formula, ("label-tax formula failed", body, mask, correction, formula))
    return F(correction, 14 * lcm(*body)), tree, correction, low, high


def unit_gap_edge_floor(body: tuple[int, ...], ruler: int, low: int, high: int) -> F:
    """Minimum G(d+(e_low-e_high)/L) over every integer level gap d>=1."""
    delta = F(body[low] - body[high], ruler)
    require(F(-11, 168) <= delta <= F(11, 168), ("label delta escaped", body, low, high, delta))
    if delta >= 0:
        # G(d+delta)=d/[49(d+delta)], increasing in d.
        return F(1, 49) / (1 + delta)
    # Here floor(d+delta)=d-1 and frac(d+delta)=1+delta>5/7.
    # The derivative numerator in d is nonnegative on [-11/168,0].
    derivative_numerator = (1 + delta) / 49 - (F(2, 7) + delta) ** 2 / 4
    require(derivative_numerator >= 0, ("negative d-derivative", body, low, high, delta))
    return (F(2, 7) + delta) ** 2 / (4 * (1 + delta))


def maximum_floor_tree(body: tuple[int, ...], mask: int):
    """Kruskal maximum tree, tie-broken by minimum signed label correction."""
    ruler = 14 * lcm(*body)
    low = {index for index in range(6) if mask & (1 << index)}
    require(low and len(low) < 6, ("trivial floor partition", body, mask))
    edges = []
    for i in low:
        for j in range(6):
            if j not in low:
                edges.append((unit_gap_edge_floor(body, ruler, i, j), body[i] - body[j], i, j))
    # For equal primary weights, smaller correction is preferred.  Matroid
    # greedy therefore maximizes total floor and then minimizes total tax.
    edges.sort(key=lambda row: (row[0], -row[1], -row[2], -row[3]), reverse=True)
    parent = list(range(6))

    def root(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    chosen = []
    for weight, correction, i, j in edges:
        left = root(i)
        right = root(j)
        if left != right:
            parent[left] = right
            chosen.append((i, j, weight, correction))
        if len(chosen) == 5:
            break
    require(len(chosen) == 5, ("Kruskal tree incomplete", body, mask, chosen))
    tree = tuple(canonical_edge(i, j) for i, j, _, _ in chosen)
    require(tuple(sorted(tree)) in TREES, ("Kruskal output is not a tree", body, mask, tree))
    floor_sum = sum((weight for _, _, weight, _ in chosen), F(0))
    correction_sum = sum(correction for _, _, _, correction in chosen)
    return floor_sum, F(correction_sum, ruler), tuple(chosen)


def worst_gap_floor(gap: int) -> F:
    """g_*(d)=min_{|delta|<=11/168} G(d+delta)=G(d-11/168)."""
    require(gap >= 1, ("nonpositive gap", gap))
    return F(2304 * gap - 935, 672 * (168 * gap - 11))


def gap_floor_at_delta(gap: int, delta: F) -> F:
    """Exact G(d+delta) on |delta|<=11/168."""
    require(gap >= 1 and -MAX_LABEL_DELTA <= delta <= MAX_LABEL_DELTA, (gap, delta))
    if delta >= 0:
        return F(gap, 49) / (gap + delta)
    return (F(gap - 1, 49) + (F(2, 7) + delta) ** 2 / 4) / (gap + delta)


def gap_tail_threshold(diameter: int) -> F:
    """Gap-sensitive sufficient threshold R(D) for the optimized double-star."""
    require(diameter >= 1, ("nonpositive diameter", diameter))
    return F(
        3024 * diameter * (168 * diameter - 11),
        2304 * diameter - 935,
    ) + F(330328, 17797)


def labels(scale: int) -> tuple[int, ...]:
    return tuple(a * scale * L - e for e, a in zip(E, TYPES))


def epsilon(scale: int) -> F:
    return sum(
        (F(e, 7 * z) for e, z in zip(E, labels(scale))),
        F(0),
    )


def interval_intersection(first, second):
    rows = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            rows.append((left, right))
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return tuple(rows)


def finite_profile(safe_cells: tuple[int, ...], scale: int):
    debt = epsilon(scale)
    singleton_sum = F(6, 7) + debt
    average_union = sum(
        (DILATION.direct_cell_mass(scale, j) for j in safe_cells),
        F(0),
    ) / len(safe_cells)
    tree_lower = (singleton_sum - average_union) / 3
    require(
        tree_lower > debt,
        ("finite carrier tree failed", scale, tree_lower, debt, average_union),
    )
    return scale, average_union, tree_lower, debt, tree_lower - debt


def addressed_control(safe_cells: tuple[int, ...], scale: int):
    z = labels(scale)
    clauses = {
        index: DILATION.R.merge_intervals(
            DILATION.R.direct_multiplier_arcs(L, z[index], GOOD_J)
        )
        for index in range(6)
    }
    singletons = sum(
        (DILATION.R.interval_mass(row) for row in clauses.values()),
        F(0),
    )
    union = DILATION.R.interval_mass(
        DILATION.R.merge_intervals(tuple(piece for row in clauses.values() for piece in row))
    )
    edge_weights = {
        (i, j): DILATION.R.interval_mass(interval_intersection(clauses[i], clauses[j]))
        for i, j in combinations(range(6), 2)
    }
    pair_total = sum(edge_weights.values(), F(0))
    require(singletons == F(6, 7) + epsilon(scale), ("singleton identity failed", scale))
    require(union == DILATION.closed_mass(scale) < F(1, 2), ("addressed union failed", scale, union))
    require(pair_total >= singletons - union > F(5, 14), ("pair overcount inequality failed", scale))
    fixed_tree_weight = sum((edge_weights[edge] for edge in ADDRESS_TREE), F(0))
    require(
        fixed_tree_weight == singletons - union > F(5, 14),
        ("address-tree telescoping identity failed", scale, fixed_tree_weight, singletons - union),
    )
    pointwise_best, pointwise_tree = max(
        (
            sum((edge_weights[edge] for edge in tree), F(0)),
            tree,
        )
        for tree in TREES
    )
    require(
        pointwise_best >= pair_total / 3 > F(5, 42),
        ("Cayley tree averaging failed", scale),
    )
    require(pointwise_best >= fixed_tree_weight, ("fixed tree beats maximum", scale))
    carrier_lower = fixed_tree_weight / len(safe_cells)
    require(
        carrier_lower > F(5, 14 * len(safe_cells)),
        ("address-to-carrier transfer failed", scale, carrier_lower),
    )
    return (
        scale,
        singletons,
        union,
        pair_total,
        ADDRESS_TREE,
        fixed_tree_weight,
        pointwise_tree,
        pointwise_best,
        carrier_lower,
    )


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    require(E == (1, 2, 4, 9, 11, 12), ("hostile body changed", E))
    require(TYPES == (2, 1, 1, 1, 1, 2), ("type vector changed", TYPES))
    require(L == 14 * lcm(*E) == 5544, ("hostile ruler changed", L))
    require(len(TREES) == 6**4 == 1296, ("Cayley tree count changed", len(TREES)))
    require(all(len(tree) == 5 for tree in TREES), "tree edge count changed")
    require(ADDRESS_TREE in TREES, ("address tree missing from Cayley ledger", ADDRESS_TREE))
    edge_multiplicities = tuple(
        (edge, sum(edge in tree for tree in TREES))
        for edge in combinations(range(6), 2)
    )
    require(
        all(count == 432 for _, count in edge_multiplicities),
        ("Cayley edge multiplicity changed", edge_multiplicities),
    )
    forest_subset_checks = 0
    for tree in TREES:
        for mask in range(64):
            active = {vertex for vertex in range(6) if mask & (1 << vertex)}
            induced_edges = sum(i in active and j in active for i, j in tree)
            require(
                induced_edges <= max(0, len(active) - 1),
                ("forest pointwise inequality failed", tree, mask, induced_edges),
            )
            forest_subset_checks += 1
    active_set_digest = hashlib.sha256(repr(ADDRESS_ACTIVE_SETS).encode()).hexdigest()
    require(len(ADDRESS_ACTIVE_SETS) == 13, ("good active-set count changed", ADDRESS_ACTIVE_SETS))
    for active_row in ADDRESS_ACTIVE_SETS:
        active = set(active_row)
        induced_edges = sum(i in active and j in active for i, j in ADDRESS_TREE)
        require(
            induced_edges == len(active) - 1,
            ("address tree does not span active segment", active_row, induced_edges),
        )

    expected_chain_names = (
        "ell4_to_ell11",
        "ell11_to_ell1",
        "rr1_to_ell12",
        "ell12_to_rr4",
        "rr4_to_rr11",
        "rr11_to_rr12",
        "component_A_to_B",
        "ell2_to_ell9",
        "ell9_to_ell1",
        "ell1_to_rr2",
        "rr2_to_rr9",
        "rr9_to_rr1",
        "component_B_to_C",
        "component_C_to_next",
        "first_above_zero",
        "last_below_one",
    )
    require(
        tuple(row[0] for row in DILATION.GOOD_CHAIN_FORMS) == expected_chain_names,
        "good endpoint-chain schema changed",
    )
    symbolic_endpoint_rows = DILATION.audit_symbolic_endpoint_chains()
    symbolic_endpoint_digest = hashlib.sha256(repr(symbolic_endpoint_rows).encode()).hexdigest()
    require(
        symbolic_endpoint_digest == EXPECTED_SYMBOLIC_ENDPOINT_DIGEST,
        ("symbolic endpoint digest changed", symbolic_endpoint_digest),
    )
    formal_numerator, formal_denominator, closed_limit = DILATION.audit_formal_identity()

    bodies = tuple(combinations(range(1, 15), 6))
    require(len(bodies) == 3003, ("body universe changed", len(bodies)))
    ruler_rows = sorted((14 * lcm(*body), body) for body in bodies)
    require(
        ruler_rows[0] == (168, (1, 2, 3, 4, 6, 12)),
        ("minimum ruler changed", ruler_rows[0]),
    )
    second_ruler = next(row for row in ruler_rows if row[0] > ruler_rows[0][0])
    require(second_ruler[0] == 336, ("second ruler changed", second_ruler))
    require(F(65, second_ruler[0]) < F(11, 42), "coarse second-ruler label bound failed")
    label_tax_digest = hashlib.sha256()
    joint_tree_digest = hashlib.sha256()
    label_tax_maximum = None
    label_tax_equality_rows = []
    selected_tax_maximum = None
    joint_floor_minimum = None
    joint_floor_equality_count = 0
    joint_slope_maximum = None
    joint_slope_equality_count = 0
    joint_unit_maximum = None
    joint_unit_equality_count = 0
    label_tax_rows = 0
    for body in bodies:
        body_ruler = 14 * lcm(*body)
        for mask in range(1, 63):
            ratio, tree, correction, low, high = optimized_label_double_star(body, mask)
            require(ratio <= F(11, 42), ("label tax exceeded", body, mask, ratio))
            ledger = (ratio, correction, body_ruler, body, mask, tree, low, high)
            if label_tax_maximum is None or ledger > label_tax_maximum:
                label_tax_maximum = ledger
            if ratio == F(11, 42):
                label_tax_equality_rows.append(ledger)
            label_tax_digest.update(f"{ledger}\n".encode())

            floor_sum, selected_tax, selected_edges = maximum_floor_tree(body, mask)
            slope = F(45, 2) / floor_sum
            intercept = (F(9, 2) * selected_tax + F(1, 39)) / floor_sum
            unit_threshold = slope + intercept
            require(slope < 313, ("joint slope reached 313", body, mask, slope))
            require(unit_threshold < 297, ("joint D=1 threshold reached 297", body, mask, unit_threshold))
            joint_ledger = (
                floor_sum,
                selected_tax,
                slope,
                intercept,
                unit_threshold,
                body,
                mask,
                selected_edges,
            )
            if joint_floor_minimum is None or floor_sum < joint_floor_minimum[0]:
                joint_floor_minimum = joint_ledger
                joint_floor_equality_count = 1
            elif floor_sum == joint_floor_minimum[0]:
                joint_floor_equality_count += 1
            slope_row = (slope,) + joint_ledger
            if joint_slope_maximum is None or slope > joint_slope_maximum[0]:
                joint_slope_maximum = slope_row
                joint_slope_equality_count = 1
            elif slope == joint_slope_maximum[0]:
                joint_slope_equality_count += 1
            unit_row = (unit_threshold,) + joint_ledger
            if joint_unit_maximum is None or unit_threshold > joint_unit_maximum[0]:
                joint_unit_maximum = unit_row
                joint_unit_equality_count = 1
            elif unit_threshold == joint_unit_maximum[0]:
                joint_unit_equality_count += 1
            selected_row = (selected_tax, body, mask, selected_edges)
            if selected_tax_maximum is None or selected_row > selected_tax_maximum:
                selected_tax_maximum = selected_row
            joint_tree_digest.update(f"{joint_ledger}\n".encode())
            label_tax_rows += 1
    require(label_tax_rows == 3003 * 62 == 186186, ("label-tax row count changed", label_tax_rows))
    require(
        label_tax_maximum is not None
        and label_tax_maximum[:5]
        == (F(11, 42), 44, 168, (1, 2, 3, 4, 6, 12), 32),
        ("maximum label tax changed", label_tax_maximum),
    )
    require(
        tuple(label_tax_equality_rows) == (label_tax_maximum,),
        ("label-tax equality locus changed", label_tax_equality_rows),
    )
    label_tax = label_tax_maximum[0]
    hostile_body = (1, 2, 3, 4, 6, 12)
    require(
        joint_floor_minimum is not None
        and joint_floor_minimum[:2] == (F(142918133, 1987853616), F(-11, 42))
        and joint_floor_minimum[5:7] == (hostile_body, 31)
        and joint_floor_equality_count == 1,
        ("joint floor minimum changed", joint_floor_minimum, joint_floor_equality_count),
    )
    require(
        joint_slope_maximum is not None
        and joint_slope_maximum[0] == F(44726706360, 142918133)
        and joint_slope_maximum[6:8] == (hostile_body, 31)
        and joint_slope_equality_count == 1,
        ("joint slope maximum changed", joint_slope_maximum, joint_slope_equality_count),
    )
    require(
        joint_unit_maximum is not None
        and joint_unit_maximum[0] == F(551653043364, 1857935729)
        and joint_unit_maximum[4] == F(-29794139316, 1857935729)
        and joint_unit_maximum[6:8] == (hostile_body, 31)
        and joint_unit_equality_count == 1,
        ("joint unit maximum changed", joint_unit_maximum, joint_unit_equality_count),
    )
    require(
        selected_tax_maximum is not None
        and selected_tax_maximum[:3] == (F(11, 42), hostile_body, 32),
        ("selected-tree tax maximum changed", selected_tax_maximum),
    )

    # Independent hostile-body path: brute-force every bipartite Cayley tree
    # and compare both Kruskal's primary optimum and its secondary tax tie-break.
    hostile_tree_checks = 0
    hostile_ruler = 14 * lcm(*hostile_body)
    for mask in range(1, 63):
        low = {index for index in range(6) if mask & (1 << index)}
        edge_rows = {
            canonical_edge(i, j): (
                unit_gap_edge_floor(hostile_body, hostile_ruler, i, j),
                hostile_body[i] - hostile_body[j],
            )
            for i in low
            for j in range(6)
            if j not in low
        }
        candidates = tuple(
            tree for tree in TREES if all(edge in edge_rows for edge in tree)
        )
        require(candidates, ("no hostile bipartite trees", mask))
        primary = {
            tree: sum((edge_rows[edge][0] for edge in tree), F(0))
            for tree in candidates
        }
        best_floor = max(primary.values())
        best_tax = min(
            F(sum(edge_rows[edge][1] for edge in tree), hostile_ruler)
            for tree in candidates
            if primary[tree] == best_floor
        )
        kruskal_floor, kruskal_tax, _ = maximum_floor_tree(hostile_body, mask)
        require(
            (kruskal_floor, kruskal_tax) == (best_floor, best_tax),
            ("hostile Kruskal cross-check failed", mask, kruskal_floor, best_floor, kruskal_tax, best_tax),
        )
        hostile_tree_checks += len(candidates)

    require(
        5 * TENT_FLOOR - F(1, 39) == F(9, 2) * TREE_TARIFF,
        "five-edge tariff identity changed",
    )
    coarse_cone_slope = F(9, 2) / TENT_FLOOR
    coarse_cone_intercept = (F(9, 2) * label_tax + F(1, 39)) / (5 * TENT_FLOOR)
    require(coarse_cone_slope == F(13392, 35), ("tree cone slope changed", coarse_cone_slope))
    require(
        coarse_cone_intercept == F(65224, 3185),
        ("tree cone intercept changed", coarse_cone_intercept),
    )
    require(F(383) > coarse_cone_slope, "clean tree cone slope failed")
    require(
        F(383 + 21) > coarse_cone_slope + coarse_cone_intercept,
        "clean tree cone intercept failed",
    )
    require(
        F(383 + 20) <= coarse_cone_slope + coarse_cone_intercept,
        "clean intercept is not sharp at D=1",
    )
    joint_slope_bound = joint_slope_maximum[0]
    joint_unit_bound = joint_unit_maximum[0]
    joint_intercept_witness = joint_unit_maximum[4]
    require(joint_slope_bound < 313, ("joint slope bound failed", joint_slope_bound))
    require(joint_unit_bound < 297, ("joint unit bound failed", joint_unit_bound))

    gap_constant = (
        F(9, 2) * label_tax + F(1, 39)
    ) / (5 * worst_gap_floor(1))
    require(gap_constant == F(330328, 17797), ("gap-tail constant changed", gap_constant))
    gap_digest = hashlib.sha256()
    for gap in (1, 2, 3, 7, 64, 1000):
        floor = worst_gap_floor(gap)
        require(
            floor == gap_floor_at_delta(gap, -MAX_LABEL_DELTA),
            ("worst-gap endpoint formula failed", gap, floor),
        )
        negative_derivative_left = (
            98 * gap * (-MAX_LABEL_DELTA)
            + 24 * gap
            + 49 * MAX_LABEL_DELTA**2
        )
        require(
            negative_derivative_left == F(10128 * gap + 121, 576) > 0,
            ("negative-delta monotonicity failed", gap, negative_derivative_left),
        )
        positive_endpoint_gap = (
            gap_floor_at_delta(gap, MAX_LABEL_DELTA) - floor
        )
        require(
            positive_endpoint_gap
            == F(
                11 * (9672 * gap + 935),
                672 * (168 * gap - 11) * (168 * gap + 11),
            )
            > 0,
            ("positive endpoint did not dominate", gap, positive_endpoint_gap),
        )
        for delta in (
            -MAX_LABEL_DELTA,
            -MAX_LABEL_DELTA / 2,
            F(0),
            MAX_LABEL_DELTA / 2,
            MAX_LABEL_DELTA,
        ):
            require(
                gap_floor_at_delta(gap, delta) >= floor,
                ("delta floor hostile control failed", gap, delta),
            )
        next_floor = worst_gap_floor(gap + 1)
        require(
            next_floor - floor
            == F(5489, 28 * (168 * gap - 11) * (168 * gap + 157))
            > 0,
            ("g-star monotonicity failed", gap),
        )
        ratio_forward = F(gap + 1) / next_floor - F(gap) / floor
        require(
            ratio_forward
            == F(
                672 * (387072 * gap**2 + 72912 * gap - 146795),
                (2304 * gap - 935) * (2304 * gap + 1369),
            )
            > 0,
            ("d/g-star monotonicity failed", gap, ratio_forward),
        )
        forward = (
            gap_tail_threshold(gap + 1)
            - 221 * (gap + 1)
            - (gap_tail_threshold(gap) - 221 * gap)
        )
        require(
            forward
            == -F(
                2654208 * gap**2 + 499968 * gap + 161024765,
                (2304 * gap - 935) * (2304 * gap + 1369),
            )
            < 0,
            ("R(D)-221D monotonicity failed", gap, forward),
        )
        row = (gap, floor, negative_derivative_left, positive_endpoint_gap, next_floor, ratio_forward, forward)
        gap_digest.update(f"{row}\n".encode())
    gap_D2_excess = gap_tail_threshold(2) - 442
    gap_D2_margin = 112 - gap_D2_excess
    require(
        gap_D2_excess == F(7302253542, 65368381)
        and gap_D2_margin == F(19005130, 65368381) > 0,
        ("D=2 gap-tail margin changed", gap_D2_excess, gap_D2_margin),
    )
    require(joint_unit_bound < 297, "D=1 joint branch failed")
    for diameter in (2, 3, 7, 64, 1000):
        require(
            gap_tail_threshold(diameter) < 221 * diameter + 112,
            ("large-D gap-tail cone failed", diameter, gap_tail_threshold(diameter)),
        )
    partition_digest = hashlib.sha256()
    partition_rows = 0
    for levels in product(range(1, 5), repeat=6):
        nonconstant = len(set(levels)) > 1
        require(
            unequal_level_graph_connected(levels) == nonconstant,
            ("unequal graph connectivity failed", levels),
        )
        if nonconstant:
            tree = minimum_to_higher_double_star(levels)
            minimum = min(levels)
            diameter = max(levels) - minimum
            tariff = tree_tariff_value(levels, tree)
            require(
                tariff <= F(5 * diameter, minimum) + F(5, 14 * minimum),
                ("double-star tariff bound failed", levels, tree, tariff),
            )
            partition_digest.update(f"{levels}|{tree}|{tariff}\n".encode())
        else:
            partition_digest.update(f"{levels}|constant\n".encode())
        partition_rows += 1
    positive_levels = (10000, 10000, 10001, 10001, 10002, 10002)
    positive_tree = minimum_to_higher_double_star(positive_levels)
    positive_tariff = tree_tariff_value(positive_levels, positive_tree)
    require(positive_tariff < TREE_TARIFF, ("positive tree tariff failed", positive_tariff))
    for diameter in (1, 2, 7, 64, 1000):
        coarse_minimum = 383 * diameter + 21
        require(
            5 * TENT_FLOOR
            - F(9, 2) * F(5 * diameter, coarse_minimum)
            - F(9, 2) * label_tax / coarse_minimum
            > F(1, 39 * coarse_minimum),
            ("clean coarse double-star cone failed", diameter, coarse_minimum),
        )
        require(
            joint_slope_bound * (diameter - 1) + joint_unit_bound
            < 313 * diameter - 16,
            ("clean joint tree cone failed", diameter),
        )

    direct_L, safe_ranges = DILATION.R.safe_cell_ranges(E)
    require(direct_L == L, ("safe ruler changed", direct_L))
    safe_cells = tuple(j for left, right in safe_ranges for j in range(left, right))
    require(len(safe_cells) == 2260 and GOOD_J in safe_cells, ("safe carrier changed", len(safe_cells)))
    carrier = CARRIER.carrier_for(E)
    carrier_mass = F(sum(right - left for left, right in carrier), CARRIER.RULER)
    require(carrier_mass == F(len(safe_cells), L), ("carrier/cell normalization failed", carrier_mass))

    eps1 = epsilon(1)
    threshold_product = F(14 * len(safe_cells), 5) * eps1
    require(F(5) < threshold_product < F(6), ("analytic threshold changed", threshold_product))
    require(
        all(epsilon(scale) <= eps1 / scale for scale in (1, 2, 3, 5, 6, 64)),
        "epsilon scaling control failed",
    )
    require(eps1 / 6 < F(5, 14 * len(safe_cells)), "large-Q carrier gate failed")

    finite_rows = tuple(
        finite_profile(safe_cells, scale)
        for scale in range(1, FINITE_MAX_Q + 1)
    )
    minimum = min((row[4], row[0]) for row in finite_rows)
    addressed_rows = tuple(
        addressed_control(safe_cells, scale)
        for scale in (1, 5, 6, 64)
    )

    semantic_payload = (
        E,
        TYPES,
        L,
        GOOD_J,
        len(safe_cells),
        carrier,
        carrier_mass,
        TREES,
        edge_multiplicities,
        forest_subset_checks,
        ADDRESS_TREE,
        ADDRESS_ACTIVE_SETS,
        active_set_digest,
        symbolic_endpoint_digest,
        formal_numerator,
        formal_denominator,
        closed_limit,
        len(bodies),
        ruler_rows[0],
        second_ruler,
        label_tax_rows,
        label_tax_maximum,
        tuple(label_tax_equality_rows),
        label_tax_digest.hexdigest(),
        selected_tax_maximum,
        joint_floor_minimum,
        joint_floor_equality_count,
        joint_slope_maximum,
        joint_slope_equality_count,
        joint_unit_maximum,
        joint_unit_equality_count,
        joint_tree_digest.hexdigest(),
        hostile_tree_checks,
        TREE_TARIFF,
        coarse_cone_slope,
        coarse_cone_intercept,
        joint_slope_bound,
        joint_intercept_witness,
        joint_unit_bound,
        MAX_LABEL_DELTA,
        gap_constant,
        gap_digest.hexdigest(),
        gap_D2_excess,
        gap_D2_margin,
        partition_rows,
        partition_digest.hexdigest(),
        positive_levels,
        positive_tree,
        positive_tariff,
        eps1,
        threshold_product,
        finite_rows,
        minimum,
        addressed_rows,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_CARRIER_SOURCE_SHA256 is not None:
        require(
            sha256(CARRIER_SOURCE) == EXPECTED_CARRIER_SOURCE_SHA256,
            ("carrier source changed", sha256(CARRIER_SOURCE)),
        )
    if EXPECTED_LOCAL_PAIR_SOURCE_SHA256 is not None:
        require(
            sha256(LOCAL_PAIR_SOURCE) == EXPECTED_LOCAL_PAIR_SOURCE_SHA256,
            ("local pair source changed", sha256(LOCAL_PAIR_SOURCE)),
        )
    if EXPECTED_DILATION_SOURCE_SHA256 is not None:
        require(sha256(DILATION_SOURCE) == EXPECTED_DILATION_SOURCE_SHA256, ("dilation source changed", sha256(DILATION_SOURCE)))
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic changed", semantic))

    lines = [
        "LRC14 hostile reflected ray all-Q carrier-averaged K6 spanning-tree bridge",
        f"carrier_engine_sha256={sha256(CARRIER_SOURCE)}",
        f"local_pair_referee_sha256={sha256(LOCAL_PAIR_SOURCE)}",
        f"dilation_referee_sha256={sha256(DILATION_SOURCE)}",
        f"body={E};L={L};levels=(2Q,Q,Q,Q,Q,2Q);safe_cells={len(safe_cells)};carrier_mass={qtext(carrier_mass)}",
        "source_to_target=pointwise pair overlap lower bounds w_ij average to omega_ij;carrier omega retains phase occupancy lost by translation-free minima",
        "hostile_local_behavior=all pointwise analytic w_ij weights vanish;this is NO_CERTIFICATE,not survival",
        "forest_inequality=for active S and forest T,|E(T[S])|<=|S|-1;subtract every tree edge pointwise",
        f"Cayley_ledger=1296 labelled K6 trees;5 edges each;every edge occurs 432 times;tree average=all-pair total/3;forest_subset_checks={forest_subset_checks}",
        f"tree_tariff={qtext(TREE_TARIFF)};identity=5*(35/2976)-1/39=(9/2)*tariff",
        "unequal_graph=complete multipartite;connected iff level word nonconstant;minimum-to-higher double-star supplies 5 unequal edges",
        f"label_tax=optimized_max={qtext(label_tax)};unique_witness={label_tax_maximum};rows={label_tax_rows};digest={label_tax_digest.hexdigest()};second_ruler={second_ruler}",
        f"edge_floor=minimum_over_integer_gap_d>=1_is_G(1+(e_low-e_high)/L);joint_tree_digest={joint_tree_digest.hexdigest()};hostile_bruteforce_tree_checks={hostile_tree_checks}",
        f"joint_floor_minimum={joint_floor_minimum};unique_count={joint_floor_equality_count};selected_tax_maximum={selected_tax_maximum}",
        f"joint_slope_maximum={joint_slope_maximum};unique_count={joint_slope_equality_count}",
        f"joint_D1_maximum={joint_unit_maximum};unique_count={joint_unit_equality_count}",
        f"coarse_double_star_cone=m>(13392/35)D+65224/3185;clean=m>=383D+21;uncertified_wedge=m<=383D+20;partition_checks={partition_rows};digest={partition_digest.hexdigest()}",
        "joint_tree_cone=all slopes<313 and all D=1 thresholds<297;clean m>=313D-16;uncertified_wedge=m<=313D-17",
        f"gap_floor=g*(d)=(2304d-935)/[672(168d-11)];max_label_delta={qtext(MAX_LABEL_DELTA)};gap_constant={qtext(gap_constant)};digest={gap_digest.hexdigest()}",
        f"gap_tail=R(D)=3024D(168D-11)/(2304D-935)+330328/17797;D2_excess_over_442={qtext(gap_D2_excess)};margin_to_112={qtext(gap_D2_margin)}",
        "final_tree_cone=D=1:m>=297;D>=2:m>=221D+112;uncertified_wedge=D=1:m<=296,D>=2:m<=221D+111",
        f"positive_tariff_control=levels={positive_levels};tree={positive_tree};tariff={qtext(positive_tariff)}",
        f"addressed_word=tree={ADDRESS_TREE};active_sets={ADDRESS_ACTIVE_SETS};active_set_digest={active_set_digest}",
        f"symbolic_endpoint_digest={symbolic_endpoint_digest};closed_mass_limit={qtext(closed_limit)}",
        "addressed_bridge=fixed tree spans every active segment,so its overlap telescopes exactly to singleton_sum-union>5/14",
        f"large_Q_floor=fixed carrier tree weight>5/(14*{len(safe_cells)});epsilon_Q<=epsilon_1/Q",
        f"epsilon_1={qtext(eps1)};(14|J|/5)epsilon_1={qtext(threshold_product)};analytic_tail=Q>=6",
        f"finite_prefix=1<=Q<={FINITE_MAX_Q};rows={len(finite_rows)};all_closed;minimum_margin={qtext(minimum[0])}@Q={minimum[1]}",
    ]
    for row in finite_rows:
        lines.append(
            f"FINITE;Q={row[0]};average_union={qtext(row[1])};tree_weight_lower={qtext(row[2])};epsilon={qtext(row[3])};margin={qtext(row[4])}"
        )
    for row in addressed_rows:
        lines.append(
            f"ADDRESS_CONTROL;Q={row[0]};union={qtext(row[2])};pair_total={qtext(row[3])};fixed_tree={row[4]};fixed_weight={qtext(row[5])};best_tree={row[6]};best_weight={qtext(row[7])};carrier_weight={qtext(row[8])}"
        )
    lines.extend(
        (
            "consequence=carrier-averaged spanning-tree criterion closes this hostile ray for every integer Q>=1",
            "scope_boundary=one dilation ray;not arbitrary reflected k1;criterion remains sufficient only",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
