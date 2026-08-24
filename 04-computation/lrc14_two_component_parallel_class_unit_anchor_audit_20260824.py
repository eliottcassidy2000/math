#!/usr/bin/env python3
"""Exact controls for the LRC(14) rank-eleven two-component relation lane.

This companion has four deliberately separated jobs.

1.  It checks the pending Euclidean short-vector support-two atlas
    ``p^2+q^2 <= 195`` and its intersection with THM-3910's 17 residual
    certificate types.  This is an audit of an incoming, not-yet-canonical
    strengthening; it does not promote that strengthening.
2.  It verifies the rank-two, two-parallel-class quotient matroid and the
    bridge/coloop stopping mechanism for graphic cut/cycle augmentation.
3.  It exhausts a small exact universe for the Dirichlet unit-anchor lemma
    used in the accompanying reflection.
4.  It rechecks hostile rows showing that a short relation can remain wholly
    internal, including the inherited THM-3818 rank-eleven equality hostile.

The file is assertion-free so normal and optimized Python have identical
semantics.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from json import dumps
from math import gcd


Q_CANON = 91**6
B_CANON = Q_CANON**2
EUCLIDEAN_SQUARE_CAP = 195

SCALE_ONE_RESIDUALS = (
    (1, 3),
    (1, 4),
    (1, 9),
    (1, 10),
    (2, 3),
    (2, 9),
    (3, 7),
    (3, 8),
    (4, 7),
    (5, 6),
    (5, 12),
    (6, 11),
    (7, 10),
    (8, 9),
    (8, 21),
    (9, 11),
)
SCALE_TWO_RESIDUAL = (1, 9)

requirements = 0


def require(condition, message):
    global requirements
    requirements += 1
    if not condition:
        raise RuntimeError(message)


def factor_exponents(value):
    require(value >= 1, "factorization input must be positive")
    factors = {}
    n = value
    prime = 2
    while prime * prime <= n:
        while n % prime == 0:
            factors[prime] = factors.get(prime, 0) + 1
            n //= prime
        prime = 3 if prime == 2 else prime + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def primitive_decoder_ratio(x, y):
    """THM-3818's finite all-scale primitive-ratio predicate.

    The common scale is intentionally unrestricted.  MISTAKE-460 forbids
    replacing this by the stronger all-inert table-free predicate.
    """

    if x <= 0 or y <= 0 or x == y:
        return False
    common = gcd(x, y)
    p, q = sorted((x // common, y // common))
    if p + q > 356:
        return False
    factors = factor_exponents(p + q)
    return all(prime % 3 == 2 and exponent <= 2 for prime, exponent in factors.items())


def decoder_edges(row):
    return tuple(
        (i, j)
        for i in range(len(row))
        for j in range(i + 1, len(row))
        if primitive_decoder_ratio(row[i], row[j])
    )


def graph_components(vertex_count, edges):
    adjacency = [set() for _ in range(vertex_count)]
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    unseen = set(range(vertex_count))
    answer = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        stack = [root]
        component = []
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbor in sorted(adjacency[vertex], reverse=True):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        answer.append(tuple(sorted(component)))
    return tuple(sorted(answer, key=lambda item: (len(item), item)))


def rational_rank(rows):
    if not rows:
        return 0
    matrix = [[Fraction(value) for value in row] for row in rows]
    row_count = len(matrix)
    column_count = len(matrix[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (index for index in range(pivot_row, row_count) if matrix[index][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        pivot_value = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot_value for value in matrix[pivot_row]]
        for index in range(row_count):
            if index == pivot_row or not matrix[index][column]:
                continue
            multiplier = matrix[index][column]
            matrix[index] = [
                left - multiplier * right
                for left, right in zip(matrix[index], matrix[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def incidence_rows(vertex_count, edges):
    rows = []
    for left, right in edges:
        row = [0] * vertex_count
        row[left] = 1
        row[right] = -1
        rows.append(tuple(row))
    return tuple(rows)


def safe_at_one_fourteenth(row):
    return all(min(value % 14, 14 - value % 14) >= 1 for value in row)


def support_two_norm_square(x, y):
    common = gcd(x, y)
    return (x // common) ** 2 + (y // common) ** 2


def euclidean_pair_atlas():
    pairs = tuple(
        (p, q)
        for p in range(1, 14)
        for q in range(p + 1, 14)
        if gcd(p, q) == 1 and p * p + q * q <= EUCLIDEAN_SQUARE_CAP
    )
    require(len(pairs) == 47, "pending Euclidean support-two atlas count changed")

    absorbed = tuple(pair for pair in SCALE_ONE_RESIDUALS if pair in pairs)
    outliers = tuple(pair for pair in SCALE_ONE_RESIDUALS if pair not in pairs)
    require(len(absorbed) == 14, "scale-one internal-tail absorption count changed")
    require(outliers == ((8, 21), (9, 11)), "scale-one Euclidean outlier list changed")
    require(SCALE_TWO_RESIDUAL in pairs, "scale-two (1,9) must remain internally absorbable")
    require(
        tuple((pair, pair[0] ** 2 + pair[1] ** 2) for pair in outliers)
        == (((8, 21), 505), ((9, 11), 202)),
        "outlier squared norms changed",
    )
    odd_character = tuple(pair for pair in SCALE_ONE_RESIDUALS if (pair[1] - pair[0]) % 2)
    even_character = tuple(pair for pair in SCALE_ONE_RESIDUALS if not (pair[1] - pair[0]) % 2)
    require(len(odd_character) == 12, "tail odd-character absorption count changed")
    require(
        even_character == ((1, 3), (1, 9), (3, 7), (9, 11)),
        "tail even-character residual list changed",
    )
    return {
        "scale_one_even_tail_character": even_character,
        "support_two_pair_count": len(pairs),
        "scale_one_absorbed_by_tail": absorbed,
        "scale_one_outliers": tuple(
            (pair[0], pair[1], pair[0] ** 2 + pair[1] ** 2) for pair in outliers
        ),
        "scale_one_odd_tail_character_count": len(odd_character),
        "scale_two_tail_character": "even",
        "scale_two_tail_norm_square": sum(value * value for value in SCALE_TWO_RESIDUAL),
    }


def parallel_class_and_coloop_control():
    class_left = tuple(range(11))
    class_right = (11, 12)
    class_of = {vertex: 0 for vertex in class_left}
    class_of.update({vertex: 1 for vertex in class_right})

    triples = tuple(combinations(range(13), 3))
    internal_pair_witnesses = 0
    for triple in triples:
        counts = (sum(class_of[v] == 0 for v in triple), sum(class_of[v] == 1 for v in triple))
        require(max(counts) >= 2, "a two-class triple must repeat one parallel class")
        internal_pair_witnesses += 1
    require(len(triples) == 286 and internal_pair_witnesses == 286, "triple saturation count changed")

    quotient_columns = tuple((vertex + 1, 0) for vertex in class_left) + ((0, 1), (0, 3))
    require(rational_rank(quotient_columns) == 2, "quotient column rank must be two")
    for triple in combinations(quotient_columns, 3):
        require(rational_rank(triple) <= 2, "every quotient triple must be dependent")

    bases = tuple((left, right) for left in class_left for right in class_right)
    require(len(bases) == 22, "two-parallel-class basis count changed")
    for first in bases:
        for second in bases:
            for removed in set(first) - set(second):
                replacements = set(second) - set(first)
                require(
                    any(
                        len({class_of[v] for v in (set(first) - {removed}) | {replacement}}) == 2
                        for replacement in replacements
                    ),
                    "basis exchange failed inside the supplied two-class ground set",
                )

    left_path = tuple((i, i + 1) for i in range(10))
    right_edge = ((11, 12),)
    internal_edges = left_path + right_edge
    cross_edge = (10, 11)
    internal_rank = rational_rank(incidence_rows(13, internal_edges))
    augmented_rank = rational_rank(incidence_rows(13, internal_edges + (cross_edge,)))
    require(internal_rank == 11 and augmented_rank == 12, "graphic rank increment changed")
    require(
        len(graph_components(13, internal_edges)) == 2
        and len(graph_components(13, internal_edges + (cross_edge,))) == 1,
        "cross edge must be a bridge between the components",
    )
    require(
        rational_rank(incidence_rows(13, internal_edges + (cross_edge,)))
        > rational_rank(incidence_rows(13, internal_edges)),
        "the supplied cross edge must be a coloop",
    )
    return {
        "bases": len(bases),
        "internal_graph_rank": internal_rank,
        "augmented_graph_rank": augmented_rank,
        "cross_edge_role": "bridge_and_coloop_no_cycle_equation",
        "triple_internal_pair_witnesses": internal_pair_witnesses,
    }


def dirichlet_unit_relation(q_bound, maximum, tail):
    """Return a+ b*maximum - c*tail = 0 with height at most q_bound.

    This executable chooses the best exact multiple.  The reflection gives
    the all-Q pigeonhole proof.
    """

    require(1 <= tail < maximum <= q_bound * q_bound, "Dirichlet domain violated")
    best = None
    for c_value in range(1, q_bound + 1):
        quotient, _ = divmod(c_value * tail, maximum)
        for b_value in (quotient, quotient + 1):
            if not 0 <= b_value <= q_bound:
                continue
            a_value = c_value * tail - b_value * maximum
            key = (abs(a_value), c_value, b_value)
            if best is None or key < best[0]:
                best = (key, (a_value, b_value, c_value))
    require(best is not None, "Dirichlet search produced no candidate")
    a_value, b_value, c_value = best[1]
    require(
        a_value + b_value * maximum - c_value * tail == 0,
        "Dirichlet relation arithmetic failed",
    )
    require(max(abs(a_value), abs(b_value), abs(c_value)) <= q_bound, "height bound failed")
    require(c_value > 0 and (a_value != 0 or b_value != 0), "relation must cross the cut")
    return best[1]


def dirichlet_small_universe():
    cases = 0
    boundary_cases = 0
    maximum_height = 0
    for q_bound in range(2, 26):
        for maximum in range(2, q_bound * q_bound + 1):
            for tail in range(1, maximum):
                relation = dirichlet_unit_relation(q_bound, maximum, tail)
                maximum_height = max(maximum_height, max(abs(value) for value in relation))
                cases += 1
                if maximum == q_bound * q_bound:
                    boundary_cases += 1
    require(maximum_height <= 25, "small-universe height escaped its declared range")

    actual_maximum = B_CANON
    actual_tail = B_CANON - 1
    actual_relation = (-1, 1, 1)
    require(
        actual_relation[0]
        + actual_relation[1] * actual_maximum
        - actual_relation[2] * actual_tail
        == 0,
        "actual-Q boundary control failed",
    )
    require(max(abs(value) for value in actual_relation) <= Q_CANON, "actual-Q height failed")
    return {
        "actual_Q_boundary_relation": actual_relation,
        "boundary_cases": boundary_cases,
        "exhaustive_cases_Q_2_through_25": cases,
        "maximum_observed_height": maximum_height,
    }


def inherited_rank_eleven_absorption_hostile():
    """Recheck THM-3818's exact 2+11 equality hostile at its decisive layer."""

    shape_pair = (1, 3)
    shape_body = (1, 3, 4, 9, 10, 16, 19, 22, 24, 33, 40)
    body_scale = 2**42
    row = shape_pair + tuple(body_scale * value for value in shape_body)
    edges = decoder_edges(row)
    components = graph_components(13, edges)
    require(tuple(sorted(map(len, components))) == (2, 11), "inherited hostile components changed")
    require(13 - len(components) == 11, "inherited hostile decoder rank changed")
    require(all((left < 2) == (right < 2) for left, right in edges), "decoder cut acquired an edge")
    require(body_scale > 2 * Q_CANON * max(shape_pair), "dominance no-crossing gate failed")
    require(sum(row) < B_CANON, "inherited hostile left the finite counterexample box")
    require(gcd(*row) == 1, "inherited hostile must be primitive")
    require(safe_at_one_fourteenth(row), "inherited hostile lost its explicit safe time")
    require(support_two_norm_square(*shape_pair) == 10, "internal Euclidean absorber changed")
    return {
        "body_scale": body_scale,
        "component_sizes": tuple(sorted(map(len, components))),
        "decoder_edges": len(edges),
        "decoder_rank": 11,
        "internal_absorber": (3, -1),
        "internal_absorber_norm_square": 10,
        "no_crossing_height_Q_reason": "body_scale_gt_2Q_max_small_component",
        "relative_tail_scale_t": 1,
        "primitive_body_maximum_U": max(shape_body),
        "row_sum": sum(row),
        "safe_time": "1/14",
    }


def residual_outlier_absorption_controls():
    body = tuple(3**power for power in range(11))
    tail_scale = 1009
    controls = []
    for pair in ((8, 21), (9, 11)):
        row = body + tuple(tail_scale * value for value in pair)
        edges = decoder_edges(row)
        components = graph_components(13, edges)
        require(tuple(sorted(map(len, components))) == (2, 11), "outlier control components changed")
        require(all((left < 11) == (right < 11) for left, right in edges), "outlier decoder cut acquired an edge")
        require(safe_at_one_fourteenth(row), "outlier control lost safe time 1/14")
        short_support_two = []
        for left, right in combinations(range(13), 2):
            norm_square = support_two_norm_square(row[left], row[right])
            if norm_square <= EUCLIDEAN_SQUARE_CAP:
                short_support_two.append((left, right, norm_square))
                require(left < 11 and right < 11, "short support-two row must remain body-internal")
        require(len(short_support_two) == 19, "powers-of-three short internal pair count changed")
        require(min(item[2] for item in short_support_two) == 10, "body absorber norm changed")
        controls.append(
            {
                "component_sizes": tuple(sorted(map(len, components))),
                "decoder_rank": 13 - len(components),
                "pair": pair,
                "pair_norm_square": pair[0] ** 2 + pair[1] ** 2,
                "short_support_two_count": len(short_support_two),
                "short_support_two_location": "eleven_body_only",
                "safe_time": "1/14",
                "tail_scale": tail_scale,
            }
        )
    return tuple(controls)


def main():
    pair_atlas = euclidean_pair_atlas()
    parallel_class = parallel_class_and_coloop_control()
    dirichlet = dirichlet_small_universe()
    inherited_hostile = inherited_rank_eleven_absorption_hostile()
    outlier_controls = residual_outlier_absorption_controls()

    semantic = {
        "B": B_CANON,
        "Q": Q_CANON,
        "dirichlet": dirichlet,
        "euclidean_pair_atlas": pair_atlas,
        "inherited_rank_eleven_hostile": inherited_hostile,
        "outlier_controls": outlier_controls,
        "parallel_class": parallel_class,
        "scope": "finite_exact_controls_and_proved_unit_anchor_lemma_not_LRC14",
    }
    semantic_hash = sha256(
        dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    print("LRC14 TWO-COMPONENT PARALLEL-CLASS / UNIT-ANCHOR AUDIT 20260824")
    print(f"Q={Q_CANON};B={B_CANON};incoming_l2_square_cap={EUCLIDEAN_SQUARE_CAP}")
    print("euclidean_pair_atlas=" + dumps(pair_atlas, sort_keys=True, separators=(",", ":")))
    print("parallel_class=" + dumps(parallel_class, sort_keys=True, separators=(",", ":")))
    print("dirichlet=" + dumps(dirichlet, sort_keys=True, separators=(",", ":")))
    print("inherited_rank_eleven_hostile=" + dumps(inherited_hostile, sort_keys=True, separators=(",", ":")))
    print("outlier_controls=" + dumps(outlier_controls, sort_keys=True, separators=(",", ":")))
    print(f"requirements={requirements}")
    print(f"semantic_sha256={semantic_hash}")
    print("RESULT=PASS;LRC14=OPEN")


if __name__ == "__main__":
    main()
