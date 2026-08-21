#!/usr/bin/env python3
"""Exact additive-support atlas for the exponent-two 3-by-4 sector.

This is the standard-library exact companion for THM-3603.  It works with

    A = {0, X, X+Y},  B = {0, U, U+V, U+V+W}

and primitive positive integral gaps (X,Y,U,V,W).  Output exchange is fixed by
putting the three-point support first.  Simultaneous reversal

    (X,Y;U,V,W) -> (Y,X;W,V,U)

is retained as an orientation and recorded as an involution on fibre words.

There are eighteen collision hyperplanes: each of the three positive
differences of A can equal each of the six positive differences of B.  The
script enumerates their row-space lattice over Q, takes exact positive
nullspace chambers, constructs every oriented sum-fibre word, computes
projection connectivity and fibre-cut debts, and independently checks the
result against all primitive positive gap vectors with every gap at most 16.

No floating-point arithmetic, external package, ``assert`` statement, or
optimization-sensitive gate is used.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import cmp_to_key, reduce
from hashlib import sha256
from itertools import combinations, product
from math import gcd, lcm
from pathlib import Path
from typing import Iterable, Sequence


GAP_COUNT = 5
HOSTILE_BOUND = 16

A_EDGE_LABELS = ("A01", "A02", "A12")
B_EDGE_LABELS = ("B01", "B02", "B03", "B12", "B13", "B23")
A_EDGE_PAIRS = ((0, 1), (0, 2), (1, 2))
B_EDGE_PAIRS = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))

A_DIFFERENCE_ROWS = (
    (1, 0, 0, 0, 0),       # X
    (1, 1, 0, 0, 0),       # X+Y
    (0, 1, 0, 0, 0),       # Y
)
B_DIFFERENCE_ROWS = (
    (0, 0, 1, 0, 0),       # U
    (0, 0, 1, 1, 0),       # U+V
    (0, 0, 1, 1, 1),       # U+V+W
    (0, 0, 0, 1, 0),       # V
    (0, 0, 0, 1, 1),       # V+W
    (0, 0, 0, 0, 1),       # W
)

HYPERPLANE_ROWS = tuple(
    tuple(a[index] - b[index] for index in range(GAP_COUNT))
    for a in A_DIFFERENCE_ROWS
    for b in B_DIFFERENCE_ROWS
)
HYPERPLANE_LABELS = tuple(
    f"{a_label}={b_label}"
    for a_label in A_EDGE_LABELS
    for b_label in B_EDGE_LABELS
)

EXPECTED_FLAT_COUNT = 786
EXPECTED_POSITIVE_CONNECTED_CONE_COUNT = 83
EXPECTED_GLOBAL_WORD_COUNT = 149
EXPECTED_REVERSAL_FIXED_COUNT = 3
EXPECTED_REVERSAL_ORBIT_COUNT = 76
EXPECTED_CONNECTED_GRAPH_MASK_PAIR_COUNT = 28
EXPECTED_ALL_GRAPH_MASK_PAIR_COUNT = 114
EXPECTED_MAX_MINIMAL_WITNESS_GAP = 5

EXPECTED_CONE_DISTRIBUTION = {
    (6, 4, 1): 1,
    (7, 4, 1): 7,
    (8, 3, 2): 12,
    (8, 4, 1): 11,
    (9, 3, 2): 52,
}
EXPECTED_CHAMBERS_PER_CONE = {1: 25, 2: 50, 3: 8}
EXPECTED_WORD_PROFILE_COUNTS = {
    (6, 8, 2, 4, 4): 1,
    (7, 5, 0, 5, 5): 3,
    (7, 6, 1, 4, 4): 4,
    (8, 4, 0, 4, 3): 18,
    (8, 4, 0, 4, 4): 9,
    (8, 5, 1, 3, 1): 8,
    (8, 5, 1, 3, 3): 2,
    (9, 3, 0, 3, 1): 10,
    (9, 3, 0, 3, 2): 48,
    (9, 3, 0, 3, 3): 28,
    (9, 4, 1, 2, 2): 18,
}
EXPECTED_CUT_DEBT_COUNTS = {
    (6, 2, 3, 2): 1,
    (7, 1, 3, 1): 3,
    (7, 2, 2, 2): 2,
    (7, 2, 3, 2): 2,
    (8, 1, 1, 2): 2,
    (8, 1, 2, 1): 35,
    (9, 1, 1, 1): 76,
    (9, 1, 2, 1): 28,
}
EXPECTED_ONE_CUT_SIZE_COUNTS = {
    (6, "-"): 1,
    (7, "-"): 4,
    (7, "2"): 3,
    (8, "2"): 27,
    (8, "3"): 10,
    (9, "2"): 86,
    (9, "23"): 18,
}
EXPECTED_EXPOSED_ONE_CUT_COUNTS = {
    (6, False): 1,
    (7, False): 7,
    (8, False): 37,
    (9, False): 46,
    (9, True): 58,
}
EXPECTED_REVERSAL_FIXED_BY_SUMSET = {6: 1, 7: 1, 8: 1}

EXPECTED_PRIMITIVE_COUNT = 1_012_441
EXPECTED_CONNECTED_COUNT = 4_617
EXPECTED_COMPONENT_COUNTS = {
    (1, 1): 4_617,
    (1, 2): 71_147,
    (1, 3): 11_797,
    (2, 1): 391,
    (2, 2): 32_343,
    (2, 3): 421_869,
    (3, 4): 470_277,
}
EXPECTED_BOUNDED_SUMSET_COUNTS = {6: 1, 7: 7, 8: 1_265, 9: 3_344}
EXPECTED_BOUNDED_WORD_COUNTS = {6: 1, 7: 7, 8: 37, 9: 104}
EXPECTED_BOUNDED_PROFILE_COUNTS = {
    (6, 8, 2, 4, 4): 1,
    (7, 5, 0, 5, 5): 3,
    (7, 6, 1, 4, 4): 4,
    (8, 4, 0, 4, 3): 782,
    (8, 4, 0, 4, 4): 9,
    (8, 5, 1, 3, 1): 472,
    (8, 5, 1, 3, 3): 2,
    (9, 3, 0, 3, 1): 274,
    (9, 3, 0, 3, 2): 1_516,
    (9, 3, 0, 3, 3): 772,
    (9, 4, 1, 2, 2): 782,
}

EXPECTED_FLAT_SEMANTIC_SHA256 = "2c5e4c1f448cfb96ff6f4f240d615339f26a853ea7406e7801cf3f581e5cbeb9"
EXPECTED_CONE_SEMANTIC_SHA256 = "a29b4a617a677ba6b6e1778e5a1e94eb47de8f9d4cb8a0df8fd25f1017c9e8f0"
EXPECTED_WORD_SEMANTIC_SHA256 = "1c32e4a56baaaecc56014899bf70f8f4c051ec0a31f377c2400b7c67c544eca8"
EXPECTED_BOUNDED_SEMANTIC_SHA256 = "0f8ace8fa981601cdc0a51e9d03e10db232b0c5ee23a99a6b999bbbd7d8dfc6d"


RationalVector = tuple[Fraction, ...]
RationalBasis = tuple[RationalVector, ...]
IntegerVector = tuple[int, ...]
Cell = tuple[int, int]
Fibre = tuple[Cell, ...]


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def semantic_sha256(lines: Iterable[str]) -> str:
    digest = sha256()
    for line in lines:
        digest.update(line.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def canonical_rowspace(rows: Sequence[Sequence[int | Fraction]]) -> RationalBasis:
    """Return the unique reduced-row-echelon basis over Q."""
    if not rows:
        return ()
    matrix = [[Fraction(value) for value in row] for row in rows]
    require(all(len(row) == GAP_COUNT for row in matrix), "row width")
    pivot_row = 0
    for column in range(GAP_COUNT):
        source = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if source is None:
            continue
        matrix[pivot_row], matrix[source] = matrix[source], matrix[pivot_row]
        pivot = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row:
                continue
            multiplier = matrix[row][column]
            if multiplier:
                matrix[row] = [
                    value - multiplier * pivot_value
                    for value, pivot_value in zip(matrix[row], matrix[pivot_row])
                ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    nonzero = [row for row in matrix if any(row)]
    return tuple(tuple(row) for row in nonzero)


def pivot_column(row: RationalVector) -> int:
    return next(index for index, value in enumerate(row) if value)


def row_in_span(basis: RationalBasis, row: Sequence[int | Fraction]) -> bool:
    remainder = [Fraction(value) for value in row]
    for basis_row in basis:
        pivot = pivot_column(basis_row)
        multiplier = remainder[pivot]
        if multiplier:
            remainder = [
                value - multiplier * basis_value
                for value, basis_value in zip(remainder, basis_row)
            ]
    return not any(remainder)


def nullspace_basis(basis: RationalBasis) -> tuple[RationalVector, ...]:
    pivots = tuple(pivot_column(row) for row in basis)
    free_columns = tuple(index for index in range(GAP_COUNT) if index not in pivots)
    vectors = []
    for free in free_columns:
        vector = [Fraction(0) for _ in range(GAP_COUNT)]
        vector[free] = Fraction(1)
        for row, pivot in zip(basis, pivots):
            vector[pivot] = -row[free]
        vectors.append(tuple(vector))
    return tuple(vectors)


def primitive_integer_vector(vector: Sequence[int | Fraction]) -> IntegerVector:
    denominator = 1
    for value in vector:
        denominator = lcm(denominator, Fraction(value).denominator)
    integers = [int(Fraction(value) * denominator) for value in vector]
    divisor = reduce(gcd, (abs(value) for value in integers if value), 0)
    require(divisor > 0, "nonzero primitive vector")
    return tuple(value // divisor for value in integers)


def primitive_direction(vector: tuple[int, int]) -> tuple[int, int]:
    divisor = gcd(abs(vector[0]), abs(vector[1]))
    require(divisor > 0, "nonzero direction")
    return vector[0] // divisor, vector[1] // divisor


def dot(left: Sequence[int | Fraction], right: Sequence[int | Fraction]) -> Fraction:
    return sum((Fraction(a) * Fraction(b) for a, b in zip(left, right)), Fraction(0))


def angle_half(vector: tuple[int, int]) -> int:
    x_value, y_value = vector
    return 0 if y_value > 0 or (y_value == 0 and x_value > 0) else 1


def compare_directions(left: tuple[int, int], right: tuple[int, int]) -> int:
    left_half = angle_half(left)
    right_half = angle_half(right)
    if left_half != right_half:
        return -1 if left_half < right_half else 1
    cross = left[0] * right[1] - left[1] * right[0]
    if cross:
        return -1 if cross > 0 else 1
    return (left > right) - (left < right)


def enumerate_flat_lattice() -> tuple[RationalBasis, ...]:
    seen: set[RationalBasis] = {()}
    frontier = [()]
    while frontier:
        basis = frontier.pop()
        rank = len(basis)
        for hyperplane in HYPERPLANE_ROWS:
            enlarged = canonical_rowspace(basis + (hyperplane,))
            if len(enlarged) > rank and enlarged not in seen:
                seen.add(enlarged)
                frontier.append(enlarged)
    require(len(seen) == EXPECTED_FLAT_COUNT, "flat count")
    return tuple(sorted(seen, key=lambda item: (len(item), item)))


def closure_mask(basis: RationalBasis) -> int:
    mask = 0
    for index, hyperplane in enumerate(HYPERPLANE_ROWS):
        if row_in_span(basis, hyperplane):
            mask |= 1 << index
    return mask


def difference_graph_masks_from_closure(mask: int) -> tuple[int, int]:
    a_mask = 0
    b_mask = 0
    for index in range(len(HYPERPLANE_ROWS)):
        if mask >> index & 1:
            a_index, b_index = divmod(index, len(B_DIFFERENCE_ROWS))
            a_mask |= 1 << a_index
            b_mask |= 1 << b_index
    return a_mask, b_mask


def component_count(vertex_count: int, edge_mask: int, pairs: Sequence[tuple[int, int]]) -> int:
    parent = list(range(vertex_count))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for index, (left, right) in enumerate(pairs):
        if edge_mask >> index & 1:
            left_root = find(left)
            right_root = find(right)
            if left_root != right_root:
                parent[right_root] = left_root
    return len({find(vertex) for vertex in range(vertex_count)})


A_COMPONENTS = tuple(
    component_count(3, mask, A_EDGE_PAIRS) for mask in range(1 << len(A_EDGE_PAIRS))
)
B_COMPONENTS = tuple(
    component_count(4, mask, B_EDGE_PAIRS) for mask in range(1 << len(B_EDGE_PAIRS))
)


def exact_equality_mask(gaps: IntegerVector) -> int:
    mask = 0
    for index, hyperplane in enumerate(HYPERPLANE_ROWS):
        if dot(hyperplane, gaps) == 0:
            mask |= 1 << index
    return mask


def difference_graph_masks(gaps: IntegerVector) -> tuple[int, int]:
    x_gap, y_gap, u_gap, v_gap, w_gap = gaps
    a_differences = (x_gap, x_gap + y_gap, y_gap)
    b_differences = (
        u_gap,
        u_gap + v_gap,
        u_gap + v_gap + w_gap,
        v_gap,
        v_gap + w_gap,
        w_gap,
    )
    a_mask = 0
    b_mask = 0
    for a_index, a_value in enumerate(a_differences):
        for b_index, b_value in enumerate(b_differences):
            if a_value == b_value:
                a_mask |= 1 << a_index
                b_mask |= 1 << b_index
    return a_mask, b_mask


def supports(gaps: IntegerVector) -> tuple[tuple[int, ...], tuple[int, ...]]:
    x_gap, y_gap, u_gap, v_gap, w_gap = gaps
    return (0, x_gap, x_gap + y_gap), (
        0,
        u_gap,
        u_gap + v_gap,
        u_gap + v_gap + w_gap,
    )


def sum_fibres(gaps: IntegerVector) -> tuple[Fibre, ...]:
    a_support, b_support = supports(gaps)
    by_sum: dict[int, list[Cell]] = defaultdict(list)
    for a_index, a_value in enumerate(a_support):
        for b_index, b_value in enumerate(b_support):
            by_sum[a_value + b_value].append((a_index, b_index))
    return tuple(
        tuple(sorted(by_sum[sum_value])) for sum_value in sorted(by_sum)
    )


def fibre_word(fibres: Sequence[Fibre]) -> str:
    return "|".join(
        "=".join(f"{a_index}{b_index}" for a_index, b_index in fibre)
        for fibre in fibres
    )


def projection_masks(fibres: Sequence[Fibre], omitted: frozenset[int] = frozenset()) -> tuple[int, int]:
    a_mask = 0
    b_mask = 0
    a_pair_index = {pair: index for index, pair in enumerate(A_EDGE_PAIRS)}
    b_pair_index = {pair: index for index, pair in enumerate(B_EDGE_PAIRS)}
    for fibre_index, fibre in enumerate(fibres):
        if fibre_index in omitted:
            continue
        for left, right in combinations(fibre, 2):
            a_pair = tuple(sorted((left[0], right[0])))
            b_pair = tuple(sorted((left[1], right[1])))
            require(a_pair[0] != a_pair[1] and b_pair[0] != b_pair[1], "collision rectangle")
            a_mask |= 1 << a_pair_index[a_pair]
            b_mask |= 1 << b_pair_index[b_pair]
    return a_mask, b_mask


def rectangle_exposed(gaps: IntegerVector, fibre: Fibre) -> bool:
    if len(fibre) != 2:
        return False
    (a_left, b_left), (a_right, b_right) = fibre
    a_support, b_support = supports(gaps)
    multiplicities = Counter(
        a_value + b_value for a_value in a_support for b_value in b_support
    )
    first_cross = a_support[a_left] + b_support[b_right]
    second_cross = a_support[a_right] + b_support[b_left]
    return multiplicities[first_cross] == 1 and multiplicities[second_cross] == 1


def fibre_cut_data(gaps: IntegerVector, fibres: Sequence[Fibre]) -> tuple[int, int, int, str, bool]:
    candidates = tuple(index for index, fibre in enumerate(fibres) if len(fibre) >= 2)
    cut_a = len(candidates) + 1
    cut_b = len(candidates) + 1
    cut_either = len(candidates) + 1
    one_cut_sizes: set[int] = set()
    exposed_one_cut = False
    for count in range(1, len(candidates) + 1):
        for selection in combinations(candidates, count):
            a_mask, b_mask = projection_masks(fibres, frozenset(selection))
            a_disconnected = A_COMPONENTS[a_mask] > 1
            b_disconnected = B_COMPONENTS[b_mask] > 1
            if a_disconnected:
                cut_a = min(cut_a, count)
            if b_disconnected:
                cut_b = min(cut_b, count)
            if a_disconnected or b_disconnected:
                cut_either = min(cut_either, count)
                if count == 1:
                    fibre = fibres[selection[0]]
                    one_cut_sizes.add(len(fibre))
                    if rectangle_exposed(gaps, fibre):
                        exposed_one_cut = True
    require(cut_a <= len(candidates), "finite A cut debt")
    require(cut_b <= len(candidates), "finite B cut debt")
    require(cut_either == min(cut_a, cut_b), "either cut debt")
    size_word = "".join(str(size) for size in sorted(one_cut_sizes)) or "-"
    return cut_either, cut_a, cut_b, size_word, exposed_one_cut


@dataclass(frozen=True)
class WordData:
    word: str
    gaps: IntegerVector
    sumset_size: int
    collision_count: int
    triple_count: int
    scalar_candidates: int
    unexposed_candidates: int
    cut_either: int
    cut_a: int
    cut_b: int
    one_cut_sizes: str
    exposed_one_cut: bool


def analyze_word(gaps: IntegerVector) -> WordData:
    require(len(gaps) == GAP_COUNT and min(gaps) > 0, "positive five-gap vector")
    fibres = sum_fibres(gaps)
    multiplicities = tuple(len(fibre) for fibre in fibres)
    require(max(multiplicities) <= 3 and sum(multiplicities) == 12, "fibre multiplicities")
    collision_count = sum(size * (size - 1) // 2 for size in multiplicities)
    triple_count = sum(size == 3 for size in multiplicities)
    require(len(fibres) == 12 - collision_count + triple_count, "sumset identity")
    candidates = tuple(fibre for fibre in fibres if len(fibre) >= 2)
    unexposed = sum(not rectangle_exposed(gaps, fibre) for fibre in candidates)
    a_mask, b_mask = projection_masks(fibres)
    require(A_COMPONENTS[a_mask] == 1 and B_COMPONENTS[b_mask] == 1, "connected word")
    require((a_mask, b_mask) == difference_graph_masks(gaps), "projection/difference graph agreement")
    cut_either, cut_a, cut_b, cut_sizes, exposed_cut = fibre_cut_data(gaps, fibres)
    return WordData(
        word=fibre_word(fibres),
        gaps=gaps,
        sumset_size=len(fibres),
        collision_count=collision_count,
        triple_count=triple_count,
        scalar_candidates=len(candidates),
        unexposed_candidates=unexposed,
        cut_either=cut_either,
        cut_a=cut_a,
        cut_b=cut_b,
        one_cut_sizes=cut_sizes,
        exposed_one_cut=exposed_cut,
    )


def chamber_samples(basis: RationalBasis, mask: int) -> tuple[IntegerVector, ...]:
    dimension = GAP_COUNT - len(basis)
    if dimension == 0:
        return ()
    nullspace = tuple(primitive_integer_vector(vector) for vector in nullspace_basis(basis))
    require(len(nullspace) == dimension, "nullspace dimension")
    samples: set[IntegerVector] = set()
    if dimension == 1:
        for sign in (1, -1):
            candidate = tuple(sign * value for value in nullspace[0])
            if min(candidate) > 0 and exact_equality_mask(candidate) == mask:
                samples.add(primitive_integer_vector(candidate))
        return tuple(sorted(samples))

    require(dimension == 2, "connected flat dimension at most two")
    first, second = nullspace
    coordinate_normals = tuple(
        (first[index], second[index]) for index in range(GAP_COUNT)
    )
    # If a coordinate form vanishes identically on this flat, the flat cannot
    # meet the strict positive orthant.
    if any(normal == (0, 0) for normal in coordinate_normals):
        return ()
    restricted_normals: set[tuple[int, int]] = set(coordinate_normals)
    for hyperplane_index, hyperplane in enumerate(HYPERPLANE_ROWS):
        if not (mask >> hyperplane_index & 1):
            normal = (int(dot(hyperplane, first)), int(dot(hyperplane, second)))
            require(normal != (0, 0), "nonclosure hyperplane restricts nontrivially")
            restricted_normals.add(normal)

    rays: set[tuple[int, int]] = set()
    for normal_x, normal_y in restricted_normals:
        require((normal_x, normal_y) != (0, 0), "nonzero chamber normal")
        ray = primitive_direction((normal_y, -normal_x))
        rays.add(ray)
        rays.add((-ray[0], -ray[1]))
    ordered_rays = sorted(rays, key=cmp_to_key(compare_directions))
    require(len(ordered_rays) >= 4 and len(ordered_rays) % 2 == 0, "central ray arrangement")

    for index, left in enumerate(ordered_rays):
        right = ordered_rays[(index + 1) % len(ordered_rays)]
        coefficient = (left[0] + right[0], left[1] + right[1])
        if coefficient == (0, 0):
            continue
        coefficient = primitive_direction(coefficient)
        candidate = tuple(
            coefficient[0] * first[coordinate] + coefficient[1] * second[coordinate]
            for coordinate in range(GAP_COUNT)
        )
        candidate = primitive_integer_vector(candidate)
        if min(candidate) <= 0 or exact_equality_mask(candidate) != mask:
            continue
        samples.add(candidate)
    return tuple(sorted(samples))


@dataclass(frozen=True)
class ConeData:
    basis: RationalBasis
    mask: int
    rank: int
    dimension: int
    words: tuple[str, ...]
    samples: tuple[IntegerVector, ...]


def enumerate_positive_connected_cones(flats: Sequence[RationalBasis]) -> tuple[ConeData, ...]:
    cones: list[ConeData] = []
    for basis in flats:
        mask = closure_mask(basis)
        a_mask, b_mask = difference_graph_masks_from_closure(mask)
        if A_COMPONENTS[a_mask] != 1 or B_COMPONENTS[b_mask] != 1:
            continue
        dimension = GAP_COUNT - len(basis)
        require(dimension <= 2, "connected collision flat has dimension at most two")
        samples = chamber_samples(basis, mask)
        if not samples:
            continue
        words = tuple(sorted({fibre_word(sum_fibres(sample)) for sample in samples}))
        require(len(words) == len(samples), "one exact sample per oriented chamber")
        cones.append(
            ConeData(
                basis=basis,
                mask=mask,
                rank=len(basis),
                dimension=dimension,
                words=words,
                samples=samples,
            )
        )
    cones.sort(key=lambda cone: (len(sum_fibres(cone.samples[0])), cone.rank, cone.mask))
    require(len(cones) == EXPECTED_POSITIVE_CONNECTED_CONE_COUNT, "positive connected cone count")
    return tuple(cones)


@dataclass(frozen=True)
class BoundedCensus:
    primitive_count: int
    connected_count: int
    component_counts: Counter[tuple[int, int]]
    graph_mask_pairs: frozenset[tuple[int, int]]
    connected_graph_mask_pairs: frozenset[tuple[int, int]]
    sumset_counts: Counter[int]
    words_by_sumset: dict[int, frozenset[str]]
    profile_counts: Counter[tuple[int, int, int, int, int]]
    best_witness: dict[str, IntegerVector]
    semantic_hash: str


def witness_key(gaps: IntegerVector) -> tuple[int, int, IntegerVector]:
    return max(gaps), sum(gaps), gaps


def hostile_bounded_census() -> BoundedCensus:
    primitive_count = 0
    connected_count = 0
    component_counts: Counter[tuple[int, int]] = Counter()
    graph_mask_pairs: set[tuple[int, int]] = set()
    connected_graph_mask_pairs: set[tuple[int, int]] = set()
    sumset_counts: Counter[int] = Counter()
    words_by_sumset_mutable: dict[int, set[str]] = defaultdict(set)
    profile_counts: Counter[tuple[int, int, int, int, int]] = Counter()
    best_witness: dict[str, IntegerVector] = {}
    digest = sha256()

    for gaps_raw in product(range(1, HOSTILE_BOUND + 1), repeat=GAP_COUNT):
        if gcd(*gaps_raw) != 1:
            continue
        gaps = tuple(gaps_raw)
        primitive_count += 1
        a_mask, b_mask = difference_graph_masks(gaps)
        a_components = A_COMPONENTS[a_mask]
        b_components = B_COMPONENTS[b_mask]
        component_counts[a_components, b_components] += 1
        graph_mask_pairs.add((a_mask, b_mask))
        connected = a_components == 1 and b_components == 1
        sumset_size = 0
        if connected:
            connected_count += 1
            connected_graph_mask_pairs.add((a_mask, b_mask))
            data = analyze_word(gaps)
            sumset_size = data.sumset_size
            sumset_counts[sumset_size] += 1
            words_by_sumset_mutable[sumset_size].add(data.word)
            profile = (
                data.sumset_size,
                data.collision_count,
                data.triple_count,
                data.scalar_candidates,
                data.unexposed_candidates,
            )
            profile_counts[profile] += 1
            current = best_witness.get(data.word)
            if current is None or witness_key(gaps) < witness_key(current):
                best_witness[data.word] = gaps
        digest.update(bytes((*gaps, a_mask, b_mask, sumset_size)))

    words_by_sumset = {
        size: frozenset(words) for size, words in words_by_sumset_mutable.items()
    }
    return BoundedCensus(
        primitive_count=primitive_count,
        connected_count=connected_count,
        component_counts=component_counts,
        graph_mask_pairs=frozenset(graph_mask_pairs),
        connected_graph_mask_pairs=frozenset(connected_graph_mask_pairs),
        sumset_counts=sumset_counts,
        words_by_sumset=words_by_sumset,
        profile_counts=profile_counts,
        best_witness=best_witness,
        semantic_hash=digest.hexdigest(),
    )


def reverse_gaps(gaps: IntegerVector) -> IntegerVector:
    x_gap, y_gap, u_gap, v_gap, w_gap = gaps
    return y_gap, x_gap, w_gap, v_gap, u_gap


def format_counter(counter: Counter | dict) -> str:
    return ",".join(f"{key}:{counter[key]}" for key in sorted(counter, key=str))


def format_gap(gaps: IntegerVector) -> str:
    return ",".join(str(value) for value in gaps)


def format_mask(mask: int) -> str:
    return f"0x{mask:05x}"


def main() -> None:
    require(len(HYPERPLANE_ROWS) == 18 and len(set(HYPERPLANE_ROWS)) == 18, "18 hyperplanes")
    flats = enumerate_flat_lattice()
    flat_masks = tuple((len(basis), closure_mask(basis)) for basis in flats)
    require(len({mask for _, mask in flat_masks}) == len(flat_masks), "flat closure masks unique")
    flat_hash = semantic_sha256(
        f"rank={rank};closure={format_mask(mask)}" for rank, mask in flat_masks
    )

    cones = enumerate_positive_connected_cones(flats)
    global_words = {word for cone in cones for word in cone.words}
    require(sum(len(cone.words) for cone in cones) == len(global_words), "word belongs to one cone")
    require(len(global_words) == EXPECTED_GLOBAL_WORD_COUNT, "global oriented word count")

    cone_distribution: Counter[tuple[int, int, int]] = Counter()
    chamber_distribution: Counter[int] = Counter()
    for cone in cones:
        sumset_size = len(sum_fibres(cone.samples[0]))
        require(all(len(sum_fibres(sample)) == sumset_size for sample in cone.samples), "cone sumset size")
        cone_distribution[sumset_size, cone.rank, cone.dimension] += 1
        chamber_distribution[len(cone.words)] += 1
    require(dict(cone_distribution) == EXPECTED_CONE_DISTRIBUTION, "cone distribution")
    require(dict(chamber_distribution) == EXPECTED_CHAMBERS_PER_CONE, "chambers per cone")

    census = hostile_bounded_census()
    require(census.primitive_count == EXPECTED_PRIMITIVE_COUNT, "bounded primitive count")
    require(census.connected_count == EXPECTED_CONNECTED_COUNT, "bounded connected count")
    require(dict(census.component_counts) == EXPECTED_COMPONENT_COUNTS, "bounded component counts")
    require(len(census.graph_mask_pairs) == EXPECTED_ALL_GRAPH_MASK_PAIR_COUNT, "all graph-mask pairs")
    require(
        len(census.connected_graph_mask_pairs) == EXPECTED_CONNECTED_GRAPH_MASK_PAIR_COUNT,
        "connected graph-mask pairs",
    )
    require(dict(census.sumset_counts) == EXPECTED_BOUNDED_SUMSET_COUNTS, "bounded sumset counts")
    require(
        {size: len(words) for size, words in census.words_by_sumset.items()}
        == EXPECTED_BOUNDED_WORD_COUNTS,
        "bounded word counts",
    )
    require(dict(census.profile_counts) == EXPECTED_BOUNDED_PROFILE_COUNTS, "bounded profile counts")
    bounded_words = set().union(*census.words_by_sumset.values())
    require(bounded_words == global_words, "hostile census realizes every global word")
    require(
        max(max(gaps) for gaps in census.best_witness.values()) == EXPECTED_MAX_MINIMAL_WITNESS_GAP,
        "maximum minimal witness gap",
    )

    word_data = {
        word: analyze_word(census.best_witness[word]) for word in sorted(global_words)
    }
    word_profile_counts: Counter[tuple[int, int, int, int, int]] = Counter()
    cut_debt_counts: Counter[tuple[int, int, int, int]] = Counter()
    one_cut_size_counts: Counter[tuple[int, str]] = Counter()
    exposed_one_cut_counts: Counter[tuple[int, bool]] = Counter()
    for data in word_data.values():
        word_profile_counts[
            data.sumset_size,
            data.collision_count,
            data.triple_count,
            data.scalar_candidates,
            data.unexposed_candidates,
        ] += 1
        cut_debt_counts[data.sumset_size, data.cut_either, data.cut_a, data.cut_b] += 1
        one_cut_size_counts[data.sumset_size, data.one_cut_sizes] += 1
        exposed_one_cut_counts[data.sumset_size, data.exposed_one_cut] += 1
    require(dict(word_profile_counts) == EXPECTED_WORD_PROFILE_COUNTS, "word profiles")
    require(dict(cut_debt_counts) == EXPECTED_CUT_DEBT_COUNTS, "cut debts")
    require(dict(one_cut_size_counts) == EXPECTED_ONE_CUT_SIZE_COUNTS, "one-cut sizes")
    require(dict(exposed_one_cut_counts) == EXPECTED_EXPOSED_ONE_CUT_COUNTS, "exposed one-cuts")

    ordered_words = tuple(sorted(global_words, key=lambda word: (word_data[word].sumset_size, word)))
    word_ids = {word: f"W{index:03d}" for index, word in enumerate(ordered_words, 1)}
    reversal = {
        word: fibre_word(sum_fibres(reverse_gaps(word_data[word].gaps)))
        for word in ordered_words
    }
    require(all(reversal[word] in global_words for word in ordered_words), "reversal closure")
    require(all(reversal[reversal[word]] == word for word in ordered_words), "reversal involution")
    fixed_words = tuple(word for word in ordered_words if reversal[word] == word)
    require(len(fixed_words) == EXPECTED_REVERSAL_FIXED_COUNT, "reversal fixed count")
    fixed_by_sumset = Counter(word_data[word].sumset_size for word in fixed_words)
    require(dict(fixed_by_sumset) == EXPECTED_REVERSAL_FIXED_BY_SUMSET, "reversal fixed levels")
    reversal_orbits = []
    unseen = set(ordered_words)
    while unseen:
        word = min(unseen, key=lambda item: int(word_ids[item][1:]))
        partner = reversal[word]
        orbit = tuple(sorted((word, partner), key=lambda item: int(word_ids[item][1:])))
        reversal_orbits.append(tuple(dict.fromkeys(orbit)))
        unseen.discard(word)
        unseen.discard(partner)
    require(len(reversal_orbits) == EXPECTED_REVERSAL_ORBIT_COUNT, "reversal orbit count")

    ordered_cones = tuple(sorted(cones, key=lambda cone: (len(sum_fibres(cone.samples[0])), cone.rank, cone.mask)))
    cone_ids = {cone.mask: f"C{index:03d}" for index, cone in enumerate(ordered_cones, 1)}
    word_to_cone: dict[str, str] = {}
    for cone in ordered_cones:
        for word in cone.words:
            require(word not in word_to_cone, "unique word cone")
            word_to_cone[word] = cone_ids[cone.mask]

    cone_hash = semantic_sha256(
        (
            f"{cone_ids[cone.mask]};rank={cone.rank};dim={cone.dimension};"
            f"closure={format_mask(cone.mask)};"
            f"words={','.join(word_ids[word] for word in cone.words)}"
        )
        for cone in ordered_cones
    )
    word_hash = semantic_sha256(
        (
            f"{word_ids[word]};cone={word_to_cone[word]};word={word};"
            f"sample={format_gap(word_data[word].gaps)};"
            f"m={word_data[word].sumset_size};C={word_data[word].collision_count};"
            f"T={word_data[word].triple_count};raw={word_data[word].scalar_candidates};"
            f"live={word_data[word].unexposed_candidates};"
            f"cuts={word_data[word].cut_either}/{word_data[word].cut_a}/{word_data[word].cut_b};"
            f"one_cut_sizes={word_data[word].one_cut_sizes};"
            f"exposed_one_cut={int(word_data[word].exposed_one_cut)};"
            f"reverse={word_ids[reversal[word]]}"
        )
        for word in ordered_words
    )

    require(flat_hash == EXPECTED_FLAT_SEMANTIC_SHA256, "flat semantic hash")
    require(cone_hash == EXPECTED_CONE_SEMANTIC_SHA256, "cone semantic hash")
    require(word_hash == EXPECTED_WORD_SEMANTIC_SHA256, "word semantic hash")
    require(census.semantic_hash == EXPECTED_BOUNDED_SEMANTIC_SHA256, "bounded semantic hash")

    source = Path(__file__).read_bytes().replace(b"\r\n", b"\n")

    print("THM-3603 EXPONENT-TWO 3-BY-4 ADDITIVE SUPPORT ATLAS")
    print("status=FINITE-EXACT+VERIFIED-EXACT;claim_scope=additive_support_only")
    print("variables=X,Y,U,V,W;all_positive;primitive_integral_bounded_control=yes")
    print("output_exchange=three_point_support_named_A;simultaneous_reversal=retained")
    print("hyperplane_count=18;flat_count=786;positive_connected_cone_count=83")
    print("global_oriented_fibre_word_count=149;reversal_fixed=3;reversal_orbits=76")
    print("cone_distribution=" + format_counter(cone_distribution))
    print("chambers_per_cone=" + format_counter(chamber_distribution))
    print("word_profile_counts=" + format_counter(word_profile_counts))
    print("cut_debt_counts=" + format_counter(cut_debt_counts))
    print("one_cut_size_counts=" + format_counter(one_cut_size_counts))
    print("exposed_one_cut_counts=" + format_counter(exposed_one_cut_counts))
    print(f"flat_semantic_sha256={flat_hash}")
    print(f"cone_semantic_sha256={cone_hash}")
    print(f"word_semantic_sha256={word_hash}")
    print(f"bounded_semantic_sha256={census.semantic_hash}")
    print(f"source_lf_sha256={sha256(source).hexdigest()}")
    print("HYPERPLANES")
    for index, label in enumerate(HYPERPLANE_LABELS):
        row = HYPERPLANE_ROWS[index]
        print(f"H{index:02d}={label};row={','.join(str(value) for value in row)}")
    print("FLAT_LATTICE")
    for index, (rank, mask) in enumerate(flat_masks, 1):
        print(f"F{index:03d}=rank:{rank};closure:{format_mask(mask)}")
    print("POSITIVE_CONNECTED_CONES")
    for cone in ordered_cones:
        cone_id = cone_ids[cone.mask]
        sumset_size = len(sum_fibres(cone.samples[0]))
        closure = ",".join(
            f"H{index:02d}" for index in range(len(HYPERPLANE_ROWS)) if cone.mask >> index & 1
        )
        words = ",".join(word_ids[word] for word in cone.words)
        sample = min((census.best_witness[word] for word in cone.words), key=witness_key)
        print(
            f"{cone_id}=m:{sumset_size};rank:{cone.rank};dim:{cone.dimension};"
            f"closure:{closure};words:{words};sample:{format_gap(sample)}"
        )
    print("ORIENTED_FIBRE_WORDS")
    for word in ordered_words:
        data = word_data[word]
        print(
            f"{word_ids[word]}=cone:{word_to_cone[word]};m:{data.sumset_size};"
            f"C:{data.collision_count};T:{data.triple_count};"
            f"raw:{data.scalar_candidates};live:{data.unexposed_candidates};"
            f"cuts:{data.cut_either}/{data.cut_a}/{data.cut_b};"
            f"one_cut_sizes:{data.one_cut_sizes};"
            f"exposed_one_cut:{int(data.exposed_one_cut)};"
            f"reverse:{word_ids[reversal[word]]};sample:{format_gap(data.gaps)};word:{word}"
        )
    print("REVERSAL_ORBITS")
    print(
        "orbits="
        + ",".join(
            "/".join(word_ids[word] for word in orbit) for orbit in reversal_orbits
        )
    )
    print("HOSTILE_BOUNDED_CONTROL")
    print(f"bound={HOSTILE_BOUND};primitive_count={census.primitive_count};connected_count={census.connected_count}")
    print("component_counts=" + format_counter(census.component_counts))
    print(
        f"all_graph_mask_pairs={len(census.graph_mask_pairs)};"
        f"connected_graph_mask_pairs={len(census.connected_graph_mask_pairs)}"
    )
    print("connected_sumset_counts=" + format_counter(census.sumset_counts))
    print(
        "connected_word_counts="
        + format_counter({size: len(words) for size, words in census.words_by_sumset.items()})
    )
    print("bounded_profile_counts=" + format_counter(census.profile_counts))
    print(
        f"bounded_words_equal_global=1;max_minimal_witness_gap="
        f"{max(max(gaps) for gaps in census.best_witness.values())}"
    )
    print(
        "boundary=no_Darboux_or_Jacobian_nonentry_claim;"
        "a_fibre_cut_removes_collision_edges_only_not_input_cells_or_coefficients;"
        "rectangle_exposed_is_an_additive_candidate_filter_not_a_validated_algebraic_deletion;"
        "the_bound_16_scan_controls_but_does_not_supply_globality;"
        "globality_comes_from_the_exact_Q_flat_and_chamber_enumeration"
    )


if __name__ == "__main__":
    main()
