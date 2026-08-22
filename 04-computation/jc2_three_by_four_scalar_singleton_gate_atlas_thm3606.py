#!/usr/bin/env python3
"""Exact scalar-arm/singleton necessary-gate atlas for THM-3606.

The parent THM-3603 classifies 149 connected oriented additive fibre words for
three weights against four.  For every collision fibre, every address in it,
and both scalar-arm orientations, this script fixes the absolute weights by
placing (-2,1) or (1,-2) at that address.  It then applies the universal
singleton-row gate to every singleton fibre.

This is a necessary support/weight test only.  It does not solve any remaining
multi-address differential equation and does not claim a Darboux pair.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import importlib.util
from itertools import product
import json
from math import lcm
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT / "04-computation/jc2_three_by_four_additive_support_atlas_thm3603.py"
EXPECTED_PARENT_SHA256 = "85c8a5e8200e0999b0e671f2419f3b9a52f03c862370cfa5a6afc7aaf8715c74"
HOSTILE_BOUND = 16
EXPECTED_EXACT_IDS = (
    "W001", "W002", "W003", "W004", "W005", "W006", "W007", "W008",
    "W009", "W010", "W011", "W012", "W015", "W016", "W018", "W020",
    "W022", "W023", "W024", "W030", "W031", "W032", "W033", "W034",
    "W035", "W036", "W037", "W038", "W040", "W041", "W042", "W043",
    "W044", "W045", "W072", "W073", "W141", "W149",
)
EXPECTED_REPRESENTATIVE_IDS = (
    "W001", "W003", "W004", "W005", "W006", "W007", "W010", "W015",
    "W016", "W032", "W033", "W036", "W037", "W040", "W044", "W045",
    "W072", "W149",
)
EXPECTED_WORD_COUNTS = {6: 1, 7: 7, 8: 26, 9: 4}
EXPECTED_SCHEME_COUNTS = {6: 12, 7: 66, 8: 90, 9: 10}
EXPECTED_ORIENTATION_COUNTS = {"+-": 90, "-+": 88}
EXPECTED_SCALAR_SIZE_COUNTS = {2: 161, 3: 17}
EXPECTED_SCALAR_CUT_COUNTS = {
    (6, False): 12,
    (7, False): 62,
    (7, True): 4,
    (8, False): 70,
    (8, True): 20,
    (9, True): 10,
}
EXPECTED_SCALAR_EXPOSED_COUNTS = {
    (6, False): 12,
    (7, False): 66,
    (8, False): 66,
    (8, True): 24,
    (9, False): 4,
    (9, True): 6,
}
EXPECTED_UNEXPOSED_WORD_IDS = (
    "W001", "W002", "W003", "W004", "W005", "W006", "W007", "W008",
    "W009", "W010", "W011", "W012", "W015", "W016", "W022", "W030",
    "W031", "W032", "W033", "W034", "W036", "W037", "W040", "W041",
    "W042", "W043", "W045", "W072", "W073", "W141", "W149",
)
EXPECTED_MAX_SCHEME_WITNESS_GAP = 12
EXPECTED_BOUNDED = {
    "connected": 6724,
    "surviving": 1288,
    "placements": 3790,
    "by_sumset": {6: 16, 7: 54, 8: 1064, 9: 154},
}
EXPECTED_SEMANTIC_SHA256 = "c21d12829dc415dded4fcd2347c23137991904da598cbe4cd8ad4e8844f5c1f0"
EXPECTED_BOUNDED_SEMANTIC_SHA256 = "68437764f3eb632545e70224b265f1e74b01a30826d2cac5c6658e150358ff07"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def semantic_sha256(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


require(lf_sha256(PARENT) == EXPECTED_PARENT_SHA256, "THM-3603 parent hash")
spec = importlib.util.spec_from_file_location("thm3603_parent", PARENT)
require(spec is not None and spec.loader is not None, "parent import spec")
parent = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = parent
spec.loader.exec_module(parent)


@dataclass(frozen=True)
class Placement:
    scalar_index: int
    scalar_fibre: tuple[tuple[int, int], ...]
    p_weights: tuple[int, int, int]
    q_weights: tuple[int, int, int, int]
    eligible: tuple[tuple[int, int, str], ...]
    singleton_rows: tuple[tuple[int, int, int, int, str], ...]
    nontrivial_singleton_components: tuple[tuple[str, ...], ...]


def fibre_label(fibre: tuple[tuple[int, int], ...]) -> str:
    return "=".join(f"{i}{j}" for i, j in fibre)


def singleton_verdict(p_weight: int, q_weight: int) -> str:
    if p_weight == 0 and q_weight == 0:
        return "00"
    if p_weight > 0 and q_weight > 0:
        return "++"
    if p_weight < 0 and q_weight < 0:
        return "--"
    return "FAIL"


def singleton_components(
    rows: tuple[tuple[int, int, int, int, str], ...]
) -> tuple[tuple[str, ...], ...]:
    vertices = tuple([f"A{i}" for i in range(3)] + [f"B{j}" for j in range(4)])
    parent_map = {vertex: vertex for vertex in vertices}

    def find(vertex: str) -> str:
        while parent_map[vertex] != vertex:
            parent_map[vertex] = parent_map[parent_map[vertex]]
            vertex = parent_map[vertex]
        return vertex

    active: set[str] = set()
    for i, j, _, _, verdict in rows:
        if verdict not in {"++", "--"}:
            continue
        left, right = f"A{i}", f"B{j}"
        active.update((left, right))
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent_map[right_root] = left_root
    parts: dict[str, list[str]] = {}
    for vertex in sorted(active):
        parts.setdefault(find(vertex), []).append(vertex)
    return tuple(sorted(tuple(part) for part in parts.values()))


def placements_for_sample(gaps: tuple[int, ...]) -> tuple[Placement, ...]:
    a_support, b_support = parent.supports(gaps)
    fibres = parent.sum_fibres(gaps)
    placements: dict[tuple[int, int, int], Placement] = {}
    for scalar_index, scalar_fibre in enumerate(fibres):
        if len(scalar_fibre) < 2:
            continue
        candidates: set[tuple[int, int]] = set()
        for i, j in scalar_fibre:
            # Place either (-2,1) or (1,-2) at this scalar address.
            candidates.add((-2 - a_support[i], 1 - b_support[j]))
            candidates.add((1 - a_support[i], -2 - b_support[j]))
        for p0, q0 in sorted(candidates):
            p_weights = tuple(p0 + value for value in a_support)
            q_weights = tuple(q0 + value for value in b_support)
            require(
                all(p_weights[i] + q_weights[j] + 1 == 0 for i, j in scalar_fibre),
                "scalar fibre weight",
            )
            eligible = tuple(
                (i, j, "-+") if (p_weights[i], q_weights[j]) == (-2, 1)
                else (i, j, "+-")
                for i, j in scalar_fibre
                if (p_weights[i], q_weights[j]) in {(-2, 1), (1, -2)}
            )
            require(eligible, "at least one scalar-arm address")
            singleton_rows = tuple(
                (
                    fibre[0][0],
                    fibre[0][1],
                    p_weights[fibre[0][0]],
                    q_weights[fibre[0][1]],
                    singleton_verdict(
                        p_weights[fibre[0][0]], q_weights[fibre[0][1]]
                    ),
                )
                for index, fibre in enumerate(fibres)
                if index != scalar_index and len(fibre) == 1
            )
            if any(row[-1] == "FAIL" for row in singleton_rows):
                continue
            placement = Placement(
                scalar_index=scalar_index,
                scalar_fibre=scalar_fibre,
                p_weights=p_weights,
                q_weights=q_weights,
                eligible=eligible,
                singleton_rows=singleton_rows,
                nontrivial_singleton_components=singleton_components(singleton_rows),
            )
            key = scalar_index, p0, q0
            require(key not in placements, "placement deduplication")
            placements[key] = placement
    return tuple(
        placements[key]
        for key in sorted(placements)
    )


def hard_rectangle_survives(gaps: tuple[int, ...]) -> bool:
    """Ignore soft cross-quadrant exceptions; retain only absolute no-go cells."""
    fibres = parent.sum_fibres(gaps)
    singleton_cells = tuple(fibre[0] for fibre in fibres if len(fibre) == 1)
    for scalar_fibre in fibres:
        if len(scalar_fibre) < 2:
            continue
        for i, j in scalar_fibre:
            # Scalar orientation (-2,1): k<=i,l>=j is always opposite sign.
            if not any(k <= i and ell >= j for k, ell in singleton_cells):
                return True
            # Scalar orientation (1,-2): k>=i,l<=j is always opposite sign.
            if not any(k >= i and ell <= j for k, ell in singleton_cells):
                return True
    return False


A_ROWS = (
    (0, 0, 0, 0, 0),
    (1, 0, 0, 0, 0),
    (1, 1, 0, 0, 0),
)
B_ROWS = (
    (0, 0, 0, 0, 0),
    (0, 0, 1, 0, 0),
    (0, 0, 1, 1, 0),
    (0, 0, 1, 1, 1),
)


def row_add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(a + b for a, b in zip(left, right))


def row_sub(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(a - b for a, b in zip(left, right))


def cell_row(cell: tuple[int, int]) -> tuple[int, ...]:
    return row_add(A_ROWS[cell[0]], B_ROWS[cell[1]])


def distance_row(side: str, low: int, high: int) -> tuple[int, ...]:
    require(low < high, "ordered support distance")
    rows = A_ROWS if side == "A" else B_ROWS
    return row_sub(rows[high], rows[low])


def word_constraints(
    gaps: tuple[int, ...],
) -> tuple[
    tuple[tuple[tuple[int, ...], int], ...],
    tuple[tuple[tuple[int, ...], int], ...],
]:
    fibres = parent.sum_fibres(gaps)
    equalities: list[tuple[tuple[int, ...], int]] = []
    inequalities: list[tuple[tuple[int, ...], int]] = []
    for fibre in fibres:
        base = cell_row(fibre[0])
        for cell in fibre[1:]:
            equalities.append((row_sub(cell_row(cell), base), 0))
    for left, right in zip(fibres, fibres[1:]):
        inequalities.append((row_sub(cell_row(right[0]), cell_row(left[0])), 1))
    for coordinate in range(5):
        row = tuple(int(index == coordinate) for index in range(5))
        inequalities.append((row, 1))
    return tuple(equalities), tuple(inequalities)


def floor_fraction(value: Fraction) -> int:
    return value.numerator // value.denominator


def ceil_fraction(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def one_parameter_integer_point(
    equalities: tuple[tuple[tuple[int, ...], int], ...],
    inequalities: tuple[tuple[tuple[int, ...], int], ...],
) -> tuple[int, ...] | None:
    """Solve an exact affine system whose equality space has dimension <=1."""
    matrix = [
        [Fraction(value) for value in row] + [Fraction(rhs)]
        for row, rhs in equalities
    ]
    pivot_row = 0
    pivots: list[int] = []
    for column in range(5):
        source = next(
            (index for index in range(pivot_row, len(matrix)) if matrix[index][column]),
            None,
        )
        if source is None:
            continue
        matrix[pivot_row], matrix[source] = matrix[source], matrix[pivot_row]
        pivot = matrix[pivot_row][column]
        matrix[pivot_row] = [value / pivot for value in matrix[pivot_row]]
        for index in range(len(matrix)):
            if index == pivot_row:
                continue
            multiplier = matrix[index][column]
            if multiplier:
                matrix[index] = [
                    value - multiplier * pivot_value
                    for value, pivot_value in zip(matrix[index], matrix[pivot_row])
                ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    for row in matrix:
        if not any(row[:5]) and row[5]:
            return None
    matrix = [row for row in matrix if any(row[:5])]
    free = tuple(column for column in range(5) if column not in pivots)
    require(len(free) <= 1, ("affine residual dimension", len(free)))

    alpha = [Fraction(0) for _ in range(5)]
    beta = [Fraction(0) for _ in range(5)]
    if free:
        alpha[free[0]] = Fraction(1)
    for row, pivot in zip(matrix, pivots):
        beta[pivot] = row[5]
        if free:
            alpha[pivot] = -row[free[0]]

    if not free:
        if any(value.denominator != 1 for value in beta):
            return None
        point = tuple(int(value) for value in beta)
        if all(sum(a * x for a, x in zip(row, point)) >= rhs for row, rhs in inequalities):
            return point
        return None

    lower: int | None = None
    upper: int | None = None
    for row, rhs in inequalities:
        slope = sum(Fraction(value) * coefficient for value, coefficient in zip(row, alpha))
        intercept = sum(Fraction(value) * constant for value, constant in zip(row, beta))
        target = Fraction(rhs) - intercept
        if slope > 0:
            bound = ceil_fraction(target / slope)
            lower = bound if lower is None else max(lower, bound)
        elif slope < 0:
            bound = floor_fraction(target / slope)
            upper = bound if upper is None else min(upper, bound)
        elif intercept < rhs:
            return None
    if lower is not None and upper is not None and lower > upper:
        return None

    period = 1
    for value in (*alpha, *beta):
        period = lcm(period, value.denominator)
    require(period >= 1, "positive congruence period")
    start = lower if lower is not None else (upper - period + 1 if upper is not None else 0)
    stop = start + period - 1
    if upper is not None:
        stop = min(stop, upper)
    for parameter in range(start, stop + 1):
        point_fraction = tuple(
            coefficient * parameter + constant
            for coefficient, constant in zip(alpha, beta)
        )
        if any(value.denominator != 1 for value in point_fraction):
            continue
        point = tuple(int(value) for value in point_fraction)
        if all(sum(a * x for a, x in zip(row, point)) >= rhs for row, rhs in inequalities):
            return point
    return None


def candidate_survives(
    gaps: tuple[int, ...],
    scalar_index: int,
    anchor: tuple[int, int],
    orientation: str,
) -> bool:
    a_support, b_support = parent.supports(gaps)
    fibres = parent.sum_fibres(gaps)
    i, j = anchor
    if orientation == "-+":
        p0, q0 = -2 - a_support[i], 1 - b_support[j]
    else:
        p0, q0 = 1 - a_support[i], -2 - b_support[j]
    p_weights = tuple(p0 + value for value in a_support)
    q_weights = tuple(q0 + value for value in b_support)
    if anchor not in fibres[scalar_index]:
        return False
    expected = (-2, 1) if orientation == "-+" else (1, -2)
    if (p_weights[i], q_weights[j]) != expected:
        return False
    if not all(p_weights[k] + q_weights[ell] + 1 == 0 for k, ell in fibres[scalar_index]):
        return False
    for index, fibre in enumerate(fibres):
        if index == scalar_index or len(fibre) != 1:
            continue
        k, ell = fibre[0]
        if singleton_verdict(p_weights[k], q_weights[ell]) == "FAIL":
            return False
    return True


def exact_word_gate_schemes(
    sample: tuple[int, ...],
) -> tuple[tuple[int, int, int, str, tuple[int, ...]], ...]:
    """Enumerate every anchor/orientation admitting some integer realization."""
    fibres = parent.sum_fibres(sample)
    singleton_cells = tuple(fibre[0] for fibre in fibres if len(fibre) == 1)
    base_equalities, base_inequalities = word_constraints(sample)
    schemes: dict[tuple[int, int, int, str], tuple[int, ...]] = {}

    for scalar_index, scalar_fibre in enumerate(fibres):
        if len(scalar_fibre) < 2:
            continue
        for i, j in scalar_fibre:
            for orientation in ("-+", "+-"):
                ordinary = list(base_inequalities)
                soft_options: list[
                    tuple[
                        tuple[tuple[tuple[int, ...], int], ...],
                        tuple[tuple[tuple[int, ...], int], ...],
                    ]
                ] = []
                impossible = False
                for k, ell in singleton_cells:
                    if orientation == "-+":
                        if k <= i and ell < j:
                            ordinary.append((distance_row("B", ell, j), 2))
                        elif k <= i and ell >= j:
                            impossible = True
                            break
                        elif k > i and ell >= j:
                            ordinary.append((distance_row("A", i, k), 3))
                        else:
                            a_distance = distance_row("A", i, k)
                            b_distance = distance_row("B", ell, j)
                            soft_options.append((
                                (((a_distance, 1),), ((b_distance, 2),)),
                                (((a_distance, 2), (b_distance, 1)), ()),
                            ))
                    else:
                        if k < i and ell <= j:
                            ordinary.append((distance_row("A", k, i), 2))
                        elif k >= i and ell <= j:
                            impossible = True
                            break
                        elif k >= i and ell > j:
                            ordinary.append((distance_row("B", j, ell), 3))
                        else:
                            a_distance = distance_row("A", k, i)
                            b_distance = distance_row("B", j, ell)
                            soft_options.append((
                                (((b_distance, 1),), ((a_distance, 2),)),
                                (((a_distance, 1), (b_distance, 2)), ()),
                            ))
                if impossible:
                    continue

                key = (scalar_index, i, j, orientation)
                witnesses: list[tuple[int, ...]] = []
                if not soft_options:
                    witness = tuple(3 * value for value in sample)
                    if all(
                        sum(a * x for a, x in zip(row, witness)) >= rhs
                        for row, rhs in ordinary
                    ) and candidate_survives(witness, scalar_index, (i, j), orientation):
                        witnesses.append(witness)
                else:
                    for choices in product((0, 1), repeat=len(soft_options)):
                        equalities = list(base_equalities)
                        inequalities = list(ordinary)
                        for options, choice in zip(soft_options, choices):
                            chosen_equalities, chosen_inequalities = options[choice]
                            equalities.extend(chosen_equalities)
                            inequalities.extend(chosen_inequalities)
                        witness = one_parameter_integer_point(
                            tuple(equalities), tuple(inequalities)
                        )
                        if witness is None:
                            continue
                        require(
                            parent.fibre_word(parent.sum_fibres(witness))
                            == parent.fibre_word(fibres),
                            "word constraint reconstruction",
                        )
                        require(
                            candidate_survives(witness, scalar_index, (i, j), orientation),
                            "symbolic gate witness",
                        )
                        witnesses.append(witness)
                if witnesses:
                    schemes[key] = min(witnesses, key=lambda value: (max(value), sum(value), value))

    return tuple(
        (*key, schemes[key]) for key in sorted(schemes)
    )


def format_counter(counter: Counter) -> str:
    return ",".join(f"{key}:{counter[key]}" for key in sorted(counter, key=str))


def format_placement(placement: Placement) -> str:
    eligible = ",".join(f"{i}{j}{orientation}" for i, j, orientation in placement.eligible)
    rows = ",".join(
        f"{i}{j}:{p}/{q}/{verdict}"
        for i, j, p, q, verdict in placement.singleton_rows
    ) or "-"
    components = "/".join(".".join(part) for part in placement.nontrivial_singleton_components) or "-"
    return (
        f"scalar={placement.scalar_index}:{fibre_label(placement.scalar_fibre)};"
        f"P={','.join(map(str, placement.p_weights))};"
        f"Q={','.join(map(str, placement.q_weights))};"
        f"eligible={eligible};singletons={rows};components={components}"
    )


def main() -> None:
    flats = parent.enumerate_flat_lattice()
    cones = parent.enumerate_positive_connected_cones(flats)
    word_sample: dict[str, tuple[int, ...]] = {}
    for cone in cones:
        for sample in cone.samples:
            word = parent.fibre_word(parent.sum_fibres(sample))
            require(word not in word_sample, "one chamber sample per word")
            word_sample[word] = sample
    require(len(word_sample) == 149, "parent word count")

    ordered_words = tuple(
        sorted(word_sample, key=lambda word: (len(parent.sum_fibres(word_sample[word])), word))
    )
    word_ids = {word: f"W{index:03d}" for index, word in enumerate(ordered_words, 1)}

    # Exact global branch: every soft cross-quadrant repair fixes an interval
    # distance to 1 or 2 and lowers the already at-most-two-dimensional word
    # cone to an affine line or point.  The one-parameter solver is exhaustive.
    exact_schemes = {
        word: exact_word_gate_schemes(sample) for word, sample in word_sample.items()
    }
    exact_words = {word for word, schemes in exact_schemes.items() if schemes}
    exact_word_counts: Counter[int] = Counter(
        len(parent.sum_fibres(word_sample[word])) for word in exact_words
    )
    exact_scheme_counts: Counter[int] = Counter()
    orientation_counts: Counter[str] = Counter()
    scalar_size_counts: Counter[int] = Counter()
    scalar_cut_counts: Counter[tuple[int, bool]] = Counter()
    scalar_exposed_counts: Counter[tuple[int, bool]] = Counter()
    words_with_scalar_cut: set[str] = set()
    words_with_unexposed_scalar: set[str] = set()
    for word, schemes in exact_schemes.items():
        sumset_size = len(parent.sum_fibres(word_sample[word]))
        exact_scheme_counts[sumset_size] += len(schemes)
        sample = word_sample[word]
        fibres = parent.sum_fibres(sample)
        for scalar_index, _, _, orientation, witness in schemes:
            orientation_counts[orientation] += 1
            scalar_size_counts[len(fibres[scalar_index])] += 1
            a_mask, b_mask = parent.projection_masks(
                fibres, frozenset((scalar_index,))
            )
            scalar_cut = (
                parent.A_COMPONENTS[a_mask] > 1 or parent.B_COMPONENTS[b_mask] > 1
            )
            scalar_cut_counts[sumset_size, scalar_cut] += 1
            scalar_exposed = parent.rectangle_exposed(
                witness, parent.sum_fibres(witness)[scalar_index]
            )
            scalar_exposed_counts[sumset_size, scalar_exposed] += 1
            if scalar_cut:
                words_with_scalar_cut.add(word)
            if not scalar_exposed:
                words_with_unexposed_scalar.add(word)
    exact_max_witness = max(
        max(scheme[-1]) for schemes in exact_schemes.values() for scheme in schemes
    )

    # The projective chamber sample is intentionally tested and rejected as a
    # sufficient affine classifier.  Common dilation is not a symmetry once
    # the scalar-arm weights -2 and 1 have been fixed.
    representative_words = {
        word for word, gaps in word_sample.items() if placements_for_sample(gaps)
    }
    generic_words = {
        word
        for word, gaps in word_sample.items()
        if placements_for_sample(tuple(3 * value for value in gaps))
    }

    # Independent hostile control on actual (not primitive-quotiented) gaps.
    bounded_vectors = 0
    bounded_surviving_vectors = 0
    bounded_placement_count = 0
    bounded_surviving_words: set[str] = set()
    bounded_minimal_witness: dict[str, tuple[int, ...]] = {}
    bounded_by_sumset: Counter[int] = Counter()
    for raw in product(range(1, HOSTILE_BOUND + 1), repeat=5):
        gaps = tuple(raw)
        a_mask, b_mask = parent.difference_graph_masks(gaps)
        if parent.A_COMPONENTS[a_mask] != 1 or parent.B_COMPONENTS[b_mask] != 1:
            continue
        bounded_vectors += 1
        placements = placements_for_sample(gaps)
        if not placements:
            continue
        bounded_surviving_vectors += 1
        bounded_placement_count += len(placements)
        word = parent.fibre_word(parent.sum_fibres(gaps))
        bounded_surviving_words.add(word)
        bounded_by_sumset[len(parent.sum_fibres(gaps))] += 1
        old = bounded_minimal_witness.get(word)
        key = (max(gaps), sum(gaps), gaps)
        if old is None or key < (max(old), sum(old), old):
            bounded_minimal_witness[word] = gaps

    exact_ids = tuple(sorted(word_ids[word] for word in exact_words))
    generic_ids = tuple(sorted(word_ids[word] for word in generic_words))
    representative_ids = tuple(sorted(word_ids[word] for word in representative_words))
    bounded_ids = tuple(sorted(word_ids[word] for word in bounded_surviving_words))
    require(exact_words == generic_words == bounded_surviving_words, "three-path word set")
    require(len(representative_words) == 18 and representative_words < exact_words,
            "projective-sample hostile")

    exact_ledger = tuple(
        (
            word_ids[word],
            len(parent.sum_fibres(word_sample[word])),
            word,
            tuple(
                (scalar_index, i, j, orientation, witness)
                for scalar_index, i, j, orientation, witness in exact_schemes[word]
            ),
        )
        for word in ordered_words
        if exact_schemes[word]
    )
    semantic = semantic_sha256(exact_ledger)
    bounded_semantic = semantic_sha256((
        bounded_vectors,
        bounded_surviving_vectors,
        bounded_placement_count,
        tuple(sorted(bounded_by_sumset.items())),
        bounded_ids,
        tuple(sorted((word_ids[word], gaps) for word, gaps in bounded_minimal_witness.items())),
    ))
    require(exact_ids == EXPECTED_EXACT_IDS, "exact word IDs")
    require(representative_ids == EXPECTED_REPRESENTATIVE_IDS, "projective hostile IDs")
    require(dict(exact_word_counts) == EXPECTED_WORD_COUNTS, "word counts")
    require(dict(exact_scheme_counts) == EXPECTED_SCHEME_COUNTS, "scheme counts")
    require(dict(orientation_counts) == EXPECTED_ORIENTATION_COUNTS, "orientation counts")
    require(dict(scalar_size_counts) == EXPECTED_SCALAR_SIZE_COUNTS, "scalar sizes")
    require(dict(scalar_cut_counts) == EXPECTED_SCALAR_CUT_COUNTS, "scalar cuts")
    require(
        dict(scalar_exposed_counts) == EXPECTED_SCALAR_EXPOSED_COUNTS,
        "scalar exposure",
    )
    require(
        tuple(sorted(word_ids[word] for word in words_with_unexposed_scalar))
        == EXPECTED_UNEXPOSED_WORD_IDS,
        "unexposed residual words",
    )
    require(exact_max_witness == EXPECTED_MAX_SCHEME_WITNESS_GAP, "scheme witness bound")
    require(
        bounded_vectors == EXPECTED_BOUNDED["connected"]
        and bounded_surviving_vectors == EXPECTED_BOUNDED["surviving"]
        and bounded_placement_count == EXPECTED_BOUNDED["placements"]
        and dict(bounded_by_sumset) == EXPECTED_BOUNDED["by_sumset"],
        "bounded hostile census",
    )
    require(semantic == EXPECTED_SEMANTIC_SHA256, "semantic digest")
    require(bounded_semantic == EXPECTED_BOUNDED_SEMANTIC_SHA256, "bounded digest")

    print("THM-3606 EXPONENT-TWO 3-BY-4 SCALAR/SINGLETON GATE ATLAS")
    print(f"parent_sha256_lf={EXPECTED_PARENT_SHA256}")
    print(f"parent_words=149;exact_surviving_words={len(exact_words)}")
    print("exact_word_counts=" + format_counter(exact_word_counts))
    print("exact_scheme_counts=" + format_counter(exact_scheme_counts))
    print("orientation_counts=" + format_counter(orientation_counts))
    print("scalar_fibre_size_counts=" + format_counter(scalar_size_counts))
    print("scalar_cut_scheme_counts=" + format_counter(scalar_cut_counts))
    print("scalar_exposed_scheme_counts=" + format_counter(scalar_exposed_counts))
    print(
        f"words_with_some_scalar_cut={len(words_with_scalar_cut)};ids="
        + ",".join(sorted(word_ids[word] for word in words_with_scalar_cut))
    )
    print(
        f"words_with_unexposed_scalar_scheme={len(words_with_unexposed_scalar)};ids="
        + ",".join(sorted(word_ids[word] for word in words_with_unexposed_scalar))
    )
    print(f"exact_max_scheme_witness_gap={exact_max_witness}")
    print("exact_word_ids=" + ",".join(exact_ids))
    print(
        f"projective_sample_words={len(representative_words)};ids="
        + ",".join(representative_ids)
    )
    print(f"generic_scaled_words={len(generic_words)};ids=" + ",".join(generic_ids))
    print(
        f"hostile_bound={HOSTILE_BOUND};connected_vectors={bounded_vectors};"
        f"surviving_vectors={bounded_surviving_vectors};"
        f"surviving_placements={bounded_placement_count};"
        f"surviving_words={len(bounded_surviving_words)}"
    )
    print("bounded_surviving_vectors_by_sumset=" + format_counter(bounded_by_sumset))
    print("three_path_word_sets_equal=1")
    print(f"semantic_sha256={semantic}")
    print(f"bounded_semantic_sha256={bounded_semantic}")
    print(f"script_sha256_lf={lf_sha256(Path(__file__).resolve())}")
    print("EXACT_GLOBAL_WORDS")
    for word_id, sumset_size, word, schemes in exact_ledger:
        best = min((scheme[-1] for scheme in schemes), key=lambda value: (max(value), sum(value), value))
        residual = tuple(
            scheme
            for scheme in schemes
            if not parent.rectangle_exposed(
                scheme[-1], parent.sum_fibres(scheme[-1])[scheme[0]]
            )
        )
        residual_anchor = "-"
        if residual:
            chosen = min(residual, key=lambda scheme: (max(scheme[-1]), sum(scheme[-1]), scheme))
            residual_anchor = (
                f"{chosen[0]}:{chosen[1]}{chosen[2]}{chosen[3]}@"
                + ",".join(map(str, chosen[-1]))
            )
        print(
            f"{word_id}=m:{sumset_size};schemes:{len(schemes)};"
            f"unexposed:{len(residual)};residual_anchor:{residual_anchor};"
            f"best_witness:{','.join(map(str,best))};"
            f"scheme_sha256:{semantic_sha256(schemes)};word:{word}"
        )
    print(
        "boundary=necessary scalar-arm and singleton-row gates only;"
        "no multi-address differential cancellation, arm nonalternation,"
        "Hermite-Pade regularity, Darboux pair, or Jacobian conclusion"
    )
    print("status=PROVED+FINITE-EXACT+VERIFIED-EXACT;global affine gate atlas")
    print("PASS")


if __name__ == "__main__":
    main()
