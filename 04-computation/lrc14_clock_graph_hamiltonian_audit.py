#!/usr/bin/env python3
"""Exact Hamiltonian audit of the seven-clock incidence skeleton.

This companion consumes the full guard-cospan incidence bank built by
``lrc14_guard_cospan_successor_private_clock_collapse.py``.  For a fixed
source residue ``s`` it regards the seven numerical clock labels as vertices
and puts a directed edge ``ell4 -> ell5`` precisely when the corresponding
13-by-13 common-x incidence matrix has positive support.  The matrix itself is
retained as an edge label.

This graph is intrinsic to the computed incidence bank, but it is **not** a
chronology graph.  Boolean multiplication below is only the formal operation
obtained after identifying consecutive printed digit labels.  No claim glues
the common-x witnesses in different cells or identifies ``h`` in one cell
with ``j`` in the next.  The Hamiltonian tests therefore classify the exact
matrix-labelled clock skeleton and provide hostile controls for broader
transition/no-go readings; they do not prove an LRC orbit statement.

The audit does not force the graph into a tournament.  Instead it records the
honest pair signature: bidirected, singly directed, and absent unordered clock
pairs.  This preserves the ties and missing comparisons that an orientation
gauge would discard.
"""

from collections import Counter
from hashlib import sha256
from itertools import permutations

import lrc14_guard_cospan_successor_private_clock_collapse as guard
import lrc14_successor_private_sharp_graph_clock_collapse as sharp


P = 13
Q7 = 7
SOURCES = tuple(range(1, P))
GENERIC_SOURCES = (1, 2, 3, 4, 5, 8, 9, 10, 11, 12)
SOURCE_CLASSES = (GENERIC_SOURCES, (6,), (7,))
REPRESENTATIVES = (1, 6, 7)

EXPECTED_EDGE_COUNTS_BY_STEP = {
    1: (6, 4, 6, 6, 4, 6),
    6: (2, 3, 4, 4, 2, 4),
    7: (4, 2, 4, 4, 3, 2),
}
EXPECTED_PAIR_SIGNATURES = {
    1: (13, 6, 2),
    6: (3, 13, 5),
    7: (3, 13, 5),
}
EXPECTED_OUT_DEGREES = {
    1: (4, 5, 4, 5, 5, 4, 5),
    6: (3, 4, 3, 4, 2, 1, 2),
    7: (3, 2, 1, 2, 4, 3, 4),
}
EXPECTED_IN_DEGREES = {
    1: (6, 3, 4, 6, 6, 4, 3),
    6: (4, 0, 0, 2, 6, 4, 3),
    7: (4, 3, 4, 6, 2, 0, 0),
}
EXPECTED_HAMILTONIAN_COUNTS = {
    1: (92, 876),
    6: (0, 0),
    7: (0, 0),
}
EXPECTED_LONGEST_SCALAR_PATH = {
    6: (6, 14),
    7: (6, 14),
}
EXPECTED_CYCLE_SUPPORT_HISTOGRAMS = {
    "safe": Counter({
        0: 628, 49: 1, 52: 2, 54: 2, 57: 2, 58: 1, 59: 2,
        60: 1, 61: 2, 62: 6, 63: 2, 64: 3, 65: 1, 66: 4,
        68: 6, 70: 14, 71: 1, 72: 5, 73: 1, 74: 6, 75: 9,
        76: 6, 77: 7, 78: 1, 80: 1, 83: 1, 86: 3, 87: 1,
        88: 1,
    }),
    "danger": Counter({0: 720}),
    "guard_free": Counter({
        0: 628, 87: 1, 88: 1, 92: 1, 95: 1, 96: 1, 101: 1,
        104: 1, 105: 2, 106: 4, 107: 3, 109: 1, 110: 2,
        111: 1, 112: 2, 113: 1, 114: 3, 115: 4, 116: 6,
        117: 1, 118: 5, 119: 9, 120: 6, 122: 1, 123: 2,
        124: 1, 125: 2, 126: 2, 127: 6, 128: 2, 129: 1,
        130: 9, 131: 1, 132: 1, 134: 1, 138: 2, 139: 1,
        140: 1, 141: 1, 143: 1,
    }),
}
EXPECTED_GENERIC_PATH_SUMMARIES = {
    "safe": (876, 39, 108),
    "danger": (0, 0, 0),
    "guard_free": (876, 73, 152),
}
EXPECTED_EXCEPTIONAL_SIX_PATH_SUMMARIES = {
    "safe": {6: (14, 44, 69), 7: (14, 46, 81)},
    "danger": {6: (0, 0, 0), 7: (0, 0, 0)},
    "guard_free": {6: (14, 73, 107), 7: (14, 69, 131)},
}
EXPECTED_DANGER_LONGEST_POSITIVE_PATH = {
    1: (3, 92, 2, 2),
    6: (3, 26, 2, 2),
    7: (3, 17, 2, 2),
}
EXPECTED_DANGER_UNION_POWER_SUPPORTS = (5, 2, 0, 0, 0)
EXPECTED_DANGER_UNION_ENTRIES = (
    (0, 11), (5, 1), (6, 0), (12, 0), (12, 1),
)
EXPECTED_EDGE_INCREMENT_SUPPORT_HISTOGRAMS = {
    "safe": Counter({9: 68, 10: 101, 11: 189}),
    "danger": Counter({1: 24, 3: 57, 4: 44, 5: 233}),
    "guard_free": Counter({10: 24, 11: 101, 12: 233}),
}
EXPECTED_EDGE_HULL_DEFECT_HISTOGRAMS = {
    "safe": Counter({
        101: 11, 102: 11, 103: 46, 109: 22, 110: 46,
        111: 22, 112: 11, 121: 90, 123: 55, 125: 44,
    }),
    "danger": Counter({12: 24, 36: 57, 48: 44, 60: 233}),
    "guard_free": Counter({
        110: 24, 120: 22, 121: 46, 122: 22, 123: 11,
        132: 66, 134: 55, 135: 11, 136: 55, 137: 46,
    }),
}
EXPECTED_CYCLE_MINIMUM_FIBRE_HISTOGRAMS = {
    "safe": Counter({2: 1, 3: 17, 4: 31, 5: 40, 6: 3}),
    "danger": Counter(),
    "guard_free": Counter({5: 2, 6: 3, 7: 14, 8: 29, 9: 28,
                            10: 15, 11: 1}),
}
EXPECTED_CYCLE_MAXIMUM_FIBRE_HISTOGRAMS = {
    "safe": Counter({5: 7, 6: 28, 7: 51, 8: 6}),
    "danger": Counter(),
    "guard_free": Counter({8: 5, 9: 12, 10: 63, 11: 12}),
}
EXPECTED_CYCLE_ALL_FIBRE_HISTOGRAMS = {
    "safe": Counter({2: 1, 3: 31, 4: 165, 5: 430,
                     6: 462, 7: 101, 8: 6}),
    "danger": Counter(),
    "guard_free": Counter({5: 3, 6: 11, 7: 60, 8: 202,
                           9: 452, 10: 395, 11: 73}),
}
EXPECTED_ORDINARY_STATE_PATH_SUMMARIES = {
    "safe": (
        92, (2, 21), (2, 14), (9, 21), (77, 236),
        Counter({
            0: 125, 1: 75, 2: 48, 3: 28, 4: 24, 5: 34, 6: 27,
            7: 78, 8: 117, 9: 154, 10: 170, 11: 172, 12: 144,
        }),
        0,
        "7fd6051b46ab288d1207f4aef2f86018d8bd58531b361d07cb1c2793080c0a2b",
    ),
    "danger": (
        0, (0, 0), (0, 0), (0, 0), (0, 0), Counter(), 0,
        "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d",
    ),
    "guard_free": (
        92, (14, 49), (14, 41), (22, 49), (240, 592),
        Counter({
            0: 99, 1: 89, 2: 96, 3: 92, 4: 103, 5: 82, 6: 89,
            7: 88, 8: 91, 9: 88, 10: 90, 11: 85, 12: 104,
        }),
        0,
        "72c9d2a57cb0147bdadeb1393abe117f446dcbebc6f692efa4e52e7ae60b68b8",
    ),
}
EXPECTED_DILATION_CYCLE_SUMMARIES = {
    "safe": (92, 45, 100, 92, 0),
    "danger": (0, 0, 0, 0, 0),
    "guard_free": (92, 82, 144, 92, 0),
}
EXPECTED_DILATION_PATH_SUMMARIES = {
    "safe": (876, 30, 106),
    "danger": (0, 0, 0),
    "guard_free": (876, 65, 154),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def bit_matrix(matrix):
    """Encode a Boolean matrix by one 13-bit mask per row."""
    return tuple(
        sum(int(value != 0) << column for column, value in enumerate(row))
        for row in matrix
    )


BIT_IDENTITY = tuple(1 << row for row in range(P))
INTEGER_IDENTITY = tuple(
    tuple(int(row == column) for column in range(P))
    for row in range(P)
)


def bit_product(left, right):
    """Boolean matrix product in row-bitset form."""
    result = []
    for middle_mask in left:
        row = 0
        while middle_mask:
            bit = middle_mask & -middle_mask
            row |= right[bit.bit_length() - 1]
            middle_mask -= bit
        result.append(row)
    return tuple(result)


def bit_support(matrix):
    return sum(row.bit_count() for row in matrix)


def bit_union(matrices):
    result = [0] * P
    for matrix in matrices:
        for row in range(P):
            result[row] |= matrix[row]
    return tuple(result)


def bit_entries(matrix):
    return tuple(
        (row, column)
        for row in range(P) for column in range(P)
        if matrix[row] & (1 << column)
    )


def integer_matrix(matrix):
    return tuple(tuple(int(value != 0) for value in row) for row in matrix)


def integer_product(left, right):
    """Ordinary product: entries count labelled intermediate state paths."""
    return tuple(
        tuple(
            sum(left[row][middle] * right[middle][column]
                for middle in range(P))
            for column in range(P)
        )
        for row in range(P)
    )


def freeze_source(bank, source):
    return tuple(
        bank[d, source, ell]
        for d in sharp.CLOCK_STEPS
        for ell in range(Q7)
    )


def source_classes(bank):
    classes = []
    unused = set(SOURCES)
    while unused:
        source = min(unused)
        signature = freeze_source(bank, source)
        equivalent = tuple(
            candidate for candidate in sorted(unused)
            if freeze_source(bank, candidate) == signature
        )
        classes.append(equivalent)
        unused.difference_update(equivalent)
    return tuple(classes)


def encode_bank(bank):
    return {
        key: bit_matrix(matrix)
        for key, matrix in bank.items()
    }


def edge_matrix(bank, source, left, right):
    require(left != right, "clock skeleton has no loop edges")
    return bank[(right - left) % Q7, source, left]


def edge_present(bank, source, left, right):
    return bool(bit_support(edge_matrix(bank, source, left, right)))


def formal_itinerary_product(bank, source, itinerary, close):
    result = BIT_IDENTITY
    pairs = list(zip(itinerary, itinerary[1:]))
    if close:
        pairs.append((itinerary[-1], itinerary[0]))
    for left, right in pairs:
        result = bit_product(result, edge_matrix(
            bank, source, left, right
        ))
    return result


def formal_integer_itinerary_product(bank, source, itinerary, close):
    """Count state sequences; unlike Boolean product, retain multiplicity."""
    result = INTEGER_IDENTITY
    pairs = list(zip(itinerary, itinerary[1:]))
    if close:
        pairs.append((itinerary[-1], itinerary[0]))
    for left, right in pairs:
        result = integer_product(result, edge_matrix(
            bank, source, left, right
        ))
    return result


def dilation_edge_matrix(bank, source, shallow, rail):
    """Edge label in the D(x)={13x}-compatible clock orientation.

    The stored bookkeeping edge points from the rail clock to the shallow
    clock.  Dilation sends that rail clock to the next shallow clock, so the
    chronological candidate reverses the stored arrow.
    """
    return edge_matrix(bank, source, rail, shallow)


def dilation_itinerary_product(bank, source, itinerary, close):
    """Formal Boolean product in the dilation-compatible orientation."""
    result = BIT_IDENTITY
    pairs = list(zip(itinerary, itinerary[1:]))
    if close:
        pairs.append((itinerary[-1], itinerary[0]))
    for shallow, rail in pairs:
        result = bit_product(
            result, dilation_edge_matrix(bank, source, shallow, rail)
        )
    return result


def dilation_scalar_itinerary_present(bank, source, itinerary, close):
    pairs = list(zip(itinerary, itinerary[1:]))
    if close:
        pairs.append((itinerary[-1], itinerary[0]))
    return all(edge_present(bank, source, rail, shallow)
               for shallow, rail in pairs)


def scalar_itinerary_present(bank, source, itinerary, close):
    pairs = list(zip(itinerary, itinerary[1:]))
    if close:
        pairs.append((itinerary[-1], itinerary[0]))
    return all(edge_present(bank, source, left, right)
               for left, right in pairs)


def pair_signature(bank, source):
    signature = Counter()
    for left in range(Q7):
        for right in range(left + 1, Q7):
            forward = edge_present(bank, source, left, right)
            backward = edge_present(bank, source, right, left)
            if forward and backward:
                signature["bidirected"] += 1
            elif forward or backward:
                signature["single"] += 1
            else:
                signature["absent"] += 1
    return (
        signature["bidirected"],
        signature["single"],
        signature["absent"],
    )


def degrees(bank, source):
    out_degrees = tuple(
        sum(edge_present(bank, source, left, right)
            for right in range(Q7) if right != left)
        for left in range(Q7)
    )
    in_degrees = tuple(
        sum(edge_present(bank, source, left, right)
            for left in range(Q7) if left != right)
        for right in range(Q7)
    )
    return out_degrees, in_degrees


def longest_scalar_paths(bank, source):
    for vertex_count in range(Q7, 0, -1):
        paths = tuple(
            itinerary
            for itinerary in permutations(range(Q7), vertex_count)
            if scalar_itinerary_present(bank, source, itinerary, False)
        )
        if paths:
            return vertex_count, paths
    raise RuntimeError("one-vertex scalar paths unexpectedly absent")


def longest_positive_paths(bank, source):
    for vertex_count in range(Q7, 0, -1):
        supports = tuple(
            bit_support(formal_itinerary_product(
                bank, source, itinerary, False
            ))
            for itinerary in permutations(range(Q7), vertex_count)
        )
        positives = tuple(value for value in supports if value)
        if positives:
            return (
                vertex_count, len(positives),
                min(positives), max(positives),
            )
    raise RuntimeError("identity relation unexpectedly empty")


def cyclic_diagonal_counts(matrix):
    return tuple(
        sum(bool(matrix[row] & (1 << ((row + shift) % P)))
            for row in range(P))
        for shift in range(P)
    )


def integer_cyclic_diagonal_counts(matrix):
    """Count all state paths ending at displacement ``column-row``."""
    return tuple(
        sum(matrix[row][(row + shift) % P] for row in range(P))
        for shift in range(P)
    )


def translation_profile(matrix):
    """Return pair support, increment support, and equivariant-hull defect.

    The minimal simultaneous-translation-equivariant relation containing R
    has all 13 translates of each occurring increment h-j.  Its support is
    therefore ``13*increment_support``.  The defect below vanishes exactly
    when R itself is translation equivariant.
    """
    multiplicities = cyclic_diagonal_counts(matrix)
    pair_support = sum(multiplicities)
    increment_support = sum(value > 0 for value in multiplicities)
    hull_defect = P * increment_support - pair_support
    return pair_support, increment_support, hull_defect, multiplicities


def counter_string(counter):
    return tuple(sorted(counter.items()))


def main():
    data = guard.rebuild_cospan_incidence_bank()
    raw_banks = data[0]
    banks = {
        sector: encode_bank(raw_banks[sector])
        for sector in guard.SECTORS
    }
    integer_banks = {
        sector: {
            key: integer_matrix(matrix)
            for key, matrix in raw_banks[sector].items()
        }
        for sector in guard.SECTORS
    }

    # The delayed guard labels change the 13-state edge relations but not the
    # scalar clock skeleton.  Exact matrix equality has the same three source
    # classes in each labelled sector.
    for sector in guard.SECTORS:
        require(
            source_classes(raw_banks[sector]) == SOURCE_CLASSES,
            f"{sector} exact source-matrix classes changed",
        )
    base = banks["safe"]
    for sector in guard.SECTORS:
        for source in SOURCES:
            for left in range(Q7):
                for right in range(Q7):
                    if left != right:
                        require(
                            edge_present(banks[sector], source, left, right)
                            == edge_present(base, source, left, right),
                            "guard sector changed the scalar clock skeleton",
                        )

    skeleton_summaries = {}
    normalized_cycles = tuple(
        (0,) + tail for tail in permutations(range(1, Q7))
    )
    hamiltonian_paths = tuple(permutations(range(Q7)))
    for source in REPRESENTATIVES:
        step_counts = tuple(
            sum(edge_present(base, source, ell, (ell + step) % Q7)
                for ell in range(Q7))
            for step in sharp.CLOCK_STEPS
        )
        require(
            step_counts == EXPECTED_EDGE_COUNTS_BY_STEP[source],
            f"source {source} edge-by-step census changed",
        )
        signature = pair_signature(base, source)
        require(
            signature == EXPECTED_PAIR_SIGNATURES[source],
            f"source {source} pair signature changed",
        )
        out_degrees, in_degrees = degrees(base, source)
        require(out_degrees == EXPECTED_OUT_DEGREES[source],
                f"source {source} out-degrees changed")
        require(in_degrees == EXPECTED_IN_DEGREES[source],
                f"source {source} in-degrees changed")
        cycle_count = sum(
            scalar_itinerary_present(base, source, itinerary, True)
            for itinerary in normalized_cycles
        )
        path_count = sum(
            scalar_itinerary_present(base, source, itinerary, False)
            for itinerary in hamiltonian_paths
        )
        require(
            (cycle_count, path_count)
            == EXPECTED_HAMILTONIAN_COUNTS[source],
            f"source {source} Hamiltonian census changed",
        )
        skeleton_summaries[source] = (
            step_counts, signature, out_degrees, in_degrees,
            cycle_count, path_count,
        )

    # The exceptional skeletons are clock reflections.  Each has two distinct
    # indegree-zero vertices, which is already enough to forbid a spanning
    # directed path: only one vertex can be first in such a path.
    require(
        all(
            edge_present(base, 6, left, right)
            == edge_present(base, 7, (-left) % Q7, (-right) % Q7)
            for left in range(Q7) for right in range(Q7) if left != right
        ),
        "exceptional scalar skeletons lost their clock reflection",
    )
    require(tuple(i for i, degree in enumerate(EXPECTED_IN_DEGREES[6])
                  if degree == 0) == (1, 2),
            "source 6 indegree-zero pair changed")
    require(tuple(i for i, degree in enumerate(EXPECTED_IN_DEGREES[7])
                  if degree == 0) == (5, 6),
            "source 7 indegree-zero pair changed")
    exceptional_scalar_paths = {}
    for source in (6, 7):
        vertex_count, paths = longest_scalar_paths(base, source)
        require(
            (vertex_count, len(paths))
            == EXPECTED_LONGEST_SCALAR_PATH[source],
            f"source {source} longest scalar paths changed",
        )
        exceptional_scalar_paths[source] = paths

    # Only sixteen of the 92 generic directed cycles retain their reversal.
    # Keeping this number records a genuine orientation asymmetry without
    # arbitrarily orienting bidirected or absent clock pairs.
    generic_scalar_cycles = {
        itinerary for itinerary in normalized_cycles
        if scalar_itinerary_present(base, 1, itinerary, True)
    }
    reverse_survivors = sum(
        (itinerary[0],) + tuple(reversed(itinerary[1:]))
        in generic_scalar_cycles
        for itinerary in generic_scalar_cycles
    )
    require(reverse_survivors == 16,
            "generic cycle-reversal census changed")
    constant_step_cycles = tuple(
        tuple((index * step) % Q7 for index in range(Q7))
        for step in sharp.CLOCK_STEPS
    )
    require(
        all(not scalar_itinerary_present(base, 1, itinerary, True)
            for itinerary in constant_step_cycles),
        "a generic constant-step clock cycle became admissible",
    )

    # THM-2642's convolution law requires every edge relation to be invariant
    # under simultaneous state translation (j,h)->(j+a,h+a).  Audit that
    # hypothesis directly.  Zero matrices are invariant only vacuously; every
    # positive incidence matrix has a strictly larger equivariant hull.
    edge_translation_summaries = {}
    for sector in guard.SECTORS:
        positive = 0
        zero = 0
        nonzero_equivariant = 0
        increment_histogram = Counter()
        defect_histogram = Counter()
        for matrix in banks[sector].values():
            pair_support, increment_support, defect, _ = (
                translation_profile(matrix)
            )
            if pair_support == 0:
                zero += 1
                require(defect == 0,
                        "zero relation has nonzero equivariant-hull defect")
            else:
                positive += 1
                nonzero_equivariant += int(defect == 0)
                increment_histogram[increment_support] += 1
                defect_histogram[defect] += 1
        require((zero, positive) == (146, 358),
                f"{sector} scalar edge census changed")
        require(nonzero_equivariant == 0,
                f"{sector} acquired a nonzero translation-equivariant edge")
        require(
            increment_histogram
            == EXPECTED_EDGE_INCREMENT_SUPPORT_HISTOGRAMS[sector],
            f"{sector} edge increment-support census changed",
        )
        require(
            defect_histogram == EXPECTED_EDGE_HULL_DEFECT_HISTOGRAMS[sector],
            f"{sector} equivariant-hull defect census changed",
        )
        edge_translation_summaries[sector] = (
            zero, positive, nonzero_equivariant,
            increment_histogram, defect_histogram,
        )

    sector_summaries = {}
    for sector in guard.SECTORS:
        bank = banks[sector]
        cycle_supports = []
        positive_cycle_diagonals = []
        positive_cycle_profiles = []
        positive_cycle_path_counts = []
        for itinerary in normalized_cycles:
            product = formal_itinerary_product(bank, 1, itinerary, True)
            support = bit_support(product)
            cycle_supports.append(support)
            if support:
                profile = translation_profile(product)
                positive_cycle_profiles.append(profile)
                positive_cycle_diagonals.append(profile[3])
                positive_cycle_path_counts.append(
                    integer_cyclic_diagonal_counts(
                        formal_integer_itinerary_product(
                            integer_banks[sector], 1, itinerary, True
                        )
                    )
                )
        cycle_histogram = Counter(cycle_supports)
        require(
            cycle_histogram == EXPECTED_CYCLE_SUPPORT_HISTOGRAMS[sector],
            f"{sector} normalized-cycle support histogram changed",
        )

        path_supports = tuple(
            bit_support(formal_itinerary_product(
                bank, 1, itinerary, False
            ))
            for itinerary in hamiltonian_paths
        )
        positive_paths = tuple(value for value in path_supports if value)
        path_summary = (
            len(positive_paths),
            min(positive_paths, default=0),
            max(positive_paths, default=0),
        )
        require(
            path_summary == EXPECTED_GENERIC_PATH_SUMMARIES[sector],
            f"{sector} generic Hamilton-path summary changed",
        )

        if sector in ("safe", "guard_free"):
            require(
                len(positive_cycle_diagonals) == len(generic_scalar_cycles)
                and all(all(count > 0 for count in diagonal_counts)
                        for diagonal_counts in positive_cycle_diagonals),
                f"{sector} positive cycle lost a state displacement",
            )
            require(
                all((support > 0)
                    == scalar_itinerary_present(
                        base, 1, itinerary, True
                    )
                    for support, itinerary
                    in zip(cycle_supports, normalized_cycles)),
                f"{sector} cycle composability left the scalar skeleton",
            )
            require(
                all((support > 0)
                    == scalar_itinerary_present(
                        base, 1, itinerary, False
                    )
                    for support, itinerary
                    in zip(path_supports, hamiltonian_paths)),
                f"{sector} path composability left the scalar skeleton",
            )

        positive_cycle_profiles = tuple(positive_cycle_profiles)
        positive_cycle_path_counts = tuple(positive_cycle_path_counts)
        full_pair_products = sum(
            pair_support == P * P
            for pair_support, _, _, _ in positive_cycle_profiles
        )
        full_difference_products = sum(
            increment_support == P
            for _, increment_support, _, _ in positive_cycle_profiles
        )
        minimum_fibre_histogram = Counter(
            min(multiplicities)
            for _, _, _, multiplicities in positive_cycle_profiles
        )
        maximum_fibre_histogram = Counter(
            max(multiplicities)
            for _, _, _, multiplicities in positive_cycle_profiles
        )
        all_fibre_histogram = Counter(
            value
            for _, _, _, multiplicities in positive_cycle_profiles
            for value in multiplicities
        )
        flat_path_counts = tuple(
            value for vector in positive_cycle_path_counts for value in vector
        )
        if flat_path_counts:
            ordinary_path_summary = (
                len(positive_cycle_path_counts),
                (min(flat_path_counts), max(flat_path_counts)),
                (min(min(vector) for vector in positive_cycle_path_counts),
                 max(min(vector) for vector in positive_cycle_path_counts)),
                (min(max(vector) for vector in positive_cycle_path_counts),
                 max(max(vector) for vector in positive_cycle_path_counts)),
                (min(sum(vector) for vector in positive_cycle_path_counts),
                 max(sum(vector) for vector in positive_cycle_path_counts)),
                Counter(value % P for value in flat_path_counts),
                sum(all(value % P == 0 for value in vector)
                    for vector in positive_cycle_path_counts),
                sha256(repr(positive_cycle_path_counts).encode()).hexdigest(),
            )
        else:
            ordinary_path_summary = (
                0, (0, 0), (0, 0), (0, 0), (0, 0), Counter(), 0,
                sha256(repr(positive_cycle_path_counts).encode()).hexdigest(),
            )
        require(
            ordinary_path_summary
            == EXPECTED_ORDINARY_STATE_PATH_SUMMARIES[sector],
            f"{sector} ordinary state-path multiplicities changed",
        )
        require(
            minimum_fibre_histogram
            == EXPECTED_CYCLE_MINIMUM_FIBRE_HISTOGRAMS[sector],
            f"{sector} minimum clutch-fibre census changed",
        )
        require(
            maximum_fibre_histogram
            == EXPECTED_CYCLE_MAXIMUM_FIBRE_HISTOGRAMS[sector],
            f"{sector} maximum clutch-fibre census changed",
        )
        require(
            all_fibre_histogram
            == EXPECTED_CYCLE_ALL_FIBRE_HISTOGRAMS[sector],
            f"{sector} all clutch-fibre census changed",
        )

        exceptional_six = {}
        for source in (6, 7):
            supports = tuple(
                bit_support(formal_itinerary_product(
                    bank, source, itinerary, False
                ))
                for itinerary in exceptional_scalar_paths[source]
            )
            positives = tuple(value for value in supports if value)
            summary = (
                len(positives),
                min(positives, default=0),
                max(positives, default=0),
            )
            require(
                summary
                == EXPECTED_EXCEPTIONAL_SIX_PATH_SUMMARIES[sector][source],
                f"{sector} source-{source} six-clock path summary changed",
            )
            exceptional_six[source] = summary

        sector_summaries[sector] = (
            cycle_histogram, path_summary, positive_cycle_diagonals,
            exceptional_six, full_pair_products, full_difference_products,
            minimum_fibre_histogram, maximum_fibre_histogram,
            all_fibre_histogram, ordinary_path_summary,
        )

    # The danger matrices have the same nonzero clock arrows as the other two
    # sectors, but their 13-state labels cannot compose along four distinct
    # clocks.  This is an internal-state obstruction, not a graph obstruction.
    danger_longest = {}
    for source in REPRESENTATIVES:
        summary = longest_positive_paths(banks["danger"], source)
        require(
            summary == EXPECTED_DANGER_LONGEST_POSITIVE_PATH[source],
            f"danger source-{source} longest positive path changed",
        )
        danger_longest[source] = summary

    # Stronger state-side explanation: even after forgetting *all* clock and
    # step labels, the danger union U has five arrows and U^3=0.  Thus the
    # danger obstruction survives this hostile quotient and is not caused by
    # the two exceptional clock skeletons.
    danger_union_entries = {}
    danger_union_powers = {}
    for source in REPRESENTATIVES:
        union = bit_union(
            banks["danger"][step, source, ell]
            for step in sharp.CLOCK_STEPS for ell in range(Q7)
        )
        power = BIT_IDENTITY
        supports = []
        for _ in range(1, 6):
            power = bit_product(power, union)
            supports.append(bit_support(power))
        require(
            tuple(supports) == EXPECTED_DANGER_UNION_POWER_SUPPORTS,
            f"danger source-{source} marginal union nilpotence changed",
        )
        danger_union_entries[source] = bit_entries(union)
        require(
            danger_union_entries[source] == EXPECTED_DANGER_UNION_ENTRIES,
            f"danger source-{source} marginal union relation changed",
        )
        danger_union_powers[source] = tuple(supports)

    # The canonical base-13 dilation D(x)={13x} satisfies j(Dx)=h(x) and
    # sends the present rail-clock phase to the next shallow-clock phase.
    # Therefore its candidate chronology reverses the stored clock arrows.
    # This remains a formal fixed-source matrix test: the code does not assert
    # that D transports the other carrier factors or the source shift.
    dilation_scalar_cycles = tuple(
        itinerary for itinerary in normalized_cycles
        if dilation_scalar_itinerary_present(base, 1, itinerary, True)
    )
    require(len(dilation_scalar_cycles) == 92,
            "dilation-reversed scalar Hamiltonian census changed")
    dilation_constant_step_zeros = {}
    dilation_cycle_summaries = {}
    dilation_path_summaries = {}
    for sector in guard.SECTORS:
        bank = banks[sector]
        zero_count = 0
        for source in SOURCES:
            for step in sharp.CLOCK_STEPS:
                itinerary = tuple((index * step) % Q7 for index in range(Q7))
                zero_count += int(not bit_support(
                    dilation_itinerary_product(
                        bank, source, itinerary, True
                    )
                ))
        require(zero_count == 72,
                f"{sector} dilation constant-step zero census changed")
        dilation_constant_step_zeros[sector] = zero_count

        cycle_products = tuple(
            dilation_itinerary_product(bank, 1, itinerary, True)
            for itinerary in normalized_cycles
        )
        positive_cycles = tuple(
            product for product in cycle_products if bit_support(product)
        )
        cycle_summary = (
            len(positive_cycles),
            min((bit_support(product) for product in positive_cycles),
                default=0),
            max((bit_support(product) for product in positive_cycles),
                default=0),
            sum(all(value > 0 for value in cyclic_diagonal_counts(product))
                for product in positive_cycles),
            sum(bit_support(product) == P * P
                for product in positive_cycles),
        )
        require(
            cycle_summary == EXPECTED_DILATION_CYCLE_SUMMARIES[sector],
            f"{sector} dilation-reversed cycle summary changed",
        )
        dilation_cycle_summaries[sector] = cycle_summary

        path_supports = tuple(
            bit_support(dilation_itinerary_product(
                bank, 1, itinerary, False
            ))
            for itinerary in hamiltonian_paths
        )
        positive_paths = tuple(value for value in path_supports if value)
        path_summary = (
            len(positive_paths),
            min(positive_paths, default=0),
            max(positive_paths, default=0),
        )
        require(
            path_summary == EXPECTED_DILATION_PATH_SUMMARIES[sector],
            f"{sector} dilation-reversed path summary changed",
        )
        dilation_path_summaries[sector] = path_summary

    print("LRC14 seven-clock matrix-labelled Hamiltonian audit")
    print("scope=full exact guard cospan; fixed-source common-x incidence bank")
    print("typing=clock graph and formal Boolean label products; not chronology")
    print(f"exact_source_matrix_classes={SOURCE_CLASSES}")
    print(
        "edge_counts_by_step_representative="
        + str(tuple((source, skeleton_summaries[source][0])
                    for source in REPRESENTATIVES))
    )
    print(
        "honest_pair_signatures_(bidirected,single,absent)="
        + str(tuple((source, skeleton_summaries[source][1])
                    for source in REPRESENTATIVES))
    )
    print(
        "out_in_degrees="
        + str(tuple((source, skeleton_summaries[source][2],
                     skeleton_summaries[source][3])
                    for source in REPRESENTATIVES))
    )
    print(
        "normalized_Hamilton_cycles_and_paths="
        + str(tuple((source, skeleton_summaries[source][4],
                     skeleton_summaries[source][5])
                    for source in REPRESENTATIVES))
    )
    print(
        "exceptional_zero_indegree_clocks=((6,(1,2)),(7,(5,6))) "
        "reflection=ell->-ell"
    )
    print("exceptional_longest_scalar_paths=((6,6,14),(7,6,14))")
    print(
        f"generic_cycle_reversal_survivors={reverse_survivors}/92 "
        "constant_step_cycles=0/6"
    )
    print(
        "edge_translation_(zero,positive,nonzero_equivariant)_by_sector="
        + str(tuple((sector,) + edge_translation_summaries[sector][:3]
                    for sector in guard.SECTORS))
    )
    print(
        "positive_edge_increment_support_hist_by_sector="
        + str(tuple((sector, counter_string(
            edge_translation_summaries[sector][3]
        )) for sector in guard.SECTORS))
    )
    print(
        "positive_edge_equivariant_hull_defect_hist_by_sector="
        + str(tuple((sector, counter_string(
            edge_translation_summaries[sector][4]
        )) for sector in guard.SECTORS))
    )
    print(
        "generic_cycle_support_hist_by_sector="
        + str(tuple((sector, counter_string(
            sector_summaries[sector][0]
        )) for sector in guard.SECTORS))
    )
    print(
        "generic_positive_Hamilton_path_(count,min,max)_by_sector="
        + str(tuple((sector, sector_summaries[sector][1])
                    for sector in guard.SECTORS))
    )
    print(
        "exceptional_six_clock_positive_(count,min,max)_by_sector="
        + str(tuple((sector, tuple(sorted(
            sector_summaries[sector][3].items()
        ))) for sector in guard.SECTORS))
    )
    print(
        "safe_and_guard_free_positive_cycles_with_all_13_state_displacements="
        "((safe,92/92),(guard_free,92/92))"
    )
    print(
        "generic_representative_positive_cycle_(full_pair,full_difference)="
        + str(tuple((sector, sector_summaries[sector][4],
                     sector_summaries[sector][5])
                    for sector in guard.SECTORS))
    )
    print(
        "all_source_positive_cycle_(count,full_pair,full_difference)="
        + str(tuple((
            sector,
            sum(count for support, count
                in sector_summaries[sector][0].items() if support)
            * len(GENERIC_SOURCES),
            sector_summaries[sector][4] * len(GENERIC_SOURCES),
            sector_summaries[sector][5] * len(GENERIC_SOURCES),
        ) for sector in guard.SECTORS))
    )
    print(
        "generic_representative_cycle_minimum_endpoint_pairs_per_clutch_hist="
        + str(tuple((sector, counter_string(
            sector_summaries[sector][6]
        )) for sector in guard.SECTORS))
    )
    print(
        "generic_representative_cycle_maximum_endpoint_pairs_per_clutch_hist="
        + str(tuple((sector, counter_string(
            sector_summaries[sector][7]
        )) for sector in guard.SECTORS))
    )
    print(
        "generic_representative_cycle_all_endpoint_pairs_per_clutch_hist="
        + str(tuple((sector, counter_string(
            sector_summaries[sector][8]
        )) for sector in guard.SECTORS))
    )
    print(
        "generic_representative_ordinary_state_path_summary_"
        "(cycles,all_clutch_range,cycle_min_range,cycle_max_range,total_range,"
        "residue_mod13_hist,all_clutches_divisible_cycles,sha256)="
        + str(tuple((
            sector,
            sector_summaries[sector][9][0],
            sector_summaries[sector][9][1],
            sector_summaries[sector][9][2],
            sector_summaries[sector][9][3],
            sector_summaries[sector][9][4],
            counter_string(sector_summaries[sector][9][5]),
            sector_summaries[sector][9][6],
            sector_summaries[sector][9][7],
        ) for sector in guard.SECTORS))
    )
    print(
        "danger_longest_positive_simple_clock_paths_(vertices,count,min,max)="
        + str(tuple((source, danger_longest[source])
                    for source in REPRESENTATIVES))
    )
    print(
        "danger_clock_step_union_entries="
        + str(tuple((source, danger_union_entries[source])
                    for source in REPRESENTATIVES))
    )
    print(
        "danger_clock_step_union_power_supports="
        + str(tuple((source, danger_union_powers[source])
                    for source in REPRESENTATIVES))
    )
    print(
        "dilation_reversed_constant_step_zero_candidates_by_sector="
        + str(tuple((sector, dilation_constant_step_zeros[sector])
                    for sector in guard.SECTORS))
    )
    print(
        "dilation_reversed_generic_cycle_"
        "(positive,min,max,all13displacements,fullpair)="
        + str(tuple((sector, dilation_cycle_summaries[sector])
                    for sector in guard.SECTORS))
    )
    print(
        "dilation_reversed_generic_path_(positive,min,max)="
        + str(tuple((sector, dilation_path_summaries[sector])
                    for sector in guard.SECTORS))
    )
    print(
        "dilation_boundary=fixed-source reversed matrix audit only; "
        "no D-pullback interval intersection or source transport is asserted"
    )
    print(
        "verdict=PASS: generic varying-clock skeletons are richly Hamiltonian; "
        "exceptional sources fail already by two zero-indegree clocks; danger "
        "fails by the source-independent state law U^3=0"
    )
    print(
        "THM2642_boundary=no positive edge is translation equivariant; full "
        "difference support in a formal cycle is not full pair support or a "
        "base-uniform convolution count; ordinary products separately count "
        "all intermediate state paths"
    )
    print(
        "semantics=no common-x witness gluing or physical clock transition is "
        "asserted; the union power deliberately forgets every clock label"
    )


if __name__ == "__main__":
    main()
