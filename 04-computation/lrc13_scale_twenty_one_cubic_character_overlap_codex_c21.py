#!/usr/bin/env python3
"""Exact quotient verifier for the scale-21 all-order overlap obstruction.

This script uses only the Python standard library.  It derives the six
owner-normalized sheet sets, checks the cubic-character obstruction in
Z[omega] with exact integer pairs, and exhausts the 12^6 independent unit
multipliers.  The quotient is the pure-overlap terminal isolated by THM-988;
it does not replay the preceding 77,810,217,408-context scalar reduction.
"""

from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path


P = 13
C = 21
RATIOS = (1, 3, 4, 9, 10, 12)
OPPOSITE_CLASSES = ((1, 12), (3, 10), (4, 9))
UNITS = tuple(value for value in range(C) if gcd(value, C) == 1)
FULL = (1 << C) - 1
Q = frozenset(RATIOS)
NQ = frozenset(2 * value % P for value in Q)


def require(condition, message):
    """Optimization-stable proof obligation."""
    if not condition:
        raise RuntimeError(message)


def sheet_mask(values):
    return sum(1 << value for value in values)


def derive_base_set(ratio):
    """Strict-left/closed-right safe representatives, reduced modulo 21."""
    return frozenset(
        value % C
        for value in range(-C + 1, C + 1)
        if (value - C * ratio) % P == 0
    )


BASE = {ratio: derive_base_set(ratio) for ratio in RATIOS}
EXPECTED_BASE = {
    1: frozenset((0, 3, 8, 16)),
    3: frozenset((6, 11, 19)),
    4: frozenset((1, 6, 14, 19)),
    9: frozenset((2, 7, 15, 20)),
    10: frozenset((2, 10, 15)),
    12: frozenset((5, 13, 18)),
}

UNIT_MASK = sheet_mask(UNITS)
GCD3_RESIDUES = frozenset(value for value in range(C) if gcd(value, C) == 3)
GCD7_RESIDUES = frozenset(value for value in range(C) if gcd(value, C) == 7)
GCD3_MASK = sheet_mask(GCD3_RESIDUES)
GCD7_MASK = sheet_mask(GCD7_RESIDUES)

MASK = {
    ratio: {
        unit: sheet_mask(unit * value % C for value in BASE[ratio])
        for unit in UNITS
    }
    for ratio in RATIOS
}


# Exact Eisenstein arithmetic: a + b*omega, with omega^2 + omega + 1 = 0.
ONE = (1, 0)
OMEGA = (0, 1)
OMEGA2 = (-1, -1)
ROOTS = (ONE, OMEGA, OMEGA2)
ZERO = (0, 0)


def eadd(left, right):
    return (left[0] + right[0], left[1] + right[1])


def eneg(value):
    return (-value[0], -value[1])


def emul(left, right):
    a, b = left
    c, d = right
    return (a * c - b * d, a * d + b * c - b * d)


def esum(values):
    answer = ZERO
    for value in values:
        answer = eadd(answer, value)
    return answer


def enorm(value):
    a, b = value
    return a * a - a * b + b * b


LOG7 = {}
power = 1
for exponent in range(6):
    LOG7[power] = exponent
    power = 3 * power % 7


def chi_exponent(value):
    """Cubic-character exponent in Z/3 for a nonzero residue modulo 7."""
    residue = value % 7
    require(residue in LOG7, "cubic character received zero modulo seven")
    return LOG7[residue] % 3


def chi(value):
    return ROOTS[chi_exponent(value)]


def cubic_coefficient(ratio):
    """Unit-pair character sum divided by the gcd-3 mark character."""
    unit_pair = tuple(
        value for value in BASE[ratio] if gcd(value, C) == 1
    )
    gcd3_point = next(
        value for value in BASE[ratio] if gcd(value, C) == 3
    )
    mark = gcd3_point // 3
    inverse_mark_character = ROOTS[-chi_exponent(mark) % 3]
    return emul(esum(chi(value) for value in unit_pair), inverse_mark_character)


COEFFICIENT = {ratio: cubic_coefficient(ratio) for ratio in RATIOS}
EXPECTED_COEFFICIENT = {
    1: eneg(OMEGA),
    12: eneg(OMEGA),
    3: eneg(OMEGA),
    10: eneg(OMEGA),
    4: eneg(OMEGA2),
    9: eneg(OMEGA2),
}


def tournament_fingerprint(coefficients):
    """Lexicographic coefficient tournament with coordinate-order tie path."""
    count = len(coefficients)
    out = [set() for _ in range(count)]
    ties = 0
    for left, right in combinations(range(count), 2):
        if coefficients[left] == coefficients[right]:
            ties += 1
            winner = left
        else:
            winner = max((left, right), key=lambda index: coefficients[index])
        loser = left + right - winner
        out[winner].add(loser)

    scores = tuple(len(edges) for edges in out)
    triangles = 0
    for first, second, third in combinations(range(count), 3):
        if (
            second in out[first]
            and third in out[second]
            and first in out[third]
        ) or (
            third in out[first]
            and second in out[third]
            and first in out[second]
        ):
            triangles += 1

    path_count = 0
    for ordering in permutations(range(count)):
        if all(ordering[index + 1] in out[ordering[index]]
               for index in range(count - 1)):
            path_count += 1

    reach = [set(edges) | {vertex} for vertex, edges in enumerate(out)]
    for middle in range(count):
        for source in range(count):
            if middle in reach[source]:
                reach[source] |= reach[middle]
    unused = set(range(count))
    sccs = []
    while unused:
        root = min(unused)
        component = {
            vertex for vertex in unused
            if vertex in reach[root] and root in reach[vertex]
        }
        sccs.append(len(component))
        unused -= component
    return ties, scores, triangles, tuple(sorted(sccs)), path_count


def main():
    require(len(UNITS) == 12, "unit group cardinality mismatch")
    require(Q == frozenset((1, 3, 4, 9, 10, 12)),
            "quadratic-residue support mismatch")
    require(NQ == frozenset((2, 5, 6, 7, 8, 11)),
            "nonresidue support mismatch")
    require(BASE == EXPECTED_BASE, "derived normalized base sets changed")
    require(
        all(BASE[-ratio % P] == frozenset(-value % C for value in BASE[ratio])
            for ratio in (3, 4)),
        "opposite-ratio symmetry changed",
    )
    require(
        BASE[12] == frozenset(-value % C for value in BASE[1] if value),
        "strict endpoint asymmetry at the plus/minus-one pair changed",
    )
    require(sum(map(len, BASE.values())) == C,
            "all-order scalar capacity is not exactly twenty-one")

    stratum_counts = Counter(
        gcd(value, C) for base in BASE.values() for value in base
    )
    require(
        stratum_counts == Counter({1: 12, 3: 6, 7: 2, 21: 1}),
        "multiplicative-stratum capacities changed",
    )

    require(COEFFICIENT == EXPECTED_COEFFICIENT,
            "cubic-character coefficient table changed")
    for ratio in RATIOS:
        base_mark = next(
            value // 3 for value in BASE[ratio] if gcd(value, C) == 3
        )
        for unit in UNITS:
            selected_unit_sum = esum(
                chi(unit * value)
                for value in BASE[ratio]
                if gcd(value, C) == 1
            )
            selected_mark_character = chi(unit * base_mark)
            require(
                selected_unit_sum
                == emul(COEFFICIENT[ratio], selected_mark_character),
                "multiplier-equivariant cubic identity failed",
            )

    omega_difference = eadd(OMEGA, eneg(OMEGA2))
    require(enorm(omega_difference) == 3,
            "cubic cancellation multiplier lost norm three")
    require(
        all(eadd(left, right) != ZERO for left in ROOTS for right in ROOTS),
        "two cubic roots unexpectedly cancel",
    )
    pair_sums = tuple({eadd(left, right) for left in ROOTS for right in ROOTS})
    contradictory_character_triples = []
    for first, second, third in product(pair_sums, repeat=3):
        mark_equation = esum((first, second, third))
        unit_equation = esum((
            emul(eneg(OMEGA), first),
            emul(eneg(OMEGA), second),
            emul(eneg(OMEGA2), third),
        ))
        if mark_equation == ZERO and unit_equation == ZERO:
            contradictory_character_triples.append((first, second, third))
    require(not contradictory_character_triples,
            "cubic partition equations admitted a solution")

    union_size_histogram = Counter()
    maximum_union = -1
    maximum_assignments = 0
    maximum_unions = set()
    maximum_missing_histogram = Counter()
    maximum_stratum_histogram = Counter()
    near_digest = sha256()
    unit_partitions = 0
    gcd3_partitions = 0
    gcd7_partitions = 0
    unit_and_gcd3_partitions = 0
    full_covers = 0

    for multipliers in product(UNITS, repeat=len(RATIOS)):
        union = 0
        selected_masks = []
        for ratio, unit in zip(RATIOS, multipliers):
            selected = MASK[ratio][unit]
            selected_masks.append(selected)
            union |= selected

        size = union.bit_count()
        union_size_histogram[size] += 1
        maximum_union = max(maximum_union, size)
        unit_partition = union & UNIT_MASK == UNIT_MASK
        gcd3_partition = union & GCD3_MASK == GCD3_MASK
        gcd7_partition = union & GCD7_MASK == GCD7_MASK
        unit_partitions += unit_partition
        gcd3_partitions += gcd3_partition
        gcd7_partitions += gcd7_partition
        unit_and_gcd3_partitions += unit_partition and gcd3_partition
        full_covers += union == FULL

        if size != 20:
            continue
        maximum_assignments += 1
        maximum_unions.add(union)
        missing = next(
            sheet for sheet in range(C) if not (union >> sheet) & 1
        )
        multiplicity = [0] * C
        for selected in selected_masks:
            for sheet in range(C):
                multiplicity[sheet] += (selected >> sheet) & 1
        duplicated = [
            sheet for sheet, count in enumerate(multiplicity) if count == 2
        ]
        require(
            len(duplicated) == 1
            and max(multiplicity) == 2
            and sum(max(0, count - 1) for count in multiplicity) == 1,
            "twenty-sheet union is not a single-defect near partition",
        )
        duplicate = duplicated[0]
        missing_gcd = gcd(missing, C)
        duplicate_gcd = gcd(duplicate, C)
        maximum_missing_histogram[missing] += 1
        maximum_stratum_histogram[(missing_gcd, duplicate_gcd)] += 1
        if missing_gcd == 1:
            require(duplicate == 4 * missing % C,
                    "unit near-cover defect is not x -> 4x")
        elif missing_gcd == 3:
            require(duplicate in (3 * missing % C, 4 * missing % C),
                    "gcd-three near-cover defect left its cubic pair")
        else:
            raise RuntimeError("near cover missed zero or a gcd-seven sheet")
        near_digest.update(bytes(multipliers + (missing, duplicate)))

    expected_maximum_unions = {
        FULL ^ (1 << sheet) for sheet in range(C) if sheet % 7 != 0
    }
    require(maximum_union == 20 and full_covers == 0,
            "all-order overlap maximum changed")
    require(
        union_size_histogram
        == Counter({
            7: 24,
            8: 192,
            9: 2_088,
            10: 12_648,
            11: 62_760,
            12: 205_200,
            13: 476_856,
            14: 730_632,
            15: 753_840,
            16: 492_096,
            17: 197_088,
            18: 45_984,
            19: 5_904,
            20: 672,
        }),
        "complete union-size histogram changed",
    )
    require(maximum_assignments == 672,
            "twenty-sheet multiplier-tuple count changed")
    require(maximum_unions == expected_maximum_unions,
            "twenty-sheet missing-sheet classification changed")
    require(
        maximum_missing_histogram
        == Counter({
            sheet: 24 if gcd(sheet, C) == 1 else 64
            for sheet in range(C) if sheet % 7 != 0
        }),
        "near-cover missing-sheet multiplicities changed",
    )
    require(
        maximum_stratum_histogram == Counter({(1, 1): 288, (3, 3): 384}),
        "near-cover defect-stratum census changed",
    )
    require(
        unit_partitions == 960
        and gcd3_partitions == 46_080
        and gcd7_partitions == 1_492_992
        and unit_and_gcd3_partitions == 0,
        "unit-edge partition/cubic obstruction census changed",
    )
    require(
        near_digest.hexdigest()
        == "7d0016777dbca7704602b3775ec76a2cc4c9bb9b39277e87fc7e8795e9842d98",
        "near-cover tuple/defect digest changed",
    )

    explicit = (2, 17, 2, 17, 2, 2)
    explicit_union = 0
    explicit_sets = []
    for ratio, unit in zip(RATIOS, explicit):
        chosen = frozenset(unit * value % C for value in BASE[ratio])
        explicit_sets.append(tuple(sorted(chosen)))
        explicit_union |= sheet_mask(chosen)
    require(explicit_union == FULL ^ (1 << 1),
            "explicit twenty-sheet witness changed")

    quotient_coefficients = tuple(
        COEFFICIENT[pair[0]] for pair in OPPOSITE_CLASSES
    )
    tournament = tournament_fingerprint(quotient_coefficients)
    require(tournament == (1, (1, 0, 2), 0, (1, 1, 1), 1),
            "ratio-pair coefficient tournament changed")

    source_hash = sha256(Path(__file__).read_bytes()).hexdigest()
    max_union_digest = sha256()
    for union in sorted(maximum_unions):
        max_union_digest.update(union.to_bytes(3, "little"))
    require(
        max_union_digest.hexdigest()
        == "9787273584e90fde50c2aa90a320b0a149d6d92ef8ea96e0c35e454698b0113a",
        "maximum-union bank digest changed",
    )

    print("scale-twenty-one all-order cubic-character overlap verifier")
    print(f"source SHA256 {source_hash}")
    print(f"supports Q={tuple(sorted(Q))} 2Q={tuple(sorted(NQ))}")
    print("normalized base sets " + " ".join(
        f"A{ratio}={tuple(sorted(BASE[ratio]))}" for ratio in RATIOS
    ))
    print("capacity strata zero:1 gcd7:2 gcd3:6 units:12 total:21")
    print("exact cubic coefficients " + " ".join(
        f"+-{pair[0]}:{COEFFICIENT[pair[0]]}"
        for pair in OPPOSITE_CLASSES
    ))
    print(
        "cubic contradiction exact: mark equation S1+S3+S4=0; "
        "unit equation -omega(S1+S3)-omega^2*S4=0; "
        f"norm(omega-omega^2)={enorm(omega_difference)}; "
        "no two cubic roots sum to zero"
    )
    print(
        f"multiplier tuples {len(UNITS) ** len(RATIOS)}; "
        f"union-size histogram "
        + " ".join(f"{size}:{count}"
                   for size, count in sorted(union_size_histogram.items()))
    )
    print(
        f"maximum union {maximum_union}; maximizing tuples "
        f"{maximum_assignments}; distinct maximum unions {len(maximum_unions)}; "
        f"full covers {full_covers}"
    )
    print(
        "maximum unions exactly Z/21 minus t for 7 not dividing t; "
        f"union-bank SHA256 {max_union_digest.hexdigest()}; "
        f"tuple-defect SHA256 {near_digest.hexdigest()}"
    )
    print(
        "near-cover defects unit->unit 288 (missing x, duplicate 4x); "
        "gcd3->gcd3 384 (missing x, duplicate in {3x,4x}); "
        "sheets 0,7,14 never missed"
    )
    print(
        f"stratum partitions units:{unit_partitions} "
        f"gcd3:{gcd3_partitions} gcd7:{gcd7_partitions} "
        f"units-and-gcd3:{unit_and_gcd3_partitions}"
    )
    print(f"explicit multipliers {explicit}; sets {tuple(explicit_sets)}; missing 1")
    print(
        "tournament vertices opposite-ratio classes +-1,+-3,+-4; "
        "pair observable exact Eisenstein coefficient; lexicographic switch; "
        "coordinate-order tie Hamiltonian path"
    )
    print(
        f"tournament coefficients {quotient_coefficients}; ties {tournament[0]}; "
        f"scores {tournament[1]}; triangles {tournament[2]}; "
        f"SCC sizes {tournament[3]}; Hamiltonian paths {tournament[4]}"
    )
    print(
        "preserved audit: CRT gcd strata plus six colored unit edges preserve "
        "the terminal full-cover predicate and cubic sums"
    )
    print(
        "lost audit: the completed tournament preserves coefficient order only; "
        "it discards additive character cancellation, sheet incidence, and "
        "multiplier compatibility"
    )
    print(
        "challenged vertices: owners, providers, individual residues, gaps, "
        "boundaries, wall events, and cover arcs lose coupling; the faithful "
        "vertices here are multiplicative strata with colored ratio-pair edges"
    )


if __name__ == "__main__":
    main()
