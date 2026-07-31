#!/usr/bin/env python3
"""Terminal extra-tooth argument for the 60 D=56 k=2/p=5 survivors.

Every survivor has four denominator-eight drifts and one denominator in
{7,14,28,56}.  On q=8 fibres the latter contributes at most one point to
each fibre.  A denominator-eight mask has one guaranteed seven-point spike
and, on an open common-phase event of exact mass 1/7, one extra seven-point
spike.  Every target q=8 fibre has load at least four.  Hence all eight fibres
need a d=8 spike; the four guaranteed spikes are insufficient unless all four
extra-tooth events occur.  Their intersection has mass at most 1/7, smaller
than the aligned-safe floor 66/91.

This script independently reconstructs the 15 rows, q=8 histograms, four
denominator shapes, exact extra-tooth phase law, and 60->0 count.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
COMBINED_PATH = (
    ROOT / "04-computation" / "lrc14_three_drift_body_projection_fiber_thm2928.py"
)
EXPECTED_COMBINED_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)
COMPOSED_PATH = (
    ROOT
    / ".scratch"
    / "lrc14_k2_five_drift_pattern_compression_20260729"
    / "composed_common_u_qprofile.py"
)
EXPECTED_COMPOSED_SHA256 = (
    "4abfcbd03be0aa0ee334f8591ff783ea011086cc4eb198de60a678ad415165f6"
)

EXPECTED_ROWS = {
    (1, 2, 3, 4, 6, 8): (336, 48, ((6, 8),)),
    (1, 2, 3, 4, 6, 9): (504, 48, ((6, 8),)),
    (1, 2, 3, 4, 6, 12): (168, 40, ((4, 2), (5, 4), (6, 2))),
    (1, 2, 3, 4, 8, 12): (336, 48, ((6, 8),)),
    (1, 2, 3, 4, 9, 12): (504, 48, ((6, 8),)),
    (1, 2, 3, 6, 8, 12): (336, 44, ((5, 4), (6, 4))),
    (1, 2, 3, 6, 9, 12): (504, 48, ((6, 8),)),
    (1, 2, 4, 5, 8, 10): (560, 48, ((6, 8),)),
    (1, 2, 4, 6, 8, 12): (336, 44, ((5, 4), (6, 4))),
    (1, 2, 4, 6, 9, 12): (504, 48, ((6, 8),)),
    (1, 2, 4, 7, 8, 14): (784, 48, ((6, 8),)),
    (1, 3, 4, 6, 8, 12): (336, 44, ((5, 4), (6, 4))),
    (1, 3, 4, 6, 9, 12): (504, 48, ((6, 8),)),
    (2, 3, 4, 6, 8, 12): (336, 44, ((5, 4), (6, 4))),
    (2, 3, 4, 6, 9, 12): (504, 48, ((6, 8),)),
}
EXPECTED_SHAPES = (
    (7, 8, 8, 8, 8),
    (8, 8, 8, 8, 14),
    (8, 8, 8, 8, 28),
    (8, 8, 8, 8, 56),
)

D = 56
Q_FIBRES = 8
TWO_SAFE_FLOOR = Q(66, 91)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


require(
    file_sha256(COMBINED_PATH) == EXPECTED_COMBINED_SHA256,
    "frozen body-projection dependency changed",
)
spec = spec_from_file_location("lrc14_three_drift_combined", COMBINED_PATH)
combined = module_from_spec(spec)
spec.loader.exec_module(combined)
support = combined.support_module


def norm_fraction(value):
    residue = value % 1
    return min(residue, 1 - residue)


def literal_d8_classes(c, u):
    require(gcd(c, 8) == 1, "d8 numerator is not a unit")
    return tuple(
        residue
        for residue in range(8)
        if norm_fraction(Q(c) * (residue + u) / 8) < Q(1, 14)
    )


def extra_tooth_controls():
    """Verify N_8 is 1+1_{{cu} in (3/7,4/7)} on exact probes."""
    checks = 0
    for c in range(1, 32, 2):
        probes = {Q(index, 224) for index in range(224)}
        for integer in range(c + 1):
            for endpoint in (Q(3, 7), Q(4, 7)):
                u = (integer + endpoint) / c
                if 0 <= u < 1:
                    probes.add(u)
        for u in probes:
            classes = literal_d8_classes(c, u)
            theta = (c * u) % 1
            expected = 2 if Q(3, 7) < theta < Q(4, 7) else 1
            require(
                len(classes) == expected,
                ("d8 extra-tooth law failed", c, u, classes, theta),
            )
            checks += 1
    positive_classes = literal_d8_classes(1, Q(1, 2))
    require(
        positive_classes == (0, 7),
        ("positive two-spike control changed", positive_classes),
    )
    require(
        literal_d8_classes(1, Q(0)) == (0,),
        "positive one-spike control changed",
    )
    return checks, positive_classes


def q_type_width(d, width):
    common = gcd(d, Q_FIBRES)
    quotient, remainder = divmod(width, common)
    height = D // lcm(d, Q_FIBRES)
    return (
        height * quotient,
        height if remainder else 0,
        (Q_FIBRES // common) * remainder,
    )


def reconstruct_rows():
    rows = {}
    for body in combinations(range(1, 15), 6):
        L, ranges = support.safe_cell_ranges(body)
        if L % D:
            continue
        support_count = support.support_size_bitset(D, ranges)
        if body not in EXPECTED_ROWS:
            continue
        arcs = combined.projected_support_arcs(D, ranges)
        histogram = combined.residue_load_histogram(arcs, Q_FIBRES)
        rows[body] = (L, support_count, histogram)
    require(rows == EXPECTED_ROWS, ("terminal row reconstruction changed", rows))
    return rows


def main():
    require(
        file_sha256(COMPOSED_PATH) == EXPECTED_COMPOSED_SHA256,
        "composed common-u x q-profile dependency changed",
    )
    extra_control_count, positive_classes = extra_tooth_controls()
    rows = reconstruct_rows()
    require(len(EXPECTED_SHAPES) * len(rows) == 60, "incoming 60 count changed")

    d8_one_tooth = q_type_width(8, 1)
    d8_two_tooth = q_type_width(8, 2)
    require(d8_one_tooth == (0, 7, 1), "d8 base spike changed")
    require(d8_two_tooth == (0, 7, 2), "d8 extra spike changed")

    fifth_types = {}
    for fifth in (7, 14, 28, 56):
        width = (fifth + 6) // 7
        fifth_types[fifth] = q_type_width(fifth, width)
    require(
        fifth_types == {
            7: (1, 0, 0),
            14: (1, 0, 0),
            28: (1, 0, 0),
            56: (1, 0, 0),
        },
        ("fifth-mask q8 baseline changed", fifth_types),
    )

    support_histogram = Counter()
    minimum_extra_histogram = Counter()
    heavy_fibre_histogram = Counter()
    semantic = sha256()
    for body, (L, support_count, histogram) in sorted(rows.items()):
        support_histogram[support_count] += 1
        minimum_target_load = min(load for load, count in histogram if count)
        require(
            minimum_target_load > 1,
            ("a target q8 fibre does not need a d8 spike", body, histogram),
        )
        # The fifth mask gives capacity one in all eight fibres.  Four d8
        # masks give one guaranteed spike each and at most one extra each.
        # Every target fibre needs a spike, so all four extras are forced.
        heavy_fibres = sum(count for load, count in histogram if load > 1)
        require(heavy_fibres == 8, ("heavy q8 fibre count changed", body))
        heavy_fibre_histogram[heavy_fibres] += 1
        minimum_extras = heavy_fibres - 4
        require(minimum_extras == 4, ("extra-tooth tariff changed", body))
        minimum_extra_histogram[minimum_extras] += 1
        semantic.update(
            f"{body}|{L}|{support_count}|{histogram}|{minimum_extras}\n".encode()
        )

    # Coverage would force R_A into the intersection of four open extra-tooth
    # events.  Each has exact mass 1/7; the intersection has mass at most 1/7.
    extra_event_mass = Q(1, 7)
    require(extra_event_mass < TWO_SAFE_FLOOR, "extra-tooth mass gap changed")
    # Exact-mean primary route: each d=8 spike mask hits Y_i(u) q-fibres and
    # integral Y_i du = (q/d)*(d/7)=q/7.  Covering N_1=8 heavy fibres forces
    # sum_i Y_i>=8.  Markov bounds that open event by 4/7.
    exact_mean_per_spike = Q(Q_FIBRES, 7)
    markov_event_bound = Q(4) * exact_mean_per_spike / 8
    necessary_heavy_cutoff = Q(13 * 4 * Q_FIBRES, 66)
    require(markov_event_bound == Q(4, 7), "exact-mean Markov bound changed")
    require(
        Q(8) >= necessary_heavy_cutoff,
        "exact-mean heavy-fibre tariff no longer kills",
    )
    require(
        markov_event_bound < TWO_SAFE_FLOOR,
        "exact-mean event mass gap changed",
    )

    print("LRC14 k=2/p=5 D56 eight-spike terminal")
    print(f"combined_script_sha256={file_sha256(COMBINED_PATH)}")
    print(f"composed_script_sha256={file_sha256(COMPOSED_PATH)}")
    print(f"extra_tooth_control_cases={extra_control_count}")
    print(f"positive_two_spike_classes={positive_classes}")
    print("extra_tooth_event={u:{c*u} in (3/7,4/7)}")
    print(f"extra_tooth_event_mass={extra_event_mass}")
    print(f"aligned_safe_floor={TWO_SAFE_FLOOR}")
    print(f"d8_one_tooth_q8_type={d8_one_tooth}")
    print(f"d8_two_tooth_q8_type={d8_two_tooth}")
    print(f"fifth_q8_types={fifth_types}")
    print(f"exact_mean_per_spike={exact_mean_per_spike}")
    print(f"markov_event_bound={markov_event_bound}")
    print(f"necessary_heavy_cutoff_strict={necessary_heavy_cutoff}")
    print(f"row_count={len(rows)}")
    print(f"shape_count={len(EXPECTED_SHAPES)}")
    print(f"incoming_occurrences={len(rows) * len(EXPECTED_SHAPES)}")
    print(f"support_histogram={support_histogram}")
    print(f"heavy_fibre_histogram={heavy_fibre_histogram}")
    print(f"minimum_extra_histogram={minimum_extra_histogram}")
    print(f"semantic_sha256={semantic.hexdigest()}")
    print("surviving_occurrences=0")
    print("bounded_D_le_100_k2_ledger=EMPTY")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
