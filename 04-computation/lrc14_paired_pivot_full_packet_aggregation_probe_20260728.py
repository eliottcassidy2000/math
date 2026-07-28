#!/usr/bin/env python3
"""Exact pivot-aggregation scout for THM-2701 Section 10.

This is an untracked research scout.  It compares the two THM-2698 event
cylinders with the complete THM-2707 fixed-skeleton packet orbit.  All wall
tests use Fraction arithmetic.  Fourier nonvanishing is certified through
the prime-cyclotomic criterion: a rational length-13 profile has a vanishing
primitive 13-character iff it is constant.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import permutations
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if not (ROOT / "04-computation").is_dir():
    ROOT = Path.cwd()
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_mixed_slope_long_word_probe as old


P = 13
R = P**6
A = P**3
C = 2 * P**5
PIVOTS = (27, 40, 53, 66)
EVENTS = (
    (Fraction(55_232_507, 24 * R), 7),
    (Fraction(58_313_459, 24 * R), 3),
)
PACKET_RADIUS = Fraction(1, 301_082_946_198_216)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def floor_fraction(value):
    return value.numerator // value.denominator


def frac(value):
    return value - floor_fraction(value)


def centered(value):
    value = frac(value)
    return value if value < Fraction(1, 2) else value - 1


def danger(value):
    return abs(centered(value)) < Fraction(1, 14)


def wall_distance(value):
    base = floor_fraction(value)
    return min(
        abs(value - (integer + sign * Fraction(1, 14)))
        for integer in range(base - 2, base + 3)
        for sign in (-1, 1)
    )


def support(center, root, pivot):
    require(danger(C * center - Fraction(root, P)), "deep root not dangerous")
    return tuple(
        s for s in range(P)
        if not danger(A * center - Fraction(s, P))
        and not danger(pivot * center + Fraction(s, P))
    )


def profile(support_set):
    support_set = set(support_set)
    return tuple(int(s in support_set) for s in range(P))


def add_profiles(rows):
    return tuple(sum(row[s] for row in rows) for s in range(P))


def quotient_rank(rows):
    """Rank after quotienting profiles by the constant vector."""
    # Fourier on nonzero characters is injective modulo constants.  Subtract
    # coordinate 12 from coordinates 0..11 to represent that quotient.
    qrows = [tuple(row[s] - row[-1] for s in range(P - 1)) for row in rows]
    return rational_rank(qrows)


def rational_rank(rows):
    """Exact row rank over Q for a narrow integer/rational matrix."""
    basis = {}
    for raw in rows:
        row = [Fraction(value) for value in raw]
        for pivot in sorted(basis):
            if row[pivot]:
                coefficient = row[pivot]
                old_row = basis[pivot]
                row = [left - coefficient * right
                       for left, right in zip(row, old_row)]
        pivot = next((i for i, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        scale = row[pivot]
        basis[pivot] = [value / scale for value in row]
    return len(basis)


def candidate_constant_kernel(rows_by_pivot):
    """Kernel of pivot weights whose support aggregate is constant."""
    columns = []
    for pivot in PIVOTS:
        row = rows_by_pivot[pivot]
        columns.append(tuple(row[s] - row[-1] for s in range(P - 1)))
    matrix = tuple(tuple(columns[j][i] for j in range(len(PIVOTS)))
                   for i in range(P - 1))
    rank = rational_rank(matrix)
    return rank, len(PIVOTS) - rank


def affine_image(values, a, b):
    return frozenset((a * x + b) % P for x in values)


def paired_event_symmetries(event_supports):
    """AGL(1,13) maps from event 0 to 1, allowing a pivot permutation."""
    left = event_supports[0]
    right = event_supports[1]
    answers = []
    for a in range(1, P):
        for b in range(P):
            images = {pivot: affine_image(left[pivot], a, b) for pivot in PIVOTS}
            for perm in permutations(PIVOTS):
                mapping = dict(zip(PIVOTS, perm))
                if all(images[pivot] == frozenset(right[mapping[pivot]])
                       for pivot in PIVOTS):
                    answers.append((a, b, tuple(mapping[pivot] for pivot in PIVOTS)))
    return tuple(answers)


def paired_object_automorphisms(event_supports, swap):
    answers = []
    for a in range(1, P):
        for b in range(P):
            for perm in permutations(PIVOTS):
                mapping = dict(zip(PIVOTS, perm))
                if all(
                    affine_image(event_supports[i][pivot], a, b)
                    == frozenset(event_supports[1 - i if swap else i][mapping[pivot]])
                    for i in range(2) for pivot in PIVOTS
                ):
                    answers.append((a, b, perm))
    return tuple(answers)


def dihedral_canonical(values):
    values = tuple(values)
    images = []
    for sign in (1, -1):
        for b in range(P):
            images.append(tuple(sorted((sign * x + b) % P for x in values)))
    return min(images)


def agl_canonical(values):
    values = tuple(values)
    return min(
        tuple(sorted((a * x + b) % P for x in values))
        for a in range(1, P)
        for b in range(P)
    )


def transform_profile(row, a, b):
    out = [0] * P
    for x, value in enumerate(row):
        out[(a * x + b) % P] = value
    return tuple(out)


def joint_orbit_key(rows, multipliers):
    return min(
        tuple(transform_profile(row, a, b) for row in rows)
        for a in multipliers for b in range(P)
    )


def reconstruct_thm2707_addresses():
    """Independent exact replay of the 4.8m-point midpoint scan."""
    m = old.m
    x = Fraction(649_039_434_905_733, 1_304_692_766_858_936)
    z = Fraction(46_873_542_509_301, 100_360_982_066_072)
    radius = Fraction(1, 1_304_692_766_858_936)
    module, prefixes, _, _, rails, present, starts = m.core.build_carrier_data()
    rows = m.shard((0, 1))[6][0]
    rail_support = old.merge_support(rails[0][3])
    present_support = tuple(present[1, (-6) % P])

    good_residues = []
    for residue in range(P):
        carry = (2 + 7 * residue) % P
        root = (6 + residue) % P
        if root and m.is_unit(rows[0][0][carry][1][6], root, 26):
            good_residues.append(residue)
    require(tuple(good_residues) == (0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12),
            "THM2707 residue bank changed")

    denominator = (z * m.T).denominator
    point = (z * m.T).numerator
    modulus = m.T * denominator
    step = (7 * m.T // R) * denominator
    scaled_rail = tuple((left * denominator, right * denominator)
                        for left, right in rail_support)
    scaled_present = tuple((left * denominator, right * denominator)
                           for left, right in present_support)
    rail_starts = tuple(left for left, _ in scaled_rail)
    present_starts = tuple(left for left, _ in scaled_present)

    def strict_interval_index(value, starts_, intervals):
        # Local binary search avoided to keep this scout dependency-free apart
        # from the already used SymPy rank calculation.
        from bisect import bisect_right
        index = bisect_right(starts_, value) - 1
        return index >= 0 and intervals[index][0] < value < intervals[index][1]

    good = []
    for n in range(R):
        if (n % P in good_residues
                and strict_interval_index(point, rail_starts, scaled_rail)
                and strict_interval_index(point, present_starts, scaled_present)):
            good.append(n)
        point = (point + step) % modulus
    require(len(good) == 3346 and 110 in good, "THM2707 packet census changed")
    return tuple(good), x, z, radius


def main():
    print("LRC14 PAIRED-PIVOT / FULL-PACKET AGGREGATION EXACT AUDIT")
    print("status=FINITE-EXACT moving-character and pivot-covariance result; "
          "endpoint current remains open")
    print("half_cycle_universe=2 THM2698 event cylinders x 4 pivots x 13 shifts")
    print("packet_universe=all n in Z/(13^6), filtered by THM2707 private-unit "
          "residues plus literal rail and present support; paired factors "
          "then inserted on the whole inherited open I")
    event_supports = []
    event_profiles = []
    for x, root in EVENTS:
        supports = {pivot: support(x, root, pivot) for pivot in PIVOTS}
        profiles = {pivot: profile(supports[pivot]) for pivot in PIVOTS}
        event_supports.append(supports)
        event_profiles.append(profiles)

    print("HALF_CYCLE")
    for i in range(2):
        print(f"event={i} supports={event_supports[i]}")
        rows = [event_profiles[i][pivot] for pivot in PIVOTS]
        print(f"event={i} support_rank={rational_rank(rows)} "
              f"moving_fourier_rank={quotient_rank(rows)} "
              f"pivot_constant_kernel={candidate_constant_kernel(event_profiles[i])}")
        aggregate = add_profiles(rows)
        print(f"event={i} equal_pivot_profile={aggregate} "
              f"constant={len(set(aggregate)) == 1}")

    paired_rows = [event_profiles[i][pivot]
                   for i in range(2) for pivot in PIVOTS]
    paired_concat = [event_profiles[0][pivot] + event_profiles[1][pivot]
                     for pivot in PIVOTS]
    print(f"eight_row_support_rank={rational_rank(paired_rows)} "
          f"eight_row_moving_rank={quotient_rank(paired_rows)}")
    print(f"fixed_pivot_two_event_concat_rank={rational_rank(paired_concat)}")
    joint_kernel_rows = {
        pivot: event_profiles[0][pivot] + event_profiles[1][pivot]
        for pivot in PIVOTS
    }
    # Quotient by one independent constant on each event.
    block_columns = []
    for pivot in PIVOTS:
        row = joint_kernel_rows[pivot]
        block_columns.append(tuple(row[s] - row[12] for s in range(12))
                             + tuple(row[13 + s] - row[25] for s in range(12)))
    block_matrix = tuple(tuple(block_columns[j][i] for j in range(4))
                         for i in range(24))
    block_rank = rational_rank(block_matrix)
    print(f"joint_pivot_constant_rank={block_rank} "
          f"kernel_dimension={4 - block_rank}")
    capacity = {
        pivot: (min(len(event_supports[i][pivot]) for i in range(2)),
                sum(len(event_supports[i][pivot]) for i in range(2)),
                len(event_supports[0][pivot]) * len(event_supports[1][pivot]))
        for pivot in PIVOTS
    }
    print(f"capacity=(min,sum,product)={capacity}")
    best_lex = max(PIVOTS, key=lambda pivot: capacity[pivot])
    print(f"unique_lexicographic_maximin_pivot={best_lex}")
    print(f"event_swap_AGL_symmetries={paired_event_symmetries(event_supports)}")
    print(f"paired_object_fixed_AGL_automorphisms="
          f"{paired_object_automorphisms(event_supports, False)}")
    print(f"paired_object_swap_AGL_automorphisms="
          f"{paired_object_automorphisms(event_supports, True)}")
    for kind, canonical in (("D13", dihedral_canonical), ("AGL1", agl_canonical)):
        classes = defaultdict(list)
        for i in range(2):
            for pivot in PIVOTS:
                classes[canonical(event_supports[i][pivot])].append((i, pivot))
        print(f"{kind}_support_orbit_count={len(classes)} "
              f"orbit_size_multiset={tuple(sorted(map(len, classes.values())))}")

    good, x, z, radius = reconstruct_thm2707_addresses()
    residue_counts = Counter(n % P for n in good)
    expected_residue_counts = {
        0: 304, 1: 305, 2: 304, 3: 305, 4: 304, 5: 305,
        6: 304, 9: 301, 10: 304, 11: 305, 12: 305,
    }
    require(dict(sorted(residue_counts.items())) == expected_residue_counts,
            "THM2707 packet residue census changed")
    atom_count = residue_counts[0]
    transit_count = len(good) - atom_count
    require((atom_count, transit_count, atom_count * transit_count)
            == (304, 3042, 924768),
            "following-atom macro census changed")
    module = old.m
    q_radius = 13 * radius
    address_profiles = []
    support_size_census = {pivot: Counter() for pivot in PIVOTS}
    aggregate_constant = 0
    all_rows_nonempty_proper = True
    wall_clear = True
    minimum_wall_ratio = None
    pivot_wins = Counter()
    total_by_pivot = {pivot: [0] * P for pivot in PIVOTS}
    all_total = [0] * P
    signature_census = Counter()
    group_totals = {
        "atom": {pivot: [0] * P for pivot in PIVOTS},
        "transit": {pivot: [0] * P for pivot in PIVOTS},
    }
    group_size_census = {
        "atom": {pivot: Counter() for pivot in PIVOTS},
        "transit": {pivot: Counter() for pivot in PIVOTS},
    }
    group_pivot_wins = {"atom": Counter(), "transit": Counter()}
    group_profile_census = {
        "atom": {pivot: Counter() for pivot in PIVOTS},
        "transit": {pivot: Counter() for pivot in PIVOTS},
    }
    group_joint_profile_census = {"atom": Counter(), "transit": Counter()}

    for n in good:
        center = frac(z + Fraction(7 * n, R))
        root = (6 + n) % P
        group = "atom" if n % P == 0 else "transit"
        profiles = {}
        sizes = {}
        for pivot in PIVOTS:
            values = support(center, root, pivot)
            row = profile(values)
            profiles[pivot] = row
            group_profile_census[group][pivot][row] += 1
            sizes[pivot] = len(values)
            support_size_census[pivot][len(values)] += 1
            group_size_census[group][pivot][len(values)] += 1
            all_rows_nonempty_proper &= 0 < len(values) < P
            for s in range(P):
                total_by_pivot[pivot][s] += row[s]
                group_totals[group][pivot][s] += row[s]
                all_total[s] += row[s]

            # Verify every inserted factor is constant on the inherited open
            # x-cylinder.  q_n has derivative 13 with respect to x.
            phases = [(C, C * center - Fraction(root, P))]
            phases += [(A, A * center - Fraction(s, P)) for s in range(P)]
            phases += [(pivot, pivot * center + Fraction(s, P)) for s in range(P)]
            for coefficient, phase in phases:
                ratio = wall_distance(phase) / (abs(coefficient) * q_radius)
                minimum_wall_ratio = ratio if minimum_wall_ratio is None else min(
                    minimum_wall_ratio, ratio)
                wall_clear &= ratio >= 1

        aggregate = add_profiles([profiles[pivot] for pivot in PIVOTS])
        aggregate_constant += int(len(set(aggregate)) == 1)
        max_size = max(sizes.values())
        winners = tuple(pivot for pivot in PIVOTS if sizes[pivot] == max_size)
        pivot_wins[winners] += 1
        group_pivot_wins[group][winners] += 1
        signature_census[tuple(sizes[pivot] for pivot in PIVOTS)] += 1
        group_joint_profile_census[group][tuple(profiles[pivot] for pivot in PIVOTS)] += 1
        address_profiles.extend(profiles[pivot] for pivot in PIVOTS)

    print("THM2707_FULL_PACKET_ORBIT")
    print(f"residue_counts={tuple(sorted(residue_counts.items()))} "
          f"atom_endpoints={atom_count} transit_packets={transit_count} "
          f"atom_to_transit_edges={atom_count * transit_count}")
    print(f"addresses={len(good)} inserted_factors_constant_on_I={wall_clear} "
          f"minimum_wall_radius_ratio={minimum_wall_ratio}")
    print(f"all_pivot_supports_nonempty_proper={all_rows_nonempty_proper}")
    print(f"support_size_census={support_size_census}")
    print(f"per_address_equal_pivot_constant_profiles={aggregate_constant}")
    print(f"argmax_pivot_census={pivot_wins}")
    print(f"size_signature_census={signature_census}")
    for pivot in PIVOTS:
        row = tuple(total_by_pivot[pivot])
        print(f"pivot={pivot} address_aggregate_profile={row} "
              f"constant={len(set(row)) == 1}")
    print(f"all_address_all_pivot_profile={tuple(all_total)} "
          f"constant={len(set(all_total)) == 1}")
    full_masses = {pivot: sum(total_by_pivot[pivot]) for pivot in PIVOTS}
    print(f"full_packet_support_masses={full_masses} "
          f"unique_total_mass_selector="
          f"{max(PIVOTS, key=lambda pivot: full_masses[pivot])}")
    for group in ("atom", "transit"):
        masses = {}
        for pivot in PIVOTS:
            row = tuple(group_totals[group][pivot])
            masses[pivot] = sum(row)
            print(f"{group}_pivot={pivot}_aggregate_profile={row} "
                  f"constant={len(set(row)) == 1}")
        print(f"{group}_support_masses={masses} size_census={group_size_census[group]}")
        print(f"{group}_argmax_pivot_census={group_pivot_wins[group]}")
    for pivot in PIVOTS:
        matching_edges = sum(
            count * group_profile_census["transit"][pivot][row]
            for row, count in group_profile_census["atom"][pivot].items()
        )
        print(f"atom_to_transit_exact_profile_matching_edges_pivot_{pivot}="
              f"{matching_edges}/924768")
    joint_matching_edges = sum(
        count * group_joint_profile_census["transit"][row]
        for row, count in group_joint_profile_census["atom"].items()
    )
    print(f"atom_to_transit_exact_joint_profile_matching_edges="
          f"{joint_matching_edges}/924768")
    for orbit_name, canonical in (
            ("translation", lambda row: joint_orbit_key((row,), (1,))[0]),
            ("dihedral", lambda row: joint_orbit_key((row,), (1, -1))[0]),
            ("AGL", lambda row: joint_orbit_key((row,), tuple(range(1, P)))[0])):
        for pivot in PIVOTS:
            left = Counter()
            right = Counter()
            for row, count in group_profile_census["atom"][pivot].items():
                left[canonical(row)] += count
            for row, count in group_profile_census["transit"][pivot].items():
                right[canonical(row)] += count
            matches = sum(count * right[key] for key, count in left.items())
            print(f"atom_to_transit_{orbit_name}_profile_matching_edges_"
                  f"pivot_{pivot}={matches}/924768")
        left = Counter()
        right = Counter()
        multipliers = ((1,) if orbit_name == "translation" else
                       (1, -1) if orbit_name == "dihedral" else
                       tuple(range(1, P)))
        for rows, count in group_joint_profile_census["atom"].items():
            left[joint_orbit_key(rows, multipliers)] += count
        for rows, count in group_joint_profile_census["transit"].items():
            right[joint_orbit_key(rows, multipliers)] += count
        matches = sum(count * right[key] for key, count in left.items())
        print(f"atom_to_transit_{orbit_name}_joint_profile_matching_edges="
              f"{matches}/924768")

    def joint_profiles_at(node):
        center = frac(z + Fraction(7 * node, R))
        root = (6 + node) % P
        return tuple(profile(support(center, root, pivot)) for pivot in PIVOTS)

    loop_base = joint_profiles_at(0)
    loop_positive = joint_profiles_at(1)
    loop_hostile = joint_profiles_at(106)
    require(loop_base == loop_positive,
            "canonical 0->1 atom macro lost its exact profile control")
    require(joint_orbit_key(loop_base, tuple(range(1, P)))
            != joint_orbit_key(loop_hostile, tuple(range(1, P))),
            "hostile atom macro unexpectedly gained an AGL profile map")
    print("macro_profile_positive_control=0->1->0 numerator=(7,-7) "
          "preserves the exact four-pivot profile")
    print("macro_profile_hostile_control=0->106->0 numerator=(742,-742) "
          "does not preserve the joint profile even up to one common "
          "AGL(1,13) shift reindexing")
    print(f"all_support_row_rank={rational_rank(address_profiles)} "
          f"all_moving_profile_rank={quotient_rank(address_profiles)}")

    # Simultaneous pivot weights that make every address profile constant.
    constraint_rows = []
    for index in range(0, len(address_profiles), 4):
        rows = address_profiles[index:index + 4]
        for s in range(P - 1):
            constraint_rows.append([row[s] - row[-1] for row in rows])
    constraint_rank = rational_rank(constraint_rows)
    print(f"global_pivot_constant_constraint_rank={constraint_rank} "
          f"kernel_dimension={4 - constraint_rank}")

    # Exact displayed THM2707 cycles.
    cycles = {
        "triangle": (110, 689654, 1379198),
        "cycle11": (110, 3447831, 2068743, 1379199, 4137376, 2758288,
                    2068744, 1379200, 3447832, 689655, 2758287),
    }
    good_set = set(good)
    for name, nodes in cycles.items():
        require(all(node in good_set for node in nodes), f"{name} left packet orbit")
        cycle_supports = {pivot: [] for pivot in PIVOTS}
        min_sizes = {}
        for node in nodes:
            center = frac(z + Fraction(7 * node, R))
            root = (6 + node) % P
            for pivot in PIVOTS:
                cycle_supports[pivot].append(profile(support(center, root, pivot)))
        aggregate_rows = {
            pivot: add_profiles(cycle_supports[pivot]) for pivot in PIVOTS
        }
        for pivot in PIVOTS:
            min_sizes[pivot] = min(sum(row) for row in cycle_supports[pivot])
        print(f"{name}_min_support_by_pivot={min_sizes}")
        for pivot in PIVOTS:
            row = aggregate_rows[pivot]
            print(f"{name}_pivot={pivot}_aggregate={row} "
                  f"constant={len(set(row)) == 1}")
        total = add_profiles(list(aggregate_rows.values()))
        print(f"{name}_all_pivot_aggregate={total} constant={len(set(total)) == 1}")

    # Missing-left-index hostile: the same observed fixed-left row K(t) has
    # both a frozen-left lift A(s,t)=K(t) and a simultaneously covariant lift
    # A(s,t)=K(t-s).  Their common-twist diagonals are, respectively, K(s)
    # and the constant K(0), so one has all moving colours and the other has
    # none.  This verifies information-theoretic non-identifiability and
    # pinpoints covariance of the absent left variable as the deciding datum.
    for i in range(2):
        for pivot in PIVOTS:
            K = event_profiles[i][pivot]
            frozen_left_diagonal = tuple(K[s] for s in range(P))
            covariant_diagonal = tuple(K[0] for _ in range(P))
            require(len(set(frozen_left_diagonal)) > 1
                    and len(set(covariant_diagonal)) == 1,
                    "hostile endpoint lifts lost their contrast")
    print("missing_left_hostile=same fixed-left row A(0,t)=K(t) admits "
          "A_frozen(s,t)=K(t), with nonconstant common diagonal K(s), and "
          "A_covariant(s,t)=K(t-s), with constant common diagonal K(0)")
    print("SCOPE=all Fourier claims above concern moving-shift profiles; by "
          "THM2568 a lawful danger-to-safe full-X endpoint completion is "
          "zero in every coarse target character")


if __name__ == "__main__":
    main()
