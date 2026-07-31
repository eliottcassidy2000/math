#!/usr/bin/env python3
"""Exact target-deck J probe on the minimal c2-safe ghost-transit bank.

The lawful coefficient action fixes the semantic rail and edge, sends the
delayed digit (h,kappa)=(0,1) to (12,0), shifts carry by +6 so the private
root is unchanged, and applies Frobenius z->z^13=z^(-1) on F13[C7].
No action on rails, owners, or present-factor supports is assumed.
"""

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
import hashlib
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = (ROOT / "04-computation" /
        "lrc14_central_half_odometer_full_local_cycle_thm2698.py")
BASE_SHA256 = "45cc393a856c00342fdf84875a0bc5a6d4c3df196ab35bb9ac2aad3cfc966c25"


def lf_sha256(path):
    data = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


if lf_sha256(BASE) != BASE_SHA256:
    raise RuntimeError("tracked THM-2698 dependency hash changed")

SPEC = importlib.util.spec_from_file_location("half", BASE)
half = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(half)

private = half.private
core = half.core
old = private.old
P, Q7, T = private.P, private.Q7, private.T
R, S = private.R, private.S
CONTENT = 26
FORWARD, REVERSE = (0, 1), (12, 0)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def freeze(value):
    if isinstance(value, list):
        return tuple(freeze(x) for x in value)
    return value


def transit_word(module):
    """The c2-safe changed grammar, rebuilt from tracked THM-2698 data."""
    deletion = module.make_comb(module.W[module.GUARD], 91, -13, 13)
    for index in module.UNIT_IDX:
        deletion = module.subtract_comb(
            deletion, module.W[index], 182, -13, 13
        )
    deletion = module.subtract_comb(deletion, module.C3, 182, -13, 13)
    return module.subtract_comb(deletion, module.C2, 182, -13, 13)


def transit_prefixes(module):
    """Return the complete [clock][half-digit][bit] safe prefix bank."""
    word = transit_word(module)
    result = []
    masses = []
    for ell in range(Q7):
        qell = module.subtract_comb(
            word, module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        by_h = []
        by_h_mass = []
        for h in range(P):
            pair = []
            pair_mass = []
            for kappa in range(2):
                digit = 2 * h + kappa
                interval = old.sat.intersect_interval(
                    qell, digit * T // (2 * P),
                    (digit + 1) * T // (2 * P)
                )
                prefix = module.make_prefix(interval)
                pair.append(prefix)
                pair_mass.append(prefix[2][-1])
            by_h.append(tuple(pair))
            by_h_mass.append(tuple(pair_mass))
        result.append(tuple(by_h))
        masses.append(tuple(by_h_mass))
    require(all(mass % P == 0
                for by_h in masses for pair in by_h for mass in pair),
            "changed delayed half-digit mass does not descend")
    return tuple(result), tuple(masses)


def primitive_data(values, root, content):
    require(root != 0 and content > 0,
            "primitive test needs nonzero root/content")
    require(all(value % content == 0 for value in values),
            "profile escaped global content")
    scale = pow(root, -1, P)
    normalized = tuple((value // content) * scale % P for value in values)
    reduced = tuple((normalized[i] - normalized[-1]) % P
                    for i in range(Q7 - 1))
    determinant = old.sat.multiplication_determinant_7(reduced)
    return normalized, reduced, determinant


def targeted_shard(bounds):
    """Build only the two digit slices needed by the J test."""
    module, _, _, _, rails, present, starts = core.build_carrier_data()
    prefixes, _ = transit_prefixes(module)
    caches = {(ell, h): {} for ell in range(Q7) for h, _ in (FORWARD, REVERSE)}
    metadata, rows = [], []
    for j in range(*bounds):
        source, owner, theta, pieces = rails[j]
        metadata.append((j, source, owner, theta))
        row = {}
        for h, kappa in (FORWARD, REVERSE):
            for edge in range(2):
                for carry in range(P):
                    row[h, kappa, edge, carry] = [0] * Q7
            for ell in range(Q7):
                overlap = old.intersect_weighted_union(
                    pieces, present[ell, (-h) % P], starts[ell, (-h) % P]
                )
                for root in range(1, P):
                    halves = (
                        old.intersect_weighted_comb(
                            overlap, module.C3, 182, 14 * root - 13, 14 * root
                        ),
                        old.intersect_weighted_comb(
                            overlap, module.C3, 182, 14 * root, 14 * root + 13
                        ),
                    )
                    for edge, half_tooth in enumerate(halves):
                        values = private.delayed_carry_pair(
                            half_tooth, prefixes[ell][h], caches[ell, h]
                        )
                        for carry in range(P):
                            predicted = (
                                2 * carry + (2 * h + kappa) // P
                                + (edge == 0)
                            ) % P
                            value = values[carry][kappa]
                            if root == predicted:
                                row[h, kappa, edge, carry][ell] = value
                            else:
                                require(value == 0, "targeted row lost root privacy")
        rows.append({key: tuple(value) for key, value in row.items()})
    return tuple(metadata), tuple(rows)


def canonical(values, edge, carry, h, kappa):
    if not any(values):
        return None
    root = (2 * carry + (2 * h + kappa) // P + (edge == 0)) % P
    require(root != 0, "positive profile has zero private root")
    normalized, reduced, determinant = primitive_data(values, root, CONTENT)
    return root, normalized, reduced, determinant


def recanonicalize(rep):
    return tuple((rep[i] - rep[6]) % P for i in range(6))


def frobenius_inverse(reduced):
    rep = tuple(reduced) + (0,)
    inverted = tuple(rep[(-i) % Q7] for i in range(Q7))
    return recanonicalize(inverted)


def rotate(reduced, shift):
    rep = tuple(reduced) + (0,)
    shifted = tuple(rep[(i - shift) % Q7] for i in range(Q7))
    return recanonicalize(shifted)


def scalar_match(left, right):
    if left == right:
        return 1
    pivot = next((i for i, x in enumerate(left) if x), None)
    if pivot is None or right[pivot] == 0:
        return None
    scale = right[pivot] * pow(left[pivot], -1, P) % P
    return scale if all(scale * x % P == y for x, y in zip(left, right)) else None


def cell_counter(y, clock, rail_js, banks, module, rails, present, starts,
                 safe_prefixes):
    """Count unit cells by semantic signature; no rail transformation."""
    h = half.floor_fraction(P * y)
    kappa = half.floor_fraction(2 * P * y) - 2 * h
    require(half.strict_interval_member(
                y, half.prefix_intervals(safe_prefixes[clock][h][kappa])),
            "requested ghost phase is outside the c2-safe semantic prefix")
    first_n = half.floor_fraction(F(6 * R, P) - y) + 1
    last_n = half.floor_fraction(F(7 * R, P) - y)
    rail_starts = [[left for left, _, _ in rail[3]] for rail in rails[:14]]
    result = Counter()
    for n in range(first_n, last_n + 1):
        x = F(n, R) + F(y, R)
        if half.shallow(x) != clock or half.owner(x) != clock:
            continue
        coordinate = x * half.T
        if not half.strict_interval_member(
                coordinate, present[clock, (-h) % P], starts[clock, (-h) % P]):
            continue
        carry = n % P
        for j in rail_js:
            if not half.strict_weighted_member(
                    coordinate, rails[j][3], rail_starts[j]):
                continue
            for edge in range(2):
                root = (2 * carry + (2 * h + kappa) // P
                        + (edge == 0)) % P
                if root == 0 or not half.half_membership(module, x, edge, root):
                    continue
                profile = canonical(banks[j][h, kappa, edge, carry],
                                    edge, carry, h, kappa)
                if profile is not None and profile[-1]:
                    result[j, edge, carry, root, profile[2]] += 1
    return result


def main():
    with ProcessPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(targeted_shard, core.SHARDS))
    metadata = {}
    banks = {}
    for metas, rows in results:
        for meta, row in zip(metas, rows):
            metadata[meta[0]] = meta[1:]
            banks[meta[0]] = row
    require(len(banks) == 162, "targeted bank is incomplete")

    stats = Counter()
    swap_stats = Counter()
    mismatch_weights = Counter()
    first_mismatch = None
    first_nonzero_exact = None
    nonzero_exact_matches = []
    for j, row in banks.items():
        for edge in range(2):
            for carry in range(P):
                forward = canonical(row[FORWARD + (edge, carry)],
                                    edge, carry, *FORWARD)
                reverse_carry = (carry + 6) % P
                reverse = canonical(row[REVERSE + (edge, reverse_carry)],
                                    edge, reverse_carry, *REVERSE)
                if forward is not None:
                    stats["forward_nonzero"] += 1
                    stats["forward_unit"] += bool(forward[-1])
                if reverse is not None:
                    stats["reverse_nonzero"] += 1
                    stats["reverse_unit"] += bool(reverse[-1])
                if forward is None or reverse is None:
                    stats["one_sided"] += (forward is None) != (reverse is None)
                    continue
                stats["paired_nonzero"] += 1
                require(forward[0] == reverse[0], "carry relabel did not preserve root")
                image = frobenius_inverse(forward[2])
                image_det = old.sat.multiplication_determinant_7(image)
                require(image_det == forward[-1],
                        "Frobenius changed the multiplication determinant")
                if image == reverse[2]:
                    stats["exact_frobenius"] += 1
                    stats["exact_nonzero_class"] += bool(any(image))
                    stats["exact_unit"] += bool(forward[-1] and reverse[-1])
                    if any(image):
                        stats[f"exact_owner_{metadata[j][1]}"] += 1
                    if any(image) and first_nonzero_exact is None:
                        first_nonzero_exact = (
                            j, metadata[j], edge, carry, reverse_carry,
                            forward[0], image, forward[-1], reverse[-1],
                        )
                    if any(image):
                        nonzero_exact_matches.append((
                            j, metadata[j], edge, carry, reverse_carry,
                            forward[0], image, forward[-1], reverse[-1],
                        ))
                if scalar_match(image, reverse[2]) is not None:
                    stats["scalar_frobenius"] += 1
                    stats["scalar_nonzero_class"] += bool(any(image))
                if any(scalar_match(rotate(image, shift), reverse[2]) is not None
                       for shift in range(Q7)):
                    stats["affine_clock_match"] += 1
                    stats["affine_nonzero_class"] += bool(any(image))
                weights = (sum(bool(x) for x in image),
                           sum(bool(x) for x in reverse[2]))
                mismatch_weights[weights] += 1
                if image != reverse[2] and first_mismatch is None:
                    first_mismatch = (
                        j, metadata[j], edge, carry, reverse_carry, forward[0],
                        forward[2], image, reverse[2], weights,
                    )

                # Larger non-inherited candidate: swap private half-tooth
                # edges while preserving the same root.  Left->right keeps
                # carry; right->left sends carry to carry-1.
                swapped_edge = 1 - edge
                swapped_carry = (carry + (0 if edge == 0 else -1)) % P
                swapped = canonical(
                    row[REVERSE + (swapped_edge, swapped_carry)],
                    swapped_edge, swapped_carry, *REVERSE
                )
                if swapped is not None:
                    swap_stats["paired_nonzero"] += 1
                    require(swapped[0] == forward[0],
                            "edge-swap relabel did not preserve root")
                    if image == swapped[2]:
                        swap_stats["exact_frobenius"] += 1
                        swap_stats["exact_nonzero_class"] += bool(any(image))
                        swap_stats["exact_unit"] += bool(forward[-1] and swapped[-1])
                    if scalar_match(image, swapped[2]) is not None:
                        swap_stats["scalar_frobenius"] += 1
                        swap_stats["scalar_nonzero_class"] += bool(any(image))
                    if any(scalar_match(rotate(image, shift), swapped[2]) is not None
                           for shift in range(Q7)):
                        swap_stats["affine_clock_match"] += 1
                        swap_stats["affine_nonzero_class"] += bool(any(image))

    # Exact forced witnesses and their strongest coefficient hostile.
    forced_forward = canonical(banks[9][FORWARD + (0, 7)], 0, 7, *FORWARD)
    forced_reverse = canonical(banks[2][REVERSE + (0, 0)], 0, 0, *REVERSE)
    same_rail_reverse = canonical(banks[9][REVERSE + (0, 0)], 0, 0, *REVERSE)
    swapped_same_rail = canonical(banks[9][REVERSE + (1, 7)], 1, 7, *REVERSE)
    forced_image = frobenius_inverse(forced_forward[2])
    require(forced_forward[0] == forced_reverse[0] == same_rail_reverse[0] == 2,
            "forced root labels changed")
    require(forced_image != forced_reverse[2]
            and forced_image != same_rail_reverse[2],
            "forced coefficient hostile unexpectedly vanished")
    require(not any(scalar_match(rotate(forced_image, shift), target) is not None
                    for shift in range(Q7)
                    for target in (same_rail_reverse[2], forced_reverse[2])),
            "forced affine-clock hostile unexpectedly vanished")
    require(swapped_same_rail is not None
            and not any(scalar_match(rotate(forced_image, shift),
                                     swapped_same_rail[2]) is not None
                        for shift in range(Q7)),
            "forced edge-swap hostile unexpectedly vanished")
    require((sum(bool(x) for x in forced_image),
             sum(bool(x) for x in same_rail_reverse[2]),
             sum(bool(x) for x in forced_reverse[2])) == (3, 2, 1),
            "forced clock-support mismatch changed")
    require(len(nonzero_exact_matches) == 14
            and all(item[0] == 104 and item[1] == (8, 6, 12)
                    and all(x == 0 for x in item[6][1:])
                    for item in nonzero_exact_matches),
            "exact Frobenius survivor is no longer the source-eight constant bank")

    # Compare actual safe cell supports without inventing a rail action.
    module, _, _, _, rails, present, starts = core.build_carrier_data()
    safe_prefixes = transit_prefixes(module)[0]
    f4 = cell_counter(F(1, 17), 4, (8, 9), banks,
                      module, rails, present, starts, safe_prefixes)
    j4 = cell_counter(F(16, 17), 4, (8, 9), banks,
                      module, rails, present, starts, safe_prefixes)
    r1 = cell_counter(F(16, 17), 1, (2, 3), banks,
                      module, rails, present, starts, safe_prefixes)
    j1 = cell_counter(F(1, 17), 1, (2, 3), banks,
                      module, rails, present, starts, safe_prefixes)

    def transformed(counter, forward_to_reverse=True, swap_edge=False):
        out = Counter()
        for (j, edge, carry, root, reduced), multiplicity in counter.items():
            if swap_edge:
                next_edge = 1 - edge
                if forward_to_reverse:
                    delta = 0 if edge == 0 else -1
                else:
                    delta = 0 if edge == 1 else 1
            else:
                next_edge = edge
                delta = 6 if forward_to_reverse else -6
            out[j, next_edge, (carry + delta) % P, root,
                frobenius_inverse(reduced)] += multiplicity
        return out

    cell_intersections = {
        "clock4": sum((transformed(f4) & j4).values()),
        "clock1": sum((transformed(j1) & r1).values()),
        "clock4_edge_swap": sum((transformed(f4, swap_edge=True) & j4).values()),
        "clock1_edge_swap": sum((transformed(j1, swap_edge=True) & r1).values()),
    }

    print("LRC14 C2-SAFE TARGET-DECK J EXACT PROBE")
    print(f"lawful_action=(rail_fixed,edge_fixed,digit_1_to_24,carry_plus_6,root_fixed,clock_Frobenius_inverse)")
    print(f"profile_stats={tuple(sorted(stats.items()))}")
    print(f"edge_swap_profile_stats={tuple(sorted(swap_stats.items()))}")
    print(f"weight_pair_hist={tuple(sorted(mismatch_weights.items()))}")
    print(f"first_same_support_mismatch={first_mismatch}")
    print(f"first_nonzero_exact_match={first_nonzero_exact}")
    print(f"nonzero_exact_matches={tuple(nonzero_exact_matches)}")
    print("forced_profiles="
          f"(forward={forced_forward[2]},J_image={forced_image},"
          f"same_rail_reverse={same_rail_reverse[2]},"
          f"edge_swap_reverse={swapped_same_rail[2]},"
          f"physical_reverse={forced_reverse[2]})")
    print("cell_bank_sizes="
          f"(forward_clock4={sum(f4.values())},Jghost_clock4={sum(j4.values())},"
          f"Jghost_clock1={sum(r1.values())},forwardghost_clock1={sum(j1.values())})")
    print("cell_signature_data="
          f"(signature_counts={(len(f4),len(j4),len(r1),len(j1))},"
          f"transformed_intersections={cell_intersections},"
          f"clock4_cardinality_defect={sum(f4.values())-sum(j4.values())})")
    print("scope=coefficient/support test only; no action on rails, semantic owners, amplitudes, or endpoint currents")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
