#!/usr/bin/env python3
"""Exact controls for THM-3405's common-centre Boolean twist theorem.

The universal proof is arithmetic.  This standard-library companion checks it
against the independent selected-mode formulas from THM-3398, exhausts small
zero-cochain mode tuples, reconstructs literal danger sets with exact
fractions, verifies the translation quotient, and freezes the q=8,16,23
boundary controls.  No floating point, randomness, or assertion-dependent
truth gate is used.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import gcd
from pathlib import Path


EXPECTED_SEMANTIC_DIGEST = "8cf22bc02893c04f54c4ac7d70121b26b94e71fc0a596604421a42421fb54ab6"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


class ExactDigest:
    def __init__(self):
        self.value = sha256()

    def add(self, item):
        self.value.update(repr(item).encode("ascii"))
        self.value.update(bytes((10,)))

    def hexdigest(self):
        return self.value.hexdigest()


def gcd_many(values):
    answer = 0
    for value in values:
        answer = gcd(answer, value)
    return answer


def extended_gcd(left, right):
    if right == 0:
        return left, 1, 0
    divisor, x1, y1 = extended_gcd(right, left % right)
    return divisor, y1, x1 - (left // right) * y1


def bezout_one(values):
    """Return coefficients whose dot product with coprime values is one."""
    coefficients = [1]
    current = values[0]
    for value in values[1:]:
        divisor, left, right = extended_gcd(current, value)
        coefficients = [left * coefficient for coefficient in coefficients] + [right]
        current = divisor
    require(current == 1, (values, current))
    require(sum(c * v for c, v in zip(coefficients, values)) == 1, values)
    return tuple(coefficients)


def owner_modes(q, speed):
    """Independent literal implementation of THM-3398 mode formula."""
    common = gcd(q, speed)
    phase_count = q // common
    phase_step = (speed // common) % phase_count
    modes = []
    for size in range(1, (phase_count + 6) // 7 + 1):
        width = q - 7 * common * (size - 1)
        require(width > 0, (q, speed, size, width))
        for start in range(phase_count):
            phases = {(start + offset) % phase_count for offset in range(size)}
            block = frozenset(
                sheet
                for sheet in range(q)
                if (phase_step * sheet) % phase_count in phases
            )
            centre = (-common * (2 * start + size - 1)) % (2 * q)
            modes.append((block, centre, width, start, size, common, phase_count))
    return tuple(modes)


def source_danger(q, speed, centre):
    answer = []
    for sheet in range(q):
        value = (speed * (centre + Fraction(sheet, q))) % 1
        if 14 * min(value, 1 - value) < 1:
            answer.append(sheet)
    return frozenset(answer)


def scalar_witness(q, speeds, modes):
    """Return the unique a mod 2q if all selected mode centres coincide."""
    common_speed = gcd_many(speeds)
    reduced = tuple(speed // common_speed for speed in speeds)
    coefficients = bezout_one(reduced)
    modulus = 2 * q
    candidate = sum(coefficient * mode[1] for coefficient, mode in zip(coefficients, modes)) % modulus
    if all((candidate * value - mode[1]) % modulus == 0 for value, mode in zip(reduced, modes)):
        return candidate
    return None


def pairwise_zero_compatible(q, speeds, modes):
    modulus_factor = 2 * q
    for left in range(len(speeds)):
        for right in range(left + 1, len(speeds)):
            residue = modes[left][1] * speeds[right] - modes[right][1] * speeds[left]
            modulus = modulus_factor * gcd(speeds[left], speeds[right])
            if residue % modulus:
                return False
    return True


def mode_scalar_audit(digest):
    checked = 0
    positive = 0
    cover_positive = 0
    nonboolean = 0
    for q in range(2, 11):
        speed_pool = tuple(range(1, min(8, 2 * q)))
        for rank in (1, 2, 3):
            for speeds in combinations(speed_pool, rank):
                banks = tuple(owner_modes(q, speed) for speed in speeds)
                # Full pairs, and a deterministic but broad triple slice.
                tuples = product(*banks)
                for index, modes in enumerate(tuples):
                    if rank == 3 and index >= 1200:
                        break
                    scalar = scalar_witness(q, speeds, modes)
                    pairwise = pairwise_zero_compatible(q, speeds, modes)
                    require((scalar is not None) == pairwise, (q, speeds, modes, scalar, pairwise))
                    checked += 1
                    if scalar is None:
                        continue
                    positive += 1
                    common_speed = gcd_many(speeds)
                    reduced = tuple(speed // common_speed for speed in speeds)
                    gauge_gcd = gcd(q, common_speed)
                    require(scalar % gauge_gcd == 0, (q, speeds, scalar, gauge_gcd))
                    twist = scalar % (2 * gauge_gcd)
                    require(twist in (0, gauge_gcd), (q, speeds, scalar, twist, gauge_gcd))
                    if twist not in (0, gauge_gcd):
                        nonboolean += 1
                    centre = Fraction(scalar, 2 * q * common_speed)
                    for speed, value, mode in zip(speeds, reduced, modes):
                        require((speed * centre - Fraction(mode[1], 2 * q)).denominator == 1,
                                (q, speeds, scalar, speed, mode))
                        require(mode[0] <= source_danger(q, speed, centre),
                                (q, speed, centre, mode[0], source_danger(q, speed, centre)))
                        require(mode[1] % mode[5] == 0, (q, speed, mode))
                        require((scalar * value - mode[1]) % (2 * q) == 0,
                                (q, speeds, scalar, value, mode))
                    if frozenset().union(*(mode[0] for mode in modes)) == frozenset(range(q)):
                        cover_positive += 1
                    digest.add((q, speeds, tuple((mode[1], mode[3], mode[4]) for mode in modes), scalar))
    require(nonboolean == 0, nonboolean)
    return checked, positive, cover_positive


def orbit_audit(digest):
    checked = 0
    for q in range(2, 43):
        modulus = 2 * q
        for common_speed in range(1, 2 * q + 1):
            gauge_gcd = gcd(q, common_speed)
            for scalar in range(modulus):
                orbit = frozenset((scalar + 2 * common_speed * shift) % modulus for shift in range(q))
                predicted = frozenset(
                    other
                    for other in range(modulus)
                    if other % (2 * gauge_gcd) == scalar % (2 * gauge_gcd)
                )
                require(orbit == predicted, (q, common_speed, scalar, orbit, predicted))
                checked += 1
            digest.add((q, common_speed, gauge_gcd))
    return checked


def literal_cover(q, speeds, centre):
    return frozenset().union(*(source_danger(q, speed, centre) for speed in speeds))


def canonical_layers(q, speeds):
    common_speed = gcd_many(speeds)
    gauge_gcd = gcd(q, common_speed)
    return (
        Fraction(0),
        Fraction(gauge_gcd, 2 * q * common_speed),
    )


def bounded_cover_audit(digest):
    checked = 0
    family_positives = 0
    for q in range(2, 15):
        full = frozenset(range(q))
        pool = tuple(range(1, min(11, 2 * q)))
        for rank in range(1, 5):
            for speeds in combinations(pool, rank):
                common_speed = gcd_many(speeds)
                gauge_gcd = gcd(q, common_speed)
                brute_family_scalars = []
                for scalar in range(2 * q):
                    centre = Fraction(scalar, 2 * q * common_speed)
                    blocks = tuple(source_danger(q, speed, centre) for speed in speeds)
                    if any(not block for block in blocks):
                        continue
                    # Every named owner is active in this audit.
                    centred = True
                    for speed, value, block in zip(
                        speeds,
                        (speed // common_speed for speed in speeds),
                        blocks,
                    ):
                        target_h = scalar * value % (2 * q)
                        if not any(mode[0] == block and mode[1] == target_h for mode in owner_modes(q, speed)):
                            centred = False
                            break
                    if centred and frozenset().union(*blocks) == full:
                        brute_family_scalars.append(scalar)
                layers = canonical_layers(q, speeds)
                layer_blocks = tuple(
                    tuple(source_danger(q, speed, centre) for speed in speeds)
                    for centre in layers
                )
                predicted_family = any(
                    all(blocks) and frozenset().union(*blocks) == full
                    for blocks in layer_blocks
                )
                brute_family = bool(brute_family_scalars)
                require(brute_family == predicted_family,
                        ("family", q, speeds, brute_family_scalars, layers))
                require(
                    all(scalar % (2 * gauge_gcd) in (0, gauge_gcd)
                        for scalar in brute_family_scalars),
                    (q, speeds, brute_family_scalars, gauge_gcd),
                )
                checked += 1
                family_positives += brute_family
                if brute_family:
                    digest.add((q, speeds, tuple(brute_family_scalars)))
    return checked, family_positives


def dilation_audit(digest):
    bases = (
        (8, (1, 3, 5, 7), Fraction(-1, 16)),
        (9, (1, 5, 6, 7), Fraction(5, 6)),
        (15, (1, 2, 3, 4, 5, 7), Fraction(0)),
    )
    checked = 0
    for q, speeds, centre in bases:
        base_d = gcd_many(speeds)
        base_g = gcd(q, base_d)
        base_a = 2 * q * base_d * centre
        require(base_a.denominator == 1, (q, speeds, centre, base_a))
        base_twist = int(base_a) % (2 * base_g)
        for scale in range(1, 8):
            q2 = scale * q
            speeds2 = tuple(scale * speed for speed in speeds)
            centre2 = centre / scale
            d2 = gcd_many(speeds2)
            g2 = gcd(q2, d2)
            a2 = 2 * q2 * d2 * centre2
            require(a2 == scale * base_a, (q, scale, a2, base_a))
            require(int(a2) % (2 * g2) == scale * base_twist % (2 * g2),
                    (q, scale, a2, g2, base_twist))
            old_blocks = tuple(source_danger(q, speed, centre) for speed in speeds)
            new_blocks = tuple(source_danger(q2, speed, centre2) for speed in speeds2)
            for old, new in zip(old_blocks, new_blocks):
                require(new == frozenset(sheet for sheet in range(q2) if sheet % q in old),
                        (q, scale, old, new))
            checked += 1
            digest.add((q, scale, int(a2), g2))
    return checked


def frozen_controls(digest):
    controls = (
        (8, (1, 3, 5, 7), Fraction(1, 16), 1),
        (16, (2, 6, 10, 14), Fraction(1, 32), 2),
        (23, (1, 4, 5, 7, 9, 11), Fraction(1, 46), 1),
        (15, (1, 2, 3, 4, 5, 7), Fraction(0), 0),
    )
    records = []
    for q, speeds, centre, expected_twist in controls:
        common_speed = gcd_many(speeds)
        gauge_gcd = gcd(q, common_speed)
        scalar = 2 * q * common_speed * centre
        require(scalar.denominator == 1, (q, speeds, centre, scalar))
        twist = int(scalar) % (2 * gauge_gcd)
        require(twist == expected_twist, (q, speeds, twist, expected_twist))
        blocks = tuple(tuple(sorted(source_danger(q, speed, centre))) for speed in speeds)
        require(frozenset().union(*(frozenset(block) for block in blocks)) == frozenset(range(q)),
                (q, speeds, centre, blocks))
        records.append((q, speeds, centre, common_speed, gauge_gcd, int(scalar), twist, blocks))
    digest.add(tuple(records))
    hostile_q = 8
    hostile_speeds = (1, 2, 3, 5, 7)
    hostile_centre = Fraction(1, 16)
    hostile_blocks = tuple(
        tuple(sorted(source_danger(hostile_q, speed, hostile_centre)))
        for speed in hostile_speeds
    )
    require(hostile_blocks[1] == (), hostile_blocks)
    require(
        frozenset().union(*(frozenset(block) for block in hostile_blocks))
        == frozenset(range(hostile_q)),
        hostile_blocks,
    )
    empty_owner_hostile = (hostile_q, hostile_speeds, hostile_centre, hostile_blocks)
    digest.add(empty_owner_hostile)

    subset_q = 9
    active_speeds = (2, 10, 12, 14)
    ambient_speeds = (1, 2, 10, 12, 14)
    active_centre = Fraction(1, 36)
    active_blocks = tuple(
        tuple(sorted(source_danger(subset_q, speed, active_centre)))
        for speed in active_speeds
    )
    require(
        all(active_blocks)
        and frozenset().union(*(frozenset(block) for block in active_blocks))
        == frozenset(range(subset_q)),
        active_blocks,
    )
    ambient_layers = canonical_layers(subset_q, ambient_speeds)
    ambient_unions = tuple(
        tuple(sorted(literal_cover(subset_q, ambient_speeds, centre)))
        for centre in ambient_layers
    )
    require(all(len(union) < subset_q for union in ambient_unions), ambient_unions)
    subset_gcd_hostile = (
        subset_q,
        active_speeds,
        active_centre,
        active_blocks,
        ambient_speeds,
        ambient_layers,
        ambient_unions,
    )
    digest.add(subset_gcd_hostile)
    return tuple(records), (empty_owner_hostile, subset_gcd_hostile)


def main():
    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "float literal",
    )

    digest = ExactDigest()
    mode_checks, mode_positive, mode_covers = mode_scalar_audit(digest)
    orbit_checks = orbit_audit(digest)
    cover_checks, family_positive = bounded_cover_audit(digest)
    dilation_checks = dilation_audit(digest)
    controls, scope_hostiles = frozen_controls(digest)
    semantic = digest.hexdigest()
    if EXPECTED_SEMANTIC_DIGEST != "TO_BE_FILLED":
        require(semantic == EXPECTED_SEMANTIC_DIGEST, (semantic, EXPECTED_SEMANTIC_DIGEST))

    print("COMMON-CENTRE GCD GAUGE AND BOOLEAN TWIST THM-3405")
    print(f"source_sha256_lf={lf_hash(source)}")
    print("status=PROVED_analytic_plus_VERIFIED-EXACT_controls;THM-3405")
    print(f"mode_tuple_checks={mode_checks};zero_compatible={mode_positive};block_covers={mode_covers}")
    print(f"translation_orbit_checks={orbit_checks}")
    print(f"bounded_active_family_checks={cover_checks};positive={family_positive}")
    print(f"dilation_checks={dilation_checks}")
    print(f"frozen_controls={controls}")
    print(f"scope_hostiles={scope_hostiles}")
    print("raw_gauge=Z/(2*gcd(q,gcd(U)));selected_mode_divisibility_reduces_to_two_orbits")
    print("canonical_layers=c0=0;c1=gcd(q,d)/(2*q*d);fixed_zero_vs_order_two")
    print("scope=zero_complete_cochain_only;not_arbitrary_physical_cover;no_LRC14_decrement")
    print(f"semantic_sha256={semantic}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
