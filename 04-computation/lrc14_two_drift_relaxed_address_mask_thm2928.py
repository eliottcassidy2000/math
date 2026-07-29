#!/usr/bin/env python3
"""Relaxed two-needle address obstruction after the THM-2928 support census.

Fix a support-hard literal body quotient (F,D,S_D) and denominator pair
d_1<=d_2 with lcm(d_1,d_2)=D.  At a fixed normalized phase, a d-denominator
comb selects lifts of at most ceil(d/7) consecutive residues after
multiplication by a unit modulo d.  We deliberately relax the real problem:

* the two cyclic blocks may choose their phases independently;
* each block is enlarged to ceil(d/7) residues;
* only the body address word S_D is retained, not the five-comb shape R_A.

For the first denominator we maximize its load on S_D over every unit
direction and every cyclic block shift.  The second comb is initially given
its entire cardinality capacity.  If this already underfills S_D, coverage is
impossible.

Only two rows with min(d_1,d_2)<=1000 survive that load relaxation.  Both
have d_1=2.  For them we retain each of the two parity masks and project the
residual P modulo d_2.  Containment in a unit image of a k-term cyclic block
would imply the NECESSARY (not sufficient) condition

    |P-P| <= 2k-1.

Both parity residuals violate it.  Thus all 22,813 of the 23,755
denominator-cardinality survivors with min(d_1,d_2)<=1000 are impossible.
Exactly 942 denominator-pair occurrences with both denominators >1000
remain, on two body words.  This is a numerator-free necessary screen, not
uniform two-drift closure.

The script pins and imports the exact body-support census on which it rests.
All checks are explicit RuntimeError checks and remain active under python -O.
"""

from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SUPPORT_PATH = HERE / "lrc14_two_drift_body_projection_support_thm2928.py"
SUPPORT_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_two_drift_body_projection_support_thm2928.out"
)
EXPECTED_SUPPORT_SHA256 = (
    "778842c0e8e7172835ca6ae673fb6156f212d4296e672bce4e7cc2815195bf1a"
)
EXPECTED_SUPPORT_OUTPUT_SHA256 = (
    "648327d3b9b5b9a50c7760f0afd89a7a33161f57fa98c1b9e181d6b5b791a25f"
)
SMALL_DENOMINATOR_LIMIT = 1000


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


require(
    file_sha256(SUPPORT_PATH) == EXPECTED_SUPPORT_SHA256,
    "support-script dependency hash changed",
)
require(
    file_sha256(SUPPORT_OUTPUT_PATH) == EXPECTED_SUPPORT_OUTPUT_SHA256,
    "support-output dependency hash changed",
)
spec = spec_from_file_location("lrc14_body_projection_support", SUPPORT_PATH)
support_module = module_from_spec(spec)
spec.loader.exec_module(support_module)


def projected_support_arcs(D, ranges):
    """Merge the image modulo D of half-open integer ranges."""
    pieces = []
    for left, right in ranges:
        length = right - left
        if length >= D:
            return [(0, D)]
        residue = left % D
        endpoint = residue + length
        if endpoint <= D:
            pieces.append((residue, endpoint))
        else:
            pieces.append((residue, D))
            pieces.append((0, endpoint - D))
    pieces.sort()
    merged = []
    for left, right in pieces:
        if merged and left <= merged[-1][1]:
            if right > merged[-1][1]:
                merged[-1] = (merged[-1][0], right)
        else:
            merged.append((left, right))
    return merged


def support_mask(D, arcs):
    mask = 0
    for left, right in arcs:
        mask |= ((1 << (right - left)) - 1) << left
    return mask


def residue_loads(arcs, d):
    """Counts of S_D in the d residue classes, using cyclic range additions."""
    common = 0
    difference = [0] * (d + 1)
    for left, right in arcs:
        quotient, remainder = divmod(right - left, d)
        common += quotient
        start = left % d
        endpoint = start + remainder
        if endpoint <= d:
            difference[start] += 1
            difference[endpoint] -= 1
        else:
            difference[start] += 1
            difference[d] -= 1
            difference[0] += 1
            difference[endpoint - d] -= 1
    loads = []
    running = common
    for residue in range(d):
        running += difference[residue]
        loads.append(running)
    require(sum(loads) == sum(b - a for a, b in arcs), "residue load lost mass")
    return loads


def maximum_unit_block_load(arcs, d):
    """Maximum S_D load of a ceil(d/7)-term cyclic unit-step block."""
    loads = residue_loads(arcs, d)
    width = (d + 6) // 7
    best = 0
    for multiplier in range(1, d):
        if gcd(multiplier, d) != 1:
            continue
        inverse = pow(multiplier, -1, d)
        word = [loads[(inverse * j) % d] for j in range(d)]
        current = sum(word[:width])
        local_best = current
        for start in range(1, d):
            current += word[(start + width - 1) % d] - word[start - 1]
            local_best = max(local_best, current)
        best = max(best, local_best)
    return best


def base_unit_block_masks(d):
    """All enlarged unit-direction cyclic blocks on Z/dZ."""
    width = (d + 6) // 7
    masks = set()
    for multiplier in range(1, d):
        if gcd(multiplier, d) != 1:
            continue
        inverse = pow(multiplier, -1, d)
        mask = sum(1 << ((inverse * j) % d) for j in range(width))
        for start in range(d):
            masks.add(mask)
            mask &= ~(1 << ((inverse * start) % d))
            mask |= 1 << ((inverse * ((start + width) % d)) % d)
    return tuple(sorted(masks))


def lift_mask(mask, d, D):
    repeat = ((1 << D) - 1) // ((1 << d) - 1)
    return mask * repeat


def project_mask(mask, D, d):
    low = (1 << d) - 1
    projected = 0
    for block in range(D // d):
        projected |= (mask >> (block * d)) & low
    return projected


def rotate_right(mask, shift, d):
    shift %= d
    if shift == 0:
        return mask
    low = (1 << d) - 1
    return ((mask >> shift) | (mask << (d - shift))) & low


def difference_excess_witness(mask, d, limit):
    """Return an exact partial |P-P| lower bound once it exceeds limit."""
    differences = 0
    remaining = mask
    shifts = 0
    while remaining and differences.bit_count() <= limit:
        bit = remaining & -remaining
        residue = bit.bit_length() - 1
        differences |= rotate_right(mask, residue, d)
        remaining -= bit
        shifts += 1
    return differences.bit_count(), shifts


def main():
    denominator_cache = {}
    load_cache = {}
    block_mask_cache = {}
    semantic_hash = sha256()

    raw_denominator_pairs = 0
    cardinality_survivors = 0
    small_denominator_pairs = 0
    load_killed = 0
    load_survivors = []
    difference_killed = 0
    difference_witnesses = []
    small_denominator_survivors = 0
    large_denominator_remainder = 0
    remaining_by_row = {}

    for F in combinations(range(1, 15), 6):
        L, ranges = support_module.safe_cell_ranges(F)
        for D in support_module.divisors(L):
            support_count = support_module.support_size_bitset(D, ranges)
            density = Q(support_count, D)
            if density > support_module.SUPPORT_CUTOFF:
                continue
            arcs = projected_support_arcs(D, ranges)
            require(
                sum(right - left for left, right in arcs) == support_count,
                ("support arc mismatch", F, D),
            )
            if D not in denominator_cache:
                denominator_cache[D] = support_module.denominator_pairs(D)
            for d1, d2 in denominator_cache[D]:
                raw_denominator_pairs += 1
                cardinality_cap = Q((d1 + 6) // 7, d1) + Q(
                    (d2 + 6) // 7, d2
                )
                if density > cardinality_cap:
                    continue
                cardinality_survivors += 1
                outcome = "large"
                sidecar = ""

                if d1 > SMALL_DENOMINATOR_LIMIT:
                    large_denominator_remainder += 1
                    key = (F, L, D, support_count, density)
                    remaining_by_row[key] = remaining_by_row.get(key, 0) + 1
                else:
                    small_denominator_pairs += 1
                    load_key = (F, D, d1)
                    if load_key not in load_cache:
                        load_cache[load_key] = maximum_unit_block_load(arcs, d1)
                    first_load = load_cache[load_key]
                    second_capacity = (D // d2) * ((d2 + 6) // 7)
                    sidecar = f"{first_load},{second_capacity}"
                    if support_count > first_load + second_capacity:
                        load_killed += 1
                        outcome = "load-kill"
                    else:
                        load_survivors.append(
                            (
                                F,
                                L,
                                D,
                                support_count,
                                density,
                                d1,
                                d2,
                                first_load,
                                second_capacity,
                            )
                        )
                        if d1 not in block_mask_cache:
                            block_mask_cache[d1] = base_unit_block_masks(d1)
                        S = support_mask(D, arcs)
                        pair_survives = False
                        for first_base_mask in block_mask_cache[d1]:
                            residual = S & ~lift_mask(first_base_mask, d1, D)
                            projected = project_mask(residual, D, d2)
                            second_width = (d2 + 6) // 7
                            if projected.bit_count() > second_width:
                                continue
                            difference_lower, shifts = difference_excess_witness(
                                projected, d2, 2 * second_width - 1
                            )
                            difference_witnesses.append(
                                (
                                    F,
                                    D,
                                    d1,
                                    d2,
                                    first_base_mask,
                                    projected.bit_count(),
                                    difference_lower,
                                    2 * second_width - 1,
                                    shifts,
                                )
                            )
                            if difference_lower <= 2 * second_width - 1:
                                pair_survives = True
                                break
                        if pair_survives:
                            small_denominator_survivors += 1
                            outcome = "unresolved-small"
                        else:
                            difference_killed += 1
                            outcome = "difference-kill"

                semantic_hash.update(
                    (
                        f"{','.join(map(str, F))}|{L}|{D}|{support_count}|"
                        f"{d1}|{d2}|{outcome}|{sidecar}\n"
                    ).encode()
                )

    remaining_rows = [
        (*key, count)
        for key, count in sorted(
            remaining_by_row.items(), key=lambda item: item[0]
        )
    ]
    require(raw_denominator_pairs == 3066274, "raw pair count changed")
    require(cardinality_survivors == 23755, "cardinality frontier changed")
    require(small_denominator_pairs == 22813, "small-denominator count changed")
    require(load_killed == 22811, "load-kill count changed")
    require(len(load_survivors) == 2, "load-survivor count changed")
    require(difference_killed == 2, "difference-kill count changed")
    require(small_denominator_survivors == 0, "small denominator survived")
    require(large_denominator_remainder == 942, "large remainder changed")
    require(
        remaining_rows
        == [
            (
                (1, 4, 5, 7, 9, 11),
                194040,
                194040,
                55392,
                Q(2308, 8085),
                371,
            ),
            (
                (1, 5, 7, 8, 9, 11),
                388080,
                388080,
                109044,
                Q(3029, 10780),
                571,
            ),
        ],
        ("remaining row ledger changed", remaining_rows),
    )
    require(
        difference_witnesses
        == [
            (
                (1, 4, 5, 7, 9, 11),
                194040,
                2,
                194040,
                1,
                27696,
                63905,
                55439,
                1261,
            ),
            (
                (1, 4, 5, 7, 9, 11),
                194040,
                2,
                194040,
                2,
                27696,
                63905,
                55439,
                1261,
            ),
            (
                (1, 5, 7, 8, 9, 11),
                388080,
                2,
                388080,
                1,
                54522,
                129569,
                110879,
                2521,
            ),
            (
                (1, 5, 7, 8, 9, 11),
                388080,
                2,
                388080,
                2,
                54522,
                129570,
                110879,
                2521,
            ),
        ],
        "difference witnesses changed",
    )

    print(f"support_script_sha256={EXPECTED_SUPPORT_SHA256}")
    print(f"support_output_sha256={EXPECTED_SUPPORT_OUTPUT_SHA256}")
    print(f"small_denominator_limit={SMALL_DENOMINATOR_LIMIT}")
    print(f"raw_denominator_pairs={raw_denominator_pairs}")
    print(f"cardinality_survivors={cardinality_survivors}")
    print(f"small_denominator_pairs={small_denominator_pairs}")
    print(f"load_killed={load_killed}")
    print(f"load_survivors={load_survivors}")
    print(f"difference_killed={difference_killed}")
    print(f"difference_witnesses={difference_witnesses}")
    print(f"small_denominator_survivors={small_denominator_survivors}")
    print(f"large_denominator_remainder={large_denominator_remainder}")
    print(f"remaining_rows={remaining_rows}")
    print(f"semantic_sha256={semantic_hash.hexdigest()}")


if __name__ == "__main__":
    main()
