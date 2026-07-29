#!/usr/bin/env python3
"""Exact terminal sieve for THM-2928's relaxed two-drift address frontier.

The preceding address-mask referee proves that every denominator pair with
min(d_1,d_2)<=1000 is impossible and leaves exactly 942 pairs with both
denominators greater than 1000.  This independent pass rederives that
942-pair universe from the body-support census and applies two further
necessary conditions.

For a body address word S_D and divisor d, write

    lambda_d(a)=#{s in S_D:s=a mod d},   k_d=ceil(d/7).

An actual d-mask selects at most k_d residue classes.  If those classes are
allowed to be completely arbitrary, its load is at most the sum of the k_d
largest lambda_d(a).  Giving the other mask its full cardinality capacity
(D/d_2)k_(d_2) is an additional relaxation.  This top-class screen kills
940 of the 942 pairs and leaves only d_1=d_2=D on each of two body rows.

There is a one-fiber exact obstruction for those diagonal pairs.  Write
D=7k.  Every enlarged actual mask has the form

    B(a,u)={a+u*j mod D:0<=j<k},   gcd(u,D)=1.

Reduction modulo k makes B(a,u) a transversal because gcd(u,k)=1.  Hence
two such masks contain at most two elements of each fiber of
Z/DZ -> Z/kZ.  Each of the two remaining S_D words has a fiber of
multiplicity at least three, so neither is coverable.

The support and preceding address referees are hash-pinned.  All arithmetic
is integral or Fraction-exact, and every check remains active under
optimized Python.
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
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
ADDRESS_PATH = HERE / "lrc14_two_drift_relaxed_address_mask_thm2928.py"
ADDRESS_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_two_drift_relaxed_address_mask_thm2928.out"
)
EXPECTED_SUPPORT_SHA256 = (
    "778842c0e8e7172835ca6ae673fb6156f212d4296e672bce4e7cc2815195bf1a"
)
EXPECTED_SUPPORT_OUTPUT_SHA256 = (
    "648327d3b9b5b9a50c7760f0afd89a7a33161f57fa98c1b9e181d6b5b791a25f"
)
EXPECTED_ADDRESS_SHA256 = (
    "870498c4f0a2d97a2d42bce593c44283c77a141fb08669b4a91133e39db5c276"
)
EXPECTED_ADDRESS_OUTPUT_SHA256 = (
    "74f7c270034dc40b4de3d33b9abf67481435d1e97eb6e52f2448d0a152cb68d7"
)
EXPECTED_SEMANTIC_SHA256 = (
    "eea51e491a52776bf56fb0747cf507c94f8e34fd86c15c326815de72f35a7d91"
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
require(
    file_sha256(ADDRESS_PATH) == EXPECTED_ADDRESS_SHA256,
    "address-script dependency hash changed",
)
require(
    file_sha256(ADDRESS_OUTPUT_PATH) == EXPECTED_ADDRESS_OUTPUT_SHA256,
    "address-output dependency hash changed",
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
    return tuple(merged)


def residue_loads(arcs, d):
    """Return lambda_d using cyclic range additions."""
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
    require(
        sum(loads) == sum(right - left for left, right in arcs),
        "residue load lost support mass",
    )
    return tuple(loads)


def point_range(point, ranges):
    for left, right in ranges:
        if left <= point < right:
            return (left, right)
    raise RuntimeError(("witness point is not in the support", point))


def qtext(value):
    return f"{value.numerator}/{value.denominator}"


def main():
    denominator_cache = {}
    arc_cache = {}
    load_cache = {}
    semantic_hash = sha256()

    body_count = 0
    divisor_rows = 0
    support_hard_rows = 0
    raw_denominator_pairs = 0
    cardinality_survivors = 0
    cardinality_survivor_rows = 0
    large_denominator_pairs = 0
    large_denominator_rows = set()
    arbitrary_class_kills = 0
    arbitrary_class_survivors = []
    closest_kill_by_row = {}

    for F in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support_module.safe_cell_ranges(F)
        for D in support_module.divisors(L):
            divisor_rows += 1
            support_count = support_module.support_size_bitset(D, ranges)
            density = Q(support_count, D)
            if density > support_module.SUPPORT_CUTOFF:
                continue
            support_hard_rows += 1
            if D not in denominator_cache:
                denominator_cache[D] = support_module.denominator_pairs(D)
            pairs = denominator_cache[D]
            raw_denominator_pairs += len(pairs)
            kept_on_row = False
            for d1, d2 in pairs:
                k1 = (d1 + 6) // 7
                k2 = (d2 + 6) // 7
                cardinality_cap = Q(k1, d1) + Q(k2, d2)
                if density > cardinality_cap:
                    continue
                cardinality_survivors += 1
                kept_on_row = True
                if d1 <= SMALL_DENOMINATOR_LIMIT:
                    continue

                large_denominator_pairs += 1
                row_key = (F, L, D, support_count, density)
                large_denominator_rows.add(row_key)
                arc_key = (F, D)
                if arc_key not in arc_cache:
                    arcs = projected_support_arcs(D, ranges)
                    require(
                        sum(right - left for left, right in arcs)
                        == support_count,
                        ("support arc mismatch", F, D),
                    )
                    arc_cache[arc_key] = arcs
                arcs = arc_cache[arc_key]
                load_key = (F, D, d1)
                if load_key not in load_cache:
                    loads = residue_loads(arcs, d1)
                    load_cache[load_key] = sum(
                        sorted(loads, reverse=True)[:k1]
                    )
                first_top_class_load = load_cache[load_key]
                second_full_capacity = (D // d2) * k2
                slack = (
                    first_top_class_load
                    + second_full_capacity
                    - support_count
                )
                record = (
                    F,
                    L,
                    D,
                    support_count,
                    density,
                    d1,
                    d2,
                    k1,
                    first_top_class_load,
                    second_full_capacity,
                    slack,
                )
                if slack < 0:
                    arbitrary_class_kills += 1
                    previous = closest_kill_by_row.get(row_key)
                    if previous is None or slack > previous[-1]:
                        closest_kill_by_row[row_key] = record
                    outcome = "top-class-kill"
                else:
                    arbitrary_class_survivors.append(record)
                    outcome = "top-class-survivor"
                semantic_hash.update(
                    (
                        f"{','.join(map(str, F))}|{L}|{D}|{support_count}|"
                        f"{d1}|{d2}|{k1}|{first_top_class_load}|"
                        f"{second_full_capacity}|{slack}|{outcome}\n"
                    ).encode()
                )
            cardinality_survivor_rows += kept_on_row

    expected_large_rows = {
        (
            (1, 4, 5, 7, 9, 11),
            194040,
            194040,
            55392,
            Q(2308, 8085),
        ),
        (
            (1, 5, 7, 8, 9, 11),
            388080,
            388080,
            109044,
            Q(3029, 10780),
        ),
    }
    expected_arbitrary_survivors = [
        (
            (1, 4, 5, 7, 9, 11),
            194040,
            194040,
            55392,
            Q(2308, 8085),
            194040,
            194040,
            27720,
            27720,
            27720,
            48,
        ),
        (
            (1, 5, 7, 8, 9, 11),
            388080,
            388080,
            109044,
            Q(3029, 10780),
            388080,
            388080,
            55440,
            55440,
            55440,
            1836,
        ),
    ]
    require(body_count == 3003, "body universe changed")
    require(divisor_rows == 251536, "divisor-row universe changed")
    require(support_hard_rows == 10976, "support frontier changed")
    require(raw_denominator_pairs == 3066274, "raw pair universe changed")
    require(cardinality_survivors == 23755, "cardinality frontier changed")
    require(
        cardinality_survivor_rows == 6292,
        "cardinality survivor-row count changed",
    )
    require(large_denominator_pairs == 942, "large pair universe changed")
    require(large_denominator_rows == expected_large_rows, "large rows changed")
    require(arbitrary_class_kills == 940, "top-class kill count changed")
    require(
        arbitrary_class_survivors == expected_arbitrary_survivors,
        "top-class survivors changed",
    )

    closest_kills = sorted(closest_kill_by_row.values())
    require(
        [
            (
                row[0],
                row[2],
                row[5],
                row[6],
                row[8],
                row[9],
                row[10],
            )
            for row in closest_kills
        ]
        == [
            (
                (1, 4, 5, 7, 9, 11),
                194040,
                97020,
                194040,
                21620,
                27720,
                -6052,
            ),
            (
                (1, 5, 7, 8, 9, 11),
                388080,
                194040,
                388080,
                40510,
                55440,
                -13094,
            ),
        ],
        "closest top-class hostile controls changed",
    )

    fiber_records = []
    for survivor in arbitrary_class_survivors:
        (
            F,
            L,
            D,
            support_count,
            _density,
            d1,
            d2,
            k1,
            _first_load,
            _second_capacity,
            _slack,
        ) = survivor
        require(L == D == d1 == d2 == 7 * k1, "non-diagonal survivor")
        arcs = arc_cache[(F, D)]
        body_L, body_ranges = support_module.safe_cell_ranges(F)
        require(body_L == L, "fiber body ruler changed")
        fiber_loads = residue_loads(arcs, k1)
        maximum = max(fiber_loads)
        # Move one cell into the first maximal plateau, avoiding any need to
        # discuss the harmless danger-interval boundary convention.
        witness_residue = fiber_loads.index(maximum) + 1
        require(
            fiber_loads[witness_residue] == maximum,
            "maximal fiber plateau lost its interior witness",
        )
        witness_points = tuple(
            witness_residue + lift * k1
            for lift in range(7)
            if any(
                left <= witness_residue + lift * k1 < right
                for left, right in arcs
            )
        )
        witness_ranges = tuple(
            point_range(point, body_ranges) for point in witness_points
        )
        histogram_counter = Counter(fiber_loads)
        histogram = tuple(
            (load, histogram_counter.get(load, 0))
            for load in range(maximum + 1)
        )
        require(sum(count for _load, count in histogram) == k1, "fiber loss")
        require(
            sum(load * count for load, count in histogram) == support_count,
            "fiber mass changed",
        )
        require(
            maximum > 2,
            "two transversals were not obstructed by a fiber",
        )
        fiber_records.append(
            (
                F,
                D,
                k1,
                histogram,
                witness_residue,
                witness_points,
                witness_ranges,
            )
        )

    expected_fiber_records = [
        (
            (1, 4, 5, 7, 9, 11),
            194040,
            27720,
            ((0, 3960), (1, 0), (2, 17472), (3, 4704), (4, 1584)),
            1981,
            (29701, 57421, 112861, 168301),
            (
                (29700, 34020),
                (57420, 63140),
                (112860, 113652),
                (168300, 170940),
            ),
        ),
        (
            (1, 5, 7, 8, 9, 11),
            388080,
            55440,
            ((0, 7920), (1, 0), (2, 33516), (3, 14004)),
            3961,
            (59401, 114841, 225721),
            (
                (59400, 68040),
                (114840, 126280),
                (225720, 227304),
            ),
        ),
    ]
    require(fiber_records == expected_fiber_records, "fiber witnesses changed")
    semantic_digest = semantic_hash.hexdigest()
    require(
        semantic_digest == EXPECTED_SEMANTIC_SHA256,
        "large-pair semantic ledger changed",
    )

    print("LRC14 two-drift arbitrary-class and seventh-fiber exact referee")
    print(f"support_script_sha256={file_sha256(SUPPORT_PATH)}")
    print(f"support_output_sha256={file_sha256(SUPPORT_OUTPUT_PATH)}")
    print(f"address_script_sha256={file_sha256(ADDRESS_PATH)}")
    print(f"address_output_sha256={file_sha256(ADDRESS_OUTPUT_PATH)}")
    print(f"body_count={body_count}")
    print(f"divisor_rows={divisor_rows}")
    print(f"support_hard_rows={support_hard_rows}")
    print(f"raw_denominator_pairs={raw_denominator_pairs}")
    print(f"cardinality_survivors={cardinality_survivors}")
    print(f"cardinality_survivor_rows={cardinality_survivor_rows}")
    print(f"large_denominator_limit={SMALL_DENOMINATOR_LIMIT}")
    print(f"large_denominator_pairs={large_denominator_pairs}")
    print(f"large_denominator_rows={sorted(large_denominator_rows)}")
    print(
        "top_class_formula=max_load_d<="
        "sum_(j=1)^ceil(d/7) lambda_d^[j]"
    )
    print(f"arbitrary_class_kills={arbitrary_class_kills}")
    print(f"arbitrary_class_survivors={arbitrary_class_survivors}")
    print(f"closest_top_class_kills={closest_kills}")
    print(
        "seventh_fiber_lemma=D=7k and gcd(u,D)=1 imply "
        "{a+uj:0<=j<k} is a transversal modulo k"
    )
    print(f"fiber_records={fiber_records}")
    print("fiber_kills=2")
    print("final_large_denominator_survivors=0")
    print(f"semantic_sha256={semantic_digest}")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
