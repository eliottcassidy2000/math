#!/usr/bin/env python3
"""Exact support-transfer ladder for the six-body LRC(14) aligned sectors.

If k tails are aligned to the body ruler and p=7-k are drifts, the body
projection Y_D(A) has mass

    (|S_D|/D) * mu(R_A).

The universal safe floors u_j give mu(R_A)>=u_k, while the union of the p
drift danger combs has mass at most 1-u_p.  Therefore every cover satisfies

    |S_D|/D <= (1-u_p)/u_k.

This referee counts the exact body/divisor rows surviving that necessary
condition for every k.  It also counts denominator-multiset shapes without
enumerating them: if a_p(D) is the number of nondecreasing p-tuples of
divisors >1 with lcm D, divisor-lattice Möbius inversion gives

    a_p(D) = sum_{e|D} mu(D/e) * C(tau(e)+p-2, p).
"""

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import comb, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMBINED_PATH = (
    ROOT / "04-computation" / "lrc14_three_drift_body_projection_fiber_thm2928.py"
)
EXPECTED_COMBINED_SHA256 = (
    "42dc165781148c702dfcd3c6535f4d02aee516af60b5ddf602a19cb1d87695e4"
)
COMBINED_OUTPUT_PATH = (
    ROOT
    / "05-knowledge"
    / "results"
    / "lrc14_three_drift_body_projection_fiber_thm2928.out"
)
EXPECTED_COMBINED_OUTPUT_SHA256 = (
    "2e211620ad7064ea06f7544b5fbac709d6d52d9a0e261b464ae26b595f09b669"
)

SAFE_FLOORS = {
    0: Q(1),
    1: Q(6, 7),
    2: Q(66, 91),
    3: Q(55, 91),
    4: Q(558, 1183),
    5: Q(478, 1365),
    6: Q(61, 273),
    7: Q(15, 154),
}

EXPECTED_CUTOFFS = {
    0: Q(139, 154),
    1: Q(106, 117),
    2: Q(887, 990),
    3: Q(125, 143),
    4: Q(26, 31),
    5: Q(375, 478),
    6: Q(39, 61),
    7: Q(0),
}

EXPECTED_ROWS = (27210, 27240, 27163, 26970, 13778, 10976, 6237, 0)
EXPECTED_BODIES = (3003, 3003, 3003, 3003, 3003, 3003, 3003, 0)
EXPECTED_DIVISORS = (219, 219, 219, 217, 206, 193, 171, 0)
EXPECTED_SHAPES = (
    161535777082757,
    3095010121875,
    50874159718,
    694921995,
    7483350,
    56419,
    171,
    0,
)
EXPECTED_OCCURRENCES = (
    1504842061942849,
    38954725590760,
    951545890235,
    21357714101,
    298255882,
    3066274,
    6237,
    0,
)
EXPECTED_SEMANTIC_SHA256 = (
    "e7ecc3d7aed7b692564a03734f0f240d4f9c6ca9d68d012a66f5a4a626d0a7c4",
    "931e04a95f70e4b16a3004496df981a0678f0839dbfa85bc6194d89290a4a76b",
    "2c2909d3d1e514b4e8f09166d23145ac6b577fe6e10f39d39a205094e7656dd5",
    "3738fce411911b8cfb64596a6ed595a01d54f9093b502b4a4161fcd63ac91b2a",
    "86109feddd36a3ece9068d4e50d091b8423930a226f3c50c5024578df41572b0",
    "38aa87113f5acf83273d756e1ad9c261f745d1f3ae106b6d1d1aa5d02ef03739",
    "ee5bd68fe275192c6eeb6b09081982c226251b0305268fdbb91df894da250941",
    "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def mobius(number):
    result = 1
    prime = 2
    remaining = number
    while prime * prime <= remaining:
        if remaining % prime:
            prime += 1
            continue
        remaining //= prime
        if remaining % prime == 0:
            return 0
        result = -result
        while remaining % prime == 0:
            remaining //= prime
        prime += 1
    if remaining > 1:
        result = -result
    return result


require(
    file_sha256(COMBINED_PATH) == EXPECTED_COMBINED_SHA256,
    "combined three-drift dependency changed",
)
require(
    file_sha256(COMBINED_OUTPUT_PATH) == EXPECTED_COMBINED_OUTPUT_SHA256,
    "combined three-drift output changed",
)
spec = spec_from_file_location("lrc14_three_drift_combined", COMBINED_PATH)
combined = module_from_spec(spec)
spec.loader.exec_module(combined)
support = combined.support_module


def lcm_multiset_shapes(D, arity):
    """Nondecreasing arity-tuples of divisors >1 with lcm exactly D."""
    total = 0
    for divisor in support.divisors(D):
        divisor_count = len(support.divisors(divisor)) - 1
        all_multisets = comb(divisor_count + arity - 1, arity)
        total += mobius(D // divisor) * all_multisets
    require(total >= 0, ("negative lcm-multiset count", D, arity, total))
    return total


def brute_shape_controls():
    cases = 0
    for D in range(1, 61):
        allowed = tuple(d for d in support.divisors(D) if d > 1)
        for arity in range(1, 5):
            brute = sum(
                lcm(*values) == D
                for values in combinations_with_replacement(allowed, arity)
            )
            exact = lcm_multiset_shapes(D, arity)
            require(
                exact == brute,
                ("lcm-multiset inversion failed", D, arity, exact, brute),
            )
            cases += 1
    return cases


def main():
    shape_control_cases = brute_shape_controls()
    cutoffs = {
        k: (Q(1) - SAFE_FLOORS[7 - k]) / SAFE_FLOORS[k]
        for k in range(8)
    }
    require(cutoffs == EXPECTED_CUTOFFS, ("support cutoffs changed", cutoffs))

    rows_by_k_D = {k: Counter() for k in range(8)}
    bodies_by_k = {k: set() for k in range(8)}
    semantic = {k: sha256() for k in range(8)}
    total_rows = 0
    body_count = 0
    for F in combinations(range(1, 15), 6):
        body_count += 1
        L, ranges = support.safe_cell_ranges(F)
        for D in support.divisors(L):
            total_rows += 1
            support_count = support.support_size_bitset(D, ranges)
            require(0 < support_count <= D, ("support size out of range", F, D))
            for k, cutoff in cutoffs.items():
                if Q(support_count, D) <= cutoff:
                    rows_by_k_D[k][D] += 1
                    bodies_by_k[k].add(F)
                    semantic[k].update(
                        f"{F}|{L}|{D}|{support_count}\n".encode()
                    )

    require(body_count == 3003, "body universe changed")
    require(total_rows == 251536, "body/divisor universe changed")
    require(sum(rows_by_k_D[4].values()) == 13778, "k=4 support count changed")
    require(len(rows_by_k_D[4]) == 206, "k=4 divisor alphabet changed")
    require(sum(rows_by_k_D[5].values()) == 10976, "k=5 support count changed")

    print("LRC14 aligned/drift support-transfer ladder referee")
    print(f"combined_script_sha256={file_sha256(COMBINED_PATH)}")
    print(f"combined_output_sha256={file_sha256(COMBINED_OUTPUT_PATH)}")
    print(f"shape_control_cases={shape_control_cases}")
    print(f"body_count={body_count}")
    print(f"body_divisor_rows={total_rows}")
    print(f"safe_floors={SAFE_FLOORS}")
    print(f"support_cutoffs={cutoffs}")
    actual_rows = []
    actual_bodies = []
    actual_divisors = []
    actual_shapes = []
    actual_occurrences = []
    actual_semantic = []
    for k in range(8):
        p = 7 - k
        row_count = sum(rows_by_k_D[k].values())
        D_count = len(rows_by_k_D[k])
        if p:
            shape_count = sum(
                lcm_multiset_shapes(D, p) for D in rows_by_k_D[k]
            )
            occurrence_count = sum(
                row_multiplicity * lcm_multiset_shapes(D, p)
                for D, row_multiplicity in rows_by_k_D[k].items()
            )
        else:
            shape_count = 0
            occurrence_count = 0
        actual_rows.append(row_count)
        actual_bodies.append(len(bodies_by_k[k]))
        actual_divisors.append(D_count)
        actual_shapes.append(shape_count)
        actual_occurrences.append(occurrence_count)
        actual_semantic.append(semantic[k].hexdigest())
        print(
            f"k={k},p={p},cutoff={cutoffs[k]},rows={row_count},"
            f"bodies={len(bodies_by_k[k])},divisors={D_count},"
            f"denominator_shapes={shape_count},"
            f"raw_occurrences={occurrence_count},"
            f"semantic_sha256={semantic[k].hexdigest()}"
        )
    require(tuple(actual_rows) == EXPECTED_ROWS, "support row ladder changed")
    require(tuple(actual_bodies) == EXPECTED_BODIES, "support body ladder changed")
    require(
        tuple(actual_divisors) == EXPECTED_DIVISORS,
        "support divisor ladder changed",
    )
    require(
        tuple(actual_shapes) == EXPECTED_SHAPES,
        "denominator shape ladder changed",
    )
    require(
        tuple(actual_occurrences) == EXPECTED_OCCURRENCES,
        "raw occurrence ladder changed",
    )
    require(
        tuple(actual_semantic) == EXPECTED_SEMANTIC_SHA256,
        "support semantic ladder changed",
    )
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
