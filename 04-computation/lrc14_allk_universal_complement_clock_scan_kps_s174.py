#!/usr/bin/env python3
"""Exact complement-clock terminals across every aligned/drift sector.

For a body/divisor row (F,L,D,S_D), the projected carrier after k aligned
clocks A is

    Y_D(A) = union_{r in S_D} (r + R_A)/D.

If k>=1, its complement is covered by the k integer clocks D*a (a in A)
and by any clocks that cover the unsupported open D-cells.  A hypothetical
cover by p=7-k drift clocks would then cover the whole circle with at most
7+r integer clocks.  Thus r<=5 contradicts the cited LRC theorem through
twelve nonzero speeds.  For k=0 one additionally uses clock D to own the
D-grid boundary, so r<=4 gives at most twelve clocks.

The unsupported-cell target is independent of k.  This program constructs
the exact common arrangement of D_1,...,D_pool_max, checks both equality
endpoints and open atoms, solves the resulting finite set-cover problem, and
weights each closed support row by the exact number of denominator multisets.
It is a support-level terminal only; it does not assert that every surviving
support row is physically realizable.
"""

import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement
from math import comb, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation/lrc14_k1_body_complement_clock_scan_kps_s172.py"
EXPECTED_BASE_SHA256 = "bdb2001cf22f7e92884e895b0095021e42e8f1febd9adbf779b250a2f6c53507"

DEFAULT_POOL_MAX = 14
MAX_COMPLETION = 5

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
CUTOFFS = {
    k: (Q(1) - SAFE_FLOORS[7 - k]) / SAFE_FLOORS[k]
    for k in range(8)
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
EXPECTED_ROWS = (27_210, 27_240, 27_163, 26_970, 13_778, 10_976, 6_237, 0)
EXPECTED_OCCURRENCES = (
    1_504_842_061_942_849,
    38_954_725_590_760,
    951_545_890_235,
    21_357_714_101,
    298_255_882,
    3_066_274,
    6_237,
    0,
)
EXPECTED_DEFAULT_CLOSED_ROWS = (17_561, 19_273, 19_198, 19_053, 5_964, 3_334, 570)
EXPECTED_DEFAULT_CLOSED_OCCURRENCES = (
    113_853_177_071_279,
    5_116_993_586_195,
    182_987_150_531,
    5_887_257_171,
    47_788_532,
    528_176,
    570,
)
EXPECTED_DEFAULT_SEMANTICS = (
    "b062ded580bae9def1ca5c2eef9c3b6463ed00c96a3677e50a3336f85e9db5e5",
    "922a071ee28536ca4be4ea6730058bac52389e5d47f0140acaa9b69d1de5ca86",
    "a0d2f6b417a65095ddc78136f4bc6384905b70371724bb1ef869f92fdfc3fafd",
    "ef908419fae1edbb8326f87e781afc75733520bed63be5a52feb6c16fc55b1d8",
    "b2ebf3cb4372b8d71c1099c6489f16eedddc436ff40d3f88504cdaba8c139b85",
    "7f5b1eaf4139c05c856adedb55d24bb7b52c8726d3b51467cde10414735ac776",
    "b3bc33b5ec23c74c181b05e6665a076b52d60b0b77f3b7cf032803b4f1f38690",
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def normalized_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_base():
    require(normalized_hash(BASE_PATH) == EXPECTED_BASE_SHA256, "base dependency changed")
    spec = spec_from_file_location("kps_s174_base", BASE_PATH)
    require(spec is not None and spec.loader is not None, "base import")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=None)
def denominator_shapes(base, D, arity):
    """Nondecreasing arity-tuples of divisors >1 having lcm exactly D."""
    return sum(
        base.mobius(D // divisor)
        * comb(len(base.all_divisors(divisor)) + arity - 2, arity)
        for divisor in base.all_divisors(D)
    )


def target_mask(base, D, gaps, points, atoms):
    """Bits for the unsupported open cells, with D-grid points removed."""
    target = 0
    for index, x in enumerate(points):
        if (D * x).denominator != 1 and any(left < x < right for left, right in gaps):
            target |= 1 << index
    offset = len(points)
    for index, (left, right) in enumerate(atoms):
        if any(
            max(left, gap_left) < min(right, gap_right)
            for gap_left, gap_right in gaps
        ):
            target |= 1 << (offset + index)
    return target


def brute_shape_controls(base):
    cases = 0
    for D in range(1, 41):
        allowed = tuple(d for d in base.all_divisors(D) if d > 1)
        for arity in range(1, 8):
            brute = sum(
                lcm(*values) == D
                for values in combinations_with_replacement(allowed, arity)
            )
            exact = denominator_shapes(base, D, arity)
            require(exact == brute, ("shape inversion", D, arity, exact, brute))
            cases += 1
    return cases


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pool-max", type=int, default=DEFAULT_POOL_MAX)
    args = parser.parse_args()
    require(args.pool_max >= 1, args.pool_max)
    require(CUTOFFS == EXPECTED_CUTOFFS, CUTOFFS)

    base = load_base()
    support = base.load_support()
    shape_control_cases = brute_shape_controls(base)

    pool = tuple(range(1, args.pool_max + 1))
    points = base.arrangement_points(pool)
    atoms = tuple(zip(points, points[1:]))
    sample_points = points + tuple((left + right) / 2 for left, right in atoms)
    masks = tuple(
        sum(
            (1 << index for index, x in enumerate(sample_points) if base.danger(clock, x)),
            0,
        )
        for clock in pool
    )
    union_mask = 0
    for mask in masks:
        union_mask |= mask
    covers_by_bit = tuple(
        tuple(index for index, mask in enumerate(masks) if mask & (1 << bit))
        for bit in range(len(sample_points))
    )

    @lru_cache(maxsize=None)
    def solve_exact(remaining, depth):
        if remaining == 0:
            return ()
        if depth == 0:
            return None
        bits = tuple(
            index
            for index in range(len(sample_points))
            if remaining & (1 << index)
        )
        bit = min(bits, key=lambda index: len(covers_by_bit[index]))
        for candidate in covers_by_bit[bit]:
            reduced = remaining & ~masks[candidate]
            if reduced == remaining:
                continue
            suffix = solve_exact(reduced, depth - 1)
            if suffix is not None:
                return (pool[candidate],) + suffix
        return None

    def solve(target):
        if target & ~union_mask:
            return None
        for depth in range(MAX_COMPLETION + 1):
            result = solve_exact(target, depth)
            if result is not None:
                require(len(result) == len(set(result)), ("duplicate clocks", result))
                return tuple(sorted(result))
        return None

    rows = [0] * 8
    input_occurrences = [0] * 8
    closed_rows = [0] * 8
    closed_occurrences = [0] * 8
    size_histograms = [Counter() for _ in range(8)]
    occurrence_histograms = [Counter() for _ in range(8)]
    semantics = [sha256() for _ in range(8)]
    distinct_targets = Counter()

    for body in combinations(range(1, 15), 6):
        L, ranges = support.safe_cell_ranges(body)
        for D in support.divisors(L):
            count = support.support_size_bitset(D, ranges)
            eligible = tuple(k for k in range(7) if Q(count, D) <= CUTOFFS[k])
            if not eligible:
                continue
            arcs = base.residue_arcs(D, ranges)
            require(
                sum(right - left for left, right in arcs) == count,
                ("support arcs", body, D, count, arcs),
            )
            gaps = base.unsupported_gaps(D, arcs)
            target = target_mask(base, D, gaps, points, atoms)
            distinct_targets[target] += 1
            completion = solve(target)

            for k in eligible:
                p = 7 - k
                weight = denominator_shapes(base, D, p)
                rows[k] += 1
                input_occurrences[k] += weight
                budget = 4 if k == 0 else 5
                if completion is None or len(completion) > budget:
                    continue
                closed_rows[k] += 1
                closed_occurrences[k] += weight
                size_histograms[k][len(completion)] += 1
                occurrence_histograms[k][len(completion)] += weight
                semantics[k].update(
                    f"{body}|{L}|{D}|{count}|{completion}\n".encode()
                )

    require(tuple(rows) == EXPECTED_ROWS, ("row ladder", tuple(rows)))
    require(
        tuple(input_occurrences) == EXPECTED_OCCURRENCES,
        ("occurrence ladder", tuple(input_occurrences)),
    )
    require(denominator_shapes(base, 14, 6) == 26, "D14 six-shape control")
    if args.pool_max == DEFAULT_POOL_MAX:
        require(
            tuple(closed_rows[:7]) == EXPECTED_DEFAULT_CLOSED_ROWS,
            ("default closed-row ladder", tuple(closed_rows[:7])),
        )
        require(
            tuple(closed_occurrences[:7]) == EXPECTED_DEFAULT_CLOSED_OCCURRENCES,
            ("default closed-occurrence ladder", tuple(closed_occurrences[:7])),
        )
        require(
            tuple(digest.hexdigest() for digest in semantics[:7])
            == EXPECTED_DEFAULT_SEMANTICS,
            "default semantic ladder changed",
        )

    print("LRC14 ALL-K UNIVERSAL COMPLEMENT-CLOCK TERMINALS")
    print("status=FINITE-EXACT support-level pointwise sufficient terminals")
    print(
        f"pool=1..{args.pool_max};arrangement_points={len(points)};"
        f"atoms={len(atoms)};shape_control_cases={shape_control_cases};"
        f"distinct_targets={len(distinct_targets)}"
    )
    for k in range(7):
        p = 7 - k
        budget = 4 if k == 0 else 5
        print(
            f"k={k};p={p};completion_budget={budget};cutoff={CUTOFFS[k]};"
            f"input_rows={rows[k]};closed_rows={closed_rows[k]};"
            f"remaining_rows={rows[k]-closed_rows[k]};"
            f"input_occurrences={input_occurrences[k]};"
            f"closed_occurrences={closed_occurrences[k]};"
            f"remaining_occurrences={input_occurrences[k]-closed_occurrences[k]};"
            f"size_histogram={tuple(sorted(size_histograms[k].items()))};"
            f"occurrence_size_histogram={tuple(sorted(occurrence_histograms[k].items()))};"
            f"semantic_sha256={semantics[k].hexdigest()}"
        )
    print(f"base_sha256={normalized_hash(BASE_PATH)}")
    print(f"solver_cache={solve_exact.cache_info()}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
