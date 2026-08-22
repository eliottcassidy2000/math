#!/usr/bin/env python3
"""Exact set-cover scout for arbitrary small complement clocks in k=1.

This strengthens the descended-body-clock scan S172.  On the common exact
arrangement of D_1,...,D_14, each unsupported projected cell contributes open
atoms and equality endpoints.  A bitset set-cover solver asks whether at most
five of these fourteen integer combs cover that pointwise target.  Any positive
answer closes the support row by the small-LRC complement-clock lemma.
"""

import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation/lrc14_k1_body_complement_clock_scan_kps_s172.py"
EXPECTED_BASE_SHA256 = "bdb2001cf22f7e92884e895b0095021e42e8f1febd9adbf779b250a2f6c53507"
DEFAULT_POOL_MAX = 14
DEFAULT_MAX_CLOCKS = 5
EXPECTED_DEFAULT = {
    "distinct_targets": 11_180,
    "closed_rows": 19_273,
    "size_histogram": ((1, 12_665), (2, 2_766), (3, 1_007), (4, 1_142), (5, 1_693)),
    "closed_occurrences": 5_116_993_586_195,
    "occurrence_histogram": (
        (1, 3_427_334_397_071),
        (2, 332_132_448_500),
        (3, 179_613_246_178),
        (4, 497_575_939_015),
        (5, 680_337_555_431),
    ),
    "semantic": "dfa63a2d6800c2ac67fa467bb5d576fcfcd7f141c01739846c3c6299018504c8",
}


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def normalized_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_base():
    spec = spec_from_file_location("kps_s173_base", BASE_PATH)
    require(spec is not None and spec.loader is not None, "base import")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def target_mask(base, D, gaps, points, atoms):
    target = 0
    for index, x in enumerate(points):
        if (D * x).denominator != 1 and any(left < x < right for left, right in gaps):
            target |= 1 << index
    offset = len(points)
    for index, (left, right) in enumerate(atoms):
        if any(max(left, gap_left) < min(right, gap_right) for gap_left, gap_right in gaps):
            target |= 1 << (offset + index)
    return target


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pool-max", type=int, default=DEFAULT_POOL_MAX)
    parser.add_argument("--max-clocks", type=int, default=DEFAULT_MAX_CLOCKS)
    args = parser.parse_args()
    require(args.pool_max >= 1, args.pool_max)
    require(0 <= args.max_clocks <= 5, args.max_clocks)
    pool = tuple(range(1, args.pool_max + 1))
    require(normalized_hash(BASE_PATH) == EXPECTED_BASE_SHA256, "base dependency changed")
    base = load_base()
    support = base.load_support()
    points = base.arrangement_points(pool)
    atoms = tuple(zip(points, points[1:]))
    sample_points = points + tuple((left + right) / 2 for left, right in atoms)
    masks = tuple(
        sum((1 << index for index, x in enumerate(sample_points) if base.danger(clock, x)), 0)
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
        bits = tuple(index for index in range(len(sample_points)) if remaining & (1 << index))
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
        for depth in range(args.max_clocks + 1):
            result = solve_exact(target, depth)
            if result is not None:
                return tuple(sorted(set(result)))
        return None

    rows = 0
    closed = []
    size_histogram = Counter()
    occurrence_histogram = Counter()
    target_histogram = Counter()
    semantic = sha256()
    for body in combinations(range(1, 15), 6):
        L, ranges = support.safe_cell_ranges(body)
        for D in support.divisors(L):
            count = support.support_size_bitset(D, ranges)
            if Q(count, D) > base.CUTOFF:
                continue
            rows += 1
            arcs = base.residue_arcs(D, ranges)
            gaps = base.unsupported_gaps(D, arcs)
            target = target_mask(base, D, gaps, points, atoms)
            target_histogram[target] += 1
            completion = solve(target)
            if completion is not None:
                row = (body, L, D, count, completion)
                closed.append(row)
                size_histogram[len(completion)] += 1
                occurrence_histogram[len(completion)] += base.denominator_shapes(D)
                semantic.update((repr(row) + "\n").encode())

    require(rows == 27_240, rows)
    closed_occurrences = sum(base.denominator_shapes(row[2]) for row in closed)
    if (
        args.pool_max == DEFAULT_POOL_MAX
        and args.max_clocks == DEFAULT_MAX_CLOCKS
    ):
        require(len(target_histogram) == EXPECTED_DEFAULT["distinct_targets"], "targets")
        require(len(closed) == EXPECTED_DEFAULT["closed_rows"], "closed rows")
        require(
            tuple(sorted(size_histogram.items())) == EXPECTED_DEFAULT["size_histogram"],
            "size histogram",
        )
        require(
            closed_occurrences == EXPECTED_DEFAULT["closed_occurrences"],
            "closed occurrences",
        )
        require(
            tuple(sorted(occurrence_histogram.items()))
            == EXPECTED_DEFAULT["occurrence_histogram"],
            "occurrence histogram",
        )
        require(semantic.hexdigest() == EXPECTED_DEFAULT["semantic"], "semantic")
    print("LRC14 K=1 UNIVERSAL COMPLEMENT-CLOCK SCOUT")
    print("status=FINITE-EXACT pointwise sufficient terminals over declared finite clock pool")
    print(f"pool=1..{args.pool_max};max_clocks={args.max_clocks};arrangement_points={len(points)};atoms={len(atoms)}")
    print(f"input_rows={rows};distinct_targets={len(target_histogram)};closed_rows={len(closed)};remaining_rows={rows-len(closed)}")
    print(f"completion_size_histogram={tuple(sorted(size_histogram.items()))}")
    print(f"closed_occurrences={closed_occurrences};remaining_occurrences={base.INPUT_OCCURRENCES-closed_occurrences}")
    print(f"occurrence_size_histogram={tuple(sorted(occurrence_histogram.items()))}")
    print(f"first_rows={tuple(closed[:20])}")
    print(f"last_rows={tuple(closed[-5:])}")
    print(f"base_sha256={normalized_hash(BASE_PATH)}")
    print(f"semantic_sha256={semantic.hexdigest()}")
    print(f"solver_cache={solve_exact.cache_info()}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
