#!/usr/bin/env python3
"""Compose THM-3366's pool-14 terminals with the refined k=3 ledger."""

from bisect import bisect_left
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
K3_PATH = ROOT / "04-computation/lrc14_k3_four_drift_expected_spike_gf_thm2928.py"
COMPLEMENT_BASE_PATH = (
    ROOT / "04-computation/lrc14_k1_body_complement_clock_scan_kps_s172.py"
)
EXPECTED_K3_SHA256 = "05e365a654b32e66b814dcbce9385a2d13c22a2c84a5474e0855dcab6262b055"
EXPECTED_COMPLEMENT_BASE_LF_SHA256 = (
    "bdb2001cf22f7e92884e895b0095021e42e8f1febd9adbf779b250a2f6c53507"
)
POOL = tuple(range(1, 15))
EXPECTED_COMPLEMENT_ROWS = 19_053
EXPECTED_COMPLEMENT_HISTOGRAM = ((1, 12_659), (2, 2_764), (3, 976), (4, 1_052), (5, 1_602))
EXPECTED_PRE = (398_241_574, 2_548_901_482, 1_904, 1_823, 107)
EXPECTED_PRE_BY_C = {
    1: 71_619_386,
    2: 1_351_841_956,
    3: 1_065_317_472,
    4: 60_122_668,
}
EXPECTED_PRE_SEMANTIC = "06b3e5a7a05d5c9a2f1633d74061f1605c2aab856b315026ad697245bac84964"
EXPECTED_REFINED_DELETION = (7, 7_648)
EXPECTED_POST = (398_241_574, 2_548_893_834, 1_897, 1_823, 107)
EXPECTED_POST_BY_C = {
    1: 71_619_386,
    2: 1_351_841_956,
    3: 1_065_309_824,
    4: 60_122_668,
}
EXPECTED_POST_SEMANTIC = "1dd557b3629ed4c48f2b85412b48476b29e300990f38bdfe65cef2372fbab8d4"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def file_hash(path):
    return sha256(path.read_bytes()).hexdigest()


def normalized_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("module import", path))
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
        if any(
            max(left, gap_left) < min(right, gap_right)
            for gap_left, gap_right in gaps
        ):
            target |= 1 << (offset + index)
    return target


def build_solver(base):
    points = base.arrangement_points(POOL)
    atoms = tuple(zip(points, points[1:]))
    samples = points + tuple((left + right) / 2 for left, right in atoms)
    masks = tuple(
        sum((1 << i for i, x in enumerate(samples) if base.danger(c, x)), 0)
        for c in POOL
    )
    union = 0
    for mask in masks:
        union |= mask
    by_bit = tuple(
        tuple(i for i, mask in enumerate(masks) if mask & (1 << bit))
        for bit in range(len(samples))
    )

    @lru_cache(maxsize=None)
    def exact(remaining, depth):
        if remaining == 0:
            return ()
        if depth == 0:
            return None
        bits = tuple(i for i in range(len(samples)) if remaining & (1 << i))
        pivot = min(bits, key=lambda i: len(by_bit[i]))
        for candidate in by_bit[pivot]:
            reduced = remaining & ~masks[candidate]
            if reduced == remaining:
                continue
            suffix = exact(reduced, depth - 1)
            if suffix is not None:
                return (POOL[candidate],) + suffix
        return None

    def solve(target):
        if target & ~union:
            return None
        for depth in range(6):
            result = exact(target, depth)
            if result is not None:
                return tuple(sorted(result))
        return None

    return points, atoms, solve, exact


def aggregate(k3, by_divisor, allowed_keys=None):
    occurrences = 0
    shapes = 0
    rows = set()
    bodies = set()
    divisors = set()
    row_counts = {}
    occurrence_by_c = Counter()
    shape_by_c = Counter()
    semantic = sha256()

    for D in sorted(by_divisor):
        q = D // 7
        distribution = k3.expected_distribution(D)
        c3_distribution = k3.c3_denominator_distribution(D)
        by_pattern_c = defaultdict(Counter)
        for (pattern, c, capacity), multiplicity in distribution.items():
            by_pattern_c[(pattern, c)][capacity] += multiplicity
        suffix_data = {key: k3.suffix(counter) for key, counter in by_pattern_c.items()}
        used_non_c3 = set()
        used_c3 = set()

        for support_count, body, L, arcs in by_divisor[D]:
            key = (body, D)
            if allowed_keys is not None and key not in allowed_keys:
                continue
            histogram_q = k3.combined.residue_load_histogram(arcs, q)
            N_by_c = {
                c: sum(count for load, count in histogram_q if load > c)
                for c in range(5)
            }
            small_loads = {}
            for divisor in (2, 3, 4):
                if D % divisor == 0:
                    histogram_d = k3.combined.residue_load_histogram(arcs, divisor)
                    small_loads[divisor] = k3.combined.top_class_load(histogram_d, 1)

            support_limits = {}
            for pattern, _c in by_pattern_c:
                if pattern in support_limits:
                    continue
                weights, marginals = k3.base.small_vectors(pattern, D, small_loads)
                support_limits[pattern] = k3.base.status_need_limit(weights, marginals)

            row_total = 0
            row_by_c = Counter()
            for (pattern, c), data in suffix_data.items():
                if c == 3:
                    continue
                threshold = support_count - support_limits[pattern]
                count = k3.at_least(data, threshold)
                if not count or not k3.mean_passes(N_by_c[c], c, q):
                    continue
                row_total += count
                row_by_c[c] += count
                for capacity in data[0][bisect_left(data[0], threshold):]:
                    used_non_c3.add((pattern, c, capacity))

            if k3.mean_passes(N_by_c[3], 3, q):
                for feature, multiplicity in c3_distribution.items():
                    pattern, _spike_d, capacity, allowance = feature
                    if capacity + support_limits[pattern] < support_count:
                        continue
                    if allowance < N_by_c[3]:
                        continue
                    row_total += multiplicity
                    row_by_c[3] += multiplicity
                    used_c3.add(feature)

            if not row_total:
                continue
            rows.add(key)
            bodies.add(body)
            divisors.add(D)
            row_counts[key] = row_total
            occurrences += row_total
            occurrence_by_c.update(row_by_c)
            semantic.update(f"{body}|{L}|{D}|{support_count}|{row_total}\n".encode())

        for feature in used_non_c3:
            _pattern, c, _capacity = feature
            multiplicity = distribution[feature]
            shapes += multiplicity
            shape_by_c[c] += multiplicity
        for feature in used_c3:
            multiplicity = c3_distribution[feature]
            shapes += multiplicity
            shape_by_c[3] += multiplicity

    return {
        "summary": (shapes, occurrences, len(rows), len(bodies), len(divisors)),
        "rows": rows,
        "row_counts": row_counts,
        "occurrence_by_c": occurrence_by_c,
        "shape_by_c": shape_by_c,
        "semantic": semantic.hexdigest(),
    }


def main():
    require(file_hash(K3_PATH) == EXPECTED_K3_SHA256, "k3 dependency changed")
    require(
        normalized_hash(COMPLEMENT_BASE_PATH) == EXPECTED_COMPLEMENT_BASE_LF_SHA256,
        "complement base changed",
    )
    k3 = load_module("kps_s176_k3", K3_PATH)
    complement = load_module("kps_s176_complement", COMPLEMENT_BASE_PATH)
    by_divisor, body_count, body_divisor_rows = k3.q7.build_rows()
    require((body_count, body_divisor_rows) == (3003, 251536), "body universe")

    points, atoms, solve, exact = build_solver(complement)
    closed = set()
    completion_histogram = Counter()
    completion_semantic = sha256()
    for D in sorted(by_divisor):
        for support_count, body, L, arcs in by_divisor[D]:
            gaps = complement.unsupported_gaps(D, arcs)
            target = target_mask(complement, D, gaps, points, atoms)
            completion = solve(target)
            if completion is None:
                continue
            closed.add((body, D))
            completion_histogram[len(completion)] += 1
            completion_semantic.update(
                f"{body}|{L}|{D}|{support_count}|{completion}\n".encode()
            )

    require(len(closed) == EXPECTED_COMPLEMENT_ROWS, len(closed))
    require(
        tuple(sorted(completion_histogram.items())) == EXPECTED_COMPLEMENT_HISTOGRAM,
        ("completion histogram", tuple(sorted(completion_histogram.items()))),
    )

    before = aggregate(k3, by_divisor)
    require(before["summary"] == EXPECTED_PRE, ("pre summary", before["summary"]))
    require(dict(before["occurrence_by_c"]) == EXPECTED_PRE_BY_C, "pre by-c")
    require(before["semantic"] == EXPECTED_PRE_SEMANTIC, "pre semantic")
    allowed = before["rows"] - closed
    after = aggregate(k3, by_divisor, allowed)
    killed_rows = before["rows"] & closed
    killed_occurrences = sum(before["row_counts"][key] for key in killed_rows)
    require(before["summary"][1] - after["summary"][1] == killed_occurrences, "subtraction")
    require(
        (len(killed_rows), killed_occurrences) == EXPECTED_REFINED_DELETION,
        ("refined deletion", len(killed_rows), killed_occurrences),
    )
    require(after["summary"] == EXPECTED_POST, ("post summary", after["summary"]))
    require(dict(after["occurrence_by_c"]) == EXPECTED_POST_BY_C, "post by-c")
    require(after["semantic"] == EXPECTED_POST_SEMANTIC, "post semantic")
    killed_detail = tuple(
        sorted((D, body, before["row_counts"][(body, D)]) for body, D in killed_rows)
    )

    print("LRC14 K=3 REFINED COMPLEMENT-CLOCK COMPOSITION")
    print("status=FINITE-EXACT intersection of final one-spike and pool-14 ledgers")
    print(
        f"arrangement_points={len(points)};atoms={len(atoms)};"
        f"support_closed_rows={len(closed)};"
        f"support_completion_histogram={tuple(sorted(completion_histogram.items()))}"
    )
    print(f"pre_summary={before['summary']}")
    print(f"pre_occurrence_by_c={before['occurrence_by_c']}")
    print(f"pre_shape_by_c={before['shape_by_c']}")
    print(
        f"refined_closed_rows={len(killed_rows)};"
        f"refined_closed_occurrences={killed_occurrences}"
    )
    print(f"refined_closed_detail={killed_detail}")
    print(f"post_summary={after['summary']}")
    print(f"post_occurrence_by_c={after['occurrence_by_c']}")
    print(f"post_shape_by_c={after['shape_by_c']}")
    print(f"completion_semantic_sha256={completion_semantic.hexdigest()}")
    print(f"pre_semantic_sha256={before['semantic']}")
    print(f"post_semantic_sha256={after['semantic']}")
    print(f"solver_cache={exact.cache_info()}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
