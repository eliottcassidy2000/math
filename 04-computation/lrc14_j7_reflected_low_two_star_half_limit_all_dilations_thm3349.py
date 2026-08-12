#!/usr/bin/env python3
"""Exact finite head for THM-3349's selected half-limit theorem.

The theorem supplies the analytic midpoint/shear bound.  This companion
reconstructs the frozen lexicographically-largest-edge selector on the
649-body, 11,856-ray atlas and checks exactly the finitely many physical
scales preceding the analytic tail.
"""

from __future__ import annotations

import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SELECTOR = ROOT / "04-computation/lrc14_j7_reflected_low_two_star_limit_selector_scout_20260812.py"
EXPECTED_SELECTOR_SHA256 = "cd8b08087f0f7e1e0c0c7d0be673629c7c2702c170c5c1e771e1d76df1d3cd1c"
EXPECTED_THRESHOLD_HIST = ((2, 140669), (3, 29), (4, 7), (5, 4), (6, 1), (9, 1))
EXPECTED_HEAD_DIGEST = "fd0b48f6975a38ca271106c75f97b9ea9a43d6348ede29289a2438b860d544b0"
EXPECTED_MAX_THRESHOLD = (
    9, (1, 2, 3, 4, 6, 12), 168, 90, (1, 5), 6, 1, F(29, 1470),
    2, 12, 1, 6, 1, -70,
)
EXPECTED_WEAKEST_HEAD = (
    F(553321691, 3504897823680),
    (2, 4, 6, 9, 11, 12), 5544, 2970, (0, 5), 12, 9, 2,
    F(174636, 553172005), F(1, 3168), 3,
)


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def load(path: Path):
    spec = spec_from_file_location("half_limit_selector", path)
    require(spec is not None and spec.loader is not None, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gamma(level: int) -> F:
    return F(1, 2) if level % 2 == 0 else F(level * level + 1, 2 * level * level)


def analytic_error_bound(scale, ruler, divisor, e, P, f, Q, cross):
    h = scale * divisor
    return (
        gamma(P) * F(e * P, h * ruler * P - e)
        + gamma(Q) * F(f * Q, h * ruler * Q - f)
        + F(abs(cross), 2 * h * h * ruler * P * Q)
    )


def tail_threshold(ruler, e, f, p, q, limit):
    divisor = gcd(p, q)
    P, Q = p // divisor, q // divisor
    cross = Q * e - P * f
    require(P + Q <= 7, ("not low", P, Q))
    require(abs(cross) < ruler, ("multiple transverse turns", cross, ruler))
    for scale in range(2, 100):
        if analytic_error_bound(scale, ruler, divisor, e, P, f, Q, cross) <= limit / 2:
            return scale, divisor, P, Q, cross
    raise RuntimeError(("tail threshold exceeds search", ruler, e, f, p, q, limit))


def main() -> None:
    selector_hash = hashlib.sha256(
        SELECTOR.read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    require(
        selector_hash == EXPECTED_SELECTOR_SHA256,
        ("selector hash", selector_hash, EXPECTED_SELECTOR_SHA256),
    )
    m = load(SELECTOR)
    bodies = m.MS.body_universe()
    rays, _ = m.MS.projective_low_two_star_rays()
    require((len(bodies), len(rays), len(bodies) * len(rays)) == (649, 11856, 7694544), "universe")

    threshold_hist = Counter()
    group_count = exception_groups = head_checks = 0
    weakest_head = max_threshold = None
    semantic = hashlib.sha256()

    for body, ruler in bodies:
        cell, _, _, positions, edges = m.body_geometry(body, ruler)
        words = []
        pair_keys = {edge: set() for edge in edges}
        for ray in rays:
            word_list = [0] * 6
            for position, value in zip(positions, ray):
                word_list[position] = value
            word = tuple(word_list)
            words.append(word)
            for edge in edges:
                pair_keys[edge].add((word[edge[0]], word[edge[1]]))
        limits = {
            (edge, p, q): m.homogeneous_limit(
                ruler, cell, body[edge[0]], p, body[edge[1]], q
            )
            for edge, pairs in pair_keys.items()
            for p, q in pairs
        }
        groups = set()
        for word in words:
            limit, edge = max(
                (limits[edge, word[edge[0]], word[edge[1]]], edge)
                for edge in edges
            )
            groups.add((edge, word[edge[0]], word[edge[1]], limit))

        group_count += len(groups)
        for edge, p, q, limit in sorted(groups):
            i, j = edge
            e, f = body[i], body[j]
            threshold, divisor, P, Q, cross = tail_threshold(ruler, e, f, p, q, limit)
            threshold_hist[threshold] += 1
            threshold_row = (
                threshold, body, ruler, cell, edge, p, q, limit,
                e, f, divisor, P, Q, cross,
            )
            if max_threshold is None or threshold_row > max_threshold:
                max_threshold = threshold_row
            if threshold > 2:
                exception_groups += 1
            for scale in range(2, threshold):
                mass = m.physical_mass(ruler, cell, e, scale * p, f, scale * q)
                margin = mass - limit / 2
                require(margin > 0, ("half-limit head", body, edge, p, q, scale, margin))
                head_checks += 1
                row = (
                    margin, body, ruler, cell, edge, p, q, scale,
                    mass, limit, threshold,
                )
                if weakest_head is None or row < weakest_head:
                    weakest_head = row
                semantic.update(
                    repr((body, ruler, cell, edge, p, q, limit, threshold, scale, mass, margin)).encode()
                )

    actual_hist = tuple(sorted(threshold_hist.items()))
    actual_digest = semantic.hexdigest()
    require(group_count == 140711, group_count)
    require(actual_hist == EXPECTED_THRESHOLD_HIST, actual_hist)
    require(exception_groups == 42, exception_groups)
    require(head_checks == 66, head_checks)
    require(max_threshold == EXPECTED_MAX_THRESHOLD, max_threshold)
    require(weakest_head == EXPECTED_WEAKEST_HEAD, weakest_head)
    require(actual_digest == EXPECTED_HEAD_DIGEST, actual_digest)

    print("LRC14 reflected low-two-star selected half-limit all-dilations audit")
    print("status=PROVED analytic midpoint bound + VERIFIED-EXACT finite head")
    print(f"bodies={len(bodies)};rays={len(rays)};assignments={len(bodies)*len(rays)}")
    print(f"selected_groups={group_count}")
    print(f"tail_threshold_hist={actual_hist}")
    print(f"exception_groups={exception_groups};head_checks={head_checks};failures=0")
    print(f"max_threshold_row={max_threshold}")
    print(f"weakest_head_row={weakest_head}")
    print(f"head_semantic_sha256={actual_digest}")
    print("conclusion=every selected pair has I_s>=I_infinity/2 for every integer s>=2")


if __name__ == "__main__":
    main()
