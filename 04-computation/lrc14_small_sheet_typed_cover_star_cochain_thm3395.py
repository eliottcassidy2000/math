#!/usr/bin/env python3
"""Exact companion for THM-3395: small-sheet typed covers.

This is a standard-library-only, exact-integer audit.  It constructs the
typed coset-cover/star-cochain certificate without using the event sweep,
then compares every LRC(14)-relevant transverse subset of size at most five
for q=2,...,7 with the independent THM-3387 rational event calculation.

No assertion is used: all checks remain live under ``python -O``.
"""

from __future__ import annotations

import hashlib
import importlib.util
import itertools
import math
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]

DEPENDENCY_PINS = {
    "01-canon/theorems/THM-3387-exact-cyclic-sheet-cover-atlas-and-q2-gcd-graph.md":
        "c540255a185efb54c67035d69f3fbd94f4c1ad3e30c4e31738a8800e81198613",
    "04-computation/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.py":
        "9b0b46874a569d674b937b37cf74a8985fca2b77e3e480a75fb4924ea602f25a",
    "05-knowledge/results/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.out":
        "b4d9ce439bab4501bfd5e2cf13eb0b0e3685b7364f30e43b7d5ca9138d25cb5c",
}

CONVERGENCE_PINS = {
    "04-computation/lrc14_q6_typed_cover_clutter_probe_20260814.py":
        "e16ae18816b0620dde6ae05db01846ce49a37c292c7cccc76532ec2fc9c728d6",
    "05-knowledge/results/lrc14_q6_typed_cover_clutter_probe_20260814.out":
        "3628260d408258b1f9ad385c23c4d4bdea8cf0c7bad65f02dbbff9dc0e5de889",
}

EXPECTED_PROFILES = {
    2: (1, 7, 3, 0, 0, 0, 0, 0),
    3: (1, 10, 45, 72, 38, 6, 0, 0, 0, 0, 0),
    4: (1, 11, 51, 118, 123, 44, 3, 0, 0, 0, 0, 0),
    5: (1, 12, 66, 220, 495, 561, 268, 45, 1, 0, 0, 0, 0),
    6: (1, 12, 66, 217, 441, 515, 304, 76, 5, 0, 0, 0, 0),
    7: (1, 12, 66, 220, 495, 792),
}

EXPECTED_EDGE_COUNTS = {
    2: ((2, 18),),
    3: ((3, 48),),
    4: ((2, 4), (3, 15), (4, 17)),
    5: ((5, 231),),
    6: ((3, 3), (4, 29), (5, 7)),
    7: (),
}

EXPECTED_EDGE_SHA256 = {
    5: "296a57df572b22ef4ed23db9996da2c7104091e7300beb6026473cc7e12c7ed2",
    6: "f52aa82e9179fb04a4fc1c53ff96035edb538c22c06d5ac5bf8f41a22f5f78c1",
}

EXPECTED_PROFILE_SHA256 = {
    5: "43f234455bfb468cfa6703265b05910c80b415bee184950a4602709d5c5c7524",
    6: "949a5ffe61bd034887d4d3795c312408e27a83e528e14ac4a0882949a48801af",
}

EXPECTED_GLOBAL_SAFE = {2: 252, 3: 585, 4: 619, 5: 1617, 6: 1471, 7: 2079}
EXPECTED_CORE_RESCUES = {2: 0, 3: 3, 4: 0, 5: 2, 6: 7, 7: 0}
EXPECTED_EXACT_ROWS = {2: 252, 3: 588, 4: 619, 5: 1619, 6: 1478, 7: 2079}
EXPECTED_COMPARISON_COUNT = 6540
EXPECTED_SEMANTIC_DIGEST = "a47bf86457583259aa3712bfd5ca849328a7fa7b2db8e194345a910abb6610c1"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256_lf(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def check_pins(pins: dict[str, str], label: str) -> None:
    for rel, expected in pins.items():
        actual = sha256_lf(ROOT / rel)
        require(actual == expected, f"{label} hash mismatch for {rel}: {actual}")


def import_thm3387():
    path = ROOT / "04-computation/lrc14_exact_cyclic_sheet_cover_atlas_thm3387.py"
    spec = importlib.util.spec_from_file_location("thm3387_event_sweep", path)
    require(spec is not None and spec.loader is not None, "cannot load THM-3387 module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest_repr(value: object) -> str:
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def transverse_pool(q: int) -> tuple[int, ...]:
    return tuple(u for u in range(1, 15) if u % q != 0)


def blocked_cosets(q: int, u: int) -> tuple[tuple[int, frozenset[int]], ...]:
    """Canonical labels and kernel cosets for multiplication by u on Z/q."""
    seen: set[frozenset[int]] = set()
    answer = []
    for k in range(q):
        block = frozenset(ell for ell in range(q) if u * (ell - k) % q == 0)
        if block not in seen:
            seen.add(block)
            answer.append((min(block), block))
    return tuple(answer)


def centered_residues(residue: int, modulus: int, strict_bound: int) -> tuple[int, ...]:
    first = residue % modulus
    lower = (-strict_bound - first + modulus - 1) // modulus
    upper = (strict_bound - first) // modulus
    return tuple(first + modulus * n for n in range(lower, upper + 1))


def gap_values(q: int, u: int, v: int, k_u: int, k_v: int) -> tuple[int, ...]:
    """Allowed integer p_uv=q*u*v*(x_u-x_v), with the strict overlap bound."""
    modulus = q * math.gcd(u, v)
    residue = (k_v - k_u) * u * v
    strict_bound = (q * (u + v) - 1) // 14
    return centered_residues(residue, modulus, strict_bound)


def star_cochain(
    q: int, speeds: Sequence[int], labels: Sequence[int]
) -> Optional[tuple[tuple[int, int, int], ...]]:
    """Find a complete closed gap cochain, using vertex 0 as the star."""
    n = len(speeds)
    if n <= 1:
        return ()
    star_fibres = [gap_values(q, speeds[0], speeds[j], labels[0], labels[j])
                   for j in range(1, n)]
    if any(not fibre for fibre in star_fibres):
        return None
    for star in itertools.product(*star_fibres):
        edge_values: dict[tuple[int, int], int] = {}
        for j, value in enumerate(star, start=1):
            edge_values[(0, j)] = value
        compatible = True
        for i in range(1, n):
            for j in range(i + 1, n):
                numerator = speeds[i] * star[j - 1] - speeds[j] * star[i - 1]
                if numerator % speeds[0] != 0:
                    compatible = False
                    break
                value = numerator // speeds[0]
                if value not in gap_values(q, speeds[i], speeds[j], labels[i], labels[j]):
                    compatible = False
                    break
                edge_values[(i, j)] = value
            if not compatible:
                break
        if compatible:
            return tuple((i, j, edge_values[(i, j)])
                         for i in range(n) for j in range(i + 1, n))
    return None


def typed_cover_witness(
    q: int, speeds: Sequence[int]
) -> Optional[tuple[tuple[int, frozenset[int]], ...]]:
    """Find cosets covering Z/q whose gap fibres admit one global star cochain."""
    choices = tuple(blocked_cosets(q, u) for u in speeds)
    if sum(max(len(block) for _, block in options) for options in choices) < q:
        return None
    universe = frozenset(range(q))
    for typed_blocks in itertools.product(*choices):
        if frozenset().union(*(block for _, block in typed_blocks)) != universe:
            continue
        labels = tuple(label for label, _ in typed_blocks)
        if star_cochain(q, speeds, labels) is not None:
            return typed_blocks
    return None


def minimal_edges(q: int, pool: Sequence[int], max_rank: int = 5) -> tuple[tuple[int, ...], ...]:
    edges: list[tuple[int, ...]] = []
    for rank in range(1, min(max_rank, len(pool)) + 1):
        for subset in itertools.combinations(pool, rank):
            subset_set = frozenset(subset)
            if any(frozenset(edge) <= subset_set for edge in edges):
                continue
            if typed_cover_witness(q, subset) is not None:
                edges.append(subset)
    return tuple(edges)


def is_independent(subset: Iterable[int], edges: Sequence[Sequence[int]]) -> bool:
    subset_set = frozenset(subset)
    return not any(frozenset(edge) <= subset_set for edge in edges)


def independence_profile(pool: Sequence[int], edges: Sequence[Sequence[int]]) -> tuple[int, ...]:
    return tuple(sum(is_independent(subset, edges)
                     for subset in itertools.combinations(pool, rank))
                 for rank in range(len(pool) + 1))


def edge_rank_profile(edges: Sequence[Sequence[int]]) -> tuple[tuple[int, int], ...]:
    return tuple((rank, sum(len(edge) == rank for edge in edges))
                 for rank in sorted({len(edge) for edge in edges}))


def global_safe_six_body(q: int, independent_counts: Sequence[int]) -> int:
    core_count = 14 // q
    # Atlas rows have a nonempty core and a nonempty transverse part.
    return sum(math.comb(core_count, c) * independent_counts[6 - c]
               for c in range(1, min(core_count, 5) + 1))


def verify_fixed_labels(
    q: int, speeds: Sequence[int], labels: Sequence[int], expect_star: bool
) -> tuple[tuple[int, int, tuple[int, ...]], ...]:
    pair_fibres = tuple((i, j, gap_values(q, speeds[i], speeds[j], labels[i], labels[j]))
                        for i in range(len(speeds)) for j in range(i + 1, len(speeds)))
    require(all(fibre for _, _, fibre in pair_fibres), "fixed-label pair fibre unexpectedly empty")
    has_star = star_cochain(q, speeds, labels) is not None
    require(has_star == expect_star, f"fixed-label star expectation failed for q={q}")
    return pair_fibres


def strict_boundary_checks() -> tuple[int, int, int]:
    q7_points = tuple(Fraction(1, 14) + Fraction(k, 7) for k in range(7))
    strict_hits = sum(min(point % 1, 1 - point % 1) < Fraction(1, 14)
                      for point in q7_points)
    closed_hits = sum(min(point % 1, 1 - point % 1) <= Fraction(1, 14)
                      for point in q7_points)
    q8_points = (Fraction(-1, 16), Fraction(1, 16))
    q8_strict_hits = sum(min(point % 1, 1 - point % 1) < Fraction(1, 14)
                         for point in q8_points)
    require((strict_hits, closed_hits, q8_strict_hits) == (0, 2, 2),
            "strict-boundary audit failed")
    return strict_hits, closed_hits, q8_strict_hits


def main() -> None:
    check_pins(DEPENDENCY_PINS, "dependency")
    check_pins(CONVERGENCE_PINS, "convergence")
    event_sweep = import_thm3387()

    all_edges: dict[int, tuple[tuple[int, ...], ...]] = {}
    all_profiles: dict[int, tuple[int, ...]] = {}
    comparison_count = 0
    comparison_failures = []

    for q in range(2, 8):
        pool = transverse_pool(q)
        # An inclusion-minimal cover has at most q owners (give each owner a
        # private sheet).  For q<=6 this is the full clutter; at q=7 only
        # ranks <=5 occur in a six-body row with a nonempty core.
        edge_limit = 5 if q == 7 else q
        edges = minimal_edges(q, pool, edge_limit)
        full_profile = independence_profile(pool, edges)
        reported_profile = full_profile if q <= 6 else full_profile[:6]
        ranks = edge_rank_profile(edges)
        require(reported_profile == EXPECTED_PROFILES[q], f"q={q} profile mismatch")
        require(ranks == EXPECTED_EDGE_COUNTS[q], f"q={q} edge-rank mismatch: {ranks}")
        all_edges[q] = edges
        all_profiles[q] = reported_profile

        for rank in range(0, min(5, len(pool)) + 1):
            for subset in itertools.combinations(pool, rank):
                structural_cover = not is_independent(subset, edges)
                event_cover = False if not subset else not event_sweep.global_transverse_survives(q, subset)
                comparison_count += 1
                if structural_cover != event_cover:
                    comparison_failures.append((q, subset, structural_cover, event_cover))

        safe = global_safe_six_body(q, full_profile)
        require(safe == EXPECTED_GLOBAL_SAFE[q], f"q={q} global-safe count mismatch")
        require(safe + EXPECTED_CORE_RESCUES[q] == EXPECTED_EXACT_ROWS[q],
                f"q={q} exact-row decomposition mismatch")

    require(comparison_count == EXPECTED_COMPARISON_COUNT, "comparison-count mismatch")
    require(not comparison_failures, f"event/structure discrepancies: {comparison_failures[:5]}")

    for q in (5, 6):
        require(digest_repr(all_edges[q]) == EXPECTED_EDGE_SHA256[q], f"q={q} edge digest mismatch")
        require(digest_repr(all_profiles[q]) == EXPECTED_PROFILE_SHA256[q],
                f"q={q} profile digest mismatch")

    q6_hostile = verify_fixed_labels(6, (2, 8, 14), (0, 1, 2), False)
    q6_positive = verify_fixed_labels(6, (2, 8, 10), (0, 1, 2), True)
    q5_hostile = verify_fixed_labels(5, (1, 2, 3, 4, 9), (0, 2, 3, 1, 4), False)
    q5_positive = verify_fixed_labels(5, (1, 2, 3, 4, 7), (0, 2, 3, 1, 4), True)
    boundary = strict_boundary_checks()

    semantic_payload = (
        tuple((q, transverse_pool(q), all_profiles[q], edge_rank_profile(all_edges[q]),
               EXPECTED_GLOBAL_SAFE[q], EXPECTED_CORE_RESCUES[q], EXPECTED_EXACT_ROWS[q])
              for q in range(2, 8)),
        tuple((q, digest_repr(all_edges[q])) for q in (5, 6)),
        q6_hostile,
        q6_positive,
        q5_hostile,
        q5_positive,
        boundary,
        comparison_count,
    )
    semantic_digest = digest_repr(semantic_payload)
    if EXPECTED_SEMANTIC_DIGEST:
        require(semantic_digest == EXPECTED_SEMANTIC_DIGEST, "semantic digest mismatch")

    print("THM-3395 small-sheet typed-cover/star-cochain exact audit")
    print("dependency pins: OK")
    print("independent q=6 convergence pins: OK")
    for q in range(2, 8):
        print(
            f"q={q}: V={transverse_pool(q)}; minimal-edge-ranks={edge_rank_profile(all_edges[q])}; "
            f"I={all_profiles[q]}; global-safe={EXPECTED_GLOBAL_SAFE[q]}; "
            f"core-rescues={EXPECTED_CORE_RESCUES[q]}; exact={EXPECTED_EXACT_ROWS[q]}"
        )
    print(f"q=5 edge SHA256: {digest_repr(all_edges[5])}")
    print(f"q=5 profile SHA256: {digest_repr(all_profiles[5])}")
    print(f"q=6 edge SHA256: {digest_repr(all_edges[6])}")
    print(f"q=6 profile SHA256: {digest_repr(all_profiles[6])}")
    print(f"q=6 hostile (2,8,14) pair fibres: {q6_hostile}; star=EMPTY")
    print(f"q=6 positive (2,8,10) pair fibres: {q6_positive}; star=NONEMPTY")
    print(f"q=5 hostile (1,2,3,4,9) pair fibres: {q5_hostile}; star=EMPTY")
    print(f"q=5 positive (1,2,3,4,7) pair fibres: {q5_positive}; star=NONEMPTY")
    print(f"strict boundary (q7 strict, q7 closed, q8 strict two-sheet): {boundary}")
    print(f"typed criterion versus independent event sweep: {comparison_count} subsets, 0 discrepancies")
    print(f"semantic SHA256: {semantic_digest}")
    print("all exact checks: PASS")


if __name__ == "__main__":
    main()
