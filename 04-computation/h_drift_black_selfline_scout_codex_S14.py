#!/usr/bin/env python3
"""Exact H-drift formula and first untested black-self-line census.

This is a bounded companion to HYP-6900 and the opus-S311 HYP-6890 entry.

Part A checks the elementary endpoint-coefficient identity

    sum_e (H(T^e)-H(T)) = a_{n-2}(T) - (n-1)H(T),

where a_{n-2} counts vertex orders with exactly one backward adjacency.  In
OCF coordinates this is

    E[Delta H | T] = 2/C(n,2) * (K(T)-nH(T)),

    K(T) = sum_I 2^(|I|+f_I),

with I ranging over vertex-disjoint directed odd-cycle collections and
f_I=n-1-sum_{C in I}(|C|-1).  Thus K is a length-weighted hard-core
partition function, not a function of H alone.

Part B reuses the certified A000568(8)=6880 fixed-path classifier from
metagraph_flow_n8_check_opus_S305.py and counts ordinary (not
converse-merged) complement self-lines.  This is the first exact n=8 test of
the proposed all-n identity 2*black_self_lines=SC(n).

The script deliberately does not claim a theorem or hypothesis identifier.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import itertools
import runpy
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb
from pathlib import Path

import numpy as np


LEGACY_N8 = Path("04-computation/metagraph_flow_n8_check_opus_S305.py")
ATLAS_EXPORT = Path("04-computation/mobius_cech_n8_node_atlas_export_codex_S13.py")
THM062 = Path("01-canon/theorems/THM-062-forward-edge-distribution.md")
THM830 = Path("01-canon/theorems/THM-830-b3-deletion-deck-mirror-current-calculus.md")
THM833 = Path("01-canon/theorems/THM-833-ornstein-uhlenbeck-law-of-the-wiggly-layer.md")
CERTIFIED_MERGED_ATLAS_SHA256 = (
    "30debad3387a4ea0ef51108ea132115efda2ac2fcdfcc2c5c1d4d23155095835"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fixed_path_tournament(n: int, mask: int) -> list[list[bool]]:
    tiles = [
        (x, y)
        for y in range(1, n - 1)
        for x in range(n, y + 1, -1)
        if x - y >= 2
    ]
    adjacency = [[False] * n for _ in range(n)]
    for x in range(2, n + 1):
        adjacency[x - 1][x - 2] = True
    for bit, (x, y) in enumerate(tiles):
        if (mask >> bit) & 1:
            adjacency[x - 1][y - 1] = True
        else:
            adjacency[y - 1][x - 1] = True
    return adjacency


def hamiltonian_count(adjacency: list[list[bool]]) -> int:
    n = len(adjacency)
    dp = [[0] * n for _ in range(1 << n)]
    for vertex in range(n):
        dp[1 << vertex][vertex] = 1
    for used in range(1 << n):
        for last in range(n):
            count = dp[used][last]
            if not count:
                continue
            for nxt in range(n):
                if not (used >> nxt) & 1 and adjacency[last][nxt]:
                    dp[used | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def forward_distribution(adjacency: list[list[bool]]) -> list[int]:
    n = len(adjacency)
    answer = [0] * n
    for order in itertools.permutations(range(n)):
        forward = sum(adjacency[order[i]][order[i + 1]] for i in range(n - 1))
        answer[forward] += 1
    return answer


def directed_cycles(adjacency: list[list[bool]], length: int) -> set[tuple[int, ...]]:
    n = len(adjacency)
    cycles: set[tuple[int, ...]] = set()
    for vertices in itertools.combinations(range(n), length):
        anchor = min(vertices)
        rest = [vertex for vertex in vertices if vertex != anchor]
        for tail in itertools.permutations(rest):
            cycle = (anchor,) + tail
            if all(adjacency[cycle[i]][cycle[(i + 1) % length]] for i in range(length)):
                cycles.add(cycle)
    return cycles


def small_ocf_data(adjacency: list[list[bool]]) -> tuple[int, int, int]:
    triangles = directed_cycles(adjacency, 3)
    fives = directed_cycles(adjacency, 5) if len(adjacency) >= 5 else set()
    disjoint_triangle_pairs = sum(
        set(left).isdisjoint(right)
        for left, right in itertools.combinations(triangles, 2)
    )
    return len(triangles), len(fives), disjoint_triangle_pairs


def direct_drift(adjacency: list[list[bool]], h_value: int) -> Fraction:
    n = len(adjacency)
    total = 0
    for u in range(n):
        for v in range(u + 1, n):
            adjacency[u][v], adjacency[v][u] = adjacency[v][u], adjacency[u][v]
            total += hamiltonian_count(adjacency) - h_value
            adjacency[u][v], adjacency[v][u] = adjacency[v][u], adjacency[u][v]
    return Fraction(total, comb(n, 2))


def drift_scout() -> dict[int, dict[str, object]]:
    summaries: dict[int, dict[str, object]] = {}
    for n in range(4, 7):
        dimension = comb(n - 1, 2)
        by_h: defaultdict[int, set[Fraction]] = defaultdict(set)
        rows = 0
        for mask in range(1 << dimension):
            adjacency = fixed_path_tournament(n, mask)
            h_value = hamiltonian_count(adjacency)
            distribution = forward_distribution(adjacency)
            one_defect = distribution[n - 2]
            drift_from_endpoint = Fraction(
                one_defect - (n - 1) * h_value, comb(n, 2)
            )
            drift_from_flips = direct_drift(adjacency, h_value)
            if drift_from_endpoint != drift_from_flips:
                raise RuntimeError("endpoint coefficient does not equal direct flip drift")

            t3, t5, alpha33 = small_ocf_data(adjacency)
            h_ocf = 1 + 2 * t3 + 2 * t5 + 4 * alpha33
            if h_value != h_ocf:
                raise RuntimeError("small OCF formula failed")
            # For n<=6 these are all possible independent odd-cycle types.
            k_value = (
                2 ** (n - 1)
                + 2 ** (n - 2) * t3
                + (2 ** (n - 4) * t5 if n >= 5 else 0)
                + (2 ** (n - 3) * alpha33 if n >= 6 else 0)
            )
            drift_from_ocf = Fraction(2 * (k_value - n * h_value), comb(n, 2))
            if drift_from_endpoint != drift_from_ocf:
                raise RuntimeError("length-weighted OCF drift formula failed")
            if n == 4 and drift_from_endpoint != Fraction(2 * (3 - h_value), 3):
                raise RuntimeError("n=4 affine H-drift formula failed")
            if n == 5 and drift_from_endpoint != Fraction(12 - h_value - 6 * t5, 5):
                raise RuntimeError("n=5 (H,t5) drift formula failed")
            by_h[h_value].add(drift_from_endpoint)
            rows += 1

        split = {
            h_value: tuple(sorted(values))
            for h_value, values in sorted(by_h.items())
            if len(values) > 1
        }
        summaries[n] = {
            "fixed_path_rows": rows,
            "H_fibres": len(by_h),
            "split_H_fibres": split,
            "split_count": len(split),
        }
    return summaries


def reflect_masks(bits: np.ndarray, reflection: np.ndarray, dimension: int) -> np.ndarray:
    answer = np.zeros(len(bits), dtype=np.int64)
    for source in range(dimension):
        answer |= ((bits >> source) & 1) << int(reflection[source])
    return answer


def n8_selfline_scout() -> dict[str, object]:
    with contextlib.redirect_stdout(io.StringIO()):
        namespace = runpy.run_path(str(LEGACY_N8))
    bits = np.asarray(namespace["bits_all"], dtype=np.int64)
    class_of = np.asarray(namespace["cls_of"], dtype=np.int32)
    representatives = np.asarray(namespace["rep"], dtype=np.int64)
    reflection = np.asarray(namespace["refl"], dtype=np.int64)
    dimension = int(namespace["m"])
    full = (1 << dimension) - 1
    reflected = reflect_masks(bits, reflection, dimension)
    grid_symmetric = reflected == bits
    partner = bits ^ full
    quasi_fixed = class_of == class_of[partner]
    black_quasi = quasi_fixed & ~grid_symmetric
    blue_quasi = quasi_fixed & grid_symmetric

    reflected_representatives = reflect_masks(representatives, reflection, dimension)
    converse_class = class_of[reflected_representatives]
    is_self_converse = converse_class == np.arange(len(representatives))
    sc_count = int(is_self_converse.sum())
    sc_from_grid = len(set(map(int, class_of[grid_symmetric])))
    if sc_count != sc_from_grid:
        raise RuntimeError("self-converse class counts disagree")

    line_side = bits < partner
    black_selfline_side = line_side & black_quasi
    blue_selfline_side = line_side & blue_quasi
    black_self_lines = int(black_selfline_side.sum())
    blue_self_lines = int(blue_selfline_side.sum())
    black_quasi_count = int(black_quasi.sum())
    blue_quasi_count = int(blue_quasi.sum())

    # The black locus has a free Klein-four action generated by flip and
    # reflection.  Canonical representatives certify the orbit count directly.
    black_masks = bits[black_quasi]
    black_reflected = reflected[black_quasi]
    orbit_representatives = np.minimum.reduce(
        (black_masks, black_masks ^ full, black_reflected, black_reflected ^ full)
    )
    klein_orbits = len(set(map(int, orbit_representatives)))
    if 4 * klein_orbits != black_quasi_count:
        raise RuntimeError("black quasi-fixed Klein-four action is not free")

    carrier_classes = class_of[black_selfline_side]
    carrier_multiplicity = Counter(map(int, carrier_classes))
    carrier_histogram = Counter(carrier_multiplicity.values())
    sc_carriers = sum(bool(is_self_converse[carrier]) for carrier in carrier_multiplicity)
    non_sc_carriers = len(carrier_multiplicity) - sc_carriers

    # A V4 orbit has two black self-lines.  They are carried by one SC class
    # twice, or by a converse pair of non-SC classes once each.
    v4_sc_carrier_orbits = 0
    v4_non_sc_carrier_orbits = 0
    for representative in set(map(int, orbit_representatives)):
        carrier = int(class_of[representative])
        if is_self_converse[carrier]:
            v4_sc_carrier_orbits += 1
        else:
            v4_non_sc_carrier_orbits += 1
    if v4_sc_carrier_orbits + v4_non_sc_carrier_orbits != klein_orbits:
        raise RuntimeError("Klein-four carrier accounting failed")

    return {
        "fixed_path_tilings": len(bits),
        "ordinary_classes": len(representatives),
        "self_converse_classes": sc_count,
        "grid_symmetric_tilings": int(grid_symmetric.sum()),
        "quasi_fixed_tilings": int(quasi_fixed.sum()),
        "blue_quasi_fixed_tilings": blue_quasi_count,
        "black_quasi_fixed_tilings": black_quasi_count,
        "blue_self_lines": blue_self_lines,
        "black_self_lines": black_self_lines,
        "two_black_self_lines_equals_SC": 2 * black_self_lines == sc_count,
        "black_quasi_fixed_equals_SC": black_quasi_count == sc_count,
        "free_black_Klein_four_orbits": klein_orbits,
        "SC_divided_by_four": Fraction(sc_count, 4),
        "carrier_classes": len(carrier_multiplicity),
        "self_converse_carrier_classes": sc_carriers,
        "non_self_converse_carrier_classes": non_sc_carriers,
        "self_lines_per_carrier_histogram": dict(sorted(carrier_histogram.items())),
        "Klein_four_orbits_with_SC_carrier": v4_sc_carrier_orbits,
        "Klein_four_orbits_with_non_SC_converse_pair": v4_non_sc_carrier_orbits,
    }


def main() -> None:
    print("=" * 88)
    print("PART A -- exact H-drift endpoint and length-weighted OCF formula")
    print("=" * 88)
    summaries = drift_scout()
    print("identity: E[Delta H|T] = (a_(n-2)-(n-1)H)/C(n,2)")
    print("OCF form: E[Delta H|T] = 2*(K-nH)/C(n,2)")
    print("K = sum_I 2^(|I|+f_I) = 2^(n-1) sum_I prod_(C in I) 2^(2-|C|)")
    print("n=4 specialization: E[Delta H|T] = 2*(3-H)/3")
    print("n=5 specialization: E[Delta H|T] = (12-H-6*t5)/5")
    for n, summary in summaries.items():
        print(
            f"n={n}: rows={summary['fixed_path_rows']} H-fibres={summary['H_fibres']} "
            f"split-fibres={summary['split_count']} {summary['split_H_fibres']}"
        )
    print("first H-only obstruction: n=5, H=15; t5=2 gives -3, t5=3 gives -21/5")
    print("proof: each Hamiltonian path is destroyed once for each of its n-1 arcs;")
    print("orders with one backward adjacency are created by flipping their unique bad arc.")
    print("THM-062 gives a_(n-2)=sum_I 2^|I|*(2^(f_I+1)-n-1),")
    print("while H=sum_I 2^|I|; subtraction gives the displayed OCF formula.")
    print("n=6 obstruction to pointwise mean reversion: E[H]=45/2 but H=23 has drift 4/3 or -4/15.")

    print()
    print("=" * 88)
    print("PART B -- n=8 corrected black-self-line law")
    print("=" * 88)
    result = n8_selfline_scout()
    for key, value in result.items():
        print(f"{key}={value}")
    print("Klein-four target: SC classes -> black self-lines must be 2-to-1;")
    print("after quotienting black lines by reflection, four SC classes are required per free V4 orbit.")
    print()
    print("ASSUMPTION CHALLENGE / TOURNAMENT ANALYSIS")
    print("drift vertices are odd-cycle collections and endpoint coefficients, not runners or raw arcs;")
    print("self-line vertices are ordinary-class equalities, kappa-lines, and free V4 orbits.")
    print("H alone destroys cycle-length occupancy; the converse-merged atlas destroys same-vs-converse.")
    print("No joint carrier tournament is formed: drift exactness and ordinary-class equality are different predicates.")
    print("Within the drift predicate, the exactness/compression carrier order is")
    print("length-weighted OCF K > endpoint pair (H,a_(n-2)) > direct flip ledger > H-only;")
    print("fingerprint: transitive score histogram {0:1,1:1,2:1,3:1}, cycles=0, SCCs=4, HP=1.")
    print()
    print("PROVENANCE")
    print(f"legacy_n8_classifier_sha256={sha256(LEGACY_N8)}")
    print(f"n8_atlas_export_sha256={sha256(ATLAS_EXPORT)}")
    print(f"THM062_sha256={sha256(THM062)}")
    print(f"THM830_sha256={sha256(THM830)}")
    print(f"THM833_sha256={sha256(THM833)}")
    print(f"certified_merged_atlas_sha256={CERTIFIED_MERGED_ATLAS_SHA256}")
    print(f"script_sha256={sha256(Path(__file__))}")


if __name__ == "__main__":
    main()
