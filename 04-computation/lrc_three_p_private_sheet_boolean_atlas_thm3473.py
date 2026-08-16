#!/usr/bin/env python3
"""Exact deterministic companion for independently audited THM-3473.

For p=14k-1 and q=3p this program compares two independent descriptions of
the THM-3469 eight-owner cover: direct strict cyclic-distance masks and the
nearest-p chart.  It verifies every private-sheet set, the complete Boolean
support spectrum, all deletion deficits, and the three residue-dependent
count formulae on a large finite window.  The theorem supplies the universal
proof; this companion is a hostile exact audit, not an extrapolation engine.
"""

from __future__ import annotations

import ast
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from json import dumps
from math import floor
from pathlib import Path
import sys


THEOREM_ID = "THM-3473"
ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = (
    (
        "THM-3469",
        "01-canon/theorems/THM-3469-three-times-p-half-twist-eight-owner-cover-boundary.md",
        "2d86e9ecf8fc14ef58a1bd0fc7a000b1fb49b42718d10d96637a3908b90f486b",
    ),
)
K_MAX = 500
FRACTION_K_MAX = 60
EXPECTED_CELLS = 5_259_000
EXPECTED_CLASSIFICATION_SHA256 = "0ded756ca28ea42ee68990cc3bb8a314d01028da4ff1727f524194efe72aab9a"
EXPECTED_SEMANTIC_SHA256 = "1a269572f813dcfafa4d99b9dda39e6fd36514005ff671ef7eec3f6416467fcc"

# Owners are indexed in the displayed THM-3469 order.
OWNER_LABELS = (
    "1", "p-1", "p+1", "2p-1", "2p", "2p+1", "3p-6", "3p-1"
)
CORE_CHANNEL = {
    (0, 1): 0,
    (0, 5): 0,
    (1, 1): 1,
    (5, 5): 1,
    (1, 5): 2,
    (5, 1): 2,
    (2, 1): 3,
    (4, 5): 3,
    (2, 5): 5,
    (4, 1): 5,
    (3, 1): 7,
    (3, 5): 7,
}
PACKET_A = (0, 3, 4, 5)
PACKET_B = (1, 2, 4, 7)
PACKET_C = (4, 6)


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def template(p: int) -> tuple[int, ...]:
    return (1, p - 1, p + 1, 2 * p - 1, 2 * p, 2 * p + 1, 3 * p - 6, 3 * p - 1)


def cyclic_distance(word: int, modulus: int) -> int:
    residue = word % modulus
    return min(residue, modulus - residue)


def direct_support(p: int, sheet: int) -> tuple[int, ...]:
    q = 3 * p
    x = 2 * sheet + 1
    return tuple(
        index
        for index, owner in enumerate(template(p))
        if 14 * cyclic_distance(owner * x, 2 * q) < 2 * q
    )


def fraction_support(p: int, sheet: int) -> tuple[int, ...]:
    q = 3 * p
    x = 2 * sheet + 1
    support = []
    for index, owner in enumerate(template(p)):
        value = Fraction(owner * x, 2 * q)
        residue = value - floor(value)
        if min(residue, 1 - residue) < Fraction(1, 14):
            support.append(index)
    return tuple(support)


def nearest_chart(p: int, x: int) -> tuple[int, int]:
    """Unique x=mp+y modulo 6p with m in Z/6 and |y|<p/2."""

    require(p % 2 == 1 and 0 <= x < 6 * p, ("chart domain", p, x))
    radius = (p - 1) // 2
    candidates = []
    for m in range(6):
        raw = x - m * p
        for y in (raw - 6 * p, raw, raw + 6 * p):
            if abs(y) <= radius:
                candidates.append((m, y))
    require(len(candidates) == 1, ("chart uniqueness", p, x, candidates))
    return candidates[0]


def chart_support(p: int, sheet: int) -> tuple[int, ...]:
    q = 3 * p
    x = 2 * sheet + 1
    m, y = nearest_chart(p, x)
    sigma = x % 6
    h = (3 * p - 1) // 7
    require(sigma in (1, 3, 5), ("odd phase", p, sheet, sigma))
    if sigma in (1, 5):
        if abs(y) <= h:
            return (CORE_CHANNEL[(m, sigma)],)
        return (6,)
    if abs(y) <= h:
        if m == 0:
            return PACKET_A
        if m == 3:
            return PACKET_B
        return (4,)
    return PACKET_C


def private_count_formula(k: int) -> tuple[int, ...]:
    residue = k % 3
    return (
        4 * k,
        4 * k - 2 * (residue == 1),
        4 * k - 2 * (residue == 0),
        4 * k,
        8 * k - 2 * (residue == 2),
        4 * k,
        4 * k - 2 * (residue == 2),
        4 * k,
    )


def support_count_formula(k: int) -> dict[tuple[int, ...], int]:
    counts = {
        (index,): count
        for index, count in enumerate(private_count_formula(k))
    }
    counts[PACKET_A] = 2 * k
    counts[PACKET_B] = 2 * k - 1
    counts[PACKET_C] = 2 * k + 2 * (k % 3 == 2)
    return counts


def residue_lemma_audit() -> dict[str, object]:
    checked = 0
    for k in range(1, 10_001):
        h = 6 * k - 1
        counts = tuple(
            (h - residue) // 6 + (h + residue) // 6 + 1
            for residue in range(6)
        )
        require(counts == (2 * k - 1,) + (2 * k,) * 5, (k, counts))
        checked += 2 * h + 1
    return {
        "k_max": 10_000,
        "integer_points": checked,
        "formula": "N_0=2k-1; N_1=...=N_5=2k",
    }


def exact_atlas_audit() -> dict[str, object]:
    classification = sha256()
    total_cells = 0
    total_owner_tests = 0
    residue_class_summaries: dict[int, Counter[tuple[int, ...]]] = {
        0: Counter(), 1: Counter(), 2: Counter()
    }
    first_rows = []

    for k in range(1, K_MAX + 1):
        p = 14 * k - 1
        q = 3 * p
        observed: Counter[tuple[int, ...]] = Counter()
        private_counts = [0] * 8
        for sheet in range(q):
            direct = direct_support(p, sheet)
            chart = chart_support(p, sheet)
            require(direct == chart, ("support mismatch", k, sheet, direct, chart))
            require(direct, ("uncovered", k, sheet))
            observed[direct] += 1
            if len(direct) == 1:
                private_counts[direct[0]] += 1
            classification.update(k.to_bytes(2, "big"))
            classification.update(sheet.to_bytes(3, "big"))
            classification.update(sum(1 << index for index in direct).to_bytes(1, "big"))

        expected = support_count_formula(k)
        require(dict(observed) == expected, ("support census", k, observed, expected))
        require(tuple(private_counts) == private_count_formula(k), (k, private_counts))
        require(all(count > 0 for count in private_counts), ("irredundancy", k, private_counts))

        # Removing coordinate i uncovers exactly its private Boolean atom.
        deletion_deficits = tuple(observed[(index,)] for index in range(8))
        require(deletion_deficits == tuple(private_counts), ("deletion", k))

        multiplicities = Counter()
        for support, count in observed.items():
            multiplicities[len(support)] += count
        epsilon = int(k % 3 == 2)
        require(
            multiplicities == Counter({
                1: 36 * k - 2 - 2 * epsilon,
                2: 2 * k + 2 * epsilon,
                4: 4 * k - 1,
            }),
            ("multiplicity", k, multiplicities),
        )
        require(sum(multiplicities.values()) == q, ("partition", k))
        residue_class_summaries[k % 3].update(observed)
        if k <= 6:
            first_rows.append((k, tuple(private_counts), tuple(sorted(multiplicities.items()))))
        total_cells += q
        total_owner_tests += 8 * q

    require(total_cells == EXPECTED_CELLS, total_cells)
    digest = classification.hexdigest()
    if EXPECTED_CLASSIFICATION_SHA256 != "UNFROZEN":
        require(digest == EXPECTED_CLASSIFICATION_SHA256, ("classification drift", digest))
    return {
        "k_max": K_MAX,
        "cells": total_cells,
        "owner_tests": total_owner_tests,
        "classification_sha256": digest,
        "first_rows": tuple(first_rows),
        "support_types": tuple(
            tuple(index + 1 for index in support)
            for support in support_count_formula(3)
        ),
        "residue_summary_digest": sha256(
            dumps(
                {
                    residue: sorted((support, count) for support, count in summary.items())
                    for residue, summary in residue_class_summaries.items()
                },
                sort_keys=True,
                separators=(",", ":"),
            ).encode("ascii")
        ).hexdigest(),
    }


def fraction_reference_audit() -> dict[str, int]:
    cells = 0
    owner_tests = 0
    for k in range(1, FRACTION_K_MAX + 1):
        p = 14 * k - 1
        q = 3 * p
        for sheet in range(q):
            require(
                fraction_support(p, sheet) == direct_support(p, sheet),
                ("fraction route", k, sheet),
            )
            cells += 1
            owner_tests += 8
    return {"k_max": FRACTION_K_MAX, "cells": cells, "owner_tests": owner_tests}


def generalized_tournament_report() -> dict[str, object]:
    """The exact bidirected 2-section; orientations are deliberately absent."""

    packets = (PACKET_A, PACKET_B, PACKET_C)
    edges = set()
    for packet in packets:
        for left in packet:
            for right in packet:
                if left < right:
                    edges.add((left, right))
    require(len(edges) == 13, edges)
    degrees = tuple(sum(vertex in edge for edge in edges) for vertex in range(8))
    require(degrees == (3, 3, 3, 3, 7, 3, 1, 3), degrees)
    require(set(PACKET_A) & set(PACKET_B) == {4}, "four-packet hub")
    require(set(PACKET_A) & set(PACKET_C) == {4}, "A/C hub")
    require(set(PACKET_B) & set(PACKET_C) == {4}, "B/C hub")
    return {
        "packets_one_based": tuple(tuple(index + 1 for index in packet) for packet in packets),
        "undirected_edges": len(edges),
        "bidirected_arcs": 2 * len(edges),
        "degrees": degrees,
        "hub_one_based": 5,
        "warning": "2-section loses sheet counts, singleton atoms, and k mod 3",
    }


def u_spine_deficit_word_report() -> dict[str, object]:
    """Ternary private-deficit state on the THM-3469 Berggren lanes."""

    def symbol(t: int) -> int:
        q = (2 * t + 1) ** 2 + 2
        if q % 42 != 39:
            return -1
        k = (q + 3) // 42
        require(q == 42 * k - 3, ("U-spine lane", t, q, k))
        return k % 3

    period = 63
    word = tuple(symbol(t) for t in range(period))
    require(all(symbol(t) == symbol(t + period) for t in range(10 * period)), "period 63")
    counts = Counter(word)
    require(counts == Counter({-1: 51, 0: 4, 1: 4, 2: 4}), counts)
    witnesses = []
    for shift in (period // 3, period // 7):
        t = next(t for t in range(period) if symbol(t) != symbol(t + shift))
        witnesses.append((shift, t, symbol(t), t + shift, symbol(t + shift)))
    return {
        "period": period,
        "counts": tuple(sorted(counts.items())),
        "minimality_witnesses": tuple(witnesses),
        "ambient_each_state": "4/63",
        "conditioned_each_state": "1/3",
    }


def security_report(source: Path) -> tuple[str, ...]:
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert present")
    forbidden = {
        "eval", "exec", "compile", "open", "system", "popen", "run", "Popen",
        "write_text", "write_bytes", "unlink", "remove", "rename",
    }
    called = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    called.update(
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    )
    require(not called & forbidden, ("forbidden calls", sorted(called & forbidden)))
    imports = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imports.add(node.module.split(".")[0])
    allowed = {
        "__future__", "ast", "collections", "fractions", "hashlib", "json",
        "math", "pathlib", "sys",
    }
    require(imports <= allowed, ("unexpected imports", sorted(imports - allowed)))
    return tuple(sorted(imports))


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    dependency_hashes = tuple(
        (label, relative, lf_sha256(ROOT / relative))
        for label, relative, _expected in DEPENDENCIES
    )
    require(dependency_hashes == DEPENDENCIES, ("dependency drift", dependency_hashes))
    security = security_report(Path(__file__))
    residues = residue_lemma_audit()
    atlas = exact_atlas_audit()
    fractions = fraction_reference_audit()
    tournament = generalized_tournament_report()
    spine_word = u_spine_deficit_word_report()

    semantic_payload = {
        "theorem": THEOREM_ID,
        "dependencies": dependency_hashes,
        "residue_lemma": residues,
        "atlas": atlas,
        "fraction_reference": fractions,
        "generalized_tournament": tournament,
        "u_spine_deficit_word": spine_word,
        "private_formula": {
            residue: private_count_formula(3 + residue)
            for residue in range(3)
        },
        "multiplicity_formula": "1:36k-2-2e2; 2:2k+2e2; 4:4k-1",
        "asymptotic_support_densities": {
            "singleton": "6/7",
            "packet_A": "1/21",
            "packet_B": "1/21",
            "packet_C": "1/21",
        },
        "deficit_lane_harmonic_coefficient": "1/3",
        "security": security,
    }
    semantic_hash = sha256(
        dumps(semantic_payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "UNFROZEN":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, ("semantic drift", semantic_hash))

    print("THM-3473 EXACT DETERMINISTIC COMPANION")
    print("STATUS: PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED")
    for label, relative, digest in dependency_hashes:
        print(f"DEPENDENCY: {label} {digest} {relative}")
    print(f"SECURITY_IMPORTS: {','.join(security)}")
    print(f"RESIDUE_LEMMA: {residues['formula']} through k={residues['k_max']}")
    print(f"DIRECT_ATLAS_K1_TO{atlas['k_max']}: cells={atlas['cells']} owner_tests={atlas['owner_tests']}")
    print(f"CLASSIFICATION_SHA256: {atlas['classification_sha256']}")
    print(f"RESIDUE_SUMMARY_SHA256: {atlas['residue_summary_digest']}")
    print(f"FRACTION_REFERENCE_K1_TO{fractions['k_max']}: cells={fractions['cells']} owner_tests={fractions['owner_tests']}")
    print(f"BOOLEAN_SUPPORT_TYPES_ONE_BASED: {atlas['support_types']}")
    print("PRIVATE_COUNTS_OWNER_ORDER: (4k,4k-2e1,4k-2e0,4k,8k-2e2,4k,4k-2e2,4k)")
    print("MULTIPLICITY_COUNTS: weight1=36k-2-2e2; weight2=2k+2e2; weight4=4k-1")
    print(f"FIRST_SIX_ROWS: {atlas['first_rows']}")
    print(f"BIDIRECTED_TWO_SECTION: {tournament}")
    print(f"U_SPINE_DEFICIT_WORD: {spine_word}")
    print("ASYMPTOTIC_SHEET_DENSITIES: singleton=6/7; each multiple packet=1/21")
    print("DEFICIT_LANES: k mod 3; each natural/harmonic coefficient=1/3")
    print(f"SEMANTIC_SHA256: {semantic_hash}")
    print("VERDICT: exact 11-support Boolean atlas and eight positive deletion deficits for every k>=1; global rank minimality and LRC(14) do not follow")


if __name__ == "__main__":
    main()
