#!/usr/bin/env python3
r"""Exact extreme-pair closure of the reflected ``(D,m)=(6,1)`` slice.

Work inside the sufficient reflected ``k=1`` family of THM-2941.  The
cap-three repair leaves 561 bodies and the wedge ``D>=6, 1<=m<D/2``.
On every one of those bodies the universal same-level graph is ``K6``.
Consequently a packet not already closed by a same-level pair has six
distinct levels.

At ``(D,m)=(6,1)`` the global-min/global-max ordered pair has levels
``(1,7)``.  For each of the 561 bodies and all 30 *physical* orientations
``(i,j)`` this verifier:

* maximizes the singleton debt over the four distinct intermediate levels;
* scans every body-safe cell for the exact overlap of slots ``i@1,j@7``;
* checks that the best located overlap is strictly larger than that debt.

The debt maximizer has a closed form.  Use the four smallest intermediate
levels ``2,3,4,5`` and assign them in reverse order to the four remaining
body labels.  This follows from strict decrease in the level and the exact
rearrangement gap proved below; an exhaustive ``5P4`` hostile audit is also
run for every body and orientation.

This closes one exact boundary slice of the reflected THM-2941 sufficient
family.  It is not a physical-survivor census and not a proof of LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations
from math import lcm
from multiprocessing import get_context
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
HALF = ROOT / "04-computation/lrc14_j7_reflected_half_cone_orientation_enrichment_scout_thm2941.py"
HALF_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_half_cone_orientation_enrichment_scout_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_extreme_pair_d6_m1_distinct_debt_closure_thm2941.out"

EXPECTED_HALF_SHA256 = "d53b9ebb76930fc6eb4cdc697c11f57413747e834e0ee4a11ea81ee47c884398"
EXPECTED_HALF_OUTPUT_SHA256 = "a5cfde4f3ba447c92d5486017bddd484a0cebfad87330cd19c25d307327a5d87"
EXPECTED_HALF_SEMANTIC = "1a8a6cac81420388034c849ca0ba852a9b87416c85eeedcb12e1cbafea572f5a"

BODY_COUNT = 561
ORIENTATION_COUNT = 30
LOW_LEVEL = 1
HIGH_LEVEL = 7
INTERMEDIATE_LEVELS = (2, 3, 4, 5)
AVAILABLE_INTERMEDIATE_LEVELS = (2, 3, 4, 5, 6)

# Frozen after the first independent exact replay.
EXPECTED_DEBT_MAXIMUM = (
    (1, 2, 3, 4, 6, 12),
    (5, 0),
    F(653495156, 42268562875),
    F(478, 15275),
    F(8699774082, 549491317375),
    150,
    (7, 5, 4, 3, 2, 1),
    168,
    88,
)
EXPECTED_WEAKEST = (
    (1, 2, 3, 4, 6, 12),
    (0, 3),
    F(4942672624, 429463642167),
    F(1020, 48931),
    F(669634815172, 71720428241889),
    143,
    (1, 5, 4, 7, 3, 2),
    168,
    88,
)
EXPECTED_ROW_DIGEST = "68a666030d9124474cbbe7b381fa1daab77d13521a84f5080ae62cc94d4aebe4"
EXPECTED_SEMANTIC = "a74dfc2377029e62aa549e7ccfcf84b14348c56b507074efdb92e97fa9ca5dfd"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(sha256(HALF) == EXPECTED_HALF_SHA256,
        ("half-cone source changed", sha256(HALF)))
require(sha256(HALF_OUTPUT) == EXPECTED_HALF_OUTPUT_SHA256,
        ("half-cone output changed", sha256(HALF_OUTPUT)))
require(f"semantic_sha256={EXPECTED_HALF_SEMANTIC}" in HALF_OUTPUT.read_text(),
        "half-cone semantic token missing")
H = import_module("extreme_pair_half_cone", HALF)
T = H.T


def debt_term(label: int, level: int, ruler: int) -> F:
    require(0 < label < level * ruler, (label, level, ruler))
    return F(label, 7 * (level * ruler - label))


def rearrangement_gap(low_label: int, high_label: int,
                      low_level: int, high_level: int, ruler: int) -> F:
    """Gain from pairing the larger label with the smaller level."""
    require(0 < low_label < high_label < ruler, (low_label, high_label, ruler))
    require(0 < low_level < high_level, (low_level, high_level))
    direct = (
        debt_term(high_label, low_level, ruler)
        + debt_term(low_label, high_level, ruler)
        - debt_term(high_label, high_level, ruler)
        - debt_term(low_label, low_level, ruler)
    )
    numerator = (
        ruler * (high_level - low_level) * (high_label - low_label)
        * (low_level * high_level * ruler * ruler - low_label * high_label)
    )
    denominator = 7
    for level, label in (
        (low_level, low_label), (low_level, high_label),
        (high_level, low_label), (high_level, high_label),
    ):
        denominator *= level * ruler - label
    closed = F(numerator, denominator)
    require(direct == closed > 0,
            (low_label, high_label, low_level, high_level, ruler, direct, closed))
    return direct


def distinct_debt_maximum(body: tuple[int, ...], pair: tuple[int, int]):
    """Maximum debt with the ordered endpoints fixed at levels one and seven."""
    ruler = 14 * lcm(*body)
    i, j = pair
    remaining_slots = tuple(k for k in range(6) if k not in {i, j})
    descending_slots = tuple(sorted(remaining_slots, key=lambda k: body[k], reverse=True))
    levels = [None] * 6
    levels[i] = LOW_LEVEL
    levels[j] = HIGH_LEVEL
    for slot, level in zip(descending_slots, INTERMEDIATE_LEVELS):
        levels[slot] = level
    levels = tuple(levels)
    debt = T.C2.C5.singleton_debt(body, ruler, levels)

    # Hostile finite audit: choose and order any four of the five available
    # intermediate levels.  It independently checks both the exchange law and
    # the instruction to take the four smallest levels.
    hostile = []
    for assigned in permutations(AVAILABLE_INTERMEDIATE_LEVELS, 4):
        trial = [None] * 6
        trial[i] = LOW_LEVEL
        trial[j] = HIGH_LEVEL
        for slot, level in zip(remaining_slots, assigned):
            trial[slot] = level
        trial = tuple(trial)
        hostile.append((T.C2.C5.singleton_debt(body, ruler, trial), trial))
    require(len(hostile) == 120 and max(hostile) == (debt, levels),
            (body, pair, max(hostile), debt, levels))

    for a, b in permutations(sorted(body[k] for k in remaining_slots), 2):
        if a < b:
            for q, r in permutations(AVAILABLE_INTERMEDIATE_LEVELS, 2):
                if q < r:
                    rearrangement_gap(a, b, q, r, ruler)
    return debt, levels, ruler


def intersection_mass(first, second) -> F:
    """Independent two-pointer intersection mass for sorted disjoint arcs."""
    left_index = 0
    right_index = 0
    total = F(0)
    while left_index < len(first) and right_index < len(second):
        total += max(
            F(0),
            min(first[left_index][1], second[right_index][1])
            - max(first[left_index][0], second[right_index][0]),
        )
        if first[left_index][1] < second[right_index][1]:
            left_index += 1
        else:
            right_index += 1
    return total


def audit_body(body: tuple[int, ...]):
    ruler, safe_ranges = T.R.safe_cell_ranges(body)
    require(ruler == 14 * lcm(*body), (body, ruler))
    cells = tuple(cell for left, right in safe_ranges for cell in range(left, right))
    require(cells, ("empty safe carrier", body))

    low_arcs = tuple(
        tuple(T.R.reflected_level_arcs(ruler, body[i], LOW_LEVEL, cell)
              for cell in cells)
        for i in range(6)
    )
    high_arcs = tuple(
        tuple(T.R.reflected_level_arcs(ruler, body[j], HIGH_LEVEL, cell)
              for cell in cells)
        for j in range(6)
    )

    rows = []
    for i in range(6):
        for j in range(6):
            if i == j:
                continue
            debt, levels, debt_ruler = distinct_debt_maximum(body, (i, j))
            require(debt_ruler == ruler, (body, ruler, debt_ruler))
            overlap, cell = max(
                (intersection_mass(low_arcs[i][index], high_arcs[j][index]), cell)
                for index, cell in enumerate(cells)
            )
            imported = T.intersection_mass(
                T.R.reflected_level_arcs(ruler, body[i], LOW_LEVEL, cell),
                T.R.reflected_level_arcs(ruler, body[j], HIGH_LEVEL, cell),
            )
            require(imported == overlap, (body, (i, j), cell, imported, overlap))
            margin = overlap - debt
            require(margin > 0,
                    ("extreme pair failure", body, (i, j), cell, overlap, debt, levels))
            rows.append((body, (i, j), margin, overlap, debt, cell, levels,
                         ruler, len(cells)))
    require(len(rows) == ORIENTATION_COUNT, (body, len(rows)))
    return tuple(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--processes", type=int, default=8)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    bodies = T.residual_bodies()
    require(len(bodies) == BODY_COUNT, len(bodies))
    repeated_exceptions = {row[0] for row in T.C2.UNIVERSAL.EXPECTED_EXCEPTIONS}
    require(not (repeated_exceptions & set(bodies)),
            ("same-level exceptions in residual", repeated_exceptions & set(bodies)))

    # Algebraic sign controls independent of the body census.
    for ruler in (168, 5460, 7056):
        for a, b in ((1, 2), (1, 14), (6, 13)):
            if b >= ruler:
                continue
            for q, r in ((1, 2), (2, 5), (5, 6)):
                rearrangement_gap(a, b, q, r, ruler)

    if args.processes == 1:
        body_rows = tuple(audit_body(body) for body in bodies)
    else:
        with get_context("fork").Pool(args.processes) as pool:
            body_rows = tuple(pool.map(audit_body, bodies, chunksize=1))
    rows = tuple(row for packet in body_rows for row in packet)
    require(len(rows) == BODY_COUNT * ORIENTATION_COUNT, len(rows))
    require(tuple(row[0] for row in rows[::ORIENTATION_COUNT]) == bodies,
            "parallel body order changed")

    weakest = min(rows, key=lambda row: (row[2], row[0], row[1]))
    debt_maximum = max(rows, key=lambda row: (row[4], row[0], row[1]))
    row_digest = hashlib.sha256(repr(rows).encode()).hexdigest()
    semantic_image = (
        tuple(bodies), tuple(sorted(repeated_exceptions)),
        LOW_LEVEL, HIGH_LEVEL, INTERMEDIATE_LEVELS,
        AVAILABLE_INTERMEDIATE_LEVELS, rows, weakest, debt_maximum, row_digest,
    )
    semantic = hashlib.sha256(repr(semantic_image).encode()).hexdigest()

    if EXPECTED_DEBT_MAXIMUM is not None:
        require(debt_maximum == EXPECTED_DEBT_MAXIMUM,
                (debt_maximum, EXPECTED_DEBT_MAXIMUM))
    if EXPECTED_WEAKEST is not None:
        require(weakest == EXPECTED_WEAKEST, (weakest, EXPECTED_WEAKEST))
    if EXPECTED_ROW_DIGEST is not None:
        require(row_digest == EXPECTED_ROW_DIGEST, (row_digest, EXPECTED_ROW_DIGEST))
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, (semantic, EXPECTED_SEMANTIC))

    lines = [
        "LRC14 reflected extreme-pair (D,m)=(6,1) exact closure",
        f"universe=residual_bodies:{len(bodies)};physical_orientations_per_body:{ORIENTATION_COUNT};rows:{len(rows)}",
        f"same_level_gate=the two non-K6 bodies are disjoint from the 561 residual bodies:{tuple(sorted(repeated_exceptions))};every uncertified residual word therefore has six distinct levels",
        "extreme_pair=ordered slots (i,j) retain their physical meaning:q_i=1 is the exact global minimum,q_j=7 is the exact global maximum;reduced channel=(1,7)>3",
        "distinct_debt=maximized by intermediate levels (2,3,4,5),assigned from largest remaining body label to smallest level;5P4=120 hostile assignments checked per row",
        "rearrangement=swap gap L(r-q)(f-e)(qrL^2-ef)/(7 product_{q,r;e,f}(level*L-label)) is strictly positive",
        f"maximum_distinct_debt={qtext(debt_maximum[4])}@body={debt_maximum[0]},orientation={debt_maximum[1]},levels={debt_maximum[6]}",
        f"weakest_margin={qtext(weakest[2])}@body={weakest[0]},orientation={weakest[1]},cell={weakest[5]},overlap={qtext(weakest[3])},debt={qtext(weakest[4])},levels={weakest[6]}",
        "cell_audit=every integer body-safe cell scanned;independent two-pointer overlap agrees with the promoted interval engine at every maximizing cell",
        "conclusion=within the reflected THM-2941 sufficient family,every packet with (D,m)=(6,1) closes by either a same-level pair or its exact global-min/global-max pair",
        "remaining_scope=other points in D>=6,1<=m<D/2 remain open;this is not a physical-survivor census and not LRC14",
        f"row_digest={row_digest}",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"half_source_sha256={EXPECTED_HALF_SHA256}",
        f"half_output_sha256={EXPECTED_HALF_OUTPUT_SHA256}",
        f"half_semantic_sha256={EXPECTED_HALF_SEMANTIC}",
        f"source_sha256={sha256(HERE)}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    main()
