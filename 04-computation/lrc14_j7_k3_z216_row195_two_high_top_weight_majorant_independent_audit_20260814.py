#!/usr/bin/env python3
"""Self-contained bounded-memory audit of the row-195 two-high repair.

The scalar packet is reconstructed from canonical pinned sources, not read
from a scratch packet or from the integrated complete-wall verifier.  The
closure then uses an independent Counter-based top-h fibre majorant.  A
strict d/7 band contains at most ceil(d/7) residue classes, so the sum of the
h largest actual-cell fibre weights bounds every unit and every translate.
"""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import re
from array import array
from collections import Counter
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PREFIX_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z216_sixteen_family_prefix_translated_two_high_closure_scout_20260812.py"
AUDIT_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.py"
NATURAL_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.py"
TORSION_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z297_ray_status_torsion_closure_thm2941.py"
RAY_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
ATLAS_SOURCE = ROOT / "04-computation/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.py"
ATLAS_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
THM3391_SOURCE = ROOT / "04-computation/lrc14_weighted_common_source_cyclic_support_capacity_thm3391.py"
THM3391_OUTPUT = ROOT / "05-knowledge/results/lrc14_weighted_common_source_cyclic_support_capacity_thm3391.out"

EXPECTED_HASHES = {
    PREFIX_SOURCE: "cfb020bfc6636a52f1eaf55f82a925e70c11c90da7f87f36b0bd77ece1ec6a62",
    AUDIT_SOURCE: "d3581bc937597270558a77c0a02398a4cfd368306f5a302be01240cdb4cef7c3",
    NATURAL_SOURCE: "430dee7ba03e0d5c9ae0df72ac512500de4f7056cb4663d1c8468bfb93a49bfe",
    TORSION_SOURCE: "d062c7ac8ebf6a433c8fb1543293e941c85625e2eb40b82fcf05fc2404539b0a",
    RAY_SOURCE: "2ef5e0639354c38b13e17e41f91acb4143c7f60973295b0e2dd0f57eb8f38db2",
    ATLAS_SOURCE: "2af6d96882f336a409a8657070ed76a75c09a53b3789101b83103b051e864ded",
    ATLAS_OUTPUT: "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda",
    THM3391_SOURCE: "22c2ea187e3d43ca55dd61611a0f6d8a70cf7b1111b1f01cb7338bc1aef7e195",
    THM3391_OUTPUT: "9cc8b652ae3552f970fae1b8f46f3b6c1d4316a5170d2f9a37eb5e59495e3062",
}

LEVEL = 216
INDEX = 195
BODY = (1, 5, 8, 9, 13, 14)
RULER = 458640
HIGH_FLOOR = 45170
LOW_LABEL = 351
EXPECTED_ATLAS = (480, 447, 33)
EXPECTED_ROWS_SHA256 = "53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649"
EXPECTED_ROW = (BODY, RULER, HIGH_FLOOR, True)
EXPECTED_ROW_SHA256 = "904340a0999946e89b1c2c4a7534d9bb73075304ad9d0b30678d293446d7a938"
EXPECTED_SCREEN_TOTALS = (65728, 27198, 36965, 1565)
EXPECTED_RESIDUAL_SHA256 = "351b44f3c270c60d59905767ae19bf798739d540326e6ba4f0521eb6819bbee1"
EXPECTED_GAP = Fraction(-17651221657, 66699959142726)
EXPECTED_SCALAR_ROW_SHA256 = "a5cdcedc61c07714cacb2164e8793c655fe57ee80fd17ca82149b6710e254036"
EXPECTED_CASES_SHA256 = "5811cb585e35831a9c58b2d0af35ab6db89e4503229e6420259d6ce6d4586c3c"
EXPECTED_CELL_STATE = (100776, 37987, 420652)
EXPECTED_CELL_SHA256 = "89236b7d22a6be06afecaaaba8f2e0be1346b1cab3c0b720314baa4014d87c36"
EXPECTED_MECHANISMS = (
    ("common-modulus-support", 76),
    ("denominator-two-measure", 1),
    ("top-weight-majorant", 60),
)
EXPECTED_WEAKEST = (10456, 130, (32760, 458640), 24800, 65520)
EXPECTED_MAJORANTS_SHA256 = "46ea72f9a0128ba19d10eba82bef883b2754185482edb22fe2e00b29b5fbbc1e"
EXPECTED_CLOSURES_SHA256 = "b31f174f9a3c2877142b6b0120d1687e04e7924722389e8b4ebcef8c9b1457da"
EXPECTED_PERIOD_HISTOGRAM = (
    (1, 760),
    (2, 7944),
    (3, 9848),
    (4, 11626),
    (5, 1616),
)
EXPECTED_HOSTILE_COSETS = (
    (
        32760,
        10920,
        2340,
        (2340, 13260, 24180),
        (0, 1, 2),
        6912,
        "871030240a15e1f19e9a330e4f6d7dac2faa8e22b03977d5f5db8c5014de19ae",
        1,
    ),
    (
        57330,
        19110,
        7800,
        (26910, 46020, 7800),
        (1, 2, 0),
        12096,
        "9b9cb7346693f73e89d3c6d9e4ec5225602466dfa8be2e5550f0b33ff7f78dbf",
        1,
    ),
)
EXPECTED_SEMANTIC_SHA256 = "7a70daaa1f25fc6fb9bfd7469ce3e5a18ba650296119789c81bc70fd32016b9c"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def digest(value):
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def digest_u32(values):
    state = hashlib.sha256()
    for value in values:
        state.update(int(value).to_bytes(4, "big"))
    return state.hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def atlas_rows():
    pattern = re.compile(
        r"^row=E=([0-9,]+);h=[^;]+;r=([0-9]+);"
        r"L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    rows = []
    for line in ATLAS_OUTPUT.read_text(encoding="utf-8").splitlines():
        if not line.startswith("row="):
            continue
        match = pattern.match(line)
        require(match is not None, line)
        if int(match.group(5)) != LEVEL:
            continue
        rows.append(
            (
                tuple(map(int, match.group(1).split(","))),
                int(match.group(3)),
                int(match.group(4)),
                LEVEL < int(match.group(4)),
            )
        )
    rows = tuple(rows)
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    require(census == EXPECTED_ATLAS, census)
    require(digest(rows) == EXPECTED_ROWS_SHA256, digest(rows))
    return rows


def scalar_packet(row):
    prefix = load("row195_selfcontained_prefix", PREFIX_SOURCE)
    require(Path(prefix.AUDIT_SOURCE).resolve() == AUDIT_SOURCE, prefix.AUDIT_SOURCE)
    require(Path(prefix.NATURAL_SOURCE).resolve() == NATURAL_SOURCE, prefix.NATURAL_SOURCE)
    require(Path(prefix.TORSION_SOURCE).resolve() == TORSION_SOURCE, prefix.TORSION_SOURCE)
    returned_index, raw_screen, direct, legacy = prefix.screen_worker((INDEX, row))
    require(returned_index == INDEX, returned_index)
    screen = tuple(raw_screen)
    require(screen[:6] == (LEVEL, BODY, RULER, HIGH_FLOOR, RULER // gcd(LEVEL, RULER), True), screen[:6])
    require(screen[16] == screen[11], "status replay")
    require(direct + legacy == screen[11], "Farkas partition")
    residual = tuple(screen[13])
    require(screen[12] == len(residual), "residual length")
    require(tuple(screen[position] for position in (9, 10, 11, 12)) == EXPECTED_SCREEN_TOTALS, screen[9:13])
    require(digest(residual) == EXPECTED_RESIDUAL_SHA256, digest(residual))

    audit = load("row195_selfcontained_audit", AUDIT_SOURCE)
    natural = load("row195_selfcontained_natural", NATURAL_SOURCE)
    engine = audit.status_engine(natural, "row195_selfcontained")
    require(Path(engine.__file__).resolve() == TORSION_SOURCE, engine.__file__)
    require(Path(engine.ray.__file__).resolve() == RAY_SOURCE, engine.ray.__file__)
    engine.FIRST = LEVEL
    engine.ray.FIRST = LEVEL
    stream = engine.ray.Stream(BODY)
    require((stream.L, stream.first, stream.high_floor) == (RULER, LEVEL, HIGH_FLOOR), (stream.L, stream.first, stream.high_floor))
    needed = {
        denominator
        for divisors in residual
        for denominator in engine.suffix_slots(divisors, stream.first_d)
    }
    low, high, sign_census, recurrence_checks = engine.build_literal_tables(stream, needed)
    two_gap, two_witness = engine.duplicate_two_high_gap(stream, residual, low, high)
    require(two_gap == EXPECTED_GAP, two_gap)
    zero_high = engine.zero_high_scalar_passes(stream, residual, low)
    one_high = engine.one_high_cases(stream, residual, low, high)
    two_cases, three_high_masks = prefix.two_high_scalar_cases(engine, stream, residual, low, high)
    require(not three_high_masks, three_high_masks)
    denominator_patterns = tuple(sorted(Counter(case[1] for case in two_cases).items()))
    packet = (
        INDEX,
        stream.L,
        screen[9],
        screen[10],
        screen[11],
        len(residual),
        digest(residual),
        two_gap,
        two_witness,
        len(zero_high),
        len(one_high),
        tuple(sorted(needed)),
        len(two_cases),
        denominator_patterns,
        two_cases,
        sign_census,
        recurrence_checks,
    )
    require(digest(packet) == EXPECTED_SCALAR_ROW_SHA256, digest(packet))
    require(len(two_cases) == 137 and digest(two_cases) == EXPECTED_CASES_SHA256, (len(two_cases), digest(two_cases)))
    require(all(case[2] == LOW_LABEL and len(case[1]) == 2 for case in two_cases), "two-high typing")
    return packet, stream, engine, two_cases, (direct, legacy)


def direct_cell_clean(cell, label, ruler):
    residue = (label * cell) % ruler
    return 14 * residue >= ruler and 14 * (residue + label) <= 13 * ruler


def fixed_cells(stream):
    cells = array("I")
    require(cells.itemsize == 4, cells.itemsize)
    for left, right in stream.ranges:
        for cell in range(left, right):
            if direct_cell_clean(cell, LEVEL, stream.L) and direct_cell_clean(cell, LOW_LABEL, stream.L):
                cells.append(cell)
    require((len(cells), cells[0], cells[-1]) == EXPECTED_CELL_STATE, (len(cells), cells[0], cells[-1]))
    require(digest_u32(cells) == EXPECTED_CELL_SHA256, digest_u32(cells))
    return cells


def broken_l_over_nine_orbit(cells, stream):
    period = stream.L // 9
    counts = array("B", (0 for _ in range(period)))
    for cell in cells:
        counts[cell % period] += 1
    histogram = tuple(sorted((weight, count) for weight, count in Counter(counts).items() if weight))
    require(histogram == EXPECTED_PERIOD_HISTOGRAM, histogram)
    source = 37987
    orbit = tuple(source + offset * period for offset in range(9))
    present = set(cells)
    mask = tuple(cell in present for cell in orbit)
    target = orbit[1]
    require(mask == (True, False, False, True, True, False, False, False, False), mask)
    require(direct_cell_clean(target, LEVEL, stream.L) and direct_cell_clean(target, LOW_LABEL, stream.L), target)
    require(not any(left <= target < right for left, right in stream.ranges), target)
    return period, histogram, max(counts), source, target, mask


def hostile_three_cell_coset_control(cells, cases):
    """Exhaust the explicit order-three repair of the gap witness pair."""
    source_cells = (198900, 46020, 122460)
    cell_set = set(cells)
    require(all(cell in cell_set for cell in source_cells), source_cells)
    clean_packet = tuple(
        (
            cell,
            tuple(
                (label, (label * cell) % RULER, direct_cell_clean(cell, label, RULER))
                for label in (LEVEL, LOW_LABEL)
            ),
        )
        for cell in source_cells
    )
    require(all(flag for _cell, labels in clean_packet for _label, _residue, flag in labels), clean_packet)
    matching = tuple(
        (case_index, case)
        for case_index, case in enumerate(cases)
        if case[1] == (32760, 57330)
    )
    require(len(matching) == 1 and matching[0][0] == 124, matching)

    records = []
    for denominator, base, expected_positions in (
        (32760, 2340, (0, 1, 2)),
        (57330, 7800, (1, 2, 0)),
    ):
        step = denominator // 3
        residues = tuple(cell % denominator for cell in source_cells)
        positions = tuple(((residue - base) % denominator) // step for residue in residues)
        require(positions == expected_positions, (denominator, residues, positions))
        require({residue % denominator for residue in residues} == {(base + offset * step) % denominator for offset in range(3)}, residues)
        require(7 * step > denominator, (denominator, step))
        state = hashlib.sha256()
        unit_count = 0
        for unit in range(1, denominator):
            if gcd(unit, denominator) != 1:
                continue
            images = tuple(sorted((unit * residue) % denominator for residue in residues))
            gaps = tuple(
                (images[(position + 1) % 3] - images[position]) % denominator
                for position in range(3)
            )
            require(gaps == (step, step, step), (denominator, unit, images, gaps))
            state.update(repr((unit, images, gaps)).encode("ascii"))
            unit_count += 1
        records.append(
            (
                denominator,
                step,
                base,
                residues,
                positions,
                unit_count,
                state.hexdigest(),
                1,
            )
        )
    records = tuple(records)
    require(records == EXPECTED_HOSTILE_COSETS, records)
    require(1 + 1 < len(source_cells), records)
    return matching[0], source_cells, clean_packet, records, (1, 1, 2, 3, 1)


def common_capacity(cells, denominators, support_cache):
    modulus = lcm(*denominators)
    if modulus not in support_cache:
        support_cache[modulus] = len({cell % modulus for cell in cells})
    support = support_cache[modulus]
    capacities = tuple(((denominator + 6) // 7) * (modulus // denominator) for denominator in denominators)
    return modulus, support, capacities, support - sum(capacities)


def top_weight_majorant(cells, denominator):
    weights = Counter(cell % denominator for cell in cells)
    width = (denominator + 6) // 7
    largest = tuple(sorted(weights.values(), reverse=True)[:width])
    upper = sum(largest)
    require(len(largest) <= width and upper <= len(cells), (denominator, width, upper))
    return (
        denominator,
        len(cells),
        len(weights),
        width,
        upper,
        len(cells) - upper,
        largest[0],
        largest[-1] if largest else 0,
        digest(largest),
    )


def close_cases(cells, cases):
    support_cache = {}
    preliminary = []
    needed_denominators = set()
    for case_index, case in enumerate(cases):
        divisors, high_ds, low_label, excess = case
        common = common_capacity(cells, high_ds, support_cache)
        if common[-1] > 0:
            mechanism = "common-modulus-support"
        elif high_ds == (2, 2):
            mechanism = "denominator-two-measure"
        else:
            mechanism = "top-weight-majorant"
            needed_denominators.update(high_ds)
        preliminary.append((case_index, divisors, high_ds, low_label, excess, common, mechanism))

    majorants = {denominator: top_weight_majorant(cells, denominator) for denominator in sorted(needed_denominators)}
    closures = []
    counts = Counter()
    weakest = None
    for record in preliminary:
        case_index, divisors, high_ds, low_label, excess, common, mechanism = record
        if mechanism == "common-modulus-support":
            certificate = (mechanism, common[-1])
        elif mechanism == "denominator-two-measure":
            certificate = (mechanism, Fraction(3, 7), Fraction(3, 91))
        else:
            first = majorants[high_ds[0]][4]
            second = majorants[high_ds[1]][4]
            lower = len(cells) - first - second
            require(lower > 0, (case_index, high_ds, lower))
            candidate = (lower, case_index, high_ds, first, second)
            weakest = candidate if weakest is None or candidate < weakest else weakest
            certificate = (mechanism, lower, majorants[high_ds[0]], majorants[high_ds[1]])
        counts[mechanism] += 1
        closures.append((case_index, divisors, high_ds, low_label, excess, common, certificate))
    mechanisms = tuple(sorted(counts.items()))
    require(mechanisms == EXPECTED_MECHANISMS, mechanisms)
    require(weakest == EXPECTED_WEAKEST, weakest)
    capacity_records = tuple(majorants[denominator] for denominator in sorted(majorants))
    closures = tuple(closures)
    require(digest(capacity_records) == EXPECTED_MAJORANTS_SHA256, digest(capacity_records))
    require(digest(closures) == EXPECTED_CLOSURES_SHA256, digest(closures))
    return (
        tuple(preliminary),
        capacity_records,
        closures,
        mechanisms,
        weakest,
        tuple(sorted(support_cache.items())),
    )


def main():
    for path, expected in EXPECTED_HASHES.items():
        require(lf_sha(path) == expected, (path, lf_sha(path), expected))
    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)), "assert node")
    require(not any(isinstance(node, ast.Constant) and isinstance(node.value, float) for node in ast.walk(syntax)), "float node")

    rows = atlas_rows()
    row = rows[INDEX]
    require(row == EXPECTED_ROW and digest(row) == EXPECTED_ROW_SHA256, (row, digest(row)))
    scalar, stream, engine, cases, farkas = scalar_packet(row)
    cells = fixed_cells(stream)
    require(engine.fixed_safe_cells(stream, (LOW_LABEL,)) == tuple(cells), "fixed-safe source parity")
    orbit = broken_l_over_nine_orbit(cells, stream)
    hostile_coset = hostile_three_cell_coset_control(cells, cases)
    preliminary, majorants, closures, mechanisms, weakest, supports = close_cases(cells, cases)
    semantic_packet = (
        tuple((path.name, expected) for path, expected in EXPECTED_HASHES.items()),
        EXPECTED_ATLAS,
        digest(rows),
        row,
        scalar,
        farkas,
        tuple(stream.ranges),
        EXPECTED_CELL_STATE,
        digest_u32(cells),
        orbit,
        hostile_coset,
        supports,
        majorants,
        preliminary,
        closures,
        mechanisms,
        weakest,
    )
    semantic = digest(semantic_packet)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    print("ROW195 SELF-CONTAINED TOP-WEIGHT MAJORANT AUDIT")
    print(f"source_sha256_lf={lf_sha(Path(__file__))}")
    print("dependency=" + ";".join(f"{path.name}:{expected}" for path, expected in EXPECTED_HASHES.items()))
    print(f"atlas={EXPECTED_ATLAS};rows_sha256={digest(rows)};row={row};row_sha256={digest(row)}")
    print(f"screen={EXPECTED_SCREEN_TOTALS};residual_sha256={scalar[6]};farkas={farkas}")
    print(f"scalar_row_sha256={digest(scalar)};gap={scalar[7]};cases={len(cases)};cases_sha256={digest(cases)};three_high=0")
    print(f"cells={EXPECTED_CELL_STATE};cell_sha256={digest_u32(cells)};ranges={len(stream.ranges)}")
    print(f"broken_L_over_9_orbit={orbit}")
    print(f"hostile_three_cell_coset={hostile_coset}")
    print(f"mechanisms={mechanisms};unique_majorant_denominators={len(majorants)};majorants_sha256={digest(majorants)};closures_sha256={digest(closures)}")
    print(f"weakest={weakest}")
    print("lemma=every unit-pullback strict d/7 band uses at most ceil(d/7) residue classes;the sum of the largest actual-cell fibre weights uniformly bounds every height and translate")
    print("d2_scope=inherited exact projected-measure terminal;safe_mass_at_least_3/7=39/91>36/91;not a pointwise weighted-cardinality claim")
    print("scope=row195_two-high_only;no_one-high_closure_claim;no_wall_or_ledger_decrement")
    print(f"semantic_sha256={semantic}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
