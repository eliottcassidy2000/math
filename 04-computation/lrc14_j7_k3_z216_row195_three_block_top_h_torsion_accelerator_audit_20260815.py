#!/usr/bin/env python3
"""Self-contained exact row-195 three-block top-h/torsion accelerator.

This self-contained verifier imports neither the active complete-wall verifier nor
any scratch scalar packet.  It reconstructs row 195, its full exact screen,
and its bounded hostile scalar terminal from pinned canonical sources.  It
then proves a disjoint accelerated terminal:

* a three-block weighted common-source bound closes every suffix completion
  on 1,408 of the 1,565 residual denominator passports;
* state-local one-high enumeration is invoked only on the 157 top-h-hostile
  passports, producing 599 cases, all with literal located-torsion
  certificates; and
* the pinned self-contained row-195 audit closes all 137 two-high cases and
  proves no three-high mask; stable filtering leaves its unique denominator-
  two measure case on the one hostile passport with no one-high case.

The historical 314,615 full one-high count and scalar digest are provenance
pins only: this accelerator deliberately never rebuilds that tuple.  All sets
are typed explicitly as states/passports or scalar cases.  This is a row-local
projected-k3 audit, not a 110-row wall run or an LRC(14) theorem.
"""

from __future__ import annotations

import ast
import hashlib
import importlib.util
import re
from array import array
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PREFIX_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z216_sixteen_family_prefix_translated_two_high_closure_scout_20260812.py"
AUDIT_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.py"
NATURAL_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.py"
SOURCE3139 = ROOT / "04-computation/lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.py"
SOURCE3114 = ROOT / "04-computation/lrc14_j7_k3_z227_screen_z226_terminal_double_layer_descent_thm3114.py"
SOURCE3113 = ROOT / "04-computation/lrc14_j7_k3_z229_terminal_z228_screen_double_layer_descent_thm3113.py"
SOURCE3111 = ROOT / "04-computation/lrc14_j7_k3_z230_exact_screen_compressed_complete_cell_descent_thm3111.py"
SOURCE3109 = ROOT / "04-computation/lrc14_j7_k3_z231_exact_screen_complete_cell_cardinality_descent_thm3109.py"
SOURCE3106 = ROOT / "04-computation/lrc14_j7_k3_z232_exact_screen_complete_cell_cardinality_descent_thm3106.py"
TORSION_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z297_ray_status_torsion_closure_thm2941.py"
RAY_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
ATLAS_SOURCE = ROOT / "04-computation/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.py"
ATLAS_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
THM3391_SOURCE = ROOT / "04-computation/lrc14_weighted_common_source_cyclic_support_capacity_thm3391.py"
THM3391_OUTPUT = ROOT / "05-knowledge/results/lrc14_weighted_common_source_cyclic_support_capacity_thm3391.out"
ROW195_SOURCE = ROOT / "04-computation/lrc14_j7_k3_z216_row195_two_high_top_weight_majorant_independent_audit_20260814.py"
ROW195_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z216_row195_two_high_top_weight_majorant_independent_audit_20260814.out"

EXPECTED_HASHES = {
    PREFIX_SOURCE: "cfb020bfc6636a52f1eaf55f82a925e70c11c90da7f87f36b0bd77ece1ec6a62",
    AUDIT_SOURCE: "d3581bc937597270558a77c0a02398a4cfd368306f5a302be01240cdb4cef7c3",
    NATURAL_SOURCE: "430dee7ba03e0d5c9ae0df72ac512500de4f7056cb4663d1c8468bfb93a49bfe",
    SOURCE3139: "92e1c22088998d37db89020ac1ebbb7d1ee17e3539152ab205f8b0cd92532e36",
    SOURCE3114: "8de1e3d03b5070a84b040ac13a173a598646107f85e7ba0defc2ca070808f162",
    SOURCE3113: "1e23ec19fa147c55fb6d38a965eedae0132f5e069b9f820bfd5c300dce4d8f89",
    SOURCE3111: "42323171481deba2371eed9947b2079976cb367dac340cf58b8f1f0c0afb5082",
    SOURCE3109: "1f74f2b2368c04f514f2c388b54c70a9ee66c9387fbc437093884b807b3eb23c",
    SOURCE3106: "f6f64ab8d8ea9b04a1a03e26fc6026efc864e44518e9cb40df4fe8471a4a7991",
    TORSION_SOURCE: "d062c7ac8ebf6a433c8fb1543293e941c85625e2eb40b82fcf05fc2404539b0a",
    RAY_SOURCE: "2ef5e0639354c38b13e17e41f91acb4143c7f60973295b0e2dd0f57eb8f38db2",
    ATLAS_SOURCE: "2af6d96882f336a409a8657070ed76a75c09a53b3789101b83103b051e864ded",
    ATLAS_OUTPUT: "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda",
    THM3391_SOURCE: "22c2ea187e3d43ca55dd61611a0f6d8a70cf7b1111b1f01cb7338bc1aef7e195",
    THM3391_OUTPUT: "9cc8b652ae3552f970fae1b8f46f3b6c1d4316a5170d2f9a37eb5e59495e3062",
    ROW195_SOURCE: "fccc10392624bbdfc2083993ad51a423e8974c135b9bc635351304a71cb0de74",
    ROW195_OUTPUT: "55ad1da385d35f5b38fdc7de2d9916f54ee874b806496d4a7b1270cb526ad30c",
}

INDEX = 195
BODY = (1, 5, 8, 9, 13, 14)
LEVEL = 216
RULER = 458640
HIGH_FLOOR = 45170
EXPECTED_ATLAS = (480, 447, 33)
EXPECTED_ROWS_SHA256 = "53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649"
EXPECTED_ROW = (BODY, RULER, HIGH_FLOOR, True)
EXPECTED_ROW_SHA256 = "904340a0999946e89b1c2c4a7534d9bb73075304ad9d0b30678d293446d7a938"
EXPECTED_SCREEN = (65728, 27198, 36965, 1565)
EXPECTED_FARKAS = (0, 36965)
EXPECTED_RESIDUAL_SHA256 = "351b44f3c270c60d59905767ae19bf798739d540326e6ba4f0521eb6819bbee1"
HISTORICAL_SCALAR_ROW_SHA256 = "a5cdcedc61c07714cacb2164e8793c655fe57ee80fd17ca82149b6710e254036"
HISTORICAL_FULL_ONE_HIGH_COUNT = 314615
EXPECTED_GAP = Fraction(-17651221657, 66699959142726)
EXPECTED_NEEDED_COUNT = 78
EXPECTED_NEEDED_SHA256 = "6b12fcbf3b4d7415cc4006d72ed57099bd635c09a084dda231e9643dad2c7803"
EXPECTED_FIRST_SOURCE = (119368, 37800, 420839)
EXPECTED_FIRST_SOURCE_SHA256 = "745925d13ce82fd6f85d5184de86a12a3cb2c370890d6f60ec3d95d354d13cae"
EXPECTED_TOP_H_PARTITION = (82160, 80687, 1473)
EXPECTED_TOP_H_SHA256 = "ce68e2ad70faf1bd0d5581ac490d84e345eea97479a42969927431c15bcc5985"
EXPECTED_TOP_H_TRIPLES_SHA256 = "9ab6d15791f021aaef022dfb27e0e231cc1826b2bb9ec1ac65747ea6fb6eab82"
EXPECTED_FAILED_TRIPLES_SHA256 = "716c0d4435ad8e7ff3b7912acf748f17591eb2c109bdea3343870aea8cb5c320"
EXPECTED_FAILED_STATE_PARTITION = (1473, 695, 295, 243, 157)
EXPECTED_ATTAINED_FAILED_SHA256 = "82edd53458fef7fafbec9e10683fa8d6799827cf6e4b49138a3b142609fae6a5"
EXPECTED_CRUDE_FAILED_SHA256 = "66bd5ec297b30d419efdb585e0df215d0de00334db89b5e539af54b7d1d89b83"
EXPECTED_STATUS_FAILED_SHA256 = "3b96905d9f42302b530dc8978de1390697de086f709cf10086c613b345c5183f"
EXPECTED_RESIDUAL_FAILED_SHA256 = "e997f387c89563b9d87fca78e049f4484f0712a21879a7f4b6088461f704e947"
EXPECTED_RESIDUAL_PARTITION = (1565, 1408, 157)
EXPECTED_EMPTY_ONE_HIGH_PASSPORTS = ((2, 2, 3920, 6370),)
EXPECTED_REDUCED_ONE_HIGH_CASES = 599
EXPECTED_REDUCED_ONE_HIGH_PASSPORTS = 156
EXPECTED_REDUCED_CASE_STATE_SHA256 = "f5d7e673e8394bf35147c3ce8f08e3198fefdc8da9f1295c78af15828a0f7cbe"
EXPECTED_REDUCED_CASES_SHA256 = "935127e640015f44655de13600eb68fe282d7c77fb7c398341d2ef5da4dcbe7e"
HISTORICAL_FULL_TWO_HIGH_COUNT = 137
EXPECTED_TWO_HIGH_SHA256 = "5811cb585e35831a9c58b2d0af35ab6db89e4503229e6420259d6ce6d4586c3c"
EXPECTED_TWO_HIGH_MECHANISMS = (
    ("common-modulus-support", 76),
    ("denominator-two-measure", 1),
    ("top-weight-majorant", 60),
)
EXPECTED_TWO_HIGH_CLOSURES_SHA256 = "b31f174f9a3c2877142b6b0120d1687e04e7924722389e8b4ebcef8c9b1457da"
EXPECTED_TORSION_SUMMARY = (
    599,
    156,
    179,
    599,
    ((2, 467), (3, 89), (4, 7), (6, 36)),
    ((2, 469), (3, 124), (4, 6)),
    (
        1,
        (2, 5733, 6370, 152880),
        2,
        (237, 320),
        (2, 2, 2, 1, 37800, 37801, 0, 1, 1, 1, 86268),
    ),
    352552,
    57330,
)
EXPECTED_TORSION_CLOSURES_SHA256 = "a67160cf481123639297f2e121776f6c17e552dbf04dee711bc9eb402a0201ac"
EXPECTED_CLOSEST_FAILED = ((10192, 65520, 458640), 0, (22696, 31152, 65520))
EXPECTED_WORST_FAILED = ((458640, 458640, 458640), -77192, (65520, 65520, 65520))
EXPECTED_FCCC_SEMANTIC = "7a70daaa1f25fc6fb9bfd7469ce3e5a18ba650296119789c81bc70fd32016b9c"

# Filled after the first self-contained pass, then required on normal/-O.
EXPECTED_STATES_SHA256 = "2a8584886131c177149055d208f5ccd3b356cc42dcdc84b6533773f6b800b935"
EXPECTED_SUFFIX_BIJECTION_SHA256 = "7d8644411a2d78ed40a6c593a7f977979803aac2e0f2427f788ee58da7c9489c"
EXPECTED_FIRST_MULTIPLICITIES = ((1, 65014), (2, 711), (3, 3))
EXPECTED_SUCCESS_RESIDUAL_SHA256 = "3ddbcf69c2c1521e1440973cdc46be0e85ca769a6ab68a247c8d00115340e372"
EXPECTED_MIN_SUCCESS = (
    4,
    (1872, 6370, 57330, 458640),
    (1872, 57330, 458640),
    (18182, 35662, 65520),
)
EXPECTED_REDUCED_ZERO_HIGH = 109
EXPECTED_REDUCED_ZERO_HIGH_SHA256 = "a1111e7538da5d8bebbc938d28db591396fb30f4504c20fe2b4a3023ccbad092"
EXPECTED_ONE_HIGH_AST_SHA256 = "001859942059725805de93d67634758a40734d8c7e55f3507ab516ddef8b94fa"
EXPECTED_REDUCED_TWO_HIGH_SHA256 = "0153f587817922eab7f9012644fe168b000233b3bded211229d3abb934eed194"
EXPECTED_TWO_HIGH_TAIL = (
    (
        (2, 2, 3920, 6370),
        (2, 2),
        351,
        Fraction(6457, 31302180),
    ),
)
EXPECTED_PARTITION_SHA256 = "a353502ccd9bcee61c6981948bbbba9af8a353d17dfc60e25f41c95efdbb6f38"
EXPECTED_SEMANTIC_SHA256 = "5ad577c0a388c0021cda9506a8147a25732f36f801e7bcf746b856a5d9ff1070"


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


def function_ast_hash(path, function_name):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    matches = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == function_name
    ]
    require(len(matches) == 1, (path, function_name, len(matches)))
    node = matches[0]
    return hashlib.sha256(
        ast.dump(node, annotate_fields=True, include_attributes=False).encode("ascii")
    ).hexdigest(), node


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


def direct_cell_clean(cell, label, ruler):
    residue = (label * cell) % ruler
    return 14 * residue >= ruler and 14 * (residue + label) <= 13 * ruler


def compact_fixed_safe_cells(stream, labels):
    require(stream.L < 2**32, stream.L)
    cells = array("I")
    require(cells.itemsize == 4, cells.itemsize)
    for left, right in stream.ranges:
        for cell in range(left, right):
            if all(
                direct_cell_clean(cell, label, stream.L)
                for label in (stream.first, *labels)
            ):
                cells.append(cell)
    return cells


def top_h_majorant(cells, denominator):
    weights = array("I", (0 for _ in range(denominator)))
    require(weights.itemsize == 4, weights.itemsize)
    for cell in cells:
        weights[cell % denominator] += 1
    require(sum(weights) == len(cells), (denominator, sum(weights), len(cells)))
    histogram = Counter(weights)
    width = (denominator + 6) // 7
    remaining = width
    upper = 0
    for weight in sorted(histogram, reverse=True):
        take = min(remaining, histogram[weight])
        upper += take * weight
        remaining -= take
        if remaining == 0:
            break
    require(remaining == 0, (denominator, width, remaining))
    return (
        denominator,
        len(cells),
        sum(weight > 0 for weight in weights),
        width,
        upper,
        max(weights),
        digest(tuple(sorted(histogram.items()))),
        weights.itemsize * len(weights),
    )


def bit_is_set(bits, residue):
    return bool(bits[residue >> 3] & (1 << (residue & 7)))


def compact_torsion_pigeonhole(cells, denominator):
    bits = bytearray((denominator + 7) // 8)
    support = 0
    for cell in cells:
        residue = cell % denominator
        byte = residue >> 3
        mask = 1 << (residue & 7)
        if not bits[byte] & mask:
            bits[byte] |= mask
            support += 1
    for qualifying_order in range(2, 8):
        if denominator % qualifying_order:
            continue
        quotient = denominator // qualifying_order
        if support <= quotient:
            continue
        selected = None
        for base in range(quotient):
            present = []
            for multiple in range(qualifying_order):
                residue = base + multiple * quotient
                if bit_is_set(bits, residue):
                    present.append(residue)
                    if len(present) == 2:
                        candidate = (present[0], present[1])
                        if selected is None or candidate[0] < selected[0]:
                            selected = candidate
                        break
        require(selected is not None, (denominator, qualifying_order, support))
        first_residue, second_residue = selected
        first_cell = None
        second_cell = None
        for cell in cells:
            residue = cell % denominator
            if residue == first_residue and first_cell is None:
                first_cell = cell
            if residue == second_residue and second_cell is None:
                second_cell = cell
            if first_cell is not None and second_cell is not None:
                break
        require(first_cell is not None and second_cell is not None, selected)
        shift = (second_residue - first_residue) % denominator
        effective = denominator // gcd(denominator, shift)
        require(2 <= effective <= qualifying_order <= 7, (denominator, selected))
        primitive = tuple(
            min(unit, effective - unit)
            for unit in range(1, effective)
            if gcd(unit, effective) == 1
        )
        require(primitive and 7 * min(primitive) >= effective, (denominator, selected))
        return (
            qualifying_order,
            effective,
            support,
            quotient,
            first_cell,
            second_cell,
            first_residue,
            second_residue,
            shift,
            min(primitive),
            len(cells),
        )
    return (None, None, support, None, None, None, None, None, None, None, len(cells))


def verify_torsion_witness(cells, denominator, witness, ruler, labels):
    qorder, effective, support, quotient, first, second, rfirst, rsecond, shift, phase, ncell = witness
    require(qorder is not None, witness)
    require(ncell == len(cells) and first in cells and second in cells, witness)
    require((first % denominator, second % denominator) == (rfirst, rsecond), witness)
    require((second - first) % denominator == shift, witness)
    require(effective == denominator // gcd(denominator, shift), witness)
    require(2 <= effective <= qorder <= 7, witness)
    primitive = tuple(
        min(unit, effective - unit)
        for unit in range(1, effective)
        if gcd(unit, effective) == 1
    )
    require(primitive and phase == min(primitive) and 7 * phase >= effective, witness)
    require(support > quotient, witness)
    require(
        all(
            direct_cell_clean(cell, label, ruler)
            for cell in (first, second)
            for label in labels
        ),
        (denominator, labels, "selected cell not fixed-safe"),
    )
    return support - quotient


def freeze(label, actual, expected):
    if expected is not None:
        require(actual == expected, (label, actual, expected))


def main():
    for path, expected in EXPECTED_HASHES.items():
        require(lf_sha(path) == expected, (path, lf_sha(path), expected))
    require(
        f"semantic_sha256={EXPECTED_FCCC_SEMANTIC}" in ROW195_OUTPUT.read_text(encoding="utf-8"),
        "row195 semantic pin",
    )
    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(syntax)
        ),
        "float literal",
    )

    rows = atlas_rows()
    row = rows[INDEX]
    require(row == EXPECTED_ROW and digest(row) == EXPECTED_ROW_SHA256, (row, digest(row)))
    prefix = load("row195_accelerator_prefix", PREFIX_SOURCE)
    require(Path(prefix.AUDIT_SOURCE).resolve() == AUDIT_SOURCE, prefix.AUDIT_SOURCE)
    require(Path(prefix.NATURAL_SOURCE).resolve() == NATURAL_SOURCE, prefix.NATURAL_SOURCE)
    require(Path(prefix.TORSION_SOURCE).resolve() == TORSION_SOURCE, prefix.TORSION_SOURCE)
    returned_index, raw_screen, direct, legacy = prefix.screen_worker((INDEX, row))
    require(returned_index == INDEX, returned_index)
    screen = tuple(raw_screen)
    require(
        screen[:6]
        == (LEVEL, BODY, RULER, HIGH_FLOOR, RULER // gcd(LEVEL, RULER), True),
        screen[:6],
    )
    require(screen[16] == screen[11] and direct + legacy == screen[11], "screen parity")
    residual = tuple(screen[13])
    require(screen[12] == len(residual), "residual length")
    require(tuple(screen[position] for position in (9, 10, 11, 12)) == EXPECTED_SCREEN, screen[9:13])
    require((direct, legacy) == EXPECTED_FARKAS, (direct, legacy))
    require(digest(residual) == EXPECTED_RESIDUAL_SHA256, digest(residual))
    print(f"phase=full_screen;screen={EXPECTED_SCREEN};residual_sha256={digest(residual)}", flush=True)

    audit = load("row195_accelerator_audit", AUDIT_SOURCE)
    natural = load("row195_accelerator_natural", NATURAL_SOURCE)
    wrapper_chain = (
        (natural, "SOURCE_3139", SOURCE3139),
        (load("row195_accelerator_3139", SOURCE3139), "SOURCE_3114", SOURCE3114),
        (load("row195_accelerator_3114", SOURCE3114), "SOURCE_3113", SOURCE3113),
        (load("row195_accelerator_3113", SOURCE3113), "SOURCE_3111", SOURCE3111),
        (load("row195_accelerator_3111", SOURCE3111), "SOURCE_3109", SOURCE3109),
        (load("row195_accelerator_3109", SOURCE3109), "SOURCE_3106", SOURCE3106),
    )
    for module, attribute, expected_path in wrapper_chain:
        require(
            Path(getattr(module, attribute)).resolve() == expected_path,
            (module.__file__, attribute, getattr(module, attribute), expected_path),
        )
    base3106 = load("row195_accelerator_3106", SOURCE3106)
    require(Path(base3106.__file__).resolve() == SOURCE3106, base3106.__file__)
    engine = audit.status_engine(natural, "row195_accelerator_engine")
    require(Path(engine.__file__).resolve() == TORSION_SOURCE, engine.__file__)
    require(Path(engine.ray.__file__).resolve() == RAY_SOURCE, engine.ray.__file__)
    engine.FIRST = LEVEL
    engine.ray.FIRST = LEVEL
    stream = engine.ray.Stream(BODY)
    require(
        (stream.L, stream.first, stream.high_floor, stream.first_d)
        == (RULER, LEVEL, HIGH_FLOOR, 6370),
        (stream.L, stream.first, stream.high_floor, stream.first_d),
    )

    trials, states, ray_checks, signs = engine.ray.ray_quotient_states(stream)
    require(len(states) == EXPECTED_SCREEN[0], len(states))
    state_keys = tuple(sorted(states))
    states_hash = digest(state_keys)
    freeze("states", states_hash, EXPECTED_STATES_SHA256)
    suffix_bijection = []
    for divisors in state_keys:
        suffix = tuple(sorted(engine.suffix_slots(divisors, stream.first_d)))
        require(tuple(sorted((stream.first_d, *suffix))) == divisors, (divisors, suffix))
        suffix_bijection.append((divisors, suffix))
    suffix_bijection = tuple(suffix_bijection)
    require(len({suffix for _divisors, suffix in suffix_bijection}) == len(states), "suffix collision")
    suffix_bijection_hash = digest(suffix_bijection)
    freeze("suffix_bijection", suffix_bijection_hash, EXPECTED_SUFFIX_BIJECTION_SHA256)
    first_multiplicities = tuple(sorted(Counter(ds.count(stream.first_d) for ds in state_keys).items()))
    freeze("first_multiplicities", first_multiplicities, EXPECTED_FIRST_MULTIPLICITIES)

    needed = tuple(
        sorted(
            {
                denominator
                for divisors in residual
                for denominator in engine.suffix_slots(divisors, stream.first_d)
            }
        )
    )
    require(len(needed) == EXPECTED_NEEDED_COUNT, len(needed))
    require(digest(needed) == EXPECTED_NEEDED_SHA256, digest(needed))
    # No residual passport repeats the fixed divisor in this row, but the
    # global suffix inverse above is valid even when it occurs repeatedly.
    require(stream.first_d not in needed, (stream.first_d, needed))
    require(all(ds.count(stream.first_d) == 1 for ds in residual), "residual repeated first divisor")

    first_cells = compact_fixed_safe_cells(stream, ())
    first_source = (len(first_cells), first_cells[0], first_cells[-1])
    require(first_source == EXPECTED_FIRST_SOURCE, first_source)
    first_source_hash = digest_u32(first_cells)
    require(first_source_hash == EXPECTED_FIRST_SOURCE_SHA256, first_source_hash)
    top_h = tuple(top_h_majorant(first_cells, denominator) for denominator in needed)
    require(digest(top_h) == EXPECTED_TOP_H_SHA256, digest(top_h))
    upper = {record[0]: record[4] for record in top_h}
    triple_records = []
    failed_triples = []
    for suffix in combinations_with_replacement(needed, 3):
        margin = len(first_cells) - sum(upper[denominator] for denominator in suffix)
        record = (suffix, margin, tuple(upper[denominator] for denominator in suffix))
        triple_records.append(record)
        if margin <= 0:
            failed_triples.append(record)
    triple_records = tuple(triple_records)
    failed_triples = tuple(failed_triples)
    triple_partition = (
        len(triple_records),
        len(triple_records) - len(failed_triples),
        len(failed_triples),
    )
    require(triple_partition == EXPECTED_TOP_H_PARTITION, triple_partition)
    require(digest(triple_records) == EXPECTED_TOP_H_TRIPLES_SHA256, digest(triple_records))
    require(digest(failed_triples) == EXPECTED_FAILED_TRIPLES_SHA256, digest(failed_triples))
    record_by_suffix = {record[0]: record for record in triple_records}
    failed_suffixes = {record[0] for record in failed_triples}
    failed_state_by_suffix = {
        tuple(sorted((stream.first_d, *suffix))): record
        for suffix, record in ((record[0], record) for record in failed_triples)
    }
    require(len(failed_state_by_suffix) == len(failed_triples), "failed state collision")
    print(f"phase=top_h;triple_partition={triple_partition}", flush=True)

    attained_failed = {
        divisors: states[divisors]
        for divisors in sorted(failed_state_by_suffix)
        if divisors in states
    }
    crude, status, screened_failed = engine.exact_common_status_screen(stream, attained_failed)
    crude_keys = set(crude)
    status_keys = set(status)
    residual_failed_keys = set(screened_failed)
    require(not (crude_keys & status_keys), "crude/status overlap")
    require(not (crude_keys & residual_failed_keys), "crude/residual overlap")
    require(not (status_keys & residual_failed_keys), "status/residual overlap")
    require(
        set(attained_failed) == crude_keys | status_keys | residual_failed_keys,
        "attained-failed partition is not exact",
    )
    residual_failed = tuple(
        divisors
        for divisors in residual
        if tuple(sorted(engine.suffix_slots(divisors, stream.first_d))) in failed_suffixes
    )
    require(tuple(screened_failed) == residual_failed, "subset screen differs from full-screen intersection")
    failed_partition = (
        len(failed_triples),
        len(attained_failed),
        len(crude),
        len(status),
        len(residual_failed),
    )
    require(failed_partition == EXPECTED_FAILED_STATE_PARTITION, failed_partition)
    attained_packet = tuple((ds, attained_failed[ds]) for ds in sorted(attained_failed))
    require(digest(attained_packet) == EXPECTED_ATTAINED_FAILED_SHA256, digest(attained_packet))
    require(digest(tuple(sorted(crude.items()))) == EXPECTED_CRUDE_FAILED_SHA256, digest(tuple(sorted(crude.items()))))
    require(digest(tuple(sorted(status.items()))) == EXPECTED_STATUS_FAILED_SHA256, digest(tuple(sorted(status.items()))))
    require(digest(residual_failed) == EXPECTED_RESIDUAL_FAILED_SHA256, digest(residual_failed))
    residual_partition = (len(residual), len(residual) - len(residual_failed), len(residual_failed))
    require(residual_partition == EXPECTED_RESIDUAL_PARTITION, residual_partition)
    success_residual = tuple(
        (
            divisors,
            tuple(sorted(engine.suffix_slots(divisors, stream.first_d))),
            record_by_suffix[tuple(sorted(engine.suffix_slots(divisors, stream.first_d)))][1],
            record_by_suffix[tuple(sorted(engine.suffix_slots(divisors, stream.first_d)))][2],
        )
        for divisors in residual
        if divisors not in set(residual_failed)
    )
    require(len(success_residual) == residual_partition[1], len(success_residual))
    require(all(record[2] > 0 for record in success_residual), "nonstrict success")
    success_residual_hash = digest(success_residual)
    freeze("success_residual", success_residual_hash, EXPECTED_SUCCESS_RESIDUAL_SHA256)
    min_success = min((record[2], record[0], record[1], record[3]) for record in success_residual)
    freeze("min_success", min_success, EXPECTED_MIN_SUCCESS)
    closest_failed = max(failed_triples, key=lambda record: (record[1], record[0]))
    worst_failed = min(failed_triples, key=lambda record: (record[1], record[0]))
    require(closest_failed == EXPECTED_CLOSEST_FAILED, closest_failed)
    require(worst_failed == EXPECTED_WORST_FAILED, worst_failed)
    print(
        f"phase=residual_partition;partition={residual_partition};"
        f"failed_sha256={digest(residual_failed)};min_success={min_success}",
        flush=True,
    )

    low, high, sign_census, recurrence_checks = engine.build_literal_tables(stream, set(needed))
    one_high_ast_hash, one_high_ast = function_ast_hash(TORSION_SOURCE, "one_high_cases")
    freeze("one_high_ast", one_high_ast_hash, EXPECTED_ONE_HIGH_AST_SHA256)
    outer_state_loops = [
        node
        for node in ast.walk(one_high_ast)
        if isinstance(node, ast.For)
        and isinstance(node.target, ast.Name)
        and node.target.id == "ds"
        and isinstance(node.iter, ast.Name)
        and node.iter.id == "residuals"
    ]
    case_adds = [
        node
        for node in ast.walk(one_high_ast)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == "cases"
        and node.func.attr == "add"
    ]
    require(len(outer_state_loops) == 1 and case_adds, "one-high locality shape")
    require(
        all(
            len(call.args) == 1
            and isinstance(call.args[0], ast.Tuple)
            and call.args[0].elts
            and isinstance(call.args[0].elts[0], ast.Name)
            and call.args[0].elts[0].id == "ds"
            for call in case_adds
        ),
        "one-high case not tagged by current state",
    )
    reduced_zero_high = engine.zero_high_scalar_passes(stream, residual_failed, low)
    reduced_one_high = engine.one_high_cases(stream, residual_failed, low, high)
    reduced_two_high, reduced_three_high_masks = prefix.two_high_scalar_cases(
        engine, stream, residual_failed, low, high
    )
    require(not reduced_three_high_masks, reduced_three_high_masks)
    require(len(reduced_one_high) == EXPECTED_REDUCED_ONE_HIGH_CASES, len(reduced_one_high))
    require(digest(reduced_one_high) == EXPECTED_REDUCED_CASES_SHA256, digest(reduced_one_high))
    case_state_counts = tuple(sorted(Counter(case[0] for case in reduced_one_high).items()))
    require(len(case_state_counts) == EXPECTED_REDUCED_ONE_HIGH_PASSPORTS, len(case_state_counts))
    require(digest(case_state_counts) == EXPECTED_REDUCED_CASE_STATE_SHA256, digest(case_state_counts))
    empty_one_high = tuple(ds for ds in residual_failed if ds not in dict(case_state_counts))
    require(empty_one_high == EXPECTED_EMPTY_ONE_HIGH_PASSPORTS, empty_one_high)

    freeze("reduced_zero_high_count", len(reduced_zero_high), EXPECTED_REDUCED_ZERO_HIGH)
    freeze("reduced_zero_high", digest(reduced_zero_high), EXPECTED_REDUCED_ZERO_HIGH_SHA256)
    reduced_two_high_hash = digest(reduced_two_high)
    freeze("reduced_two_high", reduced_two_high_hash, EXPECTED_REDUCED_TWO_HIGH_SHA256)
    freeze("two_high_tail", reduced_two_high, EXPECTED_TWO_HIGH_TAIL)
    require(len(reduced_two_high) == 1 and reduced_two_high[0][0] == empty_one_high[0], reduced_two_high)
    row195_text = ROW195_OUTPUT.read_text(encoding="utf-8")
    require(
        f"scalar_row_sha256={HISTORICAL_SCALAR_ROW_SHA256}" in row195_text
        and f"cases={HISTORICAL_FULL_TWO_HIGH_COUNT};cases_sha256={EXPECTED_TWO_HIGH_SHA256};three_high=0" in row195_text
        and f"mechanisms={EXPECTED_TWO_HIGH_MECHANISMS!r}" in row195_text
        and f"closures_sha256={EXPECTED_TWO_HIGH_CLOSURES_SHA256}" in row195_text,
        "fccc two-high bridge text",
    )
    reduced_two_high_closures = tuple(
        (
            case,
            ("denominator-two-measure", Fraction(3, 7), Fraction(3, 91)),
        )
        for case in reduced_two_high
    )
    require(
        reduced_two_high_closures[0][-1]
        == ("denominator-two-measure", Fraction(3, 7), Fraction(3, 91)),
        reduced_two_high_closures,
    )
    print(
        f"phase=hostile_scalar_tail;zero_high={len(reduced_zero_high)};"
        f"one_high={len(reduced_one_high)};one_high_passports={len(case_state_counts)};"
        f"two_high={len(reduced_two_high)};three_high_masks=0;"
        f"one_high_ast_sha256={one_high_ast_hash}",
        flush=True,
    )

    keys_by_labels = {}
    for _divisors, high_d, low_rows, _excess in reduced_one_high:
        labels = tuple(sorted(label for _denominator, label in low_rows))
        keys_by_labels.setdefault(labels, set()).add(high_d)
    cell_counts = {}
    witness_by_key = {}
    surplus_by_key = {}
    qualifying = Counter()
    effective = Counter()
    peak_cell_bytes = 0
    peak_support_bytes = 0
    for labels, denominators in sorted(keys_by_labels.items()):
        cells = compact_fixed_safe_cells(stream, labels)
        require(cells and cells.itemsize == 4, (labels, len(cells)))
        cell_counts[labels] = len(cells)
        peak_cell_bytes = max(peak_cell_bytes, cells.itemsize * len(cells))
        for denominator in sorted(denominators):
            key = (labels, denominator)
            witness = compact_torsion_pigeonhole(cells, denominator)
            require(witness[0] is not None, (key, witness))
            surplus = verify_torsion_witness(
                cells, denominator, witness, stream.L, (stream.first, *labels)
            )
            witness_by_key[key] = witness
            surplus_by_key[key] = surplus
            qualifying[witness[0]] += 1
            effective[witness[1]] += 1
            peak_support_bytes = max(peak_support_bytes, (denominator + 7) // 8)
    closures = []
    minimum = None
    for divisors, high_d, low_rows, excess in reduced_one_high:
        labels = tuple(sorted(label for _denominator, label in low_rows))
        key = (labels, high_d)
        witness = witness_by_key[key]
        surplus = surplus_by_key[key]
        closures.append(
            (
                divisors,
                high_d,
                low_rows,
                excess,
                labels,
                cell_counts[labels],
                surplus,
                witness,
            )
        )
        candidate = (surplus, divisors, high_d, labels, witness)
        if minimum is None or candidate < minimum:
            minimum = candidate
    closures = tuple(closures)
    require(len(witness_by_key) == len(reduced_one_high), (len(witness_by_key), len(reduced_one_high)))
    torsion_summary = (
        len(reduced_one_high),
        len(case_state_counts),
        len(keys_by_labels),
        len(witness_by_key),
        tuple(sorted(qualifying.items())),
        tuple(sorted(effective.items())),
        minimum,
        peak_cell_bytes,
        peak_support_bytes,
    )
    require(torsion_summary == EXPECTED_TORSION_SUMMARY, torsion_summary)
    require(digest(closures) == EXPECTED_TORSION_CLOSURES_SHA256, digest(closures))

    partition_packet = (
        residual_partition,
        (
            len(residual_failed),
            len(reduced_zero_high),
            len(reduced_one_high),
            len(case_state_counts),
            len(reduced_two_high),
            reduced_three_high_masks,
        ),
        case_state_counts,
        empty_one_high,
        digest(reduced_one_high),
        digest(reduced_zero_high),
        digest(reduced_two_high),
        digest(reduced_two_high_closures),
        one_high_ast_hash,
        (
            HISTORICAL_SCALAR_ROW_SHA256,
            HISTORICAL_FULL_ONE_HIGH_COUNT,
            HISTORICAL_FULL_TWO_HIGH_COUNT,
            EXPECTED_TWO_HIGH_SHA256,
            EXPECTED_TWO_HIGH_MECHANISMS,
            EXPECTED_TWO_HIGH_CLOSURES_SHA256,
            EXPECTED_FCCC_SEMANTIC,
        ),
    )
    partition_hash = digest(partition_packet)
    freeze("partition", partition_hash, EXPECTED_PARTITION_SHA256)
    semantic_packet = (
        tuple((path.name, expected) for path, expected in EXPECTED_HASHES.items()),
        EXPECTED_FCCC_SEMANTIC,
        EXPECTED_ATLAS,
        EXPECTED_ROWS_SHA256,
        row,
        tuple(screen[position] for position in (9, 10, 11, 12)),
        (direct, legacy),
        residual,
        states_hash,
        suffix_bijection_hash,
        first_multiplicities,
        trials,
        ray_checks,
        signs,
        needed,
        first_source,
        first_source_hash,
        top_h,
        triple_records,
        attained_packet,
        tuple(sorted(crude.items())),
        tuple(sorted(status.items())),
        residual_failed,
        success_residual,
        min_success,
        closest_failed,
        worst_failed,
        sign_census,
        recurrence_checks,
        partition_packet,
        reduced_two_high_closures,
        torsion_summary,
        digest(closures),
    )
    semantic_hash = digest(semantic_packet)
    freeze("semantic", semantic_hash, EXPECTED_SEMANTIC_SHA256)

    print("ROW195 SELF-CONTAINED THREE-BLOCK TOP-H/TORSION ACCELERATOR")
    print(f"source_sha256_lf={lf_sha(Path(__file__))}")
    print("dependency=" + ";".join(f"{path.name}:{expected}" for path, expected in EXPECTED_HASHES.items()))
    print(f"atlas={EXPECTED_ATLAS};rows_sha256={EXPECTED_ROWS_SHA256};row={row};row_sha256={digest(row)}")
    print(
        f"screen={EXPECTED_SCREEN};farkas={(direct, legacy)};"
        f"full_residual_sha256={digest(residual)};states_sha256={states_hash}"
    )
    print(
        f"suffix_bijection_sha256={suffix_bijection_hash};"
        f"first_divisor_multiplicities={first_multiplicities};"
        f"residual_repeated_first_divisor=0;needed={len(needed)};needed_sha256={digest(needed)}"
    )
    print(f"C216={first_source};C216_sha256={first_source_hash};top_h_sha256={digest(top_h)}")
    print(
        f"triple_partition={triple_partition};triple_sha256={digest(triple_records)};"
        f"failed_sha256={digest(failed_triples)};failed_state_partition={failed_partition}"
    )
    print(
        f"residual_partition={residual_partition};residual_failed_sha256={digest(residual_failed)};"
        f"success_residual_sha256={success_residual_hash};min_strict_success={min_success}"
    )
    print(f"residual_failed={residual_failed!r}")
    print(
        f"historical_scalar_only=scalar_sha256:{HISTORICAL_SCALAR_ROW_SHA256};"
        f"one_high_cases:{HISTORICAL_FULL_ONE_HIGH_COUNT};"
        f"two_gap:{EXPECTED_GAP};two_high_cases:{HISTORICAL_FULL_TWO_HIGH_COUNT};"
        f"two_high_sha256:{EXPECTED_TWO_HIGH_SHA256};full_one_high_not_rebuilt=true"
    )
    print(
        f"hostile_tail=passports:{len(residual_failed)};zero_high:{len(reduced_zero_high)};"
        f"one_high_cases:{len(reduced_one_high)};one_high_passports:{len(case_state_counts)};"
        f"two_high_cases:{len(reduced_two_high)};three_high_masks:{reduced_three_high_masks};"
        f"partition_sha256={partition_hash}"
    )
    print(
        f"reduced_one_high_passports={len(case_state_counts)};"
        f"case_state_sha256={digest(case_state_counts)};"
        f"reduced_one_high_sha256={digest(reduced_one_high)};"
        f"one_high_ast_sha256={one_high_ast_hash};"
        f"empty_one_high_passports={empty_one_high}"
    )
    print(
        f"reduced_zero_high={len(reduced_zero_high)};"
        f"reduced_zero_high_sha256={digest(reduced_zero_high)};"
        f"reduced_two_high={reduced_two_high!r};"
        f"reduced_two_high_sha256={reduced_two_high_hash}"
    )
    print(
        f"fccc_semantic={EXPECTED_FCCC_SEMANTIC};"
        f"full_two_high_mechanisms={EXPECTED_TWO_HIGH_MECHANISMS};"
        f"full_two_high_closures_sha256={EXPECTED_TWO_HIGH_CLOSURES_SHA256};"
        f"reduced_two_high_closure={reduced_two_high_closures!r}"
    )
    print(f"torsion_summary={torsion_summary!r};torsion_closures_sha256={digest(closures)}")
    print(f"closest_top_h_hostile={closest_failed!r};worst_top_h_hostile={worst_failed!r}")
    print(
        "three_block_lemma=for every suffix-label triple (z2,z3,z4), every local coordinate y, "
        "every primitive unit and every ray height, the pullback of each strict d_i/7 danger band "
        "is an arbitrary translated cyclic window using at most ceil(d_i/7) residue classes; "
        "if sum_i K_C216(d_i)<|C216| then there exists c=c(y) in the actual multiplicity-preserving "
        "source C216 safe from all three suffix blockers, while C216 is complete-safe for z1=216"
    )
    print(
        "locality_logic=the pinned one_high_cases AST initializes one set, loops independently over ds in "
        "residuals, tags every inserted case by that same ds, and sorts only after the union; therefore "
        "calling it on the exact hostile passport subset returns exactly the hostile restriction without "
        "constructing the historical 314615-case tuple"
    )
    print(
        "tail_logic=the strict high gate excludes zero-high; no three-high mask exists; 599 one-high cases "
        "close by directly verified located torsion; the only hostile passport without a one-high case is "
        "(2,2,3920,6370), whose unique two-high case is the fccc denominator-two measure terminal"
    )
    print("scope=row195_projected_k3_terminal_only;no_wall_ledger_or_LRC14_claim")
    print(f"semantic_sha256={semantic_hash}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
