#!/usr/bin/env python3
"""Exact third intrinsic-ruler prefix closure with modular multicovers.

The two earlier cost-prefixes closed twenty-three projected-k3, z1=216 wall
rows.  This scout reconstructs that state from the independent audit, derives
the next complete-family prefix from the intrinsic cost order, and screens its
four labelled rows.  Every common-status elimination is then recertified by a
deterministic threshold-layer modular multicover; solver-selected Farkas
vectors are checked but are neither normalized nor persisted.

This is a finite exact calculation in a necessary projected atlas.  It does
not restore a physical speed entry, endpoint, owner, phase, current, arbitrary
k<=1, the rung, or LRC(14).
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import itertools
from collections import Counter
from fractions import Fraction
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
AUDIT_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.py"
)
AUDIT_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.out"
)
MULTICOVER_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_threshold_layer_multicover_dual_anatomy_scout_20260803.py"
)
MULTICOVER_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_threshold_layer_multicover_dual_anatomy_scout_20260803.out"
)

AUDIT_SOURCE_SHA256 = (
    "d3581bc937597270558a77c0a02398a4cfd368306f5a302be01240cdb4cef7c3"
)
AUDIT_OUTPUT_SHA256 = (
    "e035531d0af15998a8a07d083042ff2486c80c7fdd729d902922179409d1725a"
)
AUDIT_SEMANTIC_SHA256 = (
    "4bb4eca9949eb4f07b618f9e49bbb5de63492364a29a4a5abfc84a6acf105e11"
)
MULTICOVER_SOURCE_SHA256 = (
    "e7462eabab773133688079141708d2377742e2d6f9a9480089666956f472020c"
)
MULTICOVER_OUTPUT_SHA256 = (
    "f83d7ca1a3e2557330ec7ed655dc756cef21d6abdb52326d3623563c2abc407d"
)
MULTICOVER_SEMANTIC_SHA256 = (
    "bd7ac011d960afe7c4566ac9c17e9df3d8c2cf710e6023a715cb6a6d1ce1408d"
)

LEVEL = 216
EXPECTED_CENSUS = (480, 447, 33, 357, 33)
EXPECTED_SELECTED_FAMILIES = (
    (
        1_729_728,
        1,
        72,
        72_072,
        (157,),
        ((157, (1, 4, 9, 11, 12, 13), 24, 1_729_728),),
    ),
    (
        2_121_504,
        3,
        24,
        25_872,
        (53, 293, 357),
        (
            (53, (1, 2, 6, 8, 11, 14), 26, 672_672),
            (293, (2, 4, 6, 8, 11, 14), 26, 672_672),
            (357, (2, 6, 8, 11, 12, 14), 30, 776_160),
        ),
    ),
)
EXPECTED_NEXT_FAMILY = (2_293_200, 1, 24, 76_440, (141,))
EXPECTED_SELECTED_INDICES = (157, 53, 293, 357)
EXPECTED_BY_ROW = (
    (157, (111, 93, 18, 0), ((1, 17), (2, 1))),
    (53, (27, 17, 10, 0), ((1, 10),)),
    (293, (14, 9, 5, 0), ((1, 5),)),
    (357, (271, 170, 101, 0), ((1, 101),)),
)
EXPECTED_TOTALS = (423, 289, 134, 0)
EXPECTED_FARKAS_COUNTS = (0, 134)
EXPECTED_LAYER_COUNTS = ((1, 133), (2, 1))
EXPECTED_OLD_TO_LAYER = (
    ("coordinate_union", ((1, 128),)),
    ("two_fan", ((1, 3),)),
    ("weighted_core", ((2, 1),)),
    ("zero_reduced_union", ((1, 2),)),
)
EXPECTED_UNIQUE_INSTANCES = 74
EXPECTED_INSTANCE_MULTIPLICITIES = (
    (1, 51), (2, 14), (3, 2), (4, 2), (5, 2), (6, 1), (11, 1), (14, 1),
)
EXPECTED_MULTI_ADDRESS = (157, (99, 616, 1001, 8008))
EXPECTED_MULTI_CERTIFICATE = ((3, 4), (2, 1, 2, 1), 510)
EXPECTED_PROPER_SUPPORTS = 3
EXPECTED_HOSTILE_COUNTS = (115, 50, 57, 8)
EXPECTED_ARITY3 = ((1, 4, 5), (1, 3, 1, 1), 11)

EXPECTED_SELECTED_PACKET_SHA256 = (
    "8a74873435634932fdf284a3fb7fee26c4558140f0b46fd611dda56aff0c47d1"
)
EXPECTED_SCREEN_SHA256 = (
    "17bd632e09975d5a84d4e5f8667023eb9c1b9de21c8ae16f639612ff3d380e74"
)
EXPECTED_STATUS_SHA256 = (
    "3cbc8e6f2888b2852ac65cf643755fdf2bdbeaac7c751d4f55fbdb72c13732a6"
)
EXPECTED_MINIMALITY_SHA256 = (
    "63b359ee5534803576136a31418d2c2152bf60e0202fb4fef07c154a6382da90"
)
EXPECTED_HOSTILE_SHA256 = (
    "afc67bd8bd482be01bf2a3b1f9e83e44b6fa244d54b42fd856eba42ef3e4b823"
)
EXPECTED_SEMANTIC_SHA256 = (
    "4d8547f53b287d82bb725cd6b3ed78c390ba03a37a3dcb3ab76ac98389e045a3"
)

LEDGER_BEFORE = 373_161
LEDGER_AFTER = 373_157
WALL_BEFORE = 357
WALL_AFTER = 353
FAMILIES_BEFORE = 33
FAMILIES_AFTER = 31


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, (path, "bare CR"))
    return hashlib.sha256(payload).hexdigest()


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sum_counts(records):
    return tuple(
        sum(record[position] for record in records) for position in range(4)
    )


def family_head(ranked, count):
    return tuple(
        (family[0], family[1], family[2], family[3], family[4])
        for family in ranked[:count]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    dependencies = (
        (AUDIT_SOURCE, AUDIT_SOURCE_SHA256, None),
        (AUDIT_OUTPUT, AUDIT_OUTPUT_SHA256, AUDIT_SEMANTIC_SHA256),
        (MULTICOVER_SOURCE, MULTICOVER_SOURCE_SHA256, None),
        (
            MULTICOVER_OUTPUT,
            MULTICOVER_OUTPUT_SHA256,
            MULTICOVER_SEMANTIC_SHA256,
        ),
    )
    for path, digest, semantic in dependencies:
        require(lf_sha(path) == digest, (path, "dependency changed"))
        if semantic is not None:
            require(
                f"semantic_sha256={semantic}"
                in path.read_text(encoding="utf-8"),
                (path, "semantic changed"),
            )

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require((assert_nodes, float_nodes) == (0, 0), "truth-gate syntax")

    audit = load("third_prefix_audit", AUDIT_SOURCE)
    multicover = load("third_prefix_multicover", MULTICOVER_SOURCE)
    natural = audit.load("third_prefix_natural", audit.NATURAL_SOURCE)
    order = audit.load("third_prefix_order", audit.ORDER_SOURCE)
    rows, components = natural.atlas_rows()
    order_rows, order_components = order.atlas_rows()
    require((rows, components) == (order_rows, order_components), "atlas parser")

    gcd8 = {
        index for index, row in enumerate(rows) if gcd(LEVEL, row[1]) == 8
    }
    gcd18 = {
        index for index, row in enumerate(rows) if gcd(LEVEL, row[1]) == 18
    }
    order_indices = {
        index for index, row in enumerate(rows) if not row[3]
    }
    natural_indices = {
        index for index, row in enumerate(rows)
        if (gcd(LEVEL, row[1]), row[1]) in audit.NATURAL_FAMILY_KEYS
    }
    prior = gcd8 | gcd18 | order_indices | natural_indices
    require(len(prior) == 100, len(prior))

    first_live = tuple(
        index for index, row in enumerate(rows)
        if row[3] and index not in prior
    )
    first_ranked = audit.ranked_families(rows, components, first_live)
    first_families = audit.through_next_nonsingleton(first_ranked)
    first_indices = tuple(
        index for family in first_families for index in family[4]
    )

    second_live = tuple(
        index for index, row in enumerate(rows)
        if row[3] and index not in prior | set(first_indices)
    )
    second_ranked = audit.ranked_families(rows, components, second_live)
    second_families = audit.through_next_nonsingleton(second_ranked)
    second_indices = tuple(
        index for family in second_families for index in family[4]
    )

    closed = prior | set(first_indices) | set(second_indices)
    live = tuple(
        index for index, row in enumerate(rows)
        if row[3] and index not in closed
    )
    ranked = audit.ranked_families(rows, components, live)
    selected_families = audit.through_next_nonsingleton(ranked)
    selected_indices = tuple(
        index for family in selected_families for index in family[4]
    )
    census = (
        len(rows),
        sum(row[3] for row in rows),
        sum(not row[3] for row in rows),
        len(live),
        len(ranked),
    )
    require(census == EXPECTED_CENSUS, census)
    require((len(first_indices), len(second_indices)) == (5, 18), "prior prefixes")
    require(tuple(selected_families) == EXPECTED_SELECTED_FAMILIES,
            selected_families)
    require(selected_indices == EXPECTED_SELECTED_INDICES, selected_indices)
    require(
        family_head(ranked[len(selected_families):], 1)[0]
        == EXPECTED_NEXT_FAMILY,
        ranked[:3],
    )
    require(
        selected_families[-1][0] < ranked[len(selected_families)][0],
        "cost tie crosses selected boundary",
    )

    selected_packet = tuple(
        (
            index,
            rows[index][0],
            gcd(LEVEL, rows[index][1]),
            rows[index][1],
            components[index],
            rows[index][1] * components[index],
            rows[index][3],
        )
        for index in selected_indices
    )
    selected_packet_sha = hashlib.sha256(
        repr(selected_packet).encode()
    ).hexdigest()
    if EXPECTED_SELECTED_PACKET_SHA256 is not None:
        require(selected_packet_sha == EXPECTED_SELECTED_PACKET_SHA256,
                selected_packet_sha)

    canonical = []
    counts = {}
    direct_total = 0
    legacy_total = 0
    for index in selected_indices:
        screened = audit.screen_worker((index, (LEVEL, *rows[index])))
        returned_index, screen_row, direct, legacy = screened
        require(returned_index == index, (index, returned_index))
        source = rows[index]
        divisor_gcd = gcd(LEVEL, source[1])
        require(
            screen_row[:6]
            == (
                LEVEL,
                source[0],
                source[1],
                source[2],
                source[1] // divisor_gcd,
                source[3],
            ),
            (index, screen_row[:6]),
        )
        require(screen_row[16] == screen_row[11], (index, "status verification"))
        require(direct + legacy == screen_row[11], (index, "Farkas count"))
        counts[index] = tuple(screen_row[position] for position in (9, 10, 11, 12))
        canonical.append((index, tuple(screen_row[:19])))
        direct_total += direct
        legacy_total += legacy
    canonical = tuple(canonical)
    screen_sha = hashlib.sha256(repr(canonical).encode()).hexdigest()
    if EXPECTED_SCREEN_SHA256 is not None:
        require(screen_sha == EXPECTED_SCREEN_SHA256, screen_sha)
    require(
        sum_counts(tuple(counts[index] for index in selected_indices))
        == EXPECTED_TOTALS,
        counts,
    )
    require((direct_total, legacy_total) == EXPECTED_FARKAS_COUNTS,
            (direct_total, legacy_total))
    residual_passports = tuple(
        (index, row[13]) for index, row in canonical if row[12]
    )
    require(not residual_passports, residual_passports)

    certificate_records = []
    layer_by_row = {}
    old_to_layer = {}
    old_to_layer_counter = {}
    for index in selected_indices:
        eng = audit.status_engine(natural, f"selected_{index}")
        eng.FIRST = LEVEL
        eng.ray.FIRST = LEVEL
        stream = eng.ray.Stream(rows[index][0])
        require(
            (stream.L, stream.high_floor) == (rows[index][1], rows[index][2]),
            (index, "atlas"),
        )
        _trials, states, _checks, _signs = eng.ray.ray_quotient_states(stream)
        crude, status, residual = eng.exact_common_status_screen(stream, states)
        require(
            (len(states), len(crude), len(status), len(residual)) == counts[index],
            (index, "screen reconstruction"),
        )
        row_layers = Counter()
        for divisors, witness in sorted(status.items()):
            q, cofactor, marginals, capacity_set, histogram, raw_dual = witness
            divisor_lcm = lcm(*divisors)
            require(divisor_lcm == q * cofactor, (index, divisors, "cofactor"))
            rebuilt_marginals, capacities = eng.ray.local.hunter_status_data(
                divisor_lcm, divisors, q
            )
            require(rebuilt_marginals == marginals, (index, divisors, "marginals"))
            require(
                tuple(sorted(set(capacities))) == capacity_set,
                (index, divisors, "capacities"),
            )
            eng.independent_farkas_check(
                q, marginals, capacities, histogram, raw_dual
            )
            support_size, certificates = multicover.minimal_binary_multicovers(
                capacities, marginals, histogram
            )
            require(support_size is not None and certificates,
                    (index, divisors, "no modular multicover"))
            selected, weights, gap = certificates[0]
            require(len(selected) == support_size, (index, divisors, selected))
            require(
                multicover.coordinate_minimal(capacities, selected, weights),
                (index, divisors, "coordinate-nonminimal"),
            )
            primitive_values = tuple(
                value for value in (*weights, *([1] * support_size)) if value
            )
            require(gcd(*primitive_values) == 1, (index, divisors, "nonprimitive"))
            tight = multicover.tight_patterns(capacities, selected, weights)
            submodular, supermodular, sub_witness, super_witness = (
                multicover.modularity_status(capacities, selected)
            )
            old_branch = multicover.old_branch(
                audit, capacities, marginals, histogram
            )
            old_to_layer_counter.setdefault(old_branch, Counter())[support_size] += 1
            row_layers[support_size] += 1
            certificate_records.append(
                (
                    index,
                    divisors,
                    q,
                    cofactor,
                    marginals,
                    histogram,
                    tuple(capacities),
                    old_branch,
                    selected,
                    weights,
                    gap,
                    tight,
                    submodular,
                    supermodular,
                    sub_witness,
                    super_witness,
                    len(certificates),
                )
            )
        layer_by_row[index] = tuple(sorted(row_layers.items()))

    certificate_records = tuple(certificate_records)
    status_sha = hashlib.sha256(repr(certificate_records).encode()).hexdigest()
    if EXPECTED_STATUS_SHA256 is not None:
        require(status_sha == EXPECTED_STATUS_SHA256, status_sha)
    layer_counts = tuple(
        sorted(Counter(len(record[8]) for record in certificate_records).items())
    )
    old_to_layer = tuple(
        (branch, tuple(sorted(counter.items())))
        for branch, counter in sorted(old_to_layer_counter.items())
    )
    by_row = tuple(
        (index, counts[index], layer_by_row[index]) for index in selected_indices
    )
    require(len(certificate_records) == EXPECTED_TOTALS[2], len(certificate_records))
    require(layer_counts == EXPECTED_LAYER_COUNTS, layer_counts)
    require(old_to_layer == EXPECTED_OLD_TO_LAYER, old_to_layer)
    require(by_row == EXPECTED_BY_ROW, by_row)

    instance_counts = Counter(
        (record[2], record[4], record[5], record[6])
        for record in certificate_records
    )
    multiplicities = tuple(sorted(Counter(instance_counts.values()).items()))
    require(len(instance_counts) == EXPECTED_UNIQUE_INSTANCES, len(instance_counts))
    require(multiplicities == EXPECTED_INSTANCE_MULTIPLICITIES, multiplicities)

    multi_records = tuple(
        record for record in certificate_records if len(record[8]) > 1
    )
    require(len(multi_records) == 1, multi_records)
    multi = multi_records[0]
    require((multi[0], multi[1]) == EXPECTED_MULTI_ADDRESS, multi[:2])
    require((multi[8], multi[9], multi[10]) == EXPECTED_MULTI_CERTIFICATE,
            multi[8:11])
    require(multi[11] == (0, 1, 2, 4, 8, 10), multi[11])

    thresholds, _demands, _good = multicover.tail_data(multi[6], multi[5])
    proper_support_tables = []
    for support_size in range(1, len(multi[8])):
        for selected in itertools.combinations(thresholds, support_size):
            table = multicover.exact_feasible_table(
                multi[2], multi[4], multi[6], multi[5], selected
            )
            require(table is not None, (multi[:2], selected, "proper support"))
            proper_support_tables.append((selected, table))
    proper_support_tables = tuple(proper_support_tables)
    require(len(proper_support_tables) == EXPECTED_PROPER_SUPPORTS,
            len(proper_support_tables))
    minimality_sha = hashlib.sha256(
        repr(proper_support_tables).encode()
    ).hexdigest()
    if EXPECTED_MINIMALITY_SHA256 is not None:
        require(minimality_sha == EXPECTED_MINIMALITY_SHA256, minimality_sha)

    # Hostile residual: row 64 has genuine status kills but also eight screen
    # survivors.  Its first survivor has an explicit exact common table and no
    # modular-cover contradiction.
    hostile_index = multicover.HOSTILE_INDEX
    eng = audit.status_engine(natural, "selected_hostile64")
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    stream = eng.ray.Stream(rows[hostile_index][0])
    _trials, states, _checks, _signs = eng.ray.ray_quotient_states(stream)
    crude, status, residual = eng.exact_common_status_screen(stream, states)
    hostile_counts = (len(states), len(crude), len(status), len(residual))
    require(hostile_counts == EXPECTED_HOSTILE_COUNTS, hostile_counts)
    hostile_status = []
    for divisors, witness in sorted(status.items()):
        q, _cofactor, marginals, _capacity_set, histogram, raw_dual = witness
        divisor_lcm = lcm(*divisors)
        _rebuilt, capacities = eng.ray.local.hunter_status_data(
            divisor_lcm, divisors, q
        )
        eng.independent_farkas_check(
            q, marginals, capacities, histogram, raw_dual
        )
        support_size, certificates = multicover.minimal_binary_multicovers(
            capacities, marginals, histogram
        )
        require(support_size == 1 and certificates, (hostile_index, divisors))
        hostile_status.append((divisors, certificates[0]))
    hostile_divisors = multicover.HOSTILE_DIVISORS
    require(hostile_divisors in residual, residual)
    divisor_lcm = lcm(*hostile_divisors)
    arcs = eng.ray.fibre.projected_support_arcs(divisor_lcm, stream.ranges)
    hostile_marginals, hostile_capacities = eng.ray.local.hunter_status_data(
        divisor_lcm, hostile_divisors, multicover.HOSTILE_Q
    )
    hostile_histogram = eng.ray.fibre.residue_load_histogram(
        arcs, multicover.HOSTILE_Q
    )
    require(hostile_marginals == multicover.HOSTILE_MARGINALS,
            hostile_marginals)
    require(tuple(hostile_capacities) == multicover.HOSTILE_CAPACITIES,
            hostile_capacities)
    require(hostile_histogram == multicover.HOSTILE_HISTOGRAM,
            hostile_histogram)
    hostile_thresholds, hostile_demands, hostile_good = multicover.tail_data(
        hostile_capacities, hostile_histogram
    )
    require(
        multicover.verify_table(
            multicover.HOSTILE_Q,
            hostile_marginals,
            hostile_good,
            hostile_demands,
            hostile_thresholds,
            tuple(map(Fraction, multicover.HOSTILE_TABLE)),
        ),
        "hostile common table",
    )
    hostile_support, hostile_certs = multicover.minimal_binary_multicovers(
        hostile_capacities, hostile_marginals, hostile_histogram
    )
    require(hostile_support is None and not hostile_certs,
            "hostile falsely killed")
    hostile_packet = (
        hostile_counts,
        tuple(hostile_status),
        hostile_divisors,
        hostile_marginals,
        hostile_histogram,
        multicover.HOSTILE_TABLE,
    )
    hostile_sha = hashlib.sha256(repr(hostile_packet).encode()).hexdigest()
    if EXPECTED_HOSTILE_SHA256 is not None:
        require(hostile_sha == EXPECTED_HOSTILE_SHA256, hostile_sha)

    # Positive arity control: the inherited row-138 packet genuinely needs
    # three layers, so the compiler is not hard-wired to arity at most two.
    arity3_index = 138
    arity3_divisors = (18, 49, 196, 882)
    eng = audit.status_engine(natural, "selected_arity3")
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    stream = eng.ray.Stream(rows[arity3_index][0])
    _trials, states, _checks, _signs = eng.ray.ray_quotient_states(stream)
    _crude, status, _residual = eng.exact_common_status_screen(stream, states)
    witness = status[arity3_divisors]
    q, _cofactor, marginals, _capacity_set, histogram, _raw_dual = witness
    divisor_lcm = lcm(*arity3_divisors)
    _rebuilt, capacities = eng.ray.local.hunter_status_data(
        divisor_lcm, arity3_divisors, q
    )
    arity3_support, arity3_certificates = (
        multicover.minimal_binary_multicovers(capacities, marginals, histogram)
    )
    require(arity3_support == 3 and arity3_certificates, arity3_certificates)
    require(arity3_certificates[0] == EXPECTED_ARITY3, arity3_certificates[0])

    remaining = tuple(index for index in live if index not in set(selected_indices))
    remaining_ranked = audit.ranked_families(rows, components, remaining)
    require((len(remaining), len(remaining_ranked)) == (WALL_AFTER, FAMILIES_AFTER),
            (len(remaining), len(remaining_ranked)))
    require(LEDGER_BEFORE - len(selected_indices) == LEDGER_AFTER, "ledger")
    require(WALL_BEFORE - len(selected_indices) == WALL_AFTER, "wall")
    require(FAMILIES_BEFORE - len(selected_families) == FAMILIES_AFTER,
            "families")

    semantic_packet = (
        AUDIT_SOURCE_SHA256,
        AUDIT_OUTPUT_SHA256,
        AUDIT_SEMANTIC_SHA256,
        MULTICOVER_SOURCE_SHA256,
        MULTICOVER_OUTPUT_SHA256,
        MULTICOVER_SEMANTIC_SHA256,
        census,
        first_indices,
        second_indices,
        ranked,
        selected_packet,
        canonical,
        by_row,
        certificate_records,
        multiplicities,
        proper_support_tables,
        hostile_packet,
        arity3_certificates[0],
        remaining_ranked,
        LEDGER_BEFORE,
        LEDGER_AFTER,
        WALL_BEFORE,
        WALL_AFTER,
        FAMILIES_BEFORE,
        FAMILIES_AFTER,
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    family_counts = {}
    for family in selected_families:
        family_counts[(family[2], family[3])] = sum_counts(
            tuple(counts[index] for index in family[4])
        )
    multi_summary = (
        multi[0], multi[1], multi[2], multi[4], multi[5], multi[6],
        multi[8], multi[9], multi[10], multi[11], multi[12], multi[14],
    )
    lines = [
        "LRC14 projected-k3 z216 third intrinsic-ruler prefix multicover closure",
        (
            f"dependency=two_prefix_audit:{AUDIT_SOURCE_SHA256}/"
            f"{AUDIT_OUTPUT_SHA256}/{AUDIT_SEMANTIC_SHA256};"
            f"multicover_compiler:{MULTICOVER_SOURCE_SHA256}/"
            f"{MULTICOVER_OUTPUT_SHA256}/{MULTICOVER_SEMANTIC_SHA256}"
        ),
        (
            f"universe=z1:{LEVEL};atlas_rows:{census[0]};wall:{census[1]};"
            f"order:{census[2]};prior_prefix_rows:{len(first_indices) + len(second_indices)};"
            f"live_wall:{census[3]};live_families:{census[4]}"
        ),
        (
            f"queue_head:{family_head(ranked, 3)};membership_and_complete_"
            "family_reconstruction:PASS;strict_cost_boundary:PASS"
        ),
        (
            f"selection=families:{len(selected_families)};rows:"
            f"{len(selected_indices)};indices:{selected_indices};cost:"
            f"{sum(family[0] for family in selected_families)};"
            f"packet_sha256:{selected_packet_sha}"
        ),
    ]
    for family in selected_families:
        summary = family_counts[(family[2], family[3])]
        lines.append(
            f"family=gcd{family[2]}/L{family[3]};rows:{family[1]};"
            f"indices:{family[4]};cost:{family[0]};states:{summary[0]};"
            f"crude:{summary[1]};status:{summary[2]};residual:{summary[3]}"
        )
    lines.extend(
        [
            (
                f"selected_screen=states:{EXPECTED_TOTALS[0]};"
                f"crude:{EXPECTED_TOTALS[1]};status:{EXPECTED_TOTALS[2]};"
                f"residual:{EXPECTED_TOTALS[3]};direct_farkas:{direct_total};"
                f"legacy_farkas:{legacy_total};by_row:{by_row};"
                f"screen_sha256:{screen_sha}"
            ),
            (
                f"multicover_taxonomy=layer_counts:{layer_counts};"
                f"old_to_layer:{old_to_layer};unique_instances:"
                f"{len(instance_counts)};instance_multiplicities:{multiplicities};"
                f"status_sha256:{status_sha}"
            ),
            (
                f"two_layer_circuit={multi_summary};all_singleton_tail_"
                f"supports_exactly_feasible:{len(proper_support_tables)};"
                f"minimality_sha256:{minimality_sha}"
            ),
            (
                "multicover_mechanism=sum_of_selected_nested_tail_indicators_"
                "is_pointwise_bounded_by_a_zero_constant_nonnegative_integral_"
                "coordinate_form;one_common_status_table_then_forces_tail_"
                "demand_sum_at_most_weighted_marginal_sum"
            ),
            (
                f"hostile=row64_counts:{hostile_counts};all_57_status_kills_"
                "one_layer;first_residual_has_exact_common_table:PASS;no_"
                f"multicover_contradiction:PASS;hostile_sha256:{hostile_sha}"
            ),
            (
                f"arity_control=row138_divisors:{arity3_divisors};"
                f"certificate:{arity3_certificates[0]};compiler_reaches_"
                "genuine_three_layer_support:PASS"
            ),
            (
                f"consequence=ledger:{LEDGER_BEFORE}-{len(selected_indices)}="
                f"{LEDGER_AFTER};z216_wall:{WALL_BEFORE}-"
                f"{len(selected_indices)}={WALL_AFTER};families:"
                f"{FAMILIES_BEFORE}-{len(selected_families)}="
                f"{FAMILIES_AFTER};order:0;cap:216"
            ),
            (
                f"stopping_boundary=remaining_wall:{len(remaining)};"
                f"remaining_families:{len(remaining_ranked)};next_head:"
                f"{family_head(remaining_ranked, 3)}"
            ),
            (
                "direction=empty_exact_necessary_upper_screen_excludes_each_"
                "selected_projected_wall_row;screen_feasibility_is_not_a_"
                "converse_and_intrinsic_cost_is_only_a_work_order"
            ),
            (
                "scope=projected_k3_z216_necessary_atlas_only;no_physical_"
                "speed_entry_endpoint_owner_phase_current_arbitrary_k_le_1_"
                "rung_or_LRC14_claim"
            ),
            f"truth_gates=assert_nodes:{assert_nodes};float_literals:{float_nodes}",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        ]
    )
    payload = "\n".join(lines) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
