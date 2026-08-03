#!/usr/bin/env python3
"""Independent exact audit of THM-3320 at projected k=3, z1=216.

The script reconstructs the 353-row wall after THM-3313, selects the next
whole-family cost prefix, independently screens its four rows, and compiles
every exact common-status exclusion into a deterministic threshold-layer
modular multicover.  It remains entirely inside the necessary projected
atlas and then compares the result with the separately produced canonical
companion.  No physical, arbitrary-k, rung, or LRC(14) conclusion is asserted.
"""

from __future__ import annotations

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
THIRD_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_multicover_closure_scout_20260803.py"
)
THIRD_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_multicover_closure_scout_20260803.out"
)
CANON_THEOREM = ROOT / (
    "01-canon/theorems/"
    "THM-3320-projected-k3-z216-fourth-ruler-prefix-and-affine-multicover-closure.md"
)
CANON_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.py"
)
CANON_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.out"
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
THIRD_SOURCE_SHA256 = (
    "34710b907f3274057d6cb60a0fd7bd72ae2c3879a90b830a3e376d30e5198646"
)
THIRD_OUTPUT_SHA256 = (
    "c664acf4e02b40835dfdc76ff135d17332e3d0e57b489e61b19c0f2f7b5998f5"
)
THIRD_SEMANTIC_SHA256 = (
    "4d8547f53b287d82bb725cd6b3ed78c390ba03a37a3dcb3ab76ac98389e045a3"
)
CANON_THEOREM_SHA256 = (
    "be64218cace2d1b67b9be38e1f6e4405d192406c1e0957888f284cb21573a0b1"
)
CANON_SOURCE_SHA256 = (
    "b515c70174d58ad859a08c29949cdee36a4e04122451a66237732570ab5ee213"
)
CANON_OUTPUT_SHA256 = (
    "d1845611f27d427a1d38afe349ed07bc964590fd39a3e88cd45e6ea34a86bc38"
)
CANON_SEMANTIC_SHA256 = (
    "c201276eb84f71806a7eb683e42d723b930f7cd0f79862e8d6b1e0ad07d37dd9"
)

LEVEL = 216
PRIOR_PREFIX_INDICES = (
    (210, 50, 18, 140, 299),
    (134, 121, 17, 27, 46, 56, 79, 89, 126, 138, 205, 248, 258, 286,
     298, 347, 376, 425),
    (157, 53, 293, 357),
)
EXPECTED_SELECTED_FAMILIES = (
    (
        2_293_200,
        1,
        24,
        76_440,
        (141,),
        ((141, (1, 4, 6, 10, 13, 14), 30, 2_293_200),),
    ),
    (
        2_629_536,
        3,
        24,
        30_576,
        (133, 219, 359),
        (
            (133, (1, 4, 6, 8, 13, 14), 26, 794_976),
            (219, (1, 6, 8, 12, 13, 14), 30, 917_280),
            (359, (2, 6, 8, 12, 13, 14), 30, 917_280),
        ),
    ),
)
EXPECTED_SELECTED_INDICES = (141, 133, 219, 359)
EXPECTED_BY_ROW = (
    (141, (10, 8, 2, 0), ((1, 2),)),
    (133, (10, 0, 10, 0), ((1, 7), (2, 1), (3, 2))),
    (219, (189, 147, 42, 0), ((1, 42),)),
    (359, (21, 17, 4, 0), ((1, 4),)),
)
EXPECTED_TOTALS = (230, 172, 58, 0)
EXPECTED_LAYER_COUNTS = ((1, 55), (2, 1), (3, 2))
EXPECTED_OLD_TO_LAYER = (
    ("coordinate_union", ((1, 55),)),
    ("weighted_core", ((2, 1), (3, 2))),
)
EXPECTED_NEXT_FAMILY = (3_400_992, 19, 72, 7_056)
EXPECTED_SEMANTIC_SHA256 = (
    "1ecc946ac1eaf9731b16b377ada123e3bb38437736a5897b6644fb415c587a35"
)


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


def negative_alpha_rejected(engine, q, marginals, capacities, histogram,
                            certificate):
    thresholds, alpha, z = certificate
    positive = next((i for i, value in enumerate(alpha) if value > 0), None)
    require(positive is not None, "dual has no positive threshold coefficient")
    altered = list(alpha)
    altered[positive] = -altered[positive]
    try:
        engine.independent_farkas_check(
            q,
            marginals,
            capacities,
            histogram,
            (thresholds, tuple(altered), z),
        )
    except RuntimeError:
        return True
    return False


def permuted_signature(marginals, capacities, permutation):
    pulled_marginals = tuple(marginals[permutation[j]] for j in range(4))
    pulled_capacities = []
    for new_pattern in range(16):
        old_pattern = sum(
            ((new_pattern >> j) & 1) << permutation[j] for j in range(4)
        )
        pulled_capacities.append(capacities[old_pattern])
    return pulled_marginals, tuple(pulled_capacities)


def orbit_signature(q, marginals, histogram, capacities):
    return min(
        (q, *permuted_signature(marginals, capacities, permutation), histogram)
        for permutation in itertools.permutations(range(4))
    )


def main() -> None:
    dependencies = (
        (AUDIT_SOURCE, AUDIT_SOURCE_SHA256, None),
        (AUDIT_OUTPUT, AUDIT_OUTPUT_SHA256, AUDIT_SEMANTIC_SHA256),
        (MULTICOVER_SOURCE, MULTICOVER_SOURCE_SHA256, None),
        (
            MULTICOVER_OUTPUT,
            MULTICOVER_OUTPUT_SHA256,
            MULTICOVER_SEMANTIC_SHA256,
        ),
        (THIRD_SOURCE, THIRD_SOURCE_SHA256, None),
        (THIRD_OUTPUT, THIRD_OUTPUT_SHA256, THIRD_SEMANTIC_SHA256),
        (CANON_THEOREM, CANON_THEOREM_SHA256, None),
        (CANON_SOURCE, CANON_SOURCE_SHA256, None),
        (CANON_OUTPUT, CANON_OUTPUT_SHA256, CANON_SEMANTIC_SHA256),
    )
    for path, digest, semantic in dependencies:
        require(lf_sha(path) == digest, (path, "dependency changed"))
        if semantic is not None:
            require(
                f"semantic_sha256={semantic}"
                in path.read_text(encoding="utf-8"),
                (path, "dependency semantic changed"),
            )
    canonical_text = CANON_OUTPUT.read_text(encoding="utf-8")

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require((assert_nodes, float_nodes) == (0, 0), "truth-gate syntax")

    audit = load("fourth_prefix_audit", AUDIT_SOURCE)
    multicover = load("fourth_prefix_multicover", MULTICOVER_SOURCE)
    natural = audit.load("fourth_prefix_natural", audit.NATURAL_SOURCE)
    order = audit.load("fourth_prefix_order", audit.ORDER_SOURCE)
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
    closed = gcd8 | gcd18 | order_indices | natural_indices
    require(len(closed) == 100, len(closed))

    prior_packets = []
    for prefix_number, expected_indices in enumerate(PRIOR_PREFIX_INDICES, 1):
        live = tuple(
            index for index, row in enumerate(rows)
            if row[3] and index not in closed
        )
        ranked = audit.ranked_families(rows, components, live)
        selected = audit.through_next_nonsingleton(ranked)
        selected_indices = tuple(
            index for family in selected for index in family[4]
        )
        require(
            selected_indices == expected_indices,
            (prefix_number, selected_indices),
        )
        prior_packets.append((len(live), len(ranked), selected))
        closed |= set(selected_indices)
    prior_packets = tuple(prior_packets)

    live = tuple(
        index for index, row in enumerate(rows)
        if row[3] and index not in closed
    )
    ranked = audit.ranked_families(rows, components, live)
    selected_families = audit.through_next_nonsingleton(ranked)
    selected_indices = tuple(
        index for family in selected_families for index in family[4]
    )
    require((len(rows), sum(row[3] for row in rows)) == (480, 447), "atlas")
    require((len(live), len(ranked)) == (353, 31), (len(live), len(ranked)))
    require(tuple(selected_families) == EXPECTED_SELECTED_FAMILIES,
            selected_families)
    require(selected_indices == EXPECTED_SELECTED_INDICES, selected_indices)
    require(selected_families[-1][0] < ranked[len(selected_families)][0],
            "cost boundary")
    next_family = ranked[len(selected_families)]
    require(next_family[:4] == EXPECTED_NEXT_FAMILY, next_family[:5])

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

    canonical = []
    counts = {}
    direct_total = 0
    inherited_total = 0
    for index in selected_indices:
        returned_index, screen_row, direct, inherited = audit.screen_worker(
            (index, (LEVEL, *rows[index]))
        )
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
        require(screen_row[16] == screen_row[11], (index, "exact status"))
        require(direct + inherited == screen_row[11], (index, "dual count"))
        counts[index] = tuple(
            screen_row[position] for position in (9, 10, 11, 12)
        )
        canonical.append((index, tuple(screen_row[:19])))
        direct_total += direct
        inherited_total += inherited
    canonical = tuple(canonical)
    require(
        sum_counts(tuple(counts[index] for index in selected_indices))
        == EXPECTED_TOTALS,
        counts,
    )
    require((direct_total, inherited_total) == (0, 58),
            (direct_total, inherited_total))
    require(all(counts[index][3] == 0 for index in selected_indices), counts)
    screen_sha = hashlib.sha256(repr(canonical).encode()).hexdigest()

    certificate_records = []
    layer_by_row = {}
    old_to_layer_counter = {}
    hostile_mutations = 0
    for index in selected_indices:
        engine = audit.status_engine(natural, f"fourth_selected_{index}")
        engine.FIRST = LEVEL
        engine.ray.FIRST = LEVEL
        stream = engine.ray.Stream(rows[index][0])
        require(
            (stream.L, stream.high_floor) == (rows[index][1], rows[index][2]),
            (index, "atlas mismatch"),
        )
        _trials, states, _checks, _signs = engine.ray.ray_quotient_states(stream)
        crude, status, residual = engine.exact_common_status_screen(stream, states)
        require(
            (len(states), len(crude), len(status), len(residual)) == counts[index],
            (index, "independent screen reconstruction"),
        )
        row_layers = Counter()
        mutated = False
        for divisors, witness in sorted(status.items()):
            q, cofactor, marginals, capacity_set, histogram, raw_dual = witness
            divisor_lcm = lcm(*divisors)
            require(divisor_lcm == q * cofactor, (index, divisors, "cofactor"))
            rebuilt_marginals, capacities = (
                engine.ray.local.hunter_status_data(divisor_lcm, divisors, q)
            )
            require(rebuilt_marginals == marginals,
                    (index, divisors, "marginals"))
            require(tuple(sorted(set(capacities))) == capacity_set,
                    (index, divisors, "capacities"))
            engine.independent_farkas_check(
                q, marginals, capacities, histogram, raw_dual
            )
            if not mutated:
                require(
                    negative_alpha_rejected(
                        engine,
                        q,
                        marginals,
                        capacities,
                        histogram,
                        raw_dual,
                    ),
                    (index, "negative-alpha mutation accepted"),
                )
                hostile_mutations += 1
                mutated = True

            support_size, certificates = multicover.minimal_binary_multicovers(
                capacities, marginals, histogram
            )
            require(support_size is not None and certificates,
                    (index, divisors, "uncompiled status"))
            selected, weights, gap = certificates[0]
            require(len(selected) == support_size, (index, divisors, selected))
            require(
                multicover.coordinate_minimal(capacities, selected, weights),
                (index, divisors, "coordinate-nonminimal"),
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
                    len(certificates),
                )
            )
        require(mutated, (index, "no status mutation control"))
        layer_by_row[index] = tuple(sorted(row_layers.items()))
    certificate_records = tuple(certificate_records)
    require(hostile_mutations == len(selected_indices), hostile_mutations)

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
    require(by_row == EXPECTED_BY_ROW, by_row)
    require(layer_counts == EXPECTED_LAYER_COUNTS, layer_counts)
    require(old_to_layer == EXPECTED_OLD_TO_LAYER, old_to_layer)
    require(len(certificate_records) == EXPECTED_TOTALS[2],
            len(certificate_records))

    instance_counts = Counter(
        (record[2], record[4], record[5], record[6])
        for record in certificate_records
    )
    require(len(instance_counts) == 32, len(instance_counts))
    multiplicities = tuple(sorted(Counter(instance_counts.values()).items()))

    multi_records = tuple(
        record for record in certificate_records if len(record[8]) > 1
    )
    require(len(multi_records) == 3, multi_records)
    multi_orbits = Counter(
        orbit_signature(record[2], record[4], record[5], record[6])
        for record in multi_records
    )
    require(tuple(sorted(multi_orbits.values())) == (1, 2), multi_orbits)

    # Verify THM-3320's repaired affine presentations directly on all sixteen
    # Boolean patterns, without importing its affine-majorant search routine.
    expected_affine = (
        ((1, 5), (Fraction(1, 2),) * 4 + (Fraction(3, 2),), 182),
        ((1, 5), (Fraction(1, 2),) * 4 + (Fraction(3, 2),), 182),
        ((2, 5), (Fraction(0), Fraction(1), Fraction(1), Fraction(2),
                   Fraction(0)), 510),
    )
    affine_audit = []
    for record, (selected, coefficients, expected_gap) in zip(
        multi_records, expected_affine
    ):
        q, marginals, histogram, capacities = (
            record[2], record[4], record[5], record[6]
        )
        values = multicover.layer_values(capacities, selected)
        affine_values = tuple(
            coefficients[0]
            + sum(
                coefficients[coordinate + 1]
                * ((pattern >> coordinate) & 1)
                for coordinate in range(4)
            )
            for pattern in range(16)
        )
        require(all(left >= right for left, right in zip(affine_values, values)),
                (record[:2], "affine pointwise failure"))
        _thresholds, demands, _good = multicover.tail_data(
            capacities, histogram
        )
        demand = sum(demands[threshold] for threshold in selected)
        cost = coefficients[0] * q + sum(
            coefficients[coordinate + 1] * marginals[coordinate]
            for coordinate in range(4)
        )
        gap = demand - cost
        require(gap == expected_gap, (record[:2], demand, cost, gap))
        affine_audit.append((record[0], record[1], selected, coefficients,
                             demand, cost, gap))
    affine_audit = tuple(affine_audit)

    # Audit the gap between the least *binary integral multicover* support and
    # the least tail support making the rational common-table polytope empty.
    # The two three-layer binary packets already fail on the pair (1,5), even
    # though every singleton tail is feasible.  Thus their rational support
    # rank is two: this is a genuine integrality phenomenon, not a claim that
    # three tail inequalities are necessary in the full dual cone.
    proper_support_audit = []
    for record_number, record in enumerate(multi_records):
        q, marginals, histogram, capacities = (
            record[2], record[4], record[5], record[6]
        )
        selected = record[8]
        for support_size in range(1, len(selected)):
            for proper in itertools.combinations(selected, support_size):
                table = multicover.exact_feasible_table(
                    q, marginals, capacities, histogram, proper
                )
                proper_support_audit.append(
                    (record_number, proper, table is not None, table)
                )
    proper_support_audit = tuple(proper_support_audit)
    require(len(proper_support_audit) == 14, len(proper_support_audit))
    infeasible_proper = tuple(
        (record_number, proper)
        for record_number, proper, feasible, _table in proper_support_audit
        if not feasible
    )
    require(infeasible_proper == ((0, (1, 5)), (1, (1, 5))),
            infeasible_proper)
    singleton_feasible = tuple(
        all(
            feasible
            for candidate, proper, feasible, _table in proper_support_audit
            if candidate == record_number and len(proper) == 1
        )
        for record_number in range(len(multi_records))
    )
    require(singleton_feasible == (True, True, True), singleton_feasible)
    rational_tail_supports = (2, 2, 2)
    minimality_sha = hashlib.sha256(
        repr(proper_support_audit).encode()
    ).hexdigest()

    singleton_table_audit = []
    for record_number, record in enumerate(multi_records):
        q, marginals, histogram, capacities = (
            record[2], record[4], record[5], record[6]
        )
        thresholds, _demands, _good = multicover.tail_data(
            capacities, histogram
        )
        for threshold in thresholds:
            table = multicover.exact_feasible_table(
                q, marginals, capacities, histogram, (threshold,)
            )
            require(table is not None,
                    (record[:2], threshold, "singleton infeasible"))
            singleton_table_audit.append((record_number, threshold, table))
    singleton_table_audit = tuple(singleton_table_audit)
    require(len(singleton_table_audit) == 18, len(singleton_table_audit))
    singleton_sha = hashlib.sha256(
        repr(singleton_table_audit).encode()
    ).hexdigest()

    empty_support_audit = []
    for signature in sorted(instance_counts):
        q, marginals, histogram, capacities = signature
        table = multicover.exact_feasible_table(
            q, marginals, capacities, histogram, ()
        )
        require(table is not None, (signature, "empty support infeasible"))
        empty_support_audit.append((signature, table))
    empty_support_audit = tuple(empty_support_audit)
    require(len(empty_support_audit) == 32, len(empty_support_audit))
    empty_support_sha = hashlib.sha256(
        repr(empty_support_audit).encode()
    ).hexdigest()

    # A canonical negative control: row 64 retains eight passports.  Rebuild
    # its stored exact common table and verify the compiler does not kill it.
    hostile_index = multicover.HOSTILE_INDEX
    engine = audit.status_engine(natural, "fourth_hostile64")
    engine.FIRST = LEVEL
    engine.ray.FIRST = LEVEL
    stream = engine.ray.Stream(rows[hostile_index][0])
    _trials, states, _checks, _signs = engine.ray.ray_quotient_states(stream)
    crude, status, residual = engine.exact_common_status_screen(stream, states)
    hostile_counts = (len(states), len(crude), len(status), len(residual))
    require(hostile_counts == (115, 50, 57, 8), hostile_counts)
    hostile_divisors = multicover.HOSTILE_DIVISORS
    require(hostile_divisors in residual, residual)
    divisor_lcm = lcm(*hostile_divisors)
    hostile_marginals, hostile_capacities = (
        engine.ray.local.hunter_status_data(
            divisor_lcm, hostile_divisors, multicover.HOSTILE_Q
        )
    )
    arcs = engine.ray.fibre.projected_support_arcs(divisor_lcm, stream.ranges)
    hostile_histogram = engine.ray.fibre.residue_load_histogram(
        arcs, multicover.HOSTILE_Q
    )
    require(hostile_marginals == multicover.HOSTILE_MARGINALS,
            hostile_marginals)
    require(tuple(hostile_capacities) == multicover.HOSTILE_CAPACITIES,
            hostile_capacities)
    require(hostile_histogram == multicover.HOSTILE_HISTOGRAM,
            hostile_histogram)
    thresholds, demands, good_rows = multicover.tail_data(
        hostile_capacities, hostile_histogram
    )
    require(
        multicover.verify_table(
            multicover.HOSTILE_Q,
            hostile_marginals,
            good_rows,
            demands,
            thresholds,
            tuple(multicover.Fraction(value) for value in multicover.HOSTILE_TABLE),
        ),
        "hostile exact common table",
    )
    hostile_support, hostile_certificates = (
        multicover.minimal_binary_multicovers(
            hostile_capacities, hostile_marginals, hostile_histogram
        )
    )
    require(hostile_support is None and not hostile_certificates,
            "hostile falsely compiled")
    hostile_packet = (
        hostile_counts,
        hostile_divisors,
        hostile_marginals,
        hostile_histogram,
        multicover.HOSTILE_TABLE,
    )
    hostile_sha = hashlib.sha256(repr(hostile_packet).encode()).hexdigest()

    # Opposite hostile control from THM-3308: this packet really has rational
    # tail rank three.  Every singleton and pair is exactly feasible.
    arity3_index = 138
    arity3_divisors = (18, 49, 196, 882)
    engine = audit.status_engine(natural, "fourth_audit_arity3")
    engine.FIRST = LEVEL
    engine.ray.FIRST = LEVEL
    stream = engine.ray.Stream(rows[arity3_index][0])
    _trials, states, _checks, _signs = engine.ray.ray_quotient_states(stream)
    _crude, status, _residual = engine.exact_common_status_screen(stream, states)
    q, _cofactor, marginals, _capacity_set, histogram, raw_dual = (
        status[arity3_divisors]
    )
    _rebuilt, capacities = engine.ray.local.hunter_status_data(
        lcm(*arity3_divisors), arity3_divisors, q
    )
    engine.independent_farkas_check(
        q, marginals, capacities, histogram, raw_dual
    )
    arity3_support, arity3_certificates = (
        multicover.minimal_binary_multicovers(
            capacities, marginals, histogram
        )
    )
    require(
        arity3_support == 3
        and arity3_certificates
        and arity3_certificates[0] == ((1, 4, 5), (1, 3, 1, 1), 11),
        arity3_certificates,
    )
    thresholds, _demands, _good = multicover.tail_data(capacities, histogram)
    arity3_proper = []
    for support_size in (1, 2):
        for selected in itertools.combinations(thresholds, support_size):
            table = multicover.exact_feasible_table(
                q, marginals, capacities, histogram, selected
            )
            require(table is not None,
                    (arity3_divisors, selected, "arity3 proper infeasible"))
            arity3_proper.append((selected, table))
    arity3_proper = tuple(arity3_proper)
    require(len(arity3_proper) == 21, len(arity3_proper))
    arity3_sha = hashlib.sha256(repr(arity3_proper).encode()).hexdigest()

    remaining = tuple(
        index for index in live if index not in set(selected_indices)
    )
    remaining_ranked = audit.ranked_families(rows, components, remaining)
    require((len(remaining), len(remaining_ranked)) == (349, 29),
            (len(remaining), len(remaining_ranked)))
    require(373_157 - len(selected_indices) == 373_153, "ledger")

    status_sha = hashlib.sha256(
        repr(certificate_records).encode()
    ).hexdigest()
    multi_summary = tuple(
        (
            record[0],
            record[1],
            record[2],
            record[4],
            record[5],
            record[6],
            record[8],
            record[9],
            record[10],
        )
        for record in multi_records
    )
    canonical_comparisons = (
        (
            "queue",
            f"queue_head:{family_head(ranked, 3)};" in canonical_text,
        ),
        (
            "selection",
            f"indices:{selected_indices};cost:4922736;packet_sha256:"
            f"{selected_packet_sha}" in canonical_text,
        ),
        (
            "screen",
            f"screen_sha256:{screen_sha}" in canonical_text
            and "selected_screen=states:230;crude:172;status:58;residual:0;"
            in canonical_text,
        ),
        (
            "integral_taxonomy",
            "integral_modular_taxonomy=layer_counts:((1, 55), (2, 1), "
            "(3, 2))" in canonical_text,
        ),
        (
            "full_dual_taxonomy",
            "full_common_table_taxonomy=layer_counts:((1, 55), (2, 3))"
            in canonical_text,
        ),
        (
            "affine_gaps",
            "Fraction(182, 1), 3458, Fraction(3276, 1)" in canonical_text
            and "Fraction(510, 1), 3058, Fraction(2548, 1)" in canonical_text,
        ),
        (
            "hostile64",
            "hostile=row64_counts:(115, 50, 57, 8)" in canonical_text,
        ),
        (
            "arity3_control",
            "certificate:((1, 4, 5), (1, 3, 1, 1), 11);"
            "proper_support_tables:21" in canonical_text,
        ),
        (
            "ledger",
            "consequence=ledger:373157-4=373153;z216_wall:353-4=349;"
            "families:31-2=29;order:0;cap:216" in canonical_text,
        ),
        (
            "scope",
            "scope=projected_k3_z216_necessary_atlas_only" in canonical_text,
        ),
    )
    discrepancies = tuple(
        name for name, agrees in canonical_comparisons if not agrees
    )
    require(not discrepancies, discrepancies)
    semantic_packet = (
        tuple(digest for _path, digest, _semantic in dependencies),
        tuple(semantic for _path, _digest, semantic in dependencies),
        hashlib.sha256(repr(rows).encode()).hexdigest(),
        prior_packets,
        ranked,
        selected_packet,
        canonical,
        by_row,
        certificate_records,
        affine_audit,
        multiplicities,
        tuple(sorted(multi_orbits.values())),
        proper_support_audit,
        singleton_table_audit,
        empty_support_audit,
        hostile_packet,
        arity3_certificates[0],
        arity3_proper,
        remaining_ranked,
        canonical_comparisons,
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    family_counts = {}
    for family in selected_families:
        family_counts[(family[2], family[3])] = sum_counts(
            tuple(counts[index] for index in family[4])
        )
    lines = [
        "Independent VERIFIED-EXACT cross-audit of THM-3320 z216 fourth prefix",
        (
            f"canonical=theorem:{CANON_THEOREM_SHA256};source:"
            f"{CANON_SOURCE_SHA256};output:{CANON_OUTPUT_SHA256};semantic:"
            f"{CANON_SEMANTIC_SHA256}"
        ),
        (
            f"universe=z1:{LEVEL};atlas_rows:{len(rows)};wall:447;order:33;"
            f"prior_prefix_rows:{sum(map(len, PRIOR_PREFIX_INDICES))};"
            f"live_wall:{len(live)};live_families:{len(ranked)};"
            f"row_sha256:{hashlib.sha256(repr(rows).encode()).hexdigest()};"
            "dual_parser_match:PASS"
        ),
        (
            f"queue_head:{family_head(ranked, 3)};complete_family_"
            "reconstruction:PASS;strict_cost_boundary:PASS"
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
                f"inherited_farkas:{inherited_total};by_row:{by_row};"
                f"screen_sha256:{screen_sha}"
            ),
            (
                f"multicover=layer_counts:{layer_counts};old_to_layer:"
                f"{old_to_layer};unique_instances:{len(instance_counts)};"
                f"instance_multiplicities:{multiplicities};"
                f"status_sha256:{status_sha}"
            ),
            (
                f"multi_layer_records:{multi_summary};coordinate_orbits:"
                f"{len(multi_orbits)};orbit_multiplicities:"
                f"{tuple(sorted(multi_orbits.values()))};binary_supports:"
                f"{tuple(len(record[8]) for record in multi_records)};"
                f"rational_tail_supports:{rational_tail_supports};proper_"
                f"support_tests:{len(proper_support_audit)};infeasible_"
                f"proper_supports:{infeasible_proper};"
                f"minimality_sha256:{minimality_sha}"
            ),
            (
                f"affine_audit:{affine_audit};all_16_pattern_majorants:PASS;"
                f"all_18_singleton_tables:PASS;singleton_sha256:"
                f"{singleton_sha};all_32_empty_tail_tables:PASS;"
                f"empty_support_sha256:{empty_support_sha}"
            ),
            (
                "integrality_connection=the_two_S4_equivalent_three_layer_"
                "binary_multicovers_have_exact_rational_tail_rank_two;the_"
                "pair_(1,5)_is_infeasible_while_each_singleton_is_feasible"
            ),
            (
                "orbit_connection=the_two_three_layer_packets_are_exact_"
                "coordinate_permutations;cache_by_S4_orbit_can_reuse_one_"
                "capacity_signature_on_the_next_19_row_family"
            ),
            (
                f"hostile=row64_counts:{hostile_counts};stored_exact_common_"
                "table_rechecked:PASS;no_multicover_contradiction:PASS;"
                f"hostile_sha256:{hostile_sha};negative_alpha_mutations_"
                f"rejected:{hostile_mutations}"
            ),
            (
                f"arity3_control=row138_divisors:{arity3_divisors};"
                f"certificate:{arity3_certificates[0]};proper_support_tables:"
                f"{len(arity3_proper)};all_exactly_feasible:PASS;"
                f"minimality_sha256:{arity3_sha}"
            ),
            (
                "direction=empty_exact_necessary_upper_screen_excludes_each_"
                "selected_projected_wall_row;screen_feasibility_is_not_a_"
                "converse_and_intrinsic_cost_is_only_a_work_order"
            ),
            (
                f"consequence=ledger:373157-{len(selected_indices)}=373153;"
                f"z216_wall:353-{len(selected_indices)}=349;families:31-"
                f"{len(selected_families)}=29;order:0;cap:216"
            ),
            (
                f"stopping_boundary=remaining_wall:{len(remaining)};"
                f"remaining_families:{len(remaining_ranked)};next_head:"
                f"{family_head(remaining_ranked, 3)}"
            ),
            (
                "scope=projected_k3_z216_necessary_atlas_only;no_physical_"
                "speed_entry_endpoint_owner_phase_current_arbitrary_k_le_1_"
                "rung_or_LRC14_claim"
            ),
            (
                f"canonical_comparison:{canonical_comparisons};"
                f"discrepancies:{discrepancies};verdict:INDEPENDENT_"
                "VERIFIED_EXACT"
            ),
            f"truth_gates=assert_nodes:{assert_nodes};float_literals:{float_nodes}",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        ]
    )
    print("\n".join(lines))


if __name__ == "__main__":
    main()
