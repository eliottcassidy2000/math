#!/usr/bin/env python3
"""Exact fourth z216 intrinsic-ruler prefix with affine multicovers.

THM-3313 closes the third deterministic projected-k3, z1=216 cost prefix.
This scout reconstructs that queue state from canonical artifacts, selects the
next complete-family prefix, and runs the inherited exact ray/status screen.

Every status kill is first passed through THM-3308's zero-constant integral
threshold-layer modular compiler.  Proper-support primal controls then audit
the actual common-table dual arity.  When the integral template overstates
that arity, an exact affine Boolean majorant

    sum_(t in T) 1[c(P)>=t] <= a + sum_(i in P) w_i

is derived by rational vertex enumeration.  Solver-selected inherited Farkas
vectors are checked but never persisted.

Scope is the necessary projected-k3, z1=216 atlas only.  Nothing here restores
physical entry, endpoint origin, owner, phase, current, arbitrary k<=1, the
rung, or LRC(14).
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
THM3313 = ROOT / (
    "01-canon/theorems/"
    "THM-3313-projected-k3-z216-third-ruler-cost-prefix-multicover-closure.md"
)
THIRD_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_"
    "multicover_closure_scout_20260803.py"
)
THIRD_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_"
    "multicover_closure_scout_20260803.out"
)
MULTICOVER_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_threshold_layer_multicover_dual_anatomy_scout_20260803.py"
)
MULTICOVER_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_threshold_layer_multicover_dual_anatomy_scout_20260803.out"
)

THM3313_SHA256 = (
    "88f62f7ecdaa5f897cd8ec78bdf3a352c964ca0a7ccfbfde05e28cd011f1cfff"
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
EXPECTED_CENSUS = (480, 447, 33, 353, 31)
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
EXPECTED_NEXT_FAMILY = (
    3_400_992,
    19,
    72,
    7_056,
    (21, 49, 62, 92, 129, 148, 202, 203, 212, 254, 262, 290, 303, 343,
     344, 352, 387, 428, 465),
)
EXPECTED_SELECTED_INDICES = (141, 133, 219, 359)
EXPECTED_BY_ROW = (
    (141, (10, 8, 2, 0), ((1, 2),), ((1, 2),)),
    (133, (10, 0, 10, 0), ((1, 7), (2, 1), (3, 2)), ((1, 7), (2, 3))),
    (219, (189, 147, 42, 0), ((1, 42),), ((1, 42),)),
    (359, (21, 17, 4, 0), ((1, 4),), ((1, 4),)),
)
EXPECTED_TOTALS = (230, 172, 58, 0)
EXPECTED_FARKAS_COUNTS = (0, 58)
EXPECTED_MODULAR_LAYER_COUNTS = ((1, 55), (2, 1), (3, 2))
EXPECTED_FULL_LAYER_COUNTS = ((1, 55), (2, 3))
EXPECTED_OLD_TO_MODULAR = (
    ("coordinate_union", ((1, 55),)),
    ("weighted_core", ((2, 1), (3, 2))),
)
EXPECTED_UNIQUE_INSTANCES = 32
EXPECTED_INSTANCE_MULTIPLICITIES = (
    (1, 24), (2, 2), (3, 3), (5, 1), (6, 1), (10, 1),
)
EXPECTED_MULTI_ADDRESSES = (
    (133, (84, 588, 728, 1274)),
    (133, (84, 728, 1092, 1274)),
    (133, (84, 728, 1274, 2184)),
)
EXPECTED_FULL_CERTIFICATES = (
    ((1, 5), (Fraction(1, 2), Fraction(1, 2), Fraction(1, 2),
              Fraction(1, 2), Fraction(3, 2)), Fraction(182)),
    ((1, 5), (Fraction(1, 2), Fraction(1, 2), Fraction(1, 2),
              Fraction(1, 2), Fraction(3, 2)), Fraction(182)),
    ((2, 5), (Fraction(0), Fraction(1), Fraction(1), Fraction(2),
              Fraction(0)), Fraction(510)),
)
EXPECTED_PROPER_SUPPORTS = 18
EXPECTED_EMPTY_SUPPORTS = 32
EXPECTED_HOSTILE_COUNTS = (115, 50, 57, 8)
EXPECTED_ARITY3_CERTIFICATE = ((1, 4, 5), (1, 3, 1, 1), 11)
EXPECTED_ARITY3_PROPER_SUPPORTS = 21

EXPECTED_SELECTED_PACKET_SHA256 = (
    "b8c63cd3f0b78cd83048158873dfdf8078f26ebd3ac6bab6c643b7292f7ae27d"
)
EXPECTED_SCREEN_SHA256 = (
    "ea3f2619819dba921e3d874d1592d1e07b9792c14bb5d42c5b29bb7f991387ce"
)
EXPECTED_STATUS_SHA256 = (
    "7375be6547b92ef1211f98e483ae0c8b08a0c68efb78f4e914c18cae1f63ff4d"
)
EXPECTED_FULL_CIRCUIT_SHA256 = (
    "1027a04c9718f649ff4d017e2e31c93f8197ceace52381f5427a91065de505bb"
)
EXPECTED_MINIMALITY_SHA256 = (
    "4632c5274c29f8b3304e30ada2c2c948325a1910e088048d20a95335ed9a3138"
)
EXPECTED_EMPTY_SUPPORT_SHA256 = (
    "f55ea43849bce0d13249bd77ef3619dfd5b6fb282f1ef6eedffdfbc8bc7275c2"
)
EXPECTED_HOSTILE_SHA256 = (
    "afc67bd8bd482be01bf2a3b1f9e83e44b6fa244d54b42fd856eba42ef3e4b823"
)
EXPECTED_ARITY3_MINIMALITY_SHA256 = (
    "7607d0f94fb1338fcfbd0b84730a7ebacdaf08d6921360d193c34144eb731cae"
)
EXPECTED_SEMANTIC_SHA256 = (
    "c201276eb84f71806a7eb683e42d723b930f7cd0f79862e8d6b1e0ad07d37dd9"
)

LEDGER_BEFORE = 373_157
LEDGER_AFTER = 373_153
WALL_BEFORE = 353
WALL_AFTER = 349
FAMILIES_BEFORE = 31
FAMILIES_AFTER = 29


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


def affine_rows():
    return tuple(
        (1, *tuple((pattern >> coordinate) & 1 for coordinate in range(4)))
        for pattern in range(16)
    )


def affine_value(coefficients, pattern):
    return coefficients[0] + sum(
        coefficients[coordinate + 1] * ((pattern >> coordinate) & 1)
        for coordinate in range(4)
    )


def optimal_affine_majorant(multicover, capacities, selected, q, marginals):
    """Lexicographically canonical exact affine majorant of a layer sum."""
    values = multicover.layer_values(capacities, selected)
    rows = affine_rows()
    candidates = {}
    for tight_seed in itertools.combinations(range(16), 5):
        coefficients = multicover.solve_square(
            tuple(rows[pattern] for pattern in tight_seed),
            tuple(values[pattern] for pattern in tight_seed),
        )
        if coefficients is None:
            continue
        if not all(
            affine_value(coefficients, pattern) >= values[pattern]
            for pattern in range(16)
        ):
            continue
        cost = coefficients[0] * q + sum(
            coefficients[coordinate + 1] * marginals[coordinate]
            for coordinate in range(4)
        )
        candidates[coefficients] = cost
    require(candidates, (selected, "no affine-majorant vertex"))
    cost, coefficients = min(
        (cost, coefficients) for coefficients, cost in candidates.items()
    )
    tight = tuple(
        pattern for pattern in range(16)
        if affine_value(coefficients, pattern) == values[pattern]
    )
    return coefficients, cost, tight


def primitive_affine_form(selected, coefficients):
    denominator = lcm(*(value.denominator for value in coefficients))
    raw = (
        *([denominator] * len(selected)),
        *(int(value * denominator) for value in coefficients),
    )
    divisor = gcd(*(abs(value) for value in raw if value))
    primitive = tuple(value // divisor for value in raw)
    tail_scale = primitive[0]
    require(
        all(value == tail_scale for value in primitive[:len(selected)]),
        primitive,
    )
    return primitive


def first_affine_pair_circuit(
    multicover, q, marginals, capacities, histogram
):
    thresholds, demands, _good = multicover.tail_data(capacities, histogram)
    for selected in itertools.combinations(thresholds, 2):
        table = multicover.exact_feasible_table(
            q, marginals, capacities, histogram, selected
        )
        if table is not None:
            continue
        coefficients, cost, tight = optimal_affine_majorant(
            multicover, capacities, selected, q, marginals
        )
        demand = sum(demands[threshold] for threshold in selected)
        require(demand > cost, (selected, demand, cost))
        primitive = primitive_affine_form(selected, coefficients)
        return selected, coefficients, demand - cost, demand, cost, tight, primitive
    return None


def exact_proper_support_tables(
    multicover, q, marginals, capacities, histogram, minimum
):
    thresholds, _demands, _good = multicover.tail_data(capacities, histogram)
    records = []
    for support_size in range(1, minimum):
        for selected in itertools.combinations(thresholds, support_size):
            table = multicover.exact_feasible_table(
                q, marginals, capacities, histogram, selected
            )
            require(table is not None, (selected, "infeasible proper support"))
            records.append((selected, table))
    return tuple(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    dependencies = (
        (THM3313, THM3313_SHA256, "id: THM-3313"),
        (THIRD_SOURCE, THIRD_SOURCE_SHA256, None),
        (THIRD_OUTPUT, THIRD_OUTPUT_SHA256,
         f"semantic_sha256={THIRD_SEMANTIC_SHA256}"),
        (MULTICOVER_SOURCE, MULTICOVER_SOURCE_SHA256, None),
        (MULTICOVER_OUTPUT, MULTICOVER_OUTPUT_SHA256,
         f"semantic_sha256={MULTICOVER_SEMANTIC_SHA256}"),
    )
    for path, digest, needle in dependencies:
        require(lf_sha(path) == digest, (path, "dependency changed"))
        if needle is not None:
            require(needle in path.read_text(encoding="utf-8"), (path, needle))

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require((assert_nodes, float_nodes) == (0, 0), "truth-gate syntax")

    third = load("fourth_prefix_third", THIRD_SOURCE)
    audit = load("fourth_prefix_audit", third.AUDIT_SOURCE)
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

    inherited_prefixes = []
    for prefix_number in range(1, 4):
        live = tuple(
            index for index, row in enumerate(rows)
            if row[3] and index not in closed
        )
        ranked = audit.ranked_families(rows, components, live)
        selected = audit.through_next_nonsingleton(ranked)
        indices = tuple(index for family in selected for index in family[4])
        inherited_prefixes.append((prefix_number, live, ranked, selected, indices))
        closed.update(indices)
    require(
        tuple(len(record[4]) for record in inherited_prefixes) == (5, 18, 4),
        inherited_prefixes,
    )
    require(
        inherited_prefixes[-1][4] == third.EXPECTED_SELECTED_INDICES,
        inherited_prefixes[-1][4],
    )

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
        returned_index, screen_row, direct, legacy = audit.screen_worker(
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

    status_records = []
    full_circuit_records = []
    proper_support_records = []
    modular_by_row = {}
    full_by_row = {}
    old_to_modular_counter = {}
    for index in selected_indices:
        eng = audit.status_engine(natural, f"fourth_selected_{index}")
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
        row_modular = Counter()
        row_full = Counter()
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
            modular_support, modular_certificates = (
                multicover.minimal_binary_multicovers(
                    capacities, marginals, histogram
                )
            )
            require(modular_support is not None and modular_certificates,
                    (index, divisors, "no modular certificate"))
            modular_selected, modular_weights, modular_gap = (
                modular_certificates[0]
            )
            require(
                multicover.coordinate_minimal(
                    capacities, modular_selected, modular_weights
                ),
                (index, divisors, "coordinate-nonminimal"),
            )
            old_branch = multicover.old_branch(
                audit, capacities, marginals, histogram
            )
            old_to_modular_counter.setdefault(
                old_branch, Counter()
            )[modular_support] += 1
            row_modular[modular_support] += 1

            if modular_support == 1:
                coefficients = (Fraction(0), *map(Fraction, modular_weights))
                full_selected = modular_selected
                full_gap = Fraction(modular_gap)
                thresholds, demands, _good = multicover.tail_data(
                    capacities, histogram
                )
                demand = sum(demands[threshold] for threshold in full_selected)
                cost = Fraction(demand) - full_gap
                tight = tuple(
                    pattern for pattern in range(16)
                    if affine_value(coefficients, pattern)
                    == multicover.layer_values(capacities, full_selected)[pattern]
                )
                primitive = primitive_affine_form(full_selected, coefficients)
            else:
                proper = exact_proper_support_tables(
                    multicover,
                    q,
                    marginals,
                    capacities,
                    histogram,
                    2,
                )
                proper_support_records.append(
                    (index, divisors, q, marginals, histogram,
                     tuple(capacities), proper)
                )
                affine = first_affine_pair_circuit(
                    multicover, q, marginals, capacities, histogram
                )
                require(affine is not None, (index, divisors, "no affine pair"))
                (
                    full_selected,
                    coefficients,
                    full_gap,
                    demand,
                    cost,
                    tight,
                    primitive,
                ) = affine
            full_support = len(full_selected)
            row_full[full_support] += 1
            require(
                all(
                    affine_value(coefficients, pattern)
                    >= multicover.layer_values(
                        capacities, full_selected
                    )[pattern]
                    for pattern in range(16)
                ),
                (index, divisors, "invalid affine majorant"),
            )
            require(Fraction(demand) - cost == full_gap and full_gap > 0,
                    (index, divisors, demand, cost, full_gap))

            status_records.append(
                (
                    index,
                    divisors,
                    q,
                    cofactor,
                    marginals,
                    histogram,
                    tuple(capacities),
                    old_branch,
                    modular_selected,
                    modular_weights,
                    modular_gap,
                    len(modular_certificates),
                    full_selected,
                    coefficients,
                    full_gap,
                    demand,
                    cost,
                    tight,
                    primitive,
                )
            )
            if full_support > 1:
                full_circuit_records.append(status_records[-1])
        modular_by_row[index] = tuple(sorted(row_modular.items()))
        full_by_row[index] = tuple(sorted(row_full.items()))

    status_records = tuple(status_records)
    full_circuit_records = tuple(full_circuit_records)
    proper_support_records = tuple(proper_support_records)
    status_sha = hashlib.sha256(repr(status_records).encode()).hexdigest()
    full_circuit_sha = hashlib.sha256(
        repr(full_circuit_records).encode()
    ).hexdigest()
    minimality_sha = hashlib.sha256(
        repr(proper_support_records).encode()
    ).hexdigest()
    if EXPECTED_STATUS_SHA256 is not None:
        require(status_sha == EXPECTED_STATUS_SHA256, status_sha)
    if EXPECTED_FULL_CIRCUIT_SHA256 is not None:
        require(full_circuit_sha == EXPECTED_FULL_CIRCUIT_SHA256,
                full_circuit_sha)
    if EXPECTED_MINIMALITY_SHA256 is not None:
        require(minimality_sha == EXPECTED_MINIMALITY_SHA256, minimality_sha)

    modular_layers = tuple(
        sorted(Counter(len(record[8]) for record in status_records).items())
    )
    full_layers = tuple(
        sorted(Counter(len(record[12]) for record in status_records).items())
    )
    old_to_modular = tuple(
        (branch, tuple(sorted(counter.items())))
        for branch, counter in sorted(old_to_modular_counter.items())
    )
    by_row = tuple(
        (
            index,
            counts[index],
            modular_by_row[index],
            full_by_row[index],
        )
        for index in selected_indices
    )
    require(len(status_records) == EXPECTED_TOTALS[2], len(status_records))
    require(modular_layers == EXPECTED_MODULAR_LAYER_COUNTS, modular_layers)
    require(full_layers == EXPECTED_FULL_LAYER_COUNTS, full_layers)
    require(old_to_modular == EXPECTED_OLD_TO_MODULAR, old_to_modular)
    require(by_row == EXPECTED_BY_ROW, by_row)

    instance_counts = Counter(
        (record[2], record[4], record[5], record[6])
        for record in status_records
    )
    multiplicities = tuple(sorted(Counter(instance_counts.values()).items()))
    require(len(instance_counts) == EXPECTED_UNIQUE_INSTANCES, len(instance_counts))
    require(multiplicities == EXPECTED_INSTANCE_MULTIPLICITIES, multiplicities)
    empty_support_records = []
    for q, marginals, histogram, capacities in sorted(
        instance_counts, key=repr
    ):
        table = multicover.exact_feasible_table(
            q, marginals, capacities, histogram, ()
        )
        require(table is not None, (q, marginals, "empty-tail infeasible"))
        empty_support_records.append(
            (q, marginals, histogram, capacities, table)
        )
    empty_support_records = tuple(empty_support_records)
    require(len(empty_support_records) == EXPECTED_EMPTY_SUPPORTS,
            len(empty_support_records))
    empty_support_sha = hashlib.sha256(
        repr(empty_support_records).encode()
    ).hexdigest()
    if EXPECTED_EMPTY_SUPPORT_SHA256 is not None:
        require(empty_support_sha == EXPECTED_EMPTY_SUPPORT_SHA256,
                empty_support_sha)
    require(
        tuple((record[0], record[1]) for record in full_circuit_records)
        == EXPECTED_MULTI_ADDRESSES,
        tuple((record[0], record[1]) for record in full_circuit_records),
    )
    require(
        tuple((record[12], record[13], record[14])
              for record in full_circuit_records)
        == EXPECTED_FULL_CERTIFICATES,
        tuple((record[12], record[13], record[14])
              for record in full_circuit_records),
    )
    require(
        sum(len(record[-1]) for record in proper_support_records)
        == EXPECTED_PROPER_SUPPORTS,
        proper_support_records,
    )

    # Hostile negative control: row 64 keeps eight screen residuals.  The first
    # has an explicit exact common table and no modular contradiction.
    hostile_index = multicover.HOSTILE_INDEX
    eng = audit.status_engine(natural, "fourth_hostile64")
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
    hostile_support, hostile_certificates = (
        multicover.minimal_binary_multicovers(
            hostile_capacities, hostile_marginals, hostile_histogram
        )
    )
    require(hostile_support is None and not hostile_certificates,
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

    # Hostile positive control: THM-3308's row-138 packet genuinely needs
    # support three.  Reconstruct every singleton and pair common table.
    arity3_index = 138
    arity3_divisors = (18, 49, 196, 882)
    eng = audit.status_engine(natural, "fourth_arity3")
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    stream = eng.ray.Stream(rows[arity3_index][0])
    _trials, states, _checks, _signs = eng.ray.ray_quotient_states(stream)
    _crude, status, _residual = eng.exact_common_status_screen(stream, states)
    q, _cofactor, marginals, _capacity_set, histogram, _raw_dual = (
        status[arity3_divisors]
    )
    _rebuilt, capacities = eng.ray.local.hunter_status_data(
        lcm(*arity3_divisors), arity3_divisors, q
    )
    arity3_support, arity3_certificates = (
        multicover.minimal_binary_multicovers(capacities, marginals, histogram)
    )
    require(arity3_support == 3 and arity3_certificates, arity3_certificates)
    require(arity3_certificates[0] == EXPECTED_ARITY3_CERTIFICATE,
            arity3_certificates[0])
    arity3_proper = exact_proper_support_tables(
        multicover, q, marginals, capacities, histogram, 3
    )
    require(len(arity3_proper) == EXPECTED_ARITY3_PROPER_SUPPORTS,
            len(arity3_proper))
    arity3_minimality_sha = hashlib.sha256(
        repr(arity3_proper).encode()
    ).hexdigest()
    if EXPECTED_ARITY3_MINIMALITY_SHA256 is not None:
        require(
            arity3_minimality_sha == EXPECTED_ARITY3_MINIMALITY_SHA256,
            arity3_minimality_sha,
        )

    remaining = tuple(index for index in live if index not in set(selected_indices))
    remaining_ranked = audit.ranked_families(rows, components, remaining)
    require((len(remaining), len(remaining_ranked)) == (WALL_AFTER, FAMILIES_AFTER),
            (len(remaining), len(remaining_ranked)))
    require(LEDGER_BEFORE - len(selected_indices) == LEDGER_AFTER, "ledger")
    require(WALL_BEFORE - len(selected_indices) == WALL_AFTER, "wall")
    require(FAMILIES_BEFORE - len(selected_families) == FAMILIES_AFTER,
            "families")

    semantic_packet = (
        THM3313_SHA256,
        THIRD_SOURCE_SHA256,
        THIRD_OUTPUT_SHA256,
        THIRD_SEMANTIC_SHA256,
        MULTICOVER_SOURCE_SHA256,
        MULTICOVER_OUTPUT_SHA256,
        MULTICOVER_SEMANTIC_SHA256,
        census,
        tuple((record[0], record[4]) for record in inherited_prefixes),
        ranked,
        selected_packet,
        canonical,
        by_row,
        status_records,
        full_circuit_records,
        proper_support_records,
        empty_support_records,
        multiplicities,
        hostile_packet,
        arity3_certificates[0],
        arity3_proper,
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
    circuit_summaries = tuple(
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
            record[12],
            record[13],
            record[14],
            record[15],
            record[16],
            record[17],
            record[18],
        )
        for record in full_circuit_records
    )
    lines = [
        "LRC14 projected-k3 z216 fourth intrinsic-ruler prefix affine-multicover closure",
        (
            f"dependency=THM3313:{THM3313_SHA256};third_prefix:"
            f"{THIRD_SOURCE_SHA256}/{THIRD_OUTPUT_SHA256}/"
            f"{THIRD_SEMANTIC_SHA256};multicover_compiler:"
            f"{MULTICOVER_SOURCE_SHA256}/{MULTICOVER_OUTPUT_SHA256}/"
            f"{MULTICOVER_SEMANTIC_SHA256}"
        ),
        (
            f"universe=z1:{LEVEL};atlas_rows:{census[0]};wall:{census[1]};"
            f"order:{census[2]};prior_prefix_rows:"
            f"{sum(len(record[4]) for record in inherited_prefixes)};"
            f"live_wall:{census[3]};live_families:{census[4]}"
        ),
        (
            f"queue_head:{family_head(ranked, 3)};canonical_prior_"
            "reconstruction:PASS;complete_family_reconstruction:PASS;"
            "strict_cost_boundary:PASS"
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
                f"integral_modular_taxonomy=layer_counts:{modular_layers};"
                f"old_to_modular:{old_to_modular};status_sha256:{status_sha}"
            ),
            (
                f"full_common_table_taxonomy=layer_counts:{full_layers};"
                f"unique_instances:{len(instance_counts)};instance_"
                f"multiplicities:{multiplicities};circuit_sha256:"
                f"{full_circuit_sha}"
            ),
            f"full_circuits:{circuit_summaries}",
            (
                f"proper_support_minimality=addressed_tables:"
                f"{sum(len(record[-1]) for record in proper_support_records)};"
                "all_exactly_feasible:PASS;minimality_sha256:"
                f"{minimality_sha}"
            ),
            (
                f"empty_tail_controls=unique_instances:"
                f"{len(empty_support_records)};all_marginal_systems_exactly_"
                f"feasible:PASS;empty_support_sha256:{empty_support_sha}"
            ),
            (
                "repair=the_zero_constant_nonnegative_integral_modular_"
                "template_reports_two_support_three_packets_but_exact_"
                "proper_support_audit_exposes_support_two_half_affine_"
                "majorants;integral_template_arity_is_not_full_dual_arity"
            ),
            (
                "affine_mechanism=sum_selected_tail_indicators_is_pointwise_"
                "bounded_by_an_exact_affine_Boolean_form;one_common_status_"
                "table_converts_its_constant_to_a_times_q_and_coordinates_"
                "to_weighted_marginals"
            ),
            (
                f"hostile=row64_counts:{hostile_counts};all_57_status_kills_"
                "one_layer;first_residual_exact_common_table:PASS;no_"
                f"multicover_contradiction:PASS;hostile_sha256:{hostile_sha}"
            ),
            (
                f"arity3_control=row138_divisors:{arity3_divisors};"
                f"certificate:{arity3_certificates[0]};proper_support_"
                f"tables:{len(arity3_proper)};all_exactly_feasible:PASS;"
                f"minimality_sha256:{arity3_minimality_sha}"
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
