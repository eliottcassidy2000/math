#!/usr/bin/env python3
"""Exact post-THM-3313 complete prefix with weighted tail multicovers.

THM-3313 leaves 353 projected-k3, z1=216 wall rows in 31 complete
``(gcd(216,L),L)`` families.  This companion reconstructs that state rather
than importing four addresses, selects the next complete cost prefix through
its first nonsingleton family, and runs the exact necessary upper screens.

For status explanations it extends THM-3308's binary threshold-set compiler
to a bounded threshold *multiset* language

    sum_t alpha_t 1[c(P)>=t] <= sum_i w_i 1[i in P],

where positive integral tail multiplicities have total at most three.  This
distinguishes tail-support size from total layer multiplicity.  Exact feasible
common tables audit minimum support in the full rational dual cone.

The result remains inside a necessary projected atlas.  It restores no
physical speed entry, endpoint, owner, phase, current, arbitrary k<=1, rung,
or LRC(14) conclusion.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import itertools
import multiprocessing as mp
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
    "lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_"
    "multicover_closure_scout_20260803.py"
)
THIRD_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_third_intrinsic_ruler_cost_prefix_"
    "multicover_closure_scout_20260803.out"
)

DEPENDENCIES = (
    (
        AUDIT_SOURCE,
        "d3581bc937597270558a77c0a02398a4cfd368306f5a302be01240cdb4cef7c3",
        None,
    ),
    (
        AUDIT_OUTPUT,
        "e035531d0af15998a8a07d083042ff2486c80c7fdd729d902922179409d1725a",
        "4bb4eca9949eb4f07b618f9e49bbb5de63492364a29a4a5abfc84a6acf105e11",
    ),
    (
        MULTICOVER_SOURCE,
        "e7462eabab773133688079141708d2377742e2d6f9a9480089666956f472020c",
        None,
    ),
    (
        MULTICOVER_OUTPUT,
        "f83d7ca1a3e2557330ec7ed655dc756cef21d6abdb52326d3623563c2abc407d",
        "bd7ac011d960afe7c4566ac9c17e9df3d8c2cf710e6023a715cb6a6d1ce1408d",
    ),
    (
        THIRD_SOURCE,
        "34710b907f3274057d6cb60a0fd7bd72ae2c3879a90b830a3e376d30e5198646",
        None,
    ),
    (
        THIRD_OUTPUT,
        "c664acf4e02b40835dfdc76ff135d17332e3d0e57b489e61b19c0f2f7b5998f5",
        "4d8547f53b287d82bb725cd6b3ed78c390ba03a37a3dcb3ab76ac98389e045a3",
    ),
)

LEVEL = 216
MAX_LAYER_MASS = 3
EXPECTED_PREFIX_LENGTHS = (5, 18, 4)
EXPECTED_CURRENT_CENSUS = (480, 447, 33, 353, 31)
EXPECTED_SELECTED_INDICES = (141, 133, 219, 359)
EXPECTED_SELECTED_FAMILIES = (
    (
        2_293_200, 1, 24, 76_440, (141,),
        ((141, (1, 4, 6, 10, 13, 14), 30, 2_293_200),),
    ),
    (
        2_629_536, 3, 24, 30_576, (133, 219, 359),
        (
            (133, (1, 4, 6, 8, 13, 14), 26, 794_976),
            (219, (1, 6, 8, 12, 13, 14), 30, 917_280),
            (359, (2, 6, 8, 12, 13, 14), 30, 917_280),
        ),
    ),
)
EXPECTED_BY_ROW = (
    (141, (10, 8, 2, 0), ((1, 2),)),
    (133, (10, 0, 10, 0), ((1, 7), (2, 3))),
    (219, (189, 147, 42, 0), ((1, 42),)),
    (359, (21, 17, 4, 0), ((1, 4),)),
)
EXPECTED_SCREEN_TOTALS = (230, 172, 58, 0)
EXPECTED_SUPPORT_COUNTS = ((1, 55), (2, 3))
EXPECTED_LAYER_MASS_COUNTS = ((1, 55), (2, 1), (3, 2))
EXPECTED_BRANCH_COUNTS = (
    (("coordinate_union", 1), 55),
    (("weighted_core", 2), 3),
)
EXPECTED_MULTI_CERTIFICATES = (
    (
        133, (84, 588, 728, 1274), (1, 5), (1, 2),
        (1, 1, 1, 3), 364,
    ),
    (
        133, (84, 728, 1092, 1274), (1, 5), (1, 2),
        (1, 1, 1, 3), 364,
    ),
    (
        133, (84, 728, 1274, 2184), (2, 5), (1, 1),
        (1, 1, 2, 0), 510,
    ),
)
EXPECTED_MINIMALITY_TABLES = 18
EXPECTED_HOSTILE_COUNTS = (115, 50, 57, 8)
EXPECTED_ARITY3 = (
    (1, 4, 5), (1, 1, 1), (1, 3, 1, 1), 11,
)
EXPECTED_LEDGER = (373_157, 373_153, 353, 349, 31, 29)
EXPECTED_NEXT_FAMILY = (
    3_400_992, 19, 72, 7056,
    (21, 49, 62, 92, 129, 148, 202, 203, 212, 254, 262, 290, 303,
     343, 344, 352, 387, 428, 465),
)

EXPECTED_SELECTED_PACKET_SHA256 = (
    "b8c63cd3f0b78cd83048158873dfdf8078f26ebd3ac6bab6c643b7292f7ae27d"
)
EXPECTED_SCREEN_SHA256 = (
    "ea3f2619819dba921e3d874d1592d1e07b9792c14bb5d42c5b29bb7f991387ce"
)
EXPECTED_CERTIFICATE_SHA256 = (
    "2bcd5f8550db3ae1f4e44f4781edb730d8bb078d0602aa54fba2b0a11f0d3e4f"
)
EXPECTED_MINIMALITY_SHA256 = (
    "b1d0840a625afa68f56a20616af718a853476f9852658fd048eafc86890dcfdf"
)
EXPECTED_CONTROL_SHA256 = (
    "281edd020b801a36971d77659e19daaa504c71f9fdd567b41c2c3af3b18414be"
)
EXPECTED_SEMANTIC_SHA256 = (
    "666458d6122a66a7615a2855a1b0dcf18fb778e1451b44c064bab04726d81b05"
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


def weighted_layer_values(capacities, thresholds, multiplicities):
    return tuple(
        sum(
            multiplicity * int(capacities[pattern] >= threshold)
            for threshold, multiplicity in zip(thresholds, multiplicities)
        )
        for pattern in range(16)
    )


def dominates(values, weights):
    return all(
        values[pattern]
        <= sum(
            weights[coordinate] * ((pattern >> coordinate) & 1)
            for coordinate in range(4)
        )
        for pattern in range(16)
    )


def minimal_weighted_multicovers(multicover, capacities, marginals, histogram):
    """Least support, then least total multiplicity, through mass three."""
    thresholds, demands, _good_rows = multicover.tail_data(
        capacities, histogram
    )
    for support_size in range(1, len(thresholds) + 1):
        for total_mass in range(support_size, MAX_LAYER_MASS + 1):
            certificates = []
            for selected in itertools.combinations(thresholds, support_size):
                for multiplicities in itertools.product(
                    range(1, total_mass + 1), repeat=support_size
                ):
                    if sum(multiplicities) != total_mass:
                        continue
                    values = weighted_layer_values(
                        capacities, selected, multiplicities
                    )
                    demand = sum(
                        multiplicity * demands[threshold]
                        for threshold, multiplicity
                        in zip(selected, multiplicities)
                    )
                    for weights in itertools.product(
                        range(total_mass + 1), repeat=4
                    ):
                        supply = sum(
                            weight * marginal
                            for weight, marginal in zip(weights, marginals)
                        )
                        if demand <= supply or not dominates(values, weights):
                            continue
                        primitive = tuple(
                            value for value in (*multiplicities, *weights)
                            if value
                        )
                        if gcd(*primitive) != 1:
                            continue
                        certificates.append(
                            (selected, multiplicities, weights, demand - supply)
                        )
            if certificates:
                return support_size, total_mass, tuple(sorted(certificates))
    return None, None, ()


def tight_patterns(capacities, thresholds, multiplicities, weights):
    values = weighted_layer_values(capacities, thresholds, multiplicities)
    return tuple(
        pattern for pattern in range(16)
        if values[pattern]
        == sum(
            weights[coordinate] * ((pattern >> coordinate) & 1)
            for coordinate in range(4)
        )
    )


def selected_worker(item):
    index, task = item
    audit = load(f"fourth_audit_{index}", AUDIT_SOURCE)
    multicover = load(f"fourth_multicover_{index}", MULTICOVER_SOURCE)
    natural = audit.load(f"fourth_natural_{index}", audit.NATURAL_SOURCE)
    rows, _components = natural.atlas_rows()

    returned_index, screen_row, direct, legacy = audit.screen_worker(item)
    require(returned_index == index, (index, returned_index))
    source = rows[index]
    divisor_gcd = gcd(LEVEL, source[1])
    require(
        screen_row[:6]
        == (
            LEVEL, source[0], source[1], source[2],
            source[1] // divisor_gcd, source[3],
        ),
        (index, screen_row[:6]),
    )
    require(screen_row[16] == screen_row[11], (index, "status verification"))
    require(direct + legacy == screen_row[11], (index, "Farkas count"))
    counts = tuple(screen_row[position] for position in (9, 10, 11, 12))
    require(not screen_row[13], (index, "screen residual"))

    eng = audit.status_engine(natural, f"fourth_detail_{index}")
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    stream = eng.ray.Stream(source[0])
    require(
        (stream.L, stream.high_floor) == (source[1], source[2]),
        (index, "atlas"),
    )
    _trials, states, _checks, _signs = eng.ray.ray_quotient_states(stream)
    crude, status, residual = eng.exact_common_status_screen(stream, states)
    require(
        (len(states), len(crude), len(status), len(residual)) == counts,
        (index, "screen reconstruction"),
    )

    records = []
    for divisors, witness in sorted(status.items()):
        q, cofactor, marginals, capacity_set, histogram, raw_dual = witness
        divisor_lcm = lcm(*divisors)
        require(divisor_lcm == q * cofactor, (index, divisors, "cofactor"))
        rebuilt, capacities = eng.ray.local.hunter_status_data(
            divisor_lcm, divisors, q
        )
        require(rebuilt == marginals, (index, divisors, "marginals"))
        require(
            tuple(sorted(set(capacities))) == capacity_set,
            (index, divisors, "capacities"),
        )
        eng.independent_farkas_check(
            q, marginals, capacities, histogram, raw_dual
        )
        support_size, total_mass, certificates = minimal_weighted_multicovers(
            multicover, capacities, marginals, histogram
        )
        require(support_size is not None and certificates,
                (index, divisors, "weighted compiler"))
        selected, multiplicities, weights, gap = certificates[0]
        branch = multicover.old_branch(
            audit, capacities, marginals, histogram
        )
        records.append(
            (
                index, divisors, q, cofactor, marginals, histogram,
                tuple(capacities), branch, selected, multiplicities, weights,
                gap,
                tight_patterns(
                    capacities, selected, multiplicities, weights
                ),
                total_mass, len(certificates),
            )
        )
    return (
        index, counts, tuple(screen_row[:19]), direct, legacy, tuple(records)
    )


def exact_support_minimality(multicover, multi_records):
    records = []
    for record in multi_records:
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
                    (record[0], record[1], threshold, "singleton"))
            records.append((record[0], record[1], threshold, table))
    return tuple(records)


def residual_control(audit, multicover, natural, rows):
    index = multicover.HOSTILE_INDEX
    eng = audit.status_engine(natural, "fourth_hostile64")
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    stream = eng.ray.Stream(rows[index][0])
    _trials, states, _checks, _signs = eng.ray.ray_quotient_states(stream)
    crude, status, residual = eng.exact_common_status_screen(stream, states)
    counts = (len(states), len(crude), len(status), len(residual))
    require(counts == EXPECTED_HOSTILE_COUNTS, counts)
    divisors = multicover.HOSTILE_DIVISORS
    require(divisors in residual, residual)
    divisor_lcm = lcm(*divisors)
    arcs = eng.ray.fibre.projected_support_arcs(divisor_lcm, stream.ranges)
    marginals, capacities = eng.ray.local.hunter_status_data(
        divisor_lcm, divisors, multicover.HOSTILE_Q
    )
    histogram = eng.ray.fibre.residue_load_histogram(
        arcs, multicover.HOSTILE_Q
    )
    thresholds, demands, good_rows = multicover.tail_data(
        capacities, histogram
    )
    table = tuple(map(Fraction, multicover.HOSTILE_TABLE))
    require(
        multicover.verify_table(
            multicover.HOSTILE_Q, marginals, good_rows, demands,
            thresholds, table,
        ),
        "hostile common table",
    )
    support, mass, certificates = minimal_weighted_multicovers(
        multicover, capacities, marginals, histogram
    )
    require(support is None and mass is None and not certificates,
            "hostile falsely killed")
    return (
        counts, divisors, multicover.HOSTILE_Q, marginals,
        tuple(capacities), histogram, table,
        hashlib.sha256(repr(residual).encode()).hexdigest(),
    )


def arity_three_control(audit, multicover, natural, rows):
    index = 138
    divisors = (18, 49, 196, 882)
    eng = audit.status_engine(natural, "fourth_arity3")
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    stream = eng.ray.Stream(rows[index][0])
    _trials, states, _checks, _signs = eng.ray.ray_quotient_states(stream)
    _crude, status, _residual = eng.exact_common_status_screen(stream, states)
    require(divisors in status, (index, divisors))
    q, _cofactor, marginals, _capacity_set, histogram, raw_dual = status[divisors]
    _rebuilt, capacities = eng.ray.local.hunter_status_data(
        lcm(*divisors), divisors, q
    )
    eng.independent_farkas_check(
        q, marginals, capacities, histogram, raw_dual
    )
    support, mass, certificates = minimal_weighted_multicovers(
        multicover, capacities, marginals, histogram
    )
    require((support, mass) == (3, 3), (support, mass))
    selected, multiplicities, weights, gap = certificates[0]
    require(
        (selected, multiplicities, weights, gap) == EXPECTED_ARITY3,
        certificates[0],
    )
    thresholds, _demands, _good = multicover.tail_data(
        capacities, histogram
    )
    proper = []
    for size in (1, 2):
        for subset in itertools.combinations(thresholds, size):
            table = multicover.exact_feasible_table(
                q, marginals, capacities, histogram, subset
            )
            require(table is not None, (index, divisors, subset))
            proper.append((subset, table))
    return (
        index, divisors, q, marginals, histogram, tuple(capacities),
        certificates[0], tuple(proper),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    for path, digest, semantic in DEPENDENCIES:
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

    audit = load("fourth_main_audit", AUDIT_SOURCE)
    multicover = load("fourth_main_multicover", MULTICOVER_SOURCE)
    natural = audit.load("fourth_main_natural", audit.NATURAL_SOURCE)
    order = audit.load("fourth_main_order", audit.ORDER_SOURCE)
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
    prefix_indices = []
    for _stage in range(3):
        live = tuple(
            index for index, row in enumerate(rows)
            if row[3] and index not in closed
        )
        selected = audit.through_next_nonsingleton(
            audit.ranked_families(rows, components, live)
        )
        indices = tuple(
            index for family in selected for index in family[4]
        )
        prefix_indices.append(indices)
        closed |= set(indices)
    prefix_indices = tuple(prefix_indices)
    require(tuple(map(len, prefix_indices)) == EXPECTED_PREFIX_LENGTHS,
            tuple(map(len, prefix_indices)))

    current_live = tuple(
        index for index, row in enumerate(rows)
        if row[3] and index not in closed
    )
    current_ranked = audit.ranked_families(rows, components, current_live)
    census = (
        len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows),
        len(current_live), len(current_ranked),
    )
    require(census == EXPECTED_CURRENT_CENSUS, census)
    selected_families = audit.through_next_nonsingleton(current_ranked)
    selected_indices = tuple(
        index for family in selected_families for index in family[4]
    )
    require(tuple(selected_families) == EXPECTED_SELECTED_FAMILIES,
            selected_families)
    require(selected_indices == EXPECTED_SELECTED_INDICES, selected_indices)
    require(
        selected_families[-1][0] < current_ranked[len(selected_families)][0],
        "cost tie crosses prefix",
    )
    selected_packet = tuple(
        (
            index, rows[index][0], gcd(LEVEL, rows[index][1]),
            rows[index][1], components[index],
            rows[index][1] * components[index], rows[index][3],
        )
        for index in selected_indices
    )
    selected_packet_sha = hashlib.sha256(
        repr(selected_packet).encode()
    ).hexdigest()
    if EXPECTED_SELECTED_PACKET_SHA256 is not None:
        require(selected_packet_sha == EXPECTED_SELECTED_PACKET_SHA256,
                selected_packet_sha)

    tasks = tuple(
        (index, (LEVEL, *rows[index])) for index in selected_indices
    )
    if args.processes == 1:
        detailed = tuple(selected_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(
            min(args.processes, len(tasks))
        ) as pool:
            detailed = tuple(pool.map(selected_worker, tasks, chunksize=1))
    require(tuple(record[0] for record in detailed) == selected_indices,
            "worker order")
    canonical = tuple((record[0], record[2]) for record in detailed)
    screen_sha = hashlib.sha256(repr(canonical).encode()).hexdigest()
    if EXPECTED_SCREEN_SHA256 is not None:
        require(screen_sha == EXPECTED_SCREEN_SHA256, screen_sha)
    screen_totals = sum_counts(tuple(record[1] for record in detailed))
    require(screen_totals == EXPECTED_SCREEN_TOTALS, screen_totals)
    require(sum(record[3] for record in detailed) == 0, "direct Farkas drift")
    require(sum(record[4] for record in detailed) == screen_totals[2],
            "legacy Farkas count")

    certificate_records = tuple(
        certificate for record in detailed for certificate in record[5]
    )
    require(len(certificate_records) == screen_totals[2],
            len(certificate_records))
    certificate_sha = hashlib.sha256(
        repr(certificate_records).encode()
    ).hexdigest()
    if EXPECTED_CERTIFICATE_SHA256 is not None:
        require(certificate_sha == EXPECTED_CERTIFICATE_SHA256,
                certificate_sha)
    support_counts = tuple(sorted(Counter(
        len(record[8]) for record in certificate_records
    ).items()))
    layer_mass_counts = tuple(sorted(Counter(
        sum(record[9]) for record in certificate_records
    ).items()))
    branch_counts = tuple(sorted(Counter(
        (record[7], len(record[8])) for record in certificate_records
    ).items()))
    require(support_counts == EXPECTED_SUPPORT_COUNTS, support_counts)
    require(layer_mass_counts == EXPECTED_LAYER_MASS_COUNTS, layer_mass_counts)
    require(branch_counts == EXPECTED_BRANCH_COUNTS, branch_counts)
    by_row = tuple(
        (
            record[0], record[1], tuple(sorted(Counter(
                len(certificate[8]) for certificate in record[5]
            ).items())),
        )
        for record in detailed
    )
    require(by_row == EXPECTED_BY_ROW, by_row)

    multi_records = tuple(
        record for record in certificate_records if len(record[8]) > 1
    )
    multi_summary = tuple(
        (
            record[0], record[1], record[8], record[9], record[10],
            record[11],
        )
        for record in multi_records
    )
    require(multi_summary == EXPECTED_MULTI_CERTIFICATES, multi_summary)
    minimality = exact_support_minimality(multicover, multi_records)
    require(len(minimality) == EXPECTED_MINIMALITY_TABLES, len(minimality))
    minimality_sha = hashlib.sha256(repr(minimality).encode()).hexdigest()
    if EXPECTED_MINIMALITY_SHA256 is not None:
        require(minimality_sha == EXPECTED_MINIMALITY_SHA256,
                minimality_sha)

    hostile = residual_control(audit, multicover, natural, rows)
    arity3 = arity_three_control(audit, multicover, natural, rows)
    controls = (hostile, arity3)
    control_sha = hashlib.sha256(repr(controls).encode()).hexdigest()
    if EXPECTED_CONTROL_SHA256 is not None:
        require(control_sha == EXPECTED_CONTROL_SHA256, control_sha)

    remaining = tuple(
        index for index in current_live if index not in set(selected_indices)
    )
    remaining_ranked = audit.ranked_families(rows, components, remaining)
    ledger_before, ledger_after, wall_before, wall_after, families_before, families_after = (
        EXPECTED_LEDGER
    )
    require(
        (
            ledger_before - len(selected_indices), wall_before - len(selected_indices),
            families_before - len(selected_families), len(remaining),
            len(remaining_ranked),
        )
        == (ledger_after, wall_after, families_after, wall_after, families_after),
        "ledger consequence",
    )
    require(current_ranked[len(selected_families)][:5] == EXPECTED_NEXT_FAMILY,
            current_ranked[:3])
    require(remaining_ranked[0][:5] == EXPECTED_NEXT_FAMILY,
            remaining_ranked[:2])

    semantic_packet = (
        tuple(digest for _path, digest, _semantic in DEPENDENCIES),
        hashlib.sha256(repr(rows).encode()).hexdigest(),
        prefix_indices,
        census,
        family_head(current_ranked, 3),
        selected_packet,
        canonical,
        certificate_records,
        minimality,
        controls,
        EXPECTED_LEDGER,
        family_head(remaining_ranked, 3),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected-k3 z216 fourth intrinsic-ruler prefix weighted-multicover probe",
        (
            "dependencies="
            f"two_prefix_audit:{DEPENDENCIES[0][1]}/{DEPENDENCIES[1][1]}/"
            f"{DEPENDENCIES[1][2]};THM3308:{DEPENDENCIES[2][1]}/"
            f"{DEPENDENCIES[3][1]}/{DEPENDENCIES[3][2]};THM3313:"
            f"{DEPENDENCIES[4][1]}/{DEPENDENCIES[5][1]}/{DEPENDENCIES[5][2]}"
        ),
        (
            f"universe=z1:{LEVEL};atlas_rows:{census[0]};wall:{census[1]};"
            f"order:{census[2]};prefix_lengths:{tuple(map(len, prefix_indices))};"
            f"current_live:{census[3]};current_families:{census[4]};"
            f"row_sha256:{hashlib.sha256(repr(rows).encode()).hexdigest()}"
        ),
        f"current_queue_head={family_head(current_ranked, 3)}",
        (
            f"selection=families:{len(selected_families)};rows:"
            f"{len(selected_indices)};indices:{selected_indices};cost:"
            f"{sum(family[0] for family in selected_families)};"
            f"packet_sha256:{selected_packet_sha};complete_family_and_"
            "strict_boundary:PASS"
        ),
        (
            f"screen=states:{screen_totals[0]};crude:{screen_totals[1]};"
            f"status:{screen_totals[2]};residual:{screen_totals[3]};"
            f"exact_Farkas_rechecked:{screen_totals[2]};"
            f"screen_sha256:{screen_sha}"
        ),
        f"by_row=index,screen,support_counts:{by_row}",
        (
            f"weighted_multicover=support_counts:{support_counts};"
            f"layer_mass_counts:{layer_mass_counts};branch_support:"
            f"{branch_counts};max_layer_mass:{MAX_LAYER_MASS};"
            f"certificate_sha256:{certificate_sha}"
        ),
        f"multi_support_certificates={multi_summary}",
        (
            f"full_dual_support_minimality=multi_addresses:{len(multi_records)};"
            f"singleton_tail_systems:{len(minimality)};all_exactly_feasible:PASS;"
            f"minimality_sha256:{minimality_sha}"
        ),
        (
            "mechanism=two_row133_states_need_tail_support_two_but_layer_mass_"
            "three_via_1_times_H1_plus_2_times_H5;the_third_needs_two_binary_"
            "layers_H2_plus_H5;tail_support_and_layer_multiplicity_are_distinct"
        ),
        (
            f"controls=row64:{hostile[0]};explicit_common_table_and_no_"
            f"weighted_multicover:PASS;row138_genuine_arity3:{arity3[6]};"
            f"all_{len(arity3[7])}_proper_supports_feasible:PASS;"
            f"control_sha256:{control_sha}"
        ),
        (
            f"consequence=ledger:{ledger_before}-{len(selected_indices)}="
            f"{ledger_after};z216_wall:{wall_before}-{len(selected_indices)}="
            f"{wall_after};families:{families_before}-{len(selected_families)}="
            f"{families_after};cap:216"
        ),
        (
            f"stopping_boundary=remaining_wall:{len(remaining)};families:"
            f"{len(remaining_ranked)};next_head:{family_head(remaining_ranked, 3)}"
        ),
        (
            "direction=empty_exact_necessary_upper_screen_excludes_each_"
            "selected_projected_wall_row;screen_feasibility_is_not_a_converse_"
            "and_intrinsic_cost_is_only_a_work_order"
        ),
        (
            "scope=projected_k3_z216_necessary_atlas_only;no_physical_speed_"
            "entry_endpoint_owner_phase_current_arbitrary_k_le_1_rung_or_"
            "LRC14_claim"
        ),
        f"truth_gates=assert_nodes:{assert_nodes};float_literals:{float_nodes}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    mp.freeze_support()
    main()
