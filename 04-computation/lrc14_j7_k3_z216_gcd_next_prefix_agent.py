#!/usr/bin/env python3
"""Exact next intrinsic-family prefix at projected-k3 z1=216.

The program reconstructs the current 357-row wall after the two independently
audited cost prefixes, derives the next whole-family prefix, runs every exact
necessary screen, and gives each status exclusion a solver-free event type.
Every solver-returned Farkas certificate is nevertheless rebuilt and checked
over the rationals.  Persisted semantics omit certificates and contradiction
magnitudes, as required by MISTAKE-331/333.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import multiprocessing as mp
import re
from collections import Counter
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

AUDIT_SOURCE_SHA256 = (
    "d3581bc937597270558a77c0a02398a4cfd368306f5a302be01240cdb4cef7c3"
)
AUDIT_OUTPUT_SHA256 = (
    "e035531d0af15998a8a07d083042ff2486c80c7fdd729d902922179409d1725a"
)
AUDIT_SEMANTIC_SHA256 = (
    "4bb4eca9949eb4f07b618f9e49bbb5de63492364a29a4a5abfc84a6acf105e11"
)
LEVEL = 216
LEDGER_BEFORE = 373_161
WALL_BEFORE = 357
EXPECTED_FAMILY_KEYS = ((72, 72_072), (24, 25_872))
EXPECTED_SEMANTIC_SHA256 = (
    "45a88091b7a1d6963807c81cde810ded1b3d0dc14f2fdc468dba6669295a5ef2"
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


def screen_worker(item):
    audit = load(f"gcd_next_prefix_screen_audit_{item[0]}", AUDIT_SOURCE)
    return audit.screen_worker(item)


def tail_demand(histogram, threshold):
    return sum(
        count for load_value, count in histogram if load_value >= threshold
    )


def two_threshold_overlap_cuts(capacities, marginals, histogram):
    """Find incidence cuts coupling two nested load tails.

    Suppose every status capable at the lower threshold is nonempty.  Outside
    a coordinate union S, suppose every status capable at the higher threshold
    has at least 1+e bits.  The total bit incidence M then obeys

        M >= T_low + e * max(0, T_high - sum_(j in S) m_j).

    A strict reversal is a solver-free contradiction.
    """
    cuts = []
    thresholds = tuple(
        threshold for threshold, _count in histogram if threshold > 0
    )
    total_incidence = sum(marginals)
    for low in thresholds:
        if capacities[0] >= low:
            continue
        low_demand = tail_demand(histogram, low)
        for high in thresholds:
            if high <= low:
                continue
            high_demand = tail_demand(histogram, high)
            for coordinate_mask in range(1, 16):
                outside = tuple(
                    pattern for pattern in range(16)
                    if capacities[pattern] >= high
                    and not (pattern & coordinate_mask)
                )
                if not outside:
                    continue
                extra_bits = min(pattern.bit_count() - 1 for pattern in outside)
                if extra_bits <= 0:
                    continue
                coordinates = tuple(
                    coordinate for coordinate in range(4)
                    if coordinate_mask & (1 << coordinate)
                )
                union_supply = sum(marginals[j] for j in coordinates)
                forced_outside = max(0, high_demand - union_supply)
                incidence_lower = low_demand + extra_bits * forced_outside
                if incidence_lower > total_incidence:
                    cuts.append(
                        (
                            low,
                            high,
                            coordinate_mask,
                            coordinates,
                            outside,
                            extra_bits,
                            low_demand,
                            high_demand,
                            union_supply,
                            forced_outside,
                            total_incidence,
                            incidence_lower,
                            incidence_lower - total_incidence,
                        )
                    )
    return tuple(cuts)


def negative_alpha_rejected(eng, q, marginals, capacities, histogram,
                            certificate):
    thresholds, alpha, z = certificate
    positive = next((i for i, value in enumerate(alpha) if value > 0), None)
    require(positive is not None, "Farkas certificate has no positive alpha")
    altered = list(alpha)
    altered[positive] = -altered[positive]
    try:
        eng.independent_farkas_check(
            q, marginals, capacities, histogram,
            (thresholds, tuple(altered), z),
        )
    except RuntimeError as error:
        return "negative Farkas alpha" in str(error)
    return False


def taxonomy_worker(item):
    index, task, expected_counts = item
    audit = load(f"gcd_next_prefix_taxonomy_audit_{index}", AUDIT_SOURCE)
    natural = audit.load(
        f"gcd_next_prefix_taxonomy_natural_{index}", audit.NATURAL_SOURCE
    )
    eng = audit.status_engine(natural, f"gcd_next_prefix_{index}")
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    _first, body, ruler, high_floor, wall = task
    stream = eng.ray.Stream(body)
    require(
        (stream.L, stream.high_floor) == (ruler, high_floor),
        (index, "atlas mismatch"),
    )
    _trials, states, _checks, _signs = eng.ray.ray_quotient_states(stream)
    crude, status, residual = eng.exact_common_status_screen(stream, states)
    counts = (len(states), len(crude), len(status), len(residual))
    require(counts == expected_counts, (index, counts, expected_counts))
    require(wall and not residual, (index, "selected screen has residual"))

    branches = Counter()
    certificate_counts = Counter()
    records = []
    overlap_records = []
    hostile_rejections = 0
    for divisors, witness in sorted(status.items()):
        q, cofactor, marginals, capacity_set, histogram, certificate = witness
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
            q, marginals, capacities, histogram, certificate
        )
        if hostile_rejections == 0:
            require(
                negative_alpha_rejected(
                    eng, q, marginals, capacities, histogram, certificate
                ),
                (index, "negative-alpha hostile accepted"),
            )
            hostile_rejections = 1

        unions = audit.union_cuts(capacities, marginals, histogram)
        reduced = () if unions else audit.zero_reduced_union_cuts(
            capacities, marginals, histogram
        )
        fans = () if (unions or reduced) else audit.two_fan_cuts(
            capacities, marginals, histogram
        )
        overlaps = () if (unions or reduced or fans) else (
            two_threshold_overlap_cuts(capacities, marginals, histogram)
        )
        if unions:
            branch = "coordinate_union"
            evidence = unions
        elif reduced:
            branch = "zero_reduced_union"
            evidence = reduced
        elif fans:
            branch = "two_fan"
            evidence = fans
        elif overlaps:
            branch = "two_threshold_overlap"
            evidence = overlaps
            overlap_records.append(
                (
                    index,
                    divisors,
                    q,
                    cofactor,
                    marginals,
                    capacity_set,
                    histogram,
                    tuple(capacities),
                    overlaps,
                )
            )
        else:
            branch = "weighted_core"
            evidence = ()
        branches[branch] += 1
        certificate_counts[branch] += len(evidence)
        records.append(
            (
                index,
                divisors,
                q,
                cofactor,
                marginals,
                capacity_set,
                histogram,
                branch,
                evidence,
            )
        )
    require(hostile_rejections == (1 if status else 0), index)
    summary = (
        branches["coordinate_union"],
        branches["zero_reduced_union"],
        branches["two_fan"],
        branches["two_threshold_overlap"],
        branches["weighted_core"],
    )
    certificate_summary = (
        certificate_counts["coordinate_union"],
        certificate_counts["zero_reduced_union"],
        certificate_counts["two_fan"],
        certificate_counts["two_threshold_overlap"],
    )
    require(sum(summary) == counts[2], (index, summary, counts))
    return (
        index,
        counts,
        summary,
        certificate_summary,
        tuple(records),
        tuple(overlap_records),
        hostile_rejections,
    )


def parse_control_indices(audit):
    low_text = audit.GCD8_OUTPUT.read_text(encoding="utf-8")
    match = re.search(r"residual=bodies:(.*);masks:[0-9]+;", low_text)
    require(match is not None, "low-gcd8 residual controls missing")
    low_controls = ast.literal_eval(match.group(1))
    require(low_controls and all(len(item) == 3 for item in low_controls),
            low_controls)

    costly_text = audit.COSTLY_AUDIT.read_text(encoding="utf-8")
    target_match = re.search(r"targets=(\([^\n]+\))", costly_text)
    require(target_match is not None, "costly targets missing")
    costly_targets = ast.literal_eval(target_match.group(1))
    costly_residuals = {
        int(index): int(residual)
        for index, residual in re.findall(
            r"row=([0-9]+);[^\n]*\n"
            r"screen=states:[0-9]+;crude:[0-9]+;status:[0-9]+;"
            r"residual:([0-9]+);",
            costly_text,
        )
    }
    require(set(costly_targets) == set(costly_residuals), costly_residuals)
    first_low = low_controls[0]
    expected = {first_low[0]: first_low[2], **costly_residuals}
    return tuple(expected), expected


def add_vectors(rows, length):
    return tuple(sum(row[position] for row in rows) for position in range(length))


def family_head(ranked, count):
    return tuple(
        (family[0], family[1], family[2], family[3], family[4])
        for family in ranked[:count]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=3)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    require(lf_sha(AUDIT_SOURCE) == AUDIT_SOURCE_SHA256, "audit source changed")
    require(lf_sha(AUDIT_OUTPUT) == AUDIT_OUTPUT_SHA256, "audit output changed")
    audit_output_text = AUDIT_OUTPUT.read_text(encoding="utf-8")
    require(
        f"semantic_sha256={AUDIT_SEMANTIC_SHA256}" in audit_output_text,
        "audit semantic changed",
    )
    require(
        "consequence=ledger:373184-5=373179-18=373161;"
        "z216_wall:380-5=375-18=357" in audit_output_text,
        "current ledger lineage missing",
    )

    audit = load("gcd_next_prefix_base_audit", AUDIT_SOURCE)
    for path, digest in audit.DEPENDENCIES:
        require(audit.lf_sha(path) == digest, (path, "dependency changed"))

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require((assert_nodes, float_nodes) == (0, 0), "truth-gate syntax")

    natural = audit.load("gcd_next_prefix_natural", audit.NATURAL_SOURCE)
    order = audit.load("gcd_next_prefix_order", audit.ORDER_SOURCE)
    rows, components = natural.atlas_rows()
    order_rows, order_components = order.atlas_rows()
    require((rows, components) == (order_rows, order_components),
            "atlas parser disagreement")

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
    prior_closed = gcd8 | gcd18 | order_indices | natural_indices
    require(len(prior_closed) == 100, len(prior_closed))

    first_live = tuple(
        index for index, row in enumerate(rows)
        if row[3] and index not in prior_closed
    )
    first_families = audit.through_next_nonsingleton(
        audit.ranked_families(rows, components, first_live)
    )
    first_indices = tuple(
        index for family in first_families for index in family[4]
    )

    second_closed = prior_closed | set(first_indices)
    second_live = tuple(
        index for index, row in enumerate(rows)
        if row[3] and index not in second_closed
    )
    second_families = audit.through_next_nonsingleton(
        audit.ranked_families(rows, components, second_live)
    )
    second_indices = tuple(
        index for family in second_families for index in family[4]
    )

    current_closed = second_closed | set(second_indices)
    current_live = tuple(
        index for index, row in enumerate(rows)
        if row[3] and index not in current_closed
    )
    current_ranked = audit.ranked_families(rows, components, current_live)
    require(
        (len(current_live), len(current_ranked)) == (WALL_BEFORE, 33),
        (len(current_live), len(current_ranked)),
    )
    selected_families = audit.through_next_nonsingleton(current_ranked)
    family_keys = tuple((family[2], family[3]) for family in selected_families)
    require(family_keys == EXPECTED_FAMILY_KEYS, family_keys)
    selected_indices = tuple(
        index for family in selected_families for index in family[4]
    )
    require(
        len(selected_indices) == 4
        and all(rows[index][3] for index in selected_indices),
        selected_indices,
    )
    selected_packet = tuple(
        (
            index,
            rows[index][0],
            gcd(LEVEL, rows[index][1]),
            rows[index][1],
            components[index],
            rows[index][1] * components[index],
        )
        for index in selected_indices
    )

    tasks = tuple(
        (index, (LEVEL, *rows[index])) for index in selected_indices
    )
    if args.processes == 1:
        screened = tuple(screen_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(
            min(args.processes, len(tasks))
        ) as pool:
            screened = tuple(pool.map(screen_worker, tasks, chunksize=1))
    require(tuple(item[0] for item in screened) == selected_indices,
            "screen order")

    canonical = {}
    counts = {}
    farkas_counts = {}
    for index, row, direct, legacy in screened:
        source = rows[index]
        divisor_gcd = gcd(LEVEL, source[1])
        require(
            row[:6] == (
                LEVEL,
                source[0],
                source[1],
                source[2],
                source[1] // divisor_gcd,
                source[3],
            ),
            (index, row[:6]),
        )
        require(row[16] == row[11], (index, "inherited exact audit"))
        require(direct + legacy == row[11], (index, "Farkas branch count"))
        require(row[12] == 0 and row[13] == (), (index, "screen residual"))
        canonical[index] = tuple(row[:19])
        counts[index] = tuple(row[position] for position in (9, 10, 11, 12))
        farkas_counts[index] = (direct, legacy)

    taxonomy_tasks = tuple(
        (index, (LEVEL, *rows[index]), counts[index])
        for index in selected_indices
    )
    if args.processes == 1:
        taxonomy = tuple(taxonomy_worker(task) for task in taxonomy_tasks)
    else:
        with mp.get_context("spawn").Pool(
            min(args.processes, len(taxonomy_tasks))
        ) as pool:
            taxonomy = tuple(pool.map(taxonomy_worker, taxonomy_tasks,
                                      chunksize=1))
    require(tuple(item[0] for item in taxonomy) == selected_indices,
            "taxonomy order")

    screen_totals = add_vectors(tuple(counts.values()), 4)
    taxonomy_totals = add_vectors(tuple(item[2] for item in taxonomy), 5)
    certificate_totals = add_vectors(tuple(item[3] for item in taxonomy), 4)
    exact_checked = sum(sum(item[2]) for item in taxonomy)
    require(screen_totals[2] == exact_checked, (screen_totals, exact_checked))
    require(taxonomy_totals[-1] == 0, taxonomy_totals)
    require(sum(item[6] for item in taxonomy) == len(selected_indices),
            "negative-alpha hostile count")
    overlap_records = tuple(
        record for item in taxonomy for record in item[5]
    )
    require(len(overlap_records) == taxonomy_totals[3] == 1,
            overlap_records)

    control_indices, expected_control_residuals = parse_control_indices(audit)
    require(set(control_indices).isdisjoint(selected_indices), "control overlap")
    control_tasks = tuple(
        (index, (LEVEL, *rows[index])) for index in control_indices
    )
    if args.processes == 1:
        controls = tuple(screen_worker(task) for task in control_tasks)
    else:
        with mp.get_context("spawn").Pool(
            min(args.processes, len(control_tasks))
        ) as pool:
            controls = tuple(pool.map(screen_worker, control_tasks, chunksize=1))
    control_records = []
    for index, row, direct, legacy in controls:
        require(row[12] == expected_control_residuals[index],
                (index, row[12], expected_control_residuals[index]))
        require(row[12] > 0 and row[16] == row[11],
                (index, "hostile control"))
        require(direct + legacy == row[11], (index, "control Farkas count"))
        quotient = Counter(lcm(*mask) // row[4] for mask in row[13])
        control_records.append(
            (
                index,
                tuple(row[position] for position in (9, 10, 11, 12)),
                direct,
                legacy,
                tuple(sorted(quotient.items())),
                hashlib.sha256(repr(row[13]).encode()).hexdigest(),
            )
        )
    control_records = tuple(control_records)

    remaining = tuple(
        index for index in current_live if index not in set(selected_indices)
    )
    remaining_ranked = audit.ranked_families(rows, components, remaining)
    ledger_after = LEDGER_BEFORE - len(selected_indices)
    wall_after = WALL_BEFORE - len(selected_indices)
    require(
        (len(remaining), len(remaining_ranked), ledger_after, wall_after)
        == (353, 31, 373_157, 353),
        "ledger consequence",
    )

    by_row = tuple(
        (
            index,
            counts[index],
            farkas_counts[index],
            next(item[2] for item in taxonomy if item[0] == index),
        )
        for index in selected_indices
    )
    canonical_packet = tuple(
        (index, canonical[index]) for index in selected_indices
    )
    taxonomy_packet = tuple(
        record for item in taxonomy for record in item[4]
    )
    semantic_packet = (
        AUDIT_SOURCE_SHA256,
        AUDIT_OUTPUT_SHA256,
        AUDIT_SEMANTIC_SHA256,
        hashlib.sha256(repr(rows).encode()).hexdigest(),
        tuple(map(len, (gcd8, gcd18, order_indices, natural_indices))),
        first_indices,
        second_indices,
        family_head(current_ranked, 4),
        selected_packet,
        canonical_packet,
        by_row,
        screen_totals,
        taxonomy_totals,
        certificate_totals,
        taxonomy_packet,
        overlap_records,
        control_records,
        (LEDGER_BEFORE, ledger_after, WALL_BEFORE, wall_after),
        family_head(remaining_ranked, 2),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    selected_cost = sum(family[0] for family in selected_families)
    selected_screen_sha = hashlib.sha256(
        repr(canonical_packet).encode()
    ).hexdigest()
    taxonomy_sha = hashlib.sha256(repr(taxonomy_packet).encode()).hexdigest()
    lines = [
        "LRC14 projected-k3 z216 gcd next-prefix exact probe",
        (
            f"atlas=rows:{len(rows)};wall:{sum(row[3] for row in rows)};"
            f"order:{sum(not row[3] for row in rows)};row_sha256:"
            f"{hashlib.sha256(repr(rows).encode()).hexdigest()};"
            "dual_parser_match:PASS"
        ),
        (
            f"current_queue=live:{len(current_live)};families:"
            f"{len(current_ranked)};head:{family_head(current_ranked, 4)}"
        ),
        (
            f"selection=families:{len(selected_families)};rows:"
            f"{len(selected_indices)};indices:{selected_indices};"
            f"keys:{family_keys};cost:{selected_cost};"
            "derived_from_complete_cost_rank:PASS;strict_boundary:PASS"
        ),
        (
            f"screen=states:{screen_totals[0]};crude:{screen_totals[1]};"
            f"status:{screen_totals[2]};residual:{screen_totals[3]};"
            f"exact_Farkas_rechecked:{exact_checked};"
            f"canonical_sha256:{selected_screen_sha}"
        ),
        (
            "taxonomy=coordinate_union,zero_reduced_union,two_fan,"
            f"two_threshold_overlap,weighted_core:{taxonomy_totals};"
            f"certificate_counts:{certificate_totals};"
            f"taxonomy_sha256:{taxonomy_sha}"
        ),
        f"by_row=index,screen,direct_inherited,taxonomy:{by_row}",
        f"two_threshold_overlap_record:{overlap_records}",
        (
            "overlap_lemma=if_low_tail_statuses_are_nonempty_and_every_"
            "high_tail_status_outside_coordinate_union_S_has_at_least_1+e_"
            "bits_then_total_incidence_ge_Tlow+e*max(0,Thigh-sum_S_mj)"
        ),
        (
            f"hostile_controls=index,screen,direct,inherited,quotient,"
            f"mask_sha256:{control_records};negative_alpha_mutations_rejected:"
            f"{sum(item[6] for item in taxonomy)}"
        ),
        (
            f"consequence=ledger:{LEDGER_BEFORE}-{len(selected_indices)}="
            f"{ledger_after};z216_wall:{WALL_BEFORE}-{len(selected_indices)}="
            f"{wall_after};remaining_families:{len(remaining_ranked)};"
            f"next_head:{family_head(remaining_ranked, 2)}"
        ),
        (
            "direction=empty_exact_necessary_upper_screen_excludes_each_"
            "selected_projected_wall_row;no_terminal_needed;no_converse_or_"
            "physical_lift"
        ),
        (
            "scope=projected_k3_z216_necessary_atlas_only;no_arbitrary_k_le_1_"
            "rung_endpoint_owner_phase_current_physical_cover_or_LRC14_claim"
        ),
        f"truth_gates=assert_nodes:{assert_nodes};float_literals:{float_nodes}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    print("\n".join(lines))


if __name__ == "__main__":
    mp.freeze_support()
    main()
