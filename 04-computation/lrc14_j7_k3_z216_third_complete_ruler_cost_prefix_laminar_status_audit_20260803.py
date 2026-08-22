#!/usr/bin/env python3
"""Exact third z216 complete-ruler prefix and laminar status audit.

This scratch companion reconstructs the 357 projected-k3 wall rows left by
the two pinned intrinsic-cost prefixes.  It screens the next complete prefix
(``gcd72/L72072`` and ``gcd24/L25872``), retains a known empty-screen control
and a known nonempty residual control, and independently rebuilds every status
instance used here.

It also revisits the eleven former weighted-core states.  After the earlier
coordinate-union, zero-reduced-union and two-fan templates, it searches a
canonical solver-free language of nested capacity thresholds covered
pointwise by a nonnegative modular function of the four status bits.  This is
a projected necessary-atlas computation, not a physical cover or LRC(14).
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import itertools
import multiprocessing as mp
from collections import Counter, defaultdict
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SECOND_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_second_nontrivial_ruler_cost_prefix_closure_scout_20260803.py"
)
SECOND_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_second_nontrivial_ruler_cost_prefix_closure_scout_20260803.out"
)
AUDIT_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.py"
)
AUDIT_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.out"
)
NATURAL_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.py"
)
ORDER_SOURCE = ROOT / (
    "04-computation/lrc14_j7_k3_z216_order_row_screen_descent_thm3270.py"
)

DEPENDENCIES = (
    (SECOND_SOURCE,
     "bd2236e4e3fe6c47984e4ecf10cd0b87f2d4b937a8d5c0ebb6da70e81884dd5e"),
    (SECOND_OUTPUT,
     "a5d68cb56dfeb3f44e49133451c3b2445ca387370a0a21998c0a6be5531ed8c9"),
    (AUDIT_SOURCE,
     "d3581bc937597270558a77c0a02398a4cfd368306f5a302be01240cdb4cef7c3"),
    (AUDIT_OUTPUT,
     "e035531d0af15998a8a07d083042ff2486c80c7fdd729d902922179409d1725a"),
    (NATURAL_SOURCE,
     "430dee7ba03e0d5c9ae0df72ac512500de4f7056cb4663d1c8468bfb93a49bfe"),
)

LEVEL = 216
EXPECTED_ROW_SHA256 = "53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649"
EXPECTED_CENSUS = (480, 447, 33, 123, 357, 33)
EXPECTED_RANK_HEAD = (
    (
        1_729_728, 1, 72, 72_072, (157,),
        ((157, (1, 4, 9, 11, 12, 13), 24, 1_729_728),),
    ),
    (
        2_121_504, 3, 24, 25_872, (53, 293, 357),
        (
            (53, (1, 2, 6, 8, 11, 14), 26, 672_672),
            (293, (2, 4, 6, 8, 11, 14), 26, 672_672),
            (357, (2, 6, 8, 11, 12, 14), 30, 776_160),
        ),
    ),
    (
        2_293_200, 1, 24, 76_440, (141,),
        ((141, (1, 4, 6, 10, 13, 14), 30, 2_293_200),),
    ),
)
EXPECTED_SELECTED_INDICES = (157, 53, 293, 357)
EXPECTED_SELECTED_COST = 3_851_232
EXPECTED_SELECTED_PACKET_SHA256 = (
    "2f57b06bff9e2149f8b33cf3af3cd2498a029edf9ba26d012ae03f63c2a8ac89"
)
EXPECTED_SCREEN_COUNTS = (
    (157, (111, 93, 18, 0)),
    (53, (27, 17, 10, 0)),
    (293, (14, 9, 5, 0)),
    (357, (271, 170, 101, 0)),
)
EXPECTED_TOTALS = (423, 289, 134, 0)
EXPECTED_FAMILY_TOTALS = (
    ((72, 72_072), (111, 93, 18, 0)),
    ((24, 25_872), (312, 196, 116, 0)),
)
EXPECTED_NEW_TAXONOMY = (128, 2, 3, 1, 0)
OLD_WEIGHTED_SCREEN_COUNTS = (
    (134, (29, 28, 1, 0)),
    (121, (95, 79, 16, 0)),
    (17, (25, 12, 13, 0)),
    (27, (36, 18, 18, 0)),
    (56, (64, 40, 24, 0)),
    (138, (234, 136, 98, 0)),
    (298, (207, 118, 89, 0)),
)
EXPECTED_OLD_LAMINAR_SUPPORT_COUNTS = ((2, 10), (3, 1))
EXPECTED_NEW_LAMINAR_SUPPORT_COUNTS = ((2, 1),)
EXPECTED_PROOF_SHAPES = (
    (((1, 2), (0, 2, 2, 1, 1)), 5),
    (((1, 4, 5), (0, 1, 3, 1, 1)), 1),
    (((2, 5), (0, 2, 1, 1, 0)), 1),
    (((3, 4), (0, 2, 1, 1, 2)), 1),
    (((3, 4), (0, 2, 1, 2, 1)), 2),
    (((3, 5), (0, 2, 1, 1, 0)), 2),
)
EXPECTED_OLD_LAMINAR_SHA256 = (
    "03b8a519871574c74fc3223babb502c0e5bbc9fe7055bfc746a6f6b1921786ab"
)
EXPECTED_NEW_LAMINAR_SHA256 = (
    "d7b68a6cd091efafc8898906f4c0b3ff7d9b05ba338f4b163dba8a841796bf83"
)
EXPECTED_SCREEN_SHA256 = "17bd632e09975d5a84d4e5f8667023eb9c1b9de21c8ae16f639612ff3d380e74"
EXPECTED_EMPTY_RESIDUAL_SHA256 = (
    "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d"
)
CONTROL_INDICES = (14, 64)
EXPECTED_CONTROL_TOTALS = (
    (14, (59, 24, 35, 0), 0),
    (64, (115, 50, 57, 8), 8),
)
EXPECTED_RESIDUAL_CONTROL_MASK_SHA256 = (
    "098c78ce3def914ca1c0ad26458558898fe624040222f429e6163b9a8e9855c4"
)
LEDGER_BEFORE = 373_161
LEDGER_AFTER = 373_157
WALL_BEFORE = 357
WALL_AFTER = 353
FAMILIES_BEFORE = 33
FAMILIES_AFTER = 31
EXPECTED_SEMANTIC_SHA256 = (
    "c97d088f75478c11c0b82d36c6f0dffe0bf99ae50d2bad9a803630609440424e"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, (path, "bare CR"))
    return hashlib.sha256(payload).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def totals(rows):
    return tuple(sum(row[position] for row in rows) for position in range(4))


def screen_worker(item):
    index, task = item
    second = load(f"third_prefix_second_screen_{index}", SECOND_SOURCE)
    return second.screen_worker((index, task))


def tail_demand(histogram, threshold):
    return sum(count for load_value, count in histogram if load_value >= threshold)


def canonical_laminar_modular_cover(q, capacities, marginals, histogram):
    """Shortest equal-weight nested-threshold cover, checked on all 16 states.

    For thresholds t in a selected chain, seek nonnegative integer c_0,...,c_4
    such that pointwise

      sum_t 1[capacity(P)>=t] <= c_0 + sum_i c_(i+1) 1[i in P].

    Summing against a hypothetical common status table gives the claimed
    modular upper bound.  Search order proves minimal threshold support inside
    this finite language; no optimizer or Farkas normalization is used.
    """
    thresholds = tuple(
        threshold for threshold, _count in histogram
        if threshold > 0
        and not all(capacity >= threshold for capacity in capacities)
    )
    candidates = []
    for support_size in range(1, min(3, len(thresholds)) + 1):
        for chosen in itertools.combinations(thresholds, support_size):
            demand = sum(tail_demand(histogram, threshold) for threshold in chosen)
            for coefficients in itertools.product(
                range(support_size + 1), repeat=5
            ):
                if not all(
                    sum(capacities[pattern] >= threshold for threshold in chosen)
                    <= coefficients[0]
                    + sum(
                        coefficients[coordinate + 1]
                        * ((pattern >> coordinate) & 1)
                        for coordinate in range(4)
                    )
                    for pattern in range(16)
                ):
                    continue
                supply = coefficients[0] * q + sum(
                    coefficients[coordinate + 1] * marginals[coordinate]
                    for coordinate in range(4)
                )
                if demand > supply:
                    candidates.append(
                        (
                            support_size, supply, sum(coefficients), chosen,
                            coefficients, demand, demand - supply,
                        )
                    )
        if candidates:
            break
    if not candidates:
        return None
    _size, supply, _weight, chosen, coefficients, demand, gap = min(candidates)
    return chosen, coefficients, demand, supply, gap


def detail_worker(item):
    index, task, expected_counts = item
    second = load(f"third_prefix_second_detail_{index}", SECOND_SOURCE)
    natural = load(f"third_prefix_natural_detail_{index}", NATURAL_SOURCE)
    eng = second.status_engine(natural, f"third_{index}")
    eng.FIRST = LEVEL
    eng.ray.FIRST = LEVEL
    _first, body, ruler, high, wall = task
    stream = eng.ray.Stream(body)
    require((stream.L, stream.high_floor) == (ruler, high), (index, "atlas"))
    _trials, states, _checks, _signs = eng.ray.ray_quotient_states(stream)
    crude, status, residual = eng.exact_common_status_screen(stream, states)
    counts = (len(states), len(crude), len(status), len(residual))
    require(counts == expected_counts, (index, counts, expected_counts))
    require(wall and not residual, (index, "unexpected residual"))

    branches = Counter()
    records = []
    laminar_records = []
    for divisors, witness in sorted(status.items()):
        q, cofactor, marginals, capacity_set, histogram, certificate = witness
        divisor_lcm = lcm(*divisors)
        require(divisor_lcm == q * cofactor, (index, divisors, "cofactor"))
        rebuilt, capacities = eng.ray.local.hunter_status_data(
            divisor_lcm, divisors, q
        )
        require(rebuilt == marginals, (index, divisors, "marginals"))
        require(tuple(sorted(set(capacities))) == capacity_set,
                (index, divisors, "capacities"))
        eng.independent_farkas_check(
            q, marginals, capacities, histogram, certificate
        )

        unions = second.coordinate_union_certificates(
            capacities, marginals, histogram
        )
        reduced = () if unions else second.zero_reduced_union_certificates(
            capacities, marginals, histogram
        )
        fans = () if (unions or reduced) else second.two_fan_certificates(
            capacities, marginals, histogram
        )
        laminar = None if (unions or reduced or fans) else (
            canonical_laminar_modular_cover(
                q, capacities, marginals, histogram
            )
        )
        if unions:
            branch, evidence = "coordinate_union", unions
        elif reduced:
            branch, evidence = "zero_reduced_union", reduced
        elif fans:
            branch, evidence = "two_fan", fans
        elif laminar is not None:
            branch, evidence = "laminar_modular", laminar
            laminar_records.append(
                (
                    index, divisors, q, marginals, capacity_set, histogram,
                    tuple(capacities), laminar,
                )
            )
        else:
            branch, evidence = "weighted_core", ()
        branches[branch] += 1
        records.append(
            (
                index, divisors, q, marginals, capacity_set, histogram,
                branch, evidence,
            )
        )
    return (
        index, counts,
        (
            branches["coordinate_union"],
            branches["zero_reduced_union"],
            branches["two_fan"],
            branches["laminar_modular"],
            branches["weighted_core"],
        ),
        tuple(records), tuple(laminar_records),
    )


def ranked_families(rows, components, live):
    families = defaultdict(list)
    for index in live:
        families[(gcd(LEVEL, rows[index][1]), rows[index][1])].append(index)
    ranked = []
    for (divisor_gcd, ruler), raw_indices in families.items():
        indices = tuple(raw_indices)
        packet = tuple(
            (index, rows[index][0], components[index], ruler * components[index])
            for index in indices
        )
        ranked.append(
            (
                sum(item[3] for item in packet), len(indices), divisor_gcd,
                ruler, indices, packet,
            )
        )
    return tuple(sorted(ranked))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    for path, digest in DEPENDENCIES:
        require(lf_sha(path) == digest, (path, "dependency changed"))
    require(
        "semantic_sha256=a73b57558c5fd74736fa9551c72f94d20b2c053ef7b8f36c27222ae993f8dcd8"
        in SECOND_OUTPUT.read_text(encoding="utf-8"),
        "second-prefix semantic pin",
    )
    require(
        "semantic_sha256=4bb4eca9949eb4f07b618f9e49bbb5de63492364a29a4a5abfc84a6acf105e11"
        in AUDIT_OUTPUT.read_text(encoding="utf-8")
        and "all_exact_controls=PASS" in AUDIT_OUTPUT.read_text(encoding="utf-8"),
        "independent-audit semantic pin",
    )
    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require((assert_nodes, float_nodes) == (0, 0), "truth-gate syntax")

    second = load("third_prefix_second_main", SECOND_SOURCE)
    first = load("third_prefix_first_main", second.FIRST_SOURCE)
    natural = load("third_prefix_natural_main", NATURAL_SOURCE)
    order = load("third_prefix_order_main", ORDER_SOURCE)
    rows, components = natural.atlas_rows()
    require(
        hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_SHA256,
        "row order",
    )
    prior_closed = (
        set(natural.EXPECTED_UNION_INDICES)
        | set(order.EXPECTED_ORDER_INDICES)
        | set(first.EXPECTED_SELECTED_INDICES)
        | set(second.EXPECTED_SELECTED_INDICES)
        | {
            index for index, row in enumerate(rows)
            if gcd(LEVEL, row[1]) in (8, 18)
        }
    )
    live = tuple(
        index for index, row in enumerate(rows)
        if row[3] and index not in prior_closed
    )
    ranked = ranked_families(rows, components, live)
    census = (
        len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows),
        len(prior_closed), len(live), len(ranked),
    )
    require(census == EXPECTED_CENSUS, census)
    require(ranked[:3] == EXPECTED_RANK_HEAD, ranked[:3])

    selected_families = ranked[:2]
    selected_indices = tuple(
        index for family in selected_families for index in family[4]
    )
    require(selected_indices == EXPECTED_SELECTED_INDICES, selected_indices)
    require(selected_families[-1][1] > 1, "prefix did not end at a family")
    require(selected_families[-1][0] < ranked[2][0], "cost tie crosses boundary")
    selected_packet = tuple(
        (
            index, rows[index][0], gcd(LEVEL, rows[index][1]), rows[index][1],
            components[index], rows[index][1] * components[index],
            rows[index][3], rank + 1,
        )
        for rank, family in enumerate(selected_families)
        for index in family[4]
    )
    require(sum(row[5] for row in selected_packet) == EXPECTED_SELECTED_COST,
            "selected cost")
    selected_packet_sha = hashlib.sha256(
        repr(selected_packet).encode()
    ).hexdigest()
    if EXPECTED_SELECTED_PACKET_SHA256 is not None:
        require(selected_packet_sha == EXPECTED_SELECTED_PACKET_SHA256,
                selected_packet_sha)

    screen_indices = selected_indices + CONTROL_INDICES
    screen_tasks = tuple((index, (LEVEL, *rows[index])) for index in screen_indices)
    if args.processes == 1:
        screened = tuple(screen_worker(task) for task in screen_tasks)
    else:
        with mp.get_context("spawn").Pool(
            min(args.processes, len(screen_tasks))
        ) as pool:
            screened = tuple(pool.map(screen_worker, screen_tasks, chunksize=1))
    require(tuple(row[0] for row in screened) == screen_indices, "screen order")

    canonical = {}
    count_by_index = {}
    farkas_by_index = {}
    for index, row, direct, legacy in screened:
        source = rows[index]
        divisor_gcd = gcd(LEVEL, source[1])
        require(
            row[:6] == (
                LEVEL, source[0], source[1], source[2],
                source[1] // divisor_gcd, source[3],
            ),
            (index, row[:6]),
        )
        require(row[16] == row[11], (index, "status verification"))
        require(direct + legacy == row[11], (index, "Farkas count"))
        canonical[index] = tuple(row[:19])
        count_by_index[index] = tuple(row[p] for p in (9, 10, 11, 12))
        farkas_by_index[index] = (direct, legacy)
    require(
        tuple((index, count_by_index[index]) for index in selected_indices)
        == EXPECTED_SCREEN_COUNTS,
        tuple((index, count_by_index[index]) for index in selected_indices),
    )
    selected_totals = totals(tuple(count_by_index[i] for i in selected_indices))
    require(selected_totals == EXPECTED_TOTALS, selected_totals)
    canonical_selected = tuple(
        (index, canonical[index]) for index in selected_indices
    )
    screen_sha = hashlib.sha256(repr(canonical_selected).encode()).hexdigest()
    require(screen_sha == EXPECTED_SCREEN_SHA256, screen_sha)
    residual = tuple(
        (index, canonical[index][13])
        for index in selected_indices if canonical[index][12]
    )
    residual_sha = hashlib.sha256(repr(residual).encode()).hexdigest()
    require(not residual and residual_sha == EXPECTED_EMPTY_RESIDUAL_SHA256,
            residual_sha)
    family_totals = tuple(
        (
            (family[2], family[3]),
            totals(tuple(count_by_index[index] for index in family[4])),
        )
        for family in selected_families
    )
    require(family_totals == EXPECTED_FAMILY_TOTALS, family_totals)

    old_counts = dict(OLD_WEIGHTED_SCREEN_COUNTS)
    detail_indices = tuple(old_counts) + selected_indices
    detail_tasks = tuple(
        (
            index, (LEVEL, *rows[index]),
            count_by_index[index] if index in count_by_index else old_counts[index],
        )
        for index in detail_indices
    )
    if args.processes == 1:
        details = tuple(detail_worker(task) for task in detail_tasks)
    else:
        with mp.get_context("spawn").Pool(
            min(args.processes, len(detail_tasks))
        ) as pool:
            details = tuple(pool.map(detail_worker, detail_tasks, chunksize=1))
    require(tuple(row[0] for row in details) == detail_indices, "detail order")
    detail_by_index = {row[0]: row for row in details}

    new_taxonomy = tuple(
        sum(detail_by_index[index][2][position] for index in selected_indices)
        for position in range(5)
    )
    require(new_taxonomy == EXPECTED_NEW_TAXONOMY, new_taxonomy)
    old_laminar = tuple(
        record for index in old_counts for record in detail_by_index[index][4]
    )
    new_laminar = tuple(
        record for index in selected_indices for record in detail_by_index[index][4]
    )
    require(len(old_laminar) == 11, len(old_laminar))
    require(len(new_laminar) == 1, len(new_laminar))
    require(
        not any(detail_by_index[index][2][4] for index in detail_indices),
        "unclassified weighted core",
    )
    old_support_counts = tuple(sorted(Counter(
        len(record[-1][0]) for record in old_laminar
    ).items()))
    new_support_counts = tuple(sorted(Counter(
        len(record[-1][0]) for record in new_laminar
    ).items()))
    require(old_support_counts == EXPECTED_OLD_LAMINAR_SUPPORT_COUNTS,
            old_support_counts)
    require(new_support_counts == EXPECTED_NEW_LAMINAR_SUPPORT_COUNTS,
            new_support_counts)
    proof_shapes = tuple(sorted(Counter(
        (record[-1][0], record[-1][1])
        for record in old_laminar + new_laminar
    ).items()))
    require(proof_shapes == EXPECTED_PROOF_SHAPES, proof_shapes)
    require(
        hashlib.sha256(repr(old_laminar).encode()).hexdigest()
        == EXPECTED_OLD_LAMINAR_SHA256,
        "old laminar packet",
    )
    require(
        hashlib.sha256(repr(new_laminar).encode()).hexdigest()
        == EXPECTED_NEW_LAMINAR_SHA256,
        "new laminar packet",
    )

    for index, expected_counts, expected_residual in EXPECTED_CONTROL_TOTALS:
        require(count_by_index[index] == expected_counts,
                (index, count_by_index[index]))
        require(canonical[index][12] == expected_residual,
                (index, canonical[index][12]))
    residual_control = canonical[64][13]
    require(
        hashlib.sha256(repr(residual_control).encode()).hexdigest()
        == EXPECTED_RESIDUAL_CONTROL_MASK_SHA256,
        "residual control hash",
    )
    residual_quotient = tuple(sorted(Counter(
        lcm(*mask) // canonical[64][4] for mask in residual_control
    ).items()))
    require(residual_quotient == ((8, 8),), residual_quotient)

    remaining = tuple(index for index in live if index not in set(selected_indices))
    remaining_ranked = ranked_families(rows, components, remaining)
    require((len(remaining), len(remaining_ranked)) == (353, 31),
            (len(remaining), len(remaining_ranked)))
    require(LEDGER_BEFORE - len(selected_indices) == LEDGER_AFTER, "ledger")
    require(WALL_BEFORE - len(selected_indices) == WALL_AFTER, "wall")
    require(FAMILIES_BEFORE - len(selected_families) == FAMILIES_AFTER,
            "families")

    taxonomy_packet = tuple(
        record for index in detail_indices for record in detail_by_index[index][3]
    )
    semantic_packet = (
        tuple(digest for _path, digest in DEPENDENCIES),
        EXPECTED_ROW_SHA256,
        census,
        tuple((family[0], family[1], family[2], family[3], family[4])
              for family in ranked[:3]),
        selected_packet,
        canonical_selected,
        selected_totals,
        family_totals,
        taxonomy_packet,
        old_laminar,
        new_laminar,
        proof_shapes,
        tuple((index, canonical[index]) for index in CONTROL_INDICES),
        residual_quotient,
        LEDGER_BEFORE, LEDGER_AFTER, WALL_BEFORE, WALL_AFTER,
        FAMILIES_BEFORE, FAMILIES_AFTER,
        tuple((family[0], family[1], family[2], family[3], family[4])
              for family in remaining_ranked[:3]),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    lines = [
        "LRC14 projected-k3 z216 third complete-ruler prefix and laminar status audit",
        (
            "dependencies="
            f"second:{DEPENDENCIES[0][1]}/{DEPENDENCIES[1][1]};"
            f"independent_audit:{DEPENDENCIES[2][1]}/{DEPENDENCIES[3][1]};"
            f"THM3281:{DEPENDENCIES[4][1]}"
        ),
        (
            f"universe=z1:{LEVEL};atlas_rows:{census[0]};wall:{census[1]};"
            f"order:{census[2]};prior_closed:{census[3]};live_wall:{census[4]};"
            f"live_families:{census[5]};row_sha256:{EXPECTED_ROW_SHA256}"
        ),
        (
            f"rank_head={tuple((f[0], f[1], f[2], f[3], f[4]) for f in ranked[:3])}"
        ),
        (
            "selection=cost_ordered_complete_family_prefix_through_next_"
            f"non_singleton;families:2;rows:{len(selected_indices)};"
            f"indices:{selected_indices};cost:{EXPECTED_SELECTED_COST};"
            f"packet_sha256:{selected_packet_sha};strict_boundary:PASS"
        ),
    ]
    for family, summary in zip(selected_families, family_totals):
        lines.append(
            f"family=gcd{family[2]}/L{family[3]};rows:{family[1]};"
            f"indices:{family[4]};cost:{family[0]};states:{summary[1][0]};"
            f"crude:{summary[1][1]};status:{summary[1][2]};"
            f"residual:{summary[1][3]}"
        )
    lines.extend(
        [
            (
                f"selected_screen=states:{selected_totals[0]};"
                f"crude:{selected_totals[1]};status:{selected_totals[2]};"
                f"residual:{selected_totals[3]};"
                f"legacy_Farkas_checked:{sum(farkas_by_index[i][1] for i in selected_indices)};"
                f"screen_sha256:{screen_sha};residual_sha256:{residual_sha}"
            ),
            (
                f"new_status_taxonomy=coordinate_union:{new_taxonomy[0]};"
                f"zero_reduced_union:{new_taxonomy[1]};two_fan:{new_taxonomy[2]};"
                f"laminar_modular:{new_taxonomy[3]};unclassified:{new_taxonomy[4]};"
                f"exact_Farkas_checked:{selected_totals[2]}"
            ),
            (
                f"old_weighted_resolution=states:{len(old_laminar)};"
                f"threshold_support_counts:{old_support_counts};"
                f"new_laminar_support_counts:{new_support_counts};"
                f"proof_shapes:{proof_shapes}"
            ),
            f"old_laminar_certificates={old_laminar}",
            f"new_laminar_certificates={new_laminar}",
            (
                "laminar_mechanism=capacity_tail_events_are_nested;the_sum_"
                "of_selected_tail_indicators_is_pointwise_bounded_on_all_16_"
                "patterns_by_c0_plus_a_nonnegative_modular_weight_of_the_four_"
                "status_bits;tail_demand_exceeds_c0q_plus_the_same_weighted_"
                "marginals;search_over_one_two_then_three_equal_weight_"
                "thresholds_is_exhaustive_inside_this_declared_language"
            ),
            (
                "controls=index14:59=24crude+35status+0residual;"
                "index64:115=50crude+57status+8residual;"
                f"residual_quotient:{residual_quotient};"
                f"residual_mask_sha256:{EXPECTED_RESIDUAL_CONTROL_MASK_SHA256};"
                "pinned_independent_audit:PASS"
            ),
            (
                f"consequence=ledger:{LEDGER_BEFORE}-{len(selected_indices)}="
                f"{LEDGER_AFTER};z216_wall:{WALL_BEFORE}-{len(selected_indices)}="
                f"{WALL_AFTER};families:{FAMILIES_BEFORE}-"
                f"{len(selected_families)}={FAMILIES_AFTER};cap:216"
            ),
            (
                "direction=empty_exact_necessary_upper_screens_exclude_the_"
                "four_selected_projected_wall_rows;no_converse_from_screen_"
                "feasibility_and_cost_is_only_a_work_order"
            ),
            (
                f"stopping_boundary=remaining_wall:{len(remaining)};"
                f"families:{len(remaining_ranked)};next_head:"
                f"{tuple((f[0], f[1], f[2], f[3], f[4]) for f in remaining_ranked[:3])}"
            ),
            (
                "scope=projected_k3_z216_necessary_atlas_only;no_terminal_"
                "probe_divisor_action_endpoint_owner_phase_current_physical_"
                "cover_arbitrary_k_le_1_rung_or_LRC14_claim"
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
    mp.freeze_support()
    main()
