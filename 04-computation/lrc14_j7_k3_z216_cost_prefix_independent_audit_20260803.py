#!/usr/bin/env python3
"""Independent exact audit of the two z216 intrinsic-ruler cost prefixes.

This audit does not import either 20260803 prefix scout or any of their
expected constants.  It rebuilds the projected-k3 z1=216 atlas through the
older pinned THM-3281/THM-3270 routes, reconstructs the prior disjoint
closures, derives both complete family prefixes from the cost order, replays
all twenty-three exact screens, and checks the solver-free status taxonomy.

The maintained atlas is a necessary quotient only.  Empty exact upper
screens exclude the selected projected rows; no physical LRC consequence or
converse from screen feasibility is asserted.
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
NATURAL_SOURCE = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.py"
)
ORDER_SOURCE = ROOT / (
    "04-computation/lrc14_j7_k3_z216_order_row_screen_descent_thm3270.py"
)
GCD18_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_unique_gcd18_terminal_descent_thm3261.out"
)
GCD8_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_low_cost_gcd8_terminal_descent_thm3264.out"
)
ORDER_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_order_row_screen_descent_thm3270.out"
)
NATURAL_OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.out"
)
COSTLY_REPORT = ROOT / (
    "05-knowledge/results/lrc14-z216-costly-gcd8-closure-opus-20260803.md"
)
COSTLY_AUDIT = ROOT / (
    ".scratch/lrc14_z216_gcd8_independent_audit_gcd8_independent_audit.out"
)

DEPENDENCIES = (
    (NATURAL_SOURCE,
     "430dee7ba03e0d5c9ae0df72ac512500de4f7056cb4663d1c8468bfb93a49bfe"),
    (ORDER_SOURCE,
     "dbcc60e3d691483e023486cdb9b935381e5172cea2b828b51ed3da99560fe2ab"),
    (GCD18_OUTPUT,
     "d7dae9cd7e8f0305824b30a2b8683abe6d7d8c9bd0509a0d84026322fd65344c"),
    (GCD8_OUTPUT,
     "2eebe5f7acb2d1c02fa126ee03166cd274ed209cc52d4b8a729a8fbc5f0a9782"),
    (ORDER_OUTPUT,
     "0fa283cc68ac831579e17b4e9d817cd5dd13a4c58363d51335e0c1911897d1f0"),
    (NATURAL_OUTPUT,
     "b9c23b21bf9766efbbc14aa97e2bb4268ddb6abb09b54a2b7424b8744ff686b2"),
    (COSTLY_REPORT,
     "2dfdaebf4ccd4b051d6839450c9039db4b2112d5839470a65ccc067f7b9d8fd5"),
    (COSTLY_AUDIT,
     "1432ffaca0025137da91ed392313c0a9da48d6d657c7b6faefcc157660c5e92a"),
)

LEVEL = 216
NATURAL_FAMILY_KEYS = ((24, 2352), (36, 8820), (72, 2520))
LOW_GCD8_COST_LIMIT = 2_000_000
EXPECTED_SEMANTIC_SHA256 = (
    "4bb4eca9949eb4f07b618f9e49bbb5de63492364a29a4a5abfc84a6acf105e11"
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


def status_engine(natural, suffix: str):
    """Reach the old exact status engine without a new-prefix dependency."""
    entry = load(f"cost_prefix_audit_entry_{suffix}", natural.SOURCE_3139)
    source = entry.load(f"cost_prefix_audit_source_{suffix}")
    driver = source.load(f"cost_prefix_audit_driver_{suffix}")
    predecessor = driver.load(f"cost_prefix_audit_predecessor_{suffix}")
    bridge = predecessor.load(f"cost_prefix_audit_bridge_{suffix}")
    base = bridge.load(f"cost_prefix_audit_base_{suffix}")
    return base.thm.screen_engine.eng


def screen_worker(item):
    index, task = item
    natural = load(f"cost_prefix_audit_screen_{index}", NATURAL_SOURCE)
    return natural.screen_worker((index, task))


def ranked_families(rows, components, live):
    families = defaultdict(list)
    for index in live:
        families[(gcd(LEVEL, rows[index][1]), rows[index][1])].append(index)
    ranked = []
    for (divisor_gcd, ruler), raw_indices in families.items():
        indices = tuple(raw_indices)
        packet = tuple(
            (
                index,
                rows[index][0],
                components[index],
                ruler * components[index],
            )
            for index in indices
        )
        ranked.append(
            (
                sum(item[3] for item in packet),
                len(indices),
                divisor_gcd,
                ruler,
                indices,
                packet,
            )
        )
    ranked = tuple(sorted(ranked))
    for family in ranked:
        key = (family[2], family[3])
        complete = tuple(
            index for index in live
            if (gcd(LEVEL, rows[index][1]), rows[index][1]) == key
        )
        require(family[4] == complete, (key, "incomplete family"))
    return ranked


def through_next_nonsingleton(ranked):
    stopping_rank = next(
        rank for rank, family in enumerate(ranked) if family[1] > 1
    )
    selected = ranked[:stopping_rank + 1]
    require(
        stopping_rank + 1 == len(ranked)
        or selected[-1][0] < ranked[stopping_rank + 1][0],
        "cost tie crosses prefix boundary",
    )
    return selected


def tail_demand(histogram, threshold):
    return sum(
        count for load_value, count in histogram if load_value >= threshold
    )


def union_cuts(capacities, marginals, histogram):
    cuts = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = tail_demand(histogram, threshold)
        for coordinate_mask in range(1, 16):
            event = tuple(
                bool(pattern & coordinate_mask) for pattern in range(16)
            )
            good = tuple(
                capacities[pattern] >= threshold for pattern in range(16)
            )
            if event != good:
                continue
            coordinates = tuple(
                coordinate for coordinate in range(4)
                if coordinate_mask & (1 << coordinate)
            )
            supply = sum(marginals[coordinate] for coordinate in coordinates)
            if demand > supply:
                cuts.append(
                    (threshold, coordinate_mask, coordinates, demand, supply,
                     demand - supply)
                )
    return tuple(cuts)


def zero_reduced_union_cuts(capacities, marginals, histogram):
    forbidden = sum(
        1 << coordinate for coordinate, marginal in enumerate(marginals)
        if marginal == 0
    )
    if not forbidden:
        return ()
    allowed_patterns = tuple(
        pattern for pattern in range(16) if not (pattern & forbidden)
    )
    cuts = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = tail_demand(histogram, threshold)
        for coordinate_mask in range(1, 16):
            if coordinate_mask & forbidden:
                continue
            event = tuple(
                bool(pattern & coordinate_mask) for pattern in allowed_patterns
            )
            good = tuple(
                capacities[pattern] >= threshold
                for pattern in allowed_patterns
            )
            if event != good:
                continue
            coordinates = tuple(
                coordinate for coordinate in range(4)
                if coordinate_mask & (1 << coordinate)
            )
            supply = sum(marginals[coordinate] for coordinate in coordinates)
            if demand > supply:
                cuts.append(
                    (threshold, forbidden, coordinate_mask, coordinates,
                     demand, supply, demand - supply)
                )
    return tuple(cuts)


def two_fan_cuts(capacities, marginals, histogram):
    cuts = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = tail_demand(histogram, threshold)
        for a, b in itertools.permutations(range(4), 2):
            c, d = tuple(sorted(set(range(4)).difference((a, b))))
            event = tuple(
                bool(pattern & (1 << b))
                or (
                    bool(pattern & (1 << a))
                    and bool(pattern & ((1 << c) | (1 << d)))
                )
                for pattern in range(16)
            )
            good = tuple(
                capacities[pattern] >= threshold for pattern in range(16)
            )
            if event != good:
                continue
            supply = marginals[b] + min(
                marginals[a], marginals[c] + marginals[d]
            )
            if demand > supply:
                cuts.append(
                    (threshold, (a, b, c, d), demand, supply,
                     demand - supply)
                )
    return tuple(cuts)


def taxonomy_worker(item):
    index, task, expected_counts = item
    natural = load(f"cost_prefix_audit_taxonomy_{index}", NATURAL_SOURCE)
    eng = status_engine(natural, str(index))
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
    certificate_counts = Counter()
    records = []
    weighted = []
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

        unions = union_cuts(capacities, marginals, histogram)
        reduced = () if unions else zero_reduced_union_cuts(
            capacities, marginals, histogram
        )
        fans = () if (unions or reduced) else two_fan_cuts(
            capacities, marginals, histogram
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
        else:
            branch = "weighted_core"
            evidence = ()
            weighted.append(
                (index, divisors, q, cofactor, marginals, capacity_set,
                 histogram, tuple(capacities))
            )
        branches[branch] += 1
        certificate_counts[branch] += len(evidence)
        records.append(
            (index, divisors, q, cofactor, marginals, capacity_set,
             histogram, branch, evidence)
        )
    summary = (
        branches["coordinate_union"],
        branches["zero_reduced_union"],
        branches["two_fan"],
        branches["weighted_core"],
    )
    certificate_summary = (
        certificate_counts["coordinate_union"],
        certificate_counts["zero_reduced_union"],
        certificate_counts["two_fan"],
    )
    require(sum(summary) == counts[2], (index, summary, counts))
    return (
        index, counts, summary, certificate_summary, tuple(records),
        tuple(weighted),
    )


def add_counts(rows):
    return tuple(sum(row[position] for row in rows) for position in range(4))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=3)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    for path, digest in DEPENDENCIES:
        require(lf_sha(path) == digest, (path, "dependency changed"))

    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
    float_nodes = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(syntax)
    )
    require((assert_nodes, float_nodes) == (0, 0), "truth-gate syntax")

    natural = load("cost_prefix_audit_natural", NATURAL_SOURCE)
    order = load("cost_prefix_audit_order", ORDER_SOURCE)
    rows, components = natural.atlas_rows()
    order_rows, order_components = order.atlas_rows()
    require((rows, components) == (order_rows, order_components),
            "independent atlas parser disagreement")
    require(len(rows) == len(components), "atlas/component mismatch")

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
        if (gcd(LEVEL, row[1]), row[1]) in NATURAL_FAMILY_KEYS
    }
    closure_parts = (gcd8, gcd18, order_indices, natural_indices)
    for left in range(len(closure_parts)):
        for right in range(left + 1, len(closure_parts)):
            require(
                closure_parts[left].isdisjoint(closure_parts[right]),
                (left, right, "prior closures overlap"),
            )
    require(
        tuple(map(len, closure_parts)) == (19, 1, 33, 47),
        tuple(map(len, closure_parts)),
    )
    require(all(rows[index][3] for index in gcd8 | gcd18 | natural_indices),
            "wall closure contains order row")

    low_gcd8 = {
        index for index in gcd8
        if rows[index][1] * components[index] <= LOW_GCD8_COST_LIMIT
    }
    costly_gcd8 = gcd8.difference(low_gcd8)
    require((len(low_gcd8), len(costly_gcd8)) == (17, 2),
            (low_gcd8, costly_gcd8))
    require(
        tuple(sorted(costly_gcd8)) == tuple(
            index for index in sorted(gcd8)
            if rows[index][1] == 560_560
            and components[index] == 34
        ),
        "costly gcd8 address reconstruction",
    )

    old_texts = {path: path.read_text(encoding="utf-8") for path, _ in DEPENDENCIES}
    lineage_needles = (
        (GCD18_OUTPUT, "consequence=ledger:373284-1=373283"),
        (GCD8_OUTPUT, "consequence=ledger:373283-17=373266"),
        (ORDER_OUTPUT, "consequence=ledger:373266-33=373233"),
        (NATURAL_OUTPUT, "consequence=ledger:373233-47=373186"),
        (COSTLY_REPORT, "projected ledger: 373186 -> 373184"),
        (COSTLY_AUDIT, "consequence=373186-2=373184"),
        (COSTLY_AUDIT, "all_exact_controls=PASS"),
    )
    for path, needle in lineage_needles:
        require(needle in old_texts[path], (path, needle))

    prior_closed = set().union(*closure_parts)
    require(len(prior_closed) == 100, len(prior_closed))
    first_live = tuple(
        index for index, row in enumerate(rows)
        if row[3] and index not in prior_closed
    )
    first_ranked = ranked_families(rows, components, first_live)
    first_families = through_next_nonsingleton(first_ranked)
    first_indices = tuple(
        index for family in first_families for index in family[4]
    )

    second_closed = prior_closed | set(first_indices)
    second_live = tuple(
        index for index, row in enumerate(rows)
        if row[3] and index not in second_closed
    )
    second_ranked = ranked_families(rows, components, second_live)
    second_families = through_next_nonsingleton(second_ranked)
    second_indices = tuple(
        index for family in second_families for index in family[4]
    )
    require(set(first_indices).isdisjoint(second_indices), "prefix overlap")
    require((len(first_indices), len(second_indices)) == (5, 18),
            (first_indices, second_indices))
    require(
        (len(first_live), len(first_ranked), len(second_live), len(second_ranked))
        == (380, 39, 375, 36),
        "live census",
    )

    selected_indices = first_indices + second_indices
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
    require(tuple(row[0] for row in screened) == selected_indices,
            "screen order")

    canonical = {}
    counts = {}
    farkas = {}
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
        require(row[12] == 0 and row[13] == (), (index, "residual"))
        canonical[index] = tuple(row[:19])
        counts[index] = tuple(row[position] for position in (9, 10, 11, 12))
        farkas[index] = (direct, legacy)

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
    require(tuple(row[0] for row in taxonomy) == selected_indices,
            "taxonomy order")
    taxonomy_by_index = {row[0]: row for row in taxonomy}

    def prefix_summary(indices):
        screen_counts = add_counts(tuple(counts[index] for index in indices))
        branches = add_counts(
            tuple(taxonomy_by_index[index][2] for index in indices)
        )
        certificate_counts = tuple(
            sum(taxonomy_by_index[index][3][position] for index in indices)
            for position in range(3)
        )
        exact_checked = sum(
            sum(taxonomy_by_index[index][2]) for index in indices
        )
        require(screen_counts[2] == exact_checked,
                (screen_counts, exact_checked))
        return screen_counts, branches, certificate_counts, exact_checked

    first_summary = prefix_summary(first_indices)
    second_summary = prefix_summary(second_indices)
    require(first_summary[0][3] == second_summary[0][3] == 0,
            "nonempty residual")

    remaining = tuple(
        index for index in second_live if index not in set(second_indices)
    )
    remaining_ranked = ranked_families(rows, components, remaining)
    require((len(remaining), len(remaining_ranked)) == (357, 33),
            (len(remaining), len(remaining_ranked)))
    ledger = (373_184, 373_184 - len(first_indices),
              373_184 - len(selected_indices))
    walls = (len(first_live), len(second_live), len(remaining))

    first_packet = tuple(
        (index, rows[index], components[index]) for index in first_indices
    )
    second_packet = tuple(
        (index, rows[index], components[index]) for index in second_indices
    )
    canonical_packet = tuple(
        (index, canonical[index]) for index in selected_indices
    )
    taxonomy_packet = tuple(
        record for row in taxonomy for record in row[4]
    )
    weighted_packet = tuple(
        record for row in taxonomy for record in row[5]
    )
    by_row = tuple(
        (index, counts[index], taxonomy_by_index[index][2])
        for index in selected_indices
    )
    semantic_packet = (
        tuple(digest for _path, digest in DEPENDENCIES),
        hashlib.sha256(repr(rows).encode()).hexdigest(),
        tuple(map(len, closure_parts)),
        tuple(sorted(costly_gcd8)),
        first_ranked,
        first_packet,
        second_ranked,
        second_packet,
        canonical_packet,
        by_row,
        first_summary,
        second_summary,
        taxonomy_packet,
        weighted_packet,
        ledger,
        walls,
        tuple((family[0], family[1], family[2], family[3], family[4])
              for family in remaining_ranked[:2]),
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    def family_head(ranked, count):
        return tuple(
            (family[0], family[1], family[2], family[3], family[4])
            for family in ranked[:count]
        )

    first_cost = sum(family[0] for family in first_families)
    second_cost = sum(family[0] for family in second_families)
    lines = [
        "LRC14 projected-k3 z216 intrinsic-ruler cost-prefix independent audit",
        (
            f"atlas=rows:{len(rows)};wall:{sum(row[3] for row in rows)};"
            f"order:{sum(not row[3] for row in rows)};row_sha256:"
            f"{hashlib.sha256(repr(rows).encode()).hexdigest()};"
            "dual_parser_match:PASS"
        ),
        (
            "prior_partition="
            f"gcd8:{len(gcd8)}=low:{len(low_gcd8)}+costly:{len(costly_gcd8)};"
            f"gcd18:{len(gcd18)};order:{len(order_indices)};"
            f"natural:{len(natural_indices)};union:{len(prior_closed)};"
            "pairwise_disjoint:PASS"
        ),
        (
            f"baseline=lineage:373284-1-17-33-47-2=373184;"
            f"z216_wall:{sum(row[3] for row in rows)}-19-1-"
            f"{len(natural_indices)}={len(first_live)};"
            f"costly_indices:{tuple(sorted(costly_gcd8))};"
            "pinned_prior_outputs:PASS"
        ),
        (
            f"first_queue=live:{len(first_live)};families:{len(first_ranked)};"
            f"head:{family_head(first_ranked, 4)}"
        ),
        (
            f"first_prefix=families:{len(first_families)};rows:"
            f"{len(first_indices)};indices:{first_indices};cost:{first_cost};"
            "complete_fibres:PASS;strict_boundary:PASS"
        ),
        (
            f"first_screen=states:{first_summary[0][0]};"
            f"crude:{first_summary[0][1]};status:{first_summary[0][2]};"
            f"residual:{first_summary[0][3]};taxonomy:{first_summary[1]};"
            f"certificate_counts:{first_summary[2]};"
            f"exact_Farkas_checked:{first_summary[3]}"
        ),
        (
            f"second_queue=live:{len(second_live)};"
            f"families:{len(second_ranked)};head:{family_head(second_ranked, 5)}"
        ),
        (
            f"second_prefix=families:{len(second_families)};rows:"
            f"{len(second_indices)};indices:{second_indices};"
            f"cost:{second_cost};complete_fibres:PASS;strict_boundary:PASS;"
            "disjoint_from_first:PASS"
        ),
        (
            f"second_screen=states:{second_summary[0][0]};"
            f"crude:{second_summary[0][1]};status:{second_summary[0][2]};"
            f"residual:{second_summary[0][3]};taxonomy:{second_summary[1]};"
            f"certificate_counts:{second_summary[2]};"
            f"exact_Farkas_checked:{second_summary[3]}"
        ),
        f"by_row=index,screen_counts,taxonomy:{by_row}",
        (
            "taxonomy_soundness=all_status_instances_rebuilt_from_divisors;"
            "all_returned_duals_exactly_checked;coordinate_union_bound:"
            "tail_le_sum_marginals;zero_reduced_bound:zero_marginal_deletes_"
            "incident_patterns;two_fan_bound:mB+min(mA,mC+mD);"
            f"weighted_core_exact_only:{len(weighted_packet)}"
        ),
        (
            f"consequence=ledger:{ledger[0]}-{len(first_indices)}="
            f"{ledger[1]}-{len(second_indices)}={ledger[2]};"
            f"z216_wall:{walls[0]}-{len(first_indices)}={walls[1]}-"
            f"{len(second_indices)}={walls[2]};remaining_families:"
            f"{len(remaining_ranked)}"
        ),
        (
            f"next_head:{family_head(remaining_ranked, 2)};"
            "cost_is_work_order_only_not_safety_invariant"
        ),
        (
            "direction=empty_exact_necessary_upper_screen_excludes_each_"
            "selected_projected_wall_row;no_converse_or_physical_lift"
        ),
        (
            "scope=projected_k3_z216_necessary_atlas_only;no_arbitrary_k_le_1_"
            "rung_endpoint_owner_phase_current_physical_cover_or_LRC14_claim"
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
