#!/usr/bin/env python3
"""Finite exact order-nine variance/image determinacy scout.

This script rebuilds the complete THM-4097 order-eight strong-ear stream and,
for each attained order-nine Hamilton count H, retains the first at most 32
distinct labelled child codes in deterministic generator/signature order.  It
then groups those parents by the exact scalar triple (H, F1, Var) and compares
their complete nonconstant ear-image sets.

Status: FINITE-EXACT EVIDENCE.  A zero-conflict result in this selected finite
universe is not an all-order theorem and is not an order-nine completeness
claim.  Every executable gate uses ``require`` and survives ``python -O``.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path


HERE = Path(__file__).resolve().parent
THM4097 = HERE / "tournament_order9_strong_ear_spectrum_thm4097.py"

THM4097_SHA256 = "610ca5850b272e0e75c574f2c1a710a0b96c75cc7191b1e1f1a03dfbdd1378d6"
THM4102_SELECTED_ORDER9_BANK_SHA256 = (
    "c03c203943e734d09bee4b8818227b8f184405ce4c5092dd56d0fdb6107d528c"
)
EXPECTED_SEMANTIC_SHA256 = "5f139bbcd2bfd615e4abd10d6aa8f52946eb4935bb2c0e6eb32aa694da9d4279"

CAP = 32
EXPECTED_CANDIDATE_COUNT = 46_314
EXPECTED_REPEATED_TRIPLE_GROUPS = 11_938
EXPECTED_EQUAL_TRIPLE_PAIRS = 88_527
EXPECTED_COUNT_HISTOGRAM = {
    1: 3,
    3: 2,
    5: 2,
    9: 14,
    10: 4,
    11: 1,
    12: 3,
    13: 2,
    14: 1,
    15: 2,
    16: 2,
    17: 2,
    18: 9,
    20: 3,
    21: 1,
    22: 2,
    24: 3,
    25: 1,
    26: 1,
    27: 5,
    28: 1,
    29: 1,
    32: 1_417,
}


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def canonical_digest(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(payload).hexdigest()


def load_module(name: str, path: Path, expected_hash: str):
    require(file_sha256(path) == expected_hash, f"frozen {name} source hash")
    specification = importlib.util.spec_from_file_location(name, path)
    require(
        specification is not None and specification.loader is not None,
        f"{name} import specification",
    )
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def update_row_digest(hasher, row: object) -> None:
    hasher.update(json.dumps(row, separators=(",", ":")).encode("ascii"))
    hasher.update(b"\n")


def cut_field(
    starts: list[int], ends: list[int], exposed: list[list[int]]
) -> tuple[list[list[int]], list[int]]:
    """Extract THM-4115's integral symmetric weights and zero-sum field."""
    order = len(starts)
    weights = [[0] * order for _ in range(order)]
    for left, right in combinations(range(order), 2):
        symmetric = exposed[left][right] + exposed[right][left]
        require(symmetric % 2 == 0, "integral symmetric cut weight")
        weights[left][right] = weights[right][left] = symmetric // 2
    field: list[int] = []
    for vertex in range(order):
        column = sum(exposed[other][vertex] for other in range(order))
        row = sum(exposed[vertex])
        require((column - row) % 2 == 0, "integral orientation field")
        field.append(starts[vertex] - ends[vertex] + (column - row) // 2)
    require(sum(field) == 0, "zero-sum orientation field")
    return weights, field


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")

    parent = load_module("thm4097_for_cap32", THM4097, THM4097_SHA256)
    engine = parent.load_engine()
    representatives, counts = engine.generate(8)
    require(counts == parent.A000568, "A000568 class counts through order eight")

    # Deterministic first-CAP distinct raw labelled codes in each H fibre.
    bank: dict[int, list[int]] = defaultdict(list)
    seen: dict[int, set[int]] = defaultdict(set)
    strong_order8 = 0
    ears_checked = 0
    for old in representatives[8]:
        if not engine.is_strong(8, old):
            continue
        strong_order8 += 1
        state = parent.boundary_state(old)
        for signature in range(1, (1 << 8) - 1):
            ears_checked += 1
            parent_h = parent.insertion_h(state, signature)
            if len(bank[parent_h]) >= CAP:
                continue
            child = engine.extend(old, 8, signature)
            code = parent.raw_code(child)
            if code not in seen[parent_h]:
                seen[parent_h].add(code)
                bank[parent_h].append(code)

    require(strong_order8 == 6_008, "complete strong order-eight parent count")
    require(ears_checked == 1_526_032, "complete order-eight nonconstant-ear stream")
    historical, all_classes, strong_classes = parent.read_historical_strong_values()
    require(set(bank) == historical, "all exact order-nine H fibres represented")
    require(
        (all_classes, strong_classes) == (191_536, 178_133),
        "historical order-nine class controls",
    )
    require(len(bank) == 1_482, "order-nine H-fibre count")

    count_histogram = Counter(len(codes) for codes in bank.values())
    require(dict(count_histogram) == EXPECTED_COUNT_HISTOGRAM, "cap-32 count histogram")
    candidate_count = sum(len(codes) for codes in bank.values())
    require(candidate_count == EXPECTED_CANDIDATE_COUNT, "cap-32 candidate count")

    first_selected_bank = sorted((parent_h, codes[0]) for parent_h, codes in bank.items())
    first_selected_digest = canonical_digest(first_selected_bank)
    require(
        first_selected_digest == THM4102_SELECTED_ORDER9_BANK_SHA256,
        "rank-one slice equals the frozen THM-4102 selected bank",
    )
    universe_rows = [
        [parent_h, rank, code]
        for parent_h in sorted(bank)
        for rank, code in enumerate(bank[parent_h], start=1)
    ]
    universe_digest = canonical_digest(universe_rows)

    edge_pairs = tuple(combinations(range(9), 2))
    profile_hasher = sha256()
    group_hasher = sha256()
    profile_count = 0
    triple_group_count = 0
    repeated_triple_groups = 0
    equal_triple_pairs = 0
    different_image_pairs = 0
    different_image_groups = 0
    maximum_triple_multiplicity = 0

    # H is part of the key, so one fibre can be released before the next.
    for parent_h in sorted(bank):
        groups: dict[tuple[int, int], dict[tuple[int, ...], int]] = defaultdict(
            lambda: defaultdict(int)
        )
        for rank, code in enumerate(bank[parent_h], start=1):
            adjacency = parent.decode_raw(code, 9)
            require(parent.raw_code(adjacency) == code, "raw-code round trip")
            require(engine.is_strong(9, adjacency), "retained parent strongness")
            state = parent.boundary_state(adjacency)
            direct_h = sum(state[1])
            require(direct_h == sum(state[0]) == parent_h, "retained parent H")

            weights, field = cut_field(*state)
            total_weight = sum(weights[left][right] for left, right in edge_pairs)
            defect_one = 2 * total_weight - 8 * parent_h
            require(defect_one >= 0, "nonnegative one-defect layer")
            four_variance = sum(entry * entry for entry in field) + sum(
                weights[left][right] ** 2 for left, right in edge_pairs
            )

            values = [parent.insertion_h(state, signature) for signature in range(1 << 9)]
            require(values[0] == values[-1] == parent_h, "constant-cut controls")
            require(all(value >= parent_h for value in values), "insertion lower support")
            require(all(value % 2 == 1 for value in values), "Redei odd parity")
            twice_mean = 2 * parent_h + total_weight
            require(
                2 * sum(values) == len(values) * twice_mean,
                "exact uniform mean identity",
            )
            require(
                sum((2 * value - twice_mean) ** 2 for value in values)
                == len(values) * four_variance,
                "exact Parseval variance identity",
            )

            image = tuple(sorted(set(values[1:-1])))
            require(bool(image), "nonempty nonconstant ear image")
            groups[(defect_one, four_variance)][image] += 1
            update_row_digest(
                profile_hasher,
                [parent_h, rank, code, defect_one, four_variance, list(image)],
            )
            profile_count += 1

        triple_group_count += len(groups)
        for (defect_one, four_variance), image_counts in sorted(groups.items()):
            multiplicity = sum(image_counts.values())
            maximum_triple_multiplicity = max(maximum_triple_multiplicity, multiplicity)
            if multiplicity > 1:
                repeated_triple_groups += 1
            equal_triple_pairs += multiplicity * (multiplicity - 1) // 2
            same_image_pairs = sum(
                count * (count - 1) // 2 for count in image_counts.values()
            )
            cross_image_pairs = multiplicity * (multiplicity - 1) // 2 - same_image_pairs
            different_image_pairs += cross_image_pairs
            different_image_groups += int(len(image_counts) > 1)
            update_row_digest(
                group_hasher,
                [
                    parent_h,
                    defect_one,
                    four_variance,
                    [[list(image), count] for image, count in sorted(image_counts.items())],
                ],
            )

    require(profile_count == candidate_count, "every retained parent profiled")
    require(
        repeated_triple_groups == EXPECTED_REPEATED_TRIPLE_GROUPS,
        "repeated exact-triple group count",
    )
    require(
        equal_triple_pairs == EXPECTED_EQUAL_TRIPLE_PAIRS,
        "equal exact-triple pair count",
    )
    require(
        different_image_groups == 0 and different_image_pairs == 0,
        "no image conflict inside the finite cap-32 universe",
    )

    ledger = {
        "status": "FINITE-EXACT EVIDENCE",
        "scope": "first at most 32 distinct labelled witnesses per exact order-nine H fibre",
        "not_claimed": "NO ALL-ORDER DETERMINACY CLAIM",
        "source_hashes": [THM4097_SHA256],
        "strong_order8": strong_order8,
        "ears_checked": ears_checked,
        "h_fibres": len(bank),
        "cap": CAP,
        "candidate_count": candidate_count,
        "count_histogram": sorted(count_histogram.items()),
        "first_selected_digest": first_selected_digest,
        "universe_digest": universe_digest,
        "profile_digest": profile_hasher.hexdigest(),
        "group_digest": group_hasher.hexdigest(),
        "triple_group_count": triple_group_count,
        "repeated_triple_groups": repeated_triple_groups,
        "equal_triple_pairs": equal_triple_pairs,
        "different_image_groups": different_image_groups,
        "different_image_pairs": different_image_pairs,
        "maximum_triple_multiplicity": maximum_triple_multiplicity,
        "exact_key": "(H,F1,4Var), where 4Var=sum_i h_i^2+sum_ij w_ij^2",
        "image": "sorted distinct H(T+x_S), 1<=S<=510",
    }
    semantic = canonical_digest(ledger)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(semantic == EXPECTED_SEMANTIC_SHA256, "frozen semantic digest")

    print("ORDER-9 EAR VARIANCE/IMAGE CAP-32 SCOUT")
    print("status=FINITE-EXACT EVIDENCE")
    print("scope=first_at_most_32_distinct_labelled_witnesses_per_H_fibre")
    print("source_stream=6008_strong_order8_parents_x_254_cuts=1526032")
    print("H_fibres=", len(bank), "candidate_parents=", candidate_count)
    print("candidate_count_histogram=", sorted(count_histogram.items()))
    print("rank_one_THM4102_bank_sha256=", first_selected_digest)
    print("cap32_universe_sha256=", universe_digest)
    print("profile_ledger_sha256=", profile_hasher.hexdigest())
    print("triple_group_ledger_sha256=", group_hasher.hexdigest())
    print(
        "triple_groups=", triple_group_count,
        "repeated_groups=", repeated_triple_groups,
        "equal_triple_pairs=", equal_triple_pairs,
        "max_multiplicity=", maximum_triple_multiplicity,
    )
    print(
        "same_H_F1_Var_different_image_groups=", different_image_groups,
        "pairs=", different_image_pairs,
    )
    print("semantic_sha256=", semantic)
    print("NO_ALL_ORDER_DETERMINACY_CLAIM")
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
