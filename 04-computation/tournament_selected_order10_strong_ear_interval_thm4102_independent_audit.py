#!/usr/bin/env python3
"""Independent exact audit for the provisional THM-4102.

This program does not import the THM-4102 primary compiler or the THM-4097
strong-ear compiler.  It loads only the hash-pinned order-eight isomorphism-
class generator because the selected bank is defined by that generator's
representative order.  All tournament arithmetic, extension, coding,
strongness, Hamiltonian DP, literal enumeration, and ear responses are rebuilt
here.

The ear evaluator is a path-cover cut-flow table.  It constructs start/end
subset path counts, contracts ordered two-path covers into a weighted cut,
and evaluates every cut signature by a subset recurrence.  This differs in
code structure from the primary Start/End/Q evaluator.  Every executable gate
uses ``require`` so ``python -O`` cannot remove a check.
"""

from __future__ import annotations

import importlib.util
import json
from hashlib import sha256
from itertools import permutations
from math import isqrt
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ENGINE = HERE / "strong_H_spectrum_m8_isoclass_monad_s5.py"
HISTOGRAM = ROOT / "05-knowledge/results/h_spectrum_n9_histogram_monad_s6.tsv"
LEGACY_VALUES = ROOT / "05-knowledge/results/strong_H_spectrum_m9_values_kps_S134.out"

ENGINE_SHA256 = "6ab922de4a8b6f6c15ee0ca7e0b036c3821b3e800dbdf961de72194e73346419"
HISTOGRAM_SHA256 = "e7d5594879d4c3af739cb94ca8cfd944879c4d586747d993dd6687e60126552f"
LEGACY_VALUES_SHA256 = "27fbef5b06fcf0369eeb602e513c3802ea171492e1292a3f6afa3efeadef9f55"
A000568 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}

# The full primary output freezes ten deterministic first witnesses.  The
# theorem table displays the eight boundary/lane rows between 249 and 15551.
EXPECTED_KEY_ROWS = {
    125: (75, 253, 7, 18_179_778_674_239),
    249: (75, 104, 3, 26_993_059_954_495),
    2887: (127, 158, 5, 25_133_469_073_343),
    2933: (119, 27, 4, 34_960_494_755_519),
    14649: (2575, 271, 5, 13_193_187_335_727),
    14653: (2517, 267, 4, 17_454_203_805_215),
    14655: (2393, 263, 4, 17_591_902_805_567),
    15055: (2741, 170, 4, 25_004_756_344_591),
    15551: (3081, 77, 4, 32_976_074_124_951),
    15621: (2903, 329, 4, 10_852_575_419_951),
}
LITERAL_CHILD_KEYS = (249, 15551)


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError("FAILED: " + label)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_order8_generator():
    require(file_sha256(ENGINE) == ENGINE_SHA256, "order-eight generator hash")
    spec = importlib.util.spec_from_file_location("thm4102_audit_engine_input", ENGINE)
    require(spec is not None and spec.loader is not None, "generator import specification")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def extend(adj: tuple[int, ...] | list[int], signature: int) -> tuple[int, ...]:
    """Add x=n, where signature bit i is one exactly when x -> i."""
    n = len(adj)
    require(0 <= signature < (1 << n), "extension signature range")
    child = list(adj) + [0]
    for old in range(n):
        if (signature >> old) & 1:
            child[n] |= 1 << old
        else:
            child[old] |= 1 << n
    return tuple(child)


def raw_code(adj: tuple[int, ...] | list[int]) -> int:
    code = 0
    bit_index = 0
    n = len(adj)
    for source in range(n):
        for target in range(source + 1, n):
            code |= ((adj[source] >> target) & 1) << bit_index
            bit_index += 1
    return code


def decode_raw(code: int, n: int) -> tuple[int, ...]:
    require(code >= 0, "nonnegative tournament code")
    adj = [0] * n
    bit_index = 0
    for source in range(n):
        for target in range(source + 1, n):
            if (code >> bit_index) & 1:
                adj[source] |= 1 << target
            else:
                adj[target] |= 1 << source
            bit_index += 1
    require(code >> bit_index == 0, "code fits stated tournament order")
    return tuple(adj)


def tournament_well_formed(adj: tuple[int, ...] | list[int]) -> bool:
    n = len(adj)
    full = (1 << n) - 1
    for vertex in range(n):
        if adj[vertex] & ~full or (adj[vertex] >> vertex) & 1:
            return False
    for left in range(n):
        for right in range(left + 1, n):
            if bool((adj[left] >> right) & 1) == bool((adj[right] >> left) & 1):
                return False
    return True


def hamilton_dp(adj: tuple[int, ...] | list[int]) -> int:
    """Independent endpoint Held--Karp count."""
    require(tournament_well_formed(adj), "Hamilton DP tournament input")
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for vertex in range(n):
        dp[1 << vertex][vertex] = 1
    for subset in range(1, 1 << n):
        row = dp[subset]
        for last, count in enumerate(row):
            if count == 0:
                continue
            available = adj[last] & ~subset
            while available:
                bit = available & -available
                nxt = bit.bit_length() - 1
                dp[subset | bit][nxt] += count
                available -= bit
    return sum(dp[-1])


def hamilton_literal(adj: tuple[int, ...] | list[int]) -> int:
    require(tournament_well_formed(adj), "literal Hamilton tournament input")
    n = len(adj)
    total = 0
    for word in permutations(range(n)):
        legal = True
        for index in range(n - 1):
            if not ((adj[word[index]] >> word[index + 1]) & 1):
                legal = False
                break
        total += legal
    return total


def reach(adj: tuple[int, ...] | list[int], reverse: bool) -> int:
    n = len(adj)
    seen = 1
    frontier = 1
    while frontier:
        nxt = 0
        vertices = frontier
        while vertices:
            bit = vertices & -vertices
            source = bit.bit_length() - 1
            vertices -= bit
            if reverse:
                for target in range(n):
                    if (adj[target] >> source) & 1:
                        nxt |= 1 << target
            else:
                nxt |= adj[source]
        nxt &= ~seen
        seen |= nxt
        frontier = nxt
    return seen


def is_strong(adj: tuple[int, ...] | list[int]) -> bool:
    require(tournament_well_formed(adj), "strongness tournament input")
    full = (1 << len(adj)) - 1
    return reach(adj, False) == full and reach(adj, True) == full


def endpoint_path_tables(adj: tuple[int, ...] | list[int]):
    """Return subset path counts indexed by last and independently by first."""
    n = len(adj)
    size = 1 << n
    incoming = [0] * n
    for source in range(n):
        targets = adj[source]
        while targets:
            bit = targets & -targets
            target = bit.bit_length() - 1
            incoming[target] |= 1 << source
            targets -= bit

    ends = [[0] * n for _ in range(size)]
    starts = [[0] * n for _ in range(size)]
    for vertex in range(n):
        ends[1 << vertex][vertex] = 1
        starts[1 << vertex][vertex] = 1
    for subset in range(1, size):
        if subset & (subset - 1) == 0:
            continue
        vertices = subset
        while vertices:
            bit = vertices & -vertices
            vertex = bit.bit_length() - 1
            vertices -= bit
            rest = subset ^ bit

            predecessors = incoming[vertex] & rest
            total = 0
            while predecessors:
                pred_bit = predecessors & -predecessors
                predecessor = pred_bit.bit_length() - 1
                total += ends[rest][predecessor]
                predecessors -= pred_bit
            ends[subset][vertex] = total

            successors = adj[vertex] & rest
            total = 0
            while successors:
                next_bit = successors & -successors
                successor = next_bit.bit_length() - 1
                total += starts[rest][successor]
                successors -= next_bit
            starts[subset][vertex] = total
    return ends, starts


def ear_truth_table(adj: tuple[int, ...] | list[int]) -> tuple[int, ...]:
    """Evaluate all new-vertex cuts via a weighted path-cover cut flow."""
    require(tournament_well_formed(adj), "ear evaluator tournament input")
    n = len(adj)
    size = 1 << n
    full = size - 1
    ends, starts = endpoint_path_tables(adj)

    pair_flow = [[0] * n for _ in range(n)]
    for left_set in range(1, full):
        right_set = full ^ left_set
        left_vertices = left_set
        while left_vertices:
            left_bit = left_vertices & -left_vertices
            left_end = left_bit.bit_length() - 1
            left_vertices -= left_bit
            left_count = ends[left_set][left_end]
            if left_count == 0:
                continue
            right_vertices = right_set
            while right_vertices:
                right_bit = right_vertices & -right_vertices
                right_start = right_bit.bit_length() - 1
                right_vertices -= right_bit
                right_count = starts[right_set][right_start]
                if right_count:
                    pair_flow[left_end][right_start] += left_count * right_count

    full_starts = starts[full]
    full_ends = ends[full]
    start_sum = [0] * size
    end_sum = [0] * size
    cut_flow = [0] * size
    responses = [0] * size
    total_ends = sum(full_ends)
    for signature in range(1, size):
        bit = signature & -signature
        added = bit.bit_length() - 1
        previous = signature ^ bit
        start_sum[signature] = start_sum[previous] + full_starts[added]
        end_sum[signature] = end_sum[previous] + full_ends[added]

        delta = 0
        for vertex in range(n):
            if (signature >> vertex) & 1:
                if vertex != added:
                    delta -= pair_flow[added][vertex]
            else:
                delta += pair_flow[vertex][added]
        cut_flow[signature] = cut_flow[previous] + delta
    for signature in range(size):
        responses[signature] = (
            start_sum[signature]
            + total_ends - end_sum[signature]
            + cut_flow[signature]
        )
    return tuple(responses)


def odd_intervals(values: set[int]) -> list[tuple[int, int, int]]:
    ordered = sorted(values)
    require(bool(ordered), "nonempty selected image")
    first = last = ordered[0]
    intervals = []
    for value in ordered[1:]:
        if value == last + 2:
            last = value
        else:
            intervals.append((first, last, (last - first) // 2 + 1))
            first = last = value
    intervals.append((first, last, (last - first) // 2 + 1))
    return intervals


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    return all(value % divisor for divisor in range(2, isqrt(value) + 1))


def first_missing_prime(values: set[int], start: int) -> int:
    candidate = start if start % 2 else start + 1
    while True:
        if is_prime(candidate) and candidate not in values:
            return candidate
        candidate += 2


def first_missing_seven_prime(values: set[int], start: int) -> int:
    candidate = start if start % 2 else start + 1
    while True:
        if is_prime(candidate) and 7 * candidate not in values:
            return candidate
        candidate += 2


def historical_order9_values():
    require(file_sha256(HISTOGRAM) == HISTOGRAM_SHA256, "historical histogram hash")
    require(file_sha256(LEGACY_VALUES) == LEGACY_VALUES_SHA256,
            "historical value-dump hash")
    histogram_values = set()
    all_classes = 0
    strong_classes = 0
    lines = HISTOGRAM.read_text(encoding="ascii").splitlines()
    require(lines and lines[0] == "H\tcount_all\tcount_strong", "histogram header")
    for line in lines[1:]:
        value_text, all_text, strong_text = line.split("\t")
        value = int(value_text)
        all_classes += int(all_text)
        strong_count = int(strong_text)
        strong_classes += strong_count
        if strong_count:
            histogram_values.add(value)
    legacy_values = {
        int(value)
        for value in LEGACY_VALUES.read_text(encoding="ascii").split()
    }
    require(histogram_values == legacy_values, "historical order-nine sets agree")
    require((all_classes, strong_classes) == (191_536, 178_133),
            "historical order-nine class totals")
    return histogram_values, all_classes, strong_classes


def small_evaluator_controls(representatives) -> int:
    """Compare every cut through old order five with DP and literal counting."""
    controls = 0
    for n in range(1, 6):
        for parent in representatives[n]:
            adj = tuple(parent)
            responses = ear_truth_table(adj)
            for signature, predicted in enumerate(responses):
                child = extend(adj, signature)
                dp_value = hamilton_dp(child)
                literal_value = hamilton_literal(child)
                require(predicted == dp_value, "small response versus independent DP")
                require(predicted == literal_value, "small response versus literal paths")
                controls += 1
    return controls


def main() -> None:
    generator = load_order8_generator()
    representatives, counts = generator.generate(8)
    require(counts == A000568, "A000568 representative counts through order eight")
    require(all(tournament_well_formed(tuple(adj))
                for adj in representatives[8]), "order-eight representative validity")

    small_controls = small_evaluator_controls(representatives)
    require(small_controls == 470, "frozen small evaluator control universe")

    historical_values, all_classes, strong_classes = historical_order9_values()
    selected: dict[int, tuple[int, ...]] = {}
    selected_sources: dict[int, tuple[int, int]] = {}
    strong_order8 = 0
    parent_cuts = 0
    for parent_index, parent_list in enumerate(representatives[8]):
        parent = tuple(parent_list)
        if not is_strong(parent):
            continue
        strong_order8 += 1
        responses = ear_truth_table(parent)
        for signature in range(1, (1 << 8) - 1):
            parent_cuts += 1
            value = responses[signature]
            if value not in selected:
                child = extend(parent, signature)
                require(hamilton_dp(child) == value,
                        "new selected order-nine parent DP value")
                require(is_strong(child), "new selected order-nine parent strongness")
                selected[value] = child
                selected_sources[value] = (parent_index, signature)

    require(strong_order8 == 6008, "strong order-eight parent count")
    require(parent_cuts == 1_526_032, "complete order-nine parent-ear universe")
    require(len(selected) == 1_482, "deterministic selected order-nine bank size")
    require(set(selected) == historical_values,
            "selected order-nine bank covers exact historical strong values")

    parent_failures = 0
    selected_parent_codes = []
    for value, parent in sorted(selected.items()):
        parent_failures += hamilton_dp(parent) != value or not is_strong(parent)
        selected_parent_codes.append((value, raw_code(parent)))
    require(parent_failures == 0, "all selected order-nine parents DP/strong")
    selected_bank_digest = sha256(
        json.dumps(selected_parent_codes, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    key_parent_values = sorted({row[0] for row in EXPECTED_KEY_ROWS.values()})
    literal_parent_rows = []
    for value in key_parent_values:
        literal_value = hamilton_literal(selected[value])
        require(literal_value == value, "key selected parent literal count")
        literal_parent_rows.append((value, raw_code(selected[value]), literal_value))

    values: set[int] = set()
    witnesses: dict[int, tuple[int, int, int]] = {}
    ears_checked = 0
    ear_strong_checks = 0
    for parent_h, parent in sorted(selected.items()):
        responses = ear_truth_table(parent)
        for signature in range(1, (1 << 9) - 1):
            ears_checked += 1
            child = extend(parent, signature)
            require(is_strong(child), "every selected nonconstant ear is strong")
            ear_strong_checks += 1
            value = responses[signature]
            values.add(value)
            if value not in witnesses:
                witnesses[value] = (parent_h, signature, raw_code(child))

    require(ears_checked == 1_482 * 510 == 755_820,
            "complete selected order-ten ear universe")
    require(ear_strong_checks == ears_checked, "all selected ears strongness checked")
    require(len(values) == 7_566, "selected order-ten image size")
    require((min(values), max(values)) == (125, 15_621),
            "selected order-ten image extrema")
    require(all(value % 2 == 1 for value in values), "selected values are odd")

    intervals = odd_intervals(values)
    require((249, 14_649, 7_201) in intervals, "primary solid interval")
    require((15_055, 15_551, 249) in intervals, "secondary solid interval")
    adjacent_holes = (247, 14_651, 15_053, 15_553)
    require(all(value not in values for value in adjacent_holes),
            "both solid intervals have the stated adjacent holes")
    expected_first_intervals = [
        (125, 125, 1),
        (135, 135, 1),
        (145, 147, 2),
        (153, 155, 2),
        (159, 161, 2),
        (165, 171, 4),
        (175, 231, 29),
        (235, 245, 6),
        (249, 14_649, 7_201),
        (14_653, 14_655, 2),
        (14_659, 14_659, 1),
        (14_663, 14_671, 5),
    ]
    require(intervals[:12] == expected_first_intervals,
            "frozen first selected-image components")

    witness_failures = 0
    retained_rows = []
    for value, (parent_h, signature, code) in sorted(witnesses.items()):
        child = decode_raw(code, 10)
        dp_value = hamilton_dp(child)
        strong = is_strong(child)
        witness_failures += dp_value != value or not strong
        retained_rows.append((value, parent_h, signature, code, dp_value, strong))
    require(witness_failures == 0, "all retained witnesses direct DP/strong")

    key_rows = {}
    literal_child_rows = []
    for value, expected in sorted(EXPECTED_KEY_ROWS.items()):
        parent_h, signature, code = witnesses[value]
        actual = (parent_h, signature, signature.bit_count(), code)
        require(actual == expected, "deterministic key witness row H=" + str(value))
        child = decode_raw(code, 10)
        dp_value = hamilton_dp(child)
        require(dp_value == value and is_strong(child),
                "key witness independent DP/strong H=" + str(value))
        key_rows[str(value)] = {
            "parent_h": parent_h,
            "signature": signature,
            "cut_weight": signature.bit_count(),
            "code": code,
            "dp_H": dp_value,
        }
        if value in LITERAL_CHILD_KEYS:
            literal_value = hamilton_literal(child)
            require(literal_value == value,
                    "key child literal count H=" + str(value))
            literal_child_rows.append((value, code, literal_value))

    require(14_651 == 49 * 13 * 23, "global prefix multiplicative bridge")
    prior_prefix = {
        value for value in range(1, 2_882, 2) if value not in {7, 21}
    }
    combined_prefix = prior_prefix | values | {14_651}
    prefix_missing = [
        value
        for value in range(1, 14_656, 2)
        if value not in {7, 21} and value not in combined_prefix
    ]
    require(not prefix_missing, "global allowed prefix through 14655")

    ordinary_missing = first_missing_prime(values, 2_887)
    require(ordinary_missing == 14_657 and 14_657 not in values,
            "ordinary prime lane cutoff")
    require(is_prime(14_653) and 14_653 in values,
            "ordinary prime lane last supplied value")

    seven_missing_prime = first_missing_seven_prime(values, 419)
    require(seven_missing_prime == 2_111, "seven-prime lane cutoff prime")
    require(7 * seven_missing_prime == 14_777 and 14_777 not in values,
            "seven-prime lane cutoff value")
    require(is_prime(2_099) and 7 * 2_099 == 14_693 and 14_693 in values,
            "seven-prime lane last supplied value")

    longest = sorted(intervals, key=lambda row: row[2], reverse=True)[:12]
    ledger = {
        "input_hashes": [ENGINE_SHA256, HISTOGRAM_SHA256, LEGACY_VALUES_SHA256],
        "A000568": A000568,
        "small_evaluator_controls": small_controls,
        "historical_class_totals": [all_classes, strong_classes],
        "strong_order8": strong_order8,
        "parent_cuts": parent_cuts,
        "selected_parent_codes": selected_parent_codes,
        "selected_sources": sorted(selected_sources.items()),
        "selected_bank_digest": selected_bank_digest,
        "literal_parent_rows": literal_parent_rows,
        "ears_checked": ears_checked,
        "ear_strong_checks": ear_strong_checks,
        "values": sorted(values),
        "intervals": intervals,
        "adjacent_holes": adjacent_holes,
        "retained_rows": retained_rows,
        "key_rows": key_rows,
        "literal_child_rows": literal_child_rows,
        "global_prefix": [14_651, len(prefix_missing)],
        "lane_cutoffs": [ordinary_missing, seven_missing_prime, 7 * seven_missing_prime],
    }
    semantic = sha256(
        json.dumps(ledger, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("THM-4102 INDEPENDENT SELECTED-EAR AUDIT")
    print("A000568_counts=", [counts[n] for n in range(1, 9)])
    print("small_response_DP_literal_controls=", small_controls, "failures= 0")
    print("strong_order8_parents=", strong_order8,
          "nonconstant_parent_ears=", parent_cuts)
    print("selected_order9_parents=", len(selected),
          "historical_sets_equal=", set(selected) == historical_values)
    print("selected_parent_bank_sha256=", selected_bank_digest)
    print("direct_parent_checks=", len(selected), "failures=", parent_failures)
    print("literal_key_parent_checks=", len(literal_parent_rows), "rows=", literal_parent_rows)
    print("nonconstant_order10_ears=", ears_checked,
          "direct_strongness_checks=", ear_strong_checks)
    print("selected_value_count=", len(values), "min=", min(values), "max=", max(values))
    print("primary_solid_interval= [249,14649] count=7201 adjacent_holes=(247,14651)")
    print("secondary_solid_interval= [15055,15551] count=249 adjacent_holes=(15053,15553)")
    print("first_intervals=", intervals[:12])
    print("longest_intervals=", longest)
    print("direct_retained_witness_checks=", len(witnesses),
          "failures=", witness_failures)
    print("key_witnesses=", key_rows)
    print("literal_key_child_checks=", literal_child_rows)
    print("global_allowed_prefix_through=14655 bridge_14651=49*13*23 missing=",
          len(prefix_missing))
    print("ordinary_prime_lane_through=14653 next_unforced=", ordinary_missing)
    print("seven_prime_lane_p_through=2099 next_unforced_value=", 7 * seven_missing_prime)
    print("semantic_sha256=", semantic)
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
