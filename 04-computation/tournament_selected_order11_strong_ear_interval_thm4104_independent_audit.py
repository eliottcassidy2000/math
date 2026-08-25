#!/usr/bin/env python3
"""Independent exact audit for the provisional THM-4104.

This program imports neither the THM-4104 primary compiler nor the THM-4102
primary compiler.  It loads only the hash-pinned order-eight isomorphism-class
generator, whose representative order is part of the deterministic selection
rule.  Tournament extension, raw coding, strongness, Hamiltonian DP, literal
enumeration, and both ear evaluators are rebuilt here.

The universal evaluator expands the ear count as an exact quadratic Boolean
cut polynomial and evaluates its entire truth table by a low-bit recurrence.
A separately coded memoized ordered-two-path-cover formula and direct child
Held--Karp counts audit all small controls, every retained key/boundary row,
and a broad deterministic sample.  Nonconstant ears over the checked strong
parents are strong by THM-4097's proved strong-ear lemma.  Every executable
gate uses ``require``; there are no floating-point computations.
"""

from __future__ import annotations

from functools import lru_cache
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

EXPECTED_N9_BANK_SHA256 = "c03c203943e734d09bee4b8818227b8f184405ce4c5092dd56d0fdb6107d528c"
EXPECTED_N10_BANK_SHA256 = "2f3fbd5d7f56de24a1f08ea08585dd029c70344ef444830915b5ea0d203e4b92"

EXPECTED_KEY_ROWS: dict[int, tuple[int, int, int, int]] = {
    225: (125, 32, 1, 34_095_368_048_213_567),
    429: (125, 116, 4, 33_813_343_248_580_159),
    14_657: (773, 898, 4, 4_980_304_472_833_599),
    14_777: (697, 283, 5, 25_105_481_806_249_279),
    80_265: (13_253, 903, 6, 2_181_427_721_976_895),
    80_387: (14_443, 47, 5, 31_503_056_742_464_543),
    80_405: (13_493, 283, 5, 22_513_037_423_950_911),
    80_633: (14_527, 665, 5, 15_198_813_324_868_767),
    80_875: (15_059, 182, 5, 33_617_973_800_943_135),
    84_259: (12_667, 15, 4, 36_024_248_539_118_623),
    93_751: (15_581, 841, 5, 3_645_901_060_717_615),
}


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError("FAILED: " + label)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def load_order8_generator():
    require(file_sha256(ENGINE) == ENGINE_SHA256, "order-eight generator hash")
    spec = importlib.util.spec_from_file_location("thm4104_audit_engine_input", ENGINE)
    require(spec is not None and spec.loader is not None,
            "order-eight generator import specification")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def extend(adj: tuple[int, ...] | list[int], signature: int) -> tuple[int, ...]:
    """Add x=n, with signature bit i one exactly when x -> i."""
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
    require(code >= 0, "nonnegative raw code")
    adj = [0] * n
    bit_index = 0
    for source in range(n):
        for target in range(source + 1, n):
            if (code >> bit_index) & 1:
                adj[source] |= 1 << target
            else:
                adj[target] |= 1 << source
            bit_index += 1
    require(code >> bit_index == 0, "raw code fits stated order")
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
    """Subset Hamilton-path counts indexed independently by end and start."""
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


def quadratic_ear_table(adj: tuple[int, ...] | list[int]) -> tuple[int, ...]:
    """All ear values from the exact quadratic Boolean cut polynomial."""
    require(tournament_well_formed(adj), "quadratic ear tournament input")
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
    constant = sum(full_ends)
    linear = [
        full_starts[vertex]
        - full_ends[vertex]
        + sum(pair_flow[source][vertex] for source in range(n))
        for vertex in range(n)
    ]

    responses = [0] * size
    responses[0] = constant
    for signature in range(1, size):
        bit = signature & -signature
        added = bit.bit_length() - 1
        previous = signature ^ bit
        correction = 0
        vertices = previous
        while vertices:
            old_bit = vertices & -vertices
            old = old_bit.bit_length() - 1
            vertices -= old_bit
            correction += pair_flow[old][added] + pair_flow[added][old]
        responses[signature] = responses[previous] + linear[added] - correction
    require(all(value >= 0 for value in responses), "nonnegative quadratic ear values")
    return tuple(responses)


def separate_ear_value(adj: tuple[int, ...] | list[int], signature: int) -> int:
    """Independent memoized ordered-two-path-cover formula for one cut."""
    require(tournament_well_formed(adj), "separate ear tournament input")
    n = len(adj)
    full = (1 << n) - 1
    require(0 <= signature <= full, "separate ear signature range")

    @lru_cache(maxsize=None)
    def end_count(subset: int, last: int) -> int:
        bit = 1 << last
        if not subset & bit:
            return 0
        if subset == bit:
            return 1
        previous = subset ^ bit
        return sum(
            end_count(previous, predecessor)
            for predecessor in range(n)
            if (previous >> predecessor) & 1
            and (adj[predecessor] >> last) & 1
        )

    @lru_cache(maxsize=None)
    def start_count(subset: int, first: int) -> int:
        bit = 1 << first
        if not subset & bit:
            return 0
        if subset == bit:
            return 1
        rest = subset ^ bit
        return sum(
            start_count(rest, successor)
            for successor in range(n)
            if (rest >> successor) & 1 and (adj[first] >> successor) & 1
        )

    total = sum(
        start_count(full, first)
        for first in range(n)
        if (signature >> first) & 1
    )
    total += sum(
        end_count(full, last)
        for last in range(n)
        if not (signature >> last) & 1
    )
    for left_set in range(1, full):
        right_set = full ^ left_set
        left_vertices = left_set & ~signature
        while left_vertices:
            left_bit = left_vertices & -left_vertices
            left_end = left_bit.bit_length() - 1
            left_vertices -= left_bit
            left_paths = end_count(left_set, left_end)
            if left_paths == 0:
                continue
            right_vertices = right_set & signature
            while right_vertices:
                right_bit = right_vertices & -right_vertices
                right_start = right_bit.bit_length() - 1
                right_vertices -= right_bit
                total += left_paths * start_count(right_set, right_start)
    return total


def odd_intervals(values: set[int]) -> list[tuple[int, int, int]]:
    ordered = sorted(values)
    require(bool(ordered), "nonempty selected order-eleven image")
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


def previous_prime(value: int) -> int:
    candidate = value - 1 - (value % 2)
    while not is_prime(candidate):
        candidate -= 2
    return candidate


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
        all_classes += int(all_text)
        strong_count = int(strong_text)
        strong_classes += strong_count
        if strong_count:
            histogram_values.add(int(value_text))
    legacy_values = {
        int(value)
        for value in LEGACY_VALUES.read_text(encoding="ascii").split()
    }
    require(histogram_values == legacy_values, "historical order-nine sets agree")
    require((all_classes, strong_classes) == (191_536, 178_133),
            "historical order-nine class totals")
    return histogram_values, all_classes, strong_classes


def small_evaluator_controls(representatives) -> int:
    controls = 0
    for n in range(1, 6):
        for parent_list in representatives[n]:
            parent = tuple(parent_list)
            fast = quadratic_ear_table(parent)
            for signature, predicted in enumerate(fast):
                separate = separate_ear_value(parent, signature)
                child = extend(parent, signature)
                dp_value = hamilton_dp(child)
                literal_value = hamilton_literal(child)
                require(predicted == separate, "small fast/separate ear agreement")
                require(predicted == dp_value, "small fast/direct-DP agreement")
                require(predicted == literal_value, "small fast/literal agreement")
                controls += 1
    return controls


def code_bank_digest(bank: dict[int, tuple[int, ...]]) -> tuple[str, list[tuple[int, int]]]:
    rows = [(value, raw_code(adj)) for value, adj in sorted(bank.items())]
    digest = sha256(
        json.dumps(rows, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    return digest, rows


def build_selected_order9(representatives):
    historical_values, all_classes, strong_classes = historical_order9_values()
    selected = {}
    sources = {}
    strong_order8 = 0
    cuts = 0
    for parent_index, parent_list in enumerate(representatives[8]):
        parent = tuple(parent_list)
        if not is_strong(parent):
            continue
        strong_order8 += 1
        responses = quadratic_ear_table(parent)
        for signature in range(1, (1 << 8) - 1):
            cuts += 1
            value = responses[signature]
            if value not in selected:
                selected[value] = extend(parent, signature)
                sources[value] = (parent_index, signature)
    require(strong_order8 == 6_008, "strong order-eight representative count")
    require(cuts == 1_526_032, "complete order-nine selected-source universe")
    require(len(selected) == 1_482, "selected order-nine bank size")
    require(set(selected) == historical_values, "selected order-nine historical value set")
    digest, rows = code_bank_digest(selected)
    require(digest == EXPECTED_N9_BANK_SHA256, "selected order-nine bank digest")
    failures = sum(
        hamilton_dp(adj) != value or not is_strong(adj)
        for value, adj in selected.items()
    )
    require(failures == 0, "all selected order-nine parents direct DP/strong")
    return selected, sources, rows, strong_order8, cuts, all_classes, strong_classes


def build_selected_order10(selected9):
    selected = {}
    sources = {}
    ears = 0
    for parent_h, parent in sorted(selected9.items()):
        responses = quadratic_ear_table(parent)
        for signature in range(1, (1 << 9) - 1):
            ears += 1
            value = responses[signature]
            if value not in selected:
                selected[value] = extend(parent, signature)
                sources[value] = (parent_h, signature)
    require(ears == 755_820, "complete selected order-ten ear universe")
    require(len(selected) == 7_566, "selected order-ten parent-bank size")
    digest, rows = code_bank_digest(selected)
    require(digest == EXPECTED_N10_BANK_SHA256, "selected order-ten bank digest")
    failures = sum(
        hamilton_dp(adj) != value or not is_strong(adj)
        for value, adj in selected.items()
    )
    require(failures == 0, "all selected order-ten parents direct DP/strong")
    return selected, sources, rows, digest, ears, failures


def multiplicative_prefix(seed: set[int], limit: int):
    attained = {value for value in seed if 1 <= value <= limit}
    product_witnesses = {}
    for value in range(1, limit + 1, 2):
        if value in attained:
            continue
        divisor = 3
        while divisor * divisor <= value:
            if value % divisor == 0:
                quotient = value // divisor
                if divisor in attained and quotient in attained:
                    attained.add(value)
                    product_witnesses[value] = (divisor, quotient)
                    break
            divisor += 2
    return attained, product_witnesses


def main() -> None:
    generator = load_order8_generator()
    representatives, counts = generator.generate(8)
    require(counts == A000568, "A000568 representative counts through order eight")
    small_controls = small_evaluator_controls(representatives)
    require(small_controls == 470, "frozen small evaluator universe")
    print("THM-4104 INDEPENDENT SELECTED-EAR AUDIT")
    print("A000568_counts=", [counts[n] for n in range(1, 9)])
    print("small_fast_separate_DP_literal_controls=", small_controls, "failures= 0")

    (
        selected9,
        selected9_sources,
        selected9_rows,
        strong_order8,
        order9_cuts,
        historical_all,
        historical_strong,
    ) = build_selected_order9(representatives)
    print("selected_order9_parents=", len(selected9),
          "bank_sha256=", EXPECTED_N9_BANK_SHA256)

    selected10, selected10_sources, selected10_rows, selected10_digest, order10_ears, parent_failures = (
        build_selected_order10(selected9)
    )
    print("selected_order10_parents=", len(selected10),
          "bank_sha256=", selected10_digest,
          "direct_DP_strong_failures=", parent_failures)

    parent_items = sorted(selected10.items())
    sample_signatures = {
        index: 1 + ((index * 313 + 17) % 1_022)
        for index in range(0, len(parent_items), 29)
    }
    sample_signatures[len(parent_items) - 1] = 1_022
    sample_fast = {}

    values = set()
    witnesses = {}
    ears_checked = 0
    for parent_index, (parent_h, parent) in enumerate(parent_items):
        responses = quadratic_ear_table(parent)
        if parent_index in sample_signatures:
            signature = sample_signatures[parent_index]
            sample_fast[(parent_index, signature)] = responses[signature]
        for signature in range(1, (1 << 10) - 1):
            ears_checked += 1
            value = responses[signature]
            values.add(value)
            if value not in witnesses:
                child = extend(parent, signature)
                witnesses[value] = (parent_h, signature, raw_code(child))

    require(ears_checked == 7_732_452, "complete selected order-eleven ear universe")
    require(len(values) == 43_251, "selected order-eleven value count")
    require(all(value % 2 == 1 for value in values), "selected order-eleven parity")

    intervals = odd_intervals(values)
    require((429, 80_265, 39_919) in intervals, "primary order-eleven solid interval")
    require((80_875, 84_259, 1_693) in intervals, "secondary order-eleven solid interval")
    adjacent_holes = (427, 80_267, 80_873, 84_261)
    require(all(value not in values for value in adjacent_holes),
            "solid-interval adjacent selected-image holes")
    require(14_657 in values and 14_777 in values,
            "previous ordinary and seven-prime targets are explicit")

    ordinary_missing = first_missing_prime(values, 14_657)
    require(ordinary_missing == 80_407 and 80_407 not in values,
            "next ordinary-prime lane target")
    ordinary_last = previous_prime(ordinary_missing)
    require(ordinary_last in values, "last supplied ordinary prime")

    seven_missing_prime = first_missing_seven_prime(values, 2_111)
    require(seven_missing_prime == 11_527,
            "next seven-prime lane target prime")
    require(7 * seven_missing_prime == 80_689 and 80_689 not in values,
            "next seven-prime lane target value")
    seven_last_prime = previous_prime(seven_missing_prime)
    require(7 * seven_last_prime in values, "last supplied seven-prime value")

    prior_prefix = {
        value for value in range(1, 14_656, 2) if value not in {7, 21}
    }
    attained, product_witnesses = multiplicative_prefix(
        prior_prefix | values | {1}, 80_407
    )
    prefix_missing = [
        value
        for value in range(1, 80_406, 2)
        if value not in {7, 21} and value not in attained
    ]
    require(not prefix_missing, "global allowed prefix through 80405")
    require(80_407 not in attained, "first global unforced ordinary prime")
    bridge_rows = sorted(
        (value, product_witnesses[value])
        for value in product_witnesses
        if value > 80_265
    )

    sample_rows = []
    for (parent_index, signature), fast_value in sorted(sample_fast.items()):
        parent_h, parent = parent_items[parent_index]
        separate_value = separate_ear_value(parent, signature)
        child = extend(parent, signature)
        dp_value = hamilton_dp(child)
        strong = is_strong(child)
        require(fast_value == separate_value, "broad fast/separate sample")
        require(fast_value == dp_value, "broad fast/direct-DP sample")
        require(strong, "broad sampled ear strongness")
        sample_rows.append(
            (parent_index, parent_h, signature, fast_value, raw_code(child), strong)
        )
    require(len(sample_rows) == len(sample_signatures) >= 260,
            "broad deterministic sample size")

    requested_targets = {
        14_657,
        14_777,
        429,
        80_265,
        80_875,
        84_259,
        min(values),
        max(values),
        ordinary_last,
        7 * seven_last_prime,
    }
    if 80_405 in values:
        requested_targets.add(80_405)
    key_rows = {}
    key_parent_cache = {}
    for value in sorted(requested_targets):
        require(value in witnesses, "requested key value has retained witness")
        parent_h, signature, code = witnesses[value]
        parent = selected10[parent_h]
        if parent_h not in key_parent_cache:
            key_parent_cache[parent_h] = quadratic_ear_table(parent)
        fast_value = key_parent_cache[parent_h][signature]
        separate_value = separate_ear_value(parent, signature)
        child = decode_raw(code, 11)
        dp_value = hamilton_dp(child)
        strong = is_strong(child)
        require(raw_code(extend(parent, signature)) == code,
                "key witness code reconstructs from parent and signature")
        require(fast_value == separate_value == dp_value == value,
                "key fast/separate/direct-DP agreement")
        require(strong, "key order-eleven witness strongness")
        row = (parent_h, signature, signature.bit_count(), code)
        key_rows[value] = row
        require(value in EXPECTED_KEY_ROWS and row == EXPECTED_KEY_ROWS[value],
                "frozen key row H=" + str(value))
    require(set(key_rows) == set(EXPECTED_KEY_ROWS), "frozen key target set")

    longest = sorted(intervals, key=lambda row: row[2], reverse=True)[:12]
    semantic_ledger = {
        "input_hashes": [ENGINE_SHA256, HISTOGRAM_SHA256, LEGACY_VALUES_SHA256],
        "A000568": A000568,
        "small_controls": small_controls,
        "historical_classes": [historical_all, historical_strong],
        "strong_order8": strong_order8,
        "order9_cuts": order9_cuts,
        "selected9_rows": selected9_rows,
        "selected9_sources": sorted(selected9_sources.items()),
        "selected10_rows": selected10_rows,
        "selected10_sources": sorted(selected10_sources.items()),
        "selected10_digest": selected10_digest,
        "order10_ears": order10_ears,
        "order11_ears": ears_checked,
        "order11_values": sorted(values),
        "intervals": intervals,
        "adjacent_holes": adjacent_holes,
        "sample_rows": sample_rows,
        "key_rows": key_rows,
        "lane_rows": [
            ordinary_last,
            ordinary_missing,
            seven_last_prime,
            seven_missing_prime,
            7 * seven_missing_prime,
        ],
        "global_prefix": [80_405, bridge_rows],
    }
    semantic = sha256(
        json.dumps(semantic_ledger, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("selected_order11_ears=", ears_checked,
          "strongness_by_nonconstant_ear_lemma=", ears_checked)
    print("selected_value_count=", len(values), "min=", min(values), "max=", max(values))
    print("primary_solid_interval= [429,80265] count=39919 adjacent_holes=(427,80267)")
    print("secondary_solid_interval= [80875,84259] count=1693 adjacent_holes=(80873,84261)")
    print("first_intervals=", intervals[:12])
    print("longest_intervals=", longest)
    print("broad_fast_separate_DP_strong_samples=", len(sample_rows), "failures= 0")
    print("key_boundary_witnesses=", key_rows)
    print("explicit_previous_targets=", (14_657, key_rows[14_657]),
          (14_777, key_rows[14_777]))
    print("global_allowed_prefix_through=80405 multiplicative_bridges_above_80265=",
          bridge_rows)
    print("ordinary_prime_lane_through=", ordinary_last,
          "next_unforced=", ordinary_missing)
    print("seven_prime_lane_p_through=", seven_last_prime,
          "next_unforced_value=", 7 * seven_missing_prime)
    print("semantic_sha256=", semantic)
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()
