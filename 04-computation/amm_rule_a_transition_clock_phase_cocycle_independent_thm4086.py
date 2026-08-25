#!/usr/bin/env python3
"""Structurally independent exact audit for THM-4086.

Unlike the primary enclosure compiler, this audit decides the load-bearing
R=32768 floor prefix directly by Fibonacci--Lucas comparisons.  It assigns
feed addresses first and reconstructs their binomial coefficients backwards
from C(R-1,R-1)=1.  Its R=16 control uses dense state arrays.  No historical
engine, primary module, floating point, or Python assert is used.
"""

from __future__ import annotations

import ast
from functools import lru_cache
from hashlib import sha256
import json
from math import comb
from pathlib import Path
import re


R = 32_768
OFFSETS = (854, 855, 856)
ENGINE_SHA256 = "8887080fc6e30760efa4a0ba76218ec97676cc717c6e76ccefbaeec6c73684ad"
CHECKPOINT_SHA256 = "a3c1c856979fd7b0b21ca12f93d08b3340345fd175187b9c724a14bc91b97710"
S192_RUNNER_SHA256 = "65b6fab453a0adb156bc79e1f9168fd46d95037fd345e2322861ec1f76f15ba0"
S193_RUNNER_SHA256 = "962577096266bdc6df7b004eb38412b85a8a2b47eb9dd000f4666e9829043207"
THM3644_RUNNER_SHA256 = "ddf248e8939bbdefa7b2544bbd3df1c23e47e53e00aae4013d09149eae2ca59c"
THM3644_OUTPUT_SHA256 = "28cedef2dc0b176b62f0633d93ea23a7c316993ed197743bd54f4466eb860c21"
STABLE_ROWS = (
    (854, "DIE", 8_246, 20_448, 25_379, 40_043, 4_394, 8_130, 19_875),
    (855, "CLOSED", 20_238, 20_449, 32_551, 40_044, 16_383, 8_132, None),
    (856, "CLOSED", 20_233, 20_450, 32_549, 40_045, 16_383, 8_134, None),
)
EXPECTED_CLOCK_ROWS = (
    (854, 7_709, 25_050, 7_710, 20_234),
    (855, 7_708, 25_060, 7_709, 20_238),
    (856, 7_708, 25_050, 7_709, 20_233),
)
EXPECTED_SMALL_TERMINALS = (
    (0, "CLOSED", 7),
    (1, "CLOSED", 8),
    (2, "CLOSED", 6),
    (3, "CLOSED", 7),
    (4, "CLOSED", 6),
    (5, "CLOSED", 6),
)
EXPECTED_COMMON_SHA256 = "079764dd469362281f1374d6f2a9db81652a02af0ec827420dae71fe774d20fd"
EXPECTED_INDEPENDENT_SHA256 = "1667d631e30f37800eabb1ad472113e931fe0f05beb233a78dc200a90f052b61"

CHECKS = 0


def require(condition: bool, payload: object) -> None:
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def canonical_digest(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(encoded).hexdigest()


def raw_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def audit_trusted_status_inputs(root: Path) -> None:
    paths = (
        (root / "04-computation/amm12592_transient_fast_junkflow_boxeph.py",
         ENGINE_SHA256),
        (root / "04-computation/amm12592_r32768_offset_boundary_probe_kps_s192.py",
         S192_RUNNER_SHA256),
        (root / "04-computation/amm12592_r32768_offset_boundary_lowmemory_probe_kps_s193.py",
         S193_RUNNER_SHA256),
        (root / "04-computation/amm12592_R16384_offset_transition_thm3644.py",
         THM3644_RUNNER_SHA256),
        (root / "05-knowledge/results/amm12592_R16384_offset_transition_thm3644.out",
         THM3644_OUTPUT_SHA256),
    )
    for path, expected in paths:
        require(raw_sha256(path) == expected, ("trusted input drift", str(path)))
    checkpoint = root / "05-knowledge/results/amm12592_r32768_offset855_live_checkpoint_kps_s192.out"
    require(raw_sha256(checkpoint) == CHECKPOINT_SHA256, "checkpoint drift")
    text = checkpoint.read_text(encoding="utf-8")
    parsed = tuple(ast.literal_eval(item) for item in
                   re.findall(r"stable_result=(\([^;]+\))", text))
    require(tuple(sorted(parsed)) == STABLE_ROWS, ("checkpoint stable rows", parsed))


@lru_cache(maxsize=None)
def fib_pair(n: int) -> tuple[int, int]:
    if n == 0:
        return 0, 1
    a, b = fib_pair(n >> 1)
    c = a * (2 * b - a)
    d = a * a + b * b
    if n & 1:
        return d, c + d
    return c, d


def five_power_le_phi_even(power: int, exponent_half: int) -> bool:
    if power < 0:
        return True
    fibonacci, fibonacci_next = fib_pair(2 * exponent_half)
    lucas = 2 * fibonacci_next - fibonacci
    delta = 2 * pow(5, power) - lucas
    return delta <= 0 or delta * delta < 5 * fibonacci * fibonacci


def direct_floor(n: int) -> int:
    lower = 0
    upper = n
    while lower < upper:
        middle = (lower + upper + 1) // 2
        if five_power_le_phi_even(middle, n):
            lower = middle
        else:
            upper = middle - 1
    return lower


def direct_feed_prefix() -> tuple[int, ...]:
    """Compile through the first row feed-free for all three offsets."""
    floors = [direct_floor(R)]
    require(floors[0] == 19_594, "direct R floor")
    for row in range(1, R - 1):
        candidate = floors[-1] + 1
        if five_power_le_phi_even(candidate, R + row):
            next_floor = candidate
        else:
            next_floor = candidate - 1
        require(next_floor - floors[-1] in (0, 1), ("direct step", row))
        floors.append(next_floor)
        all_feed_free = True
        for offset in OFFSETS:
            degree = next_floor + offset
            epsilon = next_floor - floors[-2]
            address = degree + row
            fed = address <= R - 1 or (epsilon == 1 and address - 1 <= R - 1)
            all_feed_free = all_feed_free and not fed
        if all_feed_free:
            break
    require(len(floors) == 7_711, ("direct prefix length", len(floors)))
    return tuple(floors)


def address_assignment(base_floors: tuple[int, ...], offset: int) -> tuple[dict[int, int], int, int]:
    assignment = {}
    last_feed = None
    first_feed_free = None
    previous_address = base_floors[0] + offset
    for row in range(1, len(base_floors)):
        degree = base_floors[row] + offset
        epsilon = base_floors[row] - base_floors[row - 1]
        address = degree + row
        consumed = tuple(range(previous_address + 1, min(address, R - 1) + 1))
        expected_count = 0
        if address <= R - 1:
            expected_count += 1
        if epsilon == 1 and address - 1 <= R - 1:
            expected_count += 1
        require(len(consumed) == expected_count, ("assigned block", offset, row))
        if consumed:
            for index in consumed:
                require(index not in assignment, ("duplicate address", offset, index))
                assignment[index] = row
            last_feed = row
        elif first_feed_free is None:
            first_feed_free = row
        previous_address = address
    require(last_feed is not None and first_feed_free == last_feed + 1,
            ("assigned boundary", offset))
    require(tuple(sorted(assignment)) == tuple(range(min(assignment), R)),
            ("complete feed interval", offset))
    return assignment, last_feed, first_feed_free


def reverse_feed_by_row(assignment: dict[int, int]) -> dict[int, int]:
    """Build C(R-1,k) downwards and aggregate directly into assigned rows."""
    minimum = min(assignment)
    binomial = 1
    feeds = {}
    for index in range(R - 1, minimum - 1, -1):
        if index in assignment:
            coefficient = (-binomial if index & 1 else binomial) - 1
            row = assignment[index]
            feeds[row] = feeds.get(row, 0) + coefficient
        if index > minimum:
            numerator = binomial * index
            denominator = R - index
            binomial, remainder = divmod(numerator, denominator)
            require(remainder == 0, ("reverse binomial", index))
    return feeds


def clamp(value: int, lower: int, upper: int) -> int:
    return min(upper, max(lower, value))


def independent_clock(base_floors: tuple[int, ...], offset: int) -> tuple[int, ...]:
    assignment, last_feed, first_feed_free = address_assignment(base_floors, offset)
    feeds = reverse_feed_by_row(assignment)
    degree_start = base_floors[0] + offset
    initial_binomial = comb(R - 2, degree_start)
    if degree_start & 1:
        initial_binomial = -initial_binomial
    initial_load = initial_binomial - (degree_start + 1) + 2
    junk = initial_load - clamp(initial_load, -2, 0)
    debt_after_last_feed = None
    for row in range(1, first_feed_free + 1):
        load = junk + feeds.get(row, 0)
        correction = clamp(load, -2, 0)
        require((load - correction) % 2 == 0, ("independent parity", offset, row))
        junk = load - correction
        if row == last_feed:
            debt_after_last_feed = -junk
    require(debt_after_last_feed is not None and debt_after_last_feed > 0,
            ("independent debt", offset))
    require(debt_after_last_feed % 2 == 0, ("independent even debt", offset))
    require(-junk == debt_after_last_feed - 2,
            ("independent first drain", offset))
    scalar_zero = last_feed + debt_after_last_feed // 2
    return offset, last_feed, debt_after_last_feed, first_feed_free, scalar_zero


def dense_initial_junk(r: int, degree: int) -> list[int]:
    row = []
    for cell in range(degree + 1):
        a_term = comb(r - 2 - cell, degree - cell)
        signed = -a_term if (degree - cell) & 1 else a_term
        load = signed - comb(degree + 1, cell + 1) + 2 * comb(degree, cell)
        lower = -2 * comb(degree - 1, cell)
        upper = 2 * comb(degree - 1, cell - 1) if cell else 0
        correction = clamp(load, lower, upper)
        row.append(load - correction)
    return row


def dense_small_rule_a(r: int, offset: int) -> tuple[str, int]:
    floors = tuple(direct_floor(r + row) for row in range(r))
    feed = []
    for index in range(r):
        if index == 0:
            feed.append(2)
        else:
            coefficient = comb(r - 1, index)
            feed.append((-coefficient if index & 1 else coefficient) - 1)
    degree_previous = floors[0] + offset
    junk = dense_initial_junk(r, degree_previous)
    for row in range(1, r - 1):
        degree = floors[row] + offset
        epsilon = degree - degree_previous
        kernel = (1, 1) if epsilon == 0 else (1, 2, 1)
        load = [0] * (degree + 1)
        for cell, value in enumerate(junk):
            for shift, multiplier in enumerate(kernel):
                target = cell + shift
                if target <= degree:
                    load[target] += multiplier * value
                else:
                    require(value == 0, ("dense support beyond degree", row, target))
        if degree + row <= r - 1:
            load[0] += feed[degree + row]
        if epsilon == 1 and degree - 1 + row <= r - 1:
            coefficient = feed[degree - 1 + row]
            load[0] += coefficient
            load[1] += coefficient
        next_junk = []
        for cell, value in enumerate(load):
            lower = -2 * comb(degree - 1, cell)
            upper = 2 * comb(degree - 1, cell - 1) if cell else 0
            next_junk.append(value - clamp(value, lower, upper))
        if next_junk[degree] != 0:
            return "DIE", row
        junk = next_junk
        degree_previous = degree
        if not any(junk) and degree + row > r - 1:
            return "CLOSED", row
    raise RuntimeError(("dense terminal not reached", r, offset))


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    audit_trusted_status_inputs(root)
    base_floors = direct_feed_prefix()
    clock_rows = tuple(independent_clock(base_floors, offset) for offset in OFFSETS)
    require(clock_rows == EXPECTED_CLOCK_ROWS, ("independent clock rows", clock_rows))

    front_rows = []
    for stable in STABLE_ROWS:
        offset, _status, _event, degree_start = stable[:4]
        require(degree_start == base_floors[0] + offset, ("independent start", offset))
        front = R - 2 - degree_start
        headroom = degree_start - front
        require(headroom == stable[7], ("independent headroom", offset))
        front_rows.append((offset, degree_start, front, headroom))
    front_rows = tuple(front_rows)

    small_terminals = tuple((offset,) + dense_small_rule_a(16, offset)
                            for offset in range(6))
    require(small_terminals == EXPECTED_SMALL_TERMINALS,
            ("independent small terminals", small_terminals))

    eta_witnesses = tuple((row,
                           direct_floor(2 * (16_384 + row))
                           - 2 * direct_floor(16_384 + row))
                          for row in (0, 2))
    zeta_witnesses = tuple((row,
                            direct_floor(2 * (16_384 + row) + 1)
                            - 2 * direct_floor(16_384 + row))
                           for row in (0, 1, 4))
    require(eta_witnesses == ((0, 0), (2, 1)), "independent eta witnesses")
    require(zeta_witnesses == ((0, 1), (1, 0), (4, 2)),
            "independent zeta witnesses")

    triangle = {
        "R16384": ((412, "DIE", 4_116), (413, "CLOSED", 10_116)),
        "R32768": tuple((row[0], row[1], row[2]) for row in STABLE_ROWS),
        "outer_offset_map": ((412, 854), (413, 856)),
        "even_profile_defect": "30+eta",
    }
    common_semantic = {
        "archive": (ENGINE_SHA256, CHECKPOINT_SHA256, S192_RUNNER_SHA256,
                    S193_RUNNER_SHA256, THM3644_RUNNER_SHA256,
                    THM3644_OUTPUT_SHA256),
        "stable_rows": STABLE_ROWS,
        "front_rows": front_rows,
        "clock_rows": clock_rows,
        "small_terminals": small_terminals,
        "eta_witnesses": eta_witnesses,
        "zeta_witnesses": zeta_witnesses,
        "triangle": triangle,
    }
    common_sha256 = canonical_digest(common_semantic)
    require(common_sha256 == EXPECTED_COMMON_SHA256,
            ("common semantic digest", common_sha256))
    independent_semantic = {
        "common_sha256": common_sha256,
        "direct_prefix_length": len(base_floors),
        "direct_prefix_sha256": canonical_digest(base_floors),
        "routes": ("Fibonacci-Lucas direct successor decisions",
                   "reverse binomial feed aggregation", "dense small engine"),
    }
    independent_sha256 = canonical_digest(independent_semantic)
    require(independent_sha256 == EXPECTED_INDEPENDENT_SHA256,
            ("independent semantic digest", independent_sha256))

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source LF")
    syntax = ast.parse(source_bytes.decode("utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)),
            "assert node present")
    require(not any(isinstance(node, ast.Constant) and isinstance(node.value, float)
                    for node in ast.walk(syntax)), "float constant present")

    print("== THM-4086 Rule-A transition clock and phase cocycle: independent audit ==")
    print("trusted_status_input=INHERITED-FINITE-EXACT;short_referee_rebuilds_clock_not_full_state")
    print(f"archive_pins=engine:{ENGINE_SHA256};S192:{S192_RUNNER_SHA256};S193:{S193_RUNNER_SHA256};checkpoint:{CHECKPOINT_SHA256}")
    print(f"R16384_pins=runner:{THM3644_RUNNER_SHA256};output:{THM3644_OUTPUT_SHA256}")
    print("route=direct Fibonacci-Lucas floor prefix + reverse binomial feed + dense R16 state")
    print(f"direct_prefix=rows:{len(base_floors)};sha256:{canonical_digest(base_floors)}")
    print(f"front_rows={front_rows}")
    print(f"stable_rows={STABLE_ROWS}")
    print(f"cell_zero_clocks={clock_rows}")
    print(f"small_R16_terminals={small_terminals};two_direction_hostile_control=PASS")
    print(f"eta_witnesses={eta_witnesses};zeta_witnesses={zeta_witnesses}")
    print("transition_split=854:front_death_before_clock;855,856:clock_exactly_equals_capture")
    print(f"common_semantic_sha256={common_sha256}")
    print(f"independent_semantic_sha256={independent_sha256}")
    print(f"CHECKS={CHECKS}")
    print("status=INDEPENDENTLY-VERIFIED-EXACT;load-bearing transition-clock fields")
    print("scope=fixed Rule-A policy;no alternative feasibility,global monotonicity,or C-star consequence")


if __name__ == "__main__":
    main()
