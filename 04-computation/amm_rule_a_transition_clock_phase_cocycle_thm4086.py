#!/usr/bin/env python3
"""Primary exact companion for THM-4086.

This program proves a rational enclosure for alpha=log_5(phi^2), compiles
every needed floor by integer division inside that enclosure, reconstructs
the autonomous Rule-A cell-zero clocks at R=32768, checks a small exact
nonmonotonicity control, and audits the two-coordinate dyadic profile
cocycle.  It imports no historical engine and uses no floating point.
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
HALF_R = 16_384
ALPHA_LOWER = (10_756, 17_987)
ALPHA_UPPER = (105_183, 175_895)
PROFILE_MAX = 2 * R - 1

ENGINE_SHA256 = "8887080fc6e30760efa4a0ba76218ec97676cc717c6e76ccefbaeec6c73684ad"
CHECKPOINT_SHA256 = "a3c1c856979fd7b0b21ca12f93d08b3340345fd175187b9c724a14bc91b97710"
S192_RUNNER_SHA256 = "65b6fab453a0adb156bc79e1f9168fd46d95037fd345e2322861ec1f76f15ba0"
S193_RUNNER_SHA256 = "962577096266bdc6df7b004eb38412b85a8a2b47eb9dd000f4666e9829043207"
THM3644_RUNNER_SHA256 = "ddf248e8939bbdefa7b2544bbd3df1c23e47e53e00aae4013d09149eae2ca59c"
THM3644_OUTPUT_SHA256 = "28cedef2dc0b176b62f0633d93ea23a7c316993ed197743bd54f4466eb860c21"

# D0,status,event,d_start,d_event,d_last,minus2,T6b,overflow_bits.
STABLE_ROWS = (
    (854, "DIE", 8_246, 20_448, 25_379, 40_043, 4_394, 8_130, 19_875),
    (855, "CLOSED", 20_238, 20_449, 32_551, 40_044, 16_383, 8_132, None),
    (856, "CLOSED", 20_233, 20_450, 32_549, 40_045, 16_383, 8_134, None),
)

# D0,last feed row,debt after last feed,first feed-free row,scalar zero row.
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
EXPECTED_ETA_COUNTS = ((0, 8_191), (1, 8_193))
EXPECTED_ZETA_COUNTS = ((0, 3_294), (1, 8_190), (2, 4_900))
EXPECTED_COMMON_SHA256 = "079764dd469362281f1374d6f2a9db81652a02af0ec827420dae71fe774d20fd"
EXPECTED_PRIMARY_SHA256 = "2b49218c35d0c2793dbcca51e9d1f66c2d07f70ff494a46b5266a849dd728259"

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
    """Decide 5^power <= phi^(2*exponent_half) over exact integers."""
    if power < 0:
        return True
    fibonacci, fibonacci_next = fib_pair(2 * exponent_half)
    lucas = 2 * fibonacci_next - fibonacci
    delta = 2 * pow(5, power) - lucas
    return delta <= 0 or delta * delta < 5 * fibonacci * fibonacci


def compile_floors() -> tuple[int, ...]:
    lower_p, lower_q = ALPHA_LOWER
    upper_p, upper_q = ALPHA_UPPER
    require(five_power_le_phi_even(lower_p, lower_q), "lower alpha enclosure")
    require(not five_power_le_phi_even(upper_p, upper_q), "upper alpha enclosure")
    floors = []
    for n in range(PROFILE_MAX + 1):
        lower_floor = lower_p * n // lower_q
        upper_floor = upper_p * n // upper_q
        require(lower_floor == upper_floor, ("floor enclosure", n))
        floors.append(lower_floor)
    return tuple(floors)


def clamp(value: int, lower: int, upper: int) -> int:
    return min(upper, max(lower, value))


def initial_cell_zero_junk(r: int, degree: int) -> int:
    signed_binomial = comb(r - 2, degree)
    if degree & 1:
        signed_binomial = -signed_binomial
    load = signed_binomial - (degree + 1) + 2
    correction = clamp(load, -2, 0)
    require((load - correction) % 2 == 0, "initial cell-zero parity")
    return load - correction


def cell_zero_clock(floors: tuple[int, ...], offset: int) -> tuple[object, ...]:
    """Reconstruct cell zero using forward consecutive-binomial streaming."""
    degree_previous = floors[R] + offset
    junk = initial_cell_zero_junk(R, degree_previous)
    address = degree_previous
    binomial = comb(R - 1, address)
    address_blocks = []
    last_feed = None
    debt_after_last_feed = None
    first_feed_free = None
    scalar_zero = None
    previous_feed_free_debt = None

    for row in range(1, R - 1):
        degree = floors[R + row] + offset
        epsilon = degree - degree_previous
        require(epsilon in (0, 1), ("degree step", offset, row, epsilon))
        target_address = degree + row
        feed = 0
        consumed = []
        while address < min(target_address, R - 1):
            numerator = binomial * (R - 1 - address)
            denominator = address + 1
            binomial, remainder = divmod(numerator, denominator)
            require(remainder == 0, ("forward binomial", offset, row))
            address += 1
            coefficient = (-binomial if address & 1 else binomial) - 1
            feed += coefficient
            consumed.append(address)

        expected_count = 0
        if target_address <= R - 1:
            expected_count += 1
        if epsilon == 1 and target_address - 1 <= R - 1:
            expected_count += 1
        require(len(consumed) == expected_count, ("feed address block", offset, row))
        if consumed:
            require(consumed == list(range(consumed[0], consumed[-1] + 1)),
                    ("consecutive feed", offset, row))
            address_blocks.append((row, tuple(consumed)))
            last_feed = row
        elif first_feed_free is None:
            first_feed_free = row

        load = junk + feed
        correction = clamp(load, -2, 0)
        require((load - correction) % 2 == 0, ("cell-zero parity", offset, row))
        junk = load - correction
        debt = -junk

        if consumed:
            debt_after_last_feed = debt
            previous_feed_free_debt = None
        else:
            if previous_feed_free_debt is not None and previous_feed_free_debt > 0:
                require(debt == max(0, previous_feed_free_debt - 2),
                        ("two-per-row clock", offset, row))
            previous_feed_free_debt = debt
            if junk == 0:
                scalar_zero = row
                break
        degree_previous = degree

    require(last_feed is not None, ("last feed", offset))
    require(first_feed_free == last_feed + 1, ("feed boundary", offset))
    require(debt_after_last_feed is not None and debt_after_last_feed > 0,
            ("positive terminal debt", offset))
    require(debt_after_last_feed % 2 == 0, ("even terminal debt", offset))
    require(scalar_zero == last_feed + debt_after_last_feed // 2,
            ("clock formula", offset))
    address_rows = tuple((row, addresses) for row, addresses in address_blocks)
    return (offset, last_feed, debt_after_last_feed, first_feed_free, scalar_zero,
            canonical_digest(address_rows))


def initial_junk_small(r: int, degree: int) -> dict[int, int]:
    a_term = comb(r - 2, degree)
    b_term = degree + 1
    lower_binomial = 1
    lower_previous = 0
    degree_binomial = 1
    sign = 1 if degree % 2 == 0 else -1
    junk = {}
    for cell in range(degree + 1):
        load = sign * a_term - b_term + 2 * degree_binomial
        correction = clamp(load, -2 * lower_binomial, 2 * lower_previous)
        require((load - correction) % 2 == 0, ("small initial parity", cell))
        if correction != load:
            junk[cell] = load - correction
        if cell < degree:
            a_term, rem_a = divmod(a_term * (degree - cell), r - 2 - cell)
            b_term, rem_b = divmod(b_term * (degree - cell), cell + 2)
            degree_binomial, rem_d = divmod(
                degree_binomial * (degree - cell), cell + 1)
            require(rem_a == rem_b == rem_d == 0, "small initial divisions")
            lower_previous = lower_binomial
            if cell < degree - 1:
                lower_binomial, rem_l = divmod(
                    lower_binomial * (degree - 1 - cell), cell + 1)
                require(rem_l == 0, "small lower division")
            else:
                lower_binomial = 0
            sign = -sign
    return junk


def two_g_small(r: int) -> tuple[int, ...]:
    coefficients = [2]
    binomial = 1
    for index in range(1, r):
        binomial, remainder = divmod(binomial * (r - index), index)
        require(remainder == 0, "small feed division")
        coefficients.append((-binomial if index & 1 else binomial) - 1)
    return tuple(coefficients)


def run_small_rule_a(floors: tuple[int, ...], r: int, offset: int) -> tuple[str, int]:
    feed_coefficients = two_g_small(r)
    degree_previous = floors[r] + offset
    junk = initial_junk_small(r, degree_previous)
    for row in range(1, r - 1):
        degree = floors[r + row] + offset
        epsilon = degree - degree_previous
        require(epsilon in (0, 1), ("small degree step", row))
        kernel = (1, 1) if epsilon == 0 else (1, 2, 1)
        load = {}
        for cell, value in junk.items():
            for shift, multiplier in enumerate(kernel):
                load[cell + shift] = load.get(cell + shift, 0) + multiplier * value
        if degree + row <= r - 1:
            load[0] = load.get(0, 0) + feed_coefficients[degree + row]
        if epsilon == 1 and degree - 1 + row <= r - 1:
            coefficient = feed_coefficients[degree - 1 + row]
            load[0] = load.get(0, 0) + coefficient
            load[1] = load.get(1, 0) + coefficient
        next_junk = {}
        for cell, value in load.items():
            lower = -2 * comb(degree - 1, cell)
            upper = 2 * comb(degree - 1, cell - 1) if cell else 0
            correction = clamp(value, lower, upper)
            require((value - correction) % 2 == 0, ("small row parity", row, cell))
            if correction != value:
                next_junk[cell] = value - correction
        if degree in next_junk:
            return "DIE", row
        junk = next_junk
        degree_previous = degree
        if not junk and degree + row > r - 1:
            return "CLOSED", row
    raise RuntimeError(("small terminal not reached", r, offset))


def count_values(values: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple((value, values.count(value)) for value in sorted(set(values)))


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    audit_trusted_status_inputs(root)
    floors = compile_floors()
    require(floors[R] == 19_594, "R floor")
    require(floors[2 * R - 1] == 39_189, "last floor")

    epsilon = tuple(floors[R + row] - floors[R + row - 1]
                    for row in range(1, R))
    require(set(epsilon) == {0, 1}, "binary kernel word")
    kernel_record = {
        "length": len(epsilon),
        "zero_steps": epsilon.count(0),
        "one_steps": epsilon.count(1),
        "sha256": canonical_digest(epsilon),
    }

    for stable in STABLE_ROWS:
        offset, _status, event, degree_start, degree_event, degree_last = stable[:6]
        require(degree_start == floors[R] + offset, ("stable start", offset))
        require(degree_event == floors[R + event] + offset, ("stable event", offset))
        require(degree_last == floors[2 * R - 1] + offset, ("stable last", offset))

    front_rows = []
    for stable in STABLE_ROWS:
        offset, _status, _event, degree_start = stable[:4]
        front = R - 2 - degree_start
        headroom = degree_start - front
        require(headroom == stable[7], ("T6b headroom", offset))
        front_rows.append((offset, degree_start, front, headroom))
    front_rows = tuple(front_rows)

    clock_audits = tuple(cell_zero_clock(floors, offset) for offset in (854, 855, 856))
    clock_rows = tuple(row[:5] for row in clock_audits)
    require(clock_rows == EXPECTED_CLOCK_ROWS, ("clock rows", clock_rows))
    require(clock_rows[0][4] < clock_rows[1][4] > clock_rows[2][4],
            "strict local clock peak")
    require(STABLE_ROWS[0][2] < clock_rows[0][4], "death precedes scalar zero")
    require(STABLE_ROWS[1][2] == clock_rows[1][4], "855 capture clock")
    require(STABLE_ROWS[2][2] == clock_rows[2][4], "856 capture clock")

    small_terminals = tuple((offset,) + run_small_rule_a(floors, 16, offset)
                            for offset in range(6))
    require(small_terminals == EXPECTED_SMALL_TERMINALS,
            ("small terminals", small_terminals))
    require(small_terminals[0][2] < small_terminals[1][2] > small_terminals[2][2],
            "small capture nonmonotonicity")

    eta = tuple(floors[2 * (HALF_R + row)] - 2 * floors[HALF_R + row]
                for row in range(HALF_R))
    zeta = tuple(floors[2 * (HALF_R + row) + 1] - 2 * floors[HALF_R + row]
                 for row in range(HALF_R))
    eta_counts = count_values(eta)
    zeta_counts = count_values(zeta)
    require(eta_counts == EXPECTED_ETA_COUNTS, ("eta counts", eta_counts))
    require(zeta_counts == EXPECTED_ZETA_COUNTS, ("zeta counts", zeta_counts))
    eta_witnesses = ((0, eta[0]), (2, eta[2]))
    zeta_witnesses = ((0, zeta[0]), (1, zeta[1]), (4, zeta[4]))
    require(eta_witnesses == ((0, 0), (2, 1)), "eta witnesses")
    require(zeta_witnesses == ((0, 1), (1, 0), (4, 2)), "zeta witnesses")

    triangle = {
        "R16384": ((412, "DIE", 4_116), (413, "CLOSED", 10_116)),
        "R32768": tuple((row[0], row[1], row[2]) for row in STABLE_ROWS),
        "outer_offset_map": ((412, 854), (413, 856)),
        "even_profile_defect": "30+eta",
    }
    require(triangle["outer_offset_map"] == ((412, 2 * 412 + 30),
                                               (413, 2 * 413 + 30)),
            "scale offset map")

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

    primary_semantic = {
        "common_sha256": common_sha256,
        "alpha_enclosure": (ALPHA_LOWER, ALPHA_UPPER, PROFILE_MAX),
        "kernel_record": kernel_record,
        "clock_address_sha256": tuple(row[5] for row in clock_audits),
        "eta": (eta_counts, canonical_digest(eta)),
        "zeta": (zeta_counts, canonical_digest(zeta)),
    }
    primary_sha256 = canonical_digest(primary_semantic)
    require(primary_sha256 == EXPECTED_PRIMARY_SHA256,
            ("primary semantic digest", primary_sha256))

    source = Path(__file__).resolve()
    source_bytes = source.read_bytes()
    require(b"\r\n" not in source_bytes, "source LF")
    syntax = ast.parse(source_bytes.decode("utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)),
            "assert node present")
    require(not any(isinstance(node, ast.Constant) and isinstance(node.value, float)
                    for node in ast.walk(syntax)), "float constant present")

    print("== THM-4086 Rule-A transition clock and phase cocycle: primary ==")
    print("scope=PROVED structural formulas + FINITE-EXACT R=32768 transition facts")
    print("trusted_status_input=INHERITED-FINITE-EXACT;short_companion_does_not_replay_full_state")
    print(f"archive_pins=engine:{ENGINE_SHA256};S192:{S192_RUNNER_SHA256};S193:{S193_RUNNER_SHA256};checkpoint:{CHECKPOINT_SHA256}")
    print(f"R16384_pins=runner:{THM3644_RUNNER_SHA256};output:{THM3644_OUTPUT_SHA256}")
    print("alpha_enclosure=10756/17987<=alpha<105183/175895;floors_exact:n=0..65535")
    print("kernel_record=" + json.dumps(kernel_record, sort_keys=True, separators=(",", ":")))
    print(f"front_rows={front_rows}")
    print(f"stable_rows={STABLE_ROWS}")
    print(f"cell_zero_clocks={clock_rows}")
    print("clock_shape=20234<20238>20233;855_to_856_capture_shift=-5")
    print(f"small_R16_terminals={small_terminals};neither_nondecreasing_nor_nonincreasing=PASS")
    print(f"eta_counts={eta_counts};eta_witnesses={eta_witnesses};eta_sha256={canonical_digest(eta)}")
    print(f"zeta_counts={zeta_counts};zeta_witnesses={zeta_witnesses};zeta_sha256={canonical_digest(zeta)}")
    print("scale_triangle=R16384:412D/413C;R32768:854D/855C/856C;outer_defect=30+eta")
    print("lift_obstruction=RuleA_856_not_equal_rowwise_degree_lift_of_RuleA_855;capture_rows:20233!=20238")
    print(f"common_semantic_sha256={common_sha256}")
    print(f"primary_semantic_sha256={primary_sha256}")
    print(f"CHECKS={CHECKS}")
    print("status=PROVED-STRUCTURAL+FINITE-EXACT+VERIFIED-EXACT;fixed Rule-A policy")
    print("nonconsequence=no global offset threshold,alternative infeasibility,uniform family,D0 law,or C-star bound")


if __name__ == "__main__":
    main()
