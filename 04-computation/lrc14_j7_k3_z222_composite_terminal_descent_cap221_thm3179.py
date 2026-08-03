#!/usr/bin/env python3
"""Exact projected-k3 z222 composite screen and terminal descent to cap221."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
from collections import Counter
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE_3174 = ROOT / (
    "04-computation/"
    "lrc14_j7_k3_z223_terminal_descent_cap222_thm3174.py"
)
OUTPUT_3174 = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z223_terminal_descent_cap222_thm3174.out"
)
OUTPUT = ROOT / (
    "05-knowledge/results/"
    "lrc14_j7_k3_z222_composite_terminal_descent_cap221_thm3179.out"
)

SOURCE_3174_SHA256 = "3ac95e58861078828671e606e1c97705ac4def2728019ac21bd78eba0b9f1c18"
OUTPUT_3174_SHA256 = "4edc6c7efb64ce1062aa863a5f4c17e4735f42bf7b4ac371c6bca660346abe43"
SEMANTIC_3174_SHA256 = "f6493459d78f6af58b999bd54c8b72e2a8f989c533a1e5a91d5e6bc337352ec7"

LEVEL = 222
EXPECTED_CENSUS = (219, 199, 20)
EXPECTED_ROW_ORDER_SHA256 = "68ea5b7128d3deafd4bb9a14f6d1d09ba867a5f117645a3e1f4f3d03e761de73"
EXPECTED_STRATA = (((False, 6), 20), ((True, 2), 1), ((True, 6), 198))
EXPECTED_EXCEPTION = (8, (1, 2, 4, 8, 10, 14), 3920, 387, 1960)
EXPECTED_SCREEN = (49676, 27341, 21285, 1050)
EXPECTED_FARKAS = (0, 21285)
EXPECTED_SCREEN_SHA256 = "146bae48b8a901cbcc3614dcba88cbeaa87f0361ef0b2a1d9b48d732628a09b5"
EXPECTED_RESIDUAL_SHA256 = "2ad57452ad83947b33e3245004dc44be6f7846d6870c2b5a9a3854f18c3b92b5"
EXPECTED_RESIDUAL_INDEX_SHA256 = "319c0ab23c866f9ccae57f82768e52231d18b1688c6fd8754eeab0f2f5e57939"
EXPECTED_DIVISOR_STATES = (((6, 3), 38), ((6, 6), 1012))
EXPECTED_COUPLING_STATES = (((3, False), 38), ((6, False), 700), ((6, True), 312))
EXPECTED_TERMINAL = (37, 1050, 984, 1150, 161, 1075, 75, 0, 0, 1, True, True)
EXPECTED_TERMINAL_SEMANTIC_SHA256 = "cfa99f5f95e2e241da6cd2a91a199df21c738645cd7a0689d67dfe9926b26e17"
EXPECTED_CLOSURE_SHA256 = "dc9fd29bc9c1c8525e12913294d312f38f58cdedbac4ffe013bfdddb8d20840c"
EXPECTED_CASE_VECTOR_SHA256 = "cf59598dc5966b961857567d81ea55e3bf0e494ea6af688e0b904ac8142a70c6"
EXPECTED_SEMANTIC_SHA256 = "8f098a0c1a00f242bc6f3dd2bfa4932a6aff187186d475d889509f3a1b380ad5"

LEDGER_BEFORE = 374025
LAYER_ROWS = 219
LEDGER_AFTER = 373806
NEXT_CAP = 221
NEXT_LEVEL = 221
NEXT_CENSUS = (90, 83, 7)
NEXT_ROW_ORDER_SHA256 = "ab4f5b7fe3b5e1330d6e51dfb150527fa27cac56b755ed5e2ae2522f0f27ceb4"

RESIDUAL = {
    (1, 2, 9, 10, 12, 14): (7, "02a6c8a376cc3e70711d8f63c2bff0335a71a8d7259f8f550f6b7911e1ab8cea"),
    (1, 2, 9, 11, 12, 14): (1, "66c262f9c3b1298ff92a0b97d5b23331bf64618979d1d6f58930c8b21d65cdd1"),
    (1, 2, 10, 11, 12, 14): (24, "31c83399e99d5a5d047fad3e3ce64bfc3d37f26c8661679b9b99387deb691166"),
    (1, 4, 5, 9, 12, 14): (1, "5e7c97e56e12fda1553c873c96d6957788731caaca5ab4c93b6868d9c77303e2"),
    (1, 4, 9, 10, 12, 14): (4, "c576071739d5b2d3439a0676d11f038508da2e25787504b1965bbf8d1413c570"),
    (1, 5, 6, 7, 9, 12): (8, "3766089b44ab8148df74def4d840435cf2965799c0f4368d1d289a00847dc2c5"),
    (1, 5, 6, 9, 12, 14): (1, "2ec261fe63eca5be6ef8f0bdaa3655ddac941941a7702ca7cbe20a6b2c278a39"),
    (1, 5, 7, 9, 12, 13): (24, "0d58f5d337bee87d5f5ef68ca3ed4405629ba452f46de95893e33b1048e300bd"),
    (1, 5, 9, 10, 12, 14): (3, "7eab79fd3be005950b716652839ff21358b0c58d0d2b0c6ccb5f664306b0388b"),
    (1, 5, 9, 11, 12, 14): (24, "19cf0b2d81a624614d5f3aac83148cf4a41932bb50a262e3e6d1dcfc20220e03"),
    (1, 5, 9, 12, 13, 14): (32, "30461594f399063d6e9fa3d0bdd7310a0690f62649cd4cac532538f057e51541"),
    (1, 6, 7, 9, 10, 12): (9, "04dcfdd971b258ddf967cd906ca78bc064da50fbeebd4bc17c6f7b3a147bbcb6"),
    (1, 6, 8, 10, 11, 14): (25, "9de2204ba660c8e877b419d4200930ed61f59b705009b47db5e34bf94d0f1c35"),
    (1, 6, 8, 10, 13, 14): (16, "bbeae1ac7b9fe7595d0f4896da3432209f15857486797feebdb604d544b20a43"),
    (1, 6, 9, 10, 12, 14): (1, "2ec261fe63eca5be6ef8f0bdaa3655ddac941941a7702ca7cbe20a6b2c278a39"),
    (1, 6, 10, 11, 12, 14): (45, "40acad0e4e1f136b4282eb94e168ffbfafab1a476c20535c7b1021884acf928f"),
    (1, 6, 10, 12, 13, 14): (5, "014ff898a9c78d8af3d9fd2713f2841961f636f4f96002e58cfff363107641da"),
    (1, 6, 11, 12, 13, 14): (5, "dda8cb8de1b7b6ce8ad467eed27ae8bf0cab3bdb67e0bddee2c567146729e01e"),
    (1, 7, 9, 10, 11, 12): (36, "43b02d343a3b78d593aabcd4aa995f903b287934041db32ec660f1f1f7ba9b08"),
    (1, 7, 9, 10, 12, 14): (5, "3912f642794a89a22fad059c8f5288d4c3fe84702a7556a3c18e4dd5e689c8f5"),
    (1, 7, 9, 11, 12, 13): (24, "4b273ad70f328d566803fad954c92846485b3a5a5cb241a8877919e26ec0b055"),
    (1, 8, 10, 11, 12, 14): (30, "6b3d9166f708a4cd9e775e6ab3fa2968a2ede5bb38e99bc353653ebf81acbfcf"),
    (1, 8, 10, 12, 13, 14): (4, "a825cd38a28d6ab46aebe86171ecf7fbb6ec24299263ce4f8c10ec3f507c2e16"),
    (1, 9, 10, 11, 12, 14): (196, "03b72f6ac0c2e46480e9c1ddb3b759391dd42e45d5b7ce36a5fdfbf884e0a42e"),
    (1, 9, 11, 12, 13, 14): (72, "953110447a33ed0e4dcf20e0a2009b88d10c5287d181e3454bd98ee2297076a1"),
    (2, 3, 10, 11, 12, 14): (6, "dbe514d9dde6d7054d741d8fb56605c1225b67a71721d735b29663f929d41cf1"),
    (2, 5, 7, 9, 12, 13): (24, "0d58f5d337bee87d5f5ef68ca3ed4405629ba452f46de95893e33b1048e300bd"),
    (2, 5, 8, 9, 12, 14): (2, "d7a7bf70e995a1607a13f49a718be8bdfa0c9b54fea64f09046c8611c51c5e0b"),
    (2, 5, 9, 11, 12, 14): (39, "1ee65a7b4425a07375965833a9e534cb89cd8ab786158adb737fec84ce4c706e"),
    (2, 5, 9, 12, 13, 14): (40, "202cfdfc70e5a8977c77189fca8e2b523d76aa63d9a9eca2bc63f28b4403fd81"),
    (2, 6, 8, 10, 12, 14): (40, "fe3e563c42cfaa1aabdfd3b43fa05959c95ceb5098a894de8037712f7ad60ec8"),
    (2, 6, 10, 11, 12, 14): (4, "44029420c4ae24619960bf50162f4b3a585d1d2e5d0b56a2d5577c873108e55f"),
    (2, 7, 9, 11, 12, 13): (24, "4b273ad70f328d566803fad954c92846485b3a5a5cb241a8877919e26ec0b055"),
    (2, 8, 10, 11, 12, 14): (7, "5e2098e92b3bec251fac9b8f0e51f423de4f6ccbf27875c3065ed641da76f66d"),
    (2, 9, 10, 11, 12, 14): (163, "8f5a66a1aa936db886d957bf2528c38d79d5a231a7caa05f7727294e92d92ec0"),
    (2, 9, 10, 12, 13, 14): (6, "d539846541f1f002e19e8d459347c206598a5d33508cffe4bcc9d962223a9784"),
    (2, 9, 11, 12, 13, 14): (93, "ee262756fb06124232b953f14322ce77790c19bab52e85eb01d6b7f61c1703a7"),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    require(b"\r" not in payload, ("bare CR", path))
    return hashlib.sha256(payload).hexdigest()


def load(name, path=SOURCE_3174):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def valuation(number, prime):
    value = 0
    while number % prime == 0:
        number //= prime
        value += 1
    return value


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def boolean_join_controls():
    divisors = (1, 2, 3, 6)
    bits = lambda q: (q % 2 == 0, q % 3 == 0)
    for left in divisors:
        require(lcm(left, left) == left, (left, "idempotence"))
        for right in divisors:
            require(
                bits(lcm(left, right))
                == tuple(a or b for a, b in zip(bits(left), bits(right))),
                (left, right, "Boolean join"),
            )
    units = tuple(
        left
        for left in divisors
        if any(lcm(left, right) == 1 for right in divisors)
    )
    require(units == (1,), units)
    atom_swap = {1: 1, 2: 3, 3: 2, 6: 6}
    require(
        all(
            atom_swap[lcm(left, right)]
            == lcm(atom_swap[left], atom_swap[right])
            for left in divisors
            for right in divisors
        ),
        "external atom-swap hostile",
    )
    return divisors, units, tuple(atom_swap[q] for q in divisors)


def screen_worker(task):
    return load(f"thm3179_screen_{task[1]}").screen_worker(task)


def terminal_worker(task):
    index, body, bank = task
    wrapper = load(f"thm3179_terminal_wrapper_{index}")
    entry = wrapper.load(f"thm3179_terminal_entry_{index}")
    source = entry.load(f"thm3179_terminal_source_{index}")
    driver = source.load(f"thm3179_terminal_driver_{index}")
    predecessor = driver.load(f"thm3179_terminal_predecessor_{index}")
    bridge = predecessor.load(f"thm3179_terminal_bridge_{index}")
    base = bridge.load(f"thm3179_terminal_base_{index}")
    return index, base.thm.terminal_probe((LEVEL, body, bank))


def modulus_data(task, row=None):
    level, body, ruler, high, wall = task
    require(level == LEVEL, level)
    require(ruler == 14 * lcm(*body), (body, ruler, "body ruler"))
    require(ruler % 37 != 0, (body, ruler, "unexpected 37-factor"))
    modulus_gcd = gcd(level, ruler)
    expected_gcd = 6 if lcm(*body) % 3 == 0 else 2
    require(modulus_gcd == expected_gcd, (body, modulus_gcd, expected_gcd))
    first_d = ruler // modulus_gcd
    if row is not None:
        require(row[:6] == (level, body, ruler, high, first_d, wall), (body, "header"))
        for ds in row[13]:
            divisor = lcm(*ds)
            require(divisor % first_d == 0, (body, ds, first_d, divisor))
            require(ruler % divisor == 0, (body, ds, divisor, ruler))
            quotient = divisor // first_d
            require(modulus_gcd % quotient == 0, (body, ds, quotient, modulus_gcd))
    return modulus_gcd, first_d


def divisor_square(rows):
    quotient_states = Counter()
    coupling_states = Counter()
    for row in rows:
        body, ruler, first_d = row[1], row[2], row[4]
        modulus_gcd = gcd(LEVEL, ruler)
        require(modulus_gcd == 6, (body, "residual outside gcd-six stratum"))
        top_two = valuation(ruler, 2)
        top_three = valuation(ruler, 3)
        for ds in row[13]:
            divisor = lcm(*ds)
            quotient = divisor // first_d
            signatures = tuple(
                (valuation(d, 2) == top_two, valuation(d, 3) == top_three)
                for d in ds
            )
            separate = any(a for a, _b in signatures) and any(b for _a, b in signatures)
            coupled = any(a and b for a, b in signatures)
            require((quotient == 6) == separate, (body, ds, quotient, signatures))
            if quotient == 3:
                require(body == (2, 6, 8, 10, 12, 14), (body, ds, "q=3 body"))
                require(divisor == ruler // 2, (body, ds, divisor, ruler))
            quotient_states[(modulus_gcd, quotient)] += 1
            coupling_states[(quotient, coupled)] += 1
    divisor_packet = tuple(sorted(quotient_states.items()))
    coupling_packet = tuple(sorted(coupling_states.items()))
    require(divisor_packet == EXPECTED_DIVISOR_STATES, divisor_packet)
    require(coupling_packet == EXPECTED_COUPLING_STATES, coupling_packet)
    return divisor_packet, coupling_packet


def terminal_semantic(row):
    # row[6] is a chosen duplicate-gap witness and is deliberately omitted.
    return (*row[:6], *row[7:])


def closure_row(row):
    return (
        row[1],
        row[4],
        (row[5].numerator, row[5].denominator),
        row[7],
        row[8],
        row[9],
        row[10],
        row[11],
        row[12],
        row[14],
        row[16],
        row[17],
    )


def metadata(base, rows, next_rows):
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    next_census = (
        len(next_rows),
        sum(row[3] for row in next_rows),
        sum(not row[3] for row in next_rows),
    )
    require(census == EXPECTED_CENSUS, census)
    require(next_census == NEXT_CENSUS, next_census)
    require(hashlib.sha256(repr(rows).encode()).hexdigest() == EXPECTED_ROW_ORDER_SHA256, "row order")
    require(hashlib.sha256(repr(next_rows).encode()).hexdigest() == NEXT_ROW_ORDER_SHA256, "next row order")
    strata = Counter()
    exceptions = []
    for index, (body, ruler, high, wall) in enumerate(rows):
        modulus_gcd, first_d = modulus_data((LEVEL, body, ruler, high, wall))
        strata[(wall, modulus_gcd)] += 1
        if modulus_gcd == 2:
            exceptions.append((index, body, ruler, high, first_d))
    strata_packet = tuple(sorted(strata.items()))
    require(strata_packet == EXPECTED_STRATA, strata_packet)
    require(tuple(exceptions) == (EXPECTED_EXCEPTION,), exceptions)
    packet = (
        "lrc14-k3-z222-composite-terminal-cap221-v1",
        SOURCE_3174_SHA256,
        OUTPUT_3174_SHA256,
        SEMANTIC_3174_SHA256,
        base.thm.ATLAS_SHA256,
        (
            EXPECTED_ROW_ORDER_SHA256,
            EXPECTED_CENSUS,
            EXPECTED_STRATA,
            EXPECTED_EXCEPTION,
            EXPECTED_SCREEN,
            EXPECTED_FARKAS,
            EXPECTED_SCREEN_SHA256,
            EXPECTED_RESIDUAL_SHA256,
            EXPECTED_RESIDUAL_INDEX_SHA256,
            tuple(RESIDUAL.items()),
            EXPECTED_DIVISOR_STATES,
            EXPECTED_COUPLING_STATES,
        ),
        (
            EXPECTED_TERMINAL,
            EXPECTED_TERMINAL_SEMANTIC_SHA256,
            EXPECTED_CLOSURE_SHA256,
            EXPECTED_CASE_VECTOR_SHA256,
        ),
        LEDGER_BEFORE,
        LEDGER_AFTER,
        NEXT_CAP,
        (NEXT_LEVEL, NEXT_CENSUS, NEXT_ROW_ORDER_SHA256),
    )
    return census, next_census, strata_packet, hashlib.sha256(repr(packet).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=12)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--metadata-only", action="store_true")
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    require(sha(SOURCE_3174) == SOURCE_3174_SHA256, "THM-3174 source changed")
    require(sha(OUTPUT_3174) == OUTPUT_3174_SHA256, "THM-3174 output changed")
    require(
        f"semantic_sha256={SEMANTIC_3174_SHA256}" in OUTPUT_3174.read_text(encoding="utf-8"),
        "THM-3174 semantic changed",
    )
    join_controls = boolean_join_controls()
    require(
        hashlib.sha256(repr(tuple(RESIDUAL.items())).encode()).hexdigest()
        == EXPECTED_RESIDUAL_INDEX_SHA256,
        "residual index changed",
    )

    wrapper = load("thm3179_wrapper")
    entry = wrapper.load("thm3179_entry")
    source = entry.load("thm3179_source")
    driver = source.load("thm3179_driver")
    predecessor = driver.load("thm3179_predecessor")
    bridge = predecessor.load("thm3179_bridge")
    base = bridge.load("thm3179_base")
    rows = entry.atlas_rows(base, LEVEL)
    next_rows = entry.atlas_rows(base, NEXT_LEVEL)
    census, next_census, strata, semantic = metadata(base, rows, next_rows)
    if args.metadata_only:
        print("SEMANTIC_SHA256", semantic)
        return
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, semantic)

    tasks = tuple((LEVEL, body, ruler, high, wall) for body, ruler, high, wall in rows)
    if args.processes == 1:
        pairs = tuple(screen_worker(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            pairs = tuple(pool.map(screen_worker, tasks))
    by_task = {task: row for task, row in pairs}
    require(len(by_task) == len(tasks), "lost screen task")
    screened = tuple(by_task[task] for task in tasks)
    for task, row in zip(tasks, screened):
        modulus_data(task, row)
    totals = tuple(sum(row[index] for row in screened) for index in (9, 10, 11, 12))
    require(totals == EXPECTED_SCREEN, totals)
    farkas = (sum(row[19] for row in screened), sum(row[20] for row in screened))
    require(farkas == EXPECTED_FARKAS, farkas)
    require(all(row[16] == row[11] for row in screened), "status control")
    canonical_screen = tuple(row[:19] for row in screened)
    screen_sha = hashlib.sha256(repr(canonical_screen).encode()).hexdigest()
    require(screen_sha == EXPECTED_SCREEN_SHA256, screen_sha)
    residual_rows = tuple(row for row in screened if row[12])
    residual = tuple((row[1], row[13]) for row in residual_rows)
    residual_sha = hashlib.sha256(repr(residual).encode()).hexdigest()
    require(residual_sha == EXPECTED_RESIDUAL_SHA256, residual_sha)
    require(tuple(row[1] for row in residual_rows) == tuple(RESIDUAL), "residual bodies")
    for row in residual_rows:
        count, mask_sha = RESIDUAL[row[1]]
        require(row[12] == count, (row[1], row[12]))
        require(hashlib.sha256(repr(row[13]).encode()).hexdigest() == mask_sha, (row[1], "mask sha"))
    divisor_states, coupling_states = divisor_square(residual_rows)

    terminal_tasks = tuple((index, row[1], row[13]) for index, row in enumerate(residual_rows))
    if args.processes == 1:
        terminal_pairs = tuple(terminal_worker(task) for task in terminal_tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(terminal_tasks))) as pool:
            terminal_pairs = tuple(pool.map(terminal_worker, terminal_tasks))
    terminal_by_index = {index: row for index, row in terminal_pairs}
    terminals = tuple(terminal_by_index[index] for index in range(len(terminal_tasks)))
    terminal_summary = (
        len(terminals),
        sum(row[4] for row in terminals),
        sum(row[7] for row in terminals),
        sum(row[8] for row in terminals),
        sum(row[9] for row in terminals),
        sum(row[10] for row in terminals),
        sum(row[11] for row in terminals),
        sum(row[12] for row in terminals),
        sum(row[14] for row in terminals),
        min(row[16] for row in terminals),
        all(row[19] for row in terminals),
        all(row[20] for row in terminals),
    )
    require(terminal_summary == EXPECTED_TERMINAL, terminal_summary)
    terminal_semantics = tuple(terminal_semantic(row) for row in terminals)
    terminal_semantic_sha = hashlib.sha256(repr(terminal_semantics).encode()).hexdigest()
    require(terminal_semantic_sha == EXPECTED_TERMINAL_SEMANTIC_SHA256, terminal_semantic_sha)
    closure_packet = (
        "z222-terminal-closure-solver-witness-free-v1",
        tuple(closure_row(row) for row in terminals),
    )
    closure_sha = hashlib.sha256(repr(closure_packet).encode()).hexdigest()
    require(closure_sha == EXPECTED_CLOSURE_SHA256, closure_sha)
    case_vector_sha = hashlib.sha256(
        repr(tuple((row[1], row[17]) for row in terminals)).encode()
    ).hexdigest()
    require(case_vector_sha == EXPECTED_CASE_VECTOR_SHA256, case_vector_sha)
    require(LEDGER_BEFORE - LAYER_ROWS == LEDGER_AFTER, "ledger arithmetic")

    lines = [
        "LRC14 projected k3 z222 composite terminal descent and cap221",
        f"dependency=THM3174_source:{SOURCE_3174_SHA256};output:{OUTPUT_3174_SHA256};semantic:{SEMANTIC_3174_SHA256}",
        f"atlas=sha256:{base.thm.ATLAS_SHA256};rows:6060",
        f"layer=z1:{LEVEL};rows:{census[0]};wall:{census[1]};order:{census[2]};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256};states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};direct_farkas:{farkas[0]};legacy_farkas:{farkas[1]};screen_record_sha256:{screen_sha};residual_rows:{len(residual_rows)};residual_sha256:{residual_sha};residual_index_sha256:{EXPECTED_RESIDUAL_INDEX_SHA256}",
        f"divisor_square=ambient:{join_controls[0]};internal_units:{join_controls[1]};external_atom_swap:{join_controls[2]};strata:{strata};exception:{EXPECTED_EXCEPTION};residual_quotients:{divisor_states};coupling:{coupling_states};q3_body:(2,6,8,10,12,14)",
        f"terminal=z1:{LEVEL};rows:{terminal_summary[0]};masks:{terminal_summary[1]};positive_two_high_gap:{sum(row[19] for row in terminals)};closed:{sum(row[20] for row in terminals)};zero_high_hostiles:{terminal_summary[2]};one_high_cases:{terminal_summary[3]};low_label_sets:{terminal_summary[4]};coarse:{terminal_summary[5]};exact:{terminal_summary[6]};maxgap:{terminal_summary[7]};failures:{terminal_summary[8]};minimum_slack:{terminal_summary[9]};terminal_semantic_sha256:{terminal_semantic_sha};closure_sha256:{closure_sha};case_vector_sha256:{case_vector_sha}",
    ]
    for row in terminals:
        lines.append(
            f"TERMINAL;E={row[1]};L={row[2]};high={row[3]};masks={row[4]};two_high_gap={ftext(row[5])};zero_high_hostiles={row[7]};one_high_cases={row[8]};low_label_sets={row[9]};coarse={row[10]};exact={row[11]};maxgap={row[12]};failed={row[14]};minimum_slack={row[16]};case_certificate_sha256={row[17]}"
        )
    lines.extend(
        [
            "direction_screen=the_ray_quotient_and_exact_Farkas_status_checks_close_a_superset_of_every_actual_projected_assignment;all_order_rows_and_the_unique_gcd2_wall_row_have_zero_residual",
            "direction_divisor_square=on_every_residual_gcd(222,L)=6_and_first_d=L/6;D=lcm(ds)_satisfies_first_d|D|L;the_missing_top_2_and_3_valuation_bits_combine_idempotently_by_lcm_OR;this_is_not_a_V4_torsor_or_order_sensitive_restoration",
            "direction_terminal=each_residual_is_a_wall_row;the_strictly_positive_duplicate_permitting_two_high_gap_and_wall_at_least_one_high_gate_force_exactly_one_high;the_one_high_bank_uses_high_ray_suprema_and_enlarges_the_actual_assignment_set",
            "direction_carrier=every_one_high_case_is_closed_by_complete_cell_projected_support_cardinality;no_max_gap_fallback_is_used",
            "evidence_boundary=all_raw_Farkas_certificates_are_verified_exactly_but_the_screen_digest_uses_only_row[:19];the_terminal_closure_digest_omits_the_chosen_duplicate_gap_maximizer_witness;all_37_residual_mask_banks_are_bound_individually",
            f"promotion_consequence=ledger {LEDGER_BEFORE}-{LAYER_ROWS}={LEDGER_AFTER};projected_k3_cap:z1<={NEXT_CAP}",
            f"next_layer=z1:{NEXT_LEVEL};rows:{next_census[0]};wall:{next_census[1]};order:{next_census[2]};row_order_sha256:{NEXT_ROW_ORDER_SHA256};status:occupied",
            "scope=projected_k3_necessary_atlas_only;the_divisor_square_is_a_noninvertible_lcm_semilattice_not_a_group_action;no_physical_cover_classification_outside_the_projection;no_k<=1_or_final_rung_or_LRC14_claim",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        ]
    )
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    mp.freeze_support()
    main()
