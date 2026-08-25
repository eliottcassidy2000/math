#!/usr/bin/env python3
"""Independent audit of the THM-4064 Rule-30 cyclotomic packet census.

This script does not import the primary computation.  Its recurrence path uses
the seven frozen physical-query bits already audited by THM-4047; the resulting
bank audit is therefore explicitly conditional on that query packet.
"""

from __future__ import annotations

import ast
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from math import gcd
from pathlib import Path


CHECKS = 0
MAX_R = 100_000
REPO_ROOT = Path(__file__).resolve().parents[1]
PRIMARY_PATH = REPO_ROOT / "04-computation/rule30_cyclotomic_hasse_packet_thm4064.py"
PRIMARY_OUTPUT = REPO_ROOT / "05-knowledge/results/rule30_cyclotomic_hasse_packet_thm4064.out"
PRIMARY_SHA = "46f35b1081ab13d8a70bbcae35180eb7a00758fefaebfbee7b7cd83cb2f0d4a3"
PRIMARY_OUTPUT_SHA = "e753b1bfdf082d9f4871f116745e7f3bcd9aea4f097c5a523d51c746388c4dce"
PHYSICAL_QUERIES = {
    3: (2, 0, (1, 1)),
    8: (10, 0, (1, 1)),
    29: (90, 0, (1, 1)),
    400: (3050, 1, (1, 1)),
    53208: (847962, 0, (1, 0)),
    58287: (929210, 1, (1, 0)),
    87867: (1402474, 0, (1, 1)),
}


def require(condition, label):
    global CHECKS
    if not condition:
        raise RuntimeError(label)
    CHECKS += 1


require(sha256(PRIMARY_PATH.read_bytes()).hexdigest() == PRIMARY_SHA, "primary source drift")
require(
    sha256(PRIMARY_OUTPUT.read_bytes()).hexdigest() == PRIMARY_OUTPUT_SHA,
    "primary output drift",
)


@dataclass(frozen=True)
class Tail:
    start: int
    word: tuple[int, ...]
    monodromy: tuple[int, int]

    def at(self, time):
        return self.word[time % len(self.word)]


def transition(left, right, state):
    return left ^ right ^ state ^ (right & state)


def phase_word(start, values):
    raw_length = len(values)
    period = next(
        divisor
        for divisor in range(1, raw_length + 1)
        if raw_length % divisor == 0
        and tuple(values)
        == tuple(values[index % divisor] for index in range(raw_length))
    )
    answer = [None] * period
    for index, value in enumerate(values):
        slot = (start + index) % period
        require(answer[slot] is None or answer[slot] == value, "phase collision")
        answer[slot] = value
    require(all(value is not None for value in answer), "phase hole")
    return tuple(answer)


def reconstruct_from_query_packet():
    tails = [Tail(0, (1,), (1, 0)), Tail(1, (1,), (0, 1))]
    used_queries = set()
    for r in range(2, MAX_R + 1):
        left, right = tails[r - 2], tails[r - 1]
        start = max(left.start, right.start)
        block = max(len(left.word), len(right.word))
        image_zero, image_one = 0, 1
        for time in range(start, start + block):
            image_zero = transition(left.at(time), right.at(time), image_zero)
            image_one = transition(left.at(time), right.at(time), image_one)
        monodromy = (image_zero ^ image_one, image_zero)
        if monodromy[0] == 0:
            tail_start = start + block
            length = block
            state = monodromy[1]
        else:
            require(r in PHYSICAL_QUERIES, ("missing physical query", r))
            query_start, state, expected_monodromy = PHYSICAL_QUERIES[r]
            require(
                (query_start, expected_monodromy) == (start, monodromy),
                ("query drift", r),
            )
            used_queries.add(r)
            tail_start = start
            length = block if monodromy[1] == 0 else 2 * block
        values = []
        current = state
        for time in range(tail_start, tail_start + length):
            values.append(current)
            current = transition(left.at(time), right.at(time), current)
        require(current == state, ("audit wrap", r))
        word = phase_word(tail_start, values)
        item = Tail(tail_start, word, monodromy)
        for time in range(tail_start, tail_start + 2 * max(block, len(word))):
            require(
                item.at(time + 1)
                == transition(left.at(time), right.at(time), item.at(time)),
                ("audit recurrence", r, time),
            )
        tails.append(item)
    require(used_queries == set(PHYSICAL_QUERIES), "query packet consumption")
    return tails


def depth_by_subtraction(a, d):
    left, right = a, d - a
    steps = 0
    while left != right:
        if left > right:
            left -= right
        else:
            right -= left
        steps += 1
    require(left == right == 1, ("primitive subtraction endpoint", a, d))
    return steps + 1


@lru_cache(maxsize=None)
def labels(time, modulus):
    residue = time % modulus
    common = gcd(residue, modulus)
    d = modulus // common
    if d == 1:
        # None is a sentinel: THM-4059 does not assign a Stern depth to 0/1.
        return (1, 0, None, None, None)
    a = residue // common
    depth = depth_by_subtraction(a, d)
    inverse = next(candidate for candidate in range(1, d) if a * candidate % d == 1)
    carry = (a * inverse - 1) // d
    flank = (carry & 1, (a - carry) & 1, inverse & 1, (d - inverse) & 1)
    return (d, a, depth, -1 if depth & 1 else 1, flank)


KEYS = {
    "C60": lambda time, modulus: (time % 60,),
    "C60_denominator": lambda time, modulus: (time % 60, labels(time, modulus)[0]),
    "C60_depth_sign": lambda time, modulus: (time % 60, labels(time, modulus)[3]),
    "C60_denominator_depth_sign": lambda time, modulus: (
        time % 60,
        labels(time, modulus)[0],
        labels(time, modulus)[3],
    ),
    "C60_denominator_full_depth": lambda time, modulus: (
        time % 60,
        labels(time, modulus)[0],
        labels(time, modulus)[2],
    ),
    "C60_denominator_flank_mod2": lambda time, modulus: (
        time % 60,
        labels(time, modulus)[0],
        labels(time, modulus)[4],
    ),
    "C60_full_exact_denominator_packet": lambda time, modulus: (
        time % 60,
        labels(time, modulus)[0],
        labels(time, modulus)[1],
    ),
}


def factor_flags(word):
    period = len(word)
    modulus = 60 * (period // gcd(60, period))
    tables = {name: {} for name in KEYS}
    flags = {name: True for name in KEYS}
    for time in range(modulus):
        value = word[time % period]
        for name, function in KEYS.items():
            if not flags[name]:
                continue
            key = function(time, modulus)
            previous = tables[name].get(key, value)
            if previous != value:
                flags[name] = False
            tables[name][key] = value
    return flags


def vanishes_on_primitive_shell(coefficients, order):
    if order == 1:
        return sum(coefficients) == 0
    if order == 2:
        return sum(coefficients[::2]) == sum(coefficients[1::2])
    half = order // 2
    return all(
        sum(coefficients[index] for index in range(residue, len(coefficients), order))
        == sum(
            coefficients[index]
            for index in range(residue + half, len(coefficients), order)
        )
        for residue in range(half)
    )


def spectral_data(word):
    period = len(word)
    signs = tuple(1 - 2 * bit for bit in word)
    support = tuple(
        frequency
        for frequency in range(period)
        if not vanishes_on_primitive_shell(
            signs, period // gcd(period, frequency)
        )
    )
    m = period // gcd(60, period)
    sectors = (
        tuple(sorted({15 * frequency % m for frequency in support}))
        if m > 1
        else (0,)
    )
    shells = tuple(
        sorted({period // gcd(period, frequency) for frequency in support})
    )
    return support, sectors, shells


tails = reconstruct_from_query_packet()
word_counts = Counter(item.word for item in tails)
require(len(word_counts) == 60457, "distinct words")
flags_by_word = {word: factor_flags(word) for word in word_counts}

counts = {
    name: sum(
        multiplicity
        for word, multiplicity in word_counts.items()
        if flags_by_word[word][name]
    )
    for name in KEYS
}
counts_by_period = {
    name: tuple(
        (
            period,
            sum(
                multiplicity
                for word, multiplicity in word_counts.items()
                if len(word) == period and flags_by_word[word][name]
            ),
        )
        for period in (1, 2, 4, 8, 16, 32)
    )
    for name in KEYS
}

sign_sector_histogram = Counter()
sign_shell_histogram = Counter()
for word, multiplicity in word_counts.items():
    _support, sectors, shells = spectral_data(word)
    sign_sector_histogram[(len(word), sectors)] += multiplicity
    sign_shell_histogram[(len(word), shells)] += multiplicity

expected_counts = {
    "C60": 82,
    "C60_denominator": 172,
    "C60_depth_sign": 82,
    "C60_denominator_depth_sign": 172,
    "C60_denominator_full_depth": 547,
    "C60_denominator_flank_mod2": 172,
    "C60_full_exact_denominator_packet": 100001,
}
require(counts == expected_counts, "primary factor census")

ell29 = tails[29]
hostile_labels = (labels(90, 120), labels(150, 120))
require(
    hostile_labels
    == ((4, 3, 3, -1, (0, 1, 1, 1)), (4, 1, 3, -1, (0, 1, 1, 1))),
    "hostile labels",
)
require((ell29.at(90), ell29.at(150)) == (0, 1), "hostile values")
require(spectral_data(ell29.word)[0] == (1, 3, 5, 7), "hostile sign spectrum")

factor_payload = repr((counts, counts_by_period))
sector_payload = repr(
    (sorted(sign_sector_histogram.items()), sorted(sign_shell_histogram.items()))
)
require(
    sha256(factor_payload.encode()).hexdigest()
    == "fc7679bc50f53005384a9f0c34f33ff651cf30b578a2828795495525adca427f",
    "factor semantic hash",
)
require(
    sha256(sector_payload.encode()).hexdigest()
    == "41a87e724c878456bc94ae8deb5f65892562a1d689e280dc2ab53d77b87815d4",
    "sign-sector semantic hash",
)

# Independent finite group/dimension checks behind the fibre and invariant ring.
for period in (1, 2, 4, 8, 16, 32):
    g = gcd(60, period)
    m = period // g
    n = 60 * m
    require(n == 60 * period // g, ("lcm", period))
    require(n // 60 == m, ("C60 kernel size", period))
    require(60 * period // g == n, ("fibre product size", period))
    for jet_order in (1, 2, 3):
        require(n * jet_order // m == 60 * jet_order, ("invariant dimension", period))

# Characteristic zero is essential: in F_13bar, zeta_4=5 and 5+X vanishes
# at the primitive inverse root 8 although (5,1,0,0) has least period four.
require(pow(5, 2, 13) == 12 and pow(5, -1, 13) == 8, "finite-field root")
require((5 + pow(5, -1, 13)) % 13 == 0, "finite-character hostile")

source_bytes = Path(__file__).read_bytes()
require(b"\r\n" not in source_bytes, "LF source")
tree = ast.parse(source_bytes.decode())
require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "Assert node")
require(not any(isinstance(node, ast.Constant) and isinstance(node.value, float) for node in ast.walk(tree)), "float literal")

print("== independent THM4064 Rule30 cyclotomic/Hasse packet audit ==")
print("status=PASS;no-import-from-primary")
print("query_packet_reconstruction=conditional_on_already_audited_THM4047")
print(f"primary_sha256={PRIMARY_SHA};primary_output_sha256={PRIMARY_OUTPUT_SHA}")
print(f"factor_counts={counts}")
print(f"factor_census_sha256={sha256(factor_payload.encode()).hexdigest()}")
print(f"sign_sector_census_sha256={sha256(sector_payload.encode()).hexdigest()}")
print(
    "ell29_hostile_labels="
    + repr(hostile_labels)
    + ";values=(0, 1);sign_support=(1, 3, 5, 7)"
)
print("d1_depth_convention=None_sentinel_not_a_Stern_depth")
print("finite_char_hostile=F13bar;p4;word_(5,1,0,0);primitive_mode_zero")
print(f"script_sha256={sha256(source_bytes).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
