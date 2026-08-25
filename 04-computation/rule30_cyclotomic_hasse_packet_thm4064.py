#!/usr/bin/env python3
"""Rule-30 C60/cyclotomic/divisor-depth synthesis probe for THM-4064.

Reconstructs the finite THM-4047 fixed-tail bank, then analyzes every tail
under the exact C_lcm(60,p) cyclotomic fibre and THM-4056/4059 labels.
"""

from __future__ import annotations

import ast
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from math import gcd
from pathlib import Path


MAX_R = 100_000
PARENT_SHA256 = (
    ("THM4047_script", "d454cc5b40315c02ebf486f29e227ebc9a79c78ede18fc449a7dd32f8dc21148"),
    ("THM4047_output", "bae500127999c350ebeff77c3145fbb4abf7b5b7292d2882eeecc6492fba75e3"),
    ("THM4055_script", "c6de48e956c425b21d20f4dc364e10d41e3fd7d7401c486c0aee08e416c1f535"),
    ("THM4055_output", "e6ed3dfa72f85e33d4b702f1cfcaa48965624f44e8a2f458601b5cd44d58728f"),
    ("THM4044_script", "cc49cd7024fdeaab6c0d668c9ca497ee113ea6c45b3c9735d836322781629898"),
    ("THM4044_output", "4e73f5d3fed3bec3966d914f4981777f1a3b78892321bd443ed83f2799185fc0"),
    ("THM4010_script", "f28dfba947010cb1b82891b1c1b6981fbc2f0555865840eddb75693d29c10888"),
    ("THM4010_output", "2a97b20bea13267d3a1d324975f921e9399e59861f2206fd9c3a72ce1108cc7e"),
    ("THM4056_script", "a15fb134839f84b5e1a4f07131f893af790adcc2c6b8189780c851a9e90d77e9"),
    ("THM4056_output", "bce0b8ccc22db995a0d96560aa66c7707b581b45e419a02589920a4249b28947"),
    ("THM4059_script", "4279d31e869764ea80febfdfaad90689cb5e679847929c04cb7ed5768efd56f4"),
    ("THM4059_output", "313f79dc4f21b69a62847796211a55230264d5d4ccf7abb8796c6ac702fdd3d8"),
)
CHECKS = 0


def require(condition, label):
    global CHECKS
    if not condition:
        raise RuntimeError(label)
    CHECKS += 1


@dataclass(frozen=True)
class Steady:
    start: int
    word: tuple[int, ...]
    monodromy: tuple[int, int]

    def value(self, time):
        require(time >= self.start, ("pre-tail query", time, self.start))
        return self.word[time % len(self.word)]


def one_step(left, right, state):
    return left ^ right ^ state ^ (right & state)


def absolute_word(start, values):
    raw_period = len(values)
    period = next(
        candidate
        for candidate in range(1, raw_period + 1)
        if raw_period % candidate == 0
        and all(
            values[index] == values[index % candidate]
            for index in range(raw_period)
        )
    )
    word = [None] * period
    for offset, value in enumerate(values):
        residue = (start + offset) % period
        require(
            word[residue] is None or word[residue] == value,
            "absolute phase conflict",
        )
        word[residue] = value
    require(all(value is not None for value in word), "absolute phase hole")
    return tuple(int(value) for value in word)


def reconstruct_bank():
    steady = [Steady(0, (1,), (1, 0)), Steady(1, (1,), (0, 1))]
    mask = (1 << (MAX_R + 1)) - 1
    front = 1
    simulation_time = 0

    def advance_to(target):
        nonlocal front, simulation_time
        require(target >= simulation_time, "nonmonotone physical query")
        while simulation_time < target:
            shifted = front << 1
            front = (
                (shifted << 1) ^ shifted ^ front ^ (shifted & front)
            ) & mask
            simulation_time += 1

    for r in range(2, MAX_R + 1):
        left = steady[r - 2]
        right = steady[r - 1]
        common_start = max(left.start, right.start)
        common_period = max(len(left.word), len(right.word))
        require(
            common_period % len(left.word) == 0
            and common_period % len(right.word) == 0,
            ("nondyadic periods", r),
        )
        image_zero, image_one = 0, 1
        for time in range(common_start, common_start + common_period):
            image_zero = one_step(left.value(time), right.value(time), image_zero)
            image_one = one_step(left.value(time), right.value(time), image_one)
        multiplier = image_zero ^ image_one
        constant = image_zero
        if multiplier == 0:
            start = common_start + common_period
            length = common_period
            initial = constant
        else:
            start = common_start
            length = common_period if constant == 0 else 2 * common_period
            advance_to(start)
            initial = (front >> r) & 1
        values = []
        state = initial
        for time in range(start, start + length):
            values.append(state)
            state = one_step(left.value(time), right.value(time), state)
        require(state == initial, ("period wrap", r))
        word = absolute_word(start, tuple(values))
        item = Steady(start, word, (multiplier, constant))
        for time in range(
            start, start + 2 * max(common_period, len(word))
        ):
            require(
                item.value(time + 1)
                == one_step(
                    left.value(time), right.value(time), item.value(time)
                ),
                ("wrapped recurrence", r, time),
            )
        steady.append(item)
    return steady


def stern_depth(a, d):
    require(d > 1 and 0 < a < d and gcd(a, d) == 1, ("depth domain", a, d))
    total = 0
    numerator, denominator = a, d
    while numerator:
        digit, remainder = divmod(denominator, numerator)
        total += digit
        denominator, numerator = numerator, remainder
    return total - 1


@lru_cache(maxsize=None)
def packet(time, modulus):
    residue = time % modulus
    common = gcd(residue, modulus)
    d = modulus // common
    if d == 1:
        # THM-4059 assigns no Stern depth to 0/1. None is a sentinel, not a
        # numerical depth or sign.
        return (1, 0, None, None, None)
    a = residue // common
    depth = stern_depth(a, d)
    inverse = pow(a, -1, d)
    carry = (a * inverse - 1) // d
    flank_mod_two = (
        carry & 1,
        (a - carry) & 1,
        inverse & 1,
        (d - inverse) & 1,
    )
    return (d, a, depth, -1 if depth & 1 else 1, flank_mod_two)


def all_factor_flags(word, key_functions):
    period = len(word)
    modulus = 60 * max(period // gcd(60, period), 1)
    observed = {name: {} for name in key_functions}
    valid = {name: True for name in key_functions}
    for time in range(modulus):
        response = word[time % period]
        for name, key_function in key_functions.items():
            if not valid[name]:
                continue
            key = key_function(time, modulus)
            if key in observed[name] and observed[name][key] != response:
                valid[name] = False
            observed[name][key] = response
    return valid


def cyclotomic_zero(coefficients, order):
    """Test a polynomial at a primitive root of 2-power order."""
    if order == 1:
        return sum(coefficients) == 0
    if order == 2:
        return (
            sum(value * (-1) ** index for index, value in enumerate(coefficients))
            == 0
        )
    half = order // 2
    remainder = [0] * half
    for index, value in enumerate(coefficients):
        remainder[index % half] += value * (-1) ** (index // half)
    return all(value == 0 for value in remainder)


def sign_fourier_support(word):
    """Fourier support of s=(-1)^word, not of the raw bit word."""
    period = len(word)
    coefficients = tuple(1 - 2 * value for value in word)
    support = []
    for frequency in range(period):
        order = period // gcd(period, frequency)
        if not cyclotomic_zero(coefficients, order):
            support.append(frequency)
    return tuple(support)


steady = reconstruct_bank()
require(len(steady) == 100001, "bank size")
require(
    Counter(len(item.word) for item in steady)
    == Counter({1: 16, 2: 10, 4: 56, 8: 668, 16: 87118, 32: 12133}),
    "parent period histogram",
)

key_functions = {
    "C60": lambda time, modulus: (time % 60,),
    "C60_denominator": lambda time, modulus: (
        time % 60,
        packet(time, modulus)[0],
    ),
    "C60_depth_sign": lambda time, modulus: (
        time % 60,
        packet(time, modulus)[3],
    ),
    "C60_denominator_depth_sign": lambda time, modulus: (
        time % 60,
        packet(time, modulus)[0],
        packet(time, modulus)[3],
    ),
    "C60_denominator_full_depth": lambda time, modulus: (
        time % 60,
        packet(time, modulus)[0],
        packet(time, modulus)[2],
    ),
    "C60_denominator_flank_mod2": lambda time, modulus: (
        time % 60,
        packet(time, modulus)[0],
        packet(time, modulus)[4],
    ),
    "C60_full_exact_denominator_packet": lambda time, modulus: (
        time % 60,
        packet(time, modulus)[0],
        packet(time, modulus)[1],
    ),
}

word_counts = Counter(item.word for item in steady)
word_factor_flags = {
    word: all_factor_flags(word, key_functions) for word in word_counts
}
factor_counts = {
    name: sum(
        count
        for word, count in word_counts.items()
        if word_factor_flags[word][name]
    )
    for name in key_functions
}
factor_counts_by_period = {
    name: tuple(
        (
            period,
            sum(
                count
                for word, count in word_counts.items()
                if len(word) == period and word_factor_flags[word][name]
            ),
        )
        for period in (1, 2, 4, 8, 16, 32)
    )
    for name in key_functions
}

require(factor_counts["C60"] == 82, "THM4055 factor count")
require(
    factor_counts["C60_full_exact_denominator_packet"] == 100001,
    "full compiler inversion",
)

sign_sector_histogram = Counter()
sign_support_histogram = Counter()
for word, count in word_counts.items():
    period = len(word)
    support = sign_fourier_support(word)
    m = period // gcd(60, period)
    sectors = (
        tuple(sorted({(15 * frequency) % m for frequency in support}))
        if m > 1
        else (0,)
    )
    sign_sector_histogram[(period, sectors)] += count
    order_shells = tuple(
        sorted({period // gcd(period, frequency) for frequency in support})
    )
    sign_support_histogram[(period, order_shells)] += count

ell29 = steady[29]
require(
    ell29.start == 90 and ell29.word == (1, 1, 0, 1, 0, 0, 1, 0),
    "ell29 parent word",
)
hostile_times = (90, 150)
hostile_packets = tuple(packet(time, 120) for time in hostile_times)
require(hostile_packets[0][0] == hostile_packets[1][0] == 4, "hostile denominator")
require(
    hostile_packets[0][2:] == hostile_packets[1][2:],
    "hostile depth/flank collision",
)
require(
    ell29.value(90) == 0 and ell29.value(150) == 1,
    "hostile response split",
)
require(
    tuple(1 - 2 * ell29.value(time) for time in hostile_times) == (1, -1),
    "hostile sign split",
)
ell29_sign_support = sign_fourier_support(ell29.word)
require(ell29_sign_support == (1, 3, 5, 7), "ell29 sign support")

# Characteristic zero is essential for the primitive-frequency lemma. In
# F_13bar, zeta_4=5 and the least-period word (5,1,0,0) has coefficient
# polynomial 5+X, which vanishes at zeta_4^-1=8.
require(pow(5, 2, 13) == 12 and pow(5, -1, 13) == 8, "finite-char root")
require((5 + pow(5, -1, 13)) % 13 == 0, "finite-char primitive zero")

factor_payload = repr((factor_counts, factor_counts_by_period))
sign_sector_payload = repr(
    (
        sorted(sign_sector_histogram.items()),
        sorted(sign_support_histogram.items()),
    )
)
print("== Rule30 C60 cyclotomic/divisor-depth synthesis ==")
print("status=FINITE-EXACT fixed-bank plus exact universal decompositions")
print(f"parent_sha256={PARENT_SHA256}")
print("typed_map=C_N --t|->zeta_N^t--> mu_N; exact_denominator=root_order")
print("full_fibre_product=mu_N~=mu_60 x_mu_g mu_p; preserves both phases and order")
print("C60_projection=z|->z^m; kernel=mu_m; loses q_torsor/nonzero isotypic sectors")
print("Hasse_parent=THM4044;mu_m-invariants=A_60,k via Q=P^m")
print("sign_Fourier_support=k=(N/p)h;mu_m_weight=15h_mod_m_for_p>=4")
print("d1_depth_convention=None_sentinel_not_a_Stern_depth")
print(f"factor_counts={factor_counts}")
print(f"factor_counts_by_period={factor_counts_by_period}")
print(f"distinct_tail_words={len(word_counts)}")
print(f"factor_census_sha256={sha256(factor_payload.encode()).hexdigest()}")
print(f"sign_sector_histogram={tuple(sorted(sign_sector_histogram.items()))}")
print(
    "sign_fourier_order_shell_histogram="
    + repr(tuple(sorted(sign_support_histogram.items())))
)
print(
    "sign_sector_census_sha256="
    + sha256(sign_sector_payload.encode()).hexdigest()
)
print("minimal_hostile=ell_29 at times 90,150 in C_120")
print(f"hostile_clock_phases={tuple(time % 60 for time in hostile_times)}")
print(f"hostile_C120_residues={tuple(time % 120 for time in hostile_times)}")
print(f"hostile_packets=(d,a,D,epsilon,flank_mod2)={hostile_packets}")
print(f"hostile_values={tuple(ell29.value(time) for time in hostile_times)}")
print(
    "ell29_sign_fourier_support_C8="
    + repr(ell29_sign_support)
    + ";lifted_C120_exponents=(15,45,75,105)"
)
print("loss=C60_plus_denominator/depth/sign/flank all identify hostile pair")
print("sidecars=point_torsor_or_labelled_numerator;dually_nontrivial_mu_2_sector")
print("Smith_firewall=no_THM4010_integral_Smith_transfer_by_exponentiation_alone")
print("finite_char_hostile=F13bar;p4;word_(5,1,0,0);primitive_mode_zero")
print("SCOPE fixed left-front tails only; moving center and Rule30 prizes OPEN")

source_bytes = Path(__file__).read_bytes()
require(b"\r\n" not in source_bytes, "LF source")
tree = ast.parse(source_bytes.decode())
require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "Assert node")
print(f"script_sha256={sha256(source_bytes).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
