#!/usr/bin/env python3
"""Exact causal-horizon audit for the AMM 12592 R=512,D0=0 rule-A death.

The script independently rebuilds the exact golden floor word and the
error-coordinate rule-A trajectory through the fatal input at row 107.  It
then proves the coefficient-capacity bound used in THM-3374:

  |[x^s](Delta-Delta')| <= 2^(s+1) binom(d,s)

for any two degree-d Lucas-box blocks.  No numerical solver or non-standard
package is used.
"""

from __future__ import annotations

from hashlib import sha256
from math import comb
from pathlib import Path


R = 512
DEATH_ROW = 107
EXPECTED_FATAL = (
    199181334599768561751691151979867147451295075845970943582846950031770442839710071820
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def trim(a: list[int]) -> list[int]:
    while a and a[-1] == 0:
        a.pop()
    return a


def sub(a: list[int], b: list[int]) -> list[int]:
    out = list(a) + [0] * max(0, len(b) - len(a))
    for i, value in enumerate(b):
        out[i] -= value
    return trim(out)


def qpow(n: int) -> list[int]:
    return [(-1) ** j * comb(n, j) for j in range(n + 1)]


def em(n: int) -> list[int]:
    return [-1] + [1] * n if n >= 0 else []


def fib_pair(n: int) -> tuple[int, int]:
    if n == 0:
        return 0, 1
    a, b = fib_pair(n >> 1)
    c = a * (2 * b - a)
    d = a * a + b * b
    return (d, c + d) if n & 1 else (c, d)


def compare_five_power_to_phi2(d: int, m: int) -> int:
    """Return sign(5**d-phi**(2m)) using Fibonacci/Lucas integers."""
    f, f_next = fib_pair(2 * m)
    lucas = 2 * f_next - f
    a = 2 * 5**d - lucas
    if a <= 0:
        return -1
    delta = a * a - 5 * f * f
    return (delta > 0) - (delta < 0)


def floor_gamma_star(m: int) -> int:
    d = 3 * m // 5
    while compare_five_power_to_phi2(d + 1, m) <= 0:
        d += 1
    while compare_five_power_to_phi2(d, m) > 0:
        d -= 1
    return d


def poly_to_bern(poly: list[int], d: int) -> list[int]:
    """Coordinates dual to x^(d-k)(1-x)^k, as in THM-3329."""
    require(len(poly) <= d + 1, "polynomial truncation exceeds degree")
    return [
        sum(poly[j] * comb(d - j, t) for j in range(len(poly)))
        for t in range(d + 1)
    ]


def bern_to_poly(cells: list[int], d: int) -> list[int]:
    out = [0] * (d + 1)
    for k, cell in enumerate(cells):
        if cell == 0:
            continue
        h = d - k
        for j in range(k + 1):
            out[h + j] += cell * (-1) ** j * comb(k, j)
    return trim(out)


def ebox(d: int, t: int) -> tuple[int, int]:
    return (
        -2 * comb(d - 1, t),
        2 * (comb(d - 1, t - 1) if t else 0),
    )


def replay_to_death(profile: list[int]) -> tuple[int, int, list[int]]:
    """Return fatal input, parity-fire count, and surviving row-top cells."""
    e = sub(qpow(R - 1), em(R - 1))
    parity_fires = 0
    top_cells: list[int] = []
    for i in range(DEATH_ROW + 1):
        d = profile[i]
        if i == DEATH_ROW:
            return e[0], parity_fires, top_cells
        corner = [
            comb(d - 1, t) - (comb(d - 1, t - 1) if t else 0)
            for t in range(d + 1)
        ]
        w = poly_to_bern(e[: d + 1], d)
        c: list[int] = []
        for t, value in enumerate(w):
            lo, hi = ebox(d, t)
            value = min(hi, max(lo, value))
            if (value - lo) % 2:
                parity_fires += 1
                delta_value = corner[t] + value
                value = value - 1 if delta_value > 0 else value + 1
            c.append(value)
        delta = [x + y for x, y in zip(corner, c)]
        require(
            all(
                abs(value) <= comb(d, t)
                and (value - comb(d, t)) % 2 == 0
                for t, value in enumerate(delta)
            ),
            f"row {i} is not Lucas-box admissible",
        )
        rem = sub(e, bern_to_poly(c, d))
        require(not rem or rem[0] == 0, f"unexpected earlier death at row {i}")
        e = rem[1:] if rem else []
        top_cells.append(delta[-1])
    raise RuntimeError("death row was not reached")


def coefficient_capacity(d: int, s: int) -> int:
    """Maximum triangle-inequality bound for a block-pair x^s coefficient."""
    direct = sum(
        2 * comb(d, h) * comb(d - h, s - h)
        for h in range(s + 1)
    )
    closed = 2 ** (s + 1) * comb(d, s)
    require(direct == closed, "Vandermonde coefficient capacity")
    return closed


def main() -> None:
    profile = [floor_gamma_star(R + i) for i in range(R)]
    fatal, parity_fires, tops = replay_to_death(profile)
    require(fatal == -EXPECTED_FATAL, "fatal constant")
    require(parity_fires == 0, "dyadic parity fire")
    require(len(tops) == DEATH_ROW, "surviving prefix length")
    require(all(value in (-1, 1) for value in tops), "top-cell parity gate")

    capacities: list[int] = []
    total = 0
    first_inconclusive = None
    for length in range(1, DEATH_ROW + 1):
        degree = profile[DEATH_ROW - length]
        total += coefficient_capacity(degree, length)
        capacities.append(total)
        if first_inconclusive is None and total >= EXPECTED_FATAL:
            first_inconclusive = length

    require(first_inconclusive == 58, "capacity threshold")
    require(capacities[56] < EXPECTED_FATAL <= capacities[57], "57/58 gate")
    five = capacities[4]
    fifty_seven = capacities[56]
    fifty_eight = capacities[57]

    source = Path(__file__).read_bytes().replace(b"\r\n", b"\n")
    print("AMM12592 R512 RULE-A CAUSAL REPAIR HORIZON")
    print("status=PROVED+VERIFIED-EXACT")
    print(
        f"carrier=R={R};D0=0;death_row={DEATH_ROW};"
        f"surviving_rows=0..{DEATH_ROW-1};parity_fires={parity_fires}"
    )
    print(f"fatal_input_constant={fatal};absolute_bits={EXPECTED_FATAL.bit_length()}")
    print(
        "coefficient_lemma=for_two_degree_d_Lucas_blocks;"
        "abs_coeff_x^s_difference<=2^(s+1)*C(d,s);"
        "proof=triangle_inequality_plus_Vandermonde"
    )
    print(
        f"degree_word_rows_97_107={profile[97:108]};"
        f"delta_word={''.join(str(profile[i+1]-profile[i]) for i in range(97,107))}"
    )
    print(f"five_row_capacity={five};bits={five.bit_length()}")
    print(f"length_57_capacity={fifty_seven};bits={fifty_seven.bit_length()}")
    print(f"length_58_capacity={fifty_eight};bits={fifty_eight.bit_length()}")
    print("first_capacity_not_excluding_repair=58")
    print(
        "conclusion=every_admissible_prefix_surviving_row_107_must_differ_"
        "from_rule_A_in_some_row_at_most_49"
    )
    print(
        "width_five_consequence=no_change_confined_to_rows_102..106_can_"
        "repair_the_rule_A_death_even_allowing_arbitrary_Lucas_box_blocks"
    )
    print(
        "boundary=the_length_58_capacity_is_only_an_upper_bound_and_does_"
        "not_construct_a_repair;earlier_or_repeated_local_moves_remain_open"
    )
    print(f"source_lf_sha256={sha256(source).hexdigest()}")


if __name__ == "__main__":
    main()
