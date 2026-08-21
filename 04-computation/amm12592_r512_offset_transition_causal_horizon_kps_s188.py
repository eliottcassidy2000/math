#!/usr/bin/env python3
"""Exact standard-library replay for THM-3577."""

from __future__ import annotations

from functools import lru_cache
from math import comb


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


@lru_cache(maxsize=None)
def fib_pair(n: int) -> tuple[int, int]:
    if n == 0:
        return 0, 1
    a, b = fib_pair(n >> 1)
    c = a * (2 * b - a)
    d = a * a + b * b
    return (d, c + d) if n & 1 else (c, d)


def five_pow_le_phi2m(d: int, m: int) -> bool:
    if d < 0:
        return True
    f, f_next = fib_pair(2 * m)
    lucas = 2 * f_next - f
    gap = 2 * 5**d - lucas
    return gap <= 0 or gap * gap < 5 * f * f


def floor_gamma_star(m: int) -> int:
    d = 3 * m // 5
    while five_pow_le_phi2m(d + 1, m):
        d += 1
    while not five_pow_le_phi2m(d, m):
        d -= 1
    require(five_pow_le_phi2m(d, m), f"golden floor lower bracket m={m}")
    require(not five_pow_le_phi2m(d + 1, m), f"golden floor upper bracket m={m}")
    return d


def two_g_coeffs(r: int) -> list[int]:
    values = [2]
    coefficient = 1
    for j in range(1, r):
        coefficient = coefficient * (r - j) // j
        values.append((-coefficient if j & 1 else coefficient) - 1)
    require(all(value % 2 == 0 for value in values), "feed parity")
    return values


def t4_row_load(r: int, d: int) -> list[int]:
    return [
        (-1) ** (d - t) * comb(r - 2 - t, d - t)
        - comb(d + 1, t + 1)
        + 2 * comb(d, t)
        for t in range(d + 1)
    ]


def first_feed_free_row(r: int, profile: list[int]) -> int:
    for i in range(1, r - 1):
        d, previous = profile[i], profile[i - 1]
        first_feed_absent = d + i > r - 1
        second_feed_present = d > previous and d - 1 + i <= r - 1
        if first_feed_absent and not second_feed_present:
            return i
    raise RuntimeError("no feed-free row")


def clamp(load: list[int], d: int) -> tuple[tuple[int, ...], int]:
    junk = []
    for t, value in enumerate(load):
        lower = -2 * comb(d - 1, t) if t < d else 0
        upper = 2 * comb(d - 1, t - 1) if t else 0
        correction = min(upper, max(lower, value))
        residual = value - correction
        require(residual % 2 == 0, f"residual parity d={d},t={t}")
        junk.append(residual)
    return tuple(value // 2 for value in junk[:d]), junk[d]


def rule_a_until(r: int, profile: list[int], stop: int) -> tuple[list[tuple[int, ...]], tuple[int, int] | None]:
    initial = t4_row_load(r, profile[0])
    require(all(value % 2 == 0 for value in initial), "initial parity")
    row, top = clamp(initial, profile[0])
    if top:
        return [], (0, top)
    rows = [row]
    feed = two_g_coeffs(r)
    for i in range(1, stop):
        d, previous_degree = profile[i], profile[i - 1]
        kernel = (1, 1) if d == previous_degree else (1, 2, 1)
        require(d - previous_degree in (0, 1), f"degree jump row={i}")
        previous = tuple(2 * value for value in rows[-1])
        load = [0] * (d + 1)
        for t, value in enumerate(previous):
            for q, coefficient in enumerate(kernel):
                load[t + q] += coefficient * value
        if d + i <= r - 1:
            load[0] += feed[d + i]
        if d > previous_degree and d - 1 + i <= r - 1:
            value = feed[d - 1 + i]
            load[0] += value
            load[1] += value
        row, top = clamp(load, d)
        if top:
            return rows, (i, top)
        rows.append(row)
    return rows, None


def entry_report(entry: tuple[int, ...], entry_degree: int) -> tuple[int, int, int]:
    support = [t for t, value in enumerate(entry) if value]
    support_max = max(support) if support else -1
    require(all(value <= 0 for value in entry), "F1 sign")
    require(support_max + 2 < entry_degree, "F1 support")
    f2_margin = entry_degree - 1 + 2 * entry[0]
    require(f2_margin >= 0, "F2")

    def y(t: int) -> int:
        return entry[t] if 0 <= t < len(entry) else 0

    margins = [
        comb(entry_degree - 1, t) + 2 * y(t - 1) + y(t - 2)
        for t in range(2, support_max + 3)
    ]
    f3_margin = min(margins, default=0)
    require(f3_margin >= 0, "F3")
    return support_max, f2_margin, f3_margin


def causal_capacity(profile: list[int], death: int, horizon: int) -> int:
    return sum(
        2 ** (s + 1) * comb(profile[death - s], s)
        for s in range(1, horizon + 1)
    )


R = 512
expected_failures = {
    0: (107, 58, 49),
    1: (110, 58, 52),
    2: (113, 60, 53),
    3: (116, 60, 56),
    4: (121, 62, 59),
}
rows_out: list[str] = []

for offset in range(6):
    profile = [floor_gamma_star(R + i) + offset for i in range(R)]
    entry = first_feed_free_row(R, profile)
    rows, failure = rule_a_until(R, profile, entry)
    if offset <= 4:
        require(failure is not None, f"expected death D0={offset}")
        death, residual = failure
        require(residual < 0, f"negative fatal residual D0={offset}")
        required = -residual
        horizon = 1
        while causal_capacity(profile, death, horizon) < required:
            horizon += 1
        last_agreement = death - horizon
        require(
            (death, horizon, last_agreement) == expected_failures[offset],
            f"failure atlas D0={offset}",
        )
        lower_capacity = causal_capacity(profile, death, horizon - 1)
        upper_capacity = causal_capacity(profile, death, horizon)
        require(lower_capacity < required <= upper_capacity, f"capacity bracket D0={offset}")
        rows_out.append(
            f"D0={offset} entry={entry} state_degree={profile[entry - 1]} "
            f"next_degree={profile[entry]} "
            f"death={death} fatal={required} horizon={horizon} "
            f"must_differ_by={last_agreement} C_prev={lower_capacity} C_hit={upper_capacity}"
        )
    else:
        require(failure is None and len(rows) == entry, "D0=5 survives to entry")
        report = entry_report(rows[-1], profile[entry])
        require(
            (entry, profile[entry - 1], profile[entry], *report)
            == (127, 386, 387, 18, 14, 1229),
            f"D0=5 entry atlas got={(entry, profile[entry - 1], profile[entry], *report)}",
        )
        rows_out.append(
            f"D0=5 entry={entry} state_degree={profile[entry - 1]} "
            f"next_degree={profile[entry]} survive=True "
            f"support_max={report[0]} F2_margin={report[1]} F3_min_margin={report[2]}"
        )

print("THM-3577 AMM R=512 offset transition and causal horizon")
for row in rows_out:
    print(row)
print("transition: Rule A fails for D0=0..4 and enters the certified F1-F3 cone at D0=5")
print("scope: exact policy atlas only; no alternative-prefix infeasibility and no uniform sub-two extractor")
print("all active truth gates passed")
