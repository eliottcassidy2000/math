#!/usr/bin/env python3
"""Exact standard-library replay for THM-3588.

The script reconstructs the five failing R=512 Rule-A prefixes from
THM-3577.  For every proposed cut it builds a positive truncated-adjoint
Pascal multiplier.  The weighted cell inequalities telescope exactly, so a
negative right-hand side is an integer Farkas certificate.
"""

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
    require(five_pow_le_phi2m(d, m), f"golden lower bracket m={m}")
    require(not five_pow_le_phi2m(d + 1, m), f"golden upper bracket m={m}")
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


def kernel(current_degree: int, previous_degree: int) -> tuple[int, ...]:
    require(current_degree - previous_degree in (0, 1), "degree jump")
    return (1, 1) if current_degree == previous_degree else (1, 2, 1)


def clamp_halved(load: list[int], d: int, row: int) -> tuple[int, ...]:
    junk: list[int] = []
    for t, value in enumerate(load):
        lower = -comb(d - 1, t) if t < d else 0
        upper = comb(d - 1, t - 1) if t else 0
        correction = min(upper, max(lower, value))
        junk.append(value - correction)
    require(junk[d] == 0, f"Rule A survives row {row}")
    return tuple(junk[:d])


def baseline_before_death(
    r: int, profile: list[int], death: int
) -> tuple[list[tuple[int, ...]], list[list[int]]]:
    initial_full = t4_row_load(r, profile[0])
    require(all(value % 2 == 0 for value in initial_full), "initial parity")
    initial = [value // 2 for value in initial_full]
    rows = [clamp_halved(initial, profile[0], 0)]
    loads: list[list[int]] = [[] for _ in range(death + 1)]
    feed = [value // 2 for value in two_g_coeffs(r)]

    for i in range(1, death + 1):
        d, previous_degree = profile[i], profile[i - 1]
        transition = kernel(d, previous_degree)
        load = [0] * (d + 1)
        for t, value in enumerate(rows[-1]):
            for q, coefficient in enumerate(transition):
                load[t + q] += coefficient * value
        if d + i <= r - 1:
            load[0] += feed[d + i]
        if d > previous_degree and d - 1 + i <= r - 1:
            value = feed[d - 1 + i]
            load[0] += value
            load[1] += value
        loads[i] = load
        if i < death:
            rows.append(clamp_halved(load, d, i))

    require(len(rows) == death, "prefix length")
    return rows, loads


def certificate(
    profile: list[int],
    rows: list[tuple[int, ...]],
    loads: list[list[int]],
    death: int,
    cut: int,
) -> tuple[int, int, int, int]:
    """Return RHS, distinct cell rows, multiplier mass, max multiplier.

    Rows before ``cut`` are fixed.  The death lower-top inequality is

        -(K_death z_(death-1))_d <= load_A[death,d].

    Its coefficient vector is propagated backwards by the transposed Pascal
    kernels, truncating at each actual state width.  Positive multipliers of
    the cell inequalities

        z_(i,t) - (K_i z_(i-1))_t <= u_A(i,t)-lower(i,t)

    cancel every free difference coordinate exactly.
    """

    require(1 <= cut < death, "certificate cut range")
    d = profile[death]
    transition = kernel(d, profile[death - 1])
    level: dict[int, int] = {}
    coefficients: dict[tuple[int, int], int] = {}

    # The death inequality contributes minus the top load difference.
    for q, value in enumerate(transition):
        t = d - q
        if 0 <= t < profile[death - 1]:
            level[t] = level.get(t, 0) + value
            coefficients[death - 1, t] = (
                coefficients.get((death - 1, t), 0) - value
            )

    rhs = loads[death][d]
    distinct_rows = 0
    multiplier_mass = 0
    max_multiplier = 0

    for i in range(death - 1, cut - 1, -1):
        d_i = profile[i]
        distinct_rows += len(level)
        multiplier_mass += sum(level.values())
        max_multiplier = max(max_multiplier, *level.values())

        for t, multiplier in level.items():
            base_load = loads[i][t]
            base_state = rows[i][t]
            base_correction = base_load - base_state
            lower = -comb(d_i - 1, t)
            rhs += multiplier * (base_correction - lower)
            coefficients[i, t] = coefficients.get((i, t), 0) + multiplier

        if i > cut:
            transition = kernel(d_i, profile[i - 1])
            previous: dict[int, int] = {}
            for t, multiplier in level.items():
                for q, value in enumerate(transition):
                    u = t - q
                    if 0 <= u < profile[i - 1]:
                        amount = multiplier * value
                        previous[u] = previous.get(u, 0) + amount
                        coefficients[i - 1, u] = (
                            coefficients.get((i - 1, u), 0) - amount
                        )
            level = previous

    residue = {key: value for key, value in coefficients.items() if value}
    require(not residue, f"exact Farkas cancellation cut={cut}")
    require(multiplier_mass > 0 and max_multiplier > 0, "positive multipliers")
    return rhs, distinct_rows, multiplier_mass, max_multiplier


R = 512
EXPECTED = {
    0: (
        107,
        36,
        49,
        41086809435714651554305675711968187478323158540573926649141693358061153498971915247,
        -76874861671160752918838113835252414677773118569312018596735213281104563872967120639,
        -99590667299884280875845575989933573725647537922985471791423475015885221419855035910,
    ),
    1: (
        110,
        39,
        52,
        132432154406083499404833844855611470852311846908930910441711037452812262877844835930,
        -67299131366770904358812237887063009272989246608815495819036140457947571682105501209,
        -105451592315219564135171410757075758985778463073169264157182969265171538405612667542,
    ),
    2: (
        113,
        40,
        53,
        3338985983972956972751398645882828535963003571887262875884142255942013894214739377351,
        -9349426258946131715675734429741281552564767664183924828959904668783436277847131542484,
        -11822052866155320337773489992843051576765929222034076575485134567922224189637016676006,
    ),
    3: (
        116,
        43,
        56,
        12585302820800871834254225868388244300005120240474964547720968195252363811275385411024,
        -6065349518170092926836849962397980226126874773437073169325318976224259558021484964106,
        -9677470848752527118014869617911141365390390748769739837180875885264640701214066949271,
    ),
    4: (
        121,
        45,
        59,
        7877782680146233652648282350697035268590370134970386735210385525966333032831210948230500,
        -552094518178517855711445625652899715444284968095857944791483565103784726945679420728711,
        -2233996250002381223173213859052295611843410320963231167526252801295528998327804201754981,
    ),
}

print("THM-3588 AMM R=512 truncated-adjoint Pascal repair horizons")
for offset, expected in EXPECTED.items():
    death, first_negative, old_bound, expected_previous, expected_cut, expected_fatal = expected
    profile = [floor_gamma_star(R + i) + offset for i in range(R)]
    rows, loads = baseline_before_death(R, profile, death)
    require(loads[death][profile[death]] == expected_fatal, f"fatal D0={offset}")

    ledger = [certificate(profile, rows, loads, death, cut) for cut in range(1, death)]
    negative_cuts = [cut for cut, data in enumerate(ledger, start=1) if data[0] < 0]
    require(negative_cuts == list(range(first_negative, death)), f"single sign wall D0={offset}")

    previous = ledger[first_negative - 2]
    active = ledger[first_negative - 1]
    require(previous[0] == expected_previous > 0, f"positive boundary D0={offset}")
    require(active[0] == expected_cut < 0, f"negative certificate D0={offset}")

    new_bound = first_negative - 1
    gain = old_bound - new_bound
    print(
        f"D0={offset} death={death} first_negative_cut={first_negative} "
        f"must_differ_by={new_bound} old_bound={old_bound} gain={gain} "
        f"B_prev={previous[0]} B_cut={active[0]} "
        f"rows={active[1]} mass={active[2]} max_multiplier={active[3]}"
    )

print("new must-differ-by vector: (35,38,39,42,44)")
print("scope: exact fixed-R fixed-policy prefix obstructions; no feasibility or uniform AMM conclusion")
print("all optimization-safe exact truth gates passed")
