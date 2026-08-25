#!/usr/bin/env python3
"""Exact companion for THM-4055's C_60/dyadic phase-fibre law.

The universal result is proved arithmetically in the theorem.  This companion
exhausts every clock state for the dyadic periods occurring in THM-4047,
checks the canonical-lift sidecar and carry law, imports the frozen bank census
by hash, and reconstructs the hostile ell_29 word from the physical orbit.
"""

from __future__ import annotations

from collections import defaultdict
from hashlib import sha256
from math import gcd, lcm
from pathlib import Path


PERIODS = (1, 2, 4, 8, 16, 32)
HISTOGRAM = ((1, 16), (2, 10), (4, 56), (8, 668), (16, 87118), (32, 12133))
PARENT_OUTPUT_SHA256 = "bae500127999c350ebeff77c3145fbb4abf7b5b7292d2882eeecc6492fba75e3"
ELL29_START = 90
ELL29_ABSOLUTE_WORD = (1, 1, 0, 1, 0, 0, 1, 0)


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def read_parent_output() -> str:
    root = Path(__file__).resolve().parents[1]
    path = root / "05-knowledge/results/rule30_left_front_affine_monodromy_thm4047.out"
    raw = path.read_bytes()
    require(sha256(raw).hexdigest() == PARENT_OUTPUT_SHA256, "THM-4047 output hash")
    text = raw.decode("ascii")
    require(f"period_histogram={HISTOGRAM}" in text, "parent histogram")
    require("(29, 90, 0, (1, 1))" in text, "parent ell_29 phase query")
    return text


def audit_phase_fibre(p: int) -> tuple[int, int, int, int]:
    g = gcd(60, p)
    m = p // g
    clock = lcm(60, p)
    require(clock == 60 * m, (p, "clock order"))

    fibres: dict[int, list[tuple[int, int]]] = defaultdict(list)
    seen = set()
    for t in range(clock):
        r = t % 60
        a = t % p
        q = ((t - r) // 60) % m
        require((r, a) not in seen, (p, "fibre-product injectivity", r, a))
        seen.add((r, a))
        fibres[r].append((a, q))
        require((r + 60 * q) % p == a, (p, "sidecar reconstruction", t))

    expected = {(r, a) for r in range(60) for a in range(p) if (r - a) % g == 0}
    require(seen == expected, (p, "fibre-product image"))
    for r in range(60):
        require(len(fibres[r]) == m, (p, "fibre cardinality", r))
        require({q for _, q in fibres[r]} == set(range(m)), (p, "sidecar bijection", r))

    # Exhaust the representative carry law.  It records why the torsor chart
    # is not a direct-product group chart when m>1.
    for r in range(60):
        for s in range(60):
            carry = (r + s) // 60
            for q in range(m):
                for u in range(m):
                    lhs = (r + 60 * q + s + 60 * u) % clock
                    rhs = ((r + s) % 60 + 60 * ((q + u + carry) % m)) % clock
                    require(lhs == rhs, (p, "carry law", r, s, q, u))

    # A least-p indicator word gives an exact hostile for the factorization
    # test.  The theorem proves the iff for every least-p response.
    word = tuple(int(a == 0) for a in range(p))
    phase_only = all(
        len({word[t % p] for t in range(r, clock, 60)}) == 1 for r in range(60)
    )
    require(phase_only == (p in (1, 2, 4)), (p, "phase-only iff"))
    return p, g, m, clock


def physical_ell29() -> tuple[int, ...]:
    # Bit r of front is ell_r(t).  The strip r<=29 is closed under the exact
    # left-front recurrence, so no discarded bit can influence this check.
    max_r = 29
    mask = (1 << (max_r + 1)) - 1
    front = 1
    values: dict[int, int] = {}
    states: dict[int, int] = {}
    for t in range(ELL29_START + 121):
        if t in (90, 94, 98):
            states[t] = front
        if t >= ELL29_START:
            values[t] = (front >> max_r) & 1
        shifted = front << 1
        front = ((shifted << 1) ^ shifted ^ front ^ (shifted & front)) & mask

    word = tuple(values[t] for t in range(96, 104))
    require(word == ELL29_ABSOLUTE_WORD, "ell_29 absolute word")
    require(states[90] == states[98] and states[90] != states[94], "closed-strip least period eight")
    for t in range(ELL29_START, ELL29_START + 61):
        require(values[t + 8] == values[t], ("ell_29 period eight", t))
        require(values[t + 4] == 1 - values[t], ("ell_29 anti-period four", t))
        require(values[t + 60] == 1 - values[t], ("ell_29 shift sixty", t))
    return word


def main() -> None:
    read_parent_output()
    fibres = tuple(audit_phase_fibre(p) for p in PERIODS)
    word = physical_ell29()

    histogram = dict(HISTOGRAM)
    phase_only_columns = sum(histogram[p] for p in (1, 2, 4))
    hidden_columns = sum(histogram[p] for p in (8, 16, 32))
    require(phase_only_columns == 82, "phase-only column count")
    require(hidden_columns == 99_919, "hidden-sidecar column count")
    require(phase_only_columns + hidden_columns == 100_001, "bank size")
    sidecar = tuple((p, p // gcd(p, 60), (p // gcd(p, 60)).bit_length() - 1) for p in PERIODS)

    # Least period eight alone does not make every clock fibre faithful.
    hostile = (1, 0, 0, 0, 0, 0, 0, 0)
    clock = lcm(60, 8)
    r1_values = {hostile[t % 8] for t in range(1, clock, 60)}
    r0_values = {hostile[t % 8] for t in range(0, clock, 60)}
    require(r1_values == {0} and r0_values == {0, 1}, "response-value hostile")

    semantic = (fibres, HISTOGRAM, phase_only_columns, hidden_columns, sidecar, word)
    digest = sha256(repr(semantic).encode("ascii")).hexdigest()
    print("SIXTY_DYADIC_RESPONSE_FIBRE_THM4055")
    print("universe=C_60 paired with p in (1,2,4,8,16,32);Rule30 fixed columns r=0..100000")
    print(f"phase_fibres_p_g_m_lcm={fibres}")
    print(f"period_histogram={HISTOGRAM}")
    print(f"phase_only_columns={phase_only_columns}")
    print(f"hidden_sidecar_columns={hidden_columns}")
    print(f"sidecar_states_and_bits={sidecar}")
    print("ell29_closed_strip_repeat=(time=90,period=8,least_against4=True)")
    print(f"ell29=(start={ELL29_START},absolute_word={word},shift60=complement)")
    print("hostile_period8_fibres=(r0=faithful,r1=constant)")
    print(f"semantic_sha256={digest}")
    print("scope=PROVED phase-fibre law;FINITE-EXACT verified fixed-bank consequence;moving observers OPEN")


if __name__ == "__main__":
    main()
