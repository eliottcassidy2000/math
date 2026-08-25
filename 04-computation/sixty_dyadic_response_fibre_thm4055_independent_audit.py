#!/usr/bin/env python3
"""Independent audit for THM-4055.

This path enumerates compatible CRT pairs rather than a combined time orbit,
uses modular inversion for the sidecar, reads THM-4047's independent output,
and evolves ordinary spatial Rule 30 rows rather than a packed front strip.
"""

from __future__ import annotations

from hashlib import sha256
from math import gcd, lcm
from pathlib import Path


PERIODS = (1, 2, 4, 8, 16, 32)
HISTOGRAM = ((1, 16), (2, 10), (4, 56), (8, 668), (16, 87118), (32, 12133))
PARENT_AUDIT_SHA256 = "4dd45b0073b15d18dfe6bc68c976fecf09fe7d86055c3bf673449b84bb521c1f"
ELL29_START = 90
ELL29_ABSOLUTE_WORD = (1, 1, 0, 1, 0, 0, 1, 0)


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def audit_parent() -> None:
    root = Path(__file__).resolve().parents[1]
    path = root / "05-knowledge/results/rule30_left_front_affine_monodromy_thm4047_independent_audit.out"
    raw = path.read_bytes()
    require(sha256(raw).hexdigest() == PARENT_AUDIT_SHA256, "THM-4047 audit hash")
    text = raw.decode("ascii")
    require(f"period_histogram={HISTOGRAM}" in text, "independent parent histogram")
    require("period_doubling_columns=(3, 8, 29, 400, 87867)" in text, "ell_29 doubling")


def compatible_pair_audit(p: int) -> tuple[int, int, int, int]:
    g = gcd(60, p)
    m = p // g
    clock = lcm(60, p)
    pairs = [(r, a) for r in range(60) for a in range(p) if (r - a) % g == 0]
    require(len(pairs) == clock, (p, "fibre product size"))

    for r in range(60):
        fibre = [a for rr, a in pairs if rr == r]
        require(len(fibre) == m, (p, "fixed phase size", r))
        recovered = set()
        for a in fibre:
            if m == 1:
                q = 0
            else:
                q = (((a - r) // g) * pow(60 // g, -1, m)) % m
            require((r + 60 * q) % p == a, (p, "inverse sidecar", r, a))
            recovered.add(q)
        require(recovered == set(range(m)), (p, "minimal torsor coordinate", r))

    # The direct-product carry formula is checked in integer representatives.
    for r in range(60):
        for s in range(60):
            c = int(r + s >= 60)
            for q in range(m):
                for u in range(m):
                    source = (r + 60 * q + s + 60 * u) % clock
                    chart = ((r + s) % 60 + 60 * ((q + u + c) % m)) % clock
                    require(source == chart, (p, "carry chart"))
    return p, g, m, clock


def ordinary_rule30_ell29() -> tuple[int, ...]:
    row = {0: 1}
    values: dict[int, int] = {}
    states: dict[int, tuple[int, ...]] = {}
    last_time = ELL29_START + 120
    for t in range(last_time + 1):
        if t in (90, 94, 98):
            states[t] = tuple(row.get(-t + r, 0) for r in range(30))
        if t >= ELL29_START:
            values[t] = row.get(-t + 29, 0)
        if t == last_time:
            break
        next_row: dict[int, int] = {}
        for j in range(-t - 1, t + 2):
            value = row.get(j - 1, 0) ^ (row.get(j, 0) | row.get(j + 1, 0))
            if value:
                next_row[j] = 1
        row = next_row

    word = tuple(values[t] for t in range(96, 104))
    require(word == ELL29_ABSOLUTE_WORD, "ordinary ell_29 word")
    require(states[90] == states[98] and states[90] != states[94], "closed-strip least period eight")
    require(all(values[t + 60] == 1 - values[t] for t in range(90, 151)), "shift-sixty control")
    return word


def main() -> None:
    audit_parent()
    fibres = tuple(compatible_pair_audit(p) for p in PERIODS)
    word = ordinary_rule30_ell29()
    histogram = dict(HISTOGRAM)
    phase_only_columns = sum(histogram[p] for p in PERIODS if 60 % p == 0)
    hidden_columns = sum(histogram[p] for p in PERIODS if 60 % p != 0)
    sidecar = tuple((p, p // gcd(p, 60), (p // gcd(p, 60)).bit_length() - 1) for p in PERIODS)
    require((phase_only_columns, hidden_columns) == (82, 99_919), "bank split")

    # Explicitly audit the response-value caveat with a least-period-eight word.
    hostile = (1, 0, 0, 0, 0, 0, 0, 0)
    require({hostile[t % 8] for t in range(1, 120, 60)} == {0}, "constant hostile fibre")
    require({hostile[t % 8] for t in range(0, 120, 60)} == {0, 1}, "faithful hostile fibre")

    semantic = (fibres, HISTOGRAM, phase_only_columns, hidden_columns, sidecar, word)
    digest = sha256(repr(semantic).encode("ascii")).hexdigest()
    print("SIXTY_DYADIC_RESPONSE_FIBRE_THM4055_INDEPENDENT_AUDIT")
    print("method=compatible-pair CRT inversion;ordinary spatial Rule30;independent THM4047 output")
    print(f"phase_fibres_p_g_m_lcm={fibres}")
    print(f"period_histogram={HISTOGRAM}")
    print(f"phase_only_columns={phase_only_columns}")
    print(f"hidden_sidecar_columns={hidden_columns}")
    print(f"sidecar_states_and_bits={sidecar}")
    print("ell29_closed_strip_repeat=(time=90,period=8,least_against4=True)")
    print(f"ell29=(start={ELL29_START},absolute_word={word},shift60=complement)")
    print("hostile_period8_fibres=(r0=faithful,r1=constant)")
    print(f"semantic_sha256={digest}")
    print("audit=PASS;no primary import;response values separated from phase states")


if __name__ == "__main__":
    main()
