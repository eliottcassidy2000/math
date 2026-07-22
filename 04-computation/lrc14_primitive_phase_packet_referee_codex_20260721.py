#!/usr/bin/env python3
"""Exact referee for THM-2058 primitive packets and deck/fan intervals."""

from __future__ import annotations

from math import gcd, lcm


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def least_abs(x: int, n: int) -> int:
    r = x % n
    return min(r, n - r)


def safe_packet(m: tuple[int, ...], n: int) -> tuple[int, ...]:
    return tuple(
        x
        for x in range(n)
        if 14 * min(least_abs(x * a, n) for a in m) >= n
    )


def primitive_packet(m: tuple[int, ...], n: int) -> tuple[int, ...]:
    return tuple(x for x in safe_packet(m, n) if gcd(x, n) == 1)


def divisors(n: int) -> tuple[int, ...]:
    small: list[int] = []
    large: list[int] = []
    d = 1
    while d * d <= n:
        if n % d == 0:
            small.append(d)
            if d * d != n:
                large.append(n // d)
        d += 1
    return tuple(small + large[::-1])


def mobius(n: int) -> int:
    value = 1
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            value = -value
            if n % p == 0:
                return 0
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        value = -value
    return value


def deck_numerator(m: tuple[int, ...], n: int) -> int:
    return max(min(least_abs(x * a, n) for a in m) for x in range(n))


def reduced_lift_union(m: tuple[int, ...], n: int) -> set[int]:
    return {
        ((n // d) * a) % n
        for d in divisors(n)
        for a in primitive_packet(m, d)
    }


def row_packet(m: tuple[int, ...], n: int, unit: int) -> tuple[int, ...]:
    return tuple(
        ell
        for ell in range(n)
        if 14 * min(least_abs(ell * unit * a, n) for a in m) >= n
    )


def one_tail_raw_interval(n: int, r: int) -> tuple[int, ...]:
    """Integer M satisfying positivity, strict -c_13 ownership, and gate failure."""
    values: list[int] = []
    for longitudinal in range(1, (n - 1) // r + 1):
        a = longitudinal
        b = n - (r + 12) * longitudinal
        tail = n - r * longitudinal
        if b <= 0 or tail <= 0:
            continue
        owner_value = 13 * b
        if owner_value <= abs(a - 12 * b):
            continue
        if a * a + b * b >= 91 * owner_value:
            continue
        values.append(longitudinal)
    return tuple(values)


def coprime_floor_count(n: int, lo: int, hi: int) -> int:
    return sum(
        mobius(e) * (hi // e - (lo - 1) // e)
        for e in divisors(n)
    )


def main() -> None:
    ap12 = (1, -1, *range(2, 13))
    templates = (
        (1, -1, 2, 4, 7),
        ap12,
        (1, -2, 3, 5, 8),
    )

    decomposition_checks = 0
    inversion_checks = 0
    for m in templates:
        for n in range(1, 241):
            packet = set(safe_packet(m, n))
            require(packet == reduced_lift_union(m, n), f"packet split m={m}, N={n}")
            primitive_count = len(primitive_packet(m, n))
            inverted = sum(
                mobius(n // d) * len(safe_packet(m, d)) for d in divisors(n)
            )
            require(primitive_count == inverted, f"Mobius inversion m={m}, N={n}")
            decomposition_checks += 1
            inversion_checks += 1

    m_beatty = (1, -1, 2)
    beatty_q = lcm(*(14 * abs(a) for a in m_beatty))
    beatty_shifts = {
        len(safe_packet(m_beatty, n + beatty_q))
        - len(safe_packet(m_beatty, n))
        for n in range(1, beatty_q + 1)
    }
    require(len(beatty_shifts) == 1, "Beatty shift law")
    beatty_shift = beatty_shifts.pop()

    require(primitive_packet(ap12, 13) == tuple(range(1, 13)), "AP p13")
    require(primitive_packet(ap12, 27) == (2, 25), "AP p27")
    require(primitive_packet(ap12, 351) == (), "AP p351")
    require(row_packet(ap12, 27, 1) == (2, 25), "unit transport M=1")
    require(row_packet(ap12, 27, 2) == (1, 26), "unit transport M=2")
    for x in safe_packet(ap12, 27):
        ell = (pow(2, -1, 27) * x) % 27
        template_vector = tuple(least_abs(x * a, 27) for a in ap12)
        row_vector = tuple(least_abs(ell * 2 * a, 27) for a in ap12)
        require(template_vector == row_vector, f"label transport x={x}")

    crt_template = (1, -1, 2, 4, 7)
    require(
        tuple(map(len, (primitive_packet(crt_template, 3),
                        primitive_packet(crt_template, 5),
                        primitive_packet(crt_template, 15)))) == (2, 4, 0),
        "CRT no-go",
    )

    one_tail_values = {
        "D34_r10": deck_numerator((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -10, 13), 34),
        "D37_r13": deck_numerator((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -13, 13), 37),
        "D41_r3": deck_numerator((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -3, 13), 41),
        "D48_r10": deck_numerator((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -10, 13), 48),
        "D51_r13": deck_numerator((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -13, 13), 51),
        "D113_r5": deck_numerator((1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, -5, 13), 113),
    }
    require(
        one_tail_values == {
            "D34_r10": 2,
            "D37_r13": 2,
            "D41_r3": 3,
            "D48_r10": 4,
            "D51_r13": 3,
            "D113_r5": 9,
        },
        "one-tail deck values",
    )

    interval_checks = 0
    interval_fingerprints: list[str] = []
    for n, r in ((34, 10), (37, 13), (41, 3), (48, 10), (51, 13), (113, 5)):
        raw = one_tail_raw_interval(n, r)
        if raw:
            require(raw == tuple(range(raw[0], raw[-1] + 1)), f"fan interval N={n}, r={r}")
            direct = sum(gcd(n, value) == 1 for value in raw)
            require(
                direct == coprime_floor_count(n, raw[0], raw[-1]),
                f"coprime floor count N={n}, r={r}",
            )
            interval_fingerprints.append(f"{n}:{r}:{raw[0]}-{raw[-1]}:{direct}")
        else:
            interval_fingerprints.append(f"{n}:{r}:empty:0")
        interval_checks += 1

    require(deck_numerator(ap12, 15) == 1, "AP D15")
    require(safe_packet(ap12, 15) == (), "AP N15 safe packet")

    hasse_checks = 0
    core = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13)
    for k in range(6):
        seven_power = 7**k
        u = pow(seven_power, -1, 113)
        w = 24 + 84 * seven_power * u
        require(w % 14 == 24 % 14, f"Hasse mod14 k={k}")
        require(w % 12 == 0, f"Hasse q12 k={k}")
        require(w % 113 == 108, f"Hasse q113 k={k}")
        require(((w - 24) // 14) % seven_power == 0, f"Hasse lift k={k}")
        distances = tuple(least_abs(47 * speed, 113) for speed in (*core, w))
        require(min(distances) == 9, f"Hasse exit k={k}")
        hasse_checks += 1

    for w in (24, 38, 108):
        b = w - 12
        determinant = max(13 * abs(b), abs(1 - 12 * b))
        require(determinant == 13 * b, f"one-tail owner w={w}")
        require(1 + b * b < 91 * determinant, f"one-tail gate failure w={w}")

    print("THM-2058 PRIMITIVE PHASE-PACKET REFEREE")
    print(f"packet_decomposition_checks={decomposition_checks}")
    print(f"mobius_inversion_checks={inversion_checks}")
    print(f"beatty_Q={beatty_q} beatty_shift={beatty_shift}")
    print("AP12_primitive_packets=p13:12,p27:2,p351:0")
    print("unit_transport_N27=M1:{2,25},M2:{1,26}")
    print("CRT_no_go=p3:2,p5:4,p15:0")
    print("one_tail_deck_numerators=" + ",".join(f"{k}:{v}" for k, v in one_tail_values.items()))
    print("fan_interval_fingerprints=" + ",".join(interval_fingerprints))
    print(f"fan_interval_checks={interval_checks}")
    print(f"same_sector_hasse_checks={hasse_checks} denominator=113 slack=13")
    print("spanning_tree_no_go=AP12,N15,D_numerator:1,safe_packet:empty")
    print("TOURNAMENT ANALYSIS=not applicable: phase-order packets and determinant magnitudes are load-bearing")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
