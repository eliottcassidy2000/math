#!/usr/bin/env python3
"""Integer-exact referee audit for THM-2060."""

from math import gcd, lcm


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def centered_abs(x, modulus):
    residue = x % modulus
    return min(residue, modulus - residue)


def ceil_div(numerator, denominator):
    return (numerator + denominator - 1) // denominator


def core_packet(core, clock):
    return [
        r
        for r in range(clock)
        if all(14 * centered_abs(c * r, clock) >= clock for c in core)
    ]


def tail_data(a, w, clock):
    full_clock = a * clock
    reduced_gcd = gcd(w, full_clock)
    h = full_clock // reduced_gcd
    u = w // reduced_gcd
    d = gcd(clock, h)
    beta = [0] * d
    for s in range(h):
        if 14 * centered_abs(u * s, h) >= h:
            beta[s % d] += 1
    return h, d, beta


def packet_count(core, a, w, clock):
    packet = core_packet(core, clock)
    h, d, beta = tail_data(a, w, clock)
    alpha = [0] * d
    for r in packet:
        alpha[r % d] += 1
    packet_types = sum(x * y for x, y in zip(alpha, beta))
    grid_count = gcd(a, w) * packet_types
    return packet, beta, packet_types, grid_count, h, d


def direct_grid_count(core, a, w, clock):
    full_clock = a * clock
    return sum(
        all(14 * centered_abs(c * k, clock) >= clock for c in core)
        and 14 * centered_abs(w * k, full_clock) >= full_clock
        for k in range(full_clock)
    )


def exceptional_residues(core, clock):
    packet = core_packet(core, clock)
    return [
        b
        for b in range(clock)
        if all(14 * centered_abs(b * r, clock) < clock for r in packet)
    ]


def main():
    algebra_rows = 0
    for a in range(1, 101):
        for w in range(1, 161):
            for clock in range(1, 81):
                full_clock = a * clock
                h = full_clock // gcd(w, full_clock)
                d = gcd(clock, h)
                q = a // gcd(a, w)
                require(h // d == q, ("q identity", a, w, clock))
                require(lcm(clock, h) == clock * q, ("lcm identity", a, w, clock))
                require(full_clock // lcm(clock, h) == gcd(a, w),
                        ("lift identity", a, w, clock))
                algebra_rows += 1

    histogram_bins = 0
    sharp_hits = set()
    for h in range(1, 211):
        for d in range(1, h + 1):
            if h % d:
                continue
            q = h // d
            lower = q - ceil_div(q, 7)
            for u in range(1, h + 1):
                if gcd(u, h) != 1:
                    continue
                beta = [0] * d
                for s in range(h):
                    if 14 * centered_abs(u * s, h) >= h:
                        beta[s % d] += 1
                for value in beta:
                    require(value >= lower, ("histogram lower bound", h, d, u, beta))
                    histogram_bins += 1
                    if value == lower:
                        sharp_hits.add(q)
                require((min(beta) > 0) == (q >= 2),
                        ("full support iff", h, d, u, beta))

    cores = (
        tuple(range(1, 14)),
        tuple(range(1, 12)) + (13,),
        (1, 3, 4, 7, 9),
        (2, 5, 11, 17),
    )
    packet_rows = 0
    positive_noncoprime_rows = 0
    for core in cores:
        for clock in range(2, 31):
            for a in range(1, 13):
                for w in range(1, 61):
                    packet, beta, packet_types, grid_count, _, _ = packet_count(
                        core, a, w, clock
                    )
                    q = a // gcd(a, w)
                    lower = q - ceil_div(q, 7)
                    require(packet_types >= lower * len(packet),
                            ("packet lower bound", core, clock, a, w))
                    require(grid_count == direct_grid_count(core, a, w, clock),
                            ("direct grid identity", core, clock, a, w))
                    if packet and a % gcd(a, w) != 0:
                        raise RuntimeError("impossible gcd divisibility state")
                    if packet and w % a and gcd(a, w) > 1:
                        require(packet_types > 0 and min(beta) > 0,
                                ("proper noncoprime tail", clock, a, w))
                        positive_noncoprime_rows += 1
                    packet_rows += 1

    two_tail_pairs = 0
    equality_pairs = []
    for q_1 in range(2, 501):
        for q_2 in range(2, 501):
            left_num = ceil_div(q_1, 7) * q_2 + ceil_div(q_2, 7) * q_1
            right_num = q_1 * q_2
            require(left_num <= right_num, ("two-tail capacity", q_1, q_2))
            if left_num == right_num:
                equality_pairs.append((q_1, q_2))
            two_tail_pairs += 1
    require(equality_pairs == [(2, 2)], ("dyadic equality", equality_pairs))

    small_clock_zero_sets = 0
    for clock in range(2, 15):
        core = tuple(c for c in range(1, 2 * clock + 1) if c % clock)
        require(exceptional_residues(core, clock) == [0],
                ("small-clock zero set", clock))
        small_clock_zero_sets += 1

    boundary_core = tuple(range(1, 14))
    require(core_packet(boundary_core, 14) == [1, 3, 5, 9, 11, 13],
            "AP boundary packet")
    require(exceptional_residues(boundary_core, 14) == [0],
            "AP boundary exceptional set")

    guardrail_packet, _, guardrail_types, guardrail_grid, _, _ = packet_count(
        tuple(range(1, 13)), 2, 26, 14
    )
    require(guardrail_packet and guardrail_types > 0 and guardrail_grid > 0,
            "a|w does not force a zero packet")
    require(
        all(14 * centered_abs(v, 28) >= 28 for v in range(2, 27, 2)),
        "t=1/28 guardrail witness",
    )

    require(exceptional_residues((1,), 2) == [0], "bulk sieve residue")
    require(
        min(centered_abs(1, 3), centered_abs(2, 3)) * 14 >= 3,
        "bulk sieve is not complete",
    )

    print("THM-2060 CRT TAIL-COSET SATURATION AUDIT")
    print("valuation/lcm identity rows checked:", algebra_rows)
    print("tail histogram bins checked:", histogram_bins)
    print("sharp q values observed through h<=210:", len(sharp_hits))
    print("packet/direct-grid rows checked:", packet_rows)
    print("proper noncoprime positive rows:", positive_noncoprime_rows)
    print("two-tail order pairs checked:", two_tail_pairs)
    print("two-tail equality pairs:", equality_pairs)
    print("small full-unit exceptional sets checked:", small_clock_zero_sets)
    print("boundary AP exceptional residues mod 14: [0]")
    print("guardrails: q=1 converse rejected; finite bulk clock iff rejected")
    print("carrier: common lift-sheet cover / CRT reduction histograms")
    print("tournament: rejected except binary ownership gauge at q1=q2=2")
    print("PASS")


if __name__ == "__main__":
    main()
