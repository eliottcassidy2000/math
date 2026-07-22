#!/usr/bin/env python3
"""Exact audit for THM-2059's CRT safe-phase packet formula."""

from math import gcd, lcm


def centered_abs(x, modulus):
    residue = x % modulus
    return min(residue, modulus - residue)


def packets(core, a, w, clock):
    q = clock * a
    g = gcd(w, q)
    h = q // g
    u = w // g
    d = gcd(clock, h)
    core_packet = [
        r
        for r in range(clock)
        if all(14 * centered_abs(c * r, clock) >= clock for c in core)
    ]
    tail_packet = [
        s for s in range(h) if 14 * centered_abs(u * s, h) >= h
    ]
    alpha = [0] * d
    beta = [0] * d
    for r in core_packet:
        alpha[r % d] += 1
    for s in tail_packet:
        beta[s % d] += 1
    packet_types = sum(x * y for x, y in zip(alpha, beta))
    predicted_grid_count = (q // lcm(clock, h)) * packet_types
    return predicted_grid_count, packet_types, alpha, beta


def direct_grid_count(core, a, w, clock):
    q = clock * a
    count = 0
    for k in range(q):
        core_safe = all(
            14 * centered_abs(c * k, clock) >= clock for c in core
        )
        tail_safe = 14 * centered_abs(w * k, q) >= q
        count += core_safe and tail_safe
    return count


def main():
    cores = (
        tuple(range(1, 14)),
        tuple(range(1, 12)) + (13,),
        tuple(range(1, 26, 2)),
        (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41),
    )
    rows = 0
    grid_indices = 0
    positive = 0
    zero = 0
    for clock in range(2, 26):
        for a in range(1, 8):
            for w in range(1, 81):
                for core in cores:
                    predicted, packet_types, _, _ = packets(core, a, w, clock)
                    direct = direct_grid_count(core, a, w, clock)
                    assert predicted == direct
                    rows += 1
                    grid_indices += clock * a
                    if packet_types:
                        positive += 1
                    else:
                        zero += 1

    missing_clock_rows = 0
    for clock in range(2, 15):
        core = tuple(c for c in range(1, 29) if c % clock)
        for a in range(1, 13):
            for w in range(1, 301):
                if w % (clock * a) == 0:
                    continue
                _, packet_types, _, _ = packets(core, a, w, clock)
                assert packet_types > 0
                missing_clock_rows += 1

    larger_clock_rows = 0
    examples = []
    for clock in range(15, 51):
        core = tuple(1 + j * clock for j in range(12))
        for a in range(1, 6):
            for w in range(1, 201):
                if w in {a * c for c in core}:
                    continue
                _, packet_types, alpha, beta = packets(core, a, w, clock)
                if packet_types:
                    larger_clock_rows += 1
                    if len(examples) < 5:
                        examples.append(
                            (clock, a, w, packet_types, alpha, beta)
                        )

    print("THM-2059 CRT PHASE-PACKET AUDIT")
    print("direct identity rows checked:", rows)
    print("direct grid indices checked:", grid_indices)
    print("positive/zero packet rows:", positive, zero)
    print("missing-clock specialization rows checked:", missing_clock_rows)
    print("larger-clock certificates found:", larger_clock_rows)
    print("first larger-clock examples (N,a,w,P,alpha,beta):")
    for example in examples:
        print(" ", example)
    print("carrier: bipartite CRT compatibility graph; tournament orientation rejected")
    print("destroyed by scalar packet sizes: reduction-class overlap")
    print("PASS")


if __name__ == "__main__":
    main()
