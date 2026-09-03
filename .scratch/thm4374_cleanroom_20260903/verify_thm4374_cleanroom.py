#!/usr/bin/env python3
"""Independent exact verifier for the proposed THM-4374.

This file has no repository imports.  It uses a direct centered-rotation
model, exact integer rational arithmetic, structural fibre representatives,
and an explicit inverse decoder for W_17.
"""

import sys
from math import gcd


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")


A = 3371
S = 1303
M = 14 * A
N = M // 2
TAIL = 11019
ACTIVE_MAX_C = 2 * A - 2
UNITS_14 = (1, 3, 5, 9, 11, 13)

checks = 0


def check(condition, label):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError("check failed: " + label)


def rat(num, den=1):
    if den == 0:
        raise ZeroDivisionError("zero rational denominator")
    if den < 0:
        num = -num
        den = -den
    d = gcd(abs(num), den)
    return (num // d, den // d)


def radd(x, y):
    return rat(x[0] * y[1] + y[0] * x[1], x[1] * y[1])


def rsub(x, y):
    return rat(x[0] * y[1] - y[0] * x[1], x[1] * y[1])


def rscale(x, k):
    return rat(x[0] * k, x[1])


def rdiv(x, y):
    return rat(x[0] * y[1], x[1] * y[0])


def rgt(x, y):
    return x[0] * y[1] > y[0] * x[1]


def rtext(x):
    return f"{x[0]}/{x[1]}"


BASE = rat(S, M)


def centered_rho(p):
    r = (S * p) % M
    if r > N:
        r -= M
    return r


def centered_c(p):
    return A - centered_rho(p)


def phase_x(p):
    c = centered_c(p)
    check(c % 2 == 0, "odd parameter has even c")
    return (c // 2) % N


def c_is_active(c):
    return 0 < c < 2 * A


def phase_is_active(x):
    return 1 <= x <= A - 1


def exit_value(p):
    check(p > 0 and p % 2 == 1, "exit parameter is positive odd")
    c = centered_c(p)
    if c_is_active(c):
        return rat(S * p + c, M * p)
    return BASE


def output_is_active(value):
    return rgt(value, BASE)


def word(p, horizon=17):
    return tuple(exit_value(p + 2 * t) for t in range(horizon + 1))


def only_integer(value, label):
    check(value[1] == 1, label + " is integral")
    return value[0]


def decode_consecutive(q_word, r):
    """Decode Q=P+2r from active outputs at r and r+1, before a wrap."""
    d0 = rsub(q_word[r], BASE)
    d1 = rsub(q_word[r + 1], BASE)
    top = rscale(radd(rat(S), rscale(d1, M)), 2)
    bottom = rscale(rsub(d0, d1), M)
    return only_integer(rdiv(top, bottom), "consecutive decoder")


def decode_one_wrap(q_word, t):
    """Decode P from active outputs at 0 and t after exactly one wrap."""
    d0 = rsub(q_word[0], BASE)
    dt = rsub(q_word[t], BASE)
    top = radd(rat(M - 2 * t * S), rscale(dt, -2 * t * M))
    bottom = rscale(rsub(dt, d0), M)
    return only_integer(rdiv(top, bottom), "one-wrap decoder")


def decode_w17(q_word):
    check(len(q_word) == 18, "decoder input length")
    active_times = [t for t, value in enumerate(q_word) if output_is_active(value)]
    check(active_times and active_times[0] <= 16, "activity by time 16")
    r = active_times[0]
    if r > 0:
        check(output_is_active(q_word[r + 1]), "delayed entry has active successor")
        q = decode_consecutive(q_word, r)
        return q - 2 * r
    if output_is_active(q_word[1]):
        return decode_consecutive(q_word, 0)
    t = 16 if output_is_active(q_word[16]) else 17
    check(output_is_active(q_word[t]), "low current phase returns by 17")
    return decode_one_wrap(q_word, t)


def realize_b(a, unit):
    """Realize a structural (a,kappa mod 14) type with a full tail fibre."""
    inv_s = pow(S, -1, M)
    b = ((A * unit - a) * inv_s) % M
    check(b % 2 == 1, "CRT representative b is odd")
    while b < TAIL or gcd(a, b) != 1:
        b += M
    kappa_num = a + S * b
    check(kappa_num % A == 0, "a+Sb divisible by A")
    kappa = kappa_num // A
    check(kappa > 0, "positive kappa")
    check(kappa % 14 == unit, "prescribed kappa residue")
    check(gcd(kappa, 14) == 1, "kappa is a unit modulo 14")
    check(gcd(a, b) == 1 and b % 2 == 1, "primitive odd b")
    return b, kappa


def structural_scales(a, unit):
    g0 = pow(unit, -1, 14)
    return tuple(range(g0, ACTIVE_MAX_C // a + 1, 14))


def expected_partition_key(a, g, horizon):
    c = a * g
    if horizon == 0:
        return ("all",)
    if horizon <= 15:
        return ("low",) if c <= 2606 else ("scale", g)
    if horizon == 16:
        return ("mid",) if 1244 <= c <= 2606 else ("scale", g)
    return ("scale", g)


def normalized_groups(groups):
    return tuple(sorted((tuple(values) for values in groups.values())))


def partition_sizes_for_words(gs, words, horizon):
    groups = {}
    for g, values in zip(gs, words):
        groups.setdefault(values[: horizon + 1], []).append(g)
    return sorted((len(values) for values in groups.values()), reverse=True)


def audit_constants_and_boundaries():
    check(A == 3371 and S == 1303, "constants")
    check(M == 47194 and N == 23597 and N == 7 * A, "moduli")
    check(gcd(A, S) == 1 and gcd(S, M) == 1, "coprime constants")
    check(pow(S, -1, M) * S % M == 1, "inverse of S modulo M")
    check(pow(1303, -1, 47194) == 1485, "published inverse")

    boundary_rows = []
    for target_rho in (-A, A):
        found = None
        for p in range(TAIL, TAIL + M, 2):
            if centered_rho(p) == target_rho:
                found = p
                break
        check(found is not None, "boundary phase exists")
        check(exit_value(found) == BASE, "boundary is inactive baseline")
        boundary_rows.append((found, target_rho))
    check(tuple(boundary_rows) == ((43823, -3371), (50565, 3371)),
          "least cofinite boundary representatives")
    return tuple(boundary_rows)


def audit_phase_cover():
    residue_to_phase = {}
    for p in range(1, M, 2):
        x = phase_x(p)
        check(x not in residue_to_phase, "odd residues map injectively to phase")
        residue_to_phase[x] = p
    check(len(residue_to_phase) == N, "odd residues cover every half-phase")

    covered = {
        (y + t * S) % N
        for t in range(17)
        for y in range(1, A)
    }
    first_hit_counts = [0] * 17
    for x0 in range(N):
        hits = []
        for t in range(17):
            x = (x0 - t * S) % N
            if phase_is_active(x):
                hits.append(t)
        check(hits, "every phase hits by time 16")
        r = hits[0]
        first_hit_counts[r] += 1
        if r > 0:
            y = (x0 - r * S) % N
            check(2068 <= y <= 3370, "delayed entry strip")
            check(phase_is_active((y - S) % N), "delayed entry successor active")

        p = residue_to_phase[x0]
        while p < TAIL:
            p += M
        check(phase_x(p) == x0, "tail representative preserves phase")
        for t in range(18):
            x = (x0 - t * S) % N
            check(output_is_active(exit_value(p + 2 * t)) == phase_is_active(x),
                  "output detects strict activity")

    check(len(covered) == N, "seventeen preimage intervals cover all phases")
    check(3370 + 16 * S == 24218 and 24218 >= N,
          "unwrapped interval cover endpoint")
    expected = [3370] + [1303] * 15 + [682]
    check(first_hit_counts == expected, "exact first-hit census")
    return tuple(first_hit_counts)


def audit_all_active_fibres():
    structural_types = 0
    scale_points = 0
    pair_total = 0
    split_counts = {1: 0, 16: 0, 17: 0}
    max_sizes = [0] * 18
    max_witnesses = [None] * 18

    for unit in UNITS_14:
        check((16 * unit) % 7 != 0 and (17 * unit) % 7 != 0,
              "return collision obstruction modulo 7")

    for a in range(2, ACTIVE_MAX_C + 1, 2):
        for unit in UNITS_14:
            gs = structural_scales(a, unit)
            if not gs:
                continue
            structural_types += 1
            b, kappa = realize_b(a, unit)
            values_by_g = []

            for g in gs:
                p = b * g
                c = a * g
                check(p >= TAIL and p % 2 == 1, "realized scale lies in odd tail")
                check(g * kappa % 14 == 1, "scale congruence")
                check(centered_c(p) == c, "realized active centered coordinate")
                values = word(p, 17)
                check(values[0] == rat(kappa, 14 * b), "current fibre exit formula")

                if c >= 2608:
                    check(output_is_active(values[1]), "high scale active at time one")
                    check(values[1] == rat(g * kappa, 14 * (b * g + 2)),
                          "time-one return formula")
                else:
                    for t in range(1, 16):
                        check(values[t] == BASE, "low scale inactive through time 15")
                    if c <= 1242:
                        check(centered_c(p + 32) == c + 5498,
                              "time-16 wrapped c formula")
                        check(values[16] == rat(g * kappa + 14,
                                               14 * (b * g + 32)),
                              "time-16 return output formula")
                    else:
                        check(centered_c(p + 32) == c + 5498,
                              "time-16 boundary/mid c formula")
                        check(values[16] == BASE, "middle scale inactive at time 16")
                    check(centered_c(p + 34) == c + 2892,
                          "time-17 wrapped c formula")
                    check(values[17] == rat(g * kappa + 14,
                                           14 * (b * g + 34)),
                          "time-17 return output formula")

                values_by_g.append(values)
                scale_points += 1

            for horizon in range(18):
                actual = {}
                expected = {}
                for g, values in zip(gs, values_by_g):
                    actual.setdefault(values[: horizon + 1], []).append(g)
                    key = expected_partition_key(a, g, horizon)
                    expected.setdefault(key, []).append(g)
                check(normalized_groups(actual) == normalized_groups(expected),
                      "exact W_H partition")
                largest = max(len(group) for group in actual.values())
                if largest > max_sizes[horizon]:
                    max_sizes[horizon] = largest
                    max_witnesses[horizon] = (a, b, kappa, unit)

            for i in range(len(gs)):
                g1 = gs[i]
                c1 = a * g1
                w1 = values_by_g[i]
                for j in range(i + 1, len(gs)):
                    g2 = gs[j]
                    c2 = a * g2
                    w2 = values_by_g[j]
                    first = None
                    for t in range(1, 18):
                        if w1[t] != w2[t]:
                            first = t
                            break
                    check(first in split_counts, "distinct scales split by horizon 17")
                    if c2 >= 2608:
                        expected_first = 1
                    elif c1 <= 1242:
                        expected_first = 16
                    else:
                        expected_first = 17
                    check(first == expected_first, "first-split trichotomy")

                    if c1 >= 2608:
                        lhs = g1 * kappa * (b * g2 + 2) - g2 * kappa * (b * g1 + 2)
                        check(lhs == 2 * kappa * (g1 - g2),
                              "time-one cross-product factorization")
                    elif first in (16, 17) and c2 <= 2606:
                        t = first
                        lhs = ((g1 * kappa + 14) * (b * g2 + 2 * t)
                               - (g2 * kappa + 14) * (b * g1 + 2 * t))
                        check(lhs == 2 * (g1 - g2) * (t * kappa - 7 * b),
                              "return cross-product factorization")
                        check((t * kappa - 7 * b) % 7 != 0,
                              "return cross-product nonzero modulo 7")

                    split_counts[first] += 1
                    pair_total += 1

    expected_max = [241] + [94] * 15 + [49, 1]
    check(max_sizes == expected_max, "sharp horizon maxima")
    check(pair_total == 281073, "complete structural pair count")
    check(split_counts == {1: 239955, 16: 30224, 17: 10894},
          "exact first-split census")
    check(scale_points == 13427, "complete structural scale-point count")
    return (structural_types, scale_points, pair_total, split_counts,
            tuple(max_sizes), tuple(max_witnesses))


def audit_sharp_witness():
    a, b, kappa = 2, 47595, 18397
    check(a + S * b == A * kappa, "sharp witness Diophantine identity")
    check(kappa % 14 == 1, "sharp witness residue")
    gs = structural_scales(a, kappa % 14)
    check(gs == tuple(1 + 14 * j for j in range(241)), "sharp witness full G")
    values = [word(b * g, 17) for g in gs]
    maxima = tuple(partition_sizes_for_words(gs, values, h)[0] for h in range(18))
    check(maxima == tuple([241] + [94] * 15 + [49, 1]),
          "one witness realizes every sharp maximum")

    low = tuple(g for g in gs if a * g <= 2606)
    mid = tuple(g for g in gs if 1244 <= a * g <= 2606)
    check(len(low) == 94 and low[0] == 1 and low[-1] == 1303,
          "sharp 94-scale bucket")
    check(len(mid) == 49 and mid[0] == 631 and mid[-1] == 1303,
          "sharp 49-scale bucket")
    return maxima, (gs[0], gs[-1]), (low[0], low[-1]), (mid[0], mid[-1])


def audit_global_decoder_and_minimality():
    literal_rows = 2 * N
    branch_counts = {"delayed": 0, "current-next": 0,
                     "current-return16": 0, "current-return17": 0}
    for j in range(literal_rows):
        p = TAIL + 2 * j
        values = word(p, 17)
        decoded = decode_w17(values)
        check(decoded == p, "direct W_17 inverse on two phase periods")
        r = next(t for t, value in enumerate(values) if output_is_active(value))
        if r > 0:
            branch_counts["delayed"] += 1
        elif output_is_active(values[1]):
            branch_counts["current-next"] += 1
        elif output_is_active(values[16]):
            branch_counts["current-return16"] += 1
        else:
            branch_counts["current-return17"] += 1

    p1, p2 = 253031, 258645
    w1 = word(p1, 17)
    w2 = word(p2, 17)
    check(w1[:17] == w2[:17], "global W_16 hostile equality")
    check(w1[17] != w2[17], "global W_17 hostile split")
    check(w1[0] == rat(155, 5614), "hostile current exit")
    for t in range(1, 17):
        check(w1[t] == BASE and w2[t] == BASE,
              "hostile baseline through time 16")
    check(w1[17] == rat(97819, 3542910), "first hostile return fraction")
    check(w2[17] == rat(99989, 3621506), "second hostile return fraction")
    check(decode_w17(w1) == p1 and decode_w17(w2) == p2,
          "decoder separates minimality hostile")

    inactive_seed = None
    for p in range(TAIL, TAIL + M, 2):
        if exit_value(p) == BASE:
            inactive_seed = p
            break
    check(inactive_seed is not None, "inactive baseline seed")
    inactive_ray = tuple(inactive_seed + j * M for j in range(6))
    check(all(exit_value(p) == BASE for p in inactive_ray),
          "infinite inactive residue ray positive control")
    return literal_rows, branch_counts, (p1, p2), w1[17], w2[17], inactive_ray


def main():
    boundaries = audit_constants_and_boundaries()
    first_hits = audit_phase_cover()
    (structural_types, scale_points, pair_total, split_counts,
     max_sizes, max_witnesses) = audit_all_active_fibres()
    witness_maxima, full_span, low_span, mid_span = audit_sharp_witness()
    (literal_rows, branch_counts, hostile, hostile_exit_1,
     hostile_exit_2, inactive_ray) = audit_global_decoder_and_minimality()

    print("STATUS PASS")
    print(f"CHECKS {checks}")
    print(f"CONSTANTS A={A} S={S} M={M} N={N} TAIL={TAIL}")
    print("PHASE_COVER " + ",".join(str(v) for v in first_hits))
    print("BOUNDARIES " + " ".join(f"P={p}:rho={rho}" for p, rho in boundaries))
    print(f"STRUCTURAL_TYPES {structural_types}")
    print(f"STRUCTURAL_SCALE_POINTS {scale_points}")
    print(f"ACTIVE_FIBRE_SCALE_PAIRS {pair_total}")
    print("FIRST_SPLIT " + " ".join(f"{t}:{split_counts[t]}" for t in (1, 16, 17)))
    print("MAX_HORIZON_FIBRE " + " ".join(
        f"{h}:{max_sizes[h]}" for h in range(18)))
    print("MAX_WITNESS_FIRST " + " ".join(
        f"{h}:{max_witnesses[h]}" for h in (0, 1, 16, 17)))
    print("SHARP_WITNESS a=2 b=47595 kappa=18397 "
          f"full={witness_maxima[0]}:{full_span[0]}..{full_span[1]} "
          f"low={witness_maxima[1]}:{low_span[0]}..{low_span[1]} "
          f"mid={witness_maxima[16]}:{mid_span[0]}..{mid_span[1]} "
          f"h17={witness_maxima[17]}")
    print(f"GLOBAL_DECODER_ROWS {literal_rows}")
    print("GLOBAL_DECODER_BRANCHES " + " ".join(
        f"{key}:{branch_counts[key]}" for key in
        ("delayed", "current-next", "current-return16", "current-return17")))
    print(f"MINIMALITY_HOSTILE P={hostile[0]},{hostile[1]} "
          f"E17={rtext(hostile_exit_1)},{rtext(hostile_exit_2)}")
    print("INACTIVE_RAY " + ",".join(str(p) for p in inactive_ray))
    print("CONGRUENCE_COROLLARY exact: W_17 inverse forces equality")
    print("SCOPE fixed THM-4365 h=420 odd tail; metric exit only; LRC(14) OPEN")


if __name__ == "__main__":
    main()
