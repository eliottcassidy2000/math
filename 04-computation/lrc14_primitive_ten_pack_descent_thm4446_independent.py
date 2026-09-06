#!/usr/bin/env python3
"""Independent exact checks for THM-4446.

This file imports neither the primary script nor any repository geometry
engine. It audits translated open-arc grid caps, literal 3p-label unions for
small hostile primes, and the dilation-ray normalization arithmetic.
"""

from fractions import Fraction as Q
from functools import reduce
from itertools import combinations
from math import gcd

DELTA = Q(1, 14)
CHECKS = 0


def need(test, detail):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(detail)


def floor_q(x):
    return x.numerator // x.denominator


def frac(x):
    return x - floor_q(x)


def norm(x):
    y = frac(x)
    return min(y, 1 - y)


def danger(x):
    return norm(x) < DELTA


def translated_grid_max(order):
    """Maximum points of a translated order-grid in the open danger arc."""
    events = {
        frac(sign * DELTA - Q(j, order))
        for sign in (-1, 1)
        for j in range(order)
    }
    ordered = sorted(events)
    probes = set(ordered)
    for i, left in enumerate(ordered):
        right = ordered[(i + 1) % len(ordered)]
        if i + 1 == len(ordered):
            right += 1
        probes.add(frac((left + right) / 2))
    return max(
        sum(danger(alpha + Q(j, order)) for j in range(order))
        for alpha in probes
    )


def bad_labels(p, t, y):
    return {j for j in range(3 * p) if danger(Q(t, 3 * p) * (y + j))}


def tail_event_probes(p, tails):
    """One exact point from every cell of the y-dependent arrangement."""
    events = {Q(0), Q(1)}
    for t in tails:
        for j in range(3 * p):
            lo = Q(t * j, 3 * p)
            hi = Q(t * (j + 1), 3 * p)
            for n in range(floor_q(lo) - 1, floor_q(hi) + 2):
                for sign in (-1, 1):
                    boundary = Q(n) + sign * DELTA
                    if lo <= boundary <= hi:
                        y = Q(3 * p, t) * boundary - j
                        if 0 <= y <= 1:
                            events.add(y)
    ordered = sorted(events)
    probes = set(ordered)
    for left, right in zip(ordered, ordered[1:]):
        probes.add((left + right) / 2)
    return sorted(probes)


def is_prime(n):
    if n < 2:
        return False
    p = 2
    while p * p <= n:
        if n % p == 0:
            return n == p
        p += 1
    return True


# Open-arc cap, including exact divisibility by seven.
for m in range(1, 211):
    need(
        translated_grid_max(m) == (m + 6) // 7,
        ("translated-grid cap", m, translated_grid_max(m), (m + 6) // 7),
    )

# Prime arithmetic and every admissible p-divisibility word.
prime_table = {}
for p in [n for n in range(2, 1001) if is_prime(n)]:
    B = (3 * p + 6) // 7
    need(B < p, ("B_p<p", p, B))
    words = [
        (d1, d2, d3)
        for d1 in (0, 1)
        for d2 in (0, 1)
        for d3 in (0, 1)
        if d1 + d2 + d3 <= 2 and (p != 3 or d1 + d2 + d3 == 0)
    ]
    margins = []
    for word in words:
        bad = sum(p if bit else B for bit in word)
        margin = 3 * p - bad
        need(margin > 0, ("word margin", p, word, bad, margin))
        margins.append(margin)
    prime_table[p] = (B, min(margins))

need(prime_table[2] == (1, 1), ("p=2 boundary", prime_table[2]))
need(prime_table[3] == (2, 3), ("p=3 boundary", prime_table[3]))
need(prime_table[7] == (3, 4), ("p=7 boundary", prime_table[7]))

# Exhaust all admissible distinct ternary-unit triples through 20 at the
# delicate small primes and every cell of their y-dependent arrangements.
literal_records = {}
tails = [t for t in range(1, 21) if t % 3]
for p in (2, 3, 5, 7):
    minimum = None
    witnesses = []
    rows = 0
    for T in combinations(tails, 3):
        divisible = sum(t % p == 0 for t in T)
        if divisible == 3 or (p == 3 and divisible):
            continue
        predicted = 3 * p - (
            divisible * p + (3 - divisible) * prime_table[p][0]
        )
        need(predicted > 0, ("predicted literal margin", p, T, predicted))
        for y in tail_event_probes(p, T):
            bad = set().union(*(bad_labels(p, t, y) for t in T))
            safe = 3 * p - len(bad)
            need(
                safe >= predicted,
                ("literal union failure", p, T, y, safe, predicted, bad),
            )
            if minimum is None or safe < minimum:
                minimum = safe
                witnesses = [(T, y)]
            elif safe == minimum and len(witnesses) < 4:
                witnesses.append((T, y))
            rows += 1
    literal_records[p] = (minimum, rows, witnesses)

# The p=2 margin is attained by a full primitive row, not only by abstract
# bad subsets.
sharp_Cprime = (1, 2, 3, 4, 6, 7, 8, 9, 10, 11)
sharp_T = (1, 4, 10)
sharp_y = Q(23, 56)
need(
    all(norm(c * sharp_y) >= DELTA for c in sharp_Cprime),
    ("sharp p2 body phase", sharp_Cprime, sharp_y),
)
sharp_safe_labels = [
    j
    for j in range(6)
    if all(not danger(Q(t, 6) * (sharp_y + j)) for t in sharp_T)
]
need(
    sharp_safe_labels == [3],
    ("sharp p2 one-label witness", sharp_T, sharp_y, sharp_safe_labels),
)
sharp_row = tuple(6 * c for c in sharp_Cprime) + sharp_T
need(
    len(set(sharp_row)) == 13 and reduce(gcd, sharp_row) == 1,
    ("sharp p2 full row type", sharp_row),
)

# Every ten-subset of [13] is gcd-one.
bodies = list(combinations(range(1, 14), 10))
need(len(bodies) == 286, ("body count", len(bodies)))
for C0 in bodies:
    need(reduce(gcd, C0) == 1, ("bounded body gcd", C0))

# Exhaust the normalization identities g|h and gcd of the divided row.
normalization_rows = 0
tail_pool = [t for t in range(1, 41) if t % 3]
for h in range(1, 101):
    for T in combinations(tail_pool, 3):
        g = reduce(gcd, (3 * h, *T))
        need(g % 3 != 0, ("normalized gcd ternary", h, T, g))
        need(h % g == 0, ("g divides h", h, T, g))
        k = h // g
        Tg = tuple(t // g for t in T)
        need(
            all(t % 3 for t in Tg) and len(set(Tg)) == 3,
            ("divided tail type", h, T, g, Tg),
        )
        need(
            reduce(gcd, (3 * k, *Tg)) == 1,
            ("divided full row primitive", h, T, g, k, Tg),
        )
        for C0 in (bodies[0], bodies[97], bodies[-1]):
            need(
                reduce(gcd, (x * k for x in C0)) == k,
                ("normalized body gcd", C0, h, T, g, k),
            )
        normalization_rows += 1

print("TEN_BODY_ENTRY_BRIDGE_CLEANROOM_CHECKS")
print(
    "open_grid_caps_m_1_to_210",
    "PASS",
    "m_3_6_7_9_21",
    {m: translated_grid_max(m) for m in (3, 6, 7, 9, 21)},
)
print(
    "prime_margins_p_2_3_5_7", {p: prime_table[p] for p in (2, 3, 5, 7)}
)
print("literal_small_prime_records", literal_records)
print(
    "sharp_p2_full_row",
    "Cprime",
    sharp_Cprime,
    "T",
    sharp_T,
    "y",
    sharp_y,
    "safe_labels",
    sharp_safe_labels,
)
print("bounded_body_gcd_one_rows", len(bodies))
print("normalization_rows", normalization_rows)
print("explicit_checks", CHECKS)
print("PASS")
