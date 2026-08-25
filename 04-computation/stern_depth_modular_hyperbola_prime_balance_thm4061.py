#!/usr/bin/env python3
"""Exact audit for the modular-hyperbola form of Stern depth packets.

The characteristic-zero analytic bound uses Weil's prime Kloosterman bound;
this companion audits every elementary reduction and its exact consequence.
"""

from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt


QMAX = 5000


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def is_prime(n):
    if n < 2:
        return False
    return all(n % d for d in range(2, isqrt(n) + 1))


def phi(n):
    return sum(gcd(a, n) == 1 for a in range(1, n))


def stern_depth_below_one(a, q):
    require(0 < a < q and gcd(a, q) == 1, (a, q))
    x, y = q, a
    digit_sum = 0
    while y:
        digit_sum += x // y
        x, y = y, x % y
    return digit_sum - 1


def packet_direct(q):
    return sum(
        -1 if stern_depth_below_one(a, q) & 1 else 1
        for a in range(1, q)
        if gcd(a, q) == 1
    )


def packet_inverse_parity(q):
    return sum(
        -1 if (a + pow(a, -1, q)) & 1 else 1
        for a in range(1, q)
        if gcd(a, q) == 1
    )


def quadrant_count(q):
    m = (q - 1) // 2
    return sum(
        (4 * r * s - 1) % q == 0
        for r in range(1, m + 1)
        for s in range(1, m + 1)
    )


def quadrant_count_inverse(q):
    m = (q - 1) // 2
    inv4 = pow(4, -1, q)
    return sum(
        1 <= (inv4 * pow(r, -1, q)) % q <= m
        for r in range(1, m + 1)
        if gcd(r, q) == 1
    )


def harmonic(n):
    return sum((Fraction(1, k) for k in range(1, n + 1)), Fraction())


def harmonic_cap(p):
    return 2 * harmonic(p - 1) - harmonic((p - 1) // 2)


def prime_bound_holds(p, packet):
    hp = harmonic_cap(p)
    rational_part = Fraction(p - 1, p * p) + Fraction(2, p) * hp
    radical_coefficient = 2 * hp * hp
    excess = Fraction(abs(packet)) - rational_part
    if excess <= 0:
        return True
    return excess * excess <= radical_coefficient * radical_coefficient * p


records = []
prime_records = []
zero_packets = []
for q in range(3, QMAX + 1, 2):
    packet_depth = packet_direct(q)
    packet_inverse = packet_inverse_parity(q)
    n_box = quadrant_count_inverse(q)
    if q <= 301:
        require(n_box == quadrant_count(q), ("box implementations", q))
    totient = phi(q)
    require(
        packet_depth == packet_inverse,
        ("depth/inverse", q, packet_depth, packet_inverse),
    )
    require(
        packet_depth == 4 * n_box - totient,
        ("quadrant", q, packet_depth, n_box, totient),
    )
    require(
        (totient // 2)
        == sum(gcd(a, q) == 1 and a % 2 == 0 for a in range(1, q)),
        ("half shell", q),
    )
    records.append((q, packet_depth, n_box, totient))
    if packet_depth == 0:
        zero_packets.append(q)
    if is_prime(q):
        require(prime_bound_holds(q, packet_depth), ("prime bound", q, packet_depth))
        prime_records.append((q, packet_depth, Fraction(packet_depth * packet_depth, q)))

top = sorted(prime_records, key=lambda row: row[2], reverse=True)[:12]
require(
    all(q % 4 == 1 for q in zero_packets if is_prime(q)),
    "prime zero congruence",
)
require(
    [(q, packet) for q, packet, _ in prime_records if q in (5, 13)]
    == [(5, 0), (13, 0)],
    "Paley controls",
)

serialization = ";".join(
    f"{q}:{packet}:{n_box}:{totient}"
    for q, packet, n_box, totient in records
)
print("status=PASS")
print("scope=odd_q_exact_through_5000;prime_weil_consequence_checked_exactly")
print("identity=S(q)=4*N_q-phi(q)")
print("N_q=#{1<=r,s<=(q-1)/2:4rs=1_mod_q}")
print("prime_bound=|S(p)|<=2*sqrt(p)*Hstar(p)^2+2*Hstar(p)/p+(p-1)/p^2")
print("Hstar(p)=2*H_(p-1)-H_((p-1)/2)")
print("odd_denominators_checked=" + str(len(records)))
print("odd_primes_checked=" + str(len(prime_records)))
print("zero_packets_count=" + str(len(zero_packets)))
print("zero_packets_first=" + repr(zero_packets[:40]))
print("largest_observed_S2_over_p=" + repr(top))
print("table_sha256=" + sha256(serialization.encode()).hexdigest())
print("boundary=Weil_bound_is_cited_not_reproved;no_zero_classification_or_LRC_claim")
