#!/usr/bin/env python3
"""Independent Farey-tree audit of the modular-hyperbola packet identity."""

from fractions import Fraction
from hashlib import sha256
from math import gcd, isqrt


QMAX = 3000


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def is_prime(n):
    if n < 2:
        return False
    return all(n % d for d in range(2, isqrt(n) + 1))


def harmonic(n):
    return sum((Fraction(1, k) for k in range(1, n + 1)), Fraction())


def farey_tree_packet_sums(limit):
    sums = [0] * (limit + 1)
    counts = [0] * (limit + 1)
    # An interval (a/b,c/d) contributes its mediant at the supplied depth.
    stack = [(0, 1, 1, 1, 1)]
    while stack:
        a, b, c, d, depth = stack.pop()
        p, q = a + c, b + d
        if q > limit:
            continue
        sums[q] += -1 if depth & 1 else 1
        counts[q] += 1
        stack.append((a, b, p, q, depth + 1))
        stack.append((p, q, c, d, depth + 1))
    return sums, counts


def hyperbola_count(q):
    m = (q - 1) // 2
    inv4 = pow(4, -1, q)
    total = 0
    for r in range(1, m + 1):
        if gcd(r, q) == 1:
            s = (inv4 * pow(r, -1, q)) % q
            total += s <= m
    return total


def bound_check(p, packet):
    hstar = 2 * harmonic(p - 1) - harmonic((p - 1) // 2)
    rational_part = Fraction(p - 1, p * p) + Fraction(2, p) * hstar
    radical_coefficient = 2 * hstar * hstar
    excess = Fraction(abs(packet)) - rational_part
    return (
        excess <= 0
        or excess * excess <= radical_coefficient * radical_coefficient * p
    )


sums, counts = farey_tree_packet_sums(QMAX)
rows = []
prime_rows = []
for q in range(3, QMAX + 1, 2):
    totient = sum(gcd(a, q) == 1 for a in range(1, q))
    require(counts[q] == totient, ("Farey completeness", q, counts[q], totient))
    n_box = hyperbola_count(q)
    require(
        sums[q] == 4 * n_box - totient,
        ("independent quadrant", q, sums[q], n_box, totient),
    )
    rows.append((q, sums[q], n_box, totient))
    if is_prime(q):
        require(bound_check(q, sums[q]), ("prime bound", q, sums[q]))
        prime_rows.append((q, sums[q]))

require(sums[5] == 0 and sums[13] == 0, "Paley zero controls")
require(sums[9] == -2 and sums[11] == 2, "first noncharacter controls")
serialization = ";".join(
    f"{q}:{packet}:{n_box}:{totient}"
    for q, packet, n_box, totient in rows
)
print("status=PASS")
print("implementation=iterative_Farey_tree_no_continued_fraction_or_packet_import")
print("fractions_generated=" + str(sum(counts)))
print("odd_denominators_checked=" + str(len(rows)))
print("odd_primes_checked=" + str(len(prime_rows)))
print("controls=p5_zero;p13_zero;q9_noncharacter;q11_hostile")
print("table_sha256=" + sha256(serialization.encode()).hexdigest())
print("boundary=finite_exact_audit;Weil_prime_bound_remains_cited")
