#!/usr/bin/env python3
"""Exact hostile audit for THM-4059.

The truth-bearing paths use explicit RuntimeError checks, not assert, so
``python -O`` executes the same audit.  All arithmetic is integral/rational.
"""

from collections import defaultdict, deque
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import gcd, isqrt


GATES = 0


def check(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"exact check failed: {label}")


def canonical_cf(p: int, q: int) -> list[int]:
    out: list[int] = []
    while q:
        out.append(p // q)
        p, q = q, p % q
    if len(out) > 1 and out[-1] == 1:
        out[-2] += 1
        out.pop()
    return out


def sb_depth(p: int, q: int) -> int:
    return sum(canonical_cf(p, q)) - 1


def epsilon(p: int, q: int) -> int:
    return -1 if sb_depth(p, q) & 1 else 1


def inverse_data(a: int, q: int) -> tuple[int, int]:
    u = pow(a, -1, q)
    k = (a * u - 1) // q
    return u, k


def inverse_formula(a: int, q: int) -> int:
    u, k = inverse_data(a, q)
    if q & 1:
        return -1 if (a + u) & 1 else 1
    return -1 if (k + 1) & 1 else 1


def xor_vector(v: tuple[int, int], w: tuple[int, int]) -> tuple[int, int]:
    return (v[0] ^ w[0], v[1] ^ w[1])


def permutation_sign(indices: list[int]) -> int:
    inversions = sum(
        indices[i] > indices[j]
        for i in range(len(indices))
        for j in range(i + 1, len(indices))
    )
    return -1 if inversions & 1 else 1


def flank_hexagon_character(a: int, q: int) -> int:
    """Relative S3 sign of the canonical parent flank against root J."""
    u, k = inverse_data(a, q)
    c1 = (k & 1, u & 1)
    c2 = ((a - k) & 1, (q - u) & 1)
    triple = (c1, xor_vector(c1, c2), c2)
    root_order = ((0, 1), (1, 1), (1, 0))
    indices = [root_order.index(v) for v in triple]
    return permutation_sign(indices)


def divisors(n: int) -> list[int]:
    small: list[int] = []
    large: list[int] = []
    for d in range(1, isqrt(n) + 1):
        if n % d:
            continue
        small.append(d)
        if d * d != n:
            large.append(n // d)
    return small + list(reversed(large))


def mobius(n: int) -> int:
    if n == 1:
        return 1
    sign = 1
    p = 2
    while p * p <= n:
        if n % p:
            p += 1
            continue
        n //= p
        if n % p == 0:
            return 0
        sign = -sign
        while n % p == 0:
            n //= p
        p += 1
    if n > 1:
        sign = -sign
    return sign


def phi(n: int) -> int:
    result = n
    p = 2
    while p * p <= n:
        if n % p:
            p += 1
            continue
        result -= result // p
        while n % p == 0:
            n //= p
        p += 1
    if n > 1:
        result -= result // n
    return result


def legendre(a: int, p: int) -> int:
    residue = pow(a % p, (p - 1) // 2, p)
    return 1 if residue == 1 else -1


def exceptional_mod_four(q: int) -> bool:
    """Whether the proved q>2 classification predicts S(q)=2 mod 4."""
    if q == 4:
        return True
    if q % 2 == 0:
        q //= 2
    if q % 2 == 0 or q == 1:
        return False
    prime = None
    n = q
    p = 3
    while p * p <= n:
        if n % p:
            p += 2
            continue
        if prime is not None:
            return False
        prime = p
        while n % p == 0:
            n //= p
        p += 2
    if n > 1:
        if prime is not None:
            return False
        prime = n
    return prime is not None and prime % 4 == 3


print("THM-4059 Stern--Brocot depth-packet character audit")
print()

# 1. Parent-flank, mod-two hexagon, and modular-inverse formula.
full_bound = 1200
primitive_cases = 0
for q in range(2, full_bound + 1):
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        u, k = inverse_data(a, q)
        check(1 <= u < q, f"inverse range {a}/{q}")
        check(0 <= k < a, f"flank numerator range {a}/{q}")
        check(k * (q - u) - u * (a - k) == -1,
              f"flank determinant {a}/{q}")
        check(k + (a - k) == a and u + (q - u) == q,
              f"flank mediant {a}/{q}")
        e = epsilon(a, q)
        check(flank_hexagon_character(a, q) == e,
              f"Farey hexagon character {a}/{q}")
        check(inverse_formula(a, q) == e,
              f"modular inverse formula {a}/{q}")
        check(sb_depth(u, q) == sb_depth(a, q),
              f"full-depth modular inversion {a}/{q}")
        check(inverse_formula(u, q) == e,
              f"reciprocal inverse formula {a}/{q}")
        primitive_cases += 1

print(f"canonical flank / inverse-formula cases through q={full_bound}: {primitive_cases}")
print("full-depth modular inversion and the odd/even inverse formulas: PASS")
print()

# 2. Packet sums, binary Farey recurrence, and mod-four classification.
packet_bound = 5000
packet_sum = [0] * (packet_bound + 1)
packet_sum[1] = 1  # Declared cyclic zero-phase weight; not D(0/1).
packet_units = 0
for q in range(2, packet_bound + 1):
    total = 0
    for a in range(1, q):
        if gcd(a, q) == 1:
            total += inverse_formula(a, q)
            packet_units += 1
    packet_sum[q] = total


@lru_cache(maxsize=None)
def packet_value(q: int) -> int:
    if q <= packet_bound:
        return packet_sum[q]
    return sum(
        inverse_formula(a, q)
        for a in range(1, q)
        if gcd(a, q) == 1
    )

check(packet_sum[1] == 1, "S(1) declared cyclic boundary")
check(packet_sum[2] == -1, "S(2) boundary")
recurrence_terms = 0
for q in range(3, packet_bound + 1):
    half_sum = 0
    for a in range(1, (q - 1) // 2 + 1):
        b = q - 2 * a
        if gcd(a, b) == 1:
            half_sum += epsilon(a, b)
            recurrence_terms += 1
    check(packet_sum[q] == 2 * half_sum, f"binary Farey recurrence q={q}")
    check((packet_sum[q] % 4 == 2) == exceptional_mod_four(q),
          f"mod-four classification q={q}")
    check((phi(q) % 4 == 2) == exceptional_mod_four(q),
          f"totient classification q={q}")

prefix = ",".join(str(packet_sum[q]) for q in range(2, packet_bound + 1))
prefix_hash = sha256(prefix.encode("ascii")).hexdigest()
check(prefix_hash == "9e5baaaff1df2977b70969f9108f2be47dd1313658b3c7e1ec484959941260ee",
      "S(2)..S(5000) prefix hash")

print(f"packet units through q={packet_bound}: {packet_units}")
print(f"binary-recurrence primitive terms: {recurrence_terms}")
print(f"S(2)..S({packet_bound}) sha256: {prefix_hash}")
print("S(q)=2 sum_{2a+b=q} epsilon(a/b) and exact mod-4 law: PASS")
print("S(2)..S(20):", packet_sum[2:21])
print()

# 3. Tournament lower stars, divisor convolution, and Möbius inversion.
star_bound = 1000
star_sum = [0] * (star_bound + 1)
for q in range(2, star_bound + 1):
    direct = 0
    for a in range(1, q):
        g = gcd(a, q)
        direct += epsilon(a // g, q // g)
    by_divisors = sum(packet_sum[d] for d in divisors(q) if d > 1)
    check(direct == by_divisors, f"lower-star divisor convolution q={q}")
    star_sum[q] = direct
    recovered = sum(mobius(q // d) * star_sum[d] for d in divisors(q))
    check(recovered == packet_sum[q], f"Mobius recovery q={q}")
    check((q - 1 + direct) % 2 == 0, f"star indegree parity q={q}")

print("B(q)=sum_{d|q,d>1}S(d); S(q)=sum_{d|q}mu(q/d)B(d) for q>1: PASS")
print("B(2)..B(20):", star_sum[2:21])
print()

# 4. Full weighted clock and signed Duffin--Schaeffer layer identities.
clock_rows = []
for n in (6, 60, 420, 840, 2520, 27720):
    direct_weight = Fraction(0)
    direct_ds = Fraction(0)
    direct_unweighted = 0
    for x in range(1, n):
        g = gcd(x, n)
        d = n // g
        a = x // g
        e = epsilon(a, d)
        weight = Fraction(d * d + 3, d + 1)
        psi = Fraction(d % 7 + 1, d + 3)
        direct_unweighted += e
        direct_weight += weight * e
        direct_ds += Fraction(2, d) * psi * e
    packet_weight = sum(
        Fraction(d * d + 3, d + 1) * packet_value(d)
        for d in divisors(n) if d > 1
    )
    packet_ds = sum(
        Fraction(2, d) * Fraction(d % 7 + 1, d + 3) * packet_value(d)
        for d in divisors(n) if d > 1
    )
    packet_unweighted = sum(packet_value(d) for d in divisors(n) if d > 1)
    check(direct_weight == packet_weight, f"arbitrary weighted clock N={n}")
    check(direct_ds == packet_ds, f"signed DS clock N={n}")
    check(direct_unweighted == packet_unweighted, f"unweighted clock N={n}")
    clock_rows.append((n, direct_unweighted))

print("arbitrary weighted and signed Duffin--Schaeffer clock identities: PASS")
print("unweighted master-clock characters:", clock_rows)
print()

# 5. Berggren/Pythagorean height columns and complete q<=1000 tree slice.
spinor_bound = 1000
expected_spinors: set[tuple[int, int]] = set()
expected_column = defaultdict(int)
for q in range(2, spinor_bound + 1):
    direct_column = 0
    for a in range(1, q):
        if gcd(a, q) == 1 and ((a ^ q) & 1):
            expected_spinors.add((a, q))
            e = epsilon(a, q)
            direct_column += e
            expected_column[q] += e
    formula = packet_sum[q] if q % 2 == 0 else packet_sum[q] // 2
    check(direct_column == formula, f"Berggren height column q={q}")

seen: set[tuple[int, int]] = set()
queue = deque([(1, 2, 0)])
while queue:
    m, n, a_parity = queue.popleft()
    if n > spinor_bound:
        continue
    check((m, n) not in seen, f"unique Berggren node {(m, n)}")
    seen.add((m, n))
    check(epsilon(m, n) == -(-1 if a_parity else 1),
          f"A-Walsh/depth sign {(m, n)}")
    children = (
        (n, 2 * n - m, a_parity ^ 1),
        (n, 2 * n + m, a_parity),
        (m, n + 2 * m, a_parity),
    )
    queue.extend(children)

check(seen == expected_spinors, "complete Berggren spinor slice q<=1000")
tree_column = defaultdict(int)
for m, n in seen:
    tree_column[n] += epsilon(m, n)
for q in range(2, spinor_bound + 1):
    check(tree_column[q] == expected_column[q], f"generated Berggren column q={q}")

for h in range(9):
    level = [(1, 2, 0)]
    for _ in range(h):
        next_level = []
        for m, n, a_parity in level:
            next_level.extend((
                (n, 2 * n - m, a_parity ^ 1),
                (n, 2 * n + m, a_parity),
                (m, n + 2 * m, a_parity),
            ))
        level = next_level
    signed = sum(epsilon(m, n) for m, n, _ in level)
    check(signed == -1, f"complete ternary row sign h={h}")

pythagorean_rows = []
for n in (6, 60, 420, 27720):
    odd_part = n
    while odd_part % 2 == 0:
        odd_part //= 2
    count = 0
    signed = 0
    for d in divisors(n):
        if d == 1:
            continue
        for a in range(1, d):
            if gcd(a, d) == 1 and ((a ^ d) & 1):
                count += 1
                signed += epsilon(a, d)
    check(count == n - (odd_part + 1) // 2, f"Pythagorean phase count N={n}")
    column_formula = sum(
        packet_value(d) if d % 2 == 0 else packet_value(d) // 2
        for d in divisors(n) if d > 1
    )
    check(signed == column_formula, f"Pythagorean signed mass N={n}")
    pythagorean_rows.append((n, count, signed))

print(f"complete primitive Berggren spinors with q<={spinor_bound}: {len(seen)}")
print("fixed-height column is S(q) (even q) or S(q)/2 (odd q): PASS")
print("every complete ternary-depth row through h=8 has signed mass -1: PASS")
print("N  primitive-phases  signed-depth-mass")
for n, count, signed in pythagorean_rows:
    print(f"{n:5d} {count:17d} {signed:18d}")
print()

# 6. Minimal hostile controls against false character/compression claims.
check(packet_sum[10] != packet_sum[2] * packet_sum[5],
      "S is not multiplicative")
q5_signs = [inverse_formula(a, 5) for a in range(1, 5)]
check(q5_signs == [1, -1, -1, 1], "first numerator-specific packet")
check(inverse_formula(2, 11) == 1, "q=11 character hostile first factor")
check(inverse_formula(4, 11) == -1, "q=11 character hostile product")
check(inverse_formula(4, 11) != inverse_formula(2, 11) ** 2,
      "depth sign is not a Dirichlet character")
for q in (3, 7, 15):
    check(all(inverse_formula(a, q) == 1 for a in range(1, q) if gcd(a, q) == 1),
          f"trivial packet character q={q}")
for p in (5, 13):
    check(all(inverse_formula(a, p) == legendre(a, p) for a in range(1, p)),
          f"Paley/Legendre packet p={p}")
digits_2_5 = canonical_cf(2, 5)[1:]
digits_4_5 = canonical_cf(4, 5)[1:]
check(len(digits_2_5) == len(digits_4_5) == 2,
      "q=5 equal positive-digit length hostile")
check(digits_2_5[0] * digits_2_5[1] == digits_4_5[0] * digits_4_5[1] == 4,
      "q=5 equal positive-digit product hostile")
check(epsilon(2, 5) == -1 and epsilon(4, 5) == 1,
      "q=5 digit-product does not determine depth sign")

print("hostiles: S(10)!=S(2)S(5); q=5 is numerator-specific;")
print("          epsilon(2/11)^2 != epsilon(4/11), so no Dirichlet character")
print("exact packets: epsilon is trivial at q=3,7,15 and Legendre at p=5,13")
print("Khinchin hostile: 2/5 and 4/5 have equal digit length/product, opposite signs")
print()
print(f"TOTAL EXACT GATES: {GATES}")
print("ALL CHECKS PASSED")
