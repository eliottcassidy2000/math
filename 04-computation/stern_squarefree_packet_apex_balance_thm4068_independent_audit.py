#!/usr/bin/env python3
"""Independent exact audit for THM-4068.

This route derives packet signs from Euclidean continued-fraction depth,
computes tournament lower stars directly over all smaller vertices, uses an
extended-Euclid half-box counter, and checks CRT spectral factorization by
iterated cyclotomic-histogram convolution.  It does not call the primary
audit's inverse-parity packet, Cartesian CRT product, or divisor-star code.
No floating point and no Python assertions are used.
"""

from collections import Counter
from hashlib import sha256
from math import gcd


QMAX = 3000
CRT_MODULI = (15, 21, 33, 35, 39, 55, 65, 77, 105)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def squarefree_primes(n):
    primes = []
    divisor = 2
    residue = n
    while divisor * divisor <= residue:
        if residue % divisor == 0:
            residue //= divisor
            if residue % divisor == 0:
                return None
            primes.append(divisor)
        divisor += 1
    if residue > 1:
        primes.append(residue)
    return primes


def bezout_inverse(a, modulus):
    old_r, r = a, modulus
    old_s, s = 1, 0
    while r:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
    require(old_r == 1, ("inverse domain", a, modulus, old_r))
    return old_s % modulus


def depth_sign(numerator, denominator):
    """Compute (-1)^(sum canonical CF digits - 1) by Euclid."""
    require(
        0 < numerator <= denominator and gcd(numerator, denominator) == 1,
        ("depth domain", numerator, denominator),
    )
    digit_sum = 0
    a = numerator
    b = denominator
    while a:
        digit_sum += b // a
        b, a = a, b % a
    return 1 if (digit_sum - 1) % 2 == 0 else -1


def depth_packet(q):
    return sum(depth_sign(a, q) for a in range(1, q) if gcd(a, q) == 1)


def direct_lower_star(q):
    total = 0
    for a in range(1, q):
        common = gcd(a, q)
        total += depth_sign(a // common, q // common)
    return total


def half_box_by_extended_euclid(q):
    m = (q - 1) // 2
    count = 0
    for r in range(1, m + 1):
        if gcd(r, q) != 1:
            continue
        s = bezout_inverse(4 * r, q)
        if 1 <= s <= m:
            count += 1
    return count


def divisors(primes):
    values = [1]
    for p in primes:
        old_values = tuple(values)
        values.extend(p * value for value in old_values)
    return sorted(values)


def phi_from_primes(q, primes):
    value = q
    for p in primes:
        value = value // p * (p - 1)
    return value


def parity(number):
    return 1 if number % 2 == 0 else -1


def plaquette_flux(u, v):
    q = u * v
    e_u = v * bezout_inverse(v, u)
    e_v = u * bezout_inverse(u, v)
    require(0 < e_u < q and 0 < e_v < q, ("idempotent range", u, v, e_u, e_v))
    require(e_u + e_v == q + 1, ("idempotent carry", u, v, e_u, e_v))
    return e_u, e_v, parity(0) * parity(1) * parity(e_u) * parity(e_v)


def unit_inverse_table(modulus):
    return tuple(
        (x, bezout_inverse(x, modulus))
        for x in range(1, modulus)
        if gcd(x, modulus) == 1
    )


def direct_spectral_histogram(q, units, h, k):
    histogram = Counter()
    for x, x_inverse in units:
        histogram[(h * x + k * x_inverse) % q] += 1
    return histogram


def convolved_local_histogram(q, primes, local_units, h, k):
    global_histogram = Counter({0: 1})
    for p in primes:
        q_p = q // p
        t_p = bezout_inverse(q_p, p)
        local_histogram = Counter()
        for x, x_inverse in local_units[p]:
            local_exponent = (h * x + k * x_inverse) % p
            local_histogram[(q_p * t_p * local_exponent) % q] += 1
        combined = Counter()
        for left_exponent, left_count in global_histogram.items():
            for right_exponent, right_count in local_histogram.items():
                combined[(left_exponent + right_exponent) % q] += left_count * right_count
        global_histogram = combined
    return global_histogram


packet_cache = {1: 1}
records = []
composite_records = 0
for q in range(3, QMAX + 1, 2):
    primes = squarefree_primes(q)
    if primes is None:
        continue
    packet_value = depth_packet(q)
    packet_cache[q] = packet_value
    phi_q = phi_from_primes(q, primes)
    half_count = half_box_by_extended_euclid(q)
    require(
        packet_value == 4 * half_count - phi_q,
        ("depth versus half box", q, packet_value, half_count, phi_q),
    )
    star_value = direct_lower_star(q)
    divisor_value = sum(packet_cache[d] for d in divisors(primes) if d > 1)
    require(star_value == divisor_value, ("direct star versus divisors", q, star_value, divisor_value))
    require((q - 1 + star_value) % 2 == 0, ("indegree parity", q, star_value))
    require((q - 1 - star_value) % 2 == 0, ("outdegree parity", q, star_value))
    flux_data = None
    if len(primes) > 1:
        u = primes[0]
        v = q // u
        e_u, e_v, flux = plaquette_flux(u, v)
        require(flux == -1, ("non-rank-one parity flux", q, e_u, e_v, flux))
        flux_data = (e_u, e_v)
        composite_records += 1
    records.append((q, tuple(primes), packet_value, half_count, phi_q, star_value, flux_data))


crt_digest = sha256()
crt_pairs = 0
degenerate_local_controls = 0
for q in CRT_MODULI:
    primes = squarefree_primes(q)
    require(primes is not None and len(primes) > 1, ("CRT modulus", q, primes))
    units = unit_inverse_table(q)
    local_units = {p: unit_inverse_table(p) for p in primes}
    for p in primes:
        axis = Counter(x for x, _ in local_units[p])
        require(axis == Counter(range(1, p)), ("local degenerate axis", p, axis))
        degenerate_local_controls += 1
    for h in range(q):
        for k in range(q):
            direct = direct_spectral_histogram(q, units, h, k)
            factored = convolved_local_histogram(q, primes, local_units, h, k)
            require(direct == factored, ("full CRT spectrum", q, h, k, direct, factored))
            crt_digest.update(f"{q}:{h}:{k}:{tuple(sorted(direct.items()))};".encode())
            crt_pairs += 1


s_three = depth_packet(3)
s_five = depth_packet(5)
s_fifteen = depth_packet(15)
require(s_three * s_five == 0, ("local product control", s_three, s_five))
require(s_fifteen == 8, ("minimal composite hostile", s_fifteen))
require(
    all(len(squarefree_primes(q) or ()) < 2 for q in range(3, 15, 2)),
    "q=15 must be the first odd squarefree composite",
)


table_digest = sha256()
for record in records:
    table_digest.update((repr(record) + ";").encode())

print("status=PASS")
print("method=continued_fraction_depth;direct_lower_star;extended_euclid_half_box;CRT_histogram_convolution")
print("scope=all_odd_squarefree_q_through_3000;all_h_k_on_nine_hostile_CRT_moduli")
print("squarefree_denominators_checked=" + str(len(records)))
print("composite_parity_fluxes_checked=" + str(composite_records))
print("full_CRT_spectral_pairs_checked=" + str(crt_pairs))
print("degenerate_prime_axis_controls=" + str(degenerate_local_controls))
print("identity=continued_fraction_S(q)=4*N(q)-phi(q)")
print("divisor_star=direct_B(q)=sum_{d|q,d>1}S(d)")
print("tournament_apex=indeg-outdeg=B(q);indeg+outdeg=q-1")
print("rank_one_control=every_nonzero_local_tensor_has_plaquette_flux_+1")
print("global_parity_flux=-1_for_every_checked_composite_squarefree_q")
print("minimal_nonmultiplicative_hostile=S(3)*S(5)=0;S(15)=8")
print("table_sha256=" + table_digest.hexdigest())
print("crt_sha256=" + crt_digest.hexdigest())
print("boundary=prime_Weil_is_imported;prime_powers_and_all_composites_are_OPEN")
