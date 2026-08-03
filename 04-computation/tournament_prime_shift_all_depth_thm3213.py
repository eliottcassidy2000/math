#!/usr/bin/env python3
"""Independent exact audit of the post-THM-3213 prime-shift argument."""

from collections import defaultdict
from functools import lru_cache
from math import comb


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


BASE = {
    (1, 0, 0): 1,
    (0, 1, 0): 1,
    (0, 0, 1): 1,
    (1, 1, 0): 1,
    (1, 0, 1): 1,
    (0, 1, 1): 1,
    (1, 1, 1): 3,
}


@lru_cache(maxsize=None)
def numerator_power(d):
    polynomial = {(0, 0, 0): 1}
    for _ in range(d):
        product = defaultdict(int)
        for left, left_value in polynomial.items():
            for right, right_value in BASE.items():
                exponent = tuple(left[i] + right[i] for i in range(3))
                product[exponent] += left_value * right_value
        polynomial = dict(product)
    return polynomial


def transitive_profile(n, modulus=None):
    """Return F_{T_n}(c)=c! S(n,c), c=0..n."""
    row = [1]
    for size in range(1, n + 1):
        nxt = [0] * (size + 1)
        for c in range(1, size + 1):
            same = row[c] if c < len(row) else 0
            merged = row[c - 1]
            value = c * (same + merged)
            nxt[c] = value if modulus is None else value % modulus
        row = nxt
    return row


def c3_ordered_coordinate(profile, d, modulus=None):
    """Compute F_{C3[T,T,T]}(d) from the quotient band."""
    n = len(profile) - 1
    total = 0
    for (a, b, c), coefficient in numerator_power(d).items():
        for q in range(n + 1):
            if q + max(a, b, c) > n:
                break
            term = (
                coefficient
                * comb(d + q - 1, d - 1)
                * profile[q + a]
                * profile[q + b]
                * profile[q + c]
            )
            total += term
            if modulus is not None:
                total %= modulus
    return total


def p_adic_valuation(value, p):
    valuation = 0
    while value % p == 0:
        valuation += 1
        value //= p
    return valuation


def prime_sieve(limit):
    flag = bytearray(b"\x01") * (limit + 1)
    flag[0:2] = b"\x00\x00"
    for p in range(2, int(limit**0.5) + 1):
        if flag[p]:
            start = p * p
            flag[start : limit + 1 : p] = b"\x00" * (((limit - start) // p) + 1)
    return [p for p in range(2, limit + 1) if flag[p]]


def finite_difference_power_mod(p, c, modulus):
    return sum(
        (-1) ** (c - j) * comb(c, j) * pow(j, p, modulus)
        for j in range(c + 1)
    ) % modulus


# The d=4 p-index first layer.
layer_weights = {2: 36, 3: 84, 4: 72, 5: 12}
p_index_rows = []
for p in (5, 7, 11, 13, 17, 19, 23, 29):
    modulus = p * p
    profile = transitive_profile(p, modulus)
    full = c3_ordered_coordinate(profile, 4, modulus)
    layer = sum(
        weight * finite_difference_power_mod(p, c, modulus)
        for c, weight in layer_weights.items()
    ) % modulus
    require(full == layer, ("p-index first layer", p, full, layer))
    p_index_rows.append((p, (full // p) % p))

for p in (7, 11, 13, 17, 19, 23, 29):
    quotient = lambda a: ((pow(a, p - 1, p * p) - 1) // p) % p
    fermat_layer = 12 * (
        24 * quotient(2) - 21 * quotient(3) + 5 * quotient(5)
    ) % p
    require(
        fermat_layer == dict(p_index_rows)[p],
        ("Fermat quotient layer", p, fermat_layer, dict(p_index_rows)[p]),
    )

exact_exception_valuations = []
for p in (5, 13):
    exact = c3_ordered_coordinate(transitive_profile(p), 4)
    exact_exception_valuations.append((p, p_adic_valuation(exact, p)))
require(exact_exception_valuations == [(5, 3), (13, 2)], exact_exception_valuations)


# Finite scout only: zeros of the rational-base Fermat quotient through 10^6.
scout_limit = 1_000_000
scout_zeros = []
for p in prime_sieve(scout_limit):
    if p <= 5:
        continue
    quotient = lambda a: ((pow(a, p - 1, p * p) - 1) // p) % p
    layer = 12 * (24 * quotient(2) - 21 * quotient(3) + 5 * quotient(5)) % p
    if layer == 0:
        scout_zeros.append(p)
require(scout_zeros == [13], scout_zeros)


# Exact prime-shift controls and minimal positive constants C_d=F_{U_m}(d).
control_primes = (5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
constants = []
shift_checks = 0
for d in range(1, 13):
    m = (d + 2) // 3
    shift = m - 1
    small_profile = transitive_profile(m)
    constant = c3_ordered_coordinate(small_profile, d)
    require(constant > 0, ("positive constant", d, m, constant))
    constants.append((d, m, shift, constant))
    for p in control_primes:
        if p <= shift:
            continue
        shifted_profile = transitive_profile(p + shift, p)
        padded_small = [x % p for x in small_profile] + [0] * (p + shift - m)
        require(
            shifted_profile == padded_small,
            ("profile shift", d, p, shift),
        )
        shifted_output = c3_ordered_coordinate(shifted_profile, d, p)
        require(
            shifted_output == constant % p,
            ("output shift", d, p, shifted_output, constant % p),
        )
        shift_checks += 1

require(constants[3] == (4, 2, 1, 1944), constants[3])
for p in control_primes:
    if p > 3:
        require(1944 % p != 0, ("d=4 all-prime control", p))


print("TOURNAMENT PRIME-SHIFT ALL-DEPTH AUDIT")
print("status=exact_controls;general_prime_shift_proof_is_algebraic")
print("d4_p_index_mod_p2=36F_p(2)+84F_p(3)+72F_p(4)+12F_p(5)")
print("d4_p_index_closed=12*(-4+8*2^p-7*3^p+4^p+5^p) mod p^2")
print("d4_p_index_quotient=12*(24q_p(2)-21q_p(3)+5q_p(5)) mod p")
print(f"p_index_first_quotients={tuple(p_index_rows)}")
print(f"p_index_exact_exception_valuations={tuple(exact_exception_valuations)}")
print(f"FINITE_SCOUT_ONLY_limit={scout_limit};zero_first_layers={tuple(scout_zeros)}")
print("prime_shift_identity=F_T_(p+s)(c)=F_T_(s+1)(c) mod p")
print("output_shift_identity=F_U_(p+s)(d)=F_U_(s+1)(d) mod p")
print(f"shift_controls={shift_checks}")
print("minimal_constants_d1_to12=" + str(tuple(constants)))
print("d4_shift=1;constant=1944=2^3*3^5;new_p3_denominator_at_r=p+1_for_all_p>3")
print("verdict=singly_factorial_normalized_non_C_finite_for_every_fixed_d")
