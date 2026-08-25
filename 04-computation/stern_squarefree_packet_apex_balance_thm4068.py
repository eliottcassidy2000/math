#!/usr/bin/env python3
"""Primary exact audit for THM-4068 squarefree Stern balance.

No floating point and no Python assertions are used.  The only non-elementary
input in THM-4068 is the prime Weil bound already pinned by THM-4061; this
program checks the exact CRT, divisor-star, coefficient-bound, and parity-flux
statements surrounding that imported input.
"""

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import gcd


QMAX = 5000


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime_factors(n):
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            factors.append(d)
            n //= d
            require(n % d != 0, ("not squarefree", d, n))
        d += 1
    if n > 1:
        factors.append(n)
    return factors


def is_squarefree(n):
    d = 2
    while d * d <= n:
        if n % (d * d) == 0:
            return False
        d += 1
    return True


def divisors_from_primes(primes):
    divisors = [1]
    for p in primes:
        divisors += [p * d for d in divisors]
    return sorted(divisors)


def packet(q):
    return sum(
        1 if (a + pow(a, -1, q)) % 2 == 0 else -1
        for a in range(1, q)
        if gcd(a, q) == 1
    )


def half_box(q):
    m = (q - 1) // 2
    inv4 = pow(4, -1, q)
    return sum(
        1 <= (inv4 * pow(r, -1, q)) % q <= m
        for r in range(1, m + 1)
        if gcd(r, q) == 1
    )


def totient(q):
    return sum(gcd(a, q) == 1 for a in range(1, q))


def harmonic(n):
    return sum((Fraction(1, j) for j in range(1, n + 1)), Fraction())


def harmonic_star(q):
    return 2 * harmonic(q - 1) - harmonic((q - 1) // 2)


def packet_bound_coefficient(q, primes):
    """Rational C(q) in the weaker exact check |S(q)| <= C(q)*sqrt(q)."""
    hs = harmonic_star(q)
    lbar = Fraction(1, q) + hs
    divisor_sum = sum(
        (harmonic_star(q // d) for d in divisors_from_primes(primes) if d < q),
        Fraction(),
    )
    abar = 1 + divisor_sum
    return (2 ** len(primes)) * lbar * abar


def kloosterman_exponents(q, h, k):
    return Counter(
        (h * x + k * pow(x, -1, q)) % q
        for x in range(1, q)
        if gcd(x, q) == 1
    )


def crt_product_exponents(q, primes, h, k):
    answer = Counter()
    for row in product(*(range(1, p) for p in primes)):
        exponent = 0
        for p, x in zip(primes, row):
            q_p = q // p
            t_p = pow(q_p, -1, p)
            exponent += q_p * (
                t_p * (h * x + k * pow(x, -1, p)) % p
            )
        answer[exponent % q] += 1
    return answer


def parity_tensor_flux(u, v):
    """Rank-one test on the CRT plaquette with local coordinates {0,1}."""
    q = u * v
    e_u = v * pow(v, -1, u)
    e_v = u * pow(u, -1, v)
    require(e_u + e_v == q + 1, ("CRT idempotent sum", u, v, e_u, e_v))
    lifts = (0, e_u, e_v, 1)
    return lifts, (-1) ** sum(lifts)


records = []
packet_cache = {1: 1}
tensor_obstructions = 0
for q in range(3, QMAX + 1, 2):
    if not is_squarefree(q):
        continue
    primes = prime_factors(q)
    s_q = packet(q)
    packet_cache[q] = s_q
    n_q = half_box(q)
    phi_q = totient(q)
    require(s_q == 4 * n_q - phi_q, ("half box", q, s_q, n_q, phi_q))
    coefficient = packet_bound_coefficient(q, primes)
    require(
        s_q * s_q <= q * coefficient * coefficient,
        ("squarefree packet bound", q, s_q, coefficient),
    )
    divisors = divisors_from_primes(primes)
    for d in divisors:
        if d > 1 and d not in packet_cache:
            packet_cache[d] = packet(d)
    b_q = sum(packet_cache[d] for d in divisors if d > 1)
    h_q = 2 * harmonic(q - 1)
    r = len(primes)
    b_coefficient = (
        (2**r)
        * (2**r - 1)
        * (1 + h_q)
        * (1 + (2**r - 1) * h_q)
    )
    require(
        b_q * b_q <= q * b_coefficient * b_coefficient,
        ("squarefree apex bound", q, b_q, b_coefficient),
    )
    if len(primes) > 1:
        u = primes[0]
        v = q // u
        lifts, flux = parity_tensor_flux(u, v)
        require(flux == -1, ("parity tensor obstruction", u, v, lifts, flux))
        tensor_obstructions += 1
    records.append((q, tuple(primes), s_q, n_q, phi_q, b_q))


crt_records = []
for q in (15, 21, 33, 35, 39, 55, 65, 77, 105):
    primes = prime_factors(q)
    probes = sorted({0, 1, 2, q - 1, *(p for p in primes), *(q // p for p in primes)})
    for h in probes:
        for k in probes:
            direct = kloosterman_exponents(q, h, k)
            product_counter = crt_product_exponents(q, primes, h, k)
            require(
                direct == product_counter,
                ("CRT Kloosterman", q, h, k, direct, product_counter),
            )
            crt_records.append((q, h, k, tuple(sorted(direct.items()))))


require(packet(3) * packet(5) == 0, "local product control")
require(packet(15) == 8, "minimal squarefree CRT hostile")

record_text = ";".join(
    f"{q}:{','.join(map(str, primes))}:{s}:{n}:{phi}:{b}"
    for q, primes, s, n, phi, b in records
)
crt_text = repr(crt_records)
largest_packet_ratios = sorted(
    ((Fraction(abs(s), phi), q, s, phi) for q, _, s, _, phi, _ in records),
    reverse=True,
)[:12]
largest_apex_ratios = sorted(
    ((Fraction(abs(b), q - 1), q, b) for q, _, _, _, _, b in records),
    reverse=True,
)[:12]

print("status=PASS")
print("scope=all_odd_squarefree_q_through_5000;exact_CRT_spectral_probes")
print("squarefree_denominators_checked=" + str(len(records)))
print("crt_spectral_pairs_checked=" + str(len(crt_records)))
print("parity_tensor_obstructions_checked=" + str(tensor_obstructions))
print("identity=S(q)=4*N(q)-phi(q)")
print("exact_bound=|S(q)|<=2^omega*sqrt(q)*Lbar(q)*Abar(q)")
print("weak_bound_exactly_checked=replace_Abar_radicals_by_1+sum_Hstar(q/d)")
print("squarefree_apex_corollary=|B(q)|<=q^(1/2+o(1))")
print("minimal_nonmultiplicative_hostile=S(3)*S(5)=0;S(15)=8")
print("CRT_parity_plaquette=eta(0,0)*eta(1,0)*eta(0,1)*eta(1,1)=-1")
print("largest_packet_ratios=" + repr(largest_packet_ratios))
print("largest_apex_ratios=" + repr(largest_apex_ratios))
print("table_sha256=" + sha256(record_text.encode()).hexdigest())
print("crt_sha256=" + sha256(crt_text.encode()).hexdigest())
print("boundary=prime_Weil_is_imported;no_prime_power_or_all_composite_claim")
