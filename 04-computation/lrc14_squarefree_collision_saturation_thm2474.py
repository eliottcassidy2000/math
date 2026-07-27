#!/usr/bin/env python3
"""Exact companion for THM-2474.

All algebra is dependency-free. Primitive Fourier values are tested by exact
reduction modulo cyclotomic polynomials; energy is evaluated with integral
Ramanujan sums and Fraction arithmetic.
"""

from fractions import Fraction as F
from math import gcd, isqrt


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def divisors(n):
    small, large = [], []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            small.append(d)
            if d * d != n:
                large.append(n // d)
    return small + large[::-1]


def prime_factors(n):
    factors = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            factors.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        factors.append(n)
    return factors


def mobius(n):
    sign = 1
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            sign = -sign
            if n % p == 0:
                return 0
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        sign = -sign
    return sign


def euler_phi(n):
    value = n
    for p in prime_factors(n):
        value = value // p * (p - 1)
    return value


def trim(poly):
    poly = list(poly)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def poly_mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def poly_divmod(num, den):
    num = trim([F(x) for x in num])
    den = trim([F(x) for x in den])
    require(den != [0], "zero polynomial divisor")
    if len(num) < len(den):
        return [F(0)], num
    quotient = [F(0)] * (len(num) - len(den) + 1)
    while len(num) >= len(den) and num != [0]:
        shift = len(num) - len(den)
        coefficient = num[-1] / den[-1]
        quotient[shift] += coefficient
        for j, value in enumerate(den):
            num[shift + j] -= coefficient * value
        num = trim(num)
    return trim(quotient), trim(num)


_cyclotomic_cache = {1: [F(-1), F(1)]}


def cyclotomic(n):
    if n in _cyclotomic_cache:
        return _cyclotomic_cache[n]
    numerator = [F(-1)] + [F(0)] * (n - 1) + [F(1)]
    denominator = [F(1)]
    for d in divisors(n):
        if d < n:
            denominator = poly_mul(denominator, cyclotomic(d))
    quotient, remainder = poly_divmod(numerator, denominator)
    require(remainder == [0], f"cyclotomic division failed at {n}")
    require(all(x.denominator == 1 for x in quotient), "nonintegral cyclotomic")
    _cyclotomic_cache[n] = quotient
    return quotient


def primitive_remainder(values, colour):
    n = len(values)
    raw = [F(0)] * n
    for s, value in enumerate(values):
        raw[(-colour * s) % n] += value
    _, remainder = poly_divmod(raw, cyclotomic(n))
    return tuple(trim(remainder))


def ramanujan(n, m):
    g = gcd(n, abs(m))
    return sum(d * mobius(n // d) for d in divisors(g))


def root_packet_correlation(n):
    # A=delta_0, F=unit mask gives C=unit mask (up to orientation).
    a = [F(int(r == 0)) for r in range(n)]
    f = [F(int(gcd(r, n) == 1)) for r in range(n)]
    correlation = []
    for shift in range(n):
        correlation.append(sum(a[(r + shift) % n] * f[r] for r in range(n)))
    require(all(correlation[s] == int(gcd(s, n) == 1) for s in range(n)),
            f"unit correlation at {n}")

    # Each codimension-one block product is zero.
    for p in prime_factors(n):
        for residue in range(p):
            a_block = sum(a[r] for r in range(n) if r % p == residue)
            f_block = sum(f[r] for r in range(n) if r % p == residue)
            require(a_block * f_block == 0, f"predecessor block p={p}, n={n}")
    require(sum(a) * sum(f) > 0, "positive top service")
    return correlation


def spectral_invoice(values):
    n = len(values)
    total = sum(values, F(0))
    collision = total / (n * n)
    signed = sum(
        values[s] * ramanujan(n, s) for s in range(n)
    ) / (n * n)
    energy = sum(
        values[s] * values[t] * ramanujan(n, t - s)
        for s in range(n) for t in range(n)
    ) / (n**4)
    return collision, signed, energy


def squarefree_controls():
    moduli = [6, 10, 15, 30, 42, 91]
    checked_colours = 0
    for n in moduli:
        require(mobius(n) != 0, f"{n} must be squarefree")
        base = root_packet_correlation(n)
        # Use nonuniform positive rational weights on all unit shifts to
        # hostile-test that the proof is not relying on a constant mask.
        values = [
            F(1 + ((s * s + 3 * s + 7) % 17), 19) if gcd(s, n) == 1 else F(0)
            for s in range(n)
        ]
        require(all(base[s] != 0 or values[s] == 0 for s in range(n)),
                "support containment")
        for colour in range(n):
            if gcd(colour, n) == 1:
                require(any(primitive_remainder(values, colour)),
                        f"primitive colour vanished: n={n}, k={colour}")
                checked_colours += 1
        collision, signed, energy = spectral_invoice(values)
        require(signed == mobius(n) * collision, f"signed ledger at {n}")
        require(energy >= collision * collision / euler_phi(n), f"energy at {n}")
    return tuple(moduli), checked_colours


def nonsquarefree_hostiles():
    moduli = [4, 8, 9, 12]
    checked_colours = 0
    for n in moduli:
        require(mobius(n) == 0, f"{n} must be nonsquarefree")
        values = root_packet_correlation(n)
        for colour in range(n):
            if gcd(colour, n) == 1:
                require(not any(primitive_remainder(values, colour)),
                        f"nonsquarefree hostile fired: n={n}, k={colour}")
                checked_colours += 1
        collision, signed, _ = spectral_invoice(values)
        require(collision > 0 and signed == 0, f"hostile ledger at {n}")
    return tuple(moduli), checked_colours


def fano_control():
    n = 91
    collision = F(1, 6365350748643250000000000)
    weight = F(1, 1537338666500000000000)
    values = [F(0)] * n
    values[34] = weight
    values[57] = weight
    require(gcd(34, n) == gcd(57, n) == gcd(57 - 34, n) == 1, "Fano shifts")
    require(sum(values) == n * n * collision, "corrected Fano normalization")
    for colour in range(1, n):
        if gcd(colour, n) == 1:
            require(any(primitive_remainder(values, colour)), f"Fano colour {colour}")
    i, signed, energy = spectral_invoice(values)
    require(i == collision and signed == collision, "Fano signed ledger")
    require(energy >= collision * collision / 72, "Fano energy")
    return collision


def is_prime(n):
    if n < 2:
        return False
    for p in range(2, isqrt(n) + 1):
        if n % p == 0:
            return False
    return True


def three_prime_and_incidence_control():
    ell = 65537
    require(is_prime(ell) and ell % 16 == 1 and gcd(ell, 91) == 1, "amplifier prime")
    table = [
        ((3, 1, 1), 1007893523, 1000000000000, True),
        ((2, 1, 1), 77530271, 3000000000000, False),
        ((3, 0, 1), 143984789, 1000000000000, False),
        ((3, 1, 0), 15379, 1000000000000, False),
    ]
    for cell, scale, right, expected in table:
        require((1360 * scale > right) == expected, f"corner inequality {cell}")

    epsilon = F(1, 10**12)

    def two_interval_nonzero(frequency, half_width):
        # Both center pairs in the carrier vanish exactly at n=4 mod 8.
        center_ok = frequency % 8 != 4
        width_phase = 2 * frequency * half_width
        width_ok = width_phase.denominator != 1
        require(F(0) < width_phase < 1, "width phase range")
        return center_ok and width_ok

    q0 = 13**2 * 17
    require(q0 == 2873 and q0 % 8 == 1, "q=17 predecessor frequency")
    require(two_interval_nonzero(q0, 169 * epsilon), "current q=17")
    require(two_interval_nonzero(q0, epsilon), "source q=17")

    x = 13**3 * 17
    c3 = 742586
    y = x + c3
    require(c3 == 2 * 13**5, "deep speed")
    require(x == 37349 and y == 13**3 * 355, "incident endpoints")
    require(y // 13 == 59995 and (y // 13) % 8 == 3, "second source frequency")
    require(two_interval_nonzero(y // 13, epsilon), "second source endpoint")
    require((x // 13**3) % 13 == (y // 13**3) % 13 == 4, "root four")
    require(y - x == c3 and gcd(1, 91) == 1, "unit incident edge")
    require(gcd(17, 91 * ell) == 1, "prescribed primitive class")
    return ell, x, y


def main():
    squarefree, squarefree_colours = squarefree_controls()
    nonsquarefree, hostile_colours = nonsquarefree_hostiles()
    fano_i = fano_control()
    ell, x, y = three_prime_and_incidence_control()
    print("THM-2474 exact companion: PASS")
    print(f"squarefree_moduli={squarefree}; primitive_colours_checked={squarefree_colours}")
    print(f"nonsquarefree_hostiles={nonsquarefree}; vanished_colours_checked={hostile_colours}")
    print(f"corrected_Fano_corner_I={fano_i}; all_72_colours_nonzero")
    print(f"three_prime_ell={ell}; full_corner=(3,1,1); faces_zero=True")
    print(f"explicit_root4_edge: X={x}; Y={y}; multiplier=1; unit_mod_91=True")
    print("Ramanujan signed ledgers, energy floors, and squarefree sharpness verified")


if __name__ == "__main__":
    main()
