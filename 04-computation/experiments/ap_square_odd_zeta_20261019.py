#!/usr/bin/env python3
"""Exact controls for the signed AP-square/Pythagorean Dirichlet bridge.

Standard library only. Run from any directory; no repository writes.
All checks remain active under python -O. Decimal displays are enclosed using
elementary series-tail bounds; the mathematical identities have proofs in the
companion note rather than being inferred from these finite checks.
"""

from collections import Counter
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, getcontext
from fractions import Fraction
from math import comb, factorial, gcd, isqrt

getcontext().prec = 80
CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def ppts(limit):
    ans = []
    for n in range(2, isqrt(limit) + 1):
        for m in range(1, min(n, isqrt(limit - n * n) + 1)):
            if (m + n) % 2 and gcd(m, n) == 1:
                a, b, c = n * n - m * m, 2 * m * n, n * n + m * m
                ans.append((max(a, b), min(a, b), c))
    return sorted(ans, key=lambda t: (t[2], t[0]))


def triangle(a, b, c, g=1):
    return tuple(2 * g * g * t for t in (c * c - 4 * a * b, c * c, c * c + 4 * a * b))


def sieve(limit):
    spf = list(range(limit + 1))
    for p in range(2, isqrt(limit) + 1):
        if spf[p] == p:
            for k in range(p * p, limit + 1, p):
                if spf[k] == k:
                    spf[k] = p
    return spf


def factor(n, spf):
    result = {}
    while n > 1:
        p = spf[n]
        result[p] = result.get(p, 0) + 1
        n //= p
    return result


def bernoulli(limit):
    values = [Fraction(1)]
    for n in range(1, limit + 1):
        values.append(-sum(Fraction(comb(n + 1, k)) * values[k] for k in range(n)) / (n + 1))
    return values


def euler(limit):
    values = [0] * (limit + 1)
    values[0] = 1
    for n in range(2, limit + 1, 2):
        values[n] = -sum(comb(n, k) * values[k] for k in range(0, n, 2))
    return values


def dec(r):
    return Decimal(r.numerator) / Decimal(r.denominator)


def complex_power(a, b, exponent):
    real, imag = 1, 0
    for _ in range(exponent):
        real, imag = real * a - imag * b, real * b + imag * a
    return real, imag


def angular_cos(a, b, exponent):
    # exponent is even; exact cosine of exponent * arg(a+ib).
    return Fraction(complex_power(a, b, exponent)[0], (a * a + b * b) ** (exponent // 2))


def atan_interval(q, terms=80):
    # Alternating rational series for arctan(1/q); even terms gives lower sum.
    low = sum((Fraction((-1) ** k, (2 * k + 1) * q ** (2 * k + 1)) for k in range(terms)), Fraction(0))
    return low, low + Fraction(1, (2 * terms + 1) * q ** (2 * terms + 1))


def zeta_interval(s, terms=10000):
    partial = sum((Decimal(1) / Decimal(k) ** s for k in range(1, terms + 1)), Decimal(0))
    return (partial + Decimal(terms + 1) ** (1 - s) / (s - 1),
            partial + Decimal(terms) ** (1 - s) / (s - 1))


def beta_interval(s, terms=10000):
    partial = sum((Decimal((-1) ** k) / Decimal(2 * k + 1) ** s for k in range(terms)), Decimal(0))
    return partial, partial + Decimal(1) / Decimal(2 * terms + 1) ** s


def main():
    limit, direct_height, convolution_limit = 100000, 300, 10000
    triples = ppts(limit)
    counts = Counter(c for a, b, c in triples)
    positive = [(a, b, c) for a, b, c in triples if c * c > 4 * a * b]
    pc = Counter(c for a, b, c in positive)
    spf = sieve(limit)

    for c in range(1, limit + 1):
        fac = factor(c, spf)
        predicted = 2 ** (len(fac) - 1) if fac and all(p % 4 == 1 for p in fac) else 0
        check(counts[c] == predicted, f"Gaussian multiplicity C={c}")
    check(counts[1] == 0 and counts[5] == 1 and counts[25] == 1 and counts[65] == 2, "axis/power/orientation controls")
    check(counts[1105] == 4 and counts[15] == 0 and counts[10] == 0, "split/inert/ramified controls")
    check(positive[0] == (40, 9, 41), "first positive base")
    check(triangle(40, 9, 41) == (482, 3362, 6242), "first positive labels")
    check(triangle(4, 3, 5) == (-46, 50, 146), "first signed labels")

    # Independent census from the pair-sum-root equation p^2+r^2=8h^2.
    direct = set()
    direct_counts = Counter()
    for h in range(1, direct_height + 1):
        y = 2 * h * h
        for p in range(1, 2 * h):
            r2 = 8 * h * h - p * p
            r = isqrt(r2)
            if r * r != r2:
                continue
            xyz = (p * p - y, y, 3 * y - p * p)
            check(p % 2 == r % 2 == 0, "pair roots even")
            check((r + p) % 4 == (r - p) % 4 == 0, "inverse integral")
            aa, bb = (r + p) // 4, (r - p) // 4
            g = gcd(gcd(aa, bb), h)
            a, b, c = aa // g, bb // g, h // g
            check(triangle(a, b, c, g) == xyz, "inverse classification")
            check(gcd(xyz[0], xyz[1]) == gcd(xyz[0], xyz[2]) == gcd(xyz[1], xyz[2]) == 2 * g * g, "pairwise vertex gcd")
            check(gcd(p, 2 * h) == gcd(p, r) == gcd(2 * h, r) == 2 * g, "pairwise root gcd")
            direct.add(xyz)
            direct_counts[h] += 1
    scaled = set()
    for a, b, c in triples:
        if c > direct_height:
            break
        for g in range(1, direct_height // c + 1):
            scaled.add(triangle(a, b, c, g))
    check(direct == scaled, "independent signed-AP census equals scaled-PPT census")
    for h in range(1, direct_height + 1):
        convolution = sum(counts[c] for c in range(1, h + 1) if h % c == 0)
        check(direct_counts[h] == convolution, f"root-height convolution h={h}")
    check(direct_counts[65] == 4, "primitive/nonprimitive same-height multiplicity")
    for h in range(1, direct_height + 1):
        fac = factor(h, spf)
        product = 1
        for p, e in fac.items():
            if p % 4 == 1:
                product *= 2 * e + 1
        check(direct_counts[h] == (product - 1) // 2, f"closed all-height coefficient h={h}")

    # Independent Dirichlet coefficient calculation for ζβ/((1+2^-s)ζ(2s)).
    zeta_beta = [0] * (convolution_limit + 1)
    for d in range(1, convolution_limit + 1, 2):
        chi = 1 if d % 4 == 1 else -1
        for n in range(d, convolution_limit + 1, d):
            zeta_beta[n] += chi
    squarefree_part = [0] * (convolution_limit + 1)
    for d in range(1, isqrt(convolution_limit) + 1):
        fac = factor(d, spf)
        mu = 0 if any(e > 1 for e in fac.values()) else (-1) ** len(fac)
        if mu:
            for n in range(1, convolution_limit // (d * d) + 1):
                squarefree_part[d * d * n] += mu * zeta_beta[n]
    full = [0] * (convolution_limit + 1)
    power, sign = 1, 1
    while power <= convolution_limit:
        for n in range(1, convolution_limit // power + 1):
            full[power * n] += sign * squarefree_part[n]
        power, sign = 2 * power, -sign
    for n in range(1, convolution_limit + 1):
        check(full[n] == (1 if n == 1 else 2 * counts[n]), f"Euler-product coefficient n={n}")

    # Independent Gaussian lattice enumeration tests the Epstein/Eisenstein
    # factors 4, 2, the axis terms, and the ramified norm-2 contribution.
    gaussian_limit = 1000
    gaussian_all, gaussian_primitive = Counter(), Counter()
    radius = isqrt(gaussian_limit)
    for c in range(-radius, radius + 1):
        for d in range(-radius, radius + 1):
            norm = c * c + d * d
            if 1 <= norm <= gaussian_limit:
                gaussian_all[norm] += 1
                if gcd(c, d) == 1:
                    gaussian_primitive[norm] += 1
    for n in range(1, gaussian_limit + 1):
        check(gaussian_all[n] == 4 * zeta_beta[n], f"full Gaussian Epstein coefficient n={n}")
        # E(i,s) = 2(1+2^-s)(1+2P(s)).
        expected_e = 2 * full[n] + (2 * full[n // 2] if n % 2 == 0 else 0)
        check(gaussian_primitive[n] == 2 * expected_e, f"primitive Eisenstein coefficient n={n}")
    check(gaussian_primitive[1] == gaussian_primitive[2] == 4 and gaussian_primitive[4] == 0, "Eisenstein axis/ramified controls")

    # Gaussian angular local factors, tested using exact rational cosines.
    angular_limit = 1000
    split_reps = {}
    for p in range(2, angular_limit + 1):
        if spf[p] == p and p % 4 == 1:
            for v in range(1, isqrt(p // 2) + 1):
                u = isqrt(p - v * v)
                if u * u + v * v == p:
                    split_reps[p] = (u, v)
                    break
    for k in range(4):
        weighted = [Fraction(0)] * (angular_limit + 1)
        weighted[1] = 1
        for a, b, c in triples:
            if c > angular_limit:
                break
            weighted[c] += 2 * angular_cos(a, b, 4 * k)
        for c in range(1, angular_limit + 1):
            fac = factor(c, spf)
            predicted = Fraction(1)
            for p, e in fac.items():
                if p % 4 != 1:
                    predicted = 0
                    break
                u, v = split_reps[p]
                predicted *= 2 * angular_cos(u, v, 8 * k * e)
            check(weighted[c] == predicted, f"angular Gaussian coefficient k={k}, C={c}")
    check(pc[5] == pc[13] == 0 and pc[65] == 1, "positive angular filter destroys multiplicativity")

    B, E = bernoulli(18), euler(8)
    expected = {5: Fraction(4725, 32), 7: Fraction(8243235, 5504), 9: Fraction(334084125375, 22459904)}
    a5, b5 = atan_interval(5)
    a239, b239 = atan_interval(239)
    pi_low, pi_high = dec(16 * a5 - 4 * b239), dec(16 * b5 - 4 * a239)
    padding = Decimal("1e-65")  # Dominates accumulated 80-digit rounding.

    print(f"Primitive PPT universe: C <= {limit}; count={len(triples)}; positive={len(positive)}")
    print(f"Independent signed-AP universe: root height <= {direct_height}; count={len(direct)}; positive={sum(x > 0 for x, y, z in direct)}")
    print(f"Independent Euler-product coefficient universe: n <= {convolution_limit}")
    print(f"Independent Gaussian lattice coefficient universe: norm <= {gaussian_limit}; full and primitive")
    print(f"Exact angular coefficient universe: C <= {angular_limit}; modes k=0,1,2,3")
    print("First signed base: (-46, 50, 146); first positive base: (482, 3362, 6242)")
    print("Primitive coefficient controls: C=1:0, 5:1, 10:0, 15:0, 25:1, 65:2, 1105:4")
    for s in (5, 7, 9):
        K = Fraction(abs(E[s - 1]) * factorial(2 * s), 2 ** (2 * s + 1) * (2 ** s + 1) * factorial(s - 1)) / abs(B[2 * s])
        check(K == expected[s], f"exact odd-zeta rational s={s}")
        total = sum((Decimal(count) / Decimal(c) ** s for c, count in counts.items()), Decimal(0))
        cone = sum((Decimal(count) / Decimal(c) ** s for c, count in pc.items()), Decimal(0))
        tail = Decimal(s) / Decimal(2 * (s - 1)) * Decimal(limit) ** (1 - s)
        zl, zh = zeta_interval(s)
        z2l, z2h = zeta_interval(2 * s)
        bl, bh = beta_interval(s)
        ll, lh = dec(K) * zl / pi_high ** s - Decimal("0.5"), dec(K) * zh / pi_low ** s - Decimal("0.5")
        factor2 = 1 + Decimal(2) ** (-s)
        el, eh = (zl * bl / (factor2 * z2h) - 1) / 2, (zh * bh / (factor2 * z2l) - 1) / 2
        check(total - padding <= lh and ll <= total + tail + padding, f"odd-zeta formula interval s={s}")
        check(total - padding <= eh and el <= total + tail + padding, f"Euler-product value interval s={s}")
        check(max(ll, el) <= min(lh, eh) + padding, f"two analytic enclosures overlap s={s}")
        check(cone + tail < total / 3, f"positive-cone one-third false at s={s}")
        print(f"s={s}: K={K}; P partial={total:.24E}; P+ partial={cone:.24E}; tail bound={tail:.3E}")
        rl = (cone - padding) / (total + tail + padding)
        rh = (cone + tail + padding) / (total - padding)
        rl = rl.quantize(Decimal(1).scaleb(rl.adjusted() - 18), rounding=ROUND_FLOOR)
        rh = rh.quantize(Decimal(1).scaleb(rh.adjusted() - 18), rounding=ROUND_CEILING)
        print(f"  outward-rounded cone/full interval: [{rl:.18E}, {rh:.18E}]")
    print(f"ALL CHECKS PASS: {CHECKS}")


if __name__ == "__main__":
    main()
