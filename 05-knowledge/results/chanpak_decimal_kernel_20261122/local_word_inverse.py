"""A short-word inverse for the modular binary path reading, with physical Farey controls."""
from argparse import ArgumentParser
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import gcd, prod

from fermat_receiver import A, D, EDGE, TARGET, address_height, flipped, join, p, receiver, tagged_data

SUPERSINGULAR = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71)
OGG = prod(SUPERSINGULAR)


def advance(state, word, modulus):
    x, z = state
    for bit in word:
        if bit == "0":
            z = (z + (x + 1) // 2) % modulus
        x = (2 * x + int(bit)) % (2 * modulus)
    return x, z


def digits(value, width):
    return format(value, f"0{width}b")


@dataclass(frozen=True)
class Switch:
    zero: str
    one: str
    slope: int
    const: int


def switch(left, right, width):
    a, b = digits(left, width), digits(right, width)
    return Switch(a, b, int(a[::-1], 2) - int(b[::-1], 2), (p(right) - p(left)) // 2)


def switches(modulus, use_shortcut=True):
    """The short one-switch search is optional; the three-switch fallback is proved."""
    period = 2 * modulus
    if use_shortcut:
        for multiplier in range(1, 33):
            v = multiplier * period
            s = switch(0, v, v.bit_length())
            divisor = gcd(s.slope, modulus)
            if divisor == 1 or (divisor == 3 and s.const % 3):
                return [s]

    power = period & -period
    odd = period // power
    u = (5 - odd) % 8
    v = u + odd
    width = ((v + 1) * power).bit_length()
    bank = [switch(u * power, v * power, width),
            switch((u + 1) * power, (v + 1) * power, width),
            switch(0, period, period.bit_length())]
    divisor = gcd(modulus, *(s.slope for s in bank))
    assert divisor in (1, 3)
    if divisor == 3 and sum(s.const for s in bank) % 3 == 0:
        i = next(i for i, s in enumerate(bank) if s.const % 3)
        s = bank[i]
        bank[i] = Switch(s.one, s.zero, -s.slope, -s.const)
    return bank


def coefficients(modulus, bank):
    """Return g and a Bezout expression for g in the switch slopes modulo M."""
    g = modulus
    coefficients = []
    for s in bank:
        reduced = gcd(g, s.slope)
        period = g // reduced
        u = pow(s.slope // reduced, -1, period) if period > 1 else 0
        v = (reduced - u * s.slope) // g
        coefficients = [v * c % modulus for c in coefficients] + [u]
        g = reduced
    return g, coefficients


def steer(modulus, incoming, target, use_shortcut=True):
    """Append a word steering (X mod 2M, f(X)/2 mod M) from incoming to target."""
    if modulus < 1:
        raise ValueError("positive modulus required")
    bank = switches(modulus, use_shortcut)
    g, bezout = coefficients(modulus, bank)
    const = sum(s.const for s in bank)
    period = 2 * modulus
    width = (period - 1).bit_length()
    desired_x, desired_z = target
    current = incoming

    def load(to):
        nonlocal current
        word = digits((to - (current[0] << width)) % period, width)
        current = advance(current, word, modulus)
        return word

    options = []
    negative = 0
    for j in range((modulus - 1).bit_length()):
        step = 1 << j
        if (step - const) % g:
            step = -step
            negative -= step
        assert (step - const) % g == 0
        alternatives = [[], []]
        for c, s in zip(bezout, bank):
            d = c * ((step - const) // g) % modulus
            loader = load(2 * d)
            alternatives[0].extend((loader, s.zero))
            alternatives[1].extend((loader, s.one))
            current = advance(current, s.zero, modulus)
        options.append(("".join(alternatives[0]), "".join(alternatives[1]), step < 0))

    tail = load(desired_x)
    rank = (desired_z - current[1] + negative) % modulus
    pieces = []
    for j, (a, b, neg) in enumerate(options):
        bit = bool((rank >> j) & 1) ^ neg
        pieces.append(b if bit else a)
    answer = "".join(pieces) + tail
    assert advance(incoming, answer, modulus) == (desired_x % period, desired_z % modulus)
    return answer


def farey(n, target, modulus, parity=0):
    ad = join(tagged_data(A, 12), tagged_data(flipped(D, *EDGE), 35))
    chart = join(ad, receiver(n))
    h, _, _, t, _, r, k = chart
    assert gcd(h, t) == 1
    assert (target - k - (t + r) * parity) % 2 == 0
    half = (target - k - (t + r) * parity) // 2
    u = half * pow(t, -1, h) % h
    z = (half - t * u) // h
    # Starting with a one gives a canonical nonempty left address.
    suffix = steer(modulus, (1, 0), ((2 * u + parity) % (2 * modulus), z % modulus))
    word = "1" + suffix
    h1 = address_height(chart, int(word, 2), 0)
    assert (h1 - target) % (2 * modulus) == 0
    return word


def density(chart, target):
    """Normalized density at all Ogg primes, conditional on each parity."""
    h, q, c, t, e, r, k = chart
    ratios = [Fraction(1)] * 4
    for prime in SUPERSINGULAR[1:]:
        if h % prime or q % prime == 0:
            if h % prime == q % prime == 0:
                assert gcd(t, prime) == 1
            continue
        assert h % (prime * prime)
        for a in (0, 1):
            for b in (0, 1):
                difference = (q * (target - k - r * a - e * b) + c * t) % prime
                count = 2 * prime - 1 if difference == 0 else prime - 1
                direct = sum((k + t * x + c * y + q * x * y + r * a + e * b - target) % prime == 0
                             for x in range(prime) for y in range(prime))
                assert direct == count
                ratios[2 * a + b] *= Fraction(count, prime)
    return ratios


def controls():
    for m in range(1, 129):
        for incoming in [(0, 0), (2 * m - 1, m - 1)]:
            for target in [(0, 0), (m - 1, 1 % m), (2 * m - 1, m - 1)]:
                steer(m, incoming, target, use_shortcut=False)
    ad = join(tagged_data(A, 12), tagged_data(flipped(D, *EDGE), 35))
    for n in range(4, 33, 4):
        ratios = density(join(ad, receiver(n)), TARGET)
        average = sum(ratios) / 4
        expected = Fraction(49, 34) if n == 12 else (
            Fraction(448, 493) if n in (4, 20, 28) else Fraction(28, 29))
        assert average == expected
        print("n", n, "local target density relative to uniform", average)
        for exponent in (1, 2):
            m = OGG ** exponent
            x = farey(n, TARGET, m)
            print(n, exponent, len(x), sha256(x.encode()).hexdigest()[:16])


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, choices=range(4, 33, 4))
    parser.add_argument("--modulus", type=int, default=OGG)
    parser.add_argument("--target", type=int, default=TARGET)
    parser.add_argument("--parity", type=int, choices=(0, 1), default=0)
    args = parser.parse_args()
    if args.n:
        print(farey(args.n, args.target, args.modulus, args.parity))
    else:
        controls()
