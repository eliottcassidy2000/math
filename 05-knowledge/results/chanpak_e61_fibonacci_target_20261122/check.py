"""Exact Mersenne carry inverse and local arc certificates at exponent 61."""
import argparse
from fractions import Fraction
from math import factorial
from pathlib import Path
import sys

PACKETS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PACKETS / "chanpak_gap_star_ports_20261122"))
sys.path.insert(0, str(PACKETS / "chanpak_e31_physical_20261122"))
from e31_certificates import CORE, ordinal
from replay import count_from_profile, flipped, p, profile, star_rows, strong


def extend_reader(source, length, contacts):
    """Attach a new root and transitive tail, using the given free contact mask."""
    k = len(source)
    n = 1 + k + length
    tail = ((1 << length) - 1) << (1 + k)
    rows = [2 | (contacts << (1 + k))]
    rows.extend((row << 1) | tail | int(i != 0) for i, row in enumerate(source))
    rows.extend((1 << n) - (2 << (1 + k + j)) + int(not (contacts >> j & 1))
                for j in range(length))
    return tuple(rows)


def reader_extensions(source, length):
    """All new strong parents with this rooted deletion and this many tail vertices."""
    if length == 0 and len(source) == 1:
        return
    for contacts in range(1 << max(0, length - 1)):
        yield extend_reader(source, length, contacts)


def strong_readers(limit):
    """Root is vertex 0. The source component after deleting it has root 1."""
    readers = {1: [(0,)]}
    yield (0,)
    for n in range(3, limit + 1):
        layer = []
        for k, smaller in readers.items():
            length = n - 1 - k
            for source in smaller:
                layer.extend(reader_extensions(source, length))
        readers[n] = layer
        yield from layer


def graft_chart(reader, data=None):
    g, slots, starts, ends, cap = profile(reader, 1) if data is None else data
    weight = 25 * g
    slope = 28 * slots + 25 * sum(starts)
    intercept = 43 * g + 25 * (cap - g) + 28 * sum(ends) - 10 * (ends[0] - ends[1])
    return slope, weight, intercept, intercept - 25 * (starts[0] - starts[1])


def ceildiv(a, b):
    return -((-a) // b)


def carry_candidates(exponent, width, chart, bit):
    """Yield each nonnegative quotient whose residual could be an integral carry."""
    a = 1 << (exponent - 1)
    b = 1 << width
    assert 1 <= width < exponent
    center = a * (2 * a - 1)
    c, h, k0, k1 = chart
    upper = ((center - 1 - (k1 if bit else k0)) // c - bit * (b - 1)) // b
    if upper < 0:
        return
    if bit:
        slope, constant = c, k1 + 1 - c
        lo, hi = constant, constant + h * upper
    else:
        slope, constant = c + h, k0 + 1
        lo, hi = constant - h * upper, constant
    carry = ceildiv(lo, b)
    carry += (center // b - carry) % slope
    for q in range(carry, hi // b + 1, slope):
        z = (center // b - q) // slope - bit
        if 0 <= z <= upper:
            yield z


def test_carry(exponent, width, chart):
    c, h, k0, k1 = chart
    b = 1 << width
    a = 1 << (exponent - 1)
    center = a * (2 * a - 1)
    tests = 0
    exact = []
    for bit in (0, 1):
        for z in carry_candidates(exponent, width, chart, bit):
            tests += 1
            if bit:
                error = k1 + 1 - c + h * p(z) + c * b * (z + 1) - center
            else:
                error = k0 + 1 + h * (p(z) - z + (z & 1)) + (c + h) * b * z - center
            if not error:
                exact.append((1 << width) * z + bit * ((1 << width) - 1))
    return tests, exact


def three_runs(exponent, chart):
    """Invert the leading one-run in all words 1^a 0^b 1^d; allow d=0."""
    a = 1 << (exponent - 1)
    target = a * (2 * a - 1) - 1
    c, h, k0, k1 = chart
    quotient, remainder = divmod(target - k1, c)
    if remainder == 0 and quotient > 0 and not ((quotient + 1) & quotient):
        yield quotient
    max_width = (target // c).bit_length()
    for b in range(1, max_width):
        zeros = 1 << b
        for d in range(max_width - b):
            ones = 1 << d
            k = k1 if d else k0
            dividend = target - k - c * (ones - 1) - h
            divisor = c * zeros * ones + h * (zeros - 1)
            prefix, remainder = divmod(dividend, divisor)
            if not remainder and prefix > 0 and not (prefix & (prefix + 1)):
                yield prefix * zeros * ones + ones - 1


def envelope(exponent, n, suffix_bits=0):
    a = 1 << (exponent - 1)
    qmax = factorial(n) * sum(Fraction(1, 1 << (i // 2 + (n - i) // 2)) for i in range(n + 1))
    cmax, kmax = 53 * qmax, (78 + Fraction(28, n + 1)) * qmax
    ratio = Fraction(25, 28 * (n + 1))
    hmax = 25 * qmax / (n + 1)
    scale = 1 << suffix_bits
    if suffix_bits:
        cmax, kmax = (scale * (cmax + hmax),
                       kmax + (scale - 1) * (cmax + hmax) + hmax)
    return cmax + hmax < 2*a - 1 and max(cmax, kmax + 1) + 2*a*ratio/scale < a


def verify_nearby():
    center = (1 << 60) * ((1 << 61) - 1)
    deletion = ordinal(CORE, (0,))
    tags = 1 | (1 << len(CORE))
    deletions = deletion, flipped(deletion, 1, 3)
    reads = [profile(d, tags) for d in deletions]
    addresses = (
        (29120694641716595143708981840017769, -2735453618),
        (30980385641069175129382758141803891, 380266414),
    )
    for y, shift in addresses:
        low, high = (count_from_profile(read, 0, y) for read in reads)
        assert (low, high) == (center + shift - 1, center + shift + 1)
        a, b = (star_rows(d, tags, 0, y) for d in deletions)
        assert all(map(strong, (a, b))) and flipped(a, 1, 3) == b
        print("order", len(a), "center offset", shift, "endpoints", low, high)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=9)
    parser.add_argument("--width", type=int, default=48)
    args = parser.parse_args()
    verify_nearby()
    assert all(envelope(61, n) for n in range(1, 20))
    for n, last in ((1, 52), (5, 46), (8, 39), (9, 36), (10, 33),
                    (12, 27), (14, 20), (16, 13), (18, 5), (19, 1)):
        assert all(envelope(61, n, bits) for bits in range(last + 1))
    readers = tests = hits = few_runs = 0
    for reader in strong_readers(args.order):
        readers += 1
        chart = graft_chart(reader)
        count, exact = test_carry(61, args.width, chart)
        simple = list(three_runs(61, chart))
        few_runs += len(simple)
        if simple:
            print("exact address with at most three runs", reader, simple)
        tests += count
        hits += len(exact)
        if exact:
            print("exact receiver", reader, "addresses", exact)
    print(readers, "rooted receivers;", tests, "carry candidates;", hits, "hits;", few_runs, "hits with at most three runs")


if __name__ == "__main__":
    main()
