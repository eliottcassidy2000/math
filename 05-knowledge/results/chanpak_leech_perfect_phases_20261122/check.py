"""Construct a physical unit seam with prescribed mixed-precision center."""
import argparse
from functools import cache
from math import gcd, prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'chanpak_gap_star_ports_20261122'))
from replay import LEFT, RIGHT, count_from_profile, p

WHEEL = 4095
DEPTHS = {3: 2, 5: 1, 7: 1, 13: 1}
DESC = (8, 3, 6, 9, 7, 0, 4, 2, 5, 1, 10, 11)


def ray(x, s, r):
    block = 1 << s
    reverse = int(format(r, f'0{s}b')[::-1], 2)
    mass = 129 + 320 * x
    slope = block * mass + 29 * (block - 1 - reverse)
    intercept = 112 + 162*x + 29*p(x) + mass*r + 29*(p(r) + (2*r < block)) - slope
    return intercept, slope


@cache
def seeds():
    result = {}
    for s in range(1, 16):
        for r in range(0, 1 << s, 2):
            i, q = ray(2, s, r)
            if gcd(q, WHEEL) == 1:
                for m in range(12):
                    result.setdefault((i + (q << m)) % WHEEL, (s, r, m))
    assert len(result) == WHEEL
    return result


def high_packets(precision, junction_sum):
    """Increase x=2 by odd-neutral groups; the pair cost is junction_sum modulo 2^precision."""
    junctions = [j for j in range(1, precision) if junction_sum >> j & 1]
    if not junctions:
        return 2
    x = 2
    position = precision + 3
    ascending = DESC[::-1]
    gap = precision + 1
    for index in range(1 + len(junctions)):
        if index:
            position += junctions[index - 1] + 1
        x |= 1 << position
        for old, new in zip(ascending, ascending[1:]):
            position += gap + (new - old - gap) % 12
            x |= 1 << position
    assert x % WHEEL == 2 and (p(x) - p(2)) % WHEEL == 0
    assert (x - 2) % (1 << precision) == 0
    assert (p(x) - p(2) + junction_sum) % (1 << precision) == 0
    return x


def width(i, q, m0, center, powers, minimum):
    """Lift one exponent digit at a time, then combine the four lifted exponents."""
    t, period = 0, 1
    for prime, precision in zip(DEPTHS, powers):
        depth = DEPTHS[prime]
        local_t, local_period = 0, 1
        for level in range(depth + 1, precision + 1):
            modulus = prime**level
            matches = [local_t + d * local_period for d in range(prime)
                       if (i + q*pow(2, m0 + 12*(local_t + d*local_period), modulus)
                           - center) % modulus == 0]
            assert len(matches) == 1
            local_t = matches[0]
            local_period *= prime
        if local_period > 1:
            t += period * (((local_t - t)*pow(period, -1, local_period)) % local_period)
            period *= local_period
    m = m0 + 12*t
    if m < minimum:
        m += 12*period*((minimum - m + 12*period-1)//(12*period))
    return m


def decode(center, two_depth, powers):
    assert center % 2 == 0 and two_depth > 0 and len(powers) == 4 and min(powers) >= 0
    s, r, m0 = seeds()[center % WHEEL]
    i0, _ = ray(2, s, r)
    two_mod = 1 << two_depth
    correction = (i0 - center) * pow(29, -1, two_mod) % two_mod
    assert not correction % 2
    x = high_packets(two_depth, correction)
    i, q = ray(x, s, r)
    m = width(i, q, m0, center, powers, max(1, two_depth))
    mod = two_mod * prod(ell**a for ell, a in zip(DEPTHS, powers))
    assert (i + q * pow(2, m, mod) - center) % mod == 0
    return x, s, r, m, mod


def finite_checks():
    assert sorted(DESC) == list(range(12))
    norm = sum(1 << ((a-b-1)%12) for j, a in enumerate(DESC) for b in DESC[j+1:])
    assert norm == 3 * WHEEL
    assert len(seeds()) == WHEEL
    for e in range(1, 13):
        center = ((1 << e) - 1) << (e-1)
        assert ((center*center+1) % 5 == 0) == (e%4 == 3)
        assert ((center*center+1) % 13 == 0) == (e%12 == 11)
        assert (center*center+1) % 25 and (center*center+1) % 169
    for center in (0, 42, 496, 3583 + 4095, 4033 + 4095, 10**12):
        x,s,r,m,_ = decode(center, 4, (2,1,1,1))
        y = ((1 << m)-1)*(1 << s) + r
        low, high = count_from_profile(RIGHT,x,y), count_from_profile(LEFT,x,y)
        assert (low + 1 - center) % (16 * WHEEL) == 0 and high-low == 2
    perfect19 = ((1 << 19) - 1) << 18
    perfect31 = ((1 << 31) - 1) << 30
    slope = 137468394624
    assert (553528478 >> 17) == 4223 and (553528478 % (1 << 17)) == 11422
    for m in (0,12,24):
        y = (1 << 17)*((4224 << m) - 1) + 11422
        assert 292+173*y+143*p(y) == perfect19+slope*((1 << m)-1)
    for prime, valuation in zip(DEPTHS, (3,2,1,1)):
        assert (perfect31-perfect19) % (prime**valuation) == 0
        assert (perfect31-perfect19) % (prime**(valuation+1)) != 0
    print('Covered all 4095 residues, verified neutral packet, source phases, and physical controls.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--center', type=int, default=((1 << 31)-1) << 30)
    parser.add_argument('--two-depth', type=int, default=32)
    parser.add_argument('--odd-depths', nargs=4, type=int, default=(3,2,2,2))
    args = parser.parse_args()
    finite_checks()
    x,s,r,m,mod = decode(args.center, args.two_depth, args.odd_depths)
    print('Requested center and modulus:', args.center, mod)
    print('Left address bit length, suffix width, suffix, growing width:', x.bit_length(), s, r, m)
