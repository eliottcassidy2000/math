"""Decode an odd target on four two-run families and reconstruct its marked graph."""
from argparse import ArgumentParser
from dataclasses import dataclass
from math import isqrt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "chanpak_decimal_kernel_20261122"))
from minimum_eleven import CERTIFICATES as D11_CERTIFICATES
from e15_pair import flipped, join, ordinal, literal_count, p, star_rows, strong, tagged_data


@dataclass(frozen=True)
class Chart:
    discriminant: int
    cube_root: int
    lower: tuple
    tags: int
    edge: tuple
    prefix: int
    tail: int
    first: int
    second: int

    receiver: tuple = ()
    receiver_mask: int = 0

    def blocks(self, upper):
        base = flipped(self.lower, *self.edge) if upper else self.lower
        return (base, self.receiver) if self.receiver else (base,)

    def profile(self, upper):
        blocks = self.blocks(upper)
        base = tagged_data(blocks[0], self.tags)
        return join(base, tagged_data(blocks[1], self.receiver_mask)) if self.receiver else base

    def address(self, a, b):
        u = ((self.prefix + 1) << a) - 1
        return ((u << (b + self.tail)) | ((1 << self.tail) - 1))

    def parameters(self):
        parity = int(self.tail > 0)
        h, _, c, _, e, _, k = self.profile(False)
        H, _, C, _, E, _, K = self.profile(True)
        assert (h, c, k + e*parity + 2) == (H, C, K + E*parity)
        middle = k + e*parity + 1 + c*((1 << self.tail) - 1) + h*(p(self.prefix) + 1)
        return h, (c << self.tail) + h, self.prefix + 1, middle


CHARTS = (
    Chart(19, 96, (6, 68, 24, 19, 35, 143, 61, 95), 40, (1, 2), 2, 0, 3, 5),
    Chart(43, 960, (4, 5, 56, 3, 11, 27, 63), 80, (2, 3), 152418, 1, 2, 1),
    Chart(67, 5280, (72, 37, 97, 22, 71, 153, 170, 31), 5, (0, 2), 162, 0, 1, 1,
          (1126, 4092, 4088, 4081, 4065, 4032, 3968, 3841, 3585, 3073, 2048, 1), 1),
    Chart(163, 640320, (6, 32, 10, 35, 15, 85, 31), 11, (0, 1), 4245148004708, 5, 1, 2),
)
MODULUS = 26896077
PERIOD = 36


def midpoint(params, a, b):
    h, s, prefix, constant = params
    return constant + ((s << b) - h)*((prefix << a) - 1)


def phase_table(charts, modulus, period):
    lookup = {}
    for chart, params in charts:
        for a in range(period):
            for b in range(period):
                for sign in (-1, 1):
                    phase = (midpoint(params, a, b) + sign) % modulus
                    entry = chart.discriminant, a, b, sign
                    if phase in lookup:
                        return None
                    lookup[phase] = entry
    return lookup


def divisors(n):
    small = [k for k in range(1, isqrt(n)+1) if n % k == 0]
    return sorted(set(small) | {n//k for k in small})


def decode(target, charts, table):
    entry = table.get(target % MODULUS)
    if entry is None:
        return None
    d, a_phase, b_phase, sign = entry
    chart, params = next((c, v) for c, v in charts if c.discriminant == d)
    h, s, prefix, constant = params
    z = target - sign - constant
    if z <= 0:
        return None
    # The two-run length is the unique L with SP*2^(L-1) < z < SP*2^L.
    length = (z // (s * prefix)).bit_length()
    if length < 2 or length % PERIOD != (a_phase+b_phase) % PERIOD:
        return None
    scaled = 1 << length
    linear = s*prefix*scaled + h - z
    discriminant = linear*linear - 4*h*prefix*s*scaled
    if linear <= 0 or discriminant < 0:
        return None
    root = isqrt(discriminant)
    if root*root != discriminant:
        return None
    for numerator in {linear-root, linear+root}:
        power, remainder = divmod(numerator, 2*h*prefix)
        if remainder or power <= 1 or power & (power-1):
            continue
        a = power.bit_length()-1
        b = length-a
        if b < 1 or (a % PERIOD, b % PERIOD) != (a_phase, b_phase):
            continue
        if midpoint(params, a, b)+sign == target:
            return chart, a, b, sign
    return None


def verify_literal(chart, a, b, center):
    address = chart.address(a, b)
    graphs = []
    tags = chart.tags | (chart.receiver_mask << len(chart.lower))
    for upper in (False, True):
        parts = chart.blocks(upper)
        deletion = ordinal(*parts) if chart.receiver else parts[0]
        graphs.append(star_rows(deletion, tags, 0, address))
    assert flipped(graphs[0], *chart.edge) == graphs[1]
    blocks = tuple(map(len, chart.blocks(False))) + (1,)*address.bit_length()
    assert [literal_count(g, blocks) for g in graphs] == [center-1, center+1]
    assert all(strong(g) for g in graphs)
    return address, len(graphs[0])


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("target", type=int, nargs="?")
    parser.add_argument("--minimal", action="store_true", help="check divisors of 2^36-1 in order")
    args = parser.parse_args()
    charts = [(c, c.parameters()) for c in CHARTS]
    for _, (h,s,prefix,_) in charts:
        assert (2*s-h)*(2*prefix-1) > 2*s*prefix
    table = phase_table(charts, MODULUS, PERIOD)
    assert table is not None and len(table) == 10368 and pow(2, PERIOD, MODULUS) == 1
    isolated = {center+s for center, *_ in D11_CERTIFICATES for s in (-1,1)}
    assert len(isolated)==4 and all(t % MODULUS not in table for t in isolated)
    if args.minimal:
        earliest = next(m for m in divisors((1 << PERIOD)-1) if m >= len(table)
                        and phase_table(charts, m, PERIOD) is not None)
        assert earliest == MODULUS
        print("least period-36 separating modulus", earliest)
    if args.target is not None:
        for center, deletion, tags, vertex in D11_CERTIFICATES:
            if args.target in (center-1, center+1):
                actual_tags = tags if args.target < center else tags ^ (1 << vertex)
                print("discriminant 11, midpoint, apex mask, flip vertex", center, actual_tags, vertex)
                return
        found = decode(args.target, charts, table)
        if found:
            chart, a, b, sign = found
            # The decoder itself uses one modular lookup and one integer square root.
            print("discriminant, widths, sign, address", chart.discriminant, a, b, sign, chart.address(a,b))
        else:
            print("target outside the four families and two isolated edges")
        return
    for chart, params in charts:
        center = chart.cube_root**3 + 744
        assert midpoint(params, chart.first, chart.second) == center
        for a,b in [(chart.first,chart.second), (1,1), (1,74), (79,2), (900,1200)]:
            for sign in (-1,1):
                assert decode(midpoint(params,a,b)+sign, charts, table) == (chart,a,b,sign)
        y,n = verify_literal(chart, chart.first, chart.second, center)
        print("D, midpoint, address, vertices",chart.discriminant,center,y,n)
    for exponent in range(1, PERIOD+1):
        center = (1 << (exponent-1))*((1 << exponent)-1)
        assert all((center+s) % MODULUS not in table for s in (-1,1))
    print("10368 distinct endpoint phases; every Euclid-form neighbor misses this table")


if __name__ == "__main__":
    main()
