"""Construct one-site tournaments and check the positive-gap phase theorem."""

import argparse
from functools import lru_cache
from math import lcm


def one_apex(value):
    count = 1
    prefix = 0
    for bit in bin(value)[2:]:
        if bit == '0':
            count += prefix + (prefix & 1)
        prefix = 2 * prefix + int(bit)
    return count


def profile(u, w):
    common = 2 * u * w + one_apex(u) + one_apex(w) - 1
    return common + u + (u & 1), common + w + (w & 1)


def phases(modulus, d, prefix):
    if min(modulus, d, prefix) < 1:
        raise ValueError("modulus, half-gap, and prefix must all be positive")
    period = lcm(modulus, 4)
    odd_part = period // (period & -period)
    ell = (period & -period).bit_length() - 1
    while pow(2, ell, odd_part) != 1 % odd_part:
        ell += 1
    nblocks = period // 4 - 1
    seed = ((prefix + 1) << ell) - 1
    switch_increment = ((1 << ell) - 1) * seed + 1
    tail = 8 * d + 1
    shift = tail.bit_length()

    def addresses(switches, odd_branch):
        q = seed
        for _ in range(switches):
            q = (q << (2 * ell)) + (1 << ell) - 1
        free_ones = 2 * ell * (nblocks - switches)
        q = ((q + 1) << free_ones) - 1
        s = (q << shift) + tail
        return 2 * (s - d), 2 * s - odd_branch

    bases = [profile(*addresses(0, branch))[0] for branch in (0, 1)]

    def decode(residue):
        residue %= modulus
        if residue % 2 == 0:
            if modulus % 2 == 0:
                raise ValueError("no odd integer has the requested residue")
            residue += modulus
        branch = next(b for b in (0, 1) if (residue - bases[b]) % 4 == 0)
        step = switch_increment // 2
        reduced_modulus = period // 4
        switches = 0
        if reduced_modulus > 1:
            inverse = pow(step, -1, reduced_modulus)
            switches = ((residue - bases[branch]) // 4 * inverse) % reduced_modulus
        return addresses(switches, branch)

    return decode


@lru_cache(maxsize=None)
def literal_count(word):
    # Vertex zero is the apex, vertices 1..n are the increasing spine.
    size = len(word) + 1
    predecessors = [0] * size
    for position, bit in enumerate(word, start=1):
        predecessors[position] = (1 << position) - 2
        if bit == '1':
            predecessors[position] |= 1
        else:
            predecessors[0] |= 1 << position

    @lru_cache(maxsize=None)
    def end_count(mask, end):
        rest = mask ^ (1 << end)
        if rest == 0:
            return 1
        choices = rest & predecessors[end]
        count = 0
        while choices:
            choice = choices & -choices
            count += end_count(rest, choice.bit_length() - 1)
            choices ^= choice
        return count

    return sum(end_count((1 << size) - 1, end) for end in range(size))


def words(u, w):
    left = bin(u)[2:]
    right = bin(w)[2:][::-1].translate(str.maketrans('01', '10'))
    return left + '0' + right, left + '1' + right


def finite_checks():
    literal_flips = 0
    for size in range(4, 10):
        for left_width in range(1, size - 2):
            right_width = size - 2 - left_width
            for u in range(1 << (left_width - 1), 1 << left_width):
                for w in range(1 << (right_width - 1), 1 << right_width):
                    assert tuple(map(literal_count, words(u, w))) == profile(u, w)
                    literal_flips += 1
    print(f"literal marked flips, orders 4..9: {literal_flips}")

    moduli = (1, 2, 3, 4, 5, 6, 8, 9, 12, 15, 16, 20, 21, 24, 28, 32, 40, 60)
    checks = 0
    for modulus in moduli:
        for d in (1, 2, 3, 7):
            for prefix in (1, 5, 12):
                decode = phases(modulus, d, prefix)
                widths = set()
                for residue in range(modulus):
                    if modulus % 2 == residue % 2 == 0:
                        continue
                    u, w = decode(residue)
                    lower, upper = profile(u, w)
                    assert lower % modulus == residue and upper == lower + 2 * d
                    assert u.bit_length() == w.bit_length()
                    assert bin(u).startswith(bin(prefix))
                    widths.add(u.bit_length())
                    checks += 1
                assert len(widths) == 1
    print(f"decoded phase/gap/prefix triples: {checks}; moduli {moduli}")

    boundary_checks = 0
    for k in range(3, 8):
        t = 1 << k
        for d in range(1, 1 + (t // 2 - 1) // 2):
            pairs = [(u, w) for u in range(t // 2, t) for w in range(t, 2 * t)
                     if (w + 1) // 2 - (u + 1) // 2 == d]
            assert len(pairs) == 4 * d
            if d == 1:
                actual = sorted(profile(u, w)[0] + 1 for u, w in pairs)
                offsets = (-15 * t // 4, -3 * t // 2, 3 * t // 2, 15 * t // 4)
                assert actual == [2 * t * t + offset for offset in offsets]
            boundary_checks += 1
    print(f"odd-order full adjacent-width scans: {boundary_checks}; k=3..7")

    for j in range(1, 8):
        modulus = 3 ** j
        period = 2 * 3 ** (j - 1)
        for numerator, denominator in ((-15, 4), (-3, 2), (3, 2), (15, 4)):
            coefficient = numerator * pow(denominator, -1, modulus)
            solutions = []
            for k in range(period):
                t = pow(2, k, modulus)
                if (2 * t * t + coefficient * t + 1) % modulus == 0:
                    solutions.append(k)
            assert len(solutions) == 2 and solutions[0] % 2 != solutions[1] % 2
    print("3-adic divisibility: two width classes per branch through 3^7")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--decode', type=int, nargs=4,
                        metavar=('MODULUS', 'RESIDUE', 'HALF_GAP', 'PREFIX'))
    args = parser.parse_args()
    if args.decode is None:
        finite_checks()
        return
    modulus, residue, d, prefix = args.decode
    u, w = phases(modulus, d, prefix)(residue)
    lower, upper = profile(u, w)
    assert lower % modulus == residue % modulus and upper - lower == 2 * d
    print("order", u.bit_length() + w.bit_length() + 2)
    print("flank addresses (hex)", hex(u), hex(w))
    print("lower and upper path counts (hex)", hex(lower), hex(upper))


if __name__ == '__main__':
    main()
