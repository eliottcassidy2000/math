"""Construct the marked tournaments, reproduce local certificates, and steer gaps."""

import argparse
from math import isqrt


def p(address):
    total = 0
    while address:
        high, bit = divmod(address, 2)
        if bit == 0:
            total += high + (high & 1)
        address = high
    return total


def flipped(rows, i, j):
    out = list(rows)
    out[i] ^= 1 << j
    out[j] ^= 1 << i
    return tuple(out)


def star_rows(core, tags, x, y):
    """Ordinal left flank, marked core, right flank, and last vertex as apex."""
    left = bin(x)[2:] if x else ''
    right = bin(y)[2:][::-1].translate(str.maketrans('01', '10')) if y else ''
    tape = left + ''.join(str((tags >> i) & 1) for i in range(len(core))) + right
    n = len(tape)
    mask = (1 << n) - 1
    rows = [mask ^ ((1 << (i + 1)) - 1) for i in range(n)]
    width = len(core)
    core_mask = ((1 << width) - 1) << len(left)
    for i, row in enumerate(core):
        index = len(left) + i
        rows[index] = (rows[index] & ~core_mask) | (row << len(left))
    apex = 0
    for i, bit in enumerate(tape):
        if bit == '1':
            apex |= 1 << i
        else:
            rows[i] |= 1 << n
    return tuple(rows + [apex])


def subpaths(rows, starts=None):
    n = len(rows)
    size = 1 << n
    end = [[0] * n for _ in range(size)]
    for i in range(n):
        if starts is None or starts >> i & 1:
            end[1 << i][i] = 1
    for subset in range(1, size):
        for i, count in enumerate(end[subset]):
            if count:
                choices = rows[i] & ~subset
                while choices:
                    bit = choices & -choices
                    end[subset | bit][bit.bit_length() - 1] += count
                    choices ^= bit
    return [sum(row) for row in end]


def profile(core, tags):
    n = len(core)
    full = (1 << n) - 1
    paths = subpaths(core)
    paths[0] = 1
    h = paths[-1]
    q = sum(paths[s] * paths[full ^ s] for s in range(full + 1))

    def eo(graph, marked):
        rooted = subpaths(graph, marked)
        return tuple(sum(rooted[s] * paths[full ^ s] for s in range(1, full + 1)
                         if s.bit_count() % 2 == parity) for parity in (0, 1))

    reverse = tuple(row ^ full ^ (1 << i) for i, row in enumerate(core))
    k = subpaths(star_rows(core, tags, 0, 0))[-1]
    return h, q, eo(core, tags), eo(reverse, full ^ tags), k


def count_from_profile(data, x, y):
    h, q, (e, o), (end_e, end_o), k = data
    return (k + h * (p(x) + p(y)) + q * x * y + (end_e + end_o) * x +
            (e + o) * y - (end_e - end_o) * (x & 1) - (e - o) * (y & 1))


BASE = (4, 17, 34, 7, 13, 27)
OTHER = flipped(BASE, 2, 3)
LEFT = (29, 320, (59, 70), (73, 89), 113)
RIGHT = (29, 320, (58, 71), (72, 90), 111)


def cubic_count(a, b):
    assert a > 0 and b > 0 and (a + b) % 2 == 1
    return 1280 * a * b + 353 * a + 287 * b + 29 * (p(a) + p(b)) + 140


def steer(target):
    z = isqrt(target // 1280) + 1
    if z < 64:
        raise ValueError('the explicit bracket is guaranteed from z=64')
    q = 1
    while (2 * q) ** 3 <= z:
        q *= 2
    u = (2 * z + q + 1 + 2 * q - 1) // (2 * q)
    twice_center = 2 * q * u + q - 1
    final = isqrt(64 * z * q)
    if final * final < 64 * z * q:
        final += 1
    if final % 2 == 0:
        final += 1

    def address_at(index):
        twice_d = 2 * index + 1
        a = (twice_center + twice_d) // 2
        b = (twice_center - twice_d) // 2
        return cubic_count(a, b), a, b

    def crossing(upper_threshold):
        low = 0
        high = (final - 1) // 2
        assert address_at(low)[0] > target > address_at(high)[0]
        while high - low > 1:
            mid = (low + high) // 2
            if 2 * address_at(mid)[0] > upper_threshold:
                low = mid
            else:
                high = mid
        return address_at(low), address_at(high)

    above = crossing(2 * target + 1)[0]
    below = crossing(2 * target - 1)[1]
    assert below[0] < target < above[0]
    assert max(target - below[0], above[0] - target) ** 3 <= 11200 ** 3 * z ** 2
    return below, above


def inverse_one_sided(target, k, h, start_mass):
    """All even positive addresses giving this target, via two digit halves."""
    if target <= k:
        return []
    b = max(1, isqrt((target - k) // start_mass).bit_length())
    size = 1 << b
    result = []
    for r in range(0, size, 2):
        reverse = int(format(r, f'0{b}b')[::-1], 2)
        den = size * start_mass + h * (size - 1 - reverse)
        rem = target - k - start_mass * r - h * p(r)
        if rem < 0:
            continue
        first = max(0, (rem - h) // (den + h))
        for z in range(first, rem // den + 1):
            if den * z + h * (p(z) + (r < size // 2) * (z & 1)) == rem:
                w = size * z + r
                if w:
                    assert k + start_mass * w + h * p(w) == target
                    result.append(w)
    return sorted(result)


def prime(n):
    return n >= 2 and (n == 2 or (n % 2 and
           all(n % d for d in range(3, isqrt(n) + 1, 2))))


def strong(rows):
    full = (1 << len(rows)) - 1
    for root in range(len(rows)):
        seen = 1 << root
        while True:
            new = seen
            for i, row in enumerate(rows):
                if seen >> i & 1:
                    new |= row
            if new == seen:
                break
            seen = new
        if seen != full:
            return False
    return True


TARGETS = (
    (4294967387, (4, 5, 48, 7, 11, 27), (2, 3), 1, 113404326),
    (1000000000039, (18, 24, 67, 69, 44, 15, 51), (3, 5), 126, 712585682),
)


def checks():
    assert profile(BASE, 22) == LEFT and profile(OTHER, 22) == RIGHT
    controls = 0
    for core, data in ((BASE, LEFT), (OTHER, RIGHT)):
        for x in range(5):
            for y in range(5):
                rows = star_rows(core, 22, x, y)
                assert subpaths(rows)[-1] == count_from_profile(data, x, y)
                controls += 1
    print('literal core and both-flank controls:', controls)
    for N in (10**7, 2**32, 10**12, 10**24, 10**48, 2**128):
        below, above = steer(N)
        print('bracket', N, 'offsets', below[0] - N, above[0] - N,
              'a,b below/above', below[1:], above[1:])
        for val, a, b in (below, above):
            assert count_from_profile(RIGHT, 2*a, 2*b) == val
            assert count_from_profile(LEFT, 2*a, 2*b) == val + 2

    for target, core, endpoints, tag, w in TARGETS:
        high_core = flipped(core, *endpoints)
        lower, higher = profile(core, tag), profile(high_core, tag)
        k, h, c = lower[4], lower[0], sum(lower[2])
        assert count_from_profile(lower, 0, w) == target
        assert count_from_profile(higher, 0, w) == target + 2
        assert w in inverse_one_sided(target, k, h, c)
        a, b = star_rows(core, tag, 0, w), star_rows(high_core, tag, 0, w)
        assert strong(a) and strong(b)
        assert sum((x ^ y).bit_count() for x, y in zip(a, b)) == 2
        print('exact physical pair', target, target + 2, 'order', len(a),
              'chart', (k, h, c), 'address', w)
    twin = TARGETS[0][0]
    assert prime(twin) and prime(twin+2)
    assert all(not (prime(n) and prime(n+2)) for n in range(2**32+1, twin, 2))
    first = TARGETS[1][0]
    assert prime(first) and all(not prime(n) for n in range(10**12+1, first, 2))
    for x in range(1, 32, 2):
        for m in range(1, 12):
            target = (320*x+129)*2**m + 11-158*x+29*p(x)
            assert count_from_profile(LEFT, x, 2**m-1) == target
            assert count_from_profile(RIGHT, x, 2**m-1) == target+2


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--target', type=int, help='find two strict bracketing unit seams')
    group.add_argument('--witness', type=int, choices=[t[0] for t in TARGETS],
                       help='print literal outgoing masks for the requested exact pair')
    args = parser.parse_args()
    if args.target is not None:
        print(steer(args.target))
    elif args.witness:
        _, core, endpoints, tag, w = next(t for t in TARGETS if t[0] == args.witness)
        print(star_rows(core, tag, 0, w))
        print(star_rows(flipped(core, *endpoints), tag, 0, w))
    else:
        checks()
