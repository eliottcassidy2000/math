from itertools import product

TRI = (2, 4, 1)
F4 = (8, 1, 3, 6)


def ordinal(*ts):
    n = sum(map(len, ts))
    p = 0
    out = []
    for t in ts:
        later = ((1 << n) - 1) ^ ((1 << (p + len(t))) - 1)
        out.extend((row << p) | later for row in t)
        p += len(t)
    return tuple(out)


def extend(t, new_to_old):
    n = len(t)
    return tuple(row | (((new_to_old >> v & 1) ^ 1) << n)
                 for v, row in enumerate(t)) + (new_to_old,)


def count(t):
    n = len(t)
    full = (1 << n) - 1
    f = [[0] * n for _ in range(full + 1)]
    for i in range(n):
        f[1 << i][i] = 1
    for mask in range(1, full):
        remaining = full ^ mask
        for end, multiplicity in enumerate(f[mask]):
            if multiplicity:
                available = t[end] & remaining
                while available:
                    bit = available & -available
                    f[mask | bit][bit.bit_length() - 1] += multiplicity
                    available ^= bit
    return sum(f[full])


def check(t, base, r, cases):
    baseline = count(base)
    positive = set()
    for masks in cases:
        current = t
        for m in masks:
            current = extend(current, m)
        delta = count(current) - baseline
        if delta:
            if delta < 2 * 3 ** (r - 1):
                raise RuntimeError((masks, baseline, delta))
            if r >= 2 and 2 * 3 ** (r - 1) < delta < 4 * 3 ** (r - 1):
                raise RuntimeError(('second gap', masks, baseline, delta))
            positive.add(delta)
    return min(positive)


for r in (1, 2, 3):
    base = ordinal(*((TRI,) * r))
    d = check(base, base, r, ((m,) for m in range(1 << len(base))))
    if d != 2 * 3 ** (r - 1):
        raise RuntimeError((r, d))
    print(f'{r} triangle SCCs, all one-vertex extensions: min positive increment {d}')

base = ordinal(TRI, TRI)
d = check(base, base, 2, product(range(1 << 6), range(1 << 7)))
if d != 6:
    raise RuntimeError(d)
print(f'2 triangle SCCs, all 8192 two-vertex extensions: min positive increment {d}')

for r in (4,):
    lhs = ordinal(*((TRI,) * r))
    rhs = ordinal(F4, *((TRI,) * (r - 1)))
    if count(rhs) - count(lhs) != 2 * 3 ** (r - 1):
        raise RuntimeError(r)
print('r=4 sharp literal example passes.')
print('Second-gap exclusion passes in every declared extension universe.')
