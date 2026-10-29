#!/usr/bin/env python3
"""Independent integer/native checks; no imports from incoming producers.
Incoming: atomic_tail_reset_20261029/ATOMIC.md @
201fa96604d51171909318cc482220e00df54c38;
seven_return_native_recursions_20261029/COLLATZ.md @
63ef3ca1172665f7100d06d405a3aa34c42efa57.
"""
import json
import hashlib
import sys
from pathlib import Path
B = ((1, 1, 0, 1, 0, (1,)), (4, 1, 1, 6, 63, (3, 1, 1, 1)), (6, 3, 11, 11, 679, (3, 3, 1, 1, 1, 1, 1)), (5, 9, 25, 11, 1135, (3, 1, 3, 1, 1, 1, 1)), (4, 27, 59, 11, 1439, (3, 1, 1, 3, 1, 1, 1)))

def val(n, p):
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v

def forward(x, i):
    e, a, b, d, r, word = B[i]
    if not (x > 0 and x % (1 << d) == r):
        raise RuntimeError('Explicit invariant check failed: x > 0 and x % (1 << d) == r')
    numerator = 3 ** e * (a * x + b)
    if not numerator % (1 << d) == 0:
        raise RuntimeError('Explicit invariant check failed: numerator % (1 << d) == 0')
    y = numerator >> d
    n = 2 * x - 1
    for expected in word:
        t = 3 * n + 1
        if not val(t, 2) == expected:
            raise RuntimeError('Explicit invariant check failed: val(t, 2) == expected')
        n = t >> expected
    if not (n == 2 * y - 1 and y > x):
        raise RuntimeError('Explicit invariant check failed: n == 2 * y - 1 and y > x')
    return y

def predecessors(y):
    out = []
    for i, (e, a, b, d, r, word) in enumerate(B):
        numerator = (y << d) - 3 ** e * b
        denominator = 3 ** e * a
        if numerator <= 0 or numerator % denominator:
            continue
        x = numerator // denominator
        if x % 81 or x % (1 << d) != r:
            continue
        if not forward(x, i) == y:
            raise RuntimeError('Explicit invariant check failed: forward(x, i) == y')
        out.append((x, i))
    return out

def levels(y, depth):
    row = [(y, (), 0, 0)]
    yield row
    for k in range(1, depth + 1):
        row = [(x, (i,) + w, j + len(B[i][5]), d + B[i][3]) for t, w, j, d in row for x, i in predecessors(t)]
        if not len(row) <= 3:
            raise RuntimeError((y, k, row))
        if not len({d - j for _, _, j, d in row}) == len(row):
            raise RuntimeError('Explicit invariant check failed: len({d - j for _, _, j, d in row}) == len(row)')
        if not all((5 * j - 3 * d == 2 * k for _, _, j, d in row)):
            raise RuntimeError('Explicit invariant check failed: all((5 * j - 3 * d == 2 * k for _, _, j, d in row))')
        yield row

def decode(y, J, D, marked=None):
    if (5 * J - 3 * D) % 2 or min(J, D) < 0:
        return None
    k = (5 * J - 3 * D) // 2
    if k < 0:
        return None
    matches = [r for r in list(levels(y, k))[-1] if r[2:] == (J, D)]
    if not len(matches) <= 1:
        raise RuntimeError('Explicit invariant check failed: len(matches) <= 1')
    if not matches:
        return None
    x, w, _, _ = matches[0]
    if marked is not None and marked != x:
        return None
    return (x, w)

def main():
    tested = 0
    for y in range(81, 810001, 81):
        for k, row in enumerate(levels(y, 12)):
            if row and y * 2048**k < 81 * 2187**k:
                raise RuntimeError('Reverse depth violates minimum-slope bound')
            tested += len(row)
    x, y, z = (85545639, 91351719, 97551864)
    if not (forward(x, 2) == y and forward(y, 2) == z):
        raise RuntimeError('Explicit invariant check failed: forward(x, 2) == y and forward(y, 2) == z')
    sharp = list(levels(z, 2))[-1]
    if not len(sharp) == 3:
        raise RuntimeError('Explicit invariant check failed: len(sharp) == 3')
    for source, word, J, D in sharp:
        if not decode(z, J, D, source) == (source, word):
            raise RuntimeError('Explicit invariant check failed: decode(z, J, D, source) == (source, word)')
    if not forward(4671, 1) == forward(5535, 4) == 5913:
        raise RuntimeError('Explicit invariant check failed: forward(4671, 1) == forward(5535, 4) == 5913')
    if not (4671 % 81 and 5535 % 81):
        raise RuntimeError('Explicit invariant check failed: 4671 % 81 and 5535 % 81')
    if not (forward(6399, 1) == 8100 and forward(8100, 0) == 12150):
        raise RuntimeError('Explicit invariant check failed: forward(6399, 1) == 8100 and forward(8100, 0) == 12150')
    if not decode(12150, 1, 1, 8100) == (8100, (0,)):
        raise RuntimeError('Explicit invariant check failed: decode(12150, 1, 1, 8100) == (8100, (0,))')
    if not decode(12150, 5, 7, 8100) is None:
        raise RuntimeError('Explicit invariant check failed: decode(12150, 5, 7, 8100) is None')
    if not decode(12150, 5, 7, 6399) == (6399, (1, 0)):
        raise RuntimeError('Explicit invariant check failed: decode(12150, 5, 7, 6399) == (6399, (1, 0))')
    for i in range(1, 5):
        e, a, b, d, r, word = B[i]
        for zeros in [0, 1, 2, 3, 4, 10, 30, 100, 257]:
            modulus = 1 << d + zeros
            x = 81 * (-b * pow(81 * a, -1, modulus) % modulus)
            if not x:
                x = 81 * modulus
            endpoint = forward(x, i)
            for _ in range(zeros):
                endpoint = forward(endpoint, 0)
            v = val(endpoint, 3)
            chi = pow(2, v, 81) * (endpoint // 3 ** v) % 81
            if not v == e + zeros:
                raise RuntimeError('Explicit invariant check failed: v == e + zeros')
            if not chi == {1: 61, 2: 13, 3: 70, 4: 34}[i]:
                raise RuntimeError('Explicit invariant check failed: chi == {1: 61, 2: 13, 3: 70, 4: 34}[i]')
            if not decode(endpoint, len(word) + zeros, d + zeros, x) == (x, (i,) + (0,) * zeros):
                raise RuntimeError('Explicit invariant check failed: decode(endpoint, len(word) + zeros, d + zeros, x) == (x, (i,) + (0,) * zeros)')
    receipt = {'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'status': 'PASS', 'endpoint_multiples81_checked': 10000, 'reverse_depth': 12, 'reverse_nodes_checked': tested, 'sharp_width3': sharp, 'long_zero_runs': [0, 1, 2, 3, 4, 10, 30, 100, 257], 'checks': ['all five native branches', '27Z collision', 'marked source hostile', 'width <=3', 'distinct clock defect', 'native-clock uniqueness']}
    path = Path(__file__).with_name('checks_optimized.json' if sys.flags.optimize else 'checks.json')
    path.write_text(json.dumps(receipt, indent=2) + '\n')
    print(json.dumps(receipt, indent=2))
if __name__ == '__main__':
    main()
