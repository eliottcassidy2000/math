#!/usr/bin/env python3
"""Exact F64 model of the 7/21 projective branch-law bridge."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json


def require(ok, message):
    if not ok:
        raise ArithmeticError(message)


def mul(x, y):
    z = 0
    while y:
        if y & 1:
            z ^= x
        y >>= 1
        x <<= 1
        if x & 64:
            x ^= 67  # X^6+X+1
    return z


def power(x, n):
    out = 1
    while n:
        if n & 1:
            out = mul(out, x)
        x = mul(x, x)
        n //= 2
    return out


def main():
    require({power(2, n) for n in range(63)} == set(range(1, 64)), 'primitive field model')
    f4 = {x for x in range(64) if power(x, 4) == x}
    f8 = {x for x in range(64) if power(x, 8) == x}
    require(len(f4) == 4 and len(f8) == 8 and f4 & f8 == {0, 1}, 'subfields')
    canon = lambda x: min(mul(x, s) for s in f4-{0})
    points = {canon(x) for x in range(1, 64)}
    require(len(points) == 21, 'projective points')
    blocks = {phase: {x for x in points if power(x, 21) == phase} for phase in f4-{0}}
    require(all(len(v) == 7 for v in blocks.values()), 'norm-phase block sizes')
    require(blocks[1] == {canon(x) for x in f8-{0}}, 'base Fano plane')
    for phase, block in blocks.items():
        u = next(x for x in range(1, 64) if power(x, 21) == phase)
        require(block == {canon(mul(u, x)) for x in f8-{0}}, 'multiplicative Fano translate')
        lines = set()
        for p in block:
            for q in block-{p}:
                line = {canon(mul(a, p) ^ mul(b, q)) for a in f4 for b in f4
                        if mul(a, p) ^ mul(b, q)}
                require(len(line) == 5, 'ambient F4 line')
                cut = frozenset(line & block)
                require(len(cut) == 3, 'Fano subplane line')
                lines.add(cut)
        require(len(lines) == 7, 'seven Fano lines')
    laws = []
    for q in (2, 4):
        counts = [q*q*(q-1), q*(q-1), q-1]
        require(sum(counts) == q**3-1, 'first nonzero coordinate census')
        require(F(counts[2], sum(counts)) == F(1, q*q+q+1), 'branch point probability')
        laws.append(dict(q=q, projective_size=q*q+q+1,
                         first_nonzero_position_probabilities=[str(F(c, sum(counts))) for c in counts]))
    require(F(1, 21) == F(1, 3)*F(1, 7), 'norm phase conditioning')
    print(json.dumps(dict(status='PASS; finite-field and distribution bridge only',
                          source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                          f4=sorted(f4), f8=sorted(f8),
                          norm_phase_blocks={str(k): sorted(v) for k, v in blocks.items()},
                          laws=laws), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
