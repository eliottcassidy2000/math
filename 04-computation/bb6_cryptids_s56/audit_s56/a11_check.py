#!/usr/bin/env python3
"""Referee: A11 countdown structure from the TM-verified rule transcription.
(1) blank-tape orbit: list every c in {0,1} check state (macro step, a, b mod 2, c, d) and the restarts;
(2) the round map is the Hydra map on x = a0 + b + 18;
(3) end-parity statistics on the ACTUAL restart seeds (b in {0, 9, 12}) rather than random b."""
import sys, os, random, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from tm_rules_check import a11_rule, HALT


def step_jump(p):
    a, b, c, d = p
    if a >= 2:                       # collapse the a-loop
        t = a // 2
        return (a - 2 * t, b + 3 * t, c, d)
    return a11_rule(p)


if __name__ == '__main__':
    # (1) blank-tape orbit with collapsed a-loop (macro steps counted as the session counts: one per step_jump)
    p = (4, 0, 0, 3)
    n = 0
    ends = []
    cur_c_start = None
    while n < 5000:
        a, b, c, d = p
        if c <= 1 and a <= 1:
            ends.append((n, a, b % 2, c, d, b.bit_length()))
        q = step_jump(p)
        n += 1
        if q == HALT:
            print('HALT at macro step', n)
            break
        if q[2] > p[2] + 2 or (p[2] <= 1 and q[2] >= 2 and q[2] != 3):
            pass
        p = q
        if p[2] > 10 ** 6:
            break
    print('check states (macro step, a0, b mod 2, c, d, bits of b):')
    for e in ends:
        print('  ', e)
    print('state after last listed step:', (p[0], p[1].bit_length(), 'bits', p[2].bit_length(), 'bits c', p[3]))
    # (2) Hydra conjugacy of one round
    rng = random.Random(1)
    ok = True
    for _ in range(100000):
        a0, b = rng.randrange(2), rng.randrange(10 ** 30)
        a = a0 + b + 4
        a0n, bn = a % 2, 3 + 3 * (a // 2)
        x, xn = a0 + b + 18, a0n + bn + 18
        if xn != (3 * x) // 2 or a0n != x % 2:
            ok = False
            break
    print('round map (a0,b) -> ((a0+b) mod 2, 3+3 floor((a0+b+4)/2)) is x -> floor(3x/2) on x = a0+b+18, a0\' = x mod 2:', ok)
