#!/usr/bin/env python3
"""audit_F (A0): the pair-chain recursion of HYP-9244 (Setting) against actual integer orbits.

For random integers y (uniform mod d^N) iterate u = y + e0, v = y under T, and at every step compare
  (1) u's branch with i = (M mod d) j + (e mod d) mod d, j = v's branch;
  (2) M_n = prod m_{i_k}/prod m_{j_k} and e_n from the recursion e' = (m_i e + r_i - (m_i/m_j) M r_j)/d
      with the identity u_n = M_n v_n + e_n (exact rationals);
  (3) B_n e_n is an integer, where M_n = A_n/B_n in lowest terms (hence e_n is an integer whenever M_n = 1).
Usage: python3 chain_vs_orbits.py
"""
import random
from fractions import Fraction as Fr
import sys
sys.path.insert(0, '.')
from bigint_pairs import MAPS


def modd(q, d):
    return (q.numerator * pow(q.denominator, -1, d)) % d


def check(name, d, br, nsamp=60, steps=400, e0=1, seed=1):
    rnd = random.Random(seed)
    N = steps + 40
    bad1 = bad2 = bad3 = 0; tot = 0; ret = 0
    for _ in range(nsamp):
        y = d ** N + rnd.randrange(d ** N)
        u, v = y + e0, y
        M, e = Fr(1), Fr(e0)
        for n in range(steps):
            j = v % d; i = u % d
            if i != (modd(M, d) * j + modd(e, d)) % d:
                bad1 += 1
            mi, ri = br[i]; mj, rj = br[j]
            e = (mi * e + ri - Fr(mi, mj) * M * rj) / d
            M = M * Fr(mi, mj)
            u = (mi * u + ri) // d; v = (mj * v + rj) // d
            tot += 1
            if u != M * v + e:
                bad2 += 1
            if (M.denominator * e).denominator != 1:
                bad3 += 1
            if M == 1:
                ret += 1
            if u == v:
                break
    print(f"{name:9s}: {tot} steps; branch-rule mismatches {bad1}; identity u = M v + e failures {bad2}; "
          f"B_n e_n non-integer {bad3}; steps at M = 1: {ret}")


if __name__ == '__main__':
    for name in ('3x+1', '3x+5', 'Z3_124', 'Z3_125', 'Z5_12311', 'Z5_12471', 'Z5_12371', '5x+1', 'Z3_1416', 'Z3_157'):
        d, br = MAPS[name]
        check(name, d, br)
    print("offset e0 = 7 for Z3_125 and Z5_12371:")
    check('Z3_125', *MAPS['Z3_125'], e0=7)
    check('Z5_12371', *MAPS['Z5_12371'], e0=7)
