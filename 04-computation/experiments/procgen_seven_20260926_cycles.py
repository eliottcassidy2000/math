#!/usr/bin/env python3
"""
procgen_seven_20260926_cycles.py -- rational sign-choice cycles of small shape (lane "seven", session
collatz-procgen-20260922, 2026-09-26).

T_eps(x) = x/2 (x even), (q x + eps)/2 (x odd), eps in {+1, -1}, on rationals with odd denominator (the
parity of a/D, D odd, is the parity of a).  A parity word w of length p with a ones and a sign word of
length a determine the periodic point x0 = c/(2^p - q^a), c = sum over the odd positions j (in order) of
eps * q^(number of later odd positions) * 2^j.  It is a genuine cycle iff the parities of its orbit agree
with w; the cycle is simple iff its p points are distinct.

A level-k sign strategy sigma contains such a cycle (as a closed walk of its parity graph G_sigma) iff sigma
agrees with the required sign at the residue mod 2^k of every odd point of the cycle.
"""
import itertools
from fractions import Fraction as Fr


def shape_cycles(q, p, a):
    """all simple cycles of shape (p,a) (one entry per cycle): list of (orbit, word, signs at odd points)"""
    out = {}
    for pos in itertools.combinations(range(p), a):
        w = [0] * p
        for i in pos:
            w[i] = 1
        for sg in itertools.product((1, -1), repeat=a):
            c, oa, si = 0, a, 0
            for j, bt in enumerate(w):
                if bt:
                    oa -= 1
                    c += sg[si] * q ** oa * 2 ** j
                    si += 1
            x0 = Fr(c, 2 ** p - q ** a)
            y, ok, si, orb = x0, True, 0, []
            for bt in w:
                orb.append(y)
                if y.numerator % 2 != bt:
                    ok = False
                    break
                if bt:
                    y = (q * y + sg[si]) / 2
                    si += 1
                else:
                    y = y / 2
            if not ok or y != x0 or len(set(orb)) != p:
                continue
            key = frozenset(orb)
            if key not in out:
                out[key] = (orb, w, list(sg))
    return list(out.values())


def residue(x, N):
    return (x.numerator * pow(x.denominator, -1, N)) % N


def present(k, sigbit, cyc):
    """sigbit(i) = 1 iff the sign of the odd node 2i+1 is '-'; True iff the level-k strategy contains cyc"""
    N = 1 << k
    orb, w, sg = cyc
    si = 0
    for y, bt in zip(orb, w):
        if bt:
            r = residue(y, N)
            s = -1 if sigbit((r - 1) // 2) else 1
            if s != sg[si]:
                return False
            si += 1
    return True


def sgn_keeps(cyc):
    """does the real-sign assignment sigma(x) = sgn(x) (the 'outward' sign) contain this cycle?"""
    orb, w, sg = cyc
    si = 0
    for y, bt in zip(orb, w):
        if bt:
            if (1 if y > 0 else -1) != sg[si]:
                return False
            si += 1
    return True
