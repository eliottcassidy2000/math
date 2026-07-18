#!/usr/bin/env python3
"""
CF-DESCENT structure of the crux  (boxeph-2026-07-18-S103)

Owner: prove the medium-modulus inverse theorem via continued-fraction descent.
Finding: CF descent of the maximizer t* delivers a clean result -- the far element is
lcm(13,14)=182 because it must block BOTH the CF convergent 1/13 (needs 13|v) AND the
adjacent mediant 1/14 (needs 14|v) -- but reaches only the far element's divisibility,
NOT the AP core. Third elementary tool (after maximality S101, sieve S102) to stop at the
same 13/14-divisibility layer. The AP core is the open inverse theorem (= LRC(14), S94).
"""
from math import gcd
from fractions import Fraction as Fr


def Mstar(V, QMAX=2500):
    best = (Fr(0), 1, 0)
    for q in range(2, QMAX + 1):
        for a in range(1, q):
            if gcd(a, q) != 1:
                continue
            m = min(min((a * v) % q, q - ((a * v) % q)) for v in V)
            if Fr(m, q) > best[0]:
                best = (Fr(m, q), q, a)
    return best


def cf(a, b):
    out = []
    while b:
        out.append(a // b)
        a, b = b, a % b
    return out


def convergents(digs):
    ps, qs, conv = [0, 1], [1, 0], []
    for d in digs:
        ps.append(d * ps[-1] + ps[-2])
        qs.append(d * qs[-1] + qs[-2])
        conv.append((ps[-1], qs[-1]))
    return conv[1:]


print('CF-descent of the deep-well tower maximizer:')
for d in [1, 2, 3]:
    V = [d * i for i in range(1, 13)] + [d * 182]
    M, q, a = Mstar(V)
    digs = cf(a, q)
    print(f'  d={d}: M={M}={float(M):.5f}  t*={a}/{q}=CF{digs}  convergents={convergents(digs)}')

print()
print('Primitive deep well {1..12,182}: t*=14/183=[0;13,14]. The two best rationals flanking t*:')
V = list(range(1, 13)) + [182]
for t, name in [(Fr(1, 13), 'CF convergent 1/13 (needs 13|v)'),
                (Fr(1, 14), 'mediant 1/14 (needs 14|v)')]:
    m = min((min((t.numerator * v) % t.denominator,
                 t.denominator - ((t.numerator * v) % t.denominator)), v) for v in V)
    print(f'  t={t}: min_v||v t|| = {m[0]}/{t.denominator}, killed by v={m[1]} '
          f'({t.denominator}|{m[1]}={m[1] % t.denominator == 0})   [{name}]')
print()
print('=> far element = lcm(13,14) = 182 blocks BOTH. CF descent EXPLAINS 182, but delivers')
print('   only far-element divisibility -- the AP core {1..12} is the open inverse theorem.')
