#!/usr/bin/env python3
"""Tables of the mod-4 race from the outputs of scrace2.c (out/scrace2_n.txt, n = 3..10).
Each line: sc strong sym r classes labelled   (sc = self-converse, strong, sym = |Aut| > 1, r = H mod 8).
Prints: the class-level race, the labelled race split into self-converse / symmetric non-self-converse / rigid,
the strong-only race (checked against the free-monoid inversion of the signed counts), mod-8 masses, and the number
SC_n of labelled self-converse tournaments against the involution bound."""
from fractions import Fraction
from math import comb
from collections import defaultdict
import os
here = os.path.dirname(os.path.abspath(__file__))
def load(n):
    d = {}
    for line in open(os.path.join(here, 'out', 'scrace2_%d.txt' % n)):
        a, b, c, r, cl, l = map(int, line.split()); d[(a, b, c, r)] = (cl, l)
    return d
def sgn(r): return 1 if r % 4 == 1 else -1
def dfact(k):
    v = 1
    while k > 1: v *= k; k -= 2
    return v
tchi = {0: 1, 1: 1, 2: 2}; rows = {}
for n in range(3, 11):
    d = load(n); tot = 2 ** comb(n, 2)
    assert sum(v[1] for v in d.values()) == tot
    def S(sel, idx, signed=True):
        return sum((sgn(k[3]) if signed else 1) * v[idx] for k, v in d.items() if sel(*k))
    rows[n] = dict(
        cls1=S(lambda a, b, c, r: r % 4 == 1, 0, False), cls3=S(lambda a, b, c, r: r % 4 == 3, 0, False),
        sc1=S(lambda a, b, c, r: a and r % 4 == 1, 0, False), sc3=S(lambda a, b, c, r: a and r % 4 == 3, 0, False),
        st1=S(lambda a, b, c, r: b and r % 4 == 1, 0, False), st3=S(lambda a, b, c, r: b and r % 4 == 3, 0, False),
        L=S(lambda *k: True, 1), Lsc=S(lambda a, b, c, r: a, 1), Lsym=S(lambda a, b, c, r: (not a) and c, 1),
        Lrig=S(lambda a, b, c, r: (not a) and (not c), 1), Lstr=S(lambda a, b, c, r: b, 1),
        Wsc=S(lambda a, b, c, r: a, 1, False), Wstr=S(lambda a, b, c, r: b, 1, False), Wsym=S(lambda a, b, c, r: (not a) and c, 1, False),
        m8=[S(lambda a, b, c, r, rr=rr: r == rr, 1, False) for rr in (1, 3, 5, 7)],
        sc8=[S(lambda a, b, c, r, rr=rr: a and r == rr, 0, False) for rr in (1, 3, 5, 7)], tot=tot)
    tchi[n] = rows[n]['L']
tchi[11] = 17075203538944            # THM-5127 (exact, from the order-10 vertex cubes)
print("1. Class-level race (isomorphism classes)")
print(" n | classes  =1 mod 4  =3 mod 4  lead | self-converse (=1, =3)  others (=1, =3) | strong (=1, =3) lead")
for n in range(3, 11):
    r = rows[n]
    print("%2d | %8d %8d %8d %6d | %6d %6d   %8d %8d | %8d %8d %6d" % (n, r['cls1'] + r['cls3'], r['cls1'], r['cls3'], r['cls1'] - r['cls3'],
          r['sc1'], r['sc3'], r['cls1'] - r['sc1'], r['cls3'] - r['sc3'], r['st1'], r['st3'], r['st1'] - r['st3']))
print()
print("2. Labelled race: D_n = (self-converse) + (non-self-converse with |Aut| > 1) + (rigid, non-self-converse)")
print(" n | D_n           self-converse   symmetric       rigid         | times 4^n: D      sc      sym      rigid | sc share")
for n in range(3, 11):
    r = rows[n]; t = r['tot']
    print("%2d | %.9f  %+.9f  %+.9f  %+.9f | %8.1f %8.1f %8.1f %8.1f | %6.1f%%   exact: D = %s, sc = %s" % (n, r['L'] / t, r['Lsc'] / t, r['Lsym'] / t, r['Lrig'] / t,
          r['L'] / t * 4 ** n, r['Lsc'] / t * 4 ** n, r['Lsym'] / t * 4 ** n, r['Lrig'] / t * 4 ** n, 100 * r['Lsc'] / r['L'], Fraction(r['L'], t), Fraction(r['Lsc'], t)))
print()
print("3. Labelled self-converse tournaments SC_n against the involution bound B_n = (n-1)!! 2^(n^2/4) (n even), n!! 2^((n^2-1)/4) (n odd)")
for n in range(3, 11):
    r = rows[n]; B = dfact(n - 1) * 2 ** (n * n // 4) if n % 2 == 0 else dfact(n) * 2 ** ((n * n - 1) // 4)
    print("%2d | SC_n = %-12d B_n = %-12d SC_n/B_n = %.6f ; fraction of all tournaments %.9f ; signed sum over them %d" % (n, r['Wsc'], B, r['Wsc'] / B, r['Wsc'] / r['tot'], r['Lsc']))
print()
print("4. Strong-only race: s_chi(n) from the census flags, and from the inversion T_chi = 1/(1 - S_chi) of the signed counts")
N = 11; t = {n: 2 ** comb(n, 2) for n in range(N + 1)}
def invert(tt):
    s = {}
    for n in range(1, N + 1): s[n] = tt[n] - sum(comb(n, k) * s[k] * tt[n - k] for k in range(1, n))
    return s
s = invert(t); schi = invert(tchi)
for n in range(3, N + 1):
    direct = rows[n]['Lstr'] if n in rows else None
    assert direct is None or (direct == schi[n] and rows[n]['Wstr'] == s[n]), (n, direct, schi[n])
    print("%2d | strong labelled %-20d signed %-18d D_n^str = %.9f  (D_n = %.9f) ; decomposable share of the signed sum %.4f%s" % (
        n, s[n], schi[n], schi[n] / s[n], tchi[n] / t[n], (tchi[n] - schi[n]) / tchi[n], "" if direct is None else "  [census flags agree]"))
print()
print("5. Labelled relative masses of H mod 8 (residues 1, 3, 5, 7) and self-converse classes by H mod 8")
for n in range(3, 11):
    r = rows[n]
    print("%2d | %s | self-converse classes %s" % (n, "  ".join("%.5f" % (4 * x / r['tot']) for x in r['m8']), r['sc8']))
print("11 | 1.00052  0.99948  1.00044  0.99956   (THM-5127 audit 2: 0.25013, 0.24987, 0.25011, 0.24989)")
