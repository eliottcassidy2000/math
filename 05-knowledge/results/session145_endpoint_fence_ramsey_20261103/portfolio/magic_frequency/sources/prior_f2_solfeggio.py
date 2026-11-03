#!/usr/bin/env python3
"""Audit F2, check (7): the Solfeggio merge (section 3.4 of the results note) and the 'denominator 7' claim."""
from fractions import Fraction
from math import comb
from sympy import primerange, factorint


def gib(a, b, n):
    g = [a, b]
    while len(g) < n: g.append(g[-1] + g[-2])
    return g


F = gib(0, 1, 400); G51 = gib(5, 1, 400); G71 = gib(7, 1, 400)
print('(7,1) = (F_9, F_10) mod 9:', (F[9] % 9, F[10] % 9) == (7, 1), '; G_n = F_(n+9) mod 9 for n < 380:', all((G71[n] - F[n + 9]) % 9 == 0 for n in range(380)))
print('(5,1): G_n = 4 F_(n+14) mod 9 for n < 380:', all((G51[n] - 4 * F[n + 14]) % 9 == 0 for n in range(380)),
      '; also G_n = 5 F_(n+2) = -4 F_(n+2) mod 9:', all((G51[n] - 5 * F[n + 2]) % 9 == 0 for n in range(380)))
# which unit multiples c F_(n+k) represent (5,1) mod 9
reps = sorted({(c, k) for c in range(1, 9) for k in range(24) if all((G51[n] - c * F[n + k]) % 9 == 0 for n in range(48))})
print('   all representations G = c F_(n+k) mod 9 (c unit, k mod 24):', reps)
dr = lambda m: m % 9 or 9
print('(7,1) digits 3,6,9 at n mod 4:', sorted({n % 4 for n in range(240) if dr(G71[n]) in (3, 6, 9)}),
      '; (5,1):', sorted({n % 4 for n in range(240) if dr(G51[n]) in (3, 6, 9)}))
print('(5,1) 24-cycle of digital roots:', [dr(G51[n]) for n in range(24)])
cls = lambda N: {1: '+-1', 8: '+-1', 2: '+-2', 7: '+-2', 4: '+-4', 5: '+-4'}.get(N % 9, 'non-primitive')
print('golden class (norm b^2-ab-a^2 mod 9 up to sign): Fib', cls(-1), '; (5,1)', cls(-29), '; (7,1)', cls(-55), '; 4*Fib', cls(16 * -1))
print('x -> 4x fixes each fibre of Z/9 -> Z/3:', all((4 * x - x) % 3 == 0 for x in range(9)),
      '; norm class under pair -> 4 * pair (norm x16):', {c: cls(16 * n) for c, n in (('+-1', 1), ('+-2', 2), ('+-4', 4))})
# the ROOT constants: recompute 34/7 and 5/7 exactly from the THM-4690(4) tabulation
S_all = Fraction(4, 3)                                   # sum_{L>=0} 4^-L
S_not3 = Fraction(1, 3) - Fraction(1, 63)                # sum_{L>=1, 3 not | L} 4^-L
E0 = Fraction(3, 8) * (2 * 2 * S_all + 4 * 6 * S_not3)
SIB = Fraction(1, 8) * 3 * (2 + 4) * S_not3
print('E0/N =', E0, '; SIB/N =', SIB, '; sum_{3|L, L>=1} 4^-L = 1/(4^3-1) =', Fraction(1, 63))
print('   factors multiplying S_not3: E0: 3/8*24 =', Fraction(3, 8) * 24, '; SIB: 1/8*3*6 =', Fraction(1, 8) * 18,
      ' -> the 9 of 63 = 9*7 cancels against these weight factors')
# the period-3 dependence is a mod-3 fact: B(x) = 4x + 1 = x + 1 (mod 3); the weights depend on m mod 3 only
print('B(x) = x + 1 mod 3 for all x:', all((4 * x + 1 - x - 1) % 3 == 0 for x in range(3)),
      '; B^3(x) = 64x + 21 = x + 3 mod 9:', all((64 * x + 21 - x - 3) % 9 == 0 for x in range(9)))
# Solfeggio primes and numerology
SOLF = [174, 285, 396, 417, 528, 639, 741, 852, 963]
ps = sorted({p for s in SOLF for p in factorint(s) if p > 3})
typ = lambda p: 'ramified' if p == 5 else ('split' if p % 5 in (1, 4) else 'inert')
print('Solfeggio primes > 3:', ps, {p: typ(p) for p in ps})
print('  split or ramified:', sum(typ(p) != 'inert' for p in ps), 'of', len(ps),
      '; binomial tail P(>= 5 split of 8 | 1/2) = %.3f' % (sum(comb(8, k) for k in range(5, 9)) / 256))
pt51 = [g for g in G51[:20] if g < 150 and all(g % q for q in range(2, int(g ** 0.5) + 1)) and g > 1]
print('(5,1) prime terms below 150:', pt51, '; overlap with Solfeggio primes:', sorted(set(pt51) & set(ps)))


def hyper_tail(Npop, K, n, k):
    return sum(Fraction(comb(K, i) * comb(Npop - K, n - i), comb(Npop, n)) for i in range(k, min(K, n) + 1))
pr150 = list(primerange(2, 151))
print('primes <= 150: %d' % len(pr150))
print('  hypergeometric P(overlap >= 3 | pop 35, K = 8, n = 5) = %.4f' % float(hyper_tail(35, 8, 5, 3)))
print('  same with population = primes 5..150 (2, 3 excluded as in the Solfeggio list) = %.4f' % float(hyper_tail(33, 8, 5, 3)))
print('  excluding the seed G_0 = 5 (prime terms 7, 13, 53, 139; overlap {13, 139}): P(>= 2 | 35, 8, 4) = %.4f ; pop 33: %.4f'
      % (float(hyper_tail(35, 8, 4, 2)), float(hyper_tail(33, 8, 4, 2))))
