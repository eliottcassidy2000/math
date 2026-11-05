#!/usr/bin/env python3
"""hc2_spectrum.py (klein-2026-11-05): the Hamiltonian-path counts h of every tournament with exactly one or exactly two
Hamiltonian cycles, from data/hval_{n}.txt ("hc h" per class, produced by scripts/hval.c from data/cls_{n}.txt).
Reports per order: min, max and number of distinct h-values of the uniquely Hamiltonian (hc = 1) and two-cycle (hc = 2)
families; the two-cycle floor f_2(n) against f(n) + 8 a(n) (f = Busch's strong floor 15, 25, 45, 75, ...; a = 1, 2, 3 then
a(n) = 5 a(n-3)); and the odd numbers below LIMIT missed by each family and by their union (cumulative over the orders).
Usage: hc2_spectrum.py N1 N2 [LIMIT]."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', 'data')
N1, N2 = int(sys.argv[1]), int(sys.argv[2]); LIMIT = int(sys.argv[3]) if len(sys.argv) > 3 else 2001
f = {3: 3, 4: 5, 5: 9, 6: 15, 7: 25, 8: 45}
for n in range(9, 40): f[n] = 5 * f[n - 3]
a = {6: 1, 7: 2, 8: 3}
for n in range(9, 40): a[n] = 5 * a[n - 3]
U = set(); D = set()
for n in range(N1, N2 + 1):
    u = set(); d = set()
    for l in open(os.path.join(DATA, 'hval_%d.txt' % n)):
        h, hv = map(int, l.split()); (u if h == 1 else d).add(hv)
    U |= u; D |= d
    print('order %2d: hc=1: min %5d (f(n) = %5d) max %5d distinct %4d | hc=2: min %5d (f(n) + 8a(n) = %5d) max %5d distinct %4d'
          % (n, min(u), f[n], max(u), len(u), min(d), f[n] + 8 * a[n], max(d), len(d)))
print('cumulative orders %d-%d: distinct h-values hc=1: %d, hc=2: %d, union: %d' % (N1, N2, len(U), len(D), len(U | D)))
for name, T in (('hc = 1', U), ('hc = 2', D), ('hc <= 2', U | D)):
    lo = min(T); miss = [x for x in range(lo, LIMIT, 2) if x not in T]
    print('  %-8s first value %d; odd numbers in [%d, %d) not attained: %d: %s' % (name, lo, lo, LIMIT, len(miss), miss[:60]))
