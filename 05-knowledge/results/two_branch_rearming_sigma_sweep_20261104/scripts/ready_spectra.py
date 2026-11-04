#!/usr/bin/env python3
"""ready_spectra.py (klein-2026-11-04) -- the exact digit-1 and digit-2 ready spectra from the generators.
Inputs: data/extgen_order{8..12}.txt ("d1 d2 H digraph6", exhaustive one-vertex extension of all ready classes, nauty
canonical forms) and data/gen1_order{11..N}.txt ("H digraph6", the pruned digit-1 generator).  Reports: class counts
against Garey's F(2n-6) (digit 1) and OEIS A377857(n-2) = (n-2) F(2n-9) - L(2n-10) (digit 2); the exact ready spectra
R1, R2 with their dyadic-window densities (a window [2^k, 2^(k+1)) is exact once every order whose floor is below
2^(k+1) is included); membership of the lineage's tree values; and a branching estimate for the two-branch grammar."""
import glob, os
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', 'data')
def fib(n):
    a, b = 0, 1
    for _ in range(n): a, b = b, a + b
    return a
def luc(n):
    a, b = 2, 1
    for _ in range(n): a, b = b, a + b
    return a
R1 = {}; R2 = {}; counts = {}
for f in sorted(glob.glob(os.path.join(DATA, 'extgen_order*.txt'))):
    n = int(f.split('order')[1].split('.')[0]); c1 = c2 = 0
    for l in open(f):
        d1, d2, H, c = l.split(); H = int(H)
        if d1 == '1': R1.setdefault(n, set()).add(H); c1 += 1
        if d2 == '1': R2.setdefault(n, set()).add(H); c2 += 1
    counts[n] = (c1, c2)
for f in sorted(glob.glob(os.path.join(DATA, 'gen1_order*.txt'))):
    n = int(f.split('order')[1].split('.')[0]); vals = set(); c1 = 0
    for l in open(f): vals.add(int(l.split()[0])); c1 += 1
    R1[n] = vals; counts[n] = (c1, counts.get(n, (None, None))[1])
for f in sorted(glob.glob(os.path.join(DATA, 'gen2_order*.txt'))):
    n = int(f.split('order')[1].split('.')[0]); vals = set(); c2 = 0
    for l in open(f): vals.add(int(l.split()[0])); c2 += 1
    R2[n] = vals; counts[n] = (counts.get(n, (None, None))[0], c2)
print('(1) class counts: order: digit-1 classes vs Garey F(2n-6); digit-2 classes vs A377857(n-2) = (n-2)F(2n-9) - L(2n-10)')
for n in sorted(counts):
    c1, c2 = counts[n]
    print('    order %2d: digit-1 %6s (F(2n-6) = %6d)  digit-2 %6s (A377857 = %6d)' % (n, c1, fib(2 * n - 6), c2 if c2 is not None else '-', (n - 2) * fib(2 * n - 9) - luc(2 * n - 10) if n >= 5 else 0))
allR1 = set().union(*R1.values()); allR2 = set().union(*R2.values())
fl1 = {n: min(v) for n, v in R1.items()}; fl2 = {n: min(v) for n, v in R2.items()}
print('(2) ready floors: digit 1:', dict(sorted(fl1.items())), '; digit 2:', dict(sorted(fl2.items())))
maxn1 = max(R1); maxn2 = max(R2)
nxt1 = 1.9 * fl1[maxn1]; nxt2 = 1.9 * fl2[maxn2]   # conservative estimate of the next order's floor (observed ratio 1.9-2.0 per vertex)
print('(3) dyadic-window densities of R1 (a window [a, 2a) is complete when 2a < floor1(%d) ~ %d) and R2 (2a < floor2(%d) ~ %d):' % (maxn1 + 1, nxt1, maxn2 + 1, nxt2))
for k in range(5, 15):
    a = 1 << k
    d1 = sum(1 for h in range(a + 1, 2 * a, 2) if h in allR1) / (a / 2); d2 = sum(1 for h in range(a + 1, 2 * a, 2) if h in allR2) / (a / 2)
    print('    [%6d, %6d): R1 %.3f%s   R2 %.3f%s' % (a, 2 * a, d1, '' if 2 * a < nxt1 else ' (incomplete)', d2, '' if 2 * a < nxt2 else ' (incomplete)'))
tree = [1429, 4289, 4291, 12869, 12871, 12875, 12877, 38609, 38611, 38615, 38617, 38627, 38629, 38633, 38635]
print('(4) tree values: in R1 (orders):', {t: sorted(n for n in R1 if t in R1[n]) for t in tree if t in allR1}, '; in R2:', {t: sorted(n for n in R2 if t in R2[n]) for t in tree if t in allR2})
print('    not in R1 within the computed orders:', [t for t in tree if t not in allR1], '; not in R2:', [t for t in tree if t not in allR2])
