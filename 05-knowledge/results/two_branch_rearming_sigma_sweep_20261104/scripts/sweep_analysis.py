#!/usr/bin/env python3
"""sweep_analysis.py (klein-2026-11-04) -- the C2 digit sweep from the ready pairs census.
Inputs: data/pairs_k.txt (k = 4..9; lines "H sigma" for every strong tournament of order k and every sole-exit vertex
a -> r, sigma = e_r(T - a); produced by scripts/ready_search.c mode `pairs` from nauty's gentourng) and
data/pairs10_digits.txt (order 10 aggregated as "H: sigma sigma ...").  The inherited C2 law (RESULT.md, JOINT_CLASSIFICATION)
realises 3H + 2 sigma as the Hamilton count of a strong tournament of order k + 2 for every pair.  Reports:
(1) the exact initial segment of odd numbers reached; (2) the odd numbers below 307 not reached and their predecessor
explanation; (3) digit-set sizes; (4) the ready spectra R0, R1, R2 (digits 0, 1, 2) below 3000 and their first gaps."""
import os
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', 'data')
pairs = {}
for k in range(4, 10):
    for l in open(os.path.join(DATA, 'pairs_%d.txt' % k)):
        h, s = map(int, l.split()); pairs.setdefault(h, set()).add(s)
pairs9 = {h: set(v) for h, v in pairs.items()}
for l in open(os.path.join(DATA, 'pairs10_digits.txt')):
    if l.startswith('#'): continue
    h, rest = l.split(':'); pairs.setdefault(int(h), set()).update(int(x) for x in rest.split())
def sweep(P): return set(3 * h + 2 * s for h in P for s in P[h] if s >= 1)
for name, P in (('orders <= 9', pairs9), ('orders <= 10', pairs)):
    sw = sweep(P); top = max(sw)
    miss = [v for v in range(1, top + 1, 2) if v not in sw]
    above = [v for v in miss if v > 305]
    print('(1) %s: %d ready values H (max %d); C2 children 3H+2sigma cover every odd number in [307, %d]; odd numbers below 307 not reached: %s' % (
        name, len(P), max(P), (above[0] - 2) if above else top, [v for v in miss if v < 307]))
print('(2) predecessors of the small misses under 3h, 3h+2, 3h+4 (digit 0, 1, 2): the predecessor h lacks that digit:')
for v in (33, 37, 39, 43, 47, 49, 55, 61, 63, 67, 79, 83, 103, 107, 121, 125, 127, 167, 203, 305):
    if v % 3 == 0: h, d = v // 3, 0
    elif v % 3 == 2: h, d = (v - 2) // 3, 1
    else: h, d = (v - 4) // 3, 2
    print('    %d = 3*%d + %d: digit set of %d is %s' % (v, h, 2 * d, h, sorted(pairs.get(h, set()))[:8]))
import statistics
rows = [(h, len([x for x in pairs[h] if x >= 1]), max(pairs[h])) for h in pairs if max(pairs[h]) >= 1]
print('(3) digit sets (orders <= 10):')
for lo, hi in ((100, 300), (300, 1000), (1000, 2000), (2000, 3000)):
    sub = [r for r in rows if lo <= r[0] < hi]
    print('    H in [%d,%d): %d ready values; mean number of digits %.1f; mean max digit / H %.3f' % (lo, hi, len(sub), statistics.mean(r[1] for r in sub), statistics.mean(r[2] / r[0] for r in sub)))
R = {d: {h for h in pairs if d in pairs[h]} for d in (0, 1, 2)}
print('(4) ready spectra below 3000 (orders <= 10; exact below the order-11 floors 438 / 682 / ? for digits 1 / 2):')
for X in (321, 500, 1000, 2000, 3000):
    odd = range(1, X + 1, 2)
    print('    odd h <= %d: in R1 %d, in R2 %d, in R0 %d, in R1 and R2 %d, of %d' % (X, sum(h in R[1] for h in odd), sum(h in R[2] for h in odd), sum(h in R[0] for h in odd), sum(h in R[1] and h in R[2] for h in odd), len(list(odd))))
print('    first odd h >= 23 not in R1:', [h for h in range(23, 300, 2) if h not in R[1]][:16])
print('    first odd h >= 23 not in R2:', [h for h in range(23, 300, 2) if h not in R[2]][:16])
print('    first odd h >= 23 not in R0:', [h for h in range(23, 600, 2) if h not in R[0]][:16])
