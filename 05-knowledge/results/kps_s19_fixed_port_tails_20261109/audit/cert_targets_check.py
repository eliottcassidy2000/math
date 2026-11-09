"""S19 audit: the LP certificate targets are admissible (odd / = 1 mod 3), in range and distinct.
Usage: python cert_targets_check.py CERT X LO"""
import sys
from fractions import Fraction
path, x, lo = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
ts, tot = [], Fraction(0)
for l in open(path):
    if l.startswith('#'): continue
    t, y = l.split(); ts.append(int(t)); tot += Fraction(y)
ok = all(lo <= t <= 200000 and (t % 2 == 1 if x == 2 else t % 3 == 1) for t in ts)
print(path, 'targets', len(ts), 'all admissible and in range:', ok, 'distinct:', len(set(ts)) == len(ts), 'sum y =', float(tot))
