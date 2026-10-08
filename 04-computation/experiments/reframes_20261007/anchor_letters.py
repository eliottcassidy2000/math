#!/usr/bin/env python3
"""First letter of each legal anchor (cycle point = 1 mod 4) versus barrier status in anchor_census.out.
First letter a = v_2(3c + 1) >= 2 for c = 1 mod 4; a = 2 iff c = 1 mod 8."""
import sys, re
sys.path.insert(0, '../runcompress_20261007')
from transitions import cycle_reps
from fractions import Fraction as Fr
def v2(q):
    n, d = q.numerator, q.denominator; v = 0
    while n % 2 == 0: n //= 2; v += 1
    return v
anchors = []
for w, c, cy in cycle_reps(9, 4):
    for p in cy:
        if p.numerator % 2 == 0: continue
        if (p.numerator * pow(p.denominator, -1, 4)) % 4 != 1: continue
        anchors.append(p)
bar = set()
for line in open('anchor_census.out'):
    m = re.match(r'\s+anchor\s+(\S+)\s+word .* odd density', line)
    if m: bar.add(Fr(m.group(1)))
from collections import Counter
tab = Counter((v2(3 * p + 1) == 2, p in bar) for p in anchors)
print(f"{len(anchors)} legal anchors; {len(bar)} barriers parsed")
print("(first letter = 2, barrier) counts:", dict(tab))
