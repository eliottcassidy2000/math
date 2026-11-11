#!/usr/bin/env python3
"""Sum the part files of cycrace.c (lines: K signed weight; the weight of a part is its total labelled weight).
Usage: merge_cyc.py out.txt part ..."""
import sys
from collections import defaultdict
sg = defaultdict(int); wt = 0
for fn in sys.argv[2:]:
    w0 = None
    for line in open(fn):
        K, s, w = map(int, line.split()); sg[K] += s; w0 = w
    wt += w0
with open(sys.argv[1], 'w') as f:
    for K in sorted(sg): f.write("%d %d %d\n" % (K, sg[K], wt))
