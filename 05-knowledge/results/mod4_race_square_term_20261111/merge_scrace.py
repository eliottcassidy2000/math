#!/usr/bin/env python3
"""Sum the part files of scrace2.c (lines: sc strong sym r classes labelled).  Usage: merge_scrace.py out.txt part ..."""
import sys
from collections import defaultdict
d = defaultdict(lambda: [0, 0])
for fn in sys.argv[2:]:
    for line in open(fn):
        a, b, c, r, cl, l = map(int, line.split()); d[(a, b, c, r)][0] += cl; d[(a, b, c, r)][1] += l
with open(sys.argv[1], 'w') as f:
    for k in sorted(d): f.write("%d %d %d %d %d %d\n" % (k + tuple(d[k])))
