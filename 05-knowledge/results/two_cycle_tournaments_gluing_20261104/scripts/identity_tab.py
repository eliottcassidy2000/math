#!/usr/bin/env python3
"""identity_tab.py (klein-2026-11-04): #cut arcs + #non-critical vertices - n for every tournament with at most three
Hamiltonian cycles (data/moonK3_{n}.txt from moonstat2 on the K = 3 generator output), by cycle count."""
import sys, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', 'data')
for n in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    c = Counter()
    for l in open(os.path.join(DATA, 'moonK3_%d.txt' % n)):
        p = l.split(); nn, h, od1, id1, ncut, ndel = map(int, p[:6]); c[(h, ncut + ndel - nn)] += 1
    print('order %d: (hc, cut + NC - n): %s' % (n, dict(sorted(c.items()))))
