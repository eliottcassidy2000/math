#!/usr/bin/env python3
"""moon_tab.py (klein-2026-11-04): tabulate the deletion statistics of two-cycle tournaments from moonstat2 output
(data/moon2_{n}.txt, one line per class: n hc od1 id1 ncut ndel : per non-critical vertex 'hc(Y-v),|O|,#intervals,k1,head1,tail1,se(U),ie(U),v').
Reports per order: #cut arcs + #non-critical vertices - n; (sole exits, sole entries) joint counts; the number of vertices
whose deletion leaves a uniquely Hamiltonian tournament; the range of the non-critical count; and the law check
D_n = a_{n-1} + (n-4) b_{n-1} with a_m = F(2m-8), b_m = F(2m-7) (THM-4868 (F): one / two sole exits)."""
import sys, os
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, '..', 'data')
def fib(k):
    a, b = 0, 1
    for _ in range(k): a, b = b, a + b
    return a
for n in range(int(sys.argv[1]), int(sys.argv[2]) + 1):
    c = Counter(); se = Counter(); ud = Counter(); nd = Counter(); tot = 0
    for l in open(os.path.join(DATA, 'moon2_%d.txt' % n)):
        p = l.split(); nn, h, od1, id1, ncut, ndel = map(int, p[:6])
        if h != 2: continue
        tot += 1; c[ncut + ndel - nn] += 1; se[(od1, id1)] += 1
        tuples = [tuple(map(int, t.split(','))) for t in p[7:]]
        ud[sum(1 for t in tuples if t[0] == 1)] += 1; nd[ndel] += 1
    a, b = fib(2 * (n - 1) - 8), fib(2 * (n - 1) - 7)
    print('order %d: %d two-cycle classes = a_{n-1} + (n-4) b_{n-1} = %d + %d*%d = %d; cut+NC-n: %s; (sole exits, sole entries): %s; '
          'unique-cycle deletions per class: %s; non-critical count range %d..%d'
          % (n, tot, a, n - 4, b, a + (n - 4) * b, dict(sorted(c.items())), dict(sorted(se.items())), dict(sorted(ud.items())), min(nd), max(nd)))
