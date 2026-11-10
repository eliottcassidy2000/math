#!/usr/bin/env python3
"""AUDIT 2 (auditor's own code).  The image of the map R -> beta_R on all labelled tournaments of orders 3-6:
how many of the 2^(C(m,2)-1) graphs with an odd number of edges occur, and the fibre sizes.  beta_R from (4.4)
(parity of even-order paths; the function of x7_deletion_rank.py).  Shows that the coincidence of the labelled rank
law with the odd-edge rank law at m = 4, 5 (README section 4.1b) is a coincidence of rank laws only, and that the
affine hull being the whole hyperplane (late item L2) does not make beta_R generic."""
from itertools import combinations
from collections import Counter
import os
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'x7_deletion_rank.py')).read().split("for m in map(int, sys.argv[1:]):")[0]
ns = {}; exec(src, ns)
beta_rows, rank = ns['beta_rows'], ns['rank']
for m in (3, 4, 5, 6):
    prs = list(combinations(range(m), 2)); img = Counter(); byrank = Counter()
    for code in range(1 << len(prs)):
        out = [0] * m
        for k, (i, j) in enumerate(prs):
            if code >> k & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        img[tuple(beta_rows(m, out))] += 1
    for b, c in img.items(): byrank[(rank(list(b)), sum(bin(x).count('1') for x in b) // 2)] += 1
    print("m = %d: %d labelled tournaments -> %d distinct beta_R of the %d graphs with an odd number of edges; fibre sizes (size: number of graphs) %s; distinct beta_R by (rank, edges): %s"
          % (m, 1 << len(prs), len(img), 2 ** (len(prs) - 1), dict(sorted(Counter(img.values()).items())), dict(sorted(byrank.items()))))
