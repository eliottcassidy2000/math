"""AUDIT D1: summarise aw_verify output on aw_prop33.py's graphs (realisers of X(n) and their parity).
usage: python3 aw_prop33_summary.py aw_prop33_verify.out"""
import sys
from collections import defaultdict
hits = defaultdict(list); tot = defaultdict(int)
for line in open(sys.argv[1]):
    p = line.split()
    if p[0] == 'SUMMARY': continue
    n, X, d, par = int(p[2]), int(p[3]), int(p[4][2:]), p[5]; tot[n] += 1
    if d == X: hits[n].append((p[0], par))
allowed = lambda sh: sorted(c for c in sh if c != 'P2')
for n in sorted(tot):
    ev = [h for h in hits[n] if h[1] == 'EVEN']
    le1 = [h[0] for h in hits[n] if h[0].count('P1,') + h[0].endswith('P1]') <= 1]
    ok_shapes = all(allowed(s[3:-1].split(',')) in (['P3'], ['P1', 'P4'], ['C4', 'P1'], ['P1', 'P3', 'P3']) for s in le1)
    print(f"n={n}: {tot[n]} graphs with complement of max degree <= 2; realisers of X(n) = {2**n - n - 2}: "
          f"{len(hits[n])}, even: {len(ev)}; those with <= 1 universal vertex: {len(le1)}, all of the listed shapes: {ok_shapes}")
