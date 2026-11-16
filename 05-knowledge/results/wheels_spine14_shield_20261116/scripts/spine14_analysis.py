#!/usr/bin/env python3
"""spine14_analysis.py census.txt: coverage of the odd numbers by the spine-14 two-apex census (second apex beating <= 3 spine
vertices), format 'V h count minpop orient'.  Reports misses below 2^k for k = 16..21, the first misses above 21, the overlap
distribution and orientation sets below 2^18 and 2^19, and the value count."""
import sys, gzip, os
from collections import Counter
p = sys.argv[1]
fh = gzip.open(p, 'rt') if p.endswith('.gz') else open(p)
S = {}; header = ''
for line in fh:
    if line.startswith('#'): header = line.strip()
    elif line.startswith('V '):
        a = line.split(); S[int(a[1])] = (int(a[2]), int(a[3]), int(a[4]))
print(header)
mx = max(S); print(f"distinct values {len(S)}, max {mx}")
for k in range(16, 22):
    B = 1 << k
    if B > mx: break
    miss = [n for n in range(1, B, 2) if n not in S]
    print(f"odd N < 2^{k} = {B}: misses {len(miss)}, first {miss[:8]}")
for k in (18, 19):
    top = 1 << k
    odd = [n for n in range(1, top, 2) if n in S]
    print(f"below 2^{k}: overlap distribution {sorted(Counter(S[n][1] for n in odd).items())}; orientation sets {sorted(Counter(S[n][2] for n in odd).items())}; values attained only with B -> A: {sum(1 for n in odd if S[n][2] == 2)}")
