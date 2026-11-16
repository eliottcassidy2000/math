#!/usr/bin/env python3
"""two_apex_analysis.py: coverage of the odd numbers by two-apex transitive tapes, from the census files
tapex_<m>.txt (all tapes, "V h count minpop orient") and tapexs_<m>.txt (second apex beating <= 3 spine vertices),
possibly gzip-compressed.  usage: two_apex_analysis.py data_dir"""
import sys, os, gzip
from collections import Counter
D = sys.argv[1]
def load(name):
    p = os.path.join(D, name); d = {}
    fh = gzip.open(p + '.gz', 'rt') if not os.path.exists(p) else open(p)
    for line in fh:
        if line.startswith('V '):
            a = line.split(); d[int(a[1])] = (int(a[2]), int(a[3]), int(a[4]))
    return d
# one-apex family: 1 + p(u), p the inversion statistic (apex packet)
M = 22; p = [0] * (1 << M)
for u in range(1, 1 << M):
    w = u >> 1; p[u] = p[w] if u & 1 else p[w] + w + (w & 1)
one = set(1 + x for x in p)
full = {m: load(f'tapex_{m}.txt') for m in range(4, 11)}
small = {m: load(f'tapexs_{m}.txt') for m in (11, 12, 13) if os.path.exists(os.path.join(D, f'tapexs_{m}.txt')) or os.path.exists(os.path.join(D, f'tapexs_{m}.txt.gz'))}
print("all two-apex tapes:")
for m in sorted(full):
    S = full[m]; mx = max(S)
    un = [n for n in range(1, mx, 2) if n not in S and n not in (7, 21)]
    print(f"  spine {m}: {len(S)} values, max {mx}, odd values below the max not attained (besides 7, 21): {len(un)}, first {un[:5]}")
S10 = full[10]
for B in (2**10, 2**12, 2**14):
    odd = range(1, B, 2); m1 = sum(1 for n in odd if n not in one); m2 = [n for n in odd if n not in S10]
    print(f"  odd N < {B}: one-apex misses {m1} ({100*m1/(B//2):.1f} %); two-apex spine <= 10 misses {m2}")
print("  overlap needed below 2^14 (min over realisations of min(popcount a, popcount b)):", sorted(Counter(S10[n][1] for n in range(1, 2**14, 2) if n in S10).items()))
print("  least spine attaining each odd N < 130:", [(n, next(m for m in range(4, 11) if n in full[m])) for n in range(1, 130, 2) if n in S10])
print("second apex beating at most three spine vertices:")
for m in sorted(small):
    S = small[m]; mx = max(S)
    for B in (2**14, 2**15, 2**16, 2**17, 2**18):
        if B > mx: break
        miss = [n for n in range(1, B, 2) if n not in S]
        print(f"  spine {m}: {len(S)} values, max {mx}; odd N < {B}: misses {len(miss)}, first {miss[:6]}")
    top = 2 ** (len(bin(mx)) - 3)
    odd = [n for n in range(1, top, 2) if n in S]
    print(f"    overlap distribution below {top}: {sorted(Counter(S[n][1] for n in odd).items())}; orientation sets: {sorted(Counter(S[n][2] for n in odd).items())}")
    only_rev = [n for n in range(1, top, 2) if n in S and S[n][2] == 2]
    print(f"    values below {top} attained only with B -> A: {len(only_rev)} {only_rev[:10]}")
