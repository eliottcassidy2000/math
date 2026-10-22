"""Audit S50 - a8: is primality the reason the repunits 31, 127 are rare at n = 8, 9 (note section 3)?
Class multiplicities N(n,h) from the mac-mini gentourng census (data).  For each window +-10% around 31 and 127,
list N(n,h) for every odd h, marking primes."""
import gzip, os
from collections import Counter
from sympy import isprime
MM = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'tournament_pseudoprimes_20261021')
for n in (8, 9):
    path = os.path.join(MM, f'census_n{n}.txt')
    fh = open(path) if os.path.exists(path) else gzip.open(path + '.gz', 'rt')
    N = Counter(int(ln.split()[0]) for ln in fh if ln.strip())
    for v in (31, 127):
        lo, hi = int(v * 0.9) | 1, int(v * 1.1)
        row = [f'{h}{"p" if isprime(h) else ""}:{N[h]}' for h in range(lo, hi + 1, 2)]
        print(f'n = {n}, window of {v}: ' + ' '.join(row))

# class multiplicities by residue mod 16, each N(n,h) divided by the median of its +-10% window (local baseline)
import statistics
for n in (8, 9):
    path = os.path.join(MM, f'census_n{n}.txt')
    fh = open(path) if os.path.exists(path) else gzip.open(path + '.gz', 'rt')
    N = Counter(int(ln.split()[0]) for ln in fh if ln.strip())
    top = max(N)
    rel = {r: [] for r in range(1, 16, 2)}
    for h in range(25, int(0.6 * top), 2):
        win = [N[x] for x in range(int(h * 0.9) | 1, int(h * 1.1) + 1, 2) if x != h]
        med = statistics.median(win)
        if med > 0:
            rel[h % 16].append(N[h] / med)
    print(f'n = {n}: mean of N(n,h)/local median over odd h in [25, {int(0.6 * top)}], by h mod 16: ' +
          ', '.join(f'{r}: {statistics.mean(v):.3f}' for r, v in rel.items()))
print('2^k - 1 = 15 (mod 16) and 2^k + 1 = 1 (mod 16) for k >= 4; (4^k - 1)/3 = 5 (mod 16) for k >= 2 (85 = 5, 341 = 5, 1365 = 5 mod 16)')
