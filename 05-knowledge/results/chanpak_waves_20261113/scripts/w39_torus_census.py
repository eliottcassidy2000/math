# Wave 39 (mac-mini-2026-11-15-frame). Run from a working directory; w39_lit_setup.py writes lit_base.pkl used by the others.
import os
# Wave 39: the "torus" of the span-3 ribbon semigroup.  For every word W of length <= LMAX over the 8 literal letters,
# test whether M_W has rank <= 1 in the 62-dim integer representation (span3_intrep.txt); for minimal rank-one words
# record the nonzero eigenvalue lambda = trace(M_W) (H(u W^k v) = c lambda^k).  Report the primes dividing the lambdas.
import sys, numpy as np
from collections import Counter
f = open(os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..', '04-computation', 'experiments', 'ribbon_spectrum_20261104'), 'span3_intrep.txt'))
r, L = map(int, f.readline().split())
y0 = np.array(list(map(int, f.readline().split())), dtype=np.float64)
hc = np.array(list(map(int, f.readline().split())), dtype=np.float64)
M = {}
for _ in range(L):
    mask = int(f.readline())
    M[mask] = np.array([list(map(int, f.readline().split())) for _ in range(r)], dtype=np.float64)
LMAX = int(sys.argv[1])


def factor(n):
    n = abs(int(round(n)))
    fs = Counter()
    p = 2
    while p * p <= n:
        while n % p == 0:
            fs[p] += 1
            n //= p
        p += 1
    if n > 1:
        fs[n] += 1
    return fs


lam = Counter()
primes = Counter()
stats = Counter()
examples = {}


def rec(word, P):
    if word:
        stats['words'] += 1
        s = np.linalg.svd(P, compute_uv=False)
        if s[0] == 0 or s[1] <= 1e-9 * s[0]:
            stats['minimal_rank_le1'] += 1
            t = np.trace(P)
            if abs(t) < 0.5:
                stats['nilpotent'] += 1
            else:
                ti = int(round(t))
                lam[ti] += 1
                for p in factor(ti):
                    primes[p] += 1
                examples.setdefault(ti, word)
            return
    if len(word) == LMAX:
        return
    for b in range(8):
        rec(word + (b,), M[b] @ P)


rec((), np.identity(r))
print(dict(stats))
print('distinct nonzero eigenvalues', len(lam))
print('primes dividing eigenvalues (number of minimal rank-one words):', sorted(primes.items()))
bad = sorted(k for k in lam if set(factor(k)) - {3, 5})
print('eigenvalues not of the form +-3^a 5^b:', bad[:60], '... total', len(bad))
print('eigenvalue multiset (value: count):', sorted(lam.items())[:60])
print('examples:', {k: examples[k] for k in sorted(lam)[:30]})
