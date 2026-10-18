"""S40 audit: p = x^2 + 1848 y^2 (p prime, p not dividing 2*3*7*11) iff p splits completely in Q(i, sqrt2, sqrt3, sqrt7, sqrt11),
i.e. (-1/p) = (2/p) = (3/p) = (7/p) = (11/p) = 1. Finite check for p < 10^7; also checks degree 32 = 2 h(-7392)."""
import math
from audit_common import sieve, h_forms
N = 10 ** 7
S = sieve(N)
bad = 0; rep_count = 0; split_count = 0
for p in range(13, N + 1):
    if not S[p]:
        continue
    split = all(pow(a % p, (p - 1) // 2, p) == 1 for a in (-1, 2, 3, 7, 11))
    rep = False
    y = 1
    while 1848 * y * y < p:
        r = p - 1848 * y * y
        x = math.isqrt(r)
        if x * x == r:
            rep = True
            break
        y += 1
    rep_count += rep; split_count += split
    if rep != split:
        bad += 1
print('primes 13 <= p <= %d: represented by x^2+1848y^2: %d; split completely in Q(i,sqrt2,sqrt3,sqrt7,sqrt11): %d; mismatches %d' % (N, rep_count, split_count, bad))
print('h(-7392) =', h_forms(-7392), '; h(-1848) =', h_forms(-1848), '; tower h(-1848*4^k), k=0..5:', [h_forms(-1848 * 4 ** k) for k in range(6)])
