#!/usr/bin/env python3
"""HYP-9657 (exploratory, in-sample): mass profile of the three-letter rooted system S, B, F (the O/Z/J grammar).

States (h, R, e, u) from (1, 1, 1, 1) under  S: (h, 2R, 0, u),  B: (h+R-e, 2R-u, u, u),  F: (h+R-e, h+R, h-e, h)
(rules checked against brute force in psi_sbf.py).  m(N) = number of distinct states reached right after a B or F
letter with h = N; such states have R <= 2N, so the cap R <= 4X loses none with h <= X.
Output: the odd N <= X with m(N) = 0, class means of m modulo small primes, m on primes and on Rabinowitsch-rough N.
Run:  python3 ozj_mass.py > ozj_mass.out
"""
from collections import Counter
import numpy as np
def step(st, ch):
    h, R, e, u = st
    if ch == "S": return (h, 2*R, 0, u)
    if ch == "B": return (h + R - e, 2*R - u, u, u)
    return (h + R - e, h + R, h - e, h)
X = 12000
seen = {(1, 1, 1, 1)}
frontier = [(1, 1, 1, 1)]
bf = set()
while frontier:
    nxt = []
    for st in frontier:
        for ch in "SBF":
            t = step(st, ch)
            if t[0] <= X and t[1] <= 4 * X:
                if ch != "S": bf.add(t)
                if t not in seen:
                    seen.add(t); nxt.append(t)
    frontier = nxt
m = Counter(s[0] for s in bf)
odd = np.arange(1, X + 1, 2)
mv = np.array([m.get(int(N), 0) for N in odd])
print("states:", len(seen), " B/F-sheet states:", len(bf), " misses:", odd[mv == 0].tolist())
lo = X // 4
sel = odd >= lo
print("window [%d, %d]: mean m = %.2f" % (lo, X, mv[sel].mean()))
for p in (3, 5, 7, 11, 13):
    res = odd[sel] % p
    means = [mv[sel][res == c].mean() for c in range(p)]
    print("p=%2d class means:" % p, " ".join("%.1f" % x for x in means))
def is_prime(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True
rough = np.array([is_prime(int(N)) and N % 4 == 3 and is_prime(int((N + 1) // 4)) for N in odd])
print("Rabinowitsch-rough in window: n = %d, mean m = %.2f (ratio %.3f)" % (int((rough & sel).sum()), mv[rough & sel].mean(), mv[rough & sel].mean() / mv[sel].mean()))
pr = np.array([is_prime(int(N)) for N in odd])
print("primes in window: mean m = %.2f (ratio %.3f)" % (mv[pr & sel].mean(), mv[pr & sel].mean() / mv[sel].mean()))
print("m at Heegner numbers:", {d: m.get(d, 0) for d in (3, 7, 11, 19, 43, 67, 163)})
small = [(int(N), int(v)) for N, v in zip(odd, mv) if N > 2647 and v <= 3][:30]
print("values above 2647 with m <= 3:", small)
