#!/usr/bin/env python3
"""holes_grammar.py (R-2): the near-join sub-grammar excludes 7 and 21 by the B-sets of the three smallest counts.
B-sets of the counts 1, 3, 5 over all rooted tournaments of orders <= 8 (every class, every root; with the end-count
e_t); the generators of each B-set as odd parts; and the exhaustive check that 7 and 21 are not h + 2^k B (k >= 1) nor
h + B with e = 0 for any rooted tournament of order <= 8 with h < 21."""
exec(open('nearjoin_check.py').read().split('# (1)')[0])
from collections import defaultdict
def rev(T): n = len(T); return [[T[j][i] for j in range(n)] for i in range(n)]
def classes1(n): return [[[0]]] if n == 1 else ([[[0, 1], [0, 0]]] if n == 2 else classes(n))
Bset = defaultdict(set); vals = set(); vals0 = set()
for n in range(1, 9):
    for T in classes1(n):
        N = 1 << n; dp, H = ham_all(T); dpr, _ = ham_all(rev(T)); h = H[N-1]
        if h >= 21: continue
        for t in range(n):
            B = sum(dpr[m][t] * H[(N-1) & ~m] for m in range(1, N) if (m >> t) & 1); e = dp[N-1][t]
            Bset[h].add(B)
            for k in range(1, 6): vals.add(h + (B << k))
            if e == 0: vals0.add(h + B)
def oddpart(x):
    while x % 2 == 0: x //= 2
    return x
for h in (1, 3, 5, 9, 11, 13, 15, 17, 19):
    bs = sorted(Bset[h]); print(f"B-set of h = {h} (orders <= 8): {bs[:14]}{' ...' if len(bs) > 14 else ''}; odd parts: {sorted(set(oddpart(b) for b in bs))}")
print("7 in {h + 2^k B}:", 7 in vals, "; 21 in {h + 2^k B}:", 21 in vals, "; 7 in dead-end appends:", 7 in vals0, "; 21 in dead-end appends:", 21 in vals0)
print("values h + 2^k B and dead-end appends below 30 from bases with h < 21:", sorted(v for v in vals | vals0 if v < 30))
