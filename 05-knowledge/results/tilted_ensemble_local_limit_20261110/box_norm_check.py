#!/usr/bin/env python3
"""Box-norm chain for the mod-4 discrepancy (README section 4.1, remark after (4.3)).

D_n = E (-1)^{c_odd(T)} over labelled tournaments on n vertices (= P(H = 1 mod 4) - P(H = 3 mod 4)).
For a tournament on n + k vertices with k disjoint marked pairs let N_k be the number of directed odd cycles that
meet every marked pair in exactly one vertex.  Claim (Cauchy-Schwarz on one vertex star at a time):

      |D_n|^(2^k)  <=  B(n, k) := E (-1)^{N_k}     (uniform tournament on n + k vertices),   1 <= k <= n,
      B(n, k)^2    <=  B(n, k+1),                   B(n, 0) = D_n.

This script computes B(n, k) exactly by brute force (all labelled tournaments, all directed odd cycles) for
n + k <= 6 and checks the chain.  Exact rational arithmetic."""
from fractions import Fraction
from itertools import combinations, permutations

def odd_cycles(N):
    """all directed odd cycles on vertex set range(N) as (vertex mask, list of arcs)"""
    res = []
    for L in range(3, N + 1, 2):
        for S in combinations(range(N), L):
            first = S[0]
            for rest in permutations(S[1:]):
                cyc = (first,) + rest
                res.append((sum(1 << v for v in S), [(cyc[i], cyc[(i + 1) % L]) for i in range(L)]))
    return res

def B(n, k):
    N = n + k
    pairs = [(2 * t, 2 * t + 1) for t in range(k)]          # marked pairs; vertices 2k .. N-1 unpaired
    cycs = [(mask, arcs) for mask, arcs in odd_cycles(N)
            if all(((mask >> a) & 1) + ((mask >> b) & 1) == 1 for a, b in pairs)]
    prs = list(combinations(range(N), 2)); idx = {p: t for t, p in enumerate(prs)}
    # arc (u, w) present iff bit of the pair equals [u < w]
    need = []
    for mask, arcs in cycs:
        on = 0; care = 0
        for u, w in arcs:
            t = idx[(min(u, w), max(u, w))]; care |= 1 << t
            if u < w: on |= 1 << t
        need.append((care, on))
    tot = 0
    for code in range(1 << len(prs)):
        par = 0
        for care, on in need:
            if code & care == on: par ^= 1
        tot += -1 if par else 1
    return Fraction(tot, 1 << len(prs))

ok = True
for n in (3, 4, 5):
    row = []
    for k in range(0, n + 1):
        if n + k > 6: break
        row.append(B(n, k))
    print("n = %d:" % n, "  ".join("B(%d,%d) = %s" % (n, k, v) for k, v in enumerate(row)))
    for k in range(1, len(row)):
        c1 = abs(row[0]) ** (2 ** k) <= row[k]; c2 = row[k - 1] ** 2 <= row[k]
        ok &= c1 and c2
        print("    k = %d: |D_n|^(2^k) = %s <= %s : %s ;  B(n,k-1)^2 = %s <= %s : %s" % (k, abs(row[0]) ** (2 ** k), row[k], c1, row[k - 1] ** 2, row[k], c2))
print("CHAIN HOLDS" if ok else "CHAIN FAILS")
