"""Reflection products on a Hamilton cycle plus one vertex (prescribed black endpoint lemma, THM-4819).
V = {v} u Z_{2m}.  sigma_t: v <-> t, x <-> 2t - x (x != t), the positional reversal of the path
(c_t, c_{t+1}, ..., c_{t-1}, v) read along a Hamilton cycle c_0 c_1 ... c_{2m-1} of T - v.  g_d = sigma_d sigma_0.
Prints, for each m, the gaps d whose product has ODD order and whether some m-subset of Z_{2m} has all its
pairwise gaps odd-order (a clique of size m in Cay(Z_{2m}, D_odd)); if none, every regular T on 2m+1 vertices has a
black Hamilton path ending at any prescribed vertex."""
import sys
from math import gcd
from itertools import combinations
def sigma(t, M):
    p = {'v': t, t: 'v'}
    for x in range(M):
        if x != t: p[x] = (2 * t - x) % M
    return p
def order_parity(perm):
    seen = set(); lens = []
    for s in perm:
        if s in seen: continue
        L = 0; x = s
        while x not in seen: seen.add(x); x = perm[x]; L += 1
        lens.append(L)
    return lens
def clique_free(M, Dodd, size):
    # is there a size-subset of Z_M containing 0 with all pairwise differences in Dodd?  (backtracking)
    S = set(Dodd)
    def ext(cur, cands):
        if len(cur) == size: return cur
        for i, c in enumerate(cands):
            r = ext(cur + [c], [d for d in cands[i+1:] if (d - c) % M in S])
            if r: return r
        return None
    return ext([0], [x for x in range(1, M) if x in S])
for m in range(3, int(sys.argv[1]) + 1):
    M = 2 * m; s0 = sigma(0, M); Dodd = []; types = {}
    for d in range(1, M):
        sd = sigma(d, M); g = {x: sd[s0[x]] for x in s0}
        lens = order_parity(g)
        if all(L % 2 == 1 for L in lens): Dodd.append(d)
        types[d] = sorted(lens, reverse=True)
    cl = clique_free(M, Dodd, m)
    print(f"m={m} n={2*m+1}: odd-order gaps {Dodd}; m-clique of odd gaps: {cl}")
    if m <= 6: print("   cycle types:", {d: types[d] for d in range(1, M)})

def closed_form(m, d):
    """cycle type of g_d from the formula (e = gcd(d, m), l = m/e)."""
    from math import gcd
    e = gcd(d, m); l = m // e
    if (d // e) % 2 == 1: return sorted([1 + 2 * l] + [l] * (2 * e - 2), reverse=True)
    return sorted([(l + 1) // 2] * 2 + [l] * (2 * e - 1), reverse=True)
def good(m, d):
    t = closed_form(m, d)
    return any(L % 2 == 0 for L in t) or t == [2 * m + 1]   # even order, or an n-cycle (vertex-transitive)
if __name__ == '__main__' and len(sys.argv) > 2:
    for m in range(3, int(sys.argv[2]) + 1):
        M = 2 * m; s0 = sigma(0, M)
        for d in range(1, M):   # the closed form agrees with the direct product
            sd = sigma(d, M); g = {x: sd[s0[x]] for x in s0}
            assert sorted(order_parity(g), reverse=True) == closed_form(m, d), (m, d)
        bad = [d for d in range(1, M) if not good(m, d)]
        cl = clique_free(M, bad, m)
        print(f"m={m} n={2*m+1}: bad gaps {bad}; m-clique of bad gaps: {cl}")
