"""Remark check (THM-4819 (D), "an alternative route that fails"): the m <= 40 left unsettled by the cycle-type route when the
unjustified clause "a full (2m+1)-cycle would make T ... pure blue" is NOT used, i.e. only "an automorphism of even order is
impossible" excludes a gap d = t - t'; a set U of m positions survives iff it is an m-clique of the remaining gaps."""
from math import gcd
def formula(m, d):
    e = gcd(d, m); l = m // e
    return ([1 + 2 * l] + [l] * (2 * e - 2)) if (d // e) % 2 else ([(l + 1) // 2] * 2 + [l] * (2 * e - 1))
def clique(M, D, size):
    S = set(D)
    def ext(cur, cands):
        if len(cur) == size: return cur
        if len(cur) + len(cands) < size: return None
        for i, c in enumerate(cands):
            r = ext(cur + [c], [d for d in cands[i + 1:] if (d - c) % M in S])
            if r: return r
        return None
    return ext([0], [x for x in range(1, M) if x in S])
uns = []
for m in range(3, 41):
    M = 2 * m
    bad = [d for d in range(1, M) if not any(x % 2 == 0 for x in formula(m, d))]
    if clique(M, bad, m): uns.append(m)
print("m <= 40 not settled by 'even order impossible' alone:", uns)
