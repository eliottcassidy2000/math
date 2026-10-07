"""Heuristic expected number of integer cycles of shape (L,k):
   E(L,k) = Lyn(L,k) / |2^L - 3^k|, Lyn = # binary Lyndon words of length L with k ones,
   restricted to k >= L/2 (product identity; L >= 2).  Rank shapes, mark Ellison pairs,
   and compare the ranking with Ellison's criterion |2^x - 3^y| < 2^x e^{-x/10}."""
from math import comb, exp, log
from sympy import mobius, divisors, factorint
def lyn(L, k):
    from math import gcd
    g = gcd(L, k) if k else L
    return sum(mobius(e) * comb(L // e, k // e) for e in divisors(g)) // L
ell = {(13,8),(14,9),(16,10),(19,12),(27,17)}
rows = []
for L in range(2, 121):
    for k in range((L + 1) // 2, L + 1):
        m = abs(2**L - 3**k)
        E = lyn(L, k) / m
        rows.append((E, L, k, m, lyn(L, k)))
rows.sort(reverse=True)
print("top 30 shapes by E(L,k) = Lyn/|2^L-3^k|, L<=120:")
for E, L, k, m, n in rows[:30]:
    tag = []
    if (L, k) in ell: tag.append("ELLISON")
    if abs(2**L - 3**k) < 2**L * exp(-L/10): tag.append("ell-crit")
    if (L,k) in {(1,1),(2,1),(3,2),(11,7)}: tag.append("CYCLE")
    print("  (%3d,%3d) m=%-22d Lyn=%-12d E=%.4f %s %s" % (L, k, m, n, E, " ".join(tag), dict(factorint(m)) if m < 10**12 else ""))
tot = sum(r[0] for r in rows)
print("sum of E over all shapes L<=120: %.4f" % tot)
print("sum of E over L<=120 excluding |m|=1 shapes: %.4f" % sum(r[0] for r in rows if r[3] > 1))
# Ellison criterion set for x <= 120
crit = [(L, k) for L in range(1, 121) for k in range(0, L + 1) if abs(2**L - 3**k) < 2**L * exp(-L/10)]
print("Ellison criterion |2^x-3^y| < 2^x e^{-x/10}, 1<=x<=120:", crit)
