"""opus-2026-07-19-S402: positive controls for the G-K pinning session.
(1) Fan-Sun [16 in G-K v4] family ML(8,4r+3,4r+11,4r+19) = (2r+7)/(8r+30)
    -- the literature refutation of G-K Conjecture 1.1 (k=2 rungs at n=4);
    verifying it validates both the citation and this repo's scanner.
(2) gridmax((1,...,11); 14) = 1/14 -- the window-boundary value is realized
    in S*_0-land (the terminal set of G-K Thm 1.4's proven chain), showing
    the proven chain alone cannot floor accumulation targets above 1/14."""
from fractions import Fraction
from math import gcd
def scan_max(V,qmax):
    best=(Fraction(0),Fraction(0))
    for q in range(2,qmax+1):
        for a in range(1,q//2+1):
            if gcd(a,q)!=1: continue
            m=min(Fraction(min((v*a)%q,q-(v*a)%q),q) for v in V)
            if m>best[0]: best=(m,Fraction(a,q))
    return best
print("Fan-Sun ML(8,4r+3,4r+11,4r+19) vs (2r+7)/(8r+30):")
for r in range(4):
    V=[8,4*r+3,4*r+11,4*r+19]
    M,t=scan_max(V,120); pred=Fraction(2*r+7,8*r+30)
    print(f"  r={r} V={V}: M={M} at t={t}; predicted {pred}; match={M==pred}")
g=max(min(min((i*j)%14,14-(i*j)%14) for i in range(1,12)) for j in range(14))
print("gridmax((1,...,11); 14) =", Fraction(g,14), "== 1/14:", Fraction(g,14)==Fraction(1,14))
