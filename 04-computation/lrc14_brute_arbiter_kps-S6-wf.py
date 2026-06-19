"""
DEFINITIVE arbiter for consec_9 mu_{1/7}.  A says 247/294. Refinement-engine says 4829/5880.
We settle by BRUTE exact computation that is impossible to get wrong:

mu_{1/7}(E) = measure of x in [0,1) with maxgap{frac(e x)} > 1/7.
For E integer with max difference D, the function x->(multiset of frac(e x)) is piecewise
with all breakpoints at rationals of denominator dividing L = lcm{1..D}? No -- breakpoints
where two points coincide are x=m/(e_i-e_j). And maxgap=theta crossing at x=(k+theta)/(e_i-e_j).
theta=1/7. So ALL breakpoints have denominator dividing 7*L where L=lcm of all |e_i-e_j|.
=> mu is a rational with denominator dividing 7*L, and is exactly computed by midpoint
evaluation over the partition by ALL rationals with denominator dividing 7*L in [0,1].

We do EXACTLY that: enumerate x = j/(7L), j=0..7L-1, and for each CELL (j/(7L),(j+1)/(7L))
take the midpoint (2j+1)/(14L) and test maxgap>theta. Sum cell widths = (#hits)/(7L).
This partition is a refinement of every possible breakpoint set, so it is GROUND TRUTH,
modulo: we must ensure no breakpoint lies strictly inside a cell. Breakpoints have denom | 7L,
so they are exactly the grid points j/(7L); none lies strictly inside (j/(7L),(j+1)/(7L)).
QED this is exact. Only cost: 7L cells.
"""
from fractions import Fraction as F
from math import gcd

def lcm(a,b): return a*b//gcd(a,b)

def maxgap_at(E, x):
    pts = sorted(set((F(e)*x) % 1 for e in E))
    if len(pts)==1: return F(1)
    g=F(0)
    for i in range(len(pts)):
        if i==len(pts)-1:
            d=(pts[0]+1)-pts[i]
        else:
            d=pts[i+1]-pts[i]
        if d>g: g=d
    return g

def mu_brute(E,theta):
    E=sorted(set(E))
    diffs=set()
    for i in range(len(E)):
        for j in range(len(E)):
            if i!=j and E[i]!=E[j]:
                diffs.add(abs(E[i]-E[j]))
    L=1
    for d in diffs: L=lcm(L,d)
    # theta=1/7 -> multiply by 7
    num,den = theta.numerator, theta.denominator
    M = L*den   # grid denominator; breakpoints denom | M
    hits=0
    # iterate cells (j/M,(j+1)/M), midpoint (2j+1)/(2M)
    for j in range(M):
        mid = F(2*j+1, 2*M)
        if maxgap_at(E,mid)>theta:
            hits+=1
    return F(hits, M), M

theta=F(1,7)
for k in [8,9]:
    E=list(range(k))
    val,M = mu_brute(E,theta)
    print(f"BRUTE consec_{k}: mu={val}={float(val):.6f}  (grid M={M})")
    print(f"   engine-A anchor: ", {8:F(691,735),9:F(247,294)}[k], "=", float({8:F(691,735),9:F(247,294)}[k]))
    print(f"   refine-engine  : ", {8:F(691,735),9:F(4829,5880)}[k], "=", float({8:F(691,735),9:F(4829,5880)}[k]))
