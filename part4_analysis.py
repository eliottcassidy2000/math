from fractions import Fraction as F
from math import exp
from itertools import combinations
from lrc import measure_B, measure_intersection, lonely_measure

def baseline(m,n): return (F(1)-F(2,n))**m
def moments(speeds,n):
    m=len(speeds)
    EN=sum(measure_B(v,n) for v in speeds)
    S2=sum(measure_intersection([speeds[i],speeds[j]],n) for i,j in combinations(range(m),2))
    ENN1=2*S2; EN2=ENN1+EN; Var=EN2-EN*EN
    return EN,ENN1,EN2,Var

# WHY the second moment method fails. Paley-Zygmund gives P(N>0) >= E[N]^2/E[N^2].
# That LOWER bounds the BUSY measure, not the lonely one. Useless for LRC.
# To show P(N=0)>0 we need an UPPER tail / anti-concentration, not 2nd moment.
# Demonstrate: compute the Paley-Zygmund lower bound on measure{N>=1} and note
# it never forces measure{N=0}>0.

print("Tight (mu=0) sets vs their full moment profile (n=5 and n=6):")
tights=[([1,2,3,4],5),([1,3,4,7],5),([1,2,3,4,5],6),([1,3,4,5,9],6)]
for s,n in tights:
    EN,ENN1,EN2,Var=moments(s,n)
    mu=lonely_measure(s,n)
    print(f"  {s} n={n}: E[N]={float(EN):.4f} E[N^2]={float(EN2):.4f} Var={float(Var):.4f} mu={mu}")
    # max N actually attained?
print()

# For these tight sets N>=1 everywhere (mu=0). Check distribution of N: what values does N take?
def N_distribution(speeds,n):
    arcs=[]
    from lrc import arc_endpoints
    all_arcs=[arc_endpoints(v,n) for v in speeds]
    pts=set([F(0),F(1)])
    for a in all_arcs:
        for lo,hi in a: pts.add(lo); pts.add(hi)
    pts=sorted(pts)
    from collections import defaultdict
    dist=defaultdict(F)
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2
        cnt=0
        for arcs in all_arcs:
            for lo,hi in arcs:
                if lo<mid<hi: cnt+=1; break
        dist[cnt]+=b-a
    return dict(dist)

print("Distribution of N(t) (measure of {N=k}) for tight sets:")
for s,n in tights:
    d=N_distribution(s,n)
    print(f"  {s} n={n}: " + "  ".join(f"N={k}:{float(v):.3f}" for k,v in sorted(d.items())))
print()

# CRUX: anti-concentration. For LRC we need min_t N(t)=0. The moment method bounds
# AVERAGE behavior. A clean sufficient condition for mu>0: if SUM of measures of
# bad arcs over a fundamental cell < cell, but that's just E[N]<1 fails (E[N]~2).
# Show: when does a 2nd-moment style argument give mu>0?
# Cantelli/Chebyshev one-sided: measure{N=0} = measure{N<=0} ; but N>=0 integer,
# {N=0} = {N < 1}. One-sided Chebyshev (Cantelli): for X with mean mu_N var s2,
# P(X - mu_N <= -lambda) <= s2/(s2+lambda^2). With lambda=mu_N (so X<=0):
# P(N=0) <= Var/(Var+E[N]^2). That's an UPPER bound on mu -> cannot prove mu>0.
print("Cantelli UPPER bound on mu: Var/(Var+E[N]^2). (only upper => can't prove LRC)")
for s,n in [([1,3,5,7],5),([1,2,5,7],5),([1,2,3,4],5)]:
    EN,ENN1,EN2,Var=moments(s,n)
    ub=float(Var/(Var+EN*EN))
    print(f"  {s} n={n}: Cantelli upper bound mu<={ub:.4f}, actual mu={float(lonely_measure(s,n)):.4f}")
