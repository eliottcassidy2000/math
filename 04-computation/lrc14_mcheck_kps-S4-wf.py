from fractions import Fraction as F
from math import gcd
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
# canonical counterexample
S=[1,6,8,10,15,16,18,22,24,26,28,378,581]
# M is a max over a finite set of "event" tau where some ||v tau|| is locally extremal.
# The standard candidate set = all fractions j/v and (j/v + j'/v')/... ; the cand() in
# the toolkit uses midpoints 1/(2v) and 1/d for pairwise sums/diffs d. The TRUE M is
# attained at a vertex of the arrangement: tau where two runners tie at the min, i.e.
# ||a tau|| = ||b tau||. Those are tau = j/(a+b) or j/(a-b) and j/(2a). cand() includes
# exactly these. To be safe, ALSO scan a fine rational grid as a sanity lower bound.
best=F(0); bestt=None
# toolkit cand
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
for t in cand(S):
    m=min(nrm(x*t) for x in S)
    if m>best: best=m; bestt=t
print("toolkit M =",best,"~",float(best),"at tau=",bestt)
# fine grid sanity (lower bound): denominators up to 3000
gb=F(0); gt=None
D=4000
for num in range(1,D//2+1):
    t=F(num,D)
    m=min(nrm(x*t) for x in S)
    if m>gb: gb=m; gt=t
print("grid (den 4000) best lower bound =",gb,"~",float(gb),"at",gt)
print(">=1/14?", best>=F(1,14), " (1/14=",float(F(1,14)),")")
