"""
BACK-TRANSFER (metagraph lesson -> LRC): do the LRC danger events have MIXED-SIGN covariance
(like HP-indicators, FKG fails) or all >=0 (FKG, as HYP-2823 assumed)?
Cov(danger_i,danger_j) = meas(||s_i t||<1/14 & ||s_j t||<1/14) - (1/7)^2 over t.
Also: idea 6 -- is the metagraph variance carried by the HIGH-H end (compatible pairs cluster there)?
"""
from fractions import Fraction as F
def nrm(x): f=x-(x.numerator//x.denominator); return min(f,1-f)
def covs(S,n=14,Q=2520):
    m=len(S); E=[F(0)]*m; J={}
    dang=[]
    for a in range(Q):
        t=F(a,Q); d=[nrm(s*t)<F(1,n) for s in S]; dang.append(d)
    for i in range(m): E[i]=F(sum(1 for d in dang if d[i]),Q)
    res=[]
    for i in range(m):
        for j in range(i+1,m):
            both=F(sum(1 for d in dang if d[i] and d[j]),Q)
            cov=both-E[i]*E[j]
            res.append((S[i],S[j],float(cov)))
    return res
ap=list(range(1,14))
r=covs(ap)
pos=[x for x in r if x[2]>1e-9]; neg=[x for x in r if x[2]<-1e-9]; zer=[x for x in r if abs(x[2])<=1e-9]
print(f"AP {{1..13}}: {len(r)} danger-event pairs:  positive Cov={len(pos)}  negative Cov={len(neg)}  ~zero={len(zer)}")
print(f"  => LRC danger events are MIXED-SIGN: {len(neg)>0} (FKG also fails for LRC, matching the metagraph)")
print(f"  most negative pairs: {sorted(neg,key=lambda x:x[2])[:5]}")
print(f"  most positive pairs: {sorted(pos,key=lambda x:-x[2])[:5]}")
# the negative pairs: which speed pairs? (resonance structure)
print(f"  sample positive pair speeds (resonant, gcd): {[(a,b) for a,b,c in sorted(pos,key=lambda x:-x[2])[:6]]}")
print()
# IDEA 6 check: metagraph variance carried by high-H end? compatible HP-pairs cluster at high H.
from itertools import permutations
def Hadj(n,adj):
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(1<<n):
        for last in range(n):
            c=dp[mask][last]
            if not c or not (mask>>last)&1: continue
            av=adj[last]&~mask
            while av:
                nb=av&-av; nx=nb.bit_length()-1; dp[mask|nb][nx]+=c; av^=nb
    return sum(dp[(1<<n)-1])
# contribution to E[H^2] by H-level: H^2 weighted by class mass
n=6; E=[(i,j) for i in range(n) for j in range(i+1,n)]; m=len(E)
from collections import defaultdict
hsum=defaultdict(int)
for bits in range(1<<m):
    adj=[0]*n
    for k,(i,j) in enumerate(E):
        if (bits>>k)&1: adj[i]|=1<<j
        else: adj[j]|=1<<i
    h=Hadj(n,adj); hsum[h]+=h*h    # contribution of each labeled T to sum H^2, grouped by H
tot=sum(hsum.values()); cum=0; hs=sorted(hsum)
print(f"IDEA 6: n={n} contribution to E[H^2] by H-value (is variance carried by HIGH H?):")
top=sorted(hsum.items(),key=lambda x:-x[0])[:6]
for h,c in top: print(f"   H={h:>2}: {100*c/tot:5.1f}% of sum H^2")
hi=sum(c for h,c in hsum.items() if h>=23)  # top half of H-range
print(f"   H>=23 (top half of range): {100*hi/tot:.1f}% of sum H^2  => variance concentrated at the high-H (regular) end: {100*hi/tot>50}")
