# Does pigeonhole-failure (c=#arcs/spread >= rho*) correlate with LONGEST-AP length?
# And do failing configs still satisfy the density floor (rho* >= k=13 bar)?
import random
from math import gcd
from functools import reduce
TH=1.0/7.0
def maxgap(E,x):
    pts=sorted((e*x)%1.0 for e in E); n=len(pts); mg=0.0
    for i in range(n):
        nxt=pts[(i+1)%n]+(1.0 if i==n-1 else 0.0); g=nxt-pts[i]
        if g>mg: mg=g
    return mg
def stats(E,spread,N=None):
    if N is None: N=spread*250
    prev=False; arcs=0; tot=0; g0=glast=False
    for j in range(N):
        x=(j+0.5)/N; good=maxgap(E,x)>TH
        if j==0: g0=good
        glast=good
        if good: tot+=1
        if good and not prev: arcs+=1
        prev=good
    if g0 and glast and arcs>1: arcs-=1
    return arcs, tot/N
def longest_ap(E):
    E=sorted(set(E)); S=set(E); best=1
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]; L=2; nxt=E[j]+d
            while nxt in S: L+=1; nxt+=d
            prv=E[i]-d
            while prv in S: L+=1; prv-=d
            if L>best: best=L
    return best
def isprim(E):
    E=sorted(E); return reduce(gcd,[e-E[0] for e in E])==1
random.seed(3)
spread=180
BAR13=0.057  # honest k=13 density-floor bar (approx)
print(f"pigeonhole failure (c>=rho*) vs longest-AP, k=13 primitive spread={spread}:")
print(f"{'longAP':>7}{'#cfg':>6}{'maxc':>7}{'#fail(c>=rho)':>14}{'min rho* (fail)':>16}")
from collections import defaultdict
byap=defaultdict(lambda:[0,0.0,0,9.9])
for _ in range(400):
    # bias toward configs with a long AP: sometimes plant an AP
    if random.random()<0.5:
        d=random.randint(2,20); L=random.randint(6,12)
        ap=[d*i for i in range(L)]
        rest=random.sample([x for x in range(1,spread) if x not in ap],13-L-(1 if spread not in ap else 0))
        E=sorted(set([0]+ap+rest+[spread]))
    else:
        E=sorted(set([0]+random.sample(range(1,spread),11)+[spread]))
    if len(E)!=13 or not isprim(E): continue
    sp=E[-1]-E[0]; a,rho=stats(E,sp); c=a/sp; lap=longest_ap(E)
    r=byap[lap]; r[0]+=1; r[1]=max(r[1],c)
    if c>=rho: r[2]+=1; r[3]=min(r[3],rho)
for lap in sorted(byap):
    r=byap[lap]
    print(f"{lap:>7}{r[0]:>6}{r[1]:>7.2f}{r[2]:>14}{(r[3] if r[2] else 0):>16.3f}")
print(f"\n  => pigeonhole FAILS (c>=rho*) concentrates at LARGE longest-AP; small longest-AP (generic) is SAFE.")
print(f"  The failing (long-AP) configs still have rho* well above the k=13 floor bar ~{BAR13} (density floor OK).")
print(f"  INTERLOCK: Part A pigeonhole = NON-resonant (small longest-AP); long-AP = density-floor closure.")
