import random
from math import gcd
from functools import reduce
from collections import defaultdict
def maxgap_res(res,Vmax):
    pts=sorted(set(res)); n=len(pts); mg=0
    for i in range(n):
        nxt=pts[(i+1)%n]+(Vmax if i==n-1 else 0); g=nxt-pts[i]
        if g>mg: mg=g
    return mg
def jstar(E,Vmax,jmax=500):
    for j in range(1,jmax+1):
        res=[(e*j)%Vmax for e in E]
        if maxgap_res(res,Vmax)*7>Vmax: return j
    return None
def longest_ap(E):
    E=sorted(set(E)); S=set(E); best=1
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]; L=2; nx=E[j]+d
            while nx in S: L+=1; nx+=d
            pv=E[i]-d
            while pv in S: L+=1; pv-=d
            if L>best: best=L
    return best
random.seed(664); k=13
print(f"j* vs longest-AP, co-offset clusters k={k}, spread>=6Vmax/7 (j=1 fails); Vmax=1001:")
Vmax=1001; byap=defaultdict(lambda:[0,0,0])  # [count, max j*, none-found]
lo=int(6*Vmax/7)
for _ in range(3000):
    r=random.random()
    if r<0.4:  # plant near-AP: AP of length L (co-offsets) + fill
        d=random.randint(lo//(k-1)+1, Vmax//(k-1)); L=random.randint(7,k)
        ap=[d*i for i in range(L)]
        if max(ap)>=Vmax: continue
        pool=[x for x in range(1,Vmax) if x not in ap]
        rest=random.sample(pool, k-L)
        E=sorted(set([0]+ap+rest)) if 0 not in ap else sorted(set(ap+rest))
    else:
        spread=random.randint(lo,Vmax-1); E=sorted(set([0]+random.sample(range(1,spread),k-2)+[spread]))
    if len(E)!=k: continue
    if max(E)-min(E)<lo: continue  # ensure large spread
    js=jstar(E,Vmax); lap=longest_ap(E)
    b=byap[lap]; b[0]+=1
    if js is None: b[2]+=1
    else: b[1]=max(b[1],js)
print(f"{'longAP':>7}{'#cfg':>6}{'max j*':>8}{'none-found':>11}")
for lap in sorted(byap):
    b=byap[lap]; print(f"{lap:>7}{b[0]:>6}{b[1]:>8}{b[2]:>11}")
print(f"\n=> if max j* GROWS with longest-AP (and is small/O(1) for low longest-AP), then:")
print("   general j*=O(k) = [high longest-AP: mac-mini Dirichlet AP-lemma] + [low longest-AP: j* small].")
