import random
from collections import defaultdict
def maxgap_res(res,Vmax):
    pts=sorted(set(res)); n=len(pts); mg=0
    for i in range(n):
        nxt=pts[(i+1)%n]+(Vmax if i==n-1 else 0); g=nxt-pts[i]
        if g>mg: mg=g
    return mg
def jstar(E,Vmax,jmax=1000):
    for j in range(1,jmax+1):
        if maxgap_res([(e*j)%Vmax for e in E],Vmax)*7>Vmax: return j
    return None
def longest_ap(E):
    E=sorted(set(E)); S=set(E); best=1
    for i in range(len(E)):
        for jj in range(i+1,len(E)):
            d=E[jj]-E[i]; L=2; nx=E[jj]+d
            while nx in S: L+=1; nx+=d
            pv=E[i]-d
            while pv in S: L+=1; pv-=d
            if L>best: best=L
    return best
random.seed(1); k=13
print(f"j* by longest-AP at LARGER Vmax (is small-j* an O(1) bound, not growing with Vmax?), k={k}:")
for Vmax in [1001,5003,20011]:
    lo=int(6*Vmax/7); byap=defaultdict(lambda:[0,0])
    for _ in range(1500):
        r=random.random()
        if r<0.5:
            d=random.randint(lo//(k-1)+1,Vmax//(k-1)); L=random.randint(7,k-2)  # longest-AP <= k-2
            ap=[d*i for i in range(L)]
            if max(ap)>=Vmax: continue
            pool=[x for x in range(1,Vmax) if x not in ap]
            E=sorted(set(([0] if 0 not in ap else [])+ap+random.sample(pool,k-L)))
        else:
            spread=random.randint(lo,Vmax-1); E=sorted(set([0]+random.sample(range(1,spread),k-2)+[spread]))
        if len(E)!=k or max(E)-min(E)<lo: continue
        lap=longest_ap(E)
        if lap>=k-1: continue  # only low-longest-AP here
        js=jstar(E,Vmax); b=byap[lap]; b[0]+=1
        if js: b[1]=max(b[1],js)
        elif js is None: b[1]=max(b[1],9999)
    mx=max((b[1] for b in byap.values()), default=0)
    print(f"  Vmax={Vmax:6d}: longest-AP<=k-2 -> MAX j* = {mx}  (per-AP: {dict(sorted((l,b[1]) for l,b in byap.items()))})")
print("\n=> if MAX j* stays O(1) (~4-6) as Vmax grows 1k->20k, then longest-AP<=k-2 => j*=O(1) is a")
print("   Vmax-uniform bound; combined with the near-AP Dirichlet case => general j*=O(k) DECOMPOSES.")
