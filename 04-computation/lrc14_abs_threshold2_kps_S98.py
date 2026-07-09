# Spread threshold for absolute resonant bound (kps-S98, fast): ~40 dissociated clusters/spread, H=6.
import numpy as np
from math import floor, gcd, sin, pi
from functools import reduce
from itertools import combinations, product
import random
random.seed(98099)
K=13; lead=(6.0/7.0)**K
def babs(m): return 0.0 if m%7==0 else abs(sin(pi*m/7))/(pi*abs(m))
def What_abs(nvec):
    act=[m for m in nvec if m!=0]; r=len(act)
    if any(m%7==0 for m in act): return 0.0
    sig=sum(nvec); pb=1.0
    for m in act: pb*=babs(m)
    if sig==0: sf=6.0/7.0
    elif sig%7==0: return 0.0
    else: sf=babs(sig)
    return (6.0/7.0)**((K-1)-r)*pb*sf
def prim(E):
    E=sorted(E); return reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
def longest_ap(E):
    S=set(E); best=2; E=sorted(E)
    for i in range(len(E)):
        for j in range(i+1,len(E)):
            d=E[j]-E[i]; L=2; nx=E[j]+d
            while nx in S: L+=1; nx+=d
            bk=E[i]-d
            while bk in S: L+=1; bk-=d
            best=max(best,L)
    return best
COORDS=list(range(1,K)); RNG=[m for m in range(-6,7) if m!=0 and m%7!=0]
# precompute (combo, vals, What) once -- independent of E
TERMS=[]
for s in range(1,4):
    for combo in combinations(COORDS,s):
        for vals in product(RNG,repeat=s):
            nvec=[0]*(K-1)
            for t,ci in enumerate(combo): nvec[ci-1]=vals[t]
            w=What_abs(nvec)
            if w>0: TERMS.append((combo,vals,w))
print(f"precomputed {len(TERMS)} nonzero-What terms (support<=3,|n|<=6).")
def abs_res(E,V):
    E=sorted(E); tot=0.0
    for combo,vals,w in TERMS:
        if sum(vals[t]*E[combo[t]] for t in range(len(combo)))%V==0: tot+=w
    return tot/lead
def signed_R(E,V):
    Ea=np.array(sorted(E)); j=np.arange(0,V)
    ph=np.mod(np.outer(j,Ea),V)/V; ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return abs(np.maximum(g-1.0/7,0).sum(axis=1).mean()-lead)/lead
print(f"{'s':>5}{'V':>5}{'max absRes/lead':>16}{'|R|signed at that E':>20}{'absRes<1?':>10}")
for s in (50,70,90,120,160,220):
    V=floor(7*s/6); best=(0,None)
    for _ in range(40):
        mid=sorted(random.sample(range(1,s),K-2)); E=[0]+mid+[s]
        if len(set(E))==K and prim(E) and longest_ap(E)<=6:
            a=abs_res(E,V)
            if a>best[0]: best=(a,E)
    if best[1] is None: continue
    a,E=best; r=signed_R(E,V)
    print(f"{s:>5}{V:>5}{a:>16.4f}{r:>20.4f}{'YES' if a<1 else 'NO':>10}")
print("=> S* = crossover: above, absolute (no-cancellation) bound closes |R|<lead; below, LEM-013 exhaustion.")
