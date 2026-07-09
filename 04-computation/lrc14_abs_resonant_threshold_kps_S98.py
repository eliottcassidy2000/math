# Spread threshold for the ABSOLUTE resonant bound Sum_{V|n.e}|What| < (6/7)^k (kps-S98).
# Large spread => few resonances => absolute bound holds (a-priori, no cancellation). Small spread =>
# more resonances + cancellation => absolute > 1, needs LEM-013 exhaustion. Find the crossover S*.
# H=8 (converged: increments shrink). Also report signed |R| to show the cancellation gap.
import numpy as np
from math import floor, gcd, sin, pi
from functools import reduce
from itertools import combinations, product
import random
random.seed(98041)
K=13; lead=(6.0/7.0)**K
def babs(m):
    return 0.0 if m%7==0 else abs(sin(pi*m/7))/(pi*abs(m))
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
def abs_res(E,V,H=8,supmax=3):
    E=sorted(E); coords=list(range(1,K)); rng=[m for m in range(-H,H+1) if m!=0 and m%7!=0]; tot=0.0
    for s in range(1,supmax+1):
        for combo in combinations(coords,s):
            for vals in product(rng,repeat=s):
                if sum(vals[t]*E[combo[t]] for t in range(s))%V!=0: continue
                nvec=[0]*(K-1)
                for t,ci in enumerate(combo): nvec[ci-1]=vals[t]
                tot+=What_abs(nvec)
    return tot/lead
def signed_R(E,V):
    Ea=np.array(sorted(E)); j=np.arange(0,V)
    ph=np.mod(np.outer(j,Ea),V)/V; ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return abs(np.maximum(g-1.0/7,0).sum(axis=1).mean()-lead)/lead
print("Absolute resonant bound vs spread (H=8). absRes<1 => a-priori |R|<lead (no cancellation).")
print(f"{'family':>10}{'s':>5}{'V':>5}{'absRes/lead':>12}{'|R|signed/lead':>15}{'<1?':>5}")
def worst_diss(s,Lcap):
    best=(0,None)
    for _ in range(3000):
        mid=sorted(random.sample(range(1,s),K-2)); E=[0]+mid+[s]
        if len(set(E))==K and prim(E) and longest_ap(E)<=Lcap:
            V=floor(7*s/6); a=abs_res(E,V)
            if a>best[0]: best=(a,E)
    return best
for s in (50,70,90,110,140,180,240):
    a,E=worst_diss(s,6)
    if E is None: continue
    V=floor(7*s/6); r=signed_R(E,V)
    print(f"{'dissoc':>10}{s:>5}{V:>5}{a:>12.4f}{r:>15.4f}{'YES' if a<1 else 'NO':>5}")
# 7-structured (all multiples of 7 + shift), several spreads
print("  --- 7-structured (diffs 0 mod 7) ---")
for s in (98,140,210):
    slots=s//7
    if slots<K-1: continue
    best=(0,None)
    for _ in range(2000):
        base=sorted(random.sample(range(1,slots),K-2)); E=sorted(set([0]+[7*b for b in base]+[7*slots]))
        if len(E)!=K: continue
        E2=E[:]; E2[1]+=1; E2=sorted(set(E2))
        if len(E2)!=K or not prim(E2): continue
        V=floor(7*s/6); a=abs_res(E2,V)
        if a>best[0]: best=(a,E2)
    if best[1] is None: continue
    a,E=best; V=floor(7*s/6); r=signed_R(E,V)
    print(f"{'7-struct':>10}{s:>5}{V:>5}{a:>12.4f}{r:>15.4f}{'YES' if a<1 else 'NO':>5}")
print("\n=> crossover spread S*: above it absRes<1 (a-priori, no cancellation); below => LEM-013 exhaustion.")
