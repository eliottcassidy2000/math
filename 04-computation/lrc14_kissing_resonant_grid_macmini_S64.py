"""
mac-mini-2026-07-09-S64 -- EXTEND kps-S97's kissing-number route: stress-test the "uniform in Vmax"
gap on the RESONANT grid 7|Vmax.

kps-S97: E_grid[W] = (6/7)^k + R, R = Poisson tail over L_V(e)={n: V|n.e}; |R| ~ c*kissing(L_V);
a-priori route |R| <= c*kissing(AP) = 0.61*lead < lead => existence. FLAGGED GAP: prove
kissing(L_V(e)) <= kissing(L_V(AP)) *uniformly in Vmax*.

MY ANGLE (HYP-5600): the resonant grid 7|V gives L_V(e) EXTRA short vectors (the mod-7 resonances),
so 7-structured sets could have HIGHER |R|/lead at 7|V than kps's non-resonant V=floor(7s/6) measured.
This is precisely where "AP maximizes kissing uniformly in V" could FAIL. Question:
  (Q1) sup_{7-struct, 7|V} |R|/lead -- does it exceed kps's AP value 0.61?  Does it stay < 1 (existence)?
  (Q2) how does it compare to gcd(7,V)=1 for the SAME sets?  (the mod-7 split)
  (Q3) if |R|/lead can exceed the AP's, the clean "AP-extremal" bound needs the mod-7 refinement --
       but existence (E_grid[W]>0) is SAVED by the resonance's own large positive W spike at x=m/7.
Decisive: report sup |R|/lead per grid-class + whether E_grid[W] ever <= 0.
"""
import numpy as np
from math import gcd, floor
from functools import reduce
import random
random.seed(6400)
K=13; lead=(6.0/7.0)**K

def prim(E):
    E=sorted(E); return len(E)>=2 and reduce(gcd,[E[i+1]-E[i] for i in range(len(E)-1)])==1
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
def Egrid_W(E,V):
    """mean over j=0..V-1 of W(j/V) = (6/7)^k + R exactly (float)."""
    Ea=np.array(sorted(E)); j=np.arange(0,V)
    ph=np.mod(np.outer(j,Ea),V)/V; ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    return np.maximum(g-1.0/7,0).sum(axis=1).mean()
def make_7struct(k, s, min_s7):
    sevens=[x for x in range(0,s+1,7)]
    if len(sevens)<min_s7: return None
    S7=sorted(random.sample(sevens, random.randint(min_s7, min(len(sevens),k-1))))
    rest_pool=[x for x in range(1,s) if x%7!=0]
    need=k-len(S7)
    if need<0 or len(rest_pool)<need: return None
    E=sorted(set(S7+[0,s]+random.sample(rest_pool,max(0,need))))
    return E if len(E)==k else None

# AP baseline (kps: |R|/lead=0.61, max kissing)
AP=[i*7 for i in range(K)]   # step-7 AP = the maximal-energy 7-structured object (its own longest-AP=13)
print(f"K={K}, lead=(6/7)^K={lead:.6f}")
print(f"AP step-7 {AP[:5]}...: longest-AP={longest_ap(AP)}")
for V in (91, 98, floor(7*max(AP)/6)):
    ew=Egrid_W(AP,V); print(f"   AP V={V} ({'7|V' if V%7==0 else 'gcd=1'}): E_grid[W]={ew:.5f}  |R|/lead={abs(ew-lead)/lead:.4f}")

print(f"\n{'grid-class':>14}{'#sets':>7}{'sup|R|/lead':>13}{'min E_grid[W]':>14}{'#(<=0)':>8}   worst-set")
for label, cond in (("7|Vmax", "res"), ("gcd(7,V)=1", "cop")):
    sup=0.0; mn=9.9; nzero=0; n=0; worst=None
    for _ in range(6000):
        s=random.randint(40,160); E=make_7struct(K,s,4)
        if E is None or not prim(E) or longest_ap(E)>K-6: continue
        mx=max(E)
        # sweep several V of the chosen class, V>mx
        Vs=[]
        base=mx+1
        for t in range(1,10):
            if cond=="res":
                V=7*((mx)//7 + t)
            else:
                V=base+t
                while V%7==0: V+=1
            if V>mx: Vs.append(V)
        for V in Vs[:8]:
            ew=Egrid_W(E,V); r=abs(ew-lead)/lead
            n+=1
            if ew<=1e-12: nzero+=1
            if ew<mn: mn=ew
            if r>sup: sup=r; worst=(E[:7],V,round(ew,5),round(r,4))
    print(f"{label:>14}{n:>7}{sup:>13.4f}{mn:>14.5f}{nzero:>8}   {worst}")

print("\nINTERPRETATION:")
print(" - sup|R|/lead vs AP's 0.61: if 7|V exceeds it => grid-kissing is NOT AP-extremal at resonance")
print("   (kps's flagged 'uniform in Vmax' gap is real at 7|V), but...")
print(" - min E_grid[W] > 0 and #(<=0)=0 => EXISTENCE still holds: the resonance's own W-spike at x=m/7")
print("   keeps the first moment positive even when |R|/lead is large. |R|<lead is the true threshold.")
