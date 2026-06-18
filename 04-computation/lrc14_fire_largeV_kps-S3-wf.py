"""
Independent confirmation of 'Lemma A fires => witness safe & M>=1/14' at LARGE V,
on a small set of explicitly constructed covering S3 13-sets, plus a sanity check
that the produced witness is exactly 1/14-safe for ALL 13 speeds.
Also probes large-V for any M<1/14 on a modest number of sets (Mval is slow at big V,
so few sets but each exact).
"""
import sys
from fractions import Fraction as F
from math import gcd
import random

def flush(*a): print(*a); sys.stdout.flush()
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def safe_components(A,h=F(1,14)):
    iv=[]
    for u in A:
        for j in range(0,u):
            c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
            if a<b: iv.append((a,b))
            else: iv.append((a,F(1))); iv.append((F(0),b))
    iv.sort(); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
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
def Mval(S):
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b
def is_covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def gcd_list(xs):
    g=0
    for x in xs: g=gcd(g,x)
    return g
H=F(1,14)
def lemmaA_window(K,Vmin,Vmax): return F(14*K+1,14*Vmin), F(14*K+13,14*Vmax)
def lemmaA_fires(P,L):
    Vmin=min(L); Vmax=max(L); Ps=safe_components(P,H); s=Vmax-Vmin
    if s==0: Kmax=2 if (13*Vmin-Vmax)>0 else -1
    else:
        if 13*Vmin-Vmax<=0: return None
        Kmax=(13*Vmin-Vmax-1)//(14*s)
    for K in range(0,Kmax+1):
        lo,hi=lemmaA_window(K,Vmin,Vmax)
        if not(lo<hi): continue
        for (a,b) in Ps:
            ov=(max(lo,a),min(hi,b))
            if ov[0]<ov[1]: return (ov[0]+ov[1])/2
    return None

def gen(n=120, seed=11, maxV=12000):
    rnd=random.Random(seed); out=[]; att=0
    while len(out)<n and att<n*500:
        att+=1
        psz=rnd.randint(2,9); P=sorted(rnd.sample(range(1,14),psz))
        csz=rnd.randint(2,max(2,13-len(P)))
        V0=rnd.randint(500,maxV); spread=rnd.choice([1,2,3,5,8,14,25,45])
        L=set()
        while len(L)<csz*4: L.add(V0+rnd.randint(0,spread))
        L=sorted(L)[:csz]
        if len(L)<2: continue
        S=sorted(set(P)|set(L))
        if len(S)<13: continue
        S=S[:13]
        if len(set(S))!=13 or gcd_list(S)!=1 or not is_covering(S): continue
        if sum(1 for v in S if v>13)<2: continue
        if max(S)<13*min(S): continue
        out.append(S)
    return out

if __name__=="__main__":
    sets=gen(n=120)
    flush("large-V covering S3 sets:", len(sets))
    fired=0; fire_fail=0; below=0; mn=None
    for S in sets:
        P=[v for v in S if v<=13]; L=[v for v in S if v>13]
        w=lemmaA_fires(P,L)
        m=Mval(S)
        if mn is None or m<mn: mn=m
        if m<H: below+=1; flush("  M<1/14 COUNTEREXAMPLE:", S, m)
        if w is not None:
            fired+=1
            wsafe=all(nrm(v*w)>=H for v in S)
            if (not wsafe) or m<H:
                fire_fail+=1; flush("  FIREFAIL:", S, float(w), wsafe, m)
    flush(f"fired={fired}  fireFail={fire_fail}  below1/14={below}  minM={mn} ({float(mn) if mn else None})  Vmax up to 12000")
    flush("1/14 =", float(H))
