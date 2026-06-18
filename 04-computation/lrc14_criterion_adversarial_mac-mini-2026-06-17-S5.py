#!/usr/bin/env python3
"""
lrc14_criterion_adversarial — mac-mini-2026-06-17-S5

Stress the CRITERION C(S)=[exists v: W(S\{v})>1/(7v)] => M(S)>=1/14 on the HARD regime:
covering 13-sets with MANY runners clustered LARGE (where the pigeonhole bound
W>=mu/Sum(u) weakens, since Sum(u)~12V). Two questions:
 (1) does EXACT C(S) still hold (some v works)?  -- if yes, the criterion is robust.
 (2) does pigeonhole-via-largest 7V*mu_V > Sum_{u<V} u hold, or do we need exact W?
Build clustered covering sets: in a window [N,N+win], greedily pick a distinct multiple
of each q in 2..14 (13 moduli -> 13 runners, all large). Also: scaled small-cores.
"""
from fractions import Fraction as F
from itertools import combinations
import random

C=F(1,14)
def darcs(v,c=C):
    hw=F(c,v); return [(F(k,v)-hw,F(k,v)+hw) for k in range(v)]
def wrapU(iv):
    o=[]
    for lo,hi in iv:
        s=lo-(lo%1); a=lo-s;b=hi-s
        if b<=1:o.append((a,b))
        else:o.append((a,F(1)));o.append((F(0),b-1))
    o=sorted(o);r=[];cl,ch=o[0]
    for lo,hi in o[1:]:
        if lo<=ch: ch=ch if ch>hi else hi
        else:r.append((cl,ch));cl,ch=lo,hi
    r.append((cl,ch));return r
def Wsafe(A,c=C):
    dz=[]
    for v in set(A): dz+=darcs(v,c)
    if not dz: return F(1)
    dz=wrapU(dz); best=F(0)
    for i in range(len(dz)):
        hi=dz[i][1]; lo=dz[(i+1)%len(dz)][0]+(1 if i==len(dz)-1 else 0)
        if lo-hi>best: best=lo-hi
    return best
def covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def crit(S):
    res=[]
    for v in sorted(set(S)):
        A=[u for u in S if u!=v]; W=Wsafe(A); thr=F(1,7*v)
        res.append((v,W,thr,W>thr))
    any_ok=any(r[3] for r in res)
    best=max(res,key=lambda r:float(r[1]-r[2]))
    return any_ok,best,res
def pigeon_largest(S):
    V=max(S); core=[u for u in S if u!=V]
    mu=F(1)-sum(F(1,7*u) for u in set(core)); su=sum(core)
    return (7*V)*mu> su, mu, su, V
# exact M
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r; return r if r<=F(1,2) else 1-r
def g(S,t): return min(nrm(v*t) for v in S)
def cand(S):
    S=sorted(set(S)); Cc=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): Cc.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): Cc.add(F(k,d)); k+=1
    Cc.add(F(1,2)); return Cc
def M(S): return max(g(S,t) for t in cand(S))

def clustered_covering(N, win, rng):
    used=set(); S=[]
    for q in range(2,15):
        cands=[x for x in range(N, N+win+1) if x%q==0 and x not in used]
        if not cands: return None
        x=rng.choice(cands); used.add(x); S.append(x)
    S=sorted(set(S))
    return S if len(S)==13 and covering(S) else None

print("="*78)
print("ADVERSARIAL: covering 13-sets with ALL runners clustered LARGE")
print("="*78)
rng=random.Random(5); ntot=0; cfail=0; pigfail_cfok=0; mbreak=0; worst=(F(99),None)
examples=[]
for _ in range(6000):
    N=rng.choice([200,500,1000,3000,10000]); win=rng.choice([40,80,160,300])
    S=clustered_covering(N,win,rng)
    if S is None: continue
    ntot+=1
    ok,best,_=crit(S)
    pig,mu,su,V=pigeon_largest(S)
    if not ok:
        cfail+=1
        if M(S)<C: mbreak+=1
        examples.append(('CFAIL',S))
    elif not pig:
        pigfail_cfok+=1
        if len(examples)<3: examples.append(('pigeon-fails-but-exact-C-ok',S,best[0]))
    m=best[1]-best[2]
    if ok and m<worst[0]: worst=(m,S,best[0])
    if ntot>=1500: break
print(f"  clustered covering 13-sets tested: {ntot}")
print(f"  EXACT C(S) failed: {cfail}  (LRC breaks: {mbreak})")
print(f"  pigeonhole-via-largest failed but exact C held: {pigfail_cfok}")
print(f"  tightest exact margin W-1/(7v): {float(worst[0]):.7f} at v={worst[2] if worst[1] else '-'}")
for e in examples[:4]: print("   eg:",e[0], e[1][:6],"...len",len(e[1]))
print()
print("READING: if EXACT C holds on all clustered sets, the criterion is robust (proof target).")
print("If pigeonhole-via-largest fails on some, we need a SHARPER W lower bound there")
print("(three-distance / the actual safe-arc structure), not the crude mu/Sum(u).")
