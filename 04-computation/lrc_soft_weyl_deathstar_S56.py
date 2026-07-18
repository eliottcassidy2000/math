from fractions import Fraction as F
from math import gcd, pi
from itertools import combinations
import random
random.seed(1)
def good_components(W, lvl=F(1,13)):
    ivs=[]
    for w in W:
        for k in range(0,w+1):
            lo=(F(k)-lvl)/w; hi=(F(k)+lvl)/w
            ivs.append((max(lo,F(0)),min(hi,F(1))))
    ivs=[(a,b) for a,b in ivs if b>a]; ivs.sort()
    m=[]
    for a,b in ivs:
        if m and a<=m[-1][1]: m[-1]=(m[-1][0],max(m[-1][1],b))
        else: m.append((a,b))
    comps=[]; prev=F(0)
    for a,b in m:
        if a>prev: comps.append((prev,a))
        prev=max(prev,b)
    if prev<F(1): comps.append((prev,F(1)))
    return comps
def covers(fam,S): return all(any(v%q==0 for v in fam) for q in S)
def is_AP(W): W=sorted(W); d=W[0]; return W==[d*i for i in range(1,13)]
# soft Weyl: position lemma holds (avg>=1/13) if C <= 464*k*mu  (k=1 binding).  bound C<=464*mu.
# verify for valid non-AP cores.
base=list(range(1,13)); pool=[v for v in range(1,36) if v%13 and v%14]
tested=0; ok=0; fails=[]; worst=None
def test(W):
    global tested,ok,worst
    W=sorted(set(W))
    if len(W)!=12 or is_AP(W) or not covers(W,range(2,13)) or any(v%13==0 or v%14==0 for v in W): return
    tested+=1
    comps=good_components(W); C=len(comps); mu=sum(float(b-a) for a,b in comps)
    bound=464.4*mu
    ratio=C/bound if bound>0 else 999
    if worst is None or ratio>worst[0]: worst=(ratio,tuple(W),C,mu)
    if C<=bound: ok+=1
    else: fails.append((tuple(W),C,mu,bound))
for k in (1,2):
    for pos in combinations(range(12),k):
        combos=[tuple(random.choice(pool) for _ in range(k)) for _ in range(150)] if k==2 else [(x,) for x in pool]
        for nv in combos:
            W=base[:]
            for p,x in zip(pos,nv): W[p]=x
            test(W)
print(f"valid non-AP cores tested: {tested}")
print(f"  C <= 464*mu (soft Weyl endpoint bound => position lemma PROVED for this core): {ok}/{tested}")
print(f"  fails: {len(fails)}")
r,W,C,mu=worst
print(f"  tightest ratio C/(464mu)={r:.3f} at W={list(W)}: C={C} components, mu={mu:.4f}, 464mu={464.4*mu:.2f}")
for f in fails[:6]: print(f"    FAIL W={list(f[0])}: C={f[1]} > 464mu={f[3]:.2f} (mu={f[2]:.4f})")
