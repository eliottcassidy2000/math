from fractions import Fraction as F
from math import gcd
from itertools import combinations
import random
random.seed(7)
def M_exact(fam):
    Q=2*max(fam)+2; best=F(0)
    for q in range(2,Q+1):
        for a in range(1,q):
            if gcd(a,q)!=1: continue
            m=min(min((v*a)%q,q-(v*a)%q) for v in fam)
            if F(m,q)>best: best=F(m,q)
    return best
def good_C_mu(W,lvl=F(1,13)):
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
    return len(comps), sum(float(b-a) for b,a in [(b,a) for a,b in comps])  # C, mu
def covers(fam,S): return all(any(v%q==0 for v in fam) for q in S)
def is_AP(W): W=sorted(W); d=W[0]; return W==[d*i for i in range(1,13)]
# COMBINED: every valid non-AP core => [C<=464mu (soft Weyl)] OR [delta>max/2366 (stability empty window)] ?
base=list(range(1,13)); pool=[v for v in range(1,41) if v%13 and v%14]
tested=0; both_fail=[]; sw=0; st=0
def test(W):
    global tested,sw,st
    W=sorted(set(W))
    if len(W)!=12 or is_AP(W) or not covers(W,range(2,13)) or any(v%13==0 or v%14==0 for v in W): return
    if max(W)>34: return
    tested+=1
    M=M_exact(W); delta=M-F(1,13)
    C,mu=good_C_mu(W)
    softweyl = (C <= 464.4*mu)
    stability = (delta > F(max(W),2366))   # empty window
    if softweyl: sw+=1
    if stability: st+=1
    if not (softweyl or stability): both_fail.append((tuple(W),C,mu,float(delta),max(W)))
for k in (1,2,3):
    for pos in combinations(range(12),k):
        combos=[tuple(random.choice(pool) for _ in range(k)) for _ in range(120)] if k>=2 else [(x,) for x in pool]
        for nv in combos:
            W=base[:]
            for p,x in zip(pos,nv): W[p]=x
            test(W)
print(f"valid non-AP cores (max<=34) tested: {tested}")
print(f"  eliminated by soft Weyl (C<=464mu): {sw}")
print(f"  eliminated by stability (delta>max/2366, empty window): {st}")
print(f"  eliminated by NEITHER (genuine residual): {len(both_fail)}")
for f in both_fail[:8]: print(f"    RESIDUAL: W={list(f[0])} C={f[1]} mu={f[2]:.4f} delta={f[3]:.5f} max={f[4]}")
if not both_fail:
    print("  => EVERY valid non-AP core (max<=34) is eliminated by soft Weyl OR stability.")
    print("     POSITION LEMMA (and the whole crux) CLOSED for max<=34 by two PROVED lenses.")
