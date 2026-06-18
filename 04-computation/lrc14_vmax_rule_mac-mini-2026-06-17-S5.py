#!/usr/bin/env python3
"""lrc14_vmax_rule — does v=largest satisfy W(S\{V})>1/(7V)? vs general exists-v. (W only, fast)"""
from fractions import Fraction as F
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
    dz=wrapU(dz); best=F(0)
    for i in range(len(dz)):
        hi=dz[i][1]; lo=dz[(i+1)%len(dz)][0]+(1 if i==len(dz)-1 else 0)
        if lo-hi>best: best=lo-hi
    return best
def covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def vmax_ok(S):
    V=max(S); return Wsafe([u for u in S if u!=V])>F(1,7*V)
def anyv_ok(S):
    return any(Wsafe([u for u in S if u!=v])>F(1,7*v) for v in set(S))

def clustered(N,win,rng):
    used=set(); S=[]
    for q in range(2,15):
        cs=[x for x in range(N,N+win+1) if x%q==0 and x not in used]
        if not cs: return None
        x=rng.choice(cs); used.add(x); S.append(x)
    S=sorted(set(S)); return S if len(S)==13 and covering(S) else None

rng=random.Random(11)
fams={'spread (drop+84k)':[], 'clustered':[], 'mixed (small+2big)':[]}
# spread: reliable covering
for drop in [1,2,3,4,5,6,7,12]:
    for k in range(1,30):
        S=sorted(set([v for v in range(1,14) if v!=drop]+[84*k]))
        if len(S)==13 and covering(S): fams['spread (drop+84k)'].append(S)
# clustered
for _ in range(800):
    S=clustered(rng.choice([60,150,400,1200,5000]), rng.choice([30,70,150,300]), rng)
    if S: fams['clustered'].append(S)
# mixed: small spread core + two large covering elements
for _ in range(3000):
    drop=rng.choice([1,2,3,4,5,6]); base=[v for v in range(1,14) if v!=drop]
    d2=rng.choice([8,9,10,11,13]); base=[v for v in base if v!=d2]
    bigs=set()
    while len(bigs)<2: bigs.add(rng.choice([84,168,126,210,154,182,84,252])*rng.randint(1,5))
    S=sorted(set(base+list(bigs)))
    if len(S)==13 and covering(S): fams['mixed (small+2big)'].append(S)

print("="*72)
print("Does v=LARGEST satisfy W(S\\{V})>1/(7V)?  vs general exists-v  (covering 13-sets)")
print("="*72)
gvf=0; gtot=0; worst=(F(99),None)
for name,L in fams.items():
    L=L[:1200]; vf=sum(0 if vmax_ok(S) else 1 for S in L); af=sum(0 if anyv_ok(S) else 1 for S in L)
    for S in L:
        V=max(S); m=Wsafe([u for u in S if u!=V])-F(1,7*V)
        if m<worst[0]: worst=(m,S)
    print(f"  {name:22s}: n={len(L):4d}  v=largest FAILS: {vf:3d}  exists-v FAILS: {af:3d}")
    gvf+=vf; gtot+=len(L)
print(f"\n  TOTAL: v=largest fails {gvf}/{gtot};  worst v=largest margin {float(worst[0]):+.6f}")
if worst[1]: print(f"    at S={worst[1][:7]}... (V={max(worst[1])})")
print("\n  => if v=largest fails sometimes but exists-v never, the proof needs the case-split;")
print("     if v=largest NEVER fails, the deterministic rule v=max(S) is the whole proof.")
