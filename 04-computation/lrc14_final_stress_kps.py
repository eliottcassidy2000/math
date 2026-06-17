#!/usr/bin/env python3
"""
lrc14_final_stress_kps  (kind-pasteur 2026-06-17) -- final adversarial stress.
 [A] wider tight census: exhaustive 13-subsets of [1,24] containing 1, but pruned:
     require the set to contain >=10 of {1..13} (tight needs heavy AP overlap; the
     12-core analysis shows you cannot deviate much). Counts ALL tight found.
 [B] sporadic-window: 13-subsets that are {1..11,13} ∪ {w} for w up to 400 -> tight w?
     (these are EXACTLY the configs that can be tight per the 12-core coverage logic)
 [C] two-stranger L-minimization: {1..11,13} drop one more, plus 2 strangers up to 200,
     exact L on float-screened survivors. Beat 1/1260?
 [D] does L -> 0 along ANY lcm-growing sequence? track min L per lcm-bucket.
"""
import sys, io
try: sys.stdout=io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
except Exception: pass
from fractions import Fraction as Fr
from itertools import combinations
from math import gcd
from functools import reduce

def arcs_of(v):
    w=Fr(1,14*v); A=[]
    for k in range(v+1):
        c=Fr(k,v); lo=c-w; hi=c+w
        if lo<0: A+=[(Fr(0),hi),(1+lo,Fr(1))]
        elif hi>1: A+=[(lo,Fr(1)),(Fr(0),hi-1)]
        else: A.append((lo,hi))
    return A
def L_exact(S):
    A=[]
    for v in set(S): A.extend(arcs_of(v))
    A=sorted((a,b) for a,b in A if b>a)
    tot=Fr(0); cl=ch=None
    for a,b in A:
        if ch is None: cl,ch=a,b
        elif a<=ch: ch=max(ch,b)
        else: tot+=ch-cl; cl,ch=a,b
    if ch is not None: tot+=ch-cl
    return 1-tot
def L_float(S):
    arcs=[]
    for v in set(S):
        inv=1.0/(14*v)
        for k in range(v+1):
            lo=(14*k-1)*inv; hi=(14*k+1)*inv
            if lo<0.0: arcs.append((0.0,hi)); arcs.append((1.0+lo,1.0))
            elif hi>1.0: arcs.append((lo,1.0)); arcs.append((0.0,hi-1.0))
            else: arcs.append((lo,hi))
    arcs=[(a,b) for a,b in arcs if b>a]; arcs.sort()
    tot=0.0; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch:
            if hi>ch: ch=hi
        else: tot+=ch-cl; cl,ch=lo,hi
    tot+=ch-cl
    return 1.0-tot
def lcm(S): return reduce(lambda a,b:a*b//gcd(a,b),S,1)
def is_prim(S): return reduce(gcd,S)==1

AP=set(range(1,14))
print("="*78)
print("[A] tight census: 13-subsets of [1,24] with >=10 of {1..13}, contain 1")
print("="*78)
pool=list(range(1,25)); tightA=[]; cnt=0
# choose 13 from [1,24]; prune by requiring >=10 overlap with {1..13}
for combo in combinations(range(2,25),12):
    S=(1,)+combo
    if len(set(S)&AP)<10: continue
    cnt+=1
    if L_float(S)>1e-9: continue
    if L_exact(S)==0: tightA.append(S)
print("tested %d; tight: %d"%(cnt,len(tightA)))
for S in sorted(tightA):
    print("   lcm=%d max=%d prim=%s %s"%(lcm(S),max(S),is_prim(S),S))

print("\n"+"="*78)
print("[B] single-stranger tightness: {1..11,13} + w, w=1..400 -> which w tight?")
print("="*78)
core=[1,2,3,4,5,6,7,8,9,10,11,13]; tightB=[]
for w in range(1,401):
    if w in core: continue
    S=tuple(sorted(core+[w]))
    if L_float(S)>1e-9: continue
    if L_exact(S)==0: tightB.append(w)
print("tight w (covering all 4 gaps of {1..11,13}):", tightB)
# also for each AP-drop-j core, which single w tightens it?
print("for each AP-drop-j core, the tightening strangers w<=400:")
for j in range(1,14):
    cj=[x for x in range(1,14) if x!=j]; ws=[]
    for w in range(1,401):
        if w in cj: continue
        S=tuple(sorted(cj+[w]))
        if L_float(S)>1e-9: continue
        if L_exact(S)==0: ws.append(w)
    print("   drop %2d: tight strangers = %s"%(j,ws))

print("\n"+"="*78)
print("[C] two-stranger L-minimization: {1..13} drop j1,j2 + two strangers <=200; beat 1/1260?")
print("="*78)
TARGET=1/1260.0; bestC=(1.0,None)
hitsC=[]
for j1,j2 in combinations(range(1,14),2):
    base=[x for x in range(1,14) if x not in (j1,j2)]  # 11 entries
    WR=list(range(14,201))
    for w1,w2 in combinations(WR,2):
        if w1 in base or w2 in base: continue
        S=base+[w1,w2]
        Lf=L_float(S)
        if Lf<bestC[0]: bestC=(Lf,tuple(sorted(S)))
        if Lf<TARGET-1e-12: hitsC.append((Lf,tuple(sorted(S)),(j1,j2,w1,w2)))
print("best two-stranger float L:", bestC[0], bestC[1])
print("configs with float L < 1/1260:", len(hitsC))
for Lf,S,info in sorted(hitsC)[:20]:
    Le=L_exact(list(S))
    print("   L=%s=%.9g (exact)  %s  drop/add %s"%(Le,float(Le),S,info))
if not hitsC:
    print("  NONE below 1/1260 (exact-confirm of float-best):")
    Le=L_exact(list(bestC[1]))
    print("   float-best exact L =", Le, "=", float(Le))

print("\n"+"="*78)
print("[D] L vs lcm: single-stranger {1..11,13}+w sweep, min L per lcm decade")
print("="*78)
buckets={}
for w in range(14,2001):
    if w in core: continue
    S=core+[w]
    if not is_prim(S): continue
    Le=L_exact(S)
    if Le==0: continue
    lc=lcm(S); dec=len(str(lc))
    if dec not in buckets or Le<buckets[dec][0]:
        buckets[dec]=(Le,w,lc)
for dec in sorted(buckets):
    Le,w,lc=buckets[dec]
    print("   lcm~1e%d: min L=%s=%.9g at w=%d (lcm=%d)"%(dec-1,Le,float(Le),w,lc))
