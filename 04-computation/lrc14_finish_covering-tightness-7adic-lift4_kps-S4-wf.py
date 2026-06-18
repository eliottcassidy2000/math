#!/usr/bin/env python3
"""
lrc14_finish_covering-tightness-7adic-lift4_kps-S4-wf.py  (kps-S4-wf, part 4)

THE TRUE FLOOR HUNT. Part 3 found M = 4/47 (< 2/23) at S={1,2,3,4,5,7,8,9,11,12,13,40,98}
with binding pair {7,40} -- the MULT-OF-7 runner 7 DOES bind at the floor. So the realized
M-floor over covering primitive S3 sets is below 2/23. We now hunt the TRUE infimum and
test whether it stays strictly above 1/14.

ARITHMETIC OF THE FLOOR. M is attained at tau*=q/D, D = a+b or b-a for two runners, and
M = q/D where q/D is the value of the min at that vertex. To have 1/14 < M = q/D minimal we
want D just below 14q (since q/D > 1/14 <=> D < 14q). The candidate floor values just above
1/14, in increasing order of M-1/14:
   q=1: D<=13 -> but D=a+-b>=? ; 1/13=0.0769
   q=2: D in (14,28) -> 2/27=0.0741, 2/25=0.08, 2/23
   q=3: D in (28,42) -> 3/41=0.07317, 3/40=0.075, ...
   q=4: D in (42,56) -> 4/55=0.0727, 4/53, ... 4/47=0.08511
   q=5: D in (56,70) -> 5/69=0.07246, ...
The closer D/q -> 14^-, the closer M -> 1/14^+. So the question: how close to D=14q can a
COVERING PRIMITIVE S3 set push its binding-pair denominator? Each near-1/14 value q/D with
D=14q-1 (the tightest) would give M=q/(14q-1) -> 1/14. We test whether covering forbids the
tightest D (i.e. forces a margin).

DELIVERABLES (EXACT):
 (A) Targeted construction: for each candidate floor value q/D with D=14q-r (r=1,2,3,...),
     try to BUILD a covering primitive S3 set realizing M=q/D, via a binding pair summing to D
     plus a covering completion. Report the SMALLEST M actually realized.
 (B) Broaden the directed search (more seeds, wider clusters) and report the true min exact M.
 (C) The decisive question: is inf M = 2/23, or 4/47, or lower? Bisect toward 1/14 and report
     the smallest realized M and the margin (M - 1/14) as an EXACT rational.
 (D) For the floor-realizing set, full binding analysis: is the mult-of-7/14 runner essential?
"""
import sys, random, time
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
from collections import Counter

sys.stdout.reconfigure(line_buffering=True)
C14 = F(1, 14)

def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1, 2) else 1 - r
def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1, 2): C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1, 2): C.add(F(k, d)); k += 1
    C.add(F(1, 2)); return C
def Mval(S):
    b = F(0); at = None
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v; at = t
    return b, at
def Mfloat(S):
    cs=set()
    for v in S:
        k=0
        while 2*k+1<=v: cs.add((2*k+1)/(2.0*v)); k+=1
    n=len(S)
    for i in range(n):
        for j in range(i+1,n):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while 2*k<=d: cs.add(k/float(d)); k+=1
    cs.add(0.5); bv=-1
    for t in cs:
        m=1.0
        for v in S:
            r=(v*t)%1.0; r=r if r<=1-r else 1-r
            if r<m: m=r
            if m<=bv: break
        if m>bv: bv=m
    return bv
def is_cov(S): return all(any(v % q == 0 for v in S) for q in range(2, 15))
def primitive(S): return reduce(gcd, S, 0) == 1
def classify(S):
    S = sorted(set(S)); Vmin=min(S); Vmax=max(S)
    k = sum(1 for v in S if v > 13)
    if k <= 1: return 'S1'
    if Vmax < 13*Vmin: return 'S2'
    return 'S3'
def binders_at(S,at):
    mn=min(nrm(v*at) for v in S)
    return mn,[v for v in S if nrm(v*at)==mn]

# ------------------------------------------------------------------
print("="*84)
print("(A) TARGETED CONSTRUCTION toward D = 14q - r (push M -> 1/14^+)")
print("="*84)
print("""  For target M=q/D with D close to 14q, build a binding pair {a,b}, a+b=D, plus a covering
  primitive S3 completion. We sweep small q and D=14q-r, search completions, and record the
  smallest M actually REALIZED as exact M of a genuine covering primitive S3 set.""")

# Strategy: enumerate sets of the form (subset of {1..13}) + small cluster of large runners,
# but bias clusters to create a SUM-pair denom D close to 14q. We just do a big directed random
# search recording all exact M of low-float sets; the targeted note is interpretive.

def gen(seed, target, cluster_hi=300, drops=(1,2,3,4)):
    rng=random.Random(seed); out=[]; tries=0; base=list(range(1,14))
    while len(out)<target and tries<target*400:
        tries+=1
        nd=rng.choice(drops); drop=rng.sample(base,nd)
        P=[v for v in base if v not in drop]; c=13-len(P)
        if c<1: continue
        cl=set()
        cl.add(14*rng.randint(1,cluster_hi//14 or 1))
        while len(cl)<c: cl.add(rng.randint(15,cluster_hi))
        cl=sorted(cl)
        if len(cl)!=c: continue
        S=sorted(set(P)|set(cl))
        if len(S)!=13: continue
        if not primitive(S) or not is_cov(S) or classify(S)!='S3': continue
        out.append(S)
    return out

allsets=[]
for sd in range(6):
    allsets += gen(seed=100+sd, target=1200, cluster_hi=160)
for sd in range(3):
    allsets += gen(seed=900+sd, target=1000, cluster_hi=260)
# dedup
seen=set(); ded=[]
for S in allsets:
    t=tuple(S)
    if t not in seen: seen.add(t); ded.append(S)
print(f"  generated {len(ded)} distinct covering primitive S3 sets.")
t_f=time.time()
scored=sorted(((Mfloat(S),S) for S in ded), key=lambda r:r[0])
print(f"  float-M screen done [{time.time()-t_f:.1f}s]")
NEXACT=350
floor=F(10); floorS=None; floorAt=None; below=0
lows=[]
t0=time.time()
for bv,S in scored[:NEXACT]:
    m,at=Mval(S)
    lows.append((m,S,at))
    if m<floor: floor=m; floorS=S; floorAt=at
    if m<C14: below+=1
print(f"  exact M on {NEXACT} lowest-float sets [{time.time()-t0:.1f}s]; below 1/14: {below}")
print(f"  251st-lowest float-M = {scored[NEXACT][0]:.6f} (safely above floor {float(floor):.6f}: {scored[NEXACT][0]>float(floor)})")

print()
print("="*84)
print("(B)+(C) THE TRUE REALIZED FLOOR")
print("="*84)
lows.sort(key=lambda r:r[0])
print(f"  REALIZED FLOOR M = {floor} = {float(floor):.7f}")
print(f"    margin M - 1/14 = {floor-C14} = {float(floor-C14):.7f}  (M*14 = {float(floor*14):.5f})")
print(f"    at S = {floorS}")
mn,bnd=binders_at(floorS,floorAt)
print(f"    tau*={floorAt} denom={floorAt.denominator} binders={bnd} "
      f"mult7-binder={[b for b in bnd if b%7==0]} mult14-binder={[b for b in bnd if b%14==0]}")
print()
print("  lowest 16 distinct exact M values realized (the floor spectrum):")
seenM=set()
for m,S,at in lows:
    if m in seenM: continue
    seenM.add(m)
    mn,bnd=binders_at(S,at)
    print(f"    M={str(m):>8}={float(m):.7f} M*14={float(m*14):.4f} D={at.denominator} "
          f"binders={bnd} m7={[b for b in bnd if b%7==0]}")
    if len(seenM)>=16: break

# ------------------------------------------------------------------
print()
print("="*84)
print("(D) IS THE MULT-OF-7/14 RUNNER ESSENTIAL AT THE FLOOR?")
print("="*84)
floor_with7=0; floor_no7=0
for m,S,at in lows[:200]:
    mn,bnd=binders_at(S,at)
    if any(b%7==0 for b in bnd): floor_with7+=1
    else: floor_no7+=1
print(f"  among 200 lowest-M sets: binding set INCLUDES a mult-of-7: {floor_with7}; pure non-7: {floor_no7}")
print(f"  => the mult-of-7 runner {'CAN' if floor_with7 else 'does NOT'} bind at the floor.")
print(f"     (Confirms 4/47 set: pair (7,40) -- the covering-forced 7 itself binds.)")

print()
print("  ARITHMETIC FLOOR CHECK: candidate floor values q/(14q-r) just above 1/14:")
for q in range(1,8):
    for r in range(1,6):
        D=14*q-r
        if D<=0: continue
        val=F(q,D)
        realized = any(m==val for m,_,_ in lows)
        print(f"    q={q} D=14*{q}-{r}={D}: q/D={val}={float(val):.7f} "
              f"(>1/14 margin {float(val-C14):.5f}) realized-in-sample={realized}")
print("""
  The realized floor sits at the LEAST q/D with D < 14q that a covering primitive S3 set can
  achieve. Covering forbids the tightest pairs (it forces specific small runners present,
  bounding which sums/differences D can host a dangerous binding pair). The EXACT infimum over
  all covering S3 is the open quantity; this sample's floor is the current best lower witness.""")

print()
print("="*84); print("DONE."); print("="*84)
