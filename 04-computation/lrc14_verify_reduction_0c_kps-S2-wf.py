#!/usr/bin/env python3
"""
PARTS 5-7 (consolidated, fast): adversarial verification of the
'thin comb cannot cover fat G_C' engine of the LRC(14) easy-dominates-hard reduction.

Uses clean helpers (no slow import side effects).
"""
from fractions import Fraction as F
import importlib.util, random, itertools
from math import gcd
from functools import reduce

spec = importlib.util.spec_from_file_location(
    "h", "04-computation/lrc14_lonely_helpers_kps-S2-wf.py")
h = importlib.util.module_from_spec(spec); spec.loader.exec_module(h)
danger=h.danger; complement=h.complement; intersect=h.intersect; measure=h.measure
normalize=h.normalize; lonely_set=h.lonely_set; lonely_measure=h.lonely_measure; M=h.M

C_drop6 = [1,2,3,4,5,7,8,9,10,11,12,13]
TARGET = F(7,858)

def L_exact(C, w, hh=F(1,14)):
    G=lonely_set(C,hh); Dw=danger(w,hh)
    return measure(intersect(G, complement(Dw)))
def components(C, hh=F(1,14)):
    return len(lonely_set(C,hh))
def comb_cover_of_arc(a, b, w, hh=F(1,14)):
    return measure(intersect(normalize([(a,b)]), danger(w,hh)))

# ============================================================================
print("="*78); print("PART 5: decoupling bound  L >= (6/7)meas(G_C) - r/(7w)  exact check"); print("="*78)
mGc=lonely_measure(C_drop6); r=components(C_drop6)
print(f"drop-6 core meas(G_C)={mGc}={float(mGc):.8f}, r={r}")
viol=0; minL=None
for m in range(1,60):
    w=14*m; L=L_exact(C_drop6,w); bound=F(6,7)*mGc - F(r,7*w)
    if L<bound: viol+=1
    if minL is None or L<minL[0]: minL=(L,w)
    if m<=8: print(f"   w={w:4d}: L={float(L):.8f}  bound={float(bound):.8f}  L>=bound? {L>=bound}")
print(f"   decoupling-bound violations over m=1..59: {viol}")
print(f"   MIN L over m=1..59: {minL[0]}={float(minL[0]):.8f} at w={minL[1]}  (>0? {minL[0]>0})")
thr=F(r,1)/(F(6)*mGc)
print(f"   bound>0 needs w>=r/(6 meas)={thr}={float(thr):.1f}; below that, exact-checked min L={float(minL[0]):.8f}>0")
print()

# ============================================================================
print("="*78); print("PART 7a: LOAD-BEARING per-arc inequality  cover_w(A) <= meas(A)/7 + 1/(7w)"); print("="*78)
random.seed(7); viol=0; nt=0; tight=None
for w in [14,28,42,84,98,140,17,23,50,99,200,313,7,9,11,13]:
    for _ in range(300):
        d1=random.randint(1,400); n1=random.randint(0,d1)
        d2=random.randint(1,400); n2=random.randint(0,d2)
        a,b=F(n1,d1),F(n2,d2)
        if a>b: a,b=b,a
        if a==b: continue
        nt+=1; cover=comb_cover_of_arc(a,b,w); rhs=(b-a)/7+F(1,7*w)
        if cover>rhs:
            viol+=1
            if viol<=5: print(f"   VIOLATION w={w} A=[{a},{b}] cover={cover}>rhs={rhs}")
        slack=rhs-cover
        if tight is None or slack<tight[0]: tight=(slack,w,a,b,cover,rhs)
print(f"   tested {nt} random arcs; VIOLATIONS: {viol}")
print(f"   tightest slack={float(tight[0]):.8f} at w={tight[1]} A=[{tight[2]},{tight[3]}] cover={float(tight[4]):.8f} rhs={float(tight[5]):.8f}")
# worst-case full-tooth arc
print("   full-tooth arcs (cover=meas(A), the extreme case):")
for w in [14,84,200]:
    a=F(3,w)-F(1,14*w); b=F(3,w)+F(1,14*w); cover=comb_cover_of_arc(a,b,w); rhs=(b-a)/7+F(1,7*w)
    print(f"     w={w}: meas(A)={float(b-a):.6f} cover={float(cover):.6f} rhs={float(rhs):.6f} holds? {cover<=rhs}")
print()

# ============================================================================
print("="*78); print("PART 6/7c: ADVERSARIAL meas(G_C) minimization over many families"); print("="*78)
overall=[(TARGET, C_drop6, "drop-6 baseline")]

print("6a. {1..k}+multiples of 14:")
best=(F(1),None); t=0
for ncore in [9,10,11]:
    core=list(range(1,ncore+1)); nl=12-ncore
    for combo in itertools.combinations([14*j for j in range(1,18)], nl):
        C=sorted(core+list(combo))
        if len(set(C))!=12: continue
        t+=1; m=lonely_measure(C)
        if m<best[0]: best=(m,C)
print(f"   tested {t}; MIN={best[0]}={float(best[0]):.8f} at {best[1]}; <7/858? {best[0]<TARGET}")
overall.append((best[0],best[1],"mult-14"))

print("6b. {1..k}+multiples of 7,12,13,84,91,168:")
best2=(F(1),None); t2=0
mults=sorted(set(b*j for b in [7,12,13,84,91,168] for j in range(1,8) if b*j<=600))
for ncore in [8,9,10]:
    core=list(range(1,ncore+1)); nl=12-ncore; cnt=0
    for combo in itertools.combinations(mults, nl):
        C=sorted(set(core+list(combo)))
        if len(C)!=12: continue
        cnt+=1
        if cnt>40000: break
        t2+=1; m=lonely_measure(C)
        if m<best2[0]: best2=(m,C)
print(f"   tested {t2}; MIN={best2[0]}={float(best2[0]):.8f} at {best2[1]}; <7/858? {best2[0]<TARGET}")
overall.append((best2[0],best2[1],"mult-mixed"))

print("7c. structured-pool random 12-cores (near-covering / spread):")
best3=(F(1),None); t3=0
pool=list(range(1,15))+[16,18,20,21,22,24,26,28,33,36,39,44,48,52,55,60,66,72,77,84,88,91,99]
random.seed(99)
for _ in range(40000):
    C=sorted(random.sample(pool,12))
    if len(set(C))!=12: continue
    t3+=1; m=lonely_measure(C)
    if m<best3[0]: best3=(m,C)
print(f"   tested {t3}; MIN={best3[0]}={float(best3[0]):.8f} at {best3[1]}; <7/858? {best3[0]<TARGET}")
overall.append((best3[0],best3[1],"structured-pool"))

print("7d. random primitive 12-sets, speeds up to 300:")
best4=(F(1),None); t4=0; nb=0
random.seed(424242)
for V in [40,80,150,300]:
    for _ in range(6000):
        C=sorted(random.sample(range(1,V+1),12))
        if reduce(gcd,C)!=1: continue
        t4+=1; m=lonely_measure(C)
        if m<best4[0]: best4=(m,C)
        if m<TARGET: nb+=1
print(f"   tested {t4}; #(<7/858)={nb}; MIN={best4[0]}={float(best4[0]):.8f} at {best4[1]}")
overall.append((best4[0],best4[1],"random-300"))

print("7e. greedy best-improvement descent from drop-6 (candidates 1..150):")
cur=list(C_drop6); cur_m=lonely_measure(cur)
for rnd in range(30):
    bm=None
    for idx in range(12):
        for newv in range(1,151):
            if newv in cur: continue
            cs=sorted(cur[:idx]+[newv]+cur[idx+1:])
            if len(set(cs))!=12: continue
            m=lonely_measure(cs)
            if m<cur_m and (bm is None or m<bm[0]): bm=(m,cs)
    if bm is None: break
    cur_m,cur=bm[0],bm[1]
print(f"   final={cur_m}={float(cur_m):.8f} at {cur}; <7/858? {cur_m<TARGET}")
overall.append((cur_m,cur,"greedy"))
print()

print("="*78); print("VERDICT DATA"); print("="*78)
gmin=min(overall,key=lambda x:x[0])
print(f"   global min meas(G_C) found: {gmin[0]}={float(gmin[0]):.8f} ({gmin[2]}) at {gmin[1]}")
print(f"   target 7/858={float(TARGET):.8f}")
print(f"   ANY 12-set with meas(G_C) < 7/858 found? {gmin[0] < TARGET}")
print(f"   per-arc comb inequality violations: {viol if 'viol' in dir() else 'NA'}")
