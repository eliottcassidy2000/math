#!/usr/bin/env python3
"""
lrc14_asm_verify_finitechecktoB_kps-S10-wf.py  (FAST single-pass engine)
Adversarial verification of the assembly piece "finite-check-to-B" (Regime-A finite
certificate) for LRC(14), k=8,9,10.

CLAIM under test (RESULT PARTIAL):
  Every primitive E with span<=B(k)=16,15,13 has meas(S7(E))<=L_y(E)<=cap_k, consec the
  UNIQUE maximizer of BOTH meas(S7) and L_y, 0 violations over 11432/6435/715 primitive
  shapes; span<=N*=7,8,10 additionally PROVED by THM-536-B2 subset domination.
  caps 2243/5880, 1979/4004, 55/91; L_y margins 683/29400, 10441/7567560 (tightest),
  35411/1031940.

Sections:
  V1  EXACT re-derivation of named constants (caps, consec meas, L_y, margins).
  V2  dual feasibility g(t)>=1[t=0] on t=0..6 (the PROVED per-E Bonferroni bound).
  V3  EXHAUSTIVE bounded-spread sweep at B(k): count, 0 over cap, consec UNIQUE argmax
      of BOTH meas(S7) and L_y, margin, and meas<=L_y (proved bound) for every E.

FAST engine: for each E, one breakpoint pass over the 6 sector-arc pullbacks gives, per
elementary interval, `free`=#unhit sectors among {1..6}.  Then
  meas(S7) = sum length*[free==0],   S_r = sum length*C(free,r).
This computes meas(S7) AND all moments in ONE pass (matches the canon J/measS7 exactly;
cross-checked in V1).  EXACT (Fraction).  No git.
"""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, comb
from functools import reduce
import sys

# --------------------------------------------------------------------------
# CANON reference tools (slow) -- used only for V1 cross-checks
# --------------------------------------------------------------------------
def J(A,E):
    E=sorted(set(E)); arcs=[(F(j,7),F(j+1,7)) for j in A]; bp=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for (a,b) in arcs:
            for end in (a,b):
                m=0
                while True:
                    xv=(end+m)/e
                    if xv>=1:break
                    if xv>=0:bp.add(xv)
                    m+=1
    bp=sorted(b for b in bp if 0<=b<1); tot=F(0)
    for lo,hi in zip(bp,bp[1:]+[F(1)]):
        if hi<=lo:continue
        mid=(lo+hi)/2
        if all(not any(a<((e*mid)%1)<b for (a,b) in arcs) for e in E): tot+=hi-lo
    return tot
def measS7_canon(E): return sum(((-1)**r)*J(set(A),E) for r in range(7) for A in combinations(range(1,7),r))

# --------------------------------------------------------------------------
# FAST single-pass engine.
# For sector j in {1..6}, sector j is "hit" at x iff some e in E has frac(e x) in [j/7,(j+1)/7).
# Equivalently sector j is UNHIT at x iff for all e, frac(e x) avoids arc j.
# Build all breakpoints = pullbacks of every sector boundary {j/7,(j+1)/7} under every e.
# On each elementary interval, for each sector j independently decide hit/unhit at midpoint;
# free = #unhit among 1..6.
# --------------------------------------------------------------------------
SECTORS=[ (F(j,7),F(j+1,7)) for j in range(1,7) ]   # arcs for sectors 1..6

def fast_profile(E):
    E=sorted(set(e for e in E))
    bp=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for j in range(1,7):
            for end in (F(j,7), F(j+1,7)):
                m=0
                while True:
                    xv=(end+m)/e
                    if xv>=1: break
                    if xv>=0: bp.add(xv)
                    m+=1
    bp=sorted(b for b in bp if 0<=b<1)
    meas_s7=F(0)
    # S_r accumulators r=0..4 ; S_0 is just total length=1
    Sacc=[F(0)]*5
    nz=[e for e in E if e!=0]
    for lo,hi in zip(bp, bp[1:]+[F(1)]):
        if hi<=lo: continue
        L=hi-lo
        mid=(lo+hi)/2
        fr=[ (e*mid)%1 for e in nz ]
        free=0
        for (a,b) in SECTORS:
            hit=False
            for v in fr:
                if a<v<b: hit=True; break
            if not hit: free+=1
        if free==0: meas_s7+=L
        for r in range(5):
            if free>=r:
                Sacc[r]+=L*comb(free,r)
    return meas_s7, Sacc   # Sacc[0]=1

DUAL = {
    8:  ([F(1),F(-1),F(1),F(-9,10),F(3,5)], 4),
    9:  ([F(1),F(-13,18),F(4,9),F(-1,6)], 3),
    10: ([F(1),F(-13,18),F(4,9),F(-1,6)], 3),
    11: ([F(1),F(-1,2),F(1,6)], 2),
    12: ([F(1),F(-1,2),F(1,6)], 2),
    13: ([F(1),F(-1,2),F(1,6)], 2),
}
def Ly_from_S(S,k):
    y,R=DUAL[k]; return sum(y[r]*S[r] for r in range(R+1))

def g_factor(t,k):
    y,R=DUAL[k]; return sum(y[r]*comb(t,r) for r in range(R+1))

CAP = {8:F(2243,5880), 9:F(1979,4004), 10:F(55,91), 11:F(66,91), 12:F(6,7)}
CONSEC_MEAS = {8:F(481,1470), 9:F(2447,5880), 10:F(8899,17640),
               11:F(3419,5880), 12:F(121103,194040)}
CLAIM_MARGIN = {8:F(683,29400), 9:F(10441,7567560), 10:F(35411,1031940)}
CLAIM_NPRIM = {8:11432, 9:6435, 10:715}

ALL_OK=True
def report(tag, ok, detail=""):
    global ALL_OK
    ALL_OK = ALL_OK and ok
    print(f"[{'PASS' if ok else 'FAIL'}] {tag}" + (f"  {detail}" if detail else ""))
    return ok

print("="*78); print("V1: EXACT constants + fast-engine cross-check vs canon J/measS7"); print("="*78)
# cross-check fast engine vs canon on a few sets
for E in [[0,1,2,3,4,5,6,7],[0,1,2,3,4,5,6,16],[0,2,3,5,7,11,13,16],[0,1,2,3,4,5,6,7,8],[0,1,2,3,4,5,6,7,8,9]]:
    m_fast,S=fast_profile(E)
    m_canon=measS7_canon(E)
    report(f"V1 fast==canon meas(S7) E={E}: {m_fast}", m_fast==m_canon, f"canon {m_canon}")

for k in (8,9,10,11,12):
    m,_=fast_profile(list(range(k)))
    report(f"V1 meas(S7(consec_{k}))={m}", m==CONSEC_MEAS[k], f"claim {CONSEC_MEAS[k]}")
for k in (8,9,10):
    m,S=fast_profile(list(range(k)))
    ly=Ly_from_S(S,k); margin=CAP[k]-ly
    report(f"V1 L_y(consec_{k})={ly} margin={margin}", margin==CLAIM_MARGIN[k], f"claim {CLAIM_MARGIN[k]}")
    report(f"V1 chain meas<=L_y<=cap consec_{k}", m<=ly<=CAP[k])

print(); print("="*78); print("V2: dual feasibility g(t)>=1[t=0] on t=0..6"); print("="*78)
for k in (8,9,10,11,12,13):
    gv=[g_factor(t,k) for t in range(7)]
    report(f"V2 g(t) k={k}: {[str(x) for x in gv]}", gv[0]==1 and all(x>=0 for x in gv))

print(); print("="*78)
print("V3: EXHAUSTIVE bounded-spread sweep at B(k)=16,15,13 (k=8,9,10)")
print("    E={0} u (k-1)-subset of {1..B}, primitive(gcd=1); check 0>cap, consec UNIQUE argmax")
print("="*78)
def prim(E): return reduce(gcd,[e for e in E if e>0])==1
for k,B in [(10,13),(9,15),(8,16)]:   # cheapest first
    cap=CAP[k]; consec=list(range(k))
    cm,cS=fast_profile(consec); cly=Ly_from_S(cS,k)
    n=0; over_meas=[]; over_ly=[]; bound_viol=[]
    max_m=cm; arg_m=consec; ties_m=[]
    max_l=cly; arg_l=consec; ties_l=[]
    for rest in combinations(range(1,B+1),k-1):
        E=[0]+list(rest)
        if not prim(E): continue
        n+=1
        m,S=fast_profile(E); ly=Ly_from_S(S,k)
        if m>ly: bound_viol.append((E,m,ly))
        if m>cap: over_meas.append((E,m))
        if ly>cap: over_ly.append((E,ly))
        if m>max_m: max_m=m; arg_m=E
        if ly>max_l: max_l=ly; arg_l=E
        if E!=consec and m>=cm: ties_m.append((E,m))
        if E!=consec and ly>=cly: ties_l.append((E,ly))
    report(f"V3 k={k} primitive count={n}", n==CLAIM_NPRIM[k], f"claim {CLAIM_NPRIM[k]}")
    report(f"V3 k={k} 0 meas(S7)>cap", len(over_meas)==0, f"{over_meas[:3]}")
    report(f"V3 k={k} 0 L_y>cap", len(over_ly)==0, f"{over_ly[:3]}")
    report(f"V3 k={k} meas<=L_y for ALL E (proved bound)", len(bound_viol)==0, f"{bound_viol[:2]}")
    report(f"V3 k={k} consec UNIQUE argmax meas(S7)",
           arg_m==consec and len(ties_m)==0, f"max={max_m} arg={arg_m} ties={ties_m[:2]}")
    report(f"V3 k={k} consec UNIQUE argmax L_y",
           arg_l==consec and len(ties_l)==0, f"max={max_l} arg={arg_l} ties={ties_l[:2]}")
    print(f"      consec meas={cm} L_y={cly} cap={cap} L_y-margin={cap-cly}")
    sys.stdout.flush()

print()
print("="*78)
print(f"OVERALL (V1-V3): {'ALL PASS' if ALL_OK else 'SOME FAIL'}")
print("="*78)
