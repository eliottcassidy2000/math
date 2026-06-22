#!/usr/bin/env python3
"""Creative swings at the LRC(14) nodes, post-main-sync. kind-mendel-2026-06-22-S4.
Swing 1: verify the zeta(2)=3/pi^2 floor (HYP-2856) and probe whether it survives G_P.
Swing 3: CRT 'constructive witness' sidestep on covering sets."""
from fractions import Fraction as F
from math import gcd, floor, pi
from functools import reduce
from itertools import combinations
def gl(xs): return reduce(lambda a,b:a*b//gcd(a,b),[x for x in xs if x],1)
def gall(xs): return reduce(gcd,[x for x in xs if x],0)
def nrm(y): f=y-floor(y); return min(f,1-f)

# ---------- exact GOOD = {maxgap of {frac(e_i x)} > thr} measure ----------
def good_meas(E, thr, P=None):
    "meas{x: maxgap{frac(e x):e in E}>thr  AND (P=None or x in G_P at 1/14)} via fine exact breakpoints"
    # breakpoints where a phase crosses another or a gap = thr boundary are irrational;
    # use a fine rational grid for GOOD; exact for G_P membership not needed (grid).
    pos=[e for e in E]
    N=20000
    cnt=0; cntP=0
    for i in range(N):
        x=F(2*i+1,2*N)
        ph=sorted(float((e*x)%1) for e in E)
        if len(ph)==0: continue
        mg=max([ph[0]+1-ph[-1]]+[ph[j+1]-ph[j] for j in range(len(ph)-1)])
        good = mg>float(thr)
        if good:
            cnt+=1
            if P is not None and all(nrm(p*x)>=F(1,14) for p in P): cntP+=1
    return cnt/N, (cntP/N if P is not None else None)

print("=== SWING 1: the zeta(2)=3/pi^2 floor (HYP-2856) ===")
print(f"3/pi^2 = 1/(2 zeta(2)) = {3/pi**2:.5f}")
for q in [3,5,7,9]:
    E=list(range(0,2*q-1))           # {0,1,...,2q-2}
    c,_=good_meas(E, F(1,q))
    print(f"q={q}: cluster {{0..{2*q-2}}}, c_q = meas(maxgap>1/{q}) ~ {c:.4f}")
print(f"  (HYP-2855 claims c_q ~ 0.55 actual; rate-V lower bound -> 3/pi^2={3/pi**2:.4f})")

print("\n=== SWING 1b: does the floor SURVIVE the G_P intersection? (the flagged gap) ===")
print("    meas(GOOD_E ∩ G_P) for sub-cluster E + small part P (the ACTUAL witness floor)")
m_P=14249/252252
tests=[
 (list(range(8)),[1,5,7,8,9]),       # consec_8 + cap-achieving P (k=8 binding)
 (list(range(8)),[2,3,4,5,6]),       # the HYP-2851 'worst P' that killed the conservative delta bound
 ([0,2,4,6,8,10,12],[1,3,5]),        # even cluster + odd P (resonance risk)
 (list(range(6)),[1,2,3,5,7,8,9]),
]
for E,P in tests:
    g_un,g_gp = good_meas(E, F(1,7), P)
    print(f"E={E if len(E)<8 else 'consec_'+str(len(E))} P={P}: meas(GOOD)~{g_un:.4f}  meas(GOOD∩G_P)~{g_gp:.4f}  (m_P={m_P:.4f}, ratio {g_gp/m_P:.1f}x)")

print("\n=== SWING 3: CRT 'constructive witness' tau=a/D for covering sets (sidestep analysis) ===")
print("    For covering S, search small D and a with ||s a/D||>=1/14 for all s in S (=> M>=1/14 explicitly)")
def is_covering(S):
    return all(any(s%q==0 for s in S) for q in range(2,15))
def crt_witness(S, Dmax=60):
    "find (a,D), D<=Dmax, gcd(a,D)=1, with ||s a/D||>=1/14 for all s; return (D,a) or None"
    best=None
    for D in range(2,Dmax+1):
        for a in range(1,D):
            if gcd(a,D)!=1: continue
            if all(nrm(F(s*a,D))>=F(1,14) for s in S):
                return (D,a,min(float(nrm(F(s*a,D))) for s in S))
    return None
import random
random.seed(1)
found=0; tot=0; fails=[]
# random covering 13-sets
trials=0
while tot<40 and trials<4000:
    trials+=1
    S=sorted(random.sample(range(1,200),13))
    if gall(S)!=1 or not is_covering(S): continue
    tot+=1
    w=crt_witness(S, Dmax=70)
    if w: found+=1
    else: fails.append(S)
print(f"covering 13-sets tested: {tot}; CRT witness (D<=70) found: {found}; failed: {len(fails)}")
if fails: print("  first failures (need larger D or genuinely tight):", fails[:3])
# the AP tight cases
for S in [list(range(1,14)), [1,2,3,4,5,6,7,8,9,10,11,12,26]]:
    w=crt_witness(S,Dmax=200)
    print(f"  S={S[:6]}...: covering={is_covering(S)} witness={w}")
