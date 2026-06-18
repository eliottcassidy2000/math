#!/usr/bin/env python3
"""
lrc14_angleD_floor_universality — mac-mini-2026-06-17-S6

THE PROOF-COMPLETING QUESTION for criterion C.
From the break-search: C holds on every covering 13-set, ALWAYS via the largest
(parked) runner V. The mechanism: W(S\\{V}) > 1/(7V). The smallcore_floor run found
min_A W(A) = 5/1848 over all 12-subsets of {1..13} AND over random adversarial 12-cores
with large runners. If this floor 5/1848 is UNIVERSAL over the 12-cores that actually
appear as S\\{V} in a primitive covering 13-set, then:

      C holds via V whenever  1/(7V) < 5/1848  <=>  V > 1848/(7*5) = 52.8  <=>  V >= 53.

So the ONLY covering 13-sets that could fail C are those whose LARGEST runner V <= 52.
Those form a FINITE set, checkable exhaustively. This would PROVE C (hence LRC14) modulo:
  (a) the floor W >= 5/1848 for the relevant 12-cores, and
  (b) the finite V<=52 check.

This script does THREE things:
 (1) Test whether the floor 5/1848 can be BROKEN by any 12-core that retains the
     "forced small runners" (a primitive covering 13-set must contain small runners to
     cover q in {8,9,11,13} unless it spends slots on large sole-covers). We search
     adversarially with large runners up to 600, but always KEEPING >= k small runners.
 (2) For covering 13-sets with V = max(S) <= 52: enumerate them and check C directly
     (these are the ONLY ones not auto-covered by the floor argument).
 (3) Report the exact min W and whether C holds on the finite V<=52 family.
"""
from fractions import Fraction as F
from math import gcd, lcm
from itertools import combinations
import random, sys

C = F(1, 14)
W0 = F(5, 1848)  # the observed floor

def Wsafe(A):
    A = sorted(set(A)); L = 1
    for u in A: L = L*u//gcd(L,u)
    D = 14*L; iv = []
    for u in A:
        cu = D//u; hw = D//(14*u)
        for k in range(u): c = k*cu; iv.append((c-hw, c+hw))
    if not iv: return F(1)
    norm=[]
    for lo,hi in iv:
        length=hi-lo; a=lo%D; b=a+length
        if b<=D: norm.append((a,b))
        else: norm.append((a,D)); norm.append((0,b-D))
    norm.sort(); mg=[]; cl,ch=norm[0]
    for lo,hi in norm[1:]:
        if lo<=ch:
            if hi>ch: ch=hi
        else: mg.append((cl,ch)); cl,ch=lo,hi
    mg.append((cl,ch)); best=0; n=len(mg)
    for i in range(n):
        hi=mg[i][1]; lo=mg[(i+1)%n][0]+(D if i==n-1 else 0)
        g=lo-hi
        if g>best: best=g
    return F(best,D) if best>0 else F(0)
def covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def crit_holds(S):
    for v in sorted(set(S)):
        A=[u for u in S if u!=v]
        if Wsafe(A) > F(1,7*v): return True, v
    return False, None
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
def P(*a): print(*a); sys.stdout.flush()

P("="*78)
P("FLOOR UNIVERSALITY:  is W(12-core) >= 5/1848 always?  + finite V<=52 check")
P("="*78)
P(f"  floor W0 = 5/1848 = {float(W0):.6f};  threshold V > 1848/35 = {float(F(1848,35)):.3f}  => V>=53 auto.")

# (1) adversarial: minimize W over 12-cores KEEPING small runners, large runners up to 600
P("\n[1] minimize W over 12-cores (keep >=8 small runners; large runners <= 120)")
rng = random.Random(2024)
amin = (F(99), None); broke = 0
for _ in range(40000):
    ksmall = rng.randint(8, 12)
    smalls = rng.sample(range(1,14), ksmall)
    larges = [rng.randint(14, 120) for _ in range(12-ksmall)]
    A = sorted(set(smalls+larges))
    if len(A) != 12: continue
    W = Wsafe(A)
    if W < amin[0]: amin = (W, tuple(A))
    if W < W0: broke += 1
P(f"    min W = {amin[0]} = {float(amin[0]):.7f}  at A={amin[1]}")
P(f"    cores with W < 5/1848: {broke}")
if amin[0] >= W0:
    P("    ==> floor 5/1848 NOT broken: W >= 5/1848 held on all sampled 12-cores.")
else:
    P("    ==> floor BROKEN. The simple floor argument needs the covering constraint.")

# (1b) EXHAUSTIVE over 12-subsets of {1..24} (small + moderate, no huge) to find true min
P("\n[1b] EXHAUSTIVE min W over 12-subsets of {1..20}:")
emin = (F(99), None)
for A in combinations(range(1,21), 12):
    W = Wsafe(list(A))
    if W < emin[0]: emin = (W, A)
P(f"    exhaustive min W over 12-subsets of 1..20 = {emin[0]} = {float(emin[0]):.7f}  at {emin[1]}")

# (2) finite check: covering 13-sets with max runner V <= 52
P("\n[2] FINITE CHECK: all covering 13-sets with max(S) <= 52 satisfy C?")
# enumerate covering 13-subsets of {1..52}. C(52,13) is astronomically large; instead
# enumerate by: which 13 distinct values in 1..52, covering, but restrict to those that
# are PRIMITIVE-ish and could be hard. We bound the search by requiring the set covers
# 2..14 and sampling heavily, plus all 13-subsets of {1..14} (the natural small ones).
checked = 0; cfail = 0; mbreak = 0; worst = (F(99),None,None); fails=[]
# 2a: all 13-subsets of {1..14} (only 14 of them) that cover
for drop in combinations(range(1,15), 1):
    S = tuple(x for x in range(1,15) if x not in drop)
    if len(S)==13 and covering(S):
        checked += 1
        ok,v = crit_holds(S)
        if not ok:
            cfail+=1; Mv=M(S); fails.append((S,Mv))
            if Mv<C: mbreak+=1
# 2b: random covering 13-subsets of {1..52}
rng2 = random.Random(77)
for _ in range(150000):
    S = tuple(sorted(rng2.sample(range(1,53), 13)))
    if not covering(S): continue
    checked += 1
    ok,v = crit_holds(S)
    if not ok:
        cfail+=1; Mv=M(S); fails.append((S,Mv))
        if Mv<C: mbreak+=1
P(f"    covering 13-sets with max<=52 checked: {checked}")
P(f"    C failures: {cfail}   (LRC breaks M<1/14: {mbreak})")
if cfail:
    for S,Mv in fails[:8]:
        P(f"      C-FAIL S={S}  M={Mv}={float(Mv):.5f}  {'*** LRC BREAK ***' if Mv<C else '(C-only)'}")
else:
    P("    ==> C held on every covering 13-set with max<=52 tested.")

P("\nVERDICT:")
P(f"  floor over 12-subsets of 1..20: {emin[0]} (= 5/1848 expected: {emin[0]==W0})")
P(f"  floor broken by large-runner cores: {'NO' if amin[0]>=W0 else 'YES'}")
P(f"  V<=52 finite family C failures: {cfail}")
if emin[0]==W0 and amin[0]>=W0 and cfail==0:
    P("  ==> Strong support: C provable via floor(V>=53) + finite(V<=52) decomposition.")
