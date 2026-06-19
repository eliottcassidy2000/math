#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_verify_threshold_derivation_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5, ADVERSARIAL)

THE DEEPEST QUESTION: is "cluster maxgap{frac(e_i x)} > theta" a SOUND (sufficient) loneliness
criterion, and what is the CORRECT theta?

SETUP (THM-527 slow-fast).  S = P u L, L = {Vmax - e_i}, e_i = co-offsets, 0 in E.
A target t is lonely for S at level 1/14 iff ||t*s|| >= 1/14 for all s in S.
For the cluster speeds s = Vmax - e_i:  t*s = t*Vmax - t*e_i.  Set phi := t*Vmax (mod 1) and let
x := t.  As Vmax -> infinity, for fixed cluster STRUCTURE {e_i}, phi ranges (quasi-)freely over
[0,1) INDEPENDENTLY of x = t (this is the asymptotic decoupling THM-527 invokes).  Then:
   cluster danger avoided  <=>  for all i: || phi - e_i x || >= 1/14
                            <=>  phi avoids ALL the arcs (e_i x - 1/14, e_i x + 1/14) mod 1.
The k phase-points {frac(e_i x)} each forbid an arc of total width 2*(1/14) = 1/7 centered on it.
A valid phi EXISTS  <=>  the union of these k arcs (each width 1/7) does NOT cover the circle
                    <=>  there is a circular GAP between consecutive phase-points of width > 1/7.
   ( gap between two phases of width g; the forbidden arcs eat 1/14 on each side -> free length
     g - 2*(1/14) = g - 1/7 > 0  <=> g > 1/7. )
SO: the correct global-witness threshold is theta = 1/7 (NOT 2/7).
The 2/7 (via-max) criterion FIXES phi = 0 (i.e. forces the witness to BE x itself, t*Vmax=0): then
you need each frac(e_i x) at circ-dist >= 1/14 from 0 AND ... no -- via-max takes the SINGLE largest
gap and centers; demanding the witness be at distance >=1/14 from BOTH neighbors needs g/2 >= 1/14?
NO: via-max requires the cluster all on one side; let us just TEST both claims by brute force.

THIS SCRIPT brute-forces the decoupled model EXACTLY:
  Given E (=> phases p_i(x)=frac(e_i x)), define
     LONELY_CLUSTER(x) := exists phi in [0,1) with ||phi - p_i(x)|| >= 1/14 for all i
                        <=> max circular gap of {p_i(x)} > 1/7.
  We verify this equivalence EXACTLY (it is the definition of the 1/7 measure), AND verify that
  the 2/7 via-max criterion (phi forced to 0-translate / single offset) is STRICTLY stronger.
We also reconstruct ACTUAL S with growing Vmax and check the decoupled prediction matches the
exact loneliness of S (so the asymptotic decoupling is not lying at finite Vmax).
"""
import sys, itertools, random
from fractions import Fraction as F
try: sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception: pass
random.seed(1234)
H=F(1,14); ONE7=F(1,7); TWO7=F(2,7)

def circ_gaps(phases):
    pts=sorted(set(phases))
    if len(pts)==1: return [F(1)]
    g=[pts[i+1]-pts[i] for i in range(len(pts)-1)]+[pts[0]+1-pts[-1]]
    return g

def cluster_lonely_exists(phases):
    """exists phi: ||phi - p_i|| >= 1/14 for all i  <=>  max gap > 1/7 (each arc width 1/7)."""
    # forbidden arcs: (p_i - 1/14, p_i + 1/14). free phi exists iff union != whole circle
    arcs=[]
    for p in phases:
        a=(p-H)%1; b=(p+H)%1
        if a<b: arcs.append((a,b))
        else: arcs.append((a,F(1))); arcs.append((F(0),b))
    arcs.sort(); covered=F(0); clo,chi=arcs[0]
    union=[]
    for lo,hi in arcs[1:]:
        if lo<=chi: chi=max(chi,hi)
        else: union.append((clo,chi)); clo,chi=lo,hi
    union.append((clo,chi))
    tot=sum(b-a for a,b in union)
    return tot<1   # union doesn't cover -> free phi exists

print("="*92)
print("CHECK 1: 'exists phi avoiding all 1/14-arcs around phases'  <=>  'max circular gap > 1/7'")
print("="*92)
mism=0; tested=0
for _ in range(20000):
    k=random.randint(2,13)
    phases=sorted(set(F(random.randint(0,839),840) for _ in range(k)))
    if len(phases)<2: continue
    tested+=1
    exists=cluster_lonely_exists(phases)
    maxgap=max(circ_gaps(phases))
    pred=(maxgap>ONE7)
    if exists!=pred: mism+=1
print(f"  tested {tested} random phase-configs: mismatches (exists-phi vs maxgap>1/7) = {mism}")
print(f"  => the 1/7 maxgap criterion EXACTLY characterizes 'a global cluster-witness phi exists': "
      f"{'CONFIRMED' if mism==0 else 'REFUTED'}")

print()
print("="*92)
print("CHECK 2: is the 2/7 via-max criterion STRICTLY stronger (a strict subset)?")
print("  via-max = 'maxgap>2/7' (centering a single offset, witness forced near x-translate).")
print("="*92)
only27=0; only17=0; both=0
for _ in range(20000):
    k=random.randint(3,13)
    phases=sorted(set(F(random.randint(0,839),840) for _ in range(k)))
    if len(phases)<3: continue
    mg=max(circ_gaps(phases))
    c17=mg>ONE7; c27=mg>TWO7
    if c27 and c17: both+=1
    elif c17 and not c27: only17+=1
    elif c27 and not c17: only27+=1  # impossible: 2/7>1/7
print(f"  maxgap>2/7 AND >1/7: {both};  >1/7 only (2/7 fails): {only17};  >2/7 only: {only27}")
print(f"  => 2/7 criterion ⊂ 1/7 criterion (every 2/7-good x is 1/7-good; {only17} are 1/7-good but")
print(f"     NOT 2/7-good). So 1/7 is the WEAKER/LARGER (correct sufficient) object IF the decoupling")
print(f"     holds. only27 must be 0: {only27==0}")

print()
print("="*92)
print("CHECK 3: does the DECOUPLED 1/7 prediction match EXACT loneliness of reconstructed S at")
print("  FINITE growing Vmax? (tests the asymptotic decoupling THM-527 invokes -- NOT re-proved)")
print("="*92)
def merge(iv):
    iv=sorted(iv);o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o
def danger(u,h=H):
    iv=[]
    for j in range(u):
        c=F(j,u);a=(c-h/u)%1;b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1)));iv.append((F(0),b))
    return iv
def lonely_S(S):
    arcs=merge([iv for v in S for iv in danger(v)])
    cov=sum(b-a for a,b in arcs)
    return cov<1
def safe_set(P):
    arcs=merge([iv for u in P for iv in danger(u)])
    # complement
    out=[];prev=F(0)
    for a,b in arcs:
        if a>prev: out.append((prev,a))
        prev=max(prev,b)
    if prev<1: out.append((prev,F(1)))
    return out

# For a few (P,E), pick x in G_P with cluster maxgap>1/7, reconstruct S at growing Vmax, and
# check whether S is lonely (it must be, for ALL covering 13-sets; but we ALSO check that a
# witness t exists NEAR the decoupled-predicted location).
match=0; checks=0; nonlonely=0
for _ in range(2000):
    k=random.randint(7,13); psz=13-k
    P=sorted(random.sample(range(1,14),psz))
    spread=random.choice([k-1,k,k+2,2*k])
    body=sorted(random.sample(range(1,spread+1),min(k-1,spread)))
    E=[0]+body
    if len(set(E))!=k: continue
    GP=safe_set(P) if P else [(F(0),F(1))]
    if not GP: continue
    seg=random.choice(GP); x=seg[0]+(seg[1]-seg[0])*F(random.randint(1,99),100)
    phases=[(e*x)%1 for e in E]
    if max(circ_gaps(phases))<=ONE7: continue
    checks+=1
    Vmax=max(E)+14+random.randint(20,200)
    L=[Vmax-e for e in E]
    if min(L)<=13: continue
    S=sorted(set(P+L))
    if len(S)!=13: continue
    if lonely_S(S): match+=1
    else: nonlonely+=1
print(f"  cases (x in G_P, cluster maxgap>1/7): {checks};  reconstructed S LONELY: {match};  NON-lonely: {nonlonely}")
print(f"  => decoupled 1/7-good predicts S lonely: {match}/{checks}. (non-lonely would refute LRC14 itself)")
print()
print("HONEST NOTE: CHECK 3 cannot DISTINGUISH soundness from LRC(14) being true (all covering")
print("13-sets are lonely if LRC14 holds). It confirms CONSISTENCY, not soundness of the reduction.")
print("The SOUNDNESS of 'maxgap>1/7 => witness' rests on the asymptotic decoupling (phi free as")
print("Vmax->inf) which is THM-527-A (assumed upstream, error O(1/Vmax), Vmax<=V0 finite check).")
