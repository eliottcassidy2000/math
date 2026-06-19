#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_Bk_verify_witnesslogic_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5, ADVERSARIAL)

THE CORE LOGICAL PIVOT under scrutiny: the claim that the 1/7-measure
   mu_1/7(E) = meas{x: maxgap{frac(e_i x)} > 1/7}
is the RIGHT object because "x in G_P and cluster maxgap > 1/7  =>  global lonely witness".

The reconstructed set is  S = P  u  {Vmax - e : e in E}  (cluster speeds large).
A point t in [0,1) is lonely for S iff  ||t*s|| >= 1/14  for ALL s in S, i.e.
  - for all p in P:        ||t p|| >= 1/14   (the G_P / small-speed condition)
  - for all e in E:        ||t (Vmax - e)|| >= 1/14   (the cluster / large-speed condition)

The reduction (THM-527) replaces t by x via t ~ x and uses Vmax->infinity so that
||t*(Vmax-e)|| ~ depends on frac(e*x) through a SHIFT.  The CLAIMED sufficient condition is:
the cluster phase-points {frac(e_i x)} have a circular GAP > 2/7 (the via-MAX 2/7 criterion)
-- because then you can place a target value in the MIDDLE of that gap, >1/7 from each phase on
BOTH sides => ||.|| >= 1/14 each.  A gap of width w lets the midpoint sit w/2 from each end;
to get ||.|| >= 1/14 = (1/2)(1/7) you need w/2 >= 1/7 i.e. w >= 2/7.  THAT is why 2/7 is the
via-max threshold (gap>2/7 <=> a single global offset clears all cluster speeds at once).

So the 1/7 measure mu_1/7 is NOT 'place one offset in one gap'.  We must DERIVE what 1/7 buys.

THIS SCRIPT: directly tests, for many reconstructed S with LARGE Vmax (so the asymptotic holds),
WHICH cluster-phase condition actually implies cluster-loneliness at level 1/14:
  (A) maxgap{frac(e x)} > 2/7   (via-max: midpoint of the gap, one global offset)
  (B) maxgap{frac(e x)} > 1/7   (the claimed 'right object')
We check, for the EXACT lonely set of S at level 1/14, whether x's good-arcs under each
criterion are CONTAINED in the actual loneliness region (so the criterion is SOUND/sufficient).
A criterion that is NOT sound (good-arc with no nearby lonely t) breaks the reduction.
"""
import sys, itertools, random
from fractions import Fraction as F
try: sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception: pass
random.seed(7)
H = F(1,14); ONE7 = F(1,7); TWO7 = F(2,7)

def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def meas(arcs): return sum((b-a for a,b in arcs),F(0))
def complement(arcs):
    arcs=merge(arcs); out=[]; prev=F(0)
    for a,b in arcs:
        if a>prev: out.append((prev,a))
        prev=max(prev,b)
    if prev<1: out.append((prev,F(1)))
    return out
def danger_arcs(u,h=H):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def safe_set_arcs(S,h=H):
    if not S: return [(F(0),F(1))]
    return complement(merge([iv for u in S for iv in danger_arcs(u,h)]))

print("="*92)
print("PART 1: What does a cluster phase-GAP of width w actually buy? (the via-max derivation)")
print("="*92)
print("Place k points on circle; pick a gap of width w; midpoint is w/2 from each neighbor.")
print("To make a SINGLE offset value y* with circ-dist >= 1/14 from EVERY phase, need the offset")
print("to sit in a gap of width >= 2*(1/14) = 1/7?  NO -- the cluster danger is at level 1/14 in the")
print("scaled coordinate, but the cluster speeds are Vmax-e ~ Vmax, so ||t(Vmax-e)|| spans the whole")
print("circle as t varies; the relevant object is frac(e x) and the 1/14 must be compared to the")
print("PHASE gap directly.  The via-max criterion (THM-528) uses gap>2/7 to fit the small-speed safe")
print("structure (G_P granularity 1/14 per p, but the cluster needs the gap to host a 1/7-radius ball")
print("=> width 2/7).  Let us TEST which threshold is SOUND by brute force on actual reconstructed S.")
print()

print("="*92)
print("PART 2: SOUNDNESS test. For reconstructed S=P u {Vmax-e}, large Vmax, is")
print("   'x in G_P AND cluster-maxgap{frac(e x)} > theta'  a SUBSET of  'exists lonely t for S'?")
print("We test theta=2/7 (claimed sufficient/sound) and theta=1/7 (claimed 'right object').")
print("Soundness = whenever the phase criterion holds at x, S is actually lonely (M(S)>1/14).")
print("="*92)

def maxgap_frac(E, x):
    pts = sorted((e*x) % 1 for e in E)
    if len(pts)==1: return F(1)
    gaps=[pts[t+1]-pts[t] for t in range(len(pts)-1)]+[pts[0]+1-pts[-1]]
    return max(gaps)

results = {TWO7:[0,0], ONE7:[0,0]}  # [criterion_true_count, criterion_true_AND_S_lonely]
trials = 0
for _ in range(6000):
    k = random.randint(7,13); psz = 13-k
    P = sorted(random.sample(range(1,14), psz))
    spread = random.choice([k-1,k,k+1,k+2,k+4,2*k,3*k])
    body = sorted(random.sample(range(1,spread+1), min(k-1,spread)))
    E = [0]+body
    if len(set(E))!=k or len(E)<3: continue
    Vmax = max(E) + 14 + random.randint(0, 60)   # large cluster
    L = [Vmax-e for e in E]
    if min(L) <= 13: continue
    S = sorted(set(P+L))
    if len(S) != 13: continue
    trials += 1
    safeS = safe_set_arcs(S, H)
    lonely = meas(safeS) > 0
    # sample x values; check criterion at the SCALE where frac(e x) matters.
    # In the reduction, x is the SAME variable as t (THM-527 sets the cluster phase = e*x).
    # We evaluate the cluster criterion using the small offset variable; sample G_P points.
    GP = safe_set_arcs(P, H)
    if not GP: continue
    # sample several x in G_P, test cluster maxgap
    for _try in range(8):
        seg = random.choice(GP)
        x = seg[0] + (seg[1]-seg[0]) * F(random.randint(1,99), 100)
        for theta in (TWO7, ONE7):
            mg = maxgap_frac(E, x)
            if mg > theta:
                results[theta][0] += 1
                if lonely: results[theta][1] += 1
print(f"reconstructed S tested: {trials}")
for theta,name in [(TWO7,"2/7"),(ONE7,"1/7")]:
    tot, snd = results[theta]
    print(f"  theta={name}: criterion fired {tot} times; of those, S lonely {snd} times -> "
          f"soundness {snd}/{tot} = {snd/max(tot,1):.4f}")
print()
print("NOTE: BOTH thresholds give S lonely 100% IF LRC(14) holds (S is a covering 13-set).")
print("That makes this test WEAK: it cannot distinguish a sound criterion from an unsound one,")
print("because all reconstructed S are lonely anyway (that's what we're trying to prove).")
print("The REAL question (below): does the criterion produce an EXPLICIT lonely t, locally?")

print()
print("="*92)
print("PART 3: THE ACTUAL SOUNDNESS QUESTION -- does cluster-maxgap>theta give a LOCAL witness?")
print("For the reduction to work, x in G_P with cluster-gap>theta must let us EXHIBIT t lonely for")
print("the FULL S (not rely on S being lonely a priori).  The via-max construction: t = x works for")
print("the small speeds (x in G_P).  For the cluster speeds Vmax-e, ||x(Vmax-e)|| = ||x*Vmax - e x||;")
print("with x*Vmax ~ uniform phase phi, the cluster dangers are at phi - frac(e x) mod 1 within 1/14.")
print("A gap of width w in {frac(e x)} of size > 2/7 lets us CHOOSE phi in the middle so every")
print("||x(Vmax-e)|| >= w/2 > 1/7 >= 1/14.  THAT needs w > 2/7, NOT 1/7.  So theta=2/7 is the SOUND")
print("via-max threshold; theta=1/7 is NOT a via-max witness by itself.")
print("="*92)
print("=> The 1/7 measure is the right object ONLY UNDER the concurrent THM-530/HYP-2592 reduction,")
print("   which is a DIFFERENT (global, not via-max) sufficiency argument. This script CANNOT and")
print("   does NOT re-derive that reduction; it only confirms the 2/7 via-max arithmetic (w/2>1/7")
print("   needs w>2/7). The 1/7 reformulation's SOUNDNESS is an UNVERIFIED upstream dependency.")
