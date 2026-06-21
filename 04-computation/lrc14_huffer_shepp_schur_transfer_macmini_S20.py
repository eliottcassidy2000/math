#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ANGLE 7 (web trawl -> falsifiable test): does the HUFFER-SHEPP Schur-convexity
of circle-coverage transfer to the LRC(14) LAYER-3 wall (consec Schur-maximizes
sum_a W_a)?

KNOWN RESULT (Huffer & Shepp 1987, J. Appl. Prob. 24:422-429; Proschan's conj):
  Throw n arcs of lengths l_1..l_n uniformly on the unit circle.  P(l_1..l_n) =
  P(circle covered) is SCHUR-CONVEX in (l_1..l_n): more SPREAD-OUT (unequal) arc
  lengths => MORE coverage.  Equal lengths = the Schur-MINIMUM of coverage.
  Proof: pairwise transfer l_1->l_1+e, l_2->l_2-e (Robin-Hood reversed) only
  increases coverage; reflection lemma "an interval of fixed length is more
  likely covered the closer it is to the center".

THE LRC(14) LAYER-3 WALL (HYP-2753/2756/2757):
  measS7(E) = sum_{a=1..6} W_a(E), each W_a a cycle-graph C_7 survival width.
  consec={0,7,14,...} (the AP) is the UNIQUE maximizer over the full-residue
  stratum.  Z/7* permutes the phases a->ca; consec is the dilation-symmetric
  (most-equal-spacing / balanced) config.

THE SIGN TENSION I MUST RESOLVE:
  Huffer-Shepp: coverage maximized at the SPREAD-OUT extreme (Schur-CONVEX).
  LRC wall:     coverage(=measS7) maximized at the BALANCED/equal point (consec).
  If the LRC functional is Schur-CONCAVE in the relevant variable, the two are
  the SAME theorem with flipped sign, and consec wins by being balanced.

FALSIFIABLE HYPOTHESIS (HYP-HS):
  The LRC LAYER-3 functional measS7 = sum_a W_a is SCHUR-CONCAVE in the
  per-residue GAP/spacing vector: a Robin-Hood transfer that EQUALIZES the
  spacing of the residue magnitudes does NOT decrease measS7; equivalently a
  spread-increasing transfer does NOT increase it.  consec = the perfectly
  equal-gap (Schur-minimal-spread) point => the Schur-concave maximum.

TESTS:
  (A) Build the survival-arc model for a single phase W_a from definitions and
      confirm sum_a W_a = measS7 (sanity, ties to HYP-2753 exact decomposition).
  (B) DIRECT Schur test on the LRC functional: enumerate pairs E, E' on the SAME
      full-residue stratum where E' majorizes E in the magnitude/gap vector
      (E' more spread out).  Check sign of measS7(E')-measS7(E).
      HYP-HS predicts measS7 DECREASES as spread increases (Schur-concave).
  (C) Cross-check vs the ACTUAL Huffer-Shepp arc-coverage measure on Z/7 to see
      whether the LRC measure and the HS measure have OPPOSITE Schur signs (the
      key claim), by simulating the analogous deterministic coverage.
"""
import itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
from collections import defaultdict

# ---------- occupancy law (exact, reused from THREAD A) ----------
def occupancy_pi(E):
    E=sorted(set(E))
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        for a in range(7*abs(e)+1): bps.add(F(a,7*abs(e)))
    bps=sorted(b for b in bps if 0<=b<=1)
    pi=[F(0)]*8
    for lo,hi in zip(bps,bps[1:]):
        if hi<=lo: continue
        xm=(lo+hi)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        pi[len(hit)]+=hi-lo
    return pi

def measS7(E):
    return occupancy_pi(E)[7]

# ---------- (A) per-phase survival width W_a, from definitions ----------
# Near x=a/7 each clock e sits at sector (e*a mod 7) and drifts at speed 7e.
# W_a = meas{ x in [a/7-1/14, a/7+1/14] : all 7 sectors covered }.
def W_a(E, a):
    """exact survival width of phase a: measure of full-coverage in the cell."""
    E=sorted(set(E))
    lo=F(a,7)-F(1,14); hi=F(a,7)+F(1,14)
    bps={lo,hi}
    for e in E:
        if e==0: continue
        # breakpoints where clock e crosses a sector boundary: e*x = j/7
        # x = j/(7e); collect those inside [lo,hi]
        d=7*abs(e)
        jlo=(lo*d); jhi=(hi*d)
        import math
        j0=int(jlo) if jlo==int(jlo) else math.floor(jlo)
        j1=int(jhi) if jhi==int(jhi) else math.ceil(jhi)
        for j in range(j0-1,j1+2):
            x=F(j,d)
            if lo<=x<=hi: bps.add(x)
    bps=sorted(bps)
    tot=F(0)
    for l,h in zip(bps,bps[1:]):
        if h<=l: continue
        xm=(l+h)/2; hit=set()
        for e in E:
            v=e*xm; v=v-(v.numerator//v.denominator); hit.add((v.numerator*7)//v.denominator)
        if len(hit)==7: tot+=h-l
    return tot

def sum_Wa(E):
    return sum(W_a(E,a) for a in range(1,7))

# ---------- majorization helpers ----------
def majorizes(x,y):
    """x majorizes y (x more spread): sorted-desc partial sums of x >= y, equal total."""
    xs=sorted(x,reverse=True); ys=sorted(y,reverse=True)
    if len(xs)!=len(ys): return False
    if sum(xs)!=sum(ys): return False
    px=0;py=0
    for i in range(len(xs)-1):
        px+=xs[i];py+=ys[i]
        if px<py: return False
    return True

def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return list(range(k))

if __name__=="__main__":
    print("="*78)
    print("(A) sanity: sum_{a=1..6} W_a(E) == measS7(E) for several shapes")
    print("="*78)
    for E in [consec(8),[0,1,2,3,4,5,6,8],[0,2,3,4,5,6,7,8],consec(9)]:
        s=sum_Wa(E); m=measS7(E)
        print(f"  E={E}: sum_Wa={s}={float(s):.6f}  measS7={m}={float(m):.6f}  match={s==m}")
        # per-phase
        ws=[W_a(E,a) for a in range(1,7)]
        print(f"     W_a (a=1..6) = {[str(w) for w in ws]}")

    print()
    print("="*78)
    print("(B) DIRECT SCHUR SIGN TEST on measS7 over a full-residue stratum.")
    print("    For shapes E,E' both = {0..6}+{extra}, build the MAGNITUDE vector")
    print("    (the multiset of |e|), test majorization, check sign of d measS7.")
    print("    HYP-HS: measS7 SCHUR-CONCAVE => spread-increasing transfer DECREASES it,")
    print("            so consec (least-spread balanced) is the max.")
    print("="*78)
    # stratum: k=8 full-residue {0..6} + one extra w in residues, vary magnitude
    # Use shapes {0..6,w} for w=7..28 (all double residue 0 if w%7==0 etc.)
    # To compare majorization cleanly, fix the residue MULTISET and vary magnitudes.
    # consec = {0,1,2,3,4,5,6,7}: magnitudes (0,1,2,3,4,5,6,7).
    base=list(range(7))
    fam=[]
    for w in range(7,29):
        E=base+[w]
        if primitive(E): fam.append(tuple(E))
    # For each pair, compare magnitude vectors by majorization
    print(f"  family {len(fam)} shapes {{0..6,w}}, w=7..28")
    print(f"  consec magnitudes={tuple(range(8))} measS7={float(measS7(consec(8))):.6f}")
    schur_concave_ok=True; tested=0; violations=[]
    for Ea,Eb in itertools.combinations(fam,2):
        ma=tuple(abs(e) for e in Ea); mb=tuple(abs(e) for e in Eb)
        if sum(ma)!=sum(mb): continue
        if majorizes(ma,mb):  # ma more spread
            tested+=1
            d=measS7(list(Ea))-measS7(list(Eb))  # spread - balanced
            if d>F(1,10**15):  # more spread gave MORE -> violates Schur-concave
                schur_concave_ok=False; violations.append((Ea,Eb,float(d)))
        elif majorizes(mb,ma):
            tested+=1
            d=measS7(list(Eb))-measS7(list(Ea))
            if d>F(1,10**15):
                schur_concave_ok=False; violations.append((Eb,Ea,float(d)))
    print(f"  comparable (same-sum, majorization-ordered) pairs tested = {tested}")
    print(f"  HYP-HS (measS7 Schur-CONCAVE in magnitude vector): {'HOLDS' if schur_concave_ok else 'FALSE'}")
    for v in violations[:8]:
        print(f"     VIOLATION: more-spread {v[0]} beat balanced {v[1]} by {v[2]:.2e}")

    print()
    print("="*78)
    print("(C) Broader Schur test: ALL full-residue k=8 shapes span<=14, group by")
    print("    magnitude-sum, test EVERY majorization-comparable pair's measS7 sign.")
    print("="*78)
    W=14
    bank=[(0,)+r for r in itertools.combinations(range(1,W+1),7)]
    bank=[E for E in bank if primitive(E) and set(e%7 for e in E)=={0,1,2,3,4,5,6}]
    print(f"  full-residue shapes span<= {W}: {len(bank)}")
    by_sum=defaultdict(list)
    for E in bank:
        by_sum[sum(abs(e) for e in E)].append(E)
    concave_viol=0; convex_viol=0; comp=0
    cache={}
    def m_(E):
        if E not in cache: cache[E]=measS7(list(E))
        return cache[E]
    for s,grp in by_sum.items():
        for Ea,Eb in itertools.combinations(grp,2):
            ma=sorted((abs(e) for e in Ea),reverse=True)
            mb=sorted((abs(e) for e in Eb),reverse=True)
            if majorizes(ma,mb):
                comp+=1; d=m_(Ea)-m_(Eb)
                if d>F(1,10**15): concave_viol+=1   # spread beat balanced
                if d<-F(1,10**15): convex_viol+=1   # balanced beat spread
            elif majorizes(mb,ma):
                comp+=1; d=m_(Eb)-m_(Ea)
                if d>F(1,10**15): concave_viol+=1
                if d<-F(1,10**15): convex_viol+=1
    print(f"  majorization-comparable pairs = {comp}")
    print(f"  Schur-CONCAVE violations (more spread strictly beat balanced) = {concave_viol}")
    print(f"  Schur-CONVEX  violations (balanced strictly beat more spread)  = {convex_viol}")
    if concave_viol==0 and convex_viol>0:
        print("  => measS7 is SCHUR-CONCAVE on this stratum: CONFIRMS HYP-HS sign-flip.")
    elif convex_viol==0 and concave_viol>0:
        print("  => measS7 is SCHUR-CONVEX (same sign as Huffer-Shepp): REFUTES HYP-HS.")
    else:
        print("  => NEITHER: measS7 not Schur-monotone in raw magnitude vector.")
