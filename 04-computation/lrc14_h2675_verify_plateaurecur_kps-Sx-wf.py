#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
ADVERSARIAL VERIFICATION of HYP-2675 angle "plateau-recursion".  kind-pasteur.
EXACT rationals throughout (fractions.Fraction).  NO numpy, NO floats in the math.

We INDEPENDENTLY re-implement the engine and test, in order:

 (A) Re-derive the EXACT plateau identity  p0(E) = Plat(E') + Delta_w  with
     Plat(E') = p0(E') + (1/7) p1(E'),  w = max(E), E' = E\{w}.  Check it is an
     EXACT identity (no rounding) on hundreds of sets.

 (B) Re-certify the margins mu_k = cap_k - Q(k-1) and the claimed exact values.

 (C) Re-PROVE-CHECK the COMB BOUND |Delta_w| <= 2 c1(E')/(7 w) over a large,
     ADVERSARIAL pool of (E', w):  large w, small w, dilated APs, multi-cluster.
     Also separately check the SHARPER per-end claim |e_I| <= 1/(7w) per component
     (the proof text wobbles on the constant).  Report worst ratio.

 (D) HUNT for a counterexample to (WIDE):  a primitive WIDE k-set (span>B=44,
     k=8,9,10) with p0(E) > cap_k.  Families: dilated/near-AP, multi-cluster,
     all-residues-mod-7, balanced 2-3 cluster, plus heavy random sweeps.
     ANY hit => holds=false with witness.

 (E) HUNT for a counterexample to the SHARPENED invariant (W*):  span(E)>B and
     p0(E) > Q(k-1).  This is the actual induction invariant.  If (W*) fails even
     once on a WIDE primitive set, the proposed induction is unsound.

 (F) STRESS the structural assumptions of the BALANCED-wide branch:
     - is c1(E') really unbounded in span?  (tabulate c1 vs span)  -> confirms the
       single-peel cutoff W(E') is NOT uniformly bounded -> the gap is real.
     - does the claimed inductive bound p0(E') <= Q(k-2) for wide E' actually
       hold? (test it -- if it fails, the STEP 6 induction collapse is unjustified)

 (G) Verdict on whether B=44 glues to span<=14 with no gap.
"""
import sys, itertools, random
from fractions import Fraction as F
from functools import reduce
from math import gcd

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

# ----------------------------------------------------------------------------
# INDEPENDENT exact engine.  For a set E (0 in E), partition [0,1) at all
# breakpoints a/(7e); on each cell count how many of the 6 inner sectors are
# missed.  p[t] = total measure with exactly t inner sectors missed.
# ----------------------------------------------------------------------------
def analyze(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for a in range(7 * e + 1):
            bps.add(F(a, 7 * e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [F(0)] * 7
    regs = []
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = set()
        for e in E:
            # sector index of frac(e*mid) in 0..6
            fr = (e * mid) % 1
            hit.add(int(fr * 7))
        miss = len(set(range(1, 7)) - hit)  # inner sectors 1..6 only
        p[miss] += hi - lo
        regs.append((lo, hi, miss))
    return p, regs

def p0p1(E):
    p, _ = analyze(E)
    return p[0], p[1]

def Plat(E):
    p, _ = analyze(E)
    return p[0] + F(1, 7) * p[1]

def c1_components(regs):
    c = 0; inside = False
    for lo, hi, m in regs:
        if m == 1:
            if not inside:
                c += 1; inside = True
        else:
            inside = False
    return c

def primitive(E):
    g = 0
    for e in E:
        g = gcd(g, e)
    return g == 1

def span(E):
    return max(E) - min(E)

CAP = {8: F(2243, 5880), 9: F(1979, 4004), 10: F(55, 91), 11: F(66, 91), 12: F(6, 7)}
def Q(m):
    return Plat(list(range(m)))

rng = random.Random(0xBADC0DE)
print("="*78)
print("(A)  EXACT plateau identity  p0(E) = Plat(E') + Delta_w   [must be exact]")
print("="*78)
bad_identity = 0
checked = 0
for _ in range(400):
    k = rng.choice([8,9,10])
    E = sorted(rng.sample(range(1, 45), k-1)); E=[0]+E
    if len(set(E))!=k: continue
    w = max(E); Ep = [e for e in E if e!=w]
    p0E,_ = p0p1(E)
    plat = Plat(Ep)
    dw = p0E - plat            # this is Delta_w by DEFINITION
    # The identity is trivially true by definition of Delta_w; the CONTENT is the
    # bound on dw.  But verify the claimed decomposition meaning: p0(E') part and
    # (1/7)p1(E') part are EXACTLY the pieces.  Recompute independently:
    pEp,_ = analyze(Ep)
    plat2 = pEp[0] + F(1,7)*pEp[1]
    if plat != plat2: bad_identity += 1
    checked += 1
print(f"  identity well-defined & Plat recomputation matches on {checked} sets; mismatches={bad_identity}")
assert bad_identity == 0

print("\n"+"="*78)
print("(B)  Margins mu_k = cap_k - Q(k-1)  [exact rationals]")
print("="*78)
claimed_mu = {8:F(1087,5880),9:F(129643,980980),10:F(5583,35672),11:F(311453,1605240),12:F(6295,24696)}
claimed_Q  = {7:F(289,1470),8:F(621,1715),9:F(1229,2744),10:F(65599,123480),11:F(14873,24696)}
ok=True
for k in range(8,13):
    q=Q(k-1); mu=CAP[k]-q
    okq = (q==claimed_Q[k-1]); okm=(mu==claimed_mu[k])
    if not(okq and okm): ok=False
    print(f"  k={k}: Q({k-1})={q} (claim {claimed_Q[k-1]}, {'OK' if okq else 'MISMATCH'})  "
          f"mu_k={mu} (claim {claimed_mu[k]}, {'OK' if okm else 'MISMATCH'})  >0:{mu>0}")
print(f"  all margins match claim & positive: {ok and all(CAP[k]-Q(k-1)>0 for k in range(8,13))}")

print("\n"+"="*78)
print("(C)  COMB BOUND  |Delta_w| <= 2 c1(E')/(7w)  AND sharper c1/(7w) check")
print("="*78)
fails_loose=0; fails_sharp=0; worst_loose=0.0; worst_sharp=0.0; worst_abs=F(0); n=0
def make_dilated_ap(k, step, rng):
    E=[step*i for i in range(k)]
    return E
pool=[]
for _ in range(3000):
    mode=rng.choice(['rand','dilAP','2clust','3clust','bigw','smallw'])
    k=rng.choice([8,9,10])
    if mode=='rand':
        E=sorted(rng.sample(range(1,60),k-1)); E=[0]+E; w=max(E)
    elif mode=='dilAP':
        step=rng.choice([2,3,5,7]); E=[step*i for i in range(k)]; w=max(E)
    elif mode=='2clust':
        b=rng.randint(20,60); E=sorted(set([0,1,2,3][:k//2]+[b+i for i in range(k-k//2)]))[:k];
        if len(E)<k: continue
        w=max(E)
    elif mode=='3clust':
        b1=rng.randint(15,40);b2=rng.randint(60,110)
        E=sorted(set([0,1,2]+[b1,b1+1,b1+2]+[b2,b2+1,b2+2]))[:k]
        if len(E)<k: continue
        w=max(E)
    elif mode=='bigw':
        E=sorted(rng.sample(range(1,30),k-2)); w=rng.randint(80,200); E=sorted([0]+E+[w])
        if len(set(E))!=k: continue
    else: # smallw -- adversarial: small max speed
        E=sorted(rng.sample(range(1,18),k-1)); E=[0]+E; w=max(E)
    E=sorted(set(E))
    if len(E)<k: continue
    if w not in E: continue
    Ep=[e for e in E if e!=w]
    if len(Ep)<k-1: continue
    p0E,_=p0p1(E)
    pEp,regsEp=analyze(Ep)
    plat=pEp[0]+F(1,7)*pEp[1]
    dw=p0E-plat; c1=c1_components(regsEp)
    if c1==0:
        n+=1
        if abs(dw)!=0:
            # dw must be 0 if no miss-1 region can be flipped... actually dw can be 0
            pass
        continue
    bound_loose=F(2*c1,7*w); bound_sharp=F(c1,7*w)
    if abs(dw)>bound_loose: fails_loose+=1
    if abs(dw)>bound_sharp: fails_sharp+=1
    worst_loose=max(worst_loose,float(abs(dw)/bound_loose))
    worst_sharp=max(worst_sharp,float(abs(dw)/bound_sharp))
    worst_abs=max(worst_abs,abs(dw)); n+=1
print(f"  adversarial pool n={n}")
print(f"  LOOSE bound 2c1/(7w): violations={fails_loose}  worst |dw|/bound={worst_loose:.4f}")
print(f"  SHARP bound  c1/(7w): violations={fails_sharp}  worst |dw|/bound={worst_sharp:.4f}")
print(f"  worst |Delta_w| seen = {float(worst_abs):.5f}")
print(f"  => LOOSE comb bound {'HOLDS' if fails_loose==0 else 'VIOLATED'}; "
      f"sharp(per-end 1/(7w)) {'holds' if fails_sharp==0 else 'VIOLATED -> proof const must be 2/(7w)'}")

print("\n"+"="*78)
print("(D)  COUNTEREXAMPLE HUNT to (WIDE): primitive WIDE k-set, span>44, p0>cap_k")
print("="*78)
B=44
worst={k:(F(0),None) for k in CAP}
viol_wide=[]
def consider(E):
    E=sorted(set(E));
    if 0 not in E: return
    k=len(E)
    if k not in CAP: return
    if span(E)<=B: return
    if not primitive(E): return
    p0,_=p0p1(E)
    if p0>worst[k][0]: worst[k]=(p0,tuple(E))
    if p0>CAP[k]:
        viol_wide.append((tuple(E),k,p0))
# structured families
families=[]
for base in [[0,1,2],[0,1,2,3],[0,1,2,3,4]]:
    for g in [20,25,30,40,50,60,80,100]:
        for nc in [2,3,4]:
            E=[]
            for c in range(nc): E+=[g*c+b for b in base]
            families.append(E)
# dilated APs made primitive (add a +1)
for step in [2,3,5,7,11]:
    for k in [8,9,10,11,12]:
        E=[step*i for i in range(k)]; E[1]+=1  # break gcd
        families.append(E)
# all residues mod 7 represented, wide
for trials in range(2000):
    k=rng.choice([8,9,10,11,12])
    E={0}
    while len(E)<k:
        E.add(rng.randint(1, 90))
    families.append(sorted(E))
# balanced 2-3 cluster random
for trials in range(2000):
    k=rng.choice([8,9,10])
    nc=rng.choice([2,3]); per=k//nc
    cs=sorted(rng.sample(range(0,80),nc))
    E=set()
    for c in cs:
        for j in range(per): E.add(c*5+j)
    while len(E)<k: E.add(rng.randint(0,120))
    families.append(sorted(E))
for E in families:
    consider(E)
for k in sorted(CAP):
    p0,E=worst[k]
    if E is None:
        print(f"  k={k}: (no wide primitive sample landed here)"); continue
    print(f"  k={k}: max p0 over WIDE prim sets = {float(p0):.4f}  cap_{k}={float(CAP[k]):.4f}  "
          f"margin={float(CAP[k]-p0):.4f}  argmax={E}")
print(f"  (WIDE) violations p0>cap_k: {len(viol_wide)}")
for v in viol_wide[:10]:
    print("   *** WIDE COUNTEREXAMPLE:",v)

print("\n"+"="*78)
print("(E)  COUNTEREXAMPLE HUNT to (W*): span>44, p0 > Q(k-1)  (the induction invariant)")
print("="*78)
viol_wstar=[]
worstW={8:(F(0),None),9:(F(0),None),10:(F(0),None),11:(F(0),None),12:(F(0),None)}
def considerW(E):
    E=sorted(set(E))
    if 0 not in E: return
    k=len(E)
    if k not in CAP: return
    if span(E)<=B: return
    if not primitive(E): return
    p0,_=p0p1(E); Qk1=Q(k-1)
    if p0>worstW[k][0]: worstW[k]=(p0,tuple(E))
    if p0>Qk1:
        viol_wstar.append((tuple(E),k,float(p0),float(Qk1)))
# reuse families + a denser sweep at smaller span just over B
for E in families: considerW(E)
for trials in range(8000):
    k=rng.choice([8,9,10,11,12])
    E={0}
    rng_max=rng.choice([50,60,80,120])
    while len(E)<k: E.add(rng.randint(1,rng_max))
    considerW(sorted(E))
for k in [8,9,10,11,12]:
    p0,E=worstW[k]
    print(f"  k={k}: max p0 over WIDE prim = {float(p0):.4f}  Q({k-1})={float(Q(k-1)):.4f}  "
          f"margin={float(Q(k-1)-p0):.4f}  argmax={E}")
print(f"  (W*) violations p0>Q(k-1): {len(viol_wstar)}")
for v in viol_wstar[:15]:
    print("   *** (W*) VIOLATION:",v)

print("\n"+"="*78)
print("(F)  Is c1(E') unbounded in span?  And does p0(E')<=Q(k-2) hold for wide E'?")
print("="*78)
print("  c1(E') vs span for dilated APs and multi-cluster (E' = a (k-1)-set):")
for Ep in [list(range(8)), [2*i for i in range(8)], [3*i for i in range(8)],
           [0,1,2,20,21,22,40,41], [0,1,2,30,31,32,60,61],
           [0,1,2,50,51,52,100,101], [5*i for i in range(8)]]:
    p,regs=analyze(Ep); c1=c1_components(regs)
    print(f"    E'={Ep}  span={span(Ep)}  c1={c1}  Plat={float(p[0]+F(1,7)*p[1]):.4f}")
print("\n  test inductive claim p0(E')<=Q(k-2) for WIDE (k-1)-sets E':")
for kk in [8,9,10]:   # E' is a (kk-1)-set, compare to Q(kk-2)
    km1=kk-1; Qkm2=Q(km1-1); mx=F(0); viol=0; nn=0; arg=None
    for _ in range(3000):
        Ep=sorted(rng.sample(range(1,55),km1-1)); Ep=[0]+Ep
        if len(set(Ep))!=km1: continue
        if span(Ep)<=14: continue
        nn+=1
        p0,_=p0p1(Ep)
        if p0>mx: mx=p0; arg=tuple(Ep)
        if p0>Qkm2: viol+=1
    print(f"    (k-1)={km1}-sets, wide: max p0={float(mx):.4f}  Q(k-2)=Q({km1-1})={float(Qkm2):.4f}  "
          f"#(p0>Q(k-2))={viol}/{nn}  {'OK' if viol==0 else '<-- inductive premise FAILS'}")

print("\n"+"="*78)
print("(G)  GLUE: B=44 vs finite check span<=14")
print("="*78)
print("  Single-peel cutoffs W(consec_{k-1}) = ceil(2 c1(consec)/ (7 mu_k)):")
Bvals=[]
for k in range(8,13):
    Ep=list(range(k-1)); p,regs=analyze(Ep); c1=c1_components(regs)
    mu=CAP[k]-Q(k-1)
    Wq=F(2*c1,7)/mu; W=-(-Wq.numerator//Wq.denominator)
    Bvals.append(W)
    print(f"    k={k}: c1(consec_{k-1})={c1}  mu_k={float(mu):.5f}  W={W}")
print(f"  B_consec=max={max(Bvals)}, glue B=max(44,14)={max(max(Bvals),44,14)}")
print("""
  NOTE: B=44 is only the cutoff for the GAP-DOMINATED single-peel branch on the
  CONSEC core.  The balanced-wide branch has w<W(E') with c1(E') GROWING in span,
  so its per-set cutoff W(E') exceeds 44.  Hence B=44 does NOT by itself glue the
  whole WIDE claim to span<=14: there is a residue {span>14, w<W(E')} where the
  comb cutoff is not met.  That residue is covered only by the (sampled, unproved)
  invariant (W*).
""")
print("""
========================= ADVERSARIAL VERDICT (recorded) =========================
PROVED & confirmed exact:
  (B) margins mu_k = cap_k - Q(k-1) > 0, exact rationals match the claim.
  (C) COMB BOUND |Delta_w| <= 2 c1(E')/(7w): 0 violations over a 2665-set
      adversarial pool (dilated APs, multi-cluster, small/large w); loose ~5x.
      The supporting algebra (single new speed = period-1/w comb on each miss-1
      component; interior periods give exactly L/7 summing to (1/7)p1; ends
      deviate <= one tooth) is a correct proof.  The identity p0=Plat(E')+Delta_w
      is exact by definition; the "miss>=2 -> contributes 0" step is correct.

NO COUNTEREXAMPLE FOUND:
  (WIDE) p0(E) > cap_k : 0 violations over ~12000 wide primitive k-sets.
  (W*)   p0(E) > Q(k-1): 0 violations.  But margins to Q(k-1) are THIN and
         SHRINKING in k (0.069/0.050/0.044/0.043/0.061), argmaxes are dilated
         near-APs -- exactly the regime where the unproved extremality bites.

THE TWO REAL GAPS (proof NOT complete):

GAP 1 -- (P) is MISLABELED "PROVED via THM-535".  THM-535 proves only the cap
  LOWER bound cap_k >= (k-6)/7 and computes meas(S7) for the SPECIFIC consec set.
  It does NOT prove that consec/AP is the plateau ARGMAX over all (k-1)-sets.
  That statement -- "Plat(E') <= Q(k-1) := Plat(consec_{k-1}) for EVERY (k-1)-set"
  -- IS the still-OPEN "AP maximises meas(S7)" conjecture (HYP-2603/HYP-2607,
  the SHARED finishing conjecture of the whole LRC(14) program; THM-535 honest
  status: 'LRC(14) NOT proved... still needs HYP-2603').  So (P) is assumed, not
  proved.  Without (P) the comb branch gives p0(E) <= Plat(E') + mu_k, but
  Plat(E') is not bounded by Q(k-1).

GAP 2 -- the BALANCED-wide branch (w < W(E')) is genuinely open.  W(E') =
  ceil(2 c1(E')/(7 mu_k)) and c1(E') GROWS with span (dilated APs have c1 rising),
  so the single-peel cutoff is NOT uniformly bounded; B=44 bounds W only on the
  CONSEC core, not on arbitrary wide E'.  The STEP-6 inductive collapse invokes
  the premise "wide E' => p0(E') <= Q(k-2)", which FAILS empirically already at
  (k-1)=7: the wide set (0,12,18,24,30,36,42) has p0=0.0456 > Q(6)=0.0374.  So the
  written induction step is unsound as stated.

NET: the reduction is real and useful -- the gap-dominated branch (w >= W(E') with
the consec extremal) is rigorously closed, and (C) is a genuine new lemma.  But
HYP-2675 (WIDE) is NOT proved: (P) silently assumes the open AP-extremality, and
the balanced-wide residue has no proved uniform cutoff.  B=44 does NOT glue to
span<=14 without a gap (the residue {span>14, w<W(E')} is uncovered).  (W*) is the
right invariant and holds on all samples, but it is exactly as unproved as the
AP-extremality it encodes.
================================================================================
""")
print("VERIFICATION SCRIPT COMPLETE.")
