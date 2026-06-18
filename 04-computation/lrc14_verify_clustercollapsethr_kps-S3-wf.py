#!/usr/bin/env python3
"""
lrc14_verify_clustercollapsethr_kps-S3-wf

ADVERSARIAL VERIFICATION of the claimed "cluster-collapse-threedistance" advance on
the S3 residual case of LRC(14).

The claim package (to verify / refute):
  (A) CLUSTER-COLLAPSE LEMMA (PROVED): for finite L subset [V0,Vmax] with 13*V0 > Vmax,
      J0 = (1/(14 V0), 1/Vmax - 1/(14 Vmax)) is level-1/14 safe for every u in L, and
      width(J0)*7*Vmax = (13 - Vmax/V0)/2, which > 1 iff Vmax < 11*V0.
  (B) CRITERION C(S) => M(S) >= 1/14 (arc-width lemma, given as PROVED upstream).
  (C) CARRY-PHASE LIMIT REDUCTION: margin_limit infimum = 1 (exact, tight). Does NOT close S3.
  (D) Sub-claims:
       - v=max determinism rule (REFUTED already by them; re-verify counterexample).
       - "realized global min best-margin ~ 1.336", "extremal-shape ~1.9002".
       - the M(S) for v=max counterexample = 123/989 ~ 0.124368.

This script:
  1. Independently re-implements W(A), M(A), C(S) with EXACT Fractions and a wrap-correct
     safe-arc routine, and cross-checks against the prompt's provided tools.
  2. Re-derives the cluster-collapse lemma inequalities from scratch on random L's.
  3. HUNTS for covering+primitive S3 13-sets that violate C(S) (i.e. C(S) false) or have M(S)<1/14.
  4. Re-derives the v=max counterexample and checks whether C(S) (over ALL v, not just v=max) holds.
  5. Stress-tests the "realized best-margin floor ~1.33" claim adversarially.
"""
from fractions import Fraction as F
from math import gcd
import random
import itertools

H = F(1,14)

# ---------------- EXACT tools (independent re-implementation) ----------------

def nrm(x):
    r = x - int(x)
    if r < 0: r += 1
    return r if r <= F(1,2) else 1 - r

def safe_components(A, h=H):
    """Return list of (lo,hi) maximal safe arcs where ||a*tau||>=h for all a in A. Wrap-correct."""
    iv = []
    for u in A:
        for j in range(0, u):
            c = F(j, u)
            a = (c - h/u) % 1
            b = (c + h/u) % 1
            if a < b:
                iv.append((a, b))
            else:
                iv.append((a, F(1)))
                iv.append((F(0), b))
    iv.sort()
    merged = []
    for a, b in iv:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    safe = []
    prev = F(0)
    for a, b in merged:
        if a > prev:
            safe.append((prev, a))
        prev = max(prev, b)
    if prev < 1:
        safe.append((prev, F(1)))
    return safe

def Wwidth(A):
    sc = safe_components(A)
    if not sc:
        return F(0)
    ws = [b - a for a, b in sc]
    # circular wrap: if first arc touches 0 and last touches 1, they join
    if sc[0][0] == 0 and sc[-1][1] == 1 and len(sc) > 1:
        ws.append((sc[0][1]) + (1 - sc[-1][0]))
    return max(ws)

def cand(S):
    S = sorted(set(S)); C = set()
    for v in S:
        k = 0
        while F(2*k+1, 2*v) <= F(1,2):
            C.add(F(2*k+1, 2*v)); k += 1
    for i in range(len(S)):
        for j in range(i+1, len(S)):
            for d in (S[i]+S[j], S[j]-S[i]):
                if d > 0:
                    k = 1
                    while F(k, d) <= F(1,2):
                        C.add(F(k, d)); k += 1
    C.add(F(1,2)); return C

def Mval(S):
    b = F(0)
    for t in cand(S):
        v = min(nrm(x*t) for x in S)
        if v > b: b = v
    return b

def is_covering(S):
    return all(any(v % q == 0 for v in S) for q in range(2,15))

def is_primitive(S):
    g = 0
    for v in S:
        g = gcd(g, v)
    return g == 1

def criterion_C(S):
    """C(S): exists v in S with W(S\{v}) > 1/(7v). Returns (holds, best_v, best_margin)."""
    S = sorted(set(S))
    best_margin = F(-1)
    best_v = None
    holds = False
    for v in S:
        A = [u for u in S if u != v]
        w = Wwidth(A)
        margin = w * 7 * v   # W > 1/(7v) <=> margin > 1
        if margin > best_margin:
            best_margin = margin; best_v = v
        if w > F(1, 7*v):
            holds = True
    return holds, best_v, best_margin

# =================================================================
print("="*70)
print("PART 0: Tool self-consistency checks")
print("="*70)

# Sanity: a known covering set; check M, C consistency
# Build a simple S3 set: small part {1..p} + cluster + Vmax
def make_s3(small, V0, offsets, Vmax):
    L = [V0 + d for d in offsets]
    S = sorted(set(list(small) + L + [Vmax]))
    return S

# A clustered example
S_test = [1,2,3,4,5,6,7,8,9,10,11,12,13]  # plain {1..13} (covering, not S3 but valid covering)
print("S={1..13}: covering?", is_covering(S_test), "primitive?", is_primitive(S_test))
print("  M(S)=", Mval(S_test), float(Mval(S_test)), "  >= 1/14?", Mval(S_test) >= H)
hC, bv, bm = criterion_C(S_test)
print("  C(S):", hC, " best_v=", bv, " best_margin=", bm, float(bm))

# =================================================================
print()
print("="*70)
print("PART 1: Re-derive CLUSTER-COLLAPSE LEMMA from scratch (exact)")
print("="*70)
# Claim: L subset [V0,Vmax], 13*V0>Vmax => J0=(1/(14V0), 1/Vmax - 1/(14Vmax)) safe for all u in L
# and width(J0)*7*Vmax = (13 - Vmax/V0)/2.
def check_cluster_collapse(L):
    V0 = min(L); Vmax = max(L)
    assert 13*V0 > Vmax, "precondition fails"
    lo = F(1, 14*V0)
    hi = F(1, Vmax) - F(1, 14*Vmax)   # = 13/(14 Vmax)
    # verify hi simplifies to 13/(14 Vmax)
    assert hi == F(13, 14*Vmax)
    assert lo < hi, "J0 empty"
    width = hi - lo
    # claimed width*7*Vmax = (13 - Vmax/V0)/2
    claimed = (F(13) - F(Vmax, V0)) / 2
    lhs = width * 7 * Vmax
    ok_width = (lhs == claimed)
    # verify safety: for u in L, and tau in (lo,hi), is u*tau in (1/14,13/14)?
    # The lemma argues 0<u*tau<1 (so frac=u*tau) and u*tau in (1/14,13/14).
    # Check endpoints rigorously: inf over u of u*lo, sup over u of u*hi.
    # For tau just above lo: u*lo >= V0*lo = 1/14 (with u>=V0). min at u=V0.
    # For tau just below hi: u*hi <= Vmax*hi = 13/14. max at u=Vmax.
    # But need u*tau < 1 for ALL u: max is Vmax*hi = 13/14 < 1. OK.
    # And need u*tau>1/14 strictly: V0*lo=1/14 exactly => need tau>lo strict (open interval). OK.
    # Pick a rational test point and verify all safe
    mid = (lo+hi)/2
    all_safe = all(nrm(u*mid) >= H for u in L)
    # Also verify the LOWER bound on u*tau holds for the smallest u over the WHOLE open interval:
    # smallest product = V0*lo = 1/14 (boundary), largest = Vmax*hi=13/14 (boundary). Interior strict.
    return ok_width, all_safe, width, claimed, lhs

random.seed(12345)
cc_fail = 0
for trial in range(20000):
    V0 = random.randint(20, 5000)
    spread = random.randint(1, min(12*V0-1, 60))  # keep 13V0>Vmax
    c = random.randint(2, 11)
    offs = sorted(random.sample(range(0, spread+1), min(c, spread+1)))
    if 0 not in offs: offs[0] = 0
    if spread not in offs: offs[-1] = spread
    L = sorted(set(V0 + d for d in offs))
    if len(L) < 2: continue
    if 13*min(L) <= max(L): continue
    okw, oks, w, claimed, lhs = check_cluster_collapse(L)
    if not (okw and oks):
        cc_fail += 1
        if cc_fail <= 5:
            print("  CLUSTER-COLLAPSE FAIL:", L, "okw",okw,"oks",oks,"lhs",lhs,"claimed",claimed)
print(f"Cluster-collapse lemma: {cc_fail} failures in 20000 random L's (exact)")

# Verify the threshold claim: width*7*Vmax > 1  iff  Vmax < 11*V0
print("\nThreshold (13 - Vmax/V0)/2 > 1  <=>  Vmax < 11*V0 :")
for (V0,Vmax) in [(100,1099),(100,1100),(100,1101),(100,1000),(100,1300)]:
    val = (F(13)-F(Vmax,V0))/2
    print(f"  V0={V0} Vmax={Vmax}: (13-Vmax/V0)/2={val}={float(val):.4f}  Vmax<11V0={Vmax<11*V0}  val>1={val>1}")

# =================================================================
print()
print("="*70)
print("PART 2: S3 set generation + adversarial hunt for C(S) failures / M<1/14")
print("="*70)
# S3: k>=2 large (>13) AND Vmax >= 13*Vmin.  Vmin=min(S).
# A covering 13-set must hit every q in 2..14. We need primitive + covering.

def case_of(S):
    S=sorted(set(S))
    Vmin=S[0]; Vmax=S[-1]
    k=sum(1 for v in S if v>13)
    if k<=1: return "S1"
    if Vmax < 13*Vmin: return "S2"
    return "S3"

def random_covering_13set(rng, max_speed=400):
    """Try to build a primitive covering 13-set, biased toward S3 (small part + large cluster)."""
    # strategy: pick small part subset of {1..13}, plus a large cluster near some V, plus maybe Vmax
    # We need covering: for each q in 2..14 some element divisible by q.
    for _ in range(200):
        S=set()
        # ensure small divisors via small runners or large multiples
        small_n = rng.randint(4, 9)
        small = rng.sample(range(1,14), small_n)
        S.update(small)
        V = rng.randint(20, max_speed)
        cspread = rng.randint(2, 45)
        csize = rng.randint(2, 11 - 0)
        cl = set()
        cl.add(V)
        while len(cl) < csize:
            cl.add(V + rng.randint(1, cspread))
        S.update(cl)
        # top off to 13 elements with more larges or smalls
        while len(S) < 13:
            if rng.random()<0.5:
                S.add(rng.randint(1,13))
            else:
                S.add(V + rng.randint(0, cspread+10))
        # trim to 13
        S=set(list(sorted(S))[:13])
        if len(S)!=13: continue
        S=sorted(S)
        if not is_primitive(S): continue
        if not is_covering(S): continue
        return S
    return None

rng=random.Random(999)
n_s3=0; n_cfail=0; n_Mfail=0
worst_margin=F(10**9); worst_margin_set=None
cfail_examples=[]; Mfail_examples=[]
TRIALS=60000
for t in range(TRIALS):
    S=random_covering_13set(rng)
    if S is None: continue
    if case_of(S)!="S3": continue
    n_s3+=1
    hC,bv,bm=criterion_C(S)
    if bm < worst_margin:
        worst_margin=bm; worst_margin_set=S[:]
    if not hC:
        n_cfail+=1
        if len(cfail_examples)<10: cfail_examples.append((S[:],bm))
        # If C fails, check actual M
        m=Mval(S)
        if m < H:
            n_Mfail+=1
            if len(Mfail_examples)<10: Mfail_examples.append((S[:],m))

print(f"S3 sets found: {n_s3}")
print(f"C(S) FAILURES (no v gives margin>1): {n_cfail}")
print(f"M(S) < 1/14 FAILURES (true LRC counterexamples): {n_Mfail}")
print(f"Worst (min) best-margin over all S3: {worst_margin} = {float(worst_margin):.6f}")
print(f"  achieved by: {worst_margin_set}")
if cfail_examples:
    print("C-failure examples (set, best_margin):")
    for s,b in cfail_examples[:10]:
        print("   ", s, float(b), " M=", float(Mval(s)))
if Mfail_examples:
    print("!!! M<1/14 examples:")
    for s,m in Mfail_examples: print("   ", s, m, float(m))

# =================================================================
print()
print("="*70)
print("PART 3: v=max determinism sub-claim + the claimed v=max counterexample")
print("="*70)
# Claim: there is an S3 set where v=max gives margin 242015/245024 ~0.98772 (<1, so v=max FAILS),
# but C(S) still HOLDS via some other v, and M(S)=123/989~0.124368 (>1/14).
# Re-derive: search for an S3 covering set where the v=max margin < 1 but some other v gives margin>1.
print("claimed v=max margin 242015/245024 =", float(F(242015,245024)))
print("claimed M = 123/989 =", float(F(123,989)), " >= 1/14?", F(123,989)>=H)

rng2=random.Random(7)
found_vmax_fail=0
examples=[]
for t in range(80000):
    S=random_covering_13set(rng2)
    if S is None: continue
    if case_of(S)!="S3": continue
    S=sorted(S)
    vmax=S[-1]
    A=[u for u in S if u!=vmax]
    wmax=Wwidth(A)
    margin_vmax=wmax*7*vmax
    if margin_vmax<=1:  # v=max does NOT fire
        # does some OTHER v fire?
        hC,bv,bm=criterion_C(S)
        found_vmax_fail+=1
        if len(examples)<8:
            examples.append((S[:],margin_vmax,hC,bv,bm,Mval(S)))
print(f"S3 sets where v=max margin <= 1 (v=max determinism FALSE): {found_vmax_fail}")
for s,mv,hC,bv,bm,M in examples:
    print(f"   S={s}")
    print(f"     v=max margin={float(mv):.5f}  C(S) via other v: {hC} best_v={bv} best_margin={float(bm):.5f}  M={M}={float(M):.6f}")

# =================================================================
print()
print("="*70)
print("PART 4: CARRY-PHASE LIMIT — is infimum really 1, and does C(S) hold at finite?")
print("="*70)
# The crux: claim says limit infimum margin = 1 (tight), realized floor ~1.33.
# We probe: take a FIXED shape (small part, offsets, removed gap D) and let V0 grow.
# Track best-margin = max_v W(S\{v})*7*v.  Does it converge from ABOVE to a limit?  Is limit>1?
# And critically: is there ANY finite covering+primitive S3 set with best-margin <= 1 (=> C fails)?
# (We already gather worst_margin in Part 2.)

# Build a parametric family and watch best-margin vs V0.
def build_shape(small, offsets, D, V0):
    # cluster = {V0+o for o in offsets}, plus Vmax=V0+D (D>spread => it's the "removed large")
    L=[V0+o for o in offsets]
    Vmax=V0+D
    S=sorted(set(list(small)+L+[Vmax]))
    return S

# pick a shape that is S3 and covering at multiple V0; need covering => small part must supply 2..14 hits
# small={1,..} won't be divisible by big q's unless cluster supplies them. Use small with required divisors.
# Try: small part {1,2,3,4,5,6,7} hits 2,3,4,5,6,7; need 8,9,10,11,12,13,14 from cluster/Vmax via divisibility.
# That is hard for arbitrary V0. Instead, just track best-margin numerically over many V0 for FIXED offset pattern
# ignoring covering, to test the "converges from above to limit" claim about the geometry.
print("Geometry-only convergence test (fixed offsets, growing V0):")
offsets=[0,14,28]; D=45; small=[1,2,3,4,5,6,7,8,9,10,11]
# trim to size 13 total: small(11)+cluster(3)+Vmax(1) = 15 -> use small of size 9
small=[1,2,3,4,5,6,7,8,9]
prev=None
for V0 in [50,100,200,400,800,1600,3200,6400,12800,25600]:
    S=build_shape(small,offsets,D,V0)
    if len(S)!=13:
        # dedup may shrink; pad
        pass
    hC,bv,bm=criterion_C(S)
    arrow=""
    if prev is not None:
        arrow = "down" if bm<prev else ("up" if bm>prev else "flat")
    print(f"  V0={V0:6d}: best-margin={float(bm):.6f} ({arrow})  C={hC}  case={case_of(S)} |S|={len(S)}")
    prev=bm
