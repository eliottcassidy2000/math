#!/usr/bin/env python3
"""
ADVERSARIAL VERIFICATION (kind-pasteur-2026-06-18-S4) of THM-527
"fixed-small-part-equidistribution" claimed closure of the fixed-small-part
single-tight-cluster sub-case of S3.

I re-derive every load-bearing step INDEPENDENTLY (my own helpers, not theirs),
then HUNT for a covering primitive S3 set in the claimed sub-family with M<1/14,
and stress-test the explicit threshold V0* and the witness construction.

ALL EXACT rationals.
"""
from fractions import Fraction as F
from math import gcd, ceil
from functools import reduce
from itertools import combinations
import random
random.seed(20260618)

# ---------------- exact tools (verbatim from prompt) ----------------
def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r
def safe_components(A,h=F(1,14)):
    iv=[]
    for u in A:
        for j in range(0,u):
            c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
            if a<b: iv.append((a,b))
            else: iv.append((a,F(1))); iv.append((F(0),b))
    iv.sort(); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C
def Mval(S):
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b
def is_cov(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def meas_safe(A): return sum(b-a for a,b in safe_components(A))

def cls(S):
    S=sorted(S); k=sum(1 for v in S if v>13)
    if k<=1: return 'S1'
    if S[-1]<13*S[0]: return 'S2'
    return 'S3'

# ====================================================================
# PART 0 — verify the FUNDAMENTAL algebraic identity the lemma rests on:
#   for u=V0+d_i, set theta=V0*tau.  Claim: ||u*tau||>=1/14  <=>  frac(theta+d_i*tau) in [1/14,13/14].
#   (i.e. the "fast phase" decomposition is exact.)
# ====================================================================
print("="*80)
print("PART 0 — algebraic identity: ||(V0+d)tau||>=1/14 <=> frac(V0 tau + d tau) in [1/14,13/14]")
print("="*80)
bad0=0; tot0=0
rng=random.Random(1)
for _ in range(3000):
    V0=rng.randint(14,500); d=rng.randint(0,40)
    tau=F(rng.randint(0,10**6),10**6)
    u=V0+d
    lhs = (nrm(u*tau)>=F(1,14))
    fr=(V0*tau + d*tau)%1   # = u*tau mod 1
    rhs = (F(1,14)<=fr<=F(13,14))
    tot0+=1
    if lhs!=rhs: bad0+=1
print("  random checks: %d ; identity mismatches: %d"%(tot0,bad0))
print("  (||x||>=1/14 iff frac(x) in [1/14,13/14] -- a definitional tautology; sanity only)")

# ====================================================================
# PART 1 — RE-DERIVE the witness lemma INDEPENDENTLY and build a GLOBAL WITNESS
# from scratch for a given (P, offsets, V0).  The lemma's TRUE content:
#   exists tau with (a) tau in G_P  and (b) frac((V0+d_i)tau) in [1/14,13/14] for all i.
# I implement a DIRECT, brute search for such a witness over the candidate set used by Mval,
# AND an independent "sweep" construction, then compare with Mval(S).  If Mval(S)>=1/14 always
# on the claimed sub-family, the closure conclusion holds (regardless of HOW the witness is found).
# ====================================================================

def cluster_safe_theta_arc(offsets, tau):
    """EXACT: set of theta=frac(V0 tau) making ALL cluster speeds safe, i.e.
       {theta: frac(theta+d_i*tau) in [1/14,13/14] for all i}.  Return list of (lo,hi) arcs in [0,1)."""
    # frac(theta+c) in [1/14,13/14] <=> theta in [1/14-c, 13/14-c] (mod1), an arc of width 6/7.
    arcs=[]
    for d in offsets:
        c=(d*tau)%1
        lo=(F(1,14)-c)%1; hi=(F(13,14)-c)%1
        if lo<=hi: arcs.append([(lo,hi)])
        else: arcs.append([(lo,F(1)),(F(0),hi)])
    # intersect all arc-sets
    cur=[(F(0),F(1))]
    for aset in arcs:
        nxt=[]
        for (l0,h0) in cur:
            for (l1,h1) in aset:
                lo=max(l0,l1); hi=min(h0,h1)
                if lo<hi: nxt.append((lo,hi))
        cur=nxt
        if not cur: break
    return cur

def witness_exists_direct(P, offsets, V0):
    """Independent existence test of a GLOBAL witness tau (exact) at a discrete but EXACT
       candidate grid: we use Mval's candidate set for the FULL S, which is exactly the set of
       tau where M is attained (M is a max over those candidates).  Returns (Mval(S), witness or None)."""
    S=sorted(set(list(P)+[V0+d for d in offsets]))
    best=F(0); wt=None
    for t in cand(S):
        mn=min(nrm(x*t) for x in S)
        if mn>best: best=mn; wt=t
    return best, wt

# ====================================================================
# PART 2 — HUNT.  Build covering+primitive S3 sets in the claimed sub-family and check M>=1/14.
#   Sub-family: S=P U {V0+d_i}, P subset {1..13} with M(P)>1/14, offsets fixed, single tight cluster.
#   We sweep MANY (P,offsets) and MANY V0 (incl. small V0 below any plausible threshold AND large).
# ====================================================================
print()
print("="*80)
print("PART 2 — HUNT for covering primitive S3 set in the sub-family with M<1/14")
print("="*80)

# enumerate proper P with M(P)>1/14, sizes 2..5 (so |L|=8..11, genuine large clusters).
# Keep ALL of sizes 2,3 (smaller P => larger cluster => the bad/dense regime); SAMPLE sizes 4,5
# to keep runtime feasible while still covering the space.
good_P=[]
rngP=random.Random(7)
for sz in range(2,6):
    allP=[list(P) for P in combinations(range(1,14),sz) if meas_safe(list(P))>0]
    if sz<=3:
        good_P+=allP
    else:
        good_P+=rngP.sample(allP, min(120, len(allP)))
print("  #proper P used (sizes 2,3 exhaustive; 4,5 sampled<=120 each), M(P)>1/14: %d"%len(good_P), flush=True)

fails=[]
tested=0; cov_tested=0; s3_cov_tested=0
minM_overall=None; minM_set=None
minM_s3cov=None; minM_s3cov_set=None

# offset patterns: contiguous, arithmetic, and random-bounded
def gen_offset_patterns(nL, rng):
    pats=[]
    pats.append(list(range(nL)))                       # contiguous 0..nL-1 (densest, the bad shape)
    pats.append([0]+sorted(rng.sample(range(1,3*nL),nL-1)))  # random spread (moderate)
    pats.append([2*i for i in range(nL)])              # even AP
    pats.append([0]+sorted(rng.sample(range(1,6*nL),nL-1)))  # wider spread (bounded to keep speeds modest)
    return [p for p in pats if len(set(p))==nL]

rng2=random.Random(42)
# V0 range: deliberately include SMALL V0 (where covering can hold but threshold may not be met)
# AND larger ones.  Representative (not every integer) to keep runtime feasible; cap magnitude so
# exact Mval (cand-set size ~ max speed) stays fast -- the dense-cluster bad shape lives at SMALL V0.
V0_list = list(range(14,160)) + [180,210,252,300,420,560]

pcount=0
for P in good_P:
    pcount+=1
    if pcount%200==0: print("    ...P progress %d/%d (cov_tested=%d s3cov=%d fails=%d)"%(pcount,len(good_P),cov_tested,s3_cov_tested,len(fails)), flush=True)
    nL=13-len(P)
    for offs in gen_offset_patterns(nL, rng2):
        dc=offs[-1]
        for V0 in V0_list:
            # tight cluster requires Vmax < 13*Vmin for the LARGE part to be a cluster;
            # but the S3 membership is about the FULL set. We just classify each set honestly.
            L=[V0+d for d in offs]
            S=sorted(set(list(P)+L))
            if len(S)!=13: continue
            if reduce(gcd,S)!=1: continue   # primitive only
            if not is_cov(S): continue      # covering only
            cov_tested+=1
            c=cls(S)
            if c!='S3': continue
            s3_cov_tested+=1
            M=Mval(S)
            if minM_s3cov is None or M<minM_s3cov:
                minM_s3cov=M; minM_s3cov_set=S
            if M<F(1,14):
                fails.append((P,offs,V0,M,S))
            tested+=1
print("  covering+primitive sets built: %d ; of which S3: %d"%(cov_tested,s3_cov_tested))
print("  S3 covering primitive sets with M<1/14: %d"%len(fails))
if minM_s3cov is not None:
    print("  MIN M over S3 covering primitive sub-family sets tested: %s = %.6f"
          %(minM_s3cov,float(minM_s3cov)))
    print("     attained at S=%s (class %s)"%(minM_s3cov_set,cls(minM_s3cov_set)))
for f in fails[:10]:
    print("    *** FAIL:",f[:4])

# ====================================================================
# PART 3 — Verify the EXPLICIT threshold V0* claim independently.
#   For a fixed (P,offsets), independently compute the constant-sign sub-arc and the
#   sweep threshold, then CHECK: for ALL V0>=V0* (covering or not), a witness exists AND M(S)>=1/14;
#   and probe V0 just BELOW the claimed scope to see whether sets there can violate.
# ====================================================================
print()
print("="*80)
print("PART 3 — explicit threshold V0*: does V0>=V0* really force M(S)>=1/14? (independent)")
print("="*80)

# Independent constructive sub-arc + threshold:
def my_cluster_theta_arc_over_interval(offsets, ap, bp):
    """theta arcs common to cluster_safe_theta_arc(offsets,tau) for ALL tau in [ap,bp].
       For speed u=V0+d, safe-theta arc center moves as tau varies; intersect endpoints is NOT
       enough in general (arc is monotone in tau per d since it's an affine shift of theta by -d*tau).
       The safe-theta arc for offset d is [1/14 - d*tau, 13/14 - d*tau] (mod1): a rigid arc of
       width 6/7 whose CENTER = 1/2 - d*tau moves LINEARLY in tau.  Over tau in [ap,bp] the center
       sweeps an arc of length d*(bp-ap).  theta is in ALL of them for ALL tau iff theta is within
       3/7 of EVERY center in the swept union.  => same as: 6/7 - circ_width(union of swept center-arcs).
       I compute the explicit common arc directly by intersecting the per-d 'safe for all tau in [ap,bp]'
       theta-sets, each of which = arc [1/14 - d*ap, 13/14 - d*bp] (shrunk by d*(bp-ap) on the trailing
       side)... carefully: safe-theta(d,tau) = [1/14 - d*tau, 13/14 - d*tau].
       Intersection over tau in [ap,bp] of a rigidly-translating arc = [max_tau lower, min_tau upper]
       = [1/14 - d*ap, 13/14 - d*bp] (since lower endpoint 1/14 - d*tau is maximized at tau=ap for d>0,
       upper endpoint 13/14 - d*tau is minimized at tau=bp). Width = 6/7 - d*(bp-ap). For d=0 full arc."""
    cur=[(F(0),F(1))]
    for d in offsets:
        width=F(6,7)-d*(bp-ap)
        if width<=0: return []
        lo=(F(1,14)-d*ap)%1
        hi=(lo+width)%1
        if lo<=hi: aset=[(lo,hi)]
        else: aset=[(lo,F(1)),(F(0),hi)]
        nxt=[]
        for (l0,h0) in cur:
            for (l1,h1) in aset:
                L=max(l0,l1); H=min(h0,h1)
                if L<H: nxt.append((L,H))
        cur=nxt
        if not cur: return []
    return cur

def my_threshold(P, offsets):
    """Find a sub-arc [ap,bp] of a P-safe arc on which the common cluster-theta set is nonempty,
       maximize its tau-length crudely-but-exactly, return V0*=ceil(1/(bp-ap)). Returns dict or None."""
    sc=safe_components(P)
    if not sc: return None
    # widest P arc
    alpha,beta=max(sc,key=lambda ab: ab[1]-ab[0])
    # We want largest s=bp-ap s.t. for SOME [ap,bp] subset [alpha,beta], common theta-set nonempty.
    # Common theta-set nonempty requires 6/7 - d_max*(bp-ap) > 0 AND the per-d arcs intersect.
    # Necessary: bp-ap < (6/7)/d_max. Search ap on a fine EXACT grid of breakpoints; for each ap
    # bisect bp up to the nonemptiness boundary.
    dmax=max(offsets) if max(offsets)>0 else 1
    hard_cap=F(6,7)/dmax
    best=None
    # candidate left endpoints: P-arc start + fractions
    NPTS=60
    for kk in range(NPTS):
        ap=alpha+(beta-alpha)*F(kk,NPTS)
        hi=min(beta-ap, hard_cap)
        if hi<=0: continue
        lo=F(0)
        # bisect: predicate nonempty common theta-set on [ap,ap+w]
        for _ in range(60):
            mid=(lo+hi)/2
            if my_cluster_theta_arc_over_interval(offsets,ap,ap+mid):
                lo=mid
            else: hi=mid
        s=lo
        if s>0 and (best is None or s>best[0]):
            best=(s,ap,ap+s)
    if best is None: return None
    s,ap,bp=best
    V0star=ceil(F(1)/s)
    return {'s':s,'ap':ap,'bp':bp,'V0star':V0star}

probe_patterns=[
    ([1,2,3],list(range(10))),
    ([1,2,3],[0,4,6,9,13,15,22,24,27]),
    ([1,2,3,4],[0,2,5,8,11,14,17,20,23]),
    ([1,2,3,5,7],[0,3,7,12,18,25,33,42]),
]
thr_fail=[]
for P,offs in probe_patterns:
    th=my_threshold(P,offs)
    if th is None:
        print("  P=%s offs=%s : NO sub-arc found"%(P,offs)); continue
    V0star=th['V0star']
    # check V0 from V0star up to V0star+200 (sample) : witness must exist & M>=1/14
    bad=0; checked=0; minM=None
    for V0 in list(range(V0star, V0star+60))+[V0star+150, 2*V0star]:
        S=sorted(set(list(P)+[V0+d for d in offs]))
        if len(S)!=13: continue
        M=Mval(S)
        checked+=1
        if minM is None or M<minM: minM=M
        if M<F(1,14):
            bad+=1; thr_fail.append((P,offs,V0,M))
    print("  P=%-12s |L|=%2d V0*=%d : checked %d V0>=V0*, min M=%s=%.5f, M<1/14 count=%d"
          %(str(P),len(offs),V0star,checked,minM,float(minM),bad))

print("  threshold violations (M<1/14 at V0>=V0*): %d"%len(thr_fail))
for f in thr_fail[:10]: print("    *** THRESHOLD FAIL:",f)

# ====================================================================
# PART 4 — Verify the EXACT witness construction gives a TRUE global witness.
#   Use the construction: pick a sub-arc [ap,bp] with common cluster-theta arc Omega*=[tlo,thi],
#   pick theta_c=midpoint, set tau*=(n+theta_c)/V0 for the right n landing in [ap,bp];
#   verify min_{v in S} ||v tau*|| >= 1/14 EXACTLY.
# ====================================================================
print()
print("="*80)
print("PART 4 — exact witness tau*: min_v ||v tau*|| >= 1/14 for actual S (independent build)")
print("="*80)
for P,offs in probe_patterns:
    th=my_threshold(P,offs)
    if th is None: continue
    ap,bp,V0star=th['ap'],th['bp'],th['V0star']
    arc=my_cluster_theta_arc_over_interval(offs,ap,bp)
    if not arc:
        print("  P=%s : empty Omega* (unexpected)"%P); continue
    # pick widest theta-arc
    tlo,thi=max(arc,key=lambda x:x[1]-x[0])
    theta_c=(tlo+thi)/2
    for V0 in [V0star, V0star+5, 3*V0star]:
        # need tau*=(n+theta_c)/V0 in [ap,bp]
        nlo=V0*ap-theta_c; nhi=V0*bp-theta_c
        n=ceil(nlo); taustar=None
        while n<=nhi:
            t=(n+theta_c)/V0
            if ap<=t<=bp: taustar=t; break
            n+=1
        if taustar is None:
            print("  P=%-12s V0=%5d: NO n lands (sweep len=%.3f)"%(str(P),V0,float(V0*(bp-ap)))); continue
        S=sorted(set(list(P)+[V0+d for d in offs]))
        mn=min(nrm(x*taustar) for x in S)
        print("  P=%-12s V0=%5d: min_v||v tau*||=%s=%.5f  >=1/14? %s"
              %(str(P),V0,mn,float(mn),mn>=F(1,14)))

# ====================================================================
# PART 5 — SCOPE sanity: are the sub-family sets actually S3, and does covering FORCE
#   tightness assumptions?  Also confirm the SPECIAL counterexample S* (MISTAKE-076) is NOT
#   in this sub-family (it has a fixed small part but its large part {38,42} is k=2, not a big cluster).
# ====================================================================
print()
print("="*80)
print("PART 5 — scope sanity")
print("="*80)
Sstar=[1,2,3,5,7,8,9,10,11,12,13,38,42]
print("  S* (MISTAKE-076) = %s : covering=%s primitive=%s class=%s M=%s"
      %(Sstar,is_cov(Sstar),reduce(gcd,Sstar)==1,cls(Sstar),Mval(Sstar)))
print("     S* large part {38,42}: k=2 small cluster (not the |L|>=8 big-cluster regime); M=%s>=1/14 OK"%Mval(Sstar))

# Does the sub-family ALWAYS land in S3? check a contiguous-cluster example growing V0.
print("  contiguous cluster P={1,2,3} offs=0..9, class vs V0:")
for V0 in [14,21,28,56,84,140,280,560]:
    S=sorted(set([1,2,3]+[V0+d for d in range(10)]))
    print("     V0=%4d: class=%s cov=%s prim=%s"%(V0,cls(S),is_cov(S),reduce(gcd,S)==1))

print()
print("DONE.")
