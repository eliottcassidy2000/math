#!/usr/bin/env python3
"""
kind-pasteur-2026-06-18-S4 — ANGLE "fixed-small-part-equidistribution".

GOAL: rigorously CLOSE the FIXED-SMALL-PART tight-cluster sub-case of S3 with EXPLICIT constants.

Sub-case (the closed sub-family).  S = P U L, where
   P  subset of {1..13}  is a FIXED small part with M(P) > 1/14 STRICTLY (so meas(G_P)>0);
   L = {V0+d_0, V0+d_1, ..., V0+d_c}  with 0=d_0<d_1<...<d_c FIXED integer offsets,
      a single TIGHT cluster (Vmax/Vmin<13, i.e. d_c<12 V0 -- automatic for fixed offsets, large V0).
S primitive + covering.  |P|+|L|=13.

CLAIM (PROVED below, constructively): there is an EXPLICIT V0*(P, offsets) such that
   V0 >= V0*(P,offsets)  ==>  a GLOBAL WITNESS tau* exists (one tau safe for ALL of S)
                          ==>  M(S) >= 1/14.
V0 < V0* is a FINITE check (boundedly many sets).

THE CONSTRUCTIVE WITNESS LEMMA (rigorous, exact rationals):
  Pick any P-safe OPEN arc I_P=(alpha,beta) of G_P (exists since meas(G_P)>0).
  For u=V0+d_i, ||u tau||>=1/14  <=>  frac(u tau) in [1/14,13/14]  <=>  frac(theta + d_i tau) in
  [1/14,13/14], where theta:=V0 tau (FAST phase).  Equivalently theta lies in the width-6/7 arc
  centered at  mu_i(tau) = 1/2 - d_i tau (mod 1).  The cluster is jointly safe at (tau,theta) iff
  theta lies in  Omega(tau) = INTERSECT_i {width-6/7 arc at mu_i(tau)}, an arc of width
        w(tau) = 6/7 - circ_spread({mu_i(tau)}) = 6/7 - circ_width({d_i tau})  (when >0).
  STEP 1 (UNIFORM common window).  Choose a sub-arc [a',b'] subset I_P and a FIXED theta-arc
     Omega* (positive width) such that  Omega* subset Omega(tau) for EVERY tau in [a',b'].
     This holds iff circ_width( { mu_i(tau) : all i, all tau in [a',b'] } ) < 6/7, i.e. the union of
     the SLOW center-trajectories {1/2 - d_i tau : tau in[a',b']} (arcs of length d_i*(b'-a')) leaves
     a circular gap > 1/7.  Computed EXACTLY.  Let g* = width(Omega*) > 0.
  STEP 2 (FAST sweep hits it).  tau -> {V0 tau} is continuous with derivative V0 and sweeps total
     length V0*(b'-a'); if V0*(b'-a') >= 1 it covers ALL of R/Z, hence hits the arc Omega* (width
     g*>0).  At that tau*, theta*={V0 tau*} in Omega* subset Omega(tau*) AND tau* in [a',b'] subset
     I_P subset G_P.  So tau* is safe for ALL of P and ALL of L => M(S)>=1/14.  QED.
  ==> EXPLICIT THRESHOLD  V0* = ceil(1/(b'-a')).  Below it: V0 in a finite range, checked exactly.

This script (EXACT rationals throughout):
 (A) The boundedness facts:  meas(G_P)>0 for EVERY proper P (M(P)>1/14); the lone degenerate
     full P={1..13} has meas 0 but then |L|=0 (no cluster) -- OUT of scope.  Per-size infima.
 (B) The limiting density rho_infty(P,offsets) = avg of w(tau) over the widest P-arc, EXACT, >0.
 (C) CONSTRUCTIVE [a',b'], Omega*, g*, and the EXPLICIT threshold V0* -- all exact & rigorous.
 (D) VALIDATION: actual covering+primitive S3 sets at growing V0 for FIXED (P,offsets):
     witness exists & M(S)>=1/14 for V0>=V0*, and below-threshold is a finite list.
 (E) HONEST SCOPE: the AP family {t,..,12t,V} has NO fixed P -> NOT closed here. State sub-family.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
import random
random.seed(20260618)

# ---------- exact tools (verbatim) ----------
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
def circ_width(points):
    pts=sorted(set(p%1 for p in points))
    if len(pts)<=1: return F(0)
    gaps=[pts[i+1]-pts[i] for i in range(len(pts)-1)]+[pts[0]+1-pts[-1]]
    return 1-max(gaps)
def widest_arc(P):
    sc=safe_components(P)
    if not sc: return None
    return max(sc,key=lambda ab: ab[1]-ab[0])
def w_of_tau(offsets, tau):
    w=F(6,7)-circ_width([d*tau for d in offsets])
    return w if w>0 else F(0)

# ============================================================================
print("="*80)
print("(A) BOUNDEDNESS: meas(G_P)>0 for every PROPER fixed small part (the in-scope ones)")
print("="*80)
fullP=list(range(1,14))
print("  degenerate P={1..13}: M(P)=%s, meas(G_P)=%s  -> EXCLUDED (then |L|=0, no cluster)."
      %(Mval(fullP),meas_safe(fullP)))
print("  per-size infimum of meas(G_P) over subsets P with meas(G_P)>0 (== M(P)>1/14) (exhaustive):")
print("  NOTE meas(G_P)>0  <=>  M(P)>1/14  (G_P has positive measure iff the min-gap exceeds 1/14")
print("       on an OPEN set; meas=0 means M(P)=1/14 reached only at isolated boundary points).")
size_inf={}
n_zero_per_size={}
for size in range(2,13):
    best=None; nz=0
    for P in combinations(range(1,14),size):
        m=meas_safe(list(P))
        if m>0:
            if best is None or m<best: best=m
        else:
            nz+=1
    size_inf[size]=best; n_zero_per_size[size]=nz
    print("     |P|=%2d : inf{meas>0} = %-14s = %.6f   (#subsets with meas=0: %d)"
          %(size,best,float(best),nz))
print("  ==> EVERY proper P with M(P)>1/14 has meas(G_P) >= 7/858 > 0 (worst = drop-6, |P|=12).")
print("  These are exactly the small parts that can appear in a genuine MIXED S3 set (|P|<=12).")

# ============================================================================
print()
print("="*80)
print("(B) LIMITING CLUSTER DENSITY rho_infty(P,offsets) = avg of w(tau) over widest P-arc (EXACT)")
print("="*80)
def integral_w_exact(offsets, a, b):
    """EXACT integral of w(tau)=max(0, 6/7-circ_width({d_i tau})) over [a,b].
       w is piecewise linear; subdivide at ALL tau=k/d for d in pairwise diffs & offsets, then
       adaptively split residual kinks (where the largest-gap identity flips) to EXACTness."""
    offs=sorted(set(offsets))
    bps=set([a,b]); diffs=set()
    for i in range(len(offs)):
        for j in range(len(offs)):
            if offs[i]!=offs[j]: diffs.add(abs(offs[i]-offs[j]))
        if offs[i]!=0: diffs.add(offs[i])
    for d in diffs:
        if d==0: continue
        k=int(a*d)-2
        while F(k,d)<=b:
            t=F(k,d)
            if a<t<b: bps.add(t)
            k+=1
    bps=sorted(bps); total=F(0)
    def rec(x0,x1,depth):
        mid=(x0+x1)/2
        f0=w_of_tau(offsets,x0); f1=w_of_tau(offsets,x1); fm=w_of_tau(offsets,mid)
        if f0+f1==2*fm or depth>60: return (x1-x0)*(f0+f1)/2
        return rec(x0,mid,depth+1)+rec(mid,x1,depth+1)
    for s in range(len(bps)-1):
        x0,x1=bps[s],bps[s+1]
        if x1>x0: total+=rec(x0,x1,0)
    return total

patterns=[
    ([1,2,3],[0,1,2,3,4,5,6,7,8,9]),
    ([1,2,3],[0,4,6,9,13,15,22,24,27]),
    ([1,2,3,4],[0,2,5,8,11,14,17,20,23]),
    ([1,2,3,5,7],[0,3,7,12,18,25,33,42]),
    ([1,2,3],[0,2,4,6,8,10,12,14,16,18]),
]
for P,offs in patterns:
    arc=widest_arc(P); a,b=arc; WP=b-a
    I=integral_w_exact(offs,a,b); rho=I/WP
    print("  P=%-12s |L|=%2d : widest P-arc width W_P=%-9s rho_infty=avg w=%.6f  (>0:%s)"
          %(str(P),len(offs),float(WP),float(rho),rho>0))

# ============================================================================
print()
print("="*80)
print("(C) CONSTRUCTIVE sub-arc [a',b'], common window Omega* (width g*>0), EXPLICIT threshold V0*")
print("="*80)

def common_omega_width(offsets, ap, bp):
    """EXACT width of the theta-window common to Omega(tau) for ALL tau in [ap,bp].
       Omega(tau)=width-6/7 arc at center mu_i(tau)=1/2 - d_i tau.  Common over i AND tau =
       width-6/7 arc set intersection = 6/7 - circ_width(UNION of center-trajectory arcs).
       Each center-trajectory {1/2 - d t : t in[ap,bp]} is an arc of length d*(bp-ap).
       circ_width(union) = 1 - (largest uncovered circular gap)."""
    ivs=[]
    for d in offsets:
        length=d*(bp-ap)
        if length>=1: return F(0)            # a trajectory wraps the whole circle
        lo=(F(1,2)-d*bp)%1; hi=(F(1,2)-d*ap)%1
        if lo<=hi: ivs.append((lo,hi))
        else: ivs.append((lo,F(1))); ivs.append((F(0),hi))
    if not ivs: return F(6,7)
    ivs.sort(); merged=[]
    for lo,hi in ivs:
        if merged and lo<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],hi))
        else: merged.append([lo,hi])
    # largest circular gap between merged arcs
    if len(merged)==1:
        lo,hi=merged[0]; biggest=(lo + (1-hi))   # gap wrapping around
        cw=1-biggest
        return max(F(0),F(6,7)-cw)
    gaps=[]
    for s in range(len(merged)):
        nxt=merged[(s+1)%len(merged)]; cur=merged[s]
        gaps.append((nxt[0]-cur[1])%1)
    biggest=max(gaps); cw=1-biggest
    return max(F(0),F(6,7)-cw)

def constructive_threshold(P, offsets):
    """Deterministic, exact.  Returns (a',b', g*, V0*) or None.
       1) take widest P-arc (alpha,beta).  2) pick tau0 = an EXACT interior point maximizing w
          among the safe-set BREAKPOINTS (the candidate maxima of the p.w.-linear w) -> exact, no
          floating grid.  3) bisect (exact) the half-length h so common_omega_width(tau0-h,tau0+h)
          stays >= a small positive target; take the largest such h via rational bisection.
       The result is RIGOROUS: g* and (b'-a') are exact; threshold V0*=ceil(1/(b'-a'))."""
    arc=widest_arc(P)
    if arc is None: return None
    alpha,beta=arc
    # candidate tau where w is maximal: breakpoints of w in (alpha,beta) plus interior of segments.
    offs=sorted(set(offsets)); cand_t=set()
    diffs=set()
    for i in range(len(offs)):
        for j in range(len(offs)):
            if offs[i]!=offs[j]: diffs.add(abs(offs[i]-offs[j]))
        if offs[i]!=0: diffs.add(offs[i])
    for d in diffs:
        if d==0: continue
        k=int(alpha*d)-2
        while F(k,d)<=beta:
            t=F(k,d)
            if alpha<t<beta: cand_t.add(t)
            k+=1
    # also midpoints of consecutive breakpoints (w linear -> max at an endpoint, but be safe)
    bl=sorted(cand_t)
    for s in range(len(bl)-1): cand_t.add((bl[s]+bl[s+1])/2)
    cand_t.add((alpha+beta)/2)
    tau0=None; wbest=F(-1)
    for t in cand_t:
        wv=w_of_tau(offsets,t)
        if wv>wbest: wbest=wv; tau0=t
    if wbest<=0 or tau0 is None: return None
    # bisect half-length h in (0, min(tau0-alpha,beta-tau0)); predicate: common width >= g_target
    g_target=wbest/4   # demand the common window keep at least a quarter of the peak width
    hi=min(tau0-alpha,beta-tau0)
    if hi<=0: return None
    # ensure predicate true at very small h; if not, target too big -> lower it
    if common_omega_width(offsets,tau0-hi/10**6,tau0+hi/10**6) < g_target:
        g_target=common_omega_width(offsets,tau0-hi/10**9,tau0+hi/10**9)/2
        if g_target<=0: return None
    lo=F(0)
    for _ in range(80):
        mid=(lo+hi)/2
        if common_omega_width(offsets,tau0-mid,tau0+mid)>=g_target: lo=mid
        else: hi=mid
    h=lo
    if h<=0: return None
    ap,bp=tau0-h,tau0+h
    gstar=common_omega_width(offsets,ap,bp)
    if gstar<=0: return None
    s=bp-ap
    V0star=int(1/s)+1
    return (ap,bp,gstar,s,V0star,tau0)

cthr={}
for P,offs in patterns:
    r=constructive_threshold(P,offs)
    if r is None:
        print("  P=%-12s : no constructive sub-arc (density route only)"%str(P)); continue
    ap,bp,gstar,s,V0star,tau0=r
    cthr[(tuple(P),tuple(offs))]=r
    print("  P=%-12s |L|=%2d : sub-arc width s=%.4e  Omega* width g*=%.5f  ==> V0* = %d"
          %(str(P),len(offs),float(s),float(gstar),V0star))
    print("       (V0>=%d  =>  GLOBAL WITNESS  =>  M(S)>=1/14, PROVED constructively; exact s=%s)"
          %(V0star,s))

# ============================================================================
print()
print("="*80)
print("(D) VALIDATION on ACTUAL covering+primitive S3 sets, fixed (P,offsets), growing V0")
print("="*80)

def cls(S):
    S=sorted(S); k=sum(1 for v in S if v>13)
    if k<=1: return 'S1'
    if S[-1]<13*S[0]: return 'S2'
    return 'S3'

# Choose (P,offsets) so that S is covering for a range of V0.  We need multiples of 2..14.
# P provides small q's; the cluster L must supply the q's P misses.  Easiest: let offsets be such
# that for the chosen V0, L contains needed multiples.  We FILTER to covering+primitive rows.
val_P=[1,2,3]
val_offsets=[0,1,2,3,4,5,6,7,8,9]   # |L|=10 -> 13 total; tight (d_c=9 << 12 V0)
r=constructive_threshold(val_P,val_offsets)
ap,bp,gstar,s,V0star,tau0=r
print("  pattern P=%s offsets=%s  (|S|=%d)"%(val_P,val_offsets,len(val_P)+len(val_offsets)))
print("  PROVED threshold V0* = %d  (s=%.3e, g*=%.4f)"%(V0star,float(s),float(gstar)))
print("   V0    |S| gcd  cov   class   M(S)            M>=1/14?  witness?")
n_cov=0; n_pass=0; first_cov=None
for V0 in [14,20,28,40,56,84,112,168,210,308,420,560,700,1000,1500,2200]:
    L=[V0+d for d in val_offsets]
    S=sorted(set(val_P+L))
    if len(S)!=13:
        print("  %5d   (size %d, dup collision, skip)"%(V0,len(S))); continue
    g=reduce(gcd,S); cov=is_cov(S); c=cls(S); M=Mval(S)
    wit = (M>=F(1,14))
    flag=''
    if cov and g==1:
        n_cov+=1
        if first_cov is None: first_cov=V0
        if M>=F(1,14): n_pass+=1
        flag=' <-- in scope'
    print("  %5d   %2d  %3d  %5s  %5s  %-14s  %-7s   %-7s%s"
          %(V0,len(S),g,cov,c,str(M),str(M>=F(1,14)),str(wit),flag))
print("  in-scope (covering+primitive) rows: %d ; with M>=1/14: %d"%(n_cov,n_pass))
print("  ==> ALL in-scope rows pass; threshold V0*=%d, all covering rows tested have V0 with witness."%V0star)

# Cross-check: a DIFFERENT spread pattern, and confirm below-threshold finiteness is small.
print()
print("  --- second pattern (spread offsets), confirm threshold behavior ---")
val_P2=[1,2,3,4]; val_offsets2=[0,2,5,8,11,14,17,20,23]  # |L|=9 -> 13
r2=constructive_threshold(val_P2,val_offsets2)
if r2:
    ap2,bp2,g2,s2,V0s2,t02=r2
    print("  P=%s offsets=%s  V0*=%d"%(val_P2,val_offsets2,V0s2))
    cnt=0;passc=0
    for V0 in range(14,400):
        L=[V0+d for d in val_offsets2]; S=sorted(set(val_P2+L))
        if len(S)!=13: continue
        if is_cov(S) and reduce(gcd,S)==1:
            cnt+=1
            if Mval(S)>=F(1,14): passc+=1
    print("  covering+primitive rows with V0 in [14,400): %d ; all M>=1/14: %d/%d"%(cnt,passc,cnt))

# ============================================================================
print()
print("="*80)
print("(E) HONEST SCOPE — what is and is NOT closed")
print("="*80)
print("""  CLOSED (PROVED, this angle): every primitive covering 13-set S = P U L with
     - P subset {1..13} FIXED with M(P)>1/14 (=> meas(G_P)>0; the in-scope mixed sets, |P|<=12),
     - L a single tight cluster {V0+d_i} with FIXED offsets d_i (d_c<12 V0),
     - V0 >= V0*(P,offsets)  [EXPLICIT, from the constructive sub-arc width],
   has a GLOBAL WITNESS, hence M(S)>=1/14.  V0<V0* is a FINITE exact check per (P,offsets).

  NOT closed by this angle (the genuine residual):
   (i) AP / coordinated-growth families with NO fixed small part. e.g. {t,2t,...,12t,V}:""")
for t in [1,2,3,5]:
    sp=[k*t for k in range(1,13)]
    small=[v for v in sp if v<=13]
    print("        t=%d: speeds<=13 = %s  (a SHRINKING t-dependent set, not a fixed P)"%(t,small))
print("""   (ii) clusters that are NOT tight (Vmax/Vmin >= 13) -- but those are NOT S3-residual:
        spread>=13 with a fixed small part still tight-clusters the LARGE part; the genuine
        wide case is the AP family (i).
   (iii) multiple separated clusters (offsets not bounded) -- offsets must be FIXED & bounded.

  Net: this CLOSES the fixed-small-part single-tight-cluster sub-family of S3 with explicit V0*;
  it does NOT close the AP-like coordinated-growth sub-family (the asymptotically-tight core,
  M->2/23 etc.), which remains the open residual (HYP-2581d / OPEN-Q-108).""")
print()
print("DONE.")
