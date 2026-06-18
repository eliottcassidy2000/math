#!/usr/bin/env python3
"""
lrc14_k5_exact_closure_macmini.py   (mac-mini, ANGLE E)

GOAL: close LRC(14) case S3 for cluster size k=|L|=5 (the first open size; pigeonhole
margin c=|L|-1=4 gives 7/4-1=3/4<1, so pigeonhole fails -> need rho*).

OBJECT (THM-527):  offset set E, 0 in E, |E|=k.  For x in [0,1) put the k points
{frac(e*x): e in E} on R/Z.  good(x) <=> circular max-gap > 2/7 (fit in open arc < 5/7).
  mu(E)      = meas{x in [0,1):       good(x)}
  rho*(P,E)  = meas{x in G_P:         good(x)},  G_P={tau: ||p tau||>=1/14 for all p in P}.
THM-527: rho*(P,E)>0  ==>  M(P u L) >= 1/14  (case S3).  A uniform floor c0>0 proves S3.

For k=5: |L|=5 large runners, |P|=13-5=8 small runners P subset {1..13}.

EXACT COMPUTATION (Fractions, no float in the decision):
  maxgap{frac(e_i x)} as a function of x is PIECEWISE LINEAR with RATIONAL breakpoints.
  Breakpoints come from (a) a point crossing 0 [frac jumps]: x = n/e_i, and
  (b) two points colliding / swapping order: frac(e_i x)=frac(e_j x) i.e. (e_i-e_j)x in Z
      => x = n/(e_i-e_j); and (c) a gap hitting the threshold 2/7 exactly.
  Collect ALL breakpoints of types (a),(b); on each open subinterval the ORDER of the
  points is fixed and each gap g(x)=A x + B/(something) is linear (rational slope/intercept).
  For each gap we solve g(x)>2/7 exactly (linear ineq) and intersect -> exact good-set.
  Intersect with G_P (exact union of rational arcs) -> rho* exact.

This gives EXACT rho* as a Fraction for every (P, E).  We then:
  - prove the k=5 pigeonhole-type lower bounds where they apply,
  - enumerate bounded-spread 5-shapes and tabulate mu(E) exactly, find the min,
  - compute rho*(P,E) exactly for the worst shapes over all P subset{1..13}, |P|=8,
  - handle the LARGE-SPREAD tail (Angle C: huge spread RAISES mu) rigorously for k=5.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

H = F(1,14)
THR = F(2,7)          # gap threshold: good <=> maxgap > 2/7

# ---------- exact arc / safe-set machinery (G_P) ----------
def danger_arcs(u, h=H):
    """arcs of tau with ||u tau|| < h  (open). returns list of (a,b) with 0<=a<b<=1."""
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-h/u); b=(c+h/u)
        a%=1; b%=1
        if a<b: iv.append((a,b))
        else:   iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def safe_set(A, h=H):
    """G_A = {tau: ||a tau||>=h for all a in A} as union of closed-ish arcs (a,b)."""
    if not A: return [(F(0),F(1))]
    dz=merge([iv for u in A for iv in danger_arcs(u,h)])
    safe=[]; prev=F(0)
    for a,b in dz:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<F(1): safe.append((prev,F(1)))
    return safe
def meas(arcs): return sum((b-a) for a,b in arcs)
def intersect(arcsA, arcsB):
    out=[]; i=j=0
    A=sorted(arcsA); B=sorted(arcsB)
    while i<len(A) and j<len(B):
        a0,a1=A[i]; b0,b1=B[j]
        lo=max(a0,b0); hi=min(a1,b1)
        if lo<hi: out.append((lo,hi))
        if a1<b1: i+=1
        else: j+=1
    return out

# ---------- EXACT good-set of an offset set E ----------
def good_set_exact(E):
    """Return exact union of arcs (a,b) in [0,1) where maxgap{frac(e x)}>2/7.
       E: list/iterable of nonneg ints (0 in E). Duplicate frac-values are collapsed
       (treat as a set of POINTS; only one point at a shared value)."""
    E=sorted(set(E))
    es=[e for e in E]
    # breakpoints in (0,1): n/e (point crosses 0) and n/(ei-ej) (collision)
    bps=set([F(0),F(1)])
    diffs=set()
    for e in es:
        if e>0:
            for n in range(1,e): bps.add(F(n,e))
    for i in range(len(es)):
        for j in range(i+1,len(es)):
            d=abs(es[i]-es[j])
            if d>0:
                diffs.add(d)
                for n in range(1,d): bps.add(F(n,d))
    bp=sorted(bps)
    good=[]
    for t in range(len(bp)-1):
        lo,hi=bp[t],bp[t+1]
        mid=(lo+hi)/2
        # exact positions at midpoint to FIX the cyclic order (order constant on (lo,hi))
        pos=sorted(set((e*mid)%1 for e in es))   # set: collapse coincident points
        # we need, on the open interval, the gap arcs as linear functions.
        # The order is fixed; each point p_k(x)=frac(e_{sigma(k)} x)=e x - floor(e mid).
        # Build, for the fixed cyclic order, the linear form of each consecutive gap.
        # Recover which e produces each sorted position (could be several e collapsing).
        # Represent each occupied "slot" by ONE e (the smallest) -> linear form e x - c.
        slot_e=[]; slot_c=[]
        seen={}
        # map sorted value -> representative e and its floor at mid
        valmap={}
        for e in es:
            v=(e*mid)%1
            fl=(e*mid)-v   # = floor(e*mid)
            if v not in valmap:
                valmap[v]=(e,fl)
        for v in sorted(valmap):
            e,fl=valmap[v]
            slot_e.append(e); slot_c.append(fl)
        m=len(slot_e)
        if m==1:
            # single point: whole circle is one gap of length 1 > 2/7 always
            good.append((lo,hi)); continue
        # consecutive gaps g_k(x)= p_{k+1}(x)-p_k(x), and wrap gap = p_0+1 - p_{m-1}
        # p_k(x) = slot_e[k]*x - slot_c[k]
        # we want x in (lo,hi) with SOME gap > THR.
        # For each gap, linear h(x)=alpha x + beta > THR  => solve.
        cur_good=[]  # subintervals of (lo,hi) that are good
        gaps=[]
        for k in range(m):
            if k<m-1:
                alpha=slot_e[k+1]-slot_e[k]
                beta=-(slot_c[k+1]-slot_c[k])
            else:
                # wrap: p_0+1 - p_{m-1}
                alpha=slot_e[0]-slot_e[m-1]
                beta=F(1)-(slot_c[0]-slot_c[m-1])
            gaps.append((alpha,beta))
        # good iff max_k (alpha_k x+beta_k) > THR. Union over gaps of {x: alpha x+beta>THR}.
        sub=[]
        for alpha,beta in gaps:
            # alpha x + beta > THR  on (lo,hi)
            if alpha==0:
                if beta>THR: sub.append((lo,hi))
                continue
            xb=(THR-beta)/alpha
            if alpha>0:
                a=max(lo,xb)
                if a<hi: sub.append((a,hi))
            else:
                b=min(hi,xb)
                if lo<b: sub.append((lo,b))
        if sub:
            for a,b in merge(sub):
                if a<b: good.append((a,b))
    return merge(good)

def mu_exact(E):
    return meas(good_set_exact(E))

def rho_star_exact(P,E):
    g=good_set_exact(E)
    gp=safe_set(P)
    return meas(intersect(g,gp)), meas(gp)

# =====================================================================
# SANITY: reproduce the known exact pure-cluster mu_k for consecutive E
# =====================================================================
print("="*90)
print("SANITY: pure consecutive-cluster mu_k = meas{x: maxgap{0,x,..,(k-1)x}>2/7}")
print("  (THM-527: mu_3=1, mu_4=19/21, mu_5=9/14, ... mu_13=829/4620)")
print("="*90)
expected={3:F(1),4:F(19,21),5:F(9,14),13:F(829,4620)}
for k in range(3,14):
    E=list(range(k))
    m=mu_exact(E)
    tag=""
    if k in expected:
        tag=" EXPECT "+str(expected[k])+("  OK" if m==expected[k] else "  *** MISMATCH ***")
    print(f"  k={k:2d}: mu_k = {str(m):>14s} = {float(m):.6f}{tag}")

# =====================================================================
# k=5 PART 1: pure mu(E) over ALL 5-shapes E (0 in E, |E|=5), bounded spread
#   spread = max(E). We must show mu(E) bounded away from 0 over ALL shapes,
#   and identify the EXTREMAL (min-mu) shape. Angle C: large spread RAISES mu.
# =====================================================================
print("\n"+"="*90)
print("k=5 PART 1: pure mu(E) over ALL shapes, 0 in E, |E|=5, by spread s=max(E)")
print("  (enumerate every subset {0}u T, T subset {1..s}, |T|=4)  -- find global min")
print("="*90)
glob_min=(F(2),None)
per_spread={}
for s in range(4, 41):                       # spread from 4 (=k-1, consecutive) up
    locmin=(F(2),None)
    for T in itertools.combinations(range(1,s+1),4):
        if max(T)!=s: continue               # require spread exactly s (s in E)
        E=(0,)+T
        m=mu_exact(E)
        if m<locmin[0]: locmin=(m,E)
        if m<glob_min[0]: glob_min=(m,E)
    per_spread[s]=locmin
    print(f"  spread s={s:2d}: min mu = {str(locmin[0]):>14s} = {float(locmin[0]):.6f}  at E={locmin[1]}")
print(f"\n  *** GLOBAL min mu over all 5-shapes (spread<=40) = {str(glob_min[0])} "
      f"= {float(glob_min[0]):.6f}  at E={glob_min[1]} ***")

# Confirm the large-spread monotone-rise (Angle C) for k=5: track min-mu as spread grows
print("\n  Large-spread behaviour (min mu vs spread): is it rising / bounded below?")
xs=sorted(per_spread)
print("    spread:  "+"  ".join(f"{s:>3d}" for s in xs[:18]))
print("    minmu :  "+"  ".join(f"{float(per_spread[s][0]):.2f}" for s in xs[:18]))
