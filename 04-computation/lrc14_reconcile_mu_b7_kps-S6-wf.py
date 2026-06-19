#!/usr/bin/env python3
"""
Reconcile two issues:
 (1) engine mu_theta(consec_8)=691/735=0.94014 vs my measpred-midpoint mu=0.94728.
     The engine is authoritative (it integrates the EXACT linear maxgap>theta
     locus per cell). My measpred uses a single midpoint per (order x arc) cell,
     which is WRONG for maxgap because maxgap is piecewise-LINEAR, not constant,
     so the indicator {maxgap>1/7} flips inside a cell. Confirm the engine value
     is the true mu by an independent FINE-grid lower/upper bracketing.
 (2) Given (1), re-evaluate whether B7(E) <= mu_TRUE(E) holds, using the engine's
     mu and a CORRECT B7 (B7's indicator {some fixed arc empty} IS constant on
     fixed-arc-augmented cells, so midpoint B7 is exact). Then state the true B7
     margins to thr_k.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb

# ---- engine (authoritative mu_theta) ----
def mu_theta(E, theta):
    E=sorted(set(E)); n=len(E); bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1); total=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2; order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        ks=[(E[order[t]]*mid).__floor__() for t in range(n)]; subs=[]
        for t in range(n):
            o1=order[t];o2=order[(t+1)%n];k1=ks[t];k2=ks[(t+1)%n];wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1];c=F(k1-k2+wrap)
            if s==0:
                if c>theta: subs.append((a,b))
            elif s>0:
                lo=max(a,(theta-c)/s);  subs.append((lo,b)) if lo<b else None
            else:
                hi=min(b,(theta-c)/s);  subs.append((a,hi)) if a<hi else None
        subs.sort(); cur=cb=None
        for lo,hi in subs:
            if cur is None: cur,cb=lo,hi
            elif lo<=cb: cb=max(cb,hi)
            else: total+=cb-cur; cur,cb=lo,hi
        if cur is not None: total+=cb-cur
    return total

# ---- independent fine-grid bracketing of mu = meas{maxgap>1/7} ----
def maxgap_at(E,x):
    pts=sorted(set((e*x)%1 for e in E))
    if len(pts)==1: return F(1)
    g=[]
    for i in range(len(pts)):
        nx=pts[(i+1)%len(pts)]
        g.append(nx-pts[i] if i+1<len(pts) else (1-pts[i])+pts[0])
    return max(g)

def mu_bracket(E, theta, N):
    """Lower/upper Riemann brackets of meas{maxgap>theta} on a uniform grid of
    N cells. maxgap is continuous & piecewise-linear, so on each subinterval the
    sup and inf of maxgap occur at the endpoints. We sample endpoints; a cell is
    surely-in if both endpoints have maxgap>theta+slack... simpler: count cell as
    'in-lower' if BOTH endpoints>theta (then interior>theta too since pw-linear
    -> min at an endpoint within a linear piece; but breakpoints inside cell may
    create a dip). To be safe use many points per cell for inf/sup."""
    lower=F(0); upper=F(0); h=F(1,N)
    PERCELL=4
    for c in range(N):
        a=c*h; b=(c+1)*h
        vals=[maxgap_at(E, a+ (b-a)*F(j,PERCELL)) for j in range(PERCELL+1)]
        mn=min(vals); mx=max(vals)
        if mn>theta: lower+=h
        if mx>theta: upper+=h
    return lower, upper

def some_fixed_arc_empty(E,x):
    arcs=[(F(2*i+1,14),F(2*i+3,14)) for i in range(7)]
    pts=[(e*x)%1 for e in E]
    for (lo,hi) in arcs:
        has=False
        for p in pts:
            if hi<=1:
                if lo<=p<hi: has=True;break
            else:
                if p>=lo or p<(hi-1): has=True;break
        if not has: return True
    return False

def b7_bps(E):
    E=sorted(set(E)); n=len(E); bp=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,e+1):
            for i in range(7):
                bp.add((F(m)+F(2*i+1,14))/e); bp.add((F(m)+F(2*i+3,14))/e)
    return sorted(b for b in bp if 0<=b<=1)

def B7(E):
    # B7 indicator constant on fixed-arc cells -> midpoint is EXACT
    bp=b7_bps(E); t=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        if some_fixed_arc_empty(E,(a+b)/2): t+=(b-a)
    return t

print("=== (1) Confirm engine mu is true mu via fine grid bracket ===")
for k in (8,9,10):
    E=list(range(k))
    me=mu_theta(E,F(1,7))
    lo,hi=mu_bracket(E,F(1,7),20000)
    inside = lo<=me<=hi
    print(f"consec_{k}: engine mu={me}={float(me):.6f}  grid bracket [{float(lo):.6f},{float(hi):.6f}]  engine in bracket? {inside}")

print()
print("=== (2) B7 (EXACT) vs engine-mu, and B7 margins to thr_k ===")
thr={8:F(3637,5880),9:None,10:None,11:None,12:F(1,7),13:F(0)}
for k in range(8,14):
    E=list(range(k))
    b=B7(E); m=mu_theta(E,F(1,7))
    line=f"consec_{k}: B7={float(b):.5f}  mu_engine={float(m):.5f}  B7<=mu? {b<=m}"
    if thr[k] is not None:
        line+=f"  thr={float(thr[k]):.5f}  B7-thr={float(b-thr[k]):+.5f}"
    print(line)

print()
print("=== (3) Does B7 >= thr hold at the BINDING k=8, and with what margin? ===")
b8=B7(list(range(8)));
print(f"B7(consec_8)={float(b8):.5f}  thr_8={float(thr[8]):.5f}  margin={float(b8-thr[8]):+.5f}")
print("(prompt claimed B7 stays in [0.94,1.0] with margin>=0.30 -- that is the")
print(" IID value 0.9755, NOT the structured deterministic B7=0.736; the real")
print(" minimizer-margin is ~0.118, far smaller than claimed.)")
