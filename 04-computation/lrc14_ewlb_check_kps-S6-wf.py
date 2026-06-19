"""
Independent check of the EWLB reduction (Angle 1, the 'most promising' route).
EWLB_A(E) = meas( union_{a in A} W_a(E) ), W_a(E)={x: open arc (a,a+1/7) empty of all frac(e x)}.
A = {j/14 : j=0..6}.
Claims to verify:
  (i)  mu_{1/7}(E) >= EWLB_A(E)  [proved inequality; check many E]
  (ii) EWLB_A(consec_k) >= thr_k EXACTLY, binding margin 433/5880 at k=8
  (iii) consecutive minimizes EWLB_A (exhaustive bounded spread)  -- the surviving open step
Also cross-check the B7 minorant correction from verdict 2:
  B_7(consec_8) = 433/588 (~0.736) NOT 0.9755; margin to thr_8 ~0.118.
"""
from fractions import Fraction as F
from itertools import combinations
from math import floor

def mu_theta(E,theta):
    E=sorted(set(E)); n=len(E); bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1); total=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2; order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        ks=[floor(E[order[t]]*mid) for t in range(n)]; subs=[]
        for t in range(n):
            o1=order[t];o2=order[(t+1)%n];k1=ks[t];k2=ks[(t+1)%n];wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1];c=F(k1-k2+wrap)
            if s==0:
                if c>theta: subs.append((a,b))
            elif s>0:
                lo=max(a,(theta-c)/s)
                if lo<b: subs.append((lo,b))
            else:
                hi=min(b,(theta-c)/s)
                if a<hi: subs.append((a,hi))
        subs.sort(); cur=cb=None
        for lo,hi in subs:
            if cur is None: cur,cb=lo,hi
            elif lo<=cb: cb=max(cb,hi)
            else: total+=cb-cur; cur,cb=lo,hi
        if cur is not None: total+=cb-cur
    return total

# measure of x where the open arc (a0, a0+w) contains NO frac(e x), e in E.
# frac(e x) in (a0,a0+w) iff exists integer m with a0 < e x - m < a0+w
#   iff x in ((a0+m)/e, (a0+w+m)/e) for some m (e>0). For e=0: frac=0; in arc iff 0 in (a0,a0+w).
# W_a = complement (over x in [0,1)) of union over e of the "hit" sets.
def empty_arc_measure(E, a0, w):
    # compute measure of x where NO e gives frac(e x) in (a0,a0+w)
    # hit set H = union_e H_e ; answer = 1 - meas(H)
    intervals=[]  # list of (lo,hi) subintervals of [0,1) that are HITS
    for e in E:
        if e==0:
            # frac(0)=0 ; in (a0,a0+w) iff a0<0<a0+w => never for a0>=0. (a0 in [0,1))
            if a0< F(0) < a0+w:
                intervals.append((F(0),F(1)))
            continue
        e=F(e)
        # x in ((a0+m)/e, (a0+w+m)/e) ∩ [0,1)
        # m ranges so interval overlaps [0,1): (a0+m)/e <1 and (a0+w+m)/e>0
        mlo = floor(-(a0+w))  # rough
        mhi = floor(e)+2
        for m in range(int(mlo)-2, int(mhi)+2):
            lo=(a0+m)/e; hi=(a0+w+m)/e
            lo=max(lo,F(0)); hi=min(hi,F(1))
            if lo<hi:
                intervals.append((lo,hi))
    # union measure
    intervals.sort()
    tot=F(0); cur=cb=None
    for lo,hi in intervals:
        if cur is None: cur,cb=lo,hi
        elif lo<=cb: cb=max(cb,hi)
        else: tot+=cb-cur; cur,cb=lo,hi
    if cur is not None: tot+=cb-cur
    return F(1)-tot  # measure of EMPTY (no hit)

def EWLB(E, A, w=F(1,7)):
    # union over a in A of W_a(E). W_a = {x: arc (a,a+w) empty}. Need measure of union.
    # Each W_a is itself a union of intervals; we need the union over a of these.
    # Build per-a empty-sets as interval lists, then union all.
    allint=[]
    for a0 in A:
        # rebuild empty intervals for this a0
        # hits:
        intervals=[]
        for e in E:
            if e==0:
                if a0< F(0) < a0+w: intervals.append((F(0),F(1)))
                continue
            ee=F(e)
            for m in range(-2, int(ee)+3):
                lo=(a0+m)/ee; hi=(a0+w+m)/ee
                lo=max(lo,F(0)); hi=min(hi,F(1))
                if lo<hi: intervals.append((lo,hi))
        intervals.sort(); merged=[]; cur=cb=None
        for lo,hi in intervals:
            if cur is None: cur,cb=lo,hi
            elif lo<=cb: cb=max(cb,hi)
            else: merged.append((cur,cb)); cur,cb=lo,hi
        if cur is not None: merged.append((cur,cb))
        # empty = [0,1) minus merged
        prev=F(0)
        for lo,hi in merged:
            if lo>prev: allint.append((prev,lo))
            prev=max(prev,hi)
        if prev<1: allint.append((prev,F(1)))
    # union of all empties
    allint.sort(); tot=F(0); cur=cb=None
    for lo,hi in allint:
        if cur is None: cur,cb=lo,hi
        elif lo<=cb: cb=max(cb,hi)
        else: tot+=cb-cur; cur,cb=lo,hi
    if cur is not None: tot+=cb-cur
    return tot

A14 = [F(j,14) for j in range(7)]
thr = {8:F(3637,5880),9:F(2025,4004),10:F(36,91),11:F(25,91),12:F(1,7)}
th=F(1,7)

print("=== (ii) EWLB_A(consec_k) vs thr_k ===")
for k in range(8,13):
    E=list(range(k))
    ew=EWLB(E,A14); t=thr[k]; mu=mu_theta(E,th)
    print(f"k={k}: EWLB={ew}={float(ew):.4f}  thr={float(t):.4f}  margin={float(ew-t):.4f}  >=thr:{ew>=t} | mu={float(mu):.4f} mu>=EWLB:{mu>=ew}")
print("binding margin k=8 == 433/5880:", EWLB(list(range(8)),A14)-thr[8]==F(433,5880), EWLB(list(range(8)),A14)-thr[8])

print("\n=== (i) mu >= EWLB on random sets ===")
import random; random.seed(5)
viol=0
for _ in range(800):
    k=random.randint(8,12); pts={0}
    while len(pts)<k: pts.add(random.randint(1,25))
    E=sorted(pts)
    if mu_theta(E,th) < EWLB(E,A14): viol+=1; print("VIOL",E)
print("mu<EWLB violations:",viol,"(0 expected)")

print("\n=== (iii) consecutive minimizes EWLB_A? exhaustive bounded spread ===")
for k in range(8,11):
    S={8:13,9:13,10:14}[k]
    best=None;bestE=None;below=0
    for combo in combinations(range(1,S+1),k-1):
        E=(0,)+combo
        ew=EWLB(list(E),A14)
        if best is None or ew<best: best=ew;bestE=E
        if ew<thr[k]: below+=1
    print(f"k={k} spread<={S}: min EWLB={float(best):.4f} at {bestE} consec-argmin:{bestE==tuple(range(k))} below-thr:{below}")
