#!/usr/bin/env python3
"""
Resolve the B_7 discrepancy. The prompt claims:
  - B_7(consec_8) stays in [0.94,1.0], "far above thr_8=0.6185";
  - B_7 iid floor = 1 - sum_{j=0}^7 (-1)^j C(7,j)(1-j/7)^k giving 0.9755 at k=8;
  - margin to thr_k >= 0.30.

My CORRECTED exact B7 (some fixed width-1/7 arc empty) gives only 0.736 at
consec_8. Two different objects are being conflated. Investigate which "B7"
the iid formula actually describes, and what the TRUE exact B7 floor is, and
whether B7 >= thr_k actually holds with the claimed margin.

iid model: k points placed independently uniform on the circle. 7 fixed arcs
each width 1/7 tiling the circle. P(a given arc empty) = (1-1/7)^k = (6/7)^k.
By inclusion-exclusion, P(at least one of 7 arcs empty) =
  sum_{j=1}^{7} (-1)^{j-1} C(7,j) P(j given arcs all empty).
P(j given arcs empty) = ((7-j)/7)^k = (1 - j/7)^k.
So P(union) = sum_{j=1}^7 (-1)^{j-1} C(7,j)(1-j/7)^k
            = -(sum_{j=1}^7 (-1)^j C(7,j)(1-j/7)^k)
            = (1 - 0) - sum_{j=0}^7 (-1)^j C(7,j)(1-j/7)^k   [adding/subtracting j=0 term =1]
Wait: sum_{j=0}^7 (-1)^j C(7,j)(1-j/7)^k has j=0 term = 1.
So 1 - sum_{j=0}^7 (...) = 1 - [1 + sum_{j=1}^7(-1)^j C(7,j)(1-j/7)^k]
                         = - sum_{j=1}^7 (-1)^j C(7,j)(1-j/7)^k
                         = sum_{j=1}^7 (-1)^{j-1} C(7,j)(1-j/7)^k = P(union). OK.
So the iid formula = P(at least one of 7 fixed disjoint width-1/7 arcs empty)
under k IID uniform points. That is EXACTLY the iid analogue of my B7. So both
should be the SAME object. Then why does iid give 0.9755 at k=8 but the
DETERMINISTIC consec_8 gives only 0.736??

Resolution: consec_8 = {0,..,7}. The points are frac(0*x),...,frac(7*x) =
{0, x, 2x, ..., 7x} mod 1. For x near rational p/q these are FAR from iid --
they are an arithmetic progression, which COVERS arcs very evenly (low
discrepancy), so FEWER arcs are empty than iid. So the deterministic B7 for the
EVEN ruler consec_8 is LOWER than the iid expectation, not equal. The iid value
is NOT a lower bound on the deterministic B7 for structured E. So the prompt's
claim "B_7 iid floor 0.9755 >= thr_8 with margin" does NOT bound the actual
B7(consec_8)=0.736. The relevant quantity is min_E B7(E), and for the AP/consec
minimizer it is ~0.736 (k=8), with margin to thr_8 only ~0.118.
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb

def maxgap_at(E, x):
    pts = sorted(set((e*x) % 1 for e in E))
    if len(pts) == 1: return F(1)
    g = []
    for i in range(len(pts)):
        nx = pts[(i+1)%len(pts)]
        g.append(nx-pts[i] if i+1<len(pts) else (1-pts[i])+pts[0])
    return max(g)

def some_fixed_arc_empty(E, x):
    arcs = [(F(2*i+1,14), F(2*i+3,14)) for i in range(7)]
    pts = [(e*x)%1 for e in E]
    for (lo,hi) in arcs:
        has=False
        for p in pts:
            if hi<=1:
                if lo<=p<hi: has=True;break
            else:
                if p>=lo or p<(hi-1): has=True;break
        if not has: return True
    return False

def bps(E):
    E=sorted(set(E)); n=len(E); bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    for e in E:
        if e==0: continue
        for m in range(0,e+1):
            for i in range(7):
                bp.add((F(m)+F(2*i+1,14))/e); bp.add((F(m)+F(2*i+3,14))/e)
    return sorted(b for b in bp if 0<=b<=1)

def measpred(E,pred):
    bp=bps(E); t=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        if pred(E,(a+b)/2): t+=(b-a)
    return t

def B7(E): return measpred(E, some_fixed_arc_empty)
def mu(E):  return measpred(E, lambda E,x: maxgap_at(E,x)>F(1,7))
def B7iid(k):
    s=F(0)
    for j in range(8): s+=(-1)**j*comb(7,j)*F(7-j,7)**k
    return 1-s

thr={8:F(3637,5880),9:None,10:None,11:None,12:F(1,7),13:F(0)}

print("=== TRUE exact B7 minorant vs iid 'floor' vs thr_k ===")
for k in range(8,14):
    E=list(range(k))
    b=B7(E); m=mu(E); bi=B7iid(k)
    line=f"k={k}: B7(consec)={float(b):.5f}  B7_iid={float(bi):.5f}  mu(consec)={float(m):.5f}"
    if thr[k] is not None:
        line+=f"  thr={float(thr[k]):.5f}  B7-thr={float(b-thr[k]):+.5f}"
    print(line)
    print(f"        B7(consec) <= mu(consec)? {b<=m}   (minorant valid)")
    print(f"        iid >= true B7? {bi>=b}  (iid OVER-estimates structured B7 => NOT a lower bound on min_E B7)")

print()
print("=== min_E B7 over k=8 small spread (is consec/AP the B7 MINIMIZER?) ===")
# Hunt for E minimizing B7 at k=8 (the relevant quantity for the route).
best=(tuple(range(8)),B7(list(range(8))))
cnt=0
for MAX in range(7,13):
    for combo in combinations(range(1,MAX),6):
        E=(0,)+combo+(MAX,); cnt+=1
        b=B7(list(E))
        if b<best[1]: best=(E,b)
print(f"tested {cnt} k=8 sets (spread<=12). min B7 = {float(best[1]):.5f} at {best[0]}")
print(f"consec_8 B7 = {float(B7(list(range(8)))):.5f}; thr_8 = {float(thr[8]):.5f}")
print(f"min B7 >= thr_8? {best[1] >= thr[8]}  margin {float(best[1]-thr[8]):+.5f}")
