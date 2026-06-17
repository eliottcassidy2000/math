#!/usr/bin/env python3
"""
lrc14_lowerbound_proof_kps  (kind-pasteur, PROVE side)

Toward a PROVABLE lower bound on min positive single-perturbation L.

STRUCTURE of the argument for drop e, add w (w not tight-preserving):
  L = meas(G_e \ D_w) = meas(G_e) - meas(G_e ∩ D_w).
  D_w is a union of w+1 arcs each of half-width 1/(14w), total measure 1/7.
  Within ANY arc A of G_e of length |A|, w can cover at most:
     - if a center k/w lands in A: it covers an interval of length min(|A|, 2/(14w))
       INTERSECTED with A, plus possibly partial arcs from adjacent centers.
  The residual in arc A is >= |A| - (number of w-arcs meeting A)*2/(14w).

KEY LOWER BOUND (per-arc): if arc A has length |A| and w has spacing 1/w between
centers, the number of w-centers whose arc can intersect A is at most
  floor(|A| * w) + 1, each contributing <= 2/(14w) coverage, BUT consecutive
w-arcs leave GAPS of length 1/w - 2/(14w) = (1 - 2/14)/w = (6/7)/w between them.
So if |A| > 1/w (more than one spacing), at least one full gap (6/7)/w of A is
uncovered UNLESS a smaller speed already covers it (but in G_e nothing does).
  => residual(A) >= |A| - ceil(|A|*w)*2/(14w) ... and crucially if w is large the
     comb is fine and residual -> |A|*(6/7).

THE MINIMUM is a tradeoff: small w (coarse comb) can ALIGN a center exactly at the
arc center and cover it fully (this is how w=24 achieves 0), but then OTHER arcs of
G_e may be missed. The optimum w=36 aligns all 4 centers (multiple of 12) but its
arcs (half-width 1/504) are too NARROW for the central arcs (half-length 1/490).

We make the central-arc residual exact and general:
  central arc has half-length L_c = 1/490 (i.e. arc [5/12 - 1/490 + ..]) hmm,
  let's just measure: central arc length = 1/245, so half-length = 1/490.
  A w-center at 5/12 covers half-width 1/(14w). For FULL coverage of one side need
  1/(14w) >= 1/490  <=>  w <= 35.  So w=35 would cover, w=36 does NOT (1/504<1/490).
  This is the SHARP threshold: w<=35 multiples-of-12 don't exist except 24; w=36 is
  the FIRST multiple of 12 that fails to cover the central arc, by the tiniest margin.

Verify this threshold logic exactly.
"""
from fractions import Fraction as F

def danger(v):
    out=[]; w=F(1,14*v)
    for k in range(v+1):
        lo=F(k,v)-w; hi=F(k,v)+w
        if lo<0: out += [(F(0),hi),(1+lo,F(1))]
        elif hi>1: out += [(lo,F(1)),(F(0),hi-1)]
        else: out.append((lo,hi))
    return [(x,y) for x,y in out if y>x]
def union(arcs):
    arcs=sorted((x,y) for x,y in arcs if y>x)
    res=[]; cl,ch=arcs[0]
    for lo,hi in arcs[1:]:
        if lo<=ch: ch=max(ch,hi)
        else: res.append((cl,ch)); cl,ch=lo,hi
    res.append((cl,ch)); return res
def total(a): return sum(y-x for x,y in union(a))
def complement(arcs):
    u=union(arcs); res=[]; prev=F(0)
    for lo,hi in u:
        if lo>prev: res.append((prev,lo))
        prev=max(prev,hi)
    if prev<1: res.append((prev,F(1)))
    return res

base=list(range(1,14))
G12=complement([a for v in base if v!=12 for a in danger(v)])
# central arc
ca=[a for a in G12 if F(2,5)<a[0]<F(1,2)][0]
print(f"central G_12 arc = [{ca[0]},{ca[1]}], length = {ca[1]-ca[0]} = {float(ca[1]-ca[0]):.6e}")
print(f"half-length = {(ca[1]-ca[0])/2} = {float((ca[1]-ca[0])/2):.6e}")
print()
print("A speed-w center placed at 5/12 covers half-width 1/(14w).")
print("To cover the central arc FULLY from its center, need 1/(14w) >= half-length.")
hl=(ca[1]-ca[0])/2
print(f"  need 1/(14w) >= {hl}  <=>  w <= 1/(14*{hl}) = {1/(14*hl)} = {float(1/(14*hl)):.4f}")
print(f"  so w <= 35 covers, w >= 36 does NOT (if center exactly at 5/12).")
print()
# but the arc is not symmetric about 5/12! check
print(f"central arc midpoint = {(ca[0]+ca[1])/2} = {float((ca[0]+ca[1])/2):.6f}, vs 5/12={float(F(5,12)):.6f}")
print(f"  arc is {'' if (ca[0]+ca[1])/2==F(5,12) else 'NOT '}centered at 5/12")
print(f"  left half (5/12 - lo) = {F(5,12)-ca[0]}, right half (hi - 5/12) = {ca[1]-F(5,12)}")
print()
# residual when w=36 center at 5/12 covers [5/12-1/504, 5/12+1/504]:
cov_lo=F(5,12)-F(1,504); cov_hi=F(5,12)+F(1,504)
res_left=max(F(0),cov_lo-ca[0])
res_right=max(F(0),ca[1]-cov_hi)
print(f"w=36 covers [{cov_lo},{cov_hi}] within arc.")
print(f"  residual left  [{ca[0]},{cov_lo}] = {res_left}")
print(f"  residual right [{cov_hi},{ca[1]}] = {res_right}")
print(f"  total residual this arc = {res_left+res_right} (matches 1/2520? {res_left+res_right==F(1,2520)})")
print()
print("WHY only left residual: right end of arc 41/98 = 13's left boundary; 5/12+1/504")
print(f"  41/98={float(F(41,98)):.6f}, 5/12+1/504={float(cov_hi):.6f}: cov_hi>=41/98? {cov_hi>=ca[1]}")
print(f"  so 36 covers the whole RIGHT side (up to 13's arc) but misses LEFT (5's side).")
print()

# GENERAL: prove min positive L >= 1/1260 over family {1..11,13,m} by the threshold.
# The argument: any m with L<1/1260 would need to cover both central arcs to within
# total residual <1/1260, i.e. each central residual <1/2520. Show no integer m does.
print("Scan: for which m is the SUM of central-arc residuals < 1/1260 but >0?")
left_arc=ca; right_arc=[a for a in G12 if F(1,2)<a[0]<F(3,5)][0]
import math
hits=[]
for m in range(14,2000):
    if m==12: continue
    Dm=danger(m)
    # just measure residual in the two central arcs
    def resid(arc):
        # arc minus Dm
        u=union(Dm)
        # complement within arc
        pts=[arc[0]];
        segs=[]; prev=arc[0]
        for lo,hi in u:
            lo=max(lo,arc[0]); hi=min(hi,arc[1])
            if lo<hi:
                if lo>prev: segs.append((prev,lo))
                prev=max(prev,hi)
        if prev<arc[1]: segs.append((prev,arc[1]))
        return sum(b-a for a,b in segs)
    r=resid(left_arc)+resid(right_arc)
    if 0<r<F(1,1260): hits.append((r,m))
print(f"  m with central residual in (0,1/1260): {hits[:10]}  (count {len(hits)})")
print("  (these would need OUTER arcs also covered to actually give total L<1/1260)")
