import sys
sys.path.insert(0,'/tmp')
from lrc_base import *

# QUANTIFY the measure surplus vs. coverage gap.
# Total band measure on (0,1/2]: each speed v -> measure 1/7 (on full circle), 
# on (0,1/2] it's... compute exactly the sum of band measures and the actual union.
def band_meas_half(v):
    tot=F(0); k=0
    while F(k,v)-F(1,14*v) < F(1,2):
        lo=max(F(k,v)-F(1,14*v),F(0)); hi=min(F(k,v)+F(1,14*v),F(1,2))
        if lo<hi: tot+=hi-lo
        k+=1
    return tot
def union_meas_half(S):
    pts=[]
    for v in S:
        k=0
        while F(k,v)-F(1,14*v) < F(1,2):
            lo=max(F(k,v)-F(1,14*v),F(0)); hi=min(F(k,v)+F(1,14*v),F(1,2))
            if lo<hi: pts.append((lo,1)); pts.append((hi,-1))
            k+=1
    pts.sort(); cov=F(0); depth=0; prev=F(0)
    for p,d in pts:
        if depth>0: cov+=p-prev
        prev=p; depth+=d
    return cov

S=list(range(1,14))
total_sum=sum(band_meas_half(v) for v in S)
union=union_meas_half(S)
overlap=total_sum-union
print(f"1..13 on (0,1/2]:")
print(f"  sum of band measures = {total_sum} = {float(total_sum):.6f}")
print(f"  union (actual coverage) = {union} = {float(union):.6f} (target 1/2={float(F(1,2)):.4f})")
print(f"  total overlap (redundancy) = {overlap} = {float(overlap):.6f}")
print(f"  surplus over interval = sum - 1/2 = {total_sum-F(1,2)} = {float(total_sum-F(1,2)):.6f}")
print()
print(f"  => union = exactly 1/2 (FULL measure coverage), overlap = {overlap} ({float(overlap):.4f})")
print(f"  The {float(overlap):.4f} of overlap is EXACTLY the redundancy that, if it could be")
print(f"  redistributed without leaving boundary touch-points, would give M<1/14.")
print(f"  But THM-503 forces overlap at every coprime pair (>=1/(7*max)), and the")
print(f"  symmetric pairs {{v,14-v}} create edge-to-edge band meetings -> survivor points.")
