# Measure-horn with the SHARP threshold, continuous circle.
# danger(v) = union of open arcs width 1/(7v) centered at j/v.  safe gap width 6/(7v).
# Claim (sharp): a connected safe component of length L > 1/(7k) contains a k-safe point.
# So if G=(core + 5 killers safe set) has max component L, the 6th killer k6 fails to cover
# whenever k6 > 1/(7L).  Threshold T_sharp = 1/(7L).  (kind-pasteur used 1/(3L), conservative.)
from fractions import Fraction as F
import itertools, random

def danger_arcs(v):  # list of (lo,hi) open arcs on [0,1), width 1/(7v) around j/v
    w = F(1,14*v)
    return [(F(j,v)-w, F(j,v)+w) for j in range(v)]

def subtract_arcs(safe, arcs):
    # safe: sorted disjoint list of (lo,hi) in [0,1). Subtract each danger arc (may wrap).
    # Normalize danger arcs into [0,1), splitting wrap-around.
    cuts=[]
    for (lo,hi) in arcs:
        lo%=1; hi%=1
        if lo<hi: cuts.append((lo,hi))
        else: cuts.append((lo,F(1))); cuts.append((F(0),hi))
    cuts.sort()
    for (clo,chi) in cuts:
        new=[]
        for (lo,hi) in safe:
            if chi<=lo or clo>=hi: new.append((lo,hi)); continue
            if clo>lo: new.append((lo,clo))
            if chi<hi: new.append((chi,hi))
        safe=new
    return safe

def safe_set(speeds):
    safe=[(F(0),F(1))]
    for v in speeds:
        safe=subtract_arcs(safe, danger_arcs(v))
    return safe

def max_component(safe):
    if not safe: return F(0)
    return max(hi-lo for (lo,hi) in safe)

# sanity: single core {1..12}, deep well killer 182 -> should be lonely (M=14/183)
P=list(range(1,13))
S=safe_set(P)
print("core {1..12}: max safe component =", float(max_component(S)), " #arcs=", len(S))
G=subtract_arcs(S, danger_arcs(182))
print("  after killer 182: max comp =", float(max_component(G)), " sharp T=1/(7L)=", float(1/(7*max_component(G))) if max_component(G)>0 else "cov")

# ---- reproduce kind-pasteur worst finite window: core [1,2,4,7,9,11,12], killers ~158..166 ----
def horn_T_sharp(P, removed):
    S=safe_set(P)
    G=subtract_arcs(S, [a for k in removed for a in danger_arcs(k)])
    L=max_component(G)
    if L==0: return None, F(0)   # fully covered by the 5 removed alone -> family lonely iff k6 also... actually covered
    return float(1/(7*L)), L

# r=5 worst reported: core [1,2,4,7,9,11,12], removed 5 consecutive killers near 158
for base in [158,160,162]:
    removed=[base+2*i for i in range(5)]  # 5 killers
    P=[1,2,4,7,9,11,12]
    Ts,L=horn_T_sharp(P,removed)
    print(f"core {P} removed {removed}: L={float(L):.6g} T_sharp(1/7L)={Ts:.2f}  T_cons(1/3L)={float(1/(3*L)):.2f}")
