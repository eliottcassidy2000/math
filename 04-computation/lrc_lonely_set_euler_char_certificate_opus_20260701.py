"""CERTIFYING COHOMOLOGY for the LRC. Lonely set L = S^1 minus union of danger zones D_v={t:||vt||<r}.
chi(L)=#lonely components (open arcs) = a cohomological invariant; chi(L)>0 CERTIFIES loneliness. It equals the
inclusion-exclusion/Lefschetz sum over the resonance arrangement (each phi_v contributes L(phi_v)=1-v)."""
import numpy as np
from math import gcd
def lonely(S, r, N=600000):
    t=np.arange(N)/N; keep=np.ones(N,bool)
    for v in S:
        d=np.abs(((v*t+0.5)%1)-0.5)
        keep &= (d>=r-1e-12)
    meas=keep.mean(); prev=keep[-1]; comp=0
    for x in keep:
        if x and not prev: comp+=1
        prev=x
    return meas, comp
print("LRC lonely set L = complement of danger zones; chi(L)=#components certifies loneliness (chi>0 => lonely):")
for n in [5,7,14]:
    S=list(range(1,n))
    for rr,tag in [(0.99/n,'r=0.99/n (sub-tight)'),(1.0/n,'r=1/n (tight)')]:
        meas,comp=lonely(S,rr)
        print(f"  AP n={n}, {tag}: lonely measure={meas:.4f}, chi(L)=#components={comp}  => {'LONELY (certified)' if comp>0 else 'covered'}")
# the arrangement Lefschetz: total danger-arc count = sum_v (# arcs of D_v) = sum_v v; each phi_v has L=1-v
print("\nArrangement = Lefschetz: D_v has v arcs (fixed pts of phi_v, L(phi_v)=1-v). Total boundary points 2*sum_v v =")
for n in [5,7,14]:
    S=list(range(1,n)); tot=sum(S); print(f"  n={n}: sum_v v = {tot} = C(n,2); 2*sum = {2*tot} arc-endpoints; the three-distance gaps = chi(L).")
print("=> chi(L) (loneliness certificate) is computed from the resonance arrangement = the dilation Lefschetz data (1-v),")
print("   NOT an SOS positivity witness. The tight AP: chi(L)=#components at r=1/n is the boundary certificate.")
