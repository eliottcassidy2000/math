"""
Creative observer reframes (owner): DILATE (solved covmin), TRANSLATE, COPY-TO-ALL-N, HAMILTONIAN PATHS.
Compute the LRC 'ordering' family: the n points {0, s_i t} as t varies trace a 1-parameter family of
circular ORDERINGS (Ham paths). Connect to tournaments (each ordering = a linear order = a transitive 'tour').
"""
import math
from fractions import Fraction
from itertools import permutations
def orderings(speeds, n, samples=20000):
    pts0=[0]+list(speeds)  # n points (observer 0 + n-1 runners)
    seen=set()
    for i in range(1,samples):
        t=i/samples
        pos=sorted(range(n), key=lambda j:(pts0[j]*t)%1)
        # circular ordering up to rotation: fix observer (index 0) first
        z=pos.index(0); rot=tuple(pos[z:]+pos[:z])
        seen.add(rot)
    return len(seen)
print("(4) HAM-PATH/ORDERING reframe: # distinct circular orderings of {0,s_i t} as t varies (the observer's")
print("    1-parameter family of Hamiltonian 'tours' of the n points):")
for n in [4,5,6]:
    AP=list(range(1,n)); o_ap=orderings(AP,n)
    # scaled AP (covmin)
    p=2 if n%2==0 else (3 if n%3==0 else n); scaled=[p*k for k in range(1,n)]
    o_sc=orderings(scaled,n)
    print(f"   n={n}: AP {{1..{n-1}}} -> {o_ap} orderings; scaled AP (covmin, x{p}) -> {o_sc}; (n-1)!={math.factorial(n-1)}")
print("   => DILATION (x p) PRESERVES the ordering count (M-invariant AND order-invariant): the covmin is the")
print("      AP's dilate, same orderings, same M. The observer (0) is FIXED by dilation; runners scale.")
print()
print("(3) COPY-TO-ALL-N reframe: make every point an observer. 'all lonely' = min PAIRWISE gap maximized =")
print("    the COVERING RADIUS = the equally-spaced cusp = covmin 1/n. The single-observer LRC, summed over")
print("    all n points, IS the covering-min (what I just solved). Tournament dual: sum OCF over all v =")
print("    sum_C |C| mu(C) (each odd cycle weighted by its length) -- the all-vertex odd-cycle census.")
print()
print("(2) TRANSLATE reframe: observer at c!=0 -> inhomogeneous LRC M_c(S)=max_t min_s ||s t - c||. Tournament")
print("    is VERTEX-TRANSITIVE so translating the marked vertex is trivial (relabel) -- UNLESS we read")
print("    translation as the COMPLEMENT/REVERSAL R (flip all arcs = shift the observer's signature by the")
print("    antipode). So: LRC translation <-> tournament complement (both move the observer's frame).")
print()
print("(1) SAME-CHANGE analogy: the analogy SURVIVES under DILATION (both have M(pS)=M(S)/order-invariance and")
print("    the OCF is scale-free) and under COPY-TO-ALL (covmin <-> OCF-sum); it BREAKS under raw TRANSLATION")
print("    (LRC origin special, tournament vertex-transitive) unless translation:=complement. The robust shared")
print("    symmetry is DILATION + the marked observer; the fragile one is translation.")
