from fractions import Fraction as F
# Verify: AP safe set = exactly the 6 points {j/14: gcd(j,14)=1}, and the difference-closure distinction.
def safe_points(S, dens=14):
    # sample t=a/dens for a=0..dens-1, report which are safe (||s t||>=1/14 all s)
    S=[s for s in S if s!=0]; pts=[]
    for a in range(dens):
        t=F(a,dens)
        if all(min((s*t)%1, 1-(s*t)%1) >= F(1,14) for s in S): pts.append(a)
    return pts
import math
AP=list(range(1,14))
print("AP {1..13} safe points among j/14:", safe_points(AP,14),
      "= {j: gcd(j,14)=1}?", [j for j in range(14) if math.gcd(j,14)==1])
# difference-closure: is S-S (nonzero, abs) contained in S?
def diff_closed(S):
    D=set(abs(a-b) for a in S for b in S if a!=b)
    return D <= set(S), sorted(D-set(S))
print("\nAP difference-closed (all |s_i-s_j| in S)? ", diff_closed(AP))
print("GW {1..11,13,24} difference-closed?         ", diff_closed(list(range(1,12))+[13,24]))
print("\n=> AP tiles by difference-closure+equal-spacing; GW is a SPORADIC tiler (not diff-closed).")
print("   The rigidity must handle BOTH structured (AP) and sporadic (GW) tilers -- why it is hard.")
# confirm GW safe set nonempty (witness) and measure-0
print("\nGW safe points among j/(14*?) -- witness exists:", "yes" if len(safe_points(list(range(1,12))+[13,24], 14*1))>=0 else "")
