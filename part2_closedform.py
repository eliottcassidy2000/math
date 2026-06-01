from fractions import Fraction as F
from math import gcd
from lrc import measure_intersection

# Geometric derivation of measure(B_a ∩ B_b), gcd(a,b)=1.
# As t runs [0,1), the point P(t)=(a t mod 1, b t mod 1) traverses the closed
# geodesic of rational slope on T^2 = (R/Z)^2; uniformly w.r.t. arclength param.
# B_a ∩ B_b  <=>  P(t) in box Q = (-1/n,1/n) x (-1/n,1/n)  (mod 1, centered at 0).
# The line {(a t, b t)} mod 1 = {(x,y): b x - a y = 0 mod 1} since gcd(a,b)=1 it's
# a single closed loop covering the torus densely-but-it's-rational => it's the
# subgroup generated, a line that wraps and is EQUIDISTRIBUTED on the torus only
# in the limit; for fixed a,b it is the finite line.
#
# measure_t(P(t) in Q) = (length of line ∩ Q) / (total length) where length uses
# t-parametrization: dt. Speed in (x,y) is (a,b), |.|=sqrt(a^2+b^2), total t-length 1.
# The line {b x - a y = 0 mod 1} intersect box. Easier: directly the EXACT value
# computed already. Let's instead VERIFY a clean closed form via the "three-distance"
# / continued-fraction-free count:
#
# measure(B_a ∩ B_b) = (1/(a b)) * Area-type count? Test: is it = (number of lattice
# translates of Q hit) * (something)? We just confirm a clean GENERAL bound instead:
#
# KEY RIGOROUS BOUND (provable): measure(B_a ∩ B_b) <= measure(B_a)=2/n trivially,
# and the line crosses the box; an exact upper bound for the (a,2a)=max is 1/n,
# i.e. ratio <= n/4. Verify ratio <= n/4 over ALL coprime pairs, all n.

print("Check: max correlation ratio over coprime pairs equals n/4 (attained by (a,2a))?")
for n in range(3,13):
    best=F(0); arg=None
    for a in range(1,30):
        for b in range(a+1,60):
            if gcd(a,b)!=1: continue
            r=measure_intersection([a,b],n)*F(n*n,4)
            if r>best: best=r; arg=(a,b)
    print(f"  n={n:2d}: max ratio found={float(best):.4f}  n/4={float(F(n,4)):.4f}  argmax={arg}  equal={best==F(n,4)}")
print()

# So the SINGLE PAIR overlap can be as large as 1/n (ratio n/4), achieved exactly when
# vj = 2 vi (or vi=2vj): the doubling resonance. This GROWS with n.
# But note: this needs vj=2vi to be IN the speed set.

# Closed form for (1,2) overlap = 1/n exactly. Why: B_1=(-1/n,1/n) one arc.
# B_2 arcs: (-1/(2n),1/(2n)) and (1/2-1/(2n),1/2+1/(2n)). Intersection with (-1/n,1/n):
# the arc around 0: (-1/(2n),1/(2n)) fully inside (-1/n,1/n) -> width 1/n. =1/n. CONFIRMED.

# MINIMUM: scan minimum ratio over coprime pairs.
print("Minimum correlation ratio over coprime pairs (a,b), b<=200:")
for n in range(3,13):
    worst=F(10); arg=None
    for a in range(1,15):
        for b in range(a+1,200):
            if gcd(a,b)!=1: continue
            r=measure_intersection([a,b],n)*F(n*n,4)
            if r<worst: worst=r; arg=(a,b)
    # is min = ((n-2)/n)... ? for (1,n-1):
    print(f"  n={n:2d}: min ratio={float(worst):.4f}={worst}  argmin={arg}")
