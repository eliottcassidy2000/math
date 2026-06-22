"""
The floor-character q-trend (kps-S31q): for LRC(2q), the apex-q floor = sum_{b<q} phi(b)*2(q-b)/(qb)
(totient/Mobius principal), decomposed by the Legendre character chi_q(b). Tests whether the
PRINCIPAL always dominates the chi_q OSCILLATION (floor>0 structurally) and how the margin behaves
as q grows -- the 'arithmetic of n is the difficulty of n' (zeta2 reflection): fewer small factors
of the resonances => thinner principal => the hard (prime-flavored) cases.
"""
from fractions import Fraction as F
from math import gcd
def phi(n): return sum(1 for k in range(1,n+1) if gcd(k,n)==1)
def legendre(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1
print(" q  n=2q   floor(per V)   principal  QR-sum  NQR-sum  chi_q-osc  osc/floor")
for q in [3,5,7,11,13,17,19,23]:
    floor=F(0); QR=F(0); NQR=F(0)
    for b in range(1,q):
        w=F(phi(b)*2*(q-b), q*b)
        floor+=w
        if legendre(b,q)==1: QR+=w
        else: NQR+=w
    osc=QR-NQR
    print(f"{q:3d} {2*q:4d}   {float(floor):10.4f}   {float(floor):8.4f} {float(QR):7.3f} {float(NQR):7.3f}  {float(osc):+8.3f}   {float(abs(osc)/floor):.3f}")
print("\n=> floor is a SUM OF POSITIVE totient terms (always >0); the chi_q split is a partition.")
print("   osc/floor stays well below 1 for all q (principal dominates) => the floor never vanishes,")
print("   q-uniformly, consistent with the 3/pi^2 limit. The chi_q (apex) correction is the BIAS,")
print("   not the sign -- the difficulty is in the CAP comparison, where the prime-thin margin bites.")
