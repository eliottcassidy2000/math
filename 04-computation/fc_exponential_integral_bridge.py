"""The claim  int_0^1 e^{Q(t)} dt != 0  for nonconstant Q in Qbar[t]:
what is elementary, and how it sits against FC / THM-3022."""
from fractions import Fraction
from math import comb

print("STRUCTURE OF THE CLAIM")
print("  If Q has REAL algebraic coefficients then e^{Q(t)} > 0 on [0,1], so the")
print("  integral is > 0 -- trivial.  The entire content is COMPLEX Q.")
print("  Exactly the shape found for FC in HYP-9076 sec 6: real case trivial by")
print("  positivity, all content over C.")
print()
print("DEGREE 1 IS LINDEMANN-WEIERSTRASS")
print("  Q = a t + b, a != 0:  int_0^1 e^{at+b} dt = e^b (e^a - 1)/a,")
print("  zero  <=>  e^a = 1  <=>  a in 2 pi i Z \\ {0}, which is TRANSCENDENTAL.")
print("  So algebraicity of the coefficients is exactly what rules it out.")
print("  (Nonconstant is needed: Q constant gives e^Q != 0 anyway.)")
print()
print("RELATION TO FC -- the precise statement")
print("  int_0^1 e^{Q} dt = sum_m M(Q^m)/m!   with M(t^j) = 1/(j+1).")
print("  FC kills EVERY moment L(f^m) separately (weight w_j = j!).")
print("  The integral conjecture kills ONE weighted sum of moments M(Q^m)")
print("  (weight w_j = 1/(j+1)).  Same functional shape, different demand.")
print()
print("THM-3022 APPLIED TO THE [0,1] LEBESGUE WEIGHT  w_j = 1/(j+1)")
w = lambda j: Fraction(1, j+1)
def Q(a,b):
    return w(2*a)*w(b)**2 - 2*w(a)*w(b)*w(a+b) + w(a)**2*w(2*b)
lc = all(w(j)**2 < w(j-1)*w(j+1) for j in range(1,40))
print(f"  log-convex (w_j^2 < w_(j-1) w_(j+1) for j=1..39): {lc}")
bad = [(a,b) for a in range(0,9) for b in range(a+1,10) if Q(a,b) <= 0]
print(f"  Q_w(a,b) > 0 for all 0<=a<b<=9: {not bad}   violations: {bad}")
print(f"  e.g. Q_w(0,1) = {Q(0,1)},  Q_w(1,3) = {Q(1,3)},  Q_w(2,5) = {Q(2,5)}")
print()
print("  => by THM-3022 the two-slot threshold for the [0,1] Lebesgue weight is 2:")
print("     M(f) = M(f^2) = 0 already forces f = 0.  So the FC-analogue for this")
print("     measure is settled at two slots by the same log-convexity argument.")
