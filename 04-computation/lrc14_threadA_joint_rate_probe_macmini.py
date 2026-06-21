# Adversarial probe: does a SINGLE THM-546 peel of f2 give the O(1/f1) rate?
# THM-546: |Delta_w| <= (6/49) V(E')/w, V(E') <= 7*sigma(E') = 7*sum(e in E').
# Peel f2: E' = B u {f1}. sigma(E') ~ f1 (dominated by f1). w = f2 = gamma*f1.
# bound = (6/49)*7*f1 / (gamma*f1) = (6/7)/gamma = O(1) NOT O(1/f1).
from fractions import Fraction as F
gamma = F(2,1)
for f1 in (20,40,80,160,320):
    sigma_Eprime = sum([0,2,4,6,8,10,12]) + f1   # B + f1
    V = 7*sigma_Eprime
    w = int(gamma*f1)
    bound = F(6,49)*V/w
    print(f"f1={f1}: single-peel THM-546 bound (6/49)V(B u f1)/f2 = {float(bound):.4f}")
print()
print("=> a SINGLE THM-546 peel of f2 gives an O(1) error (~0.74), NOT O(1/f1).")
print("   The O(1/f1) rate is the EMPIRICAL convergence, not what one THM-546 peel proves.")
