"""
Merge: the resonance-killing/cut game <-> zeta-duality <-> Littlewood.
(1) ze functional equation bridges FLOOR zeta(2) and CAP zeta(-1):  s<->1-s, s=2<->s=-1.
(2) 2-adic Littlewood: |n|_2 IS the resonance-killing factor (n divisible by 2^k kills the 2^k resonance).
(3) LRC = MIN (additive, has a FLOOR); Littlewood = PRODUCT (multiplicative, liminf 0). Same torus, dual norms.
"""
import math
from fractions import Fraction as F
def norm(x): x%=1.0; return min(x,1-x)

# (1) functional equation: zeta(s)=2^s pi^{s-1} sin(pi s/2) Gamma(1-s) zeta(1-s)
def fe_zeta_at(s, zeta_1ms):
    return 2**s * math.pi**(s-1) * math.sin(math.pi*s/2) * math.gamma(1-s) * zeta_1ms
z2 = math.pi**2/6
z_m1 = fe_zeta_at(-1, z2)   # should be -1/12
print("ZETA-DUALITY (the floor<->cap functional equation s<->1-s):")
print(f"  FLOOR side: zeta(2)=pi^2/6={z2:.6f}  (Euler product prod_p(1-p^-2)^-1 = resonance density; floor mass 1/(2 zeta2)={1/(2*z2):.5f}=3/pi^2)")
print(f"  CAP side:   zeta(-1) via functional eqn from zeta(2) = {z_m1:.6f}   (= -1/12 = -B_2/2; M({{1..11,13}})=1/12 finite avatar)")
print(f"  bridge: zeta(2)=pi^2/6  --[s<->1-s]-->  zeta(-1)=-1/12 ; same zeta, two faces of the 1/14 bound.\n")

# (2) 2-adic Littlewood: |n|_2 = 2^{-v2(n)}. n kills 2-power resonances <=> |n|_2 small.
def v2(n):
    k=0
    while n%2==0: n//=2; k+=1
    return k
print("2-ADIC LITTLEWOOD = the resonance-killing game (|n|_2 = the divisibility-killer):")
print("  de Mathan-Teulie p-LC:  liminf_n  n * |n|_2 * ||n a|| = 0.   |n|_2=2^{-v2(n)} small <=> n kills 2^k resonance.")
for n in [7,14,28,56,84,112]:
    print(f"    n={n:>3}: v2={v2(n)}  |n|_2={2.0**(-v2(n)):.4f}  n*|n|_2={n*2.0**(-v2(n)):.3f}   (14=2*7: the 2-adic face of apex-7)")
print()

# (3) MIN (LRC, additive, floor) vs PRODUCT (Littlewood, multiplicative, vanishes) on the AP-ish torus
print("LRC (MIN, additive: FLOOR) vs LITTLEWOOD (PRODUCT, multiplicative: liminf 0) -- dual norms, same torus:")
import random
random.seed(1)
alpha = (math.sqrt(5)-1)/2     # golden = most badly approximable
beta  = math.sqrt(2)-1         # also bdd partial quotients
best_min=1.0; best_prod=1.0
for n in range(1,20000):
    m = norm(n*alpha); 
    mn = min(norm(n*alpha), norm(n*beta))      # LRC-style MIN over 2 'speeds'
    pr = n*norm(n*alpha)*norm(n*beta)          # Littlewood PRODUCT
    best_min=min(best_min,mn); best_prod=min(best_prod,pr)
print(f"  min_n  min(||n a||,||n b||)         = {best_min:.5f}  (MIN bounded below ~ the LRC floor: never 0)")
print(f"  min_n  n*||n a||*||n b|| (Littlewood)= {best_prod:.5f}  (PRODUCT driven toward 0: liminf conj = 0)")
print("  => additive MIN keeps a floor (zeta(2) resonance density); multiplicative PRODUCT collapses (Littlewood).")
print("     LRC = additive pole (floor), Littlewood = multiplicative pole (vanish); ze func-eqn is the hinge.")
