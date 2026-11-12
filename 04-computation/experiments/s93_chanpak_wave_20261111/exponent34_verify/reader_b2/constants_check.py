# reader_b2: arithmetic of the explicit constants in the proof of (B2) of THM-5141 (pure arithmetic, no search).
# Labels: every line below is an exact/float evaluation of a closed formula appearing in the reader's re-derivation.
import math
from fractions import Fraction as Fr

print("== (A)4 constant: (5 sqrt2/3) * sum_j (2 sqrt2/3)^j  (must be < 42)")
r = 2*math.sqrt(2)/3
print("   ratio 2sqrt2/3 =", r, " value =", (5*math.sqrt(2)/3)/(1-r))

print("== van der Corput step: (8 lam L + 2)(4/sqrt(lam) + 3) <= 56 L sqrt(lam) + 8/sqrt(lam) + 6 for 0 < lam < 1")
worst = 0.0
for L in (1, 2, 7, 64, 1000, 10**6):
    for e in range(-400, 0):
        lam = 10**(e/40)
        lhs = (8*lam*L + 2)*(4/math.sqrt(lam) + 3)
        rhs = 56*L*math.sqrt(lam) + 8/math.sqrt(lam) + 6
        worst = max(worst, lhs/rhs)
print("   max lhs/rhs on grid =", worst, "(<= 1 required)")

print("== Erdos-Turan chain with H = ceil(L |a|^(-1/3)), L >= |a|^(1/3) >= 1:")
print("   count - 2 dl L <= 12 L/H + (8/pi) * sum_{h<=H} (1/h) (56 L lam_h^(1/2) + 8 lam_h^(-1/2) + 6)")
z32 = sum(h**-1.5 for h in range(1, 2000000)) + 2/math.sqrt(2000000)   # zeta(3/2) with tail bound
c1 = 12 + (8/math.pi)*56*2*math.sqrt(2)      # coefficient of |a|^(1/3)
c2 = (8/math.pi)*8*math.sqrt(2)*z32          # coefficient of L^(3/2) |a|^(-1/2)   (2^ds <= 2)
c3 = (8/math.pi)*6*(1 + math.log(2))         # constant
c4 = (8/math.pi)*6                           # coefficient of ln L
print("   zeta(3/2) <=", z32)
print("   coefficients: |a|^(1/3): %.2f   L^(3/2)|a|^(-1/2): %.2f   const: %.2f   ln L: %.2f" % (c1, c2, c3, c4))
print("   => count <= 2 dl L + 416 |a|^(1/3) + 76 L^(3/2) |a|^(-1/2) + 26 + 16 ln L")

print("== number of odd integers in (L, 2L] is <= L for real L >= 1 (grid check)")
bad = 0
for i in range(1000, 400001):
    L = i/1000.0
    lo = math.floor(L) + 1; hi = math.floor(2*L)
    n = sum(1 for k in range(lo, hi+1) if k % 2) if hi - lo < 50 else ((hi + 1)//2 - lo//2)
    if n > L: bad += 1
print("   violations:", bad)

print("== the sum: 6 * 2^(2/5) =", 6*2**0.4, " (< 8)")
q = 6*2**0.4/8
m = max(((1+s)**3)*q**s for s in range(1, 5000))
arg = max(range(1, 5000), key=lambda s: ((1+s)**3)*q**s)
print("   sup_s (1+s)^3 (6*2^(2/5)/8)^s = %.1f at s = %d  (finite: the factor w in (B2) is not needed)" % (m, arg))
m2 = max(((1+s)**2)*q**s for s in range(1, 5000))
print("   sup_s (1+s)^2 (6*2^(2/5)/8)^s = %.1f" % m2)

print("== exponents: min(A K^(3/2), B/K) <= A^(2/5) B^(3/5):  (3/2)(2/5) - (3/5) =", Fr(3,2)*Fr(2,5) - Fr(3,5))
print("   tail: sum_{i>=s2} 2^-i 4^w = 2^(2w - s2 + 1) <= 2^(3w/2 + 3/2);  8^s1 <= 2^(3w/2 + 3/2)")
for w in range(4, 12):
    s1 = (w+1)//2; s2 = w//2
    assert 2*w - s2 + 1 <= 1.5*w + 1.5 + 1e-9 and 3*s1 <= 1.5*w + 1.5 + 1e-9
print("   checked for w = 4..11 (parity pattern)")
