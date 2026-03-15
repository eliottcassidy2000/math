#!/usr/bin/env python3
"""Analyze W(n) EGF structure: f(x) = Σ W(n) x^n/n!.
Try to decompose f(x) = 1/(1-x) + corrections.
Since CV²(n) = W(n)/n! - 1, we have f(x) = 1/(1-x) + g(x) where g(x) = Σ (W(n)/n! - 1) x^n.
"""
from fractions import Fraction as F
import math

wn = {1: 1, 2: 2, 3: 8, 4: 32, 5: 158, 6: 928, 7: 6350, 8: 49752, 9: 439670,
      10: 4327904, 11: 46963358, 12: 556953448, 13: 7166360054,
      14: 99428495088, 15: 1479600188798, 16: 23506712352248,
      17: 397095175477430, 18: 7107209383674112, 19: 134345623603516190,
      20: 2674426516381764744, 21: 55925062706620208438,
      22: 1225582324497129993488, 23: 28088043261491650347134,
      24: 671901551280362054066968, 25: 16746599265666151628174198,
      26: 434187457363955400414175008, 27: 11692423738081050318010736030}

# W(n)/n! - 1 = CV²(n)
# f(x) = Σ W(n) x^n/n! = 1/(1-x) + Σ CV²(n) x^n
# The correction is: Σ (W(n)/n! - 1) x^n

# Since CV²(n) ~ 2/n, and Σ 2/n * x^n = -2 ln(1-x), we expect
# the correction to behave like -2 ln(1-x) near x=1.

# Let's check: define R(n) = W(n) - n!. Then R(n)/n! = CV²(n).
# R(n)/n! ~ 2/n, so R(n) ~ 2(n-1)!

print("=== W(n) - n! decomposition ===")
print("R(n) = W(n) - n!, R(n)/(n-1)!:")
for n in range(3, 28):
    R = wn[n] - math.factorial(n)
    ratio = F(R, math.factorial(n-1))
    print(f"  n={n:2d}: R(n)/(n-1)! = {float(ratio):.10f}")

# So R(n)/(n-1)! → 2. Let's subtract: S(n) = R(n) - 2(n-1)! = W(n) - n! - 2(n-1)!
print("\n=== S(n) = W(n) - n! - 2(n-1)!, S(n)/(n-2)!: ===")
for n in range(3, 28):
    S = wn[n] - math.factorial(n) - 2*math.factorial(n-1)
    ratio = F(S, math.factorial(n-2))
    print(f"  n={n:2d}: S(n)/(n-2)! = {float(ratio):.10f}")

# This should converge too. Let's define T(n) = S(n) - c(n-2)! and find c
print("\n=== S(n)/(n-2)! values (looking for limit): ===")
for n in range(5, 28):
    S = wn[n] - math.factorial(n) - 2*math.factorial(n-1)
    ratio = F(S, math.factorial(n-2))
    print(f"  n={n:2d}: {float(ratio):.10f}")

# It looks like S(n)/(n-2)! → some constant. Let me see if it's related to 2+something.
# Actually, the expansion W(n)/n! = 1 + 2/n + c/n² + d/n³ + ...
# means W(n) = n! + 2(n-1)! + c(n-2)! + d(n-3)! + ...
# This is an asymptotic series in falling factorials.

# Let me compute the coefficients by successive subtraction
print("\n=== Asymptotic expansion W(n) = Σ a_j (n-j)! ===")
print("a_0 = 1 (from n!), a_1 = 2 (from 2(n-1)!)")
residual = dict()
for n in range(1, 28):
    residual[n] = F(wn[n])

# Subtract a_0 * n!
for n in range(1, 28):
    residual[n] -= F(math.factorial(n))

# Subtract a_1 * (n-1)!
a1 = F(2)
for n in range(1, 28):
    residual[n] -= a1 * F(math.factorial(n-1))

# Now residual[n] ~ a2 * (n-2)! for large n
print("\nresidual[n]/(n-2)! for a2:")
for n in range(3, 28):
    if n >= 2:
        ratio = F(residual[n], math.factorial(n-2))
        print(f"  n={n:2d}: {float(ratio):.10f}")

# The limit gives a2. Let me see what it converges to.
# n=27: let me get the exact fraction
n = 27
ratio_27 = F(residual[27], math.factorial(25))
print(f"\nExact a2 estimate from n=27: {ratio_27}")
print(f"Float: {float(ratio_27):.15f}")

# Seems to converge to about 2. Let me check more carefully
# Actually, from (CV²*n - 2)*n → something:
# CV² = W/n! - 1 = 2/n + c/n² + ...
# (CV²*n - 2)*n = c + ...
# From earlier output, this was about -0.20 at n=27, still decreasing

# So actually (CV²*n - 2)*n doesn't converge to a constant — it goes to 0 slowly
# That means c in 2/n + c/n² is actually 0, and we have 2/n + d/n^α for some α < 2
# Or maybe it's 2/n + c*ln(n)/n² or similar

# Let me re-examine: (CV² - 2/n) * n²
print("\n=== (CV² - 2/n) * n² ===")
for n in range(5, 28):
    cv2 = F(wn[n], math.factorial(n)) - 1
    correction = (cv2 - F(2, n)) * n * n
    print(f"  n={n:2d}: {float(correction):.10f}")

print("\n=== (CV² - 2/n) * n² / ln(n) ===")
for n in range(5, 28):
    cv2 = F(wn[n], math.factorial(n)) - 1
    correction = float((cv2 - F(2, n)) * n * n) / math.log(n)
    print(f"  n={n:2d}: {correction:.10f}")

# Now let me try a completely different approach: fit the EGF
# f(x) = Σ W(n)/n! x^n
# We know f(x) ~ 1/(1-x) near x=1 (since W(n)/n! → 1)
# More precisely, f(x) = 1/(1-x) - 2 ln(1-x) + C + O((1-x)ln(1-x))
# because Σ (2/n) x^n = -2 ln(1-x)

# Check: at what rate does f(x) - 1/(1-x) + 2 ln(1-x) converge?
# Numerically evaluate f at x=0.5:
print("\n=== Numerical evaluation of f(x) at x=0.5 ===")
x = 0.5
f_val = sum(wn[n] / math.factorial(n) * x**n for n in range(1, 28))
pole = 1/(1-x)  # = 2
log_term = -2 * math.log(1-x)  # = 2*ln(2) ≈ 1.3863
print(f"f(0.5) ≈ {f_val:.10f}")
print(f"1/(1-x) = {pole:.10f}")
print(f"-2*ln(1-x) = {log_term:.10f}")
print(f"f - 1/(1-x) = {f_val - pole:.10f}")
print(f"f - 1/(1-x) + 2*ln(1-x) = {f_val - pole + 2*math.log(1-x):.10f}")

# The Hertzsprung-connection means N(n,0) = A002464(n).
# For the EGF of NUD perms: exp(-x)/(1-x)²
# Let's verify: for W(n), the bivariate EGF at u=2 should give W(n)
# F(x,u) = Σ_n (Σ_j N(n,j) u^j) x^n/n!
# F(x,1) = exp(-x)/(1-x)²
# F(x,0) = EGF of Hertzsprung = Σ A002464(n) x^n/n!

# Let's compute F(x,u) at x=0.5 for u=0,1,2
print("\n=== Bivariate evaluation at x=0.5 ===")
# Need the N(n,j) data - recompute
from nud_adj1_dist import compute_nud_adj1_dist
all_dist = {}
for n in range(1, 16):
    all_dist[n] = compute_nud_adj1_dist(n)

for u in [0, 0.5, 1, 1.5, 2]:
    fu = sum(sum(u**j * all_dist[n][j] for j in range(len(all_dist[n]))) / math.factorial(n) * 0.5**n for n in range(1, 16))
    print(f"  F(0.5, {u}) ≈ {fu:.10f}")

# F(x,1) should equal exp(-x)/(1-x)²
nud_egf_at_05 = math.exp(-0.5) / (1-0.5)**2
print(f"  exp(-0.5)/(0.5)² = {nud_egf_at_05:.10f}")

# Let's try F(x,u) = exp((u-1)x) / (1-x)² ?
# At u=1: exp(0)/(1-x)² = 1/(1-x)², but NUD EGF is exp(-x)/(1-x)². No.
# Try F(x,u) = exp(-(2-u)x) / (1-x)^2?
# At u=1: exp(-x)/(1-x)² ✓
# At u=2: exp(0)/(1-x)² = 1/(1-x)², gives n!+1 for coeff n... W(n) ≈ n! but not exactly.
# W(n)/n! = 1 + 2/n + ... which is larger than 1, and 1/(1-x)² gives (n+1).
# So W(n) ≈ (n+1)! for this form? No, (n+1) not (n+1)!. Hmm.
# [x^n] 1/(1-x)² = n+1, so n!*(n+1) = (n+1)!. But W(n) ≈ n!, so this is way off.

# Try F(x,u) = exp(a(u)*x) / (1-x)^{b(u)}?
# At u=1: exp(-x)/(1-x)², so a(1) = -1, b(1) = 2.
# At u=0: Hertzsprung EGF (no closed form known)
# This approach won't work for a(u), b(u) because Hertzsprung doesn't have EGF of this form.

# Let me instead try: does F(x,u) satisfy a PDE?
# ∂F/∂x = ? in terms of F and its derivatives...
# This is hard without a closed form.

# Let me try a more data-driven approach: look for a recurrence in N(n,j)
print("\n=== Looking for bivariate recurrence ===")
# Try: N(n,j) = (n-1)*N(n-1,j) + (n-2)*N(n-2,j) + α*N(n-1,j-1) + β*N(n-2,j-1) + γ*N(n-3,j) + ...
# The NUD recurrence NUD(n) = (n-1)*NUD(n-1) + (n-2)*NUD(n-2) summed over j gives the base.
# The correction must be j-dependent.

# From the correction data, the j=0 correction grows like n*correction(n-1,0).
# Let me look at the correction / (n-2)! ratio for j=0:
corr_j0 = [-1, 2, 6, 12, 36, 180, 1164, 8772, 75132, 721284, 7670460, 89471556, 1135411068]
# These are for n=3,4,...,15

print("Correction at j=0 / (n-3)!:")
for i, c in enumerate(corr_j0):
    n = i + 3
    if n >= 3:
        ratio = c / math.factorial(max(n-3, 0)) if n >= 3 else 0
        print(f"  n={n}: {c} / {math.factorial(n-3)} = {ratio:.6f}")

# Ratio of consecutive corrections at j=0:
print("\nRatio corr(n,0)/corr(n-1,0):")
for i in range(1, len(corr_j0)):
    n = i + 3
    if corr_j0[i-1] != 0:
        print(f"  n={n}: {corr_j0[i]/corr_j0[i-1]:.6f}")
