#!/usr/bin/env python3
"""audit_sec8.py -- independent check of the numbers quoted in section 8 of the note (folds and bands; table rows)."""
import math
from fractions import Fraction as F
def V(x): return math.pi ** (x / 2) / math.gamma(x / 2 + 1)
print("ball volumes V_n, n = 1..8:", [round(V(n), 4) for n in range(1, 9)], "; integer maximum at n =", max(range(1, 30), key=V))
lo, hi = 4.0, 7.0
for _ in range(200):
    a = lo + (hi - lo) / 3; b = hi - (hi - lo) / 3
    if V(a) < V(b): lo = a
    else: hi = b
print(f"real maximum of V at x = {lo:.5f}, V = {V(lo):.5f}")
print("even dimensions pi^k/k!, k = 1..5:", [round(math.pi ** k / math.factorial(k), 4) for k in range(1, 6)], "; maximum at k =", max(range(1, 20), key=lambda k: math.pi ** k / math.factorial(k)))
print("n^(1/n), n = 1..6:", [round(n ** (1 / n), 4) for n in range(1, 7)], "; integer maximum at n =", max(range(1, 30), key=lambda n: n ** (1 / n)), "; real maximum at e")
def busch(n):   # least number of Hamiltonian paths of a strong tournament on n >= 3 vertices (Busch 2006)
    if n % 3 == 0: return 3 * 5 ** ((n - 3) // 3)
    if n % 3 == 1: return 5 ** ((n - 1) // 3)
    return 9 * 5 ** ((n - 5) // 3)
maxH = {2: 1, 3: 3, 4: 5, 5: 15, 6: 45, 7: 189, 8: 661}      # A038375
lam = 2 * 5 ** (1 / 3)
print(f"lambda = 2*5^(1/3) = {lam:.4f}")
for n in range(2, 9):
    mean = F(math.factorial(n), 2 ** (n - 1)); f = busch(n + 1)
    kappa = {1: 1.0, 0: 3 / 5 ** (2 / 3), 2: 9 / 5 ** (4 / 3)}[(n + 1) % 3]
    print(f"  n={n}: f(n+1) = {f}, E[H_n] = {mean}, ratio = {float(f / mean):.3f}, (kappa/2) lambda^n/n! = {kappa / 2 * lam ** n / math.factorial(n):.3f}"
          + (f", f(n+1)/max H(n) = {F(f, maxH[n])}" if n in maxH else ""))
print("f(n+1) > max H(n):", {n: busch(n + 1) > maxH[n] for n in range(2, 8)}, "; equal at n = 5:", busch(6) == maxH[5])
print("bands [f(n), max H(n)], n = 3..7:", [(busch(n), maxH[n]) for n in range(3, 8)], "-> 7 lies between the bands of orders 4 and 5; 21 lies INSIDE the band of order 6")
def dfact(k):
    r = 1
    while k > 1: r *= k; k -= 2
    return r
print("share bound (n-1)!! 2^(-n(n-2)/4), n = 2, 4, 6, 8, 10, 12:", [str(F(dfact(n - 1), 2 ** (n * (n - 2) // 4))) for n in (2, 4, 6, 8, 10, 12)], "(monotone decreasing: no interior peak)")
print("odd n: n!! 2^(-(n-1)^2/4), n = 3, 5, 7, 9:", [str(F(dfact(n), 2 ** ((n - 1) ** 2 // 4))) for n in (3, 5, 7, 9)])
print("square term bound at n = 12 in units of 4^-12: 11!! 2^-30 4^12 =", dfact(11), "/ 64 =", dfact(11) / 64)
print("ratios D_n/D_(n+1), n = 8..11:", [round(a / b, 2) for a, b in ((1361 / 32768, 63103 / 2 ** 23), (63103 / 2 ** 23, 650471 / 2 ** 28), (650471 / 2 ** 28, 260546929 / 2 ** 39), (260546929 / 2 ** 39, 5.8917e-5))])
print("upper end of a strong window: largest n with f(n) <= h is about 1 + 3 ln h / ln 5 =", f"1 + {3 / math.log(5):.3f} ln h; at h = 223: {1 + 3 * math.log(223) / math.log(5):.1f} (window LENGTH there is at most 4)")
