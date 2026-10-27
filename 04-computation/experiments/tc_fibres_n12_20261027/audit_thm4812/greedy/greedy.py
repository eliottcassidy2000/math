"""greedy.py -- audit of THM-4812 (c), case m >= 19 (no script for this step exists in the repository packet).

For a regular (in particular vertex-transitive) tournament of odd order m:
  c3 = m(m^2-1)/24 cyclic triangles, t = (m^2-1)/8 through each vertex, a triangle meets at most D = 3(t-1) others,
  so the number i_k of k-sets of pairwise disjoint cyclic triangles satisfies i_k >= prod_{i<k}(c3 - i(D+1)) / k!
  (valid while every factor is positive), and H = I(Omega, 2) >= sum_k 2^k i_k.
Claim checked: sum_k 2^k i_k (lower bound) > 23 * 3^((m-1)/2) (Dixon) for every odd m in [19, 401] (exact rationals),
and the analytic tail for m > 401 with k = floor(m/20).
"""
from fractions import Fraction
import math

def greedy_lb(m):
    c3 = m * (m * m - 1) // 24
    assert c3 * 24 == m * (m * m - 1)
    t = (m * m - 1) // 8
    D = 3 * (t - 1)
    total = Fraction(1)       # k = 0
    prod = 1
    k = 0
    while True:
        fac = c3 - k * (D + 1)          # factor for choosing triangle number k+1
        if fac <= 0:
            break
        prod *= fac
        k += 1
        total += Fraction(2 ** k * prod, math.factorial(k))
    return total, k

worst = None
fails = []
for m in range(19, 402, 2):
    lb, kmax = greedy_lb(m)
    rhs = 23 * 3 ** ((m - 1) // 2)
    ratio = lb / rhs
    if lb <= rhs:
        fails.append(m)
    lr = math.log(lb.numerator) - math.log(lb.denominator) - math.log(rhs)
    if worst is None or lr < worst[1]:
        worst = (m, lr)
    if m in (19, 21, 23, 25, 27, 31, 41, 51, 101, 201, 401):
        print(f"m={m}: kmax={kmax}, log(greedy LB) - log(23*3^((m-1)/2)) = {lr:.3f}")
print(f"odd m in [19,401] with greedy LB <= 23*3^((m-1)/2): {fails}; minimum log-margin {worst[1]:.3f} at m = {worst[0]}")

# ---- analytic tail as stated in THM-4812: k = floor(m/20), c3 - k(D+1) >= 0.0229 m^3, 2^k i_k >= (0.0458 m^3)^k / k!
def stated_margin(m):
    k = m // 20
    return k * math.log(0.0458 * m ** 3) - math.lgamma(k + 1) - math.log(23) - (m - 1) / 2 * math.log(3)

def exact_single_term_margin(m):
    """log of the exact single-term bound 2^k prod_{i<k}(c3 - i(D+1))/k! with k = floor(m/20), minus log(23 * 3^((m-1)/2))."""
    k = m // 20
    c3 = m * (m * m - 1) // 24; t = (m * m - 1) // 8; D = 3 * (t - 1)
    s = k * math.log(2) - math.lgamma(k + 1)
    for i in range(k):
        fac = c3 - i * (D + 1)
        assert fac > 0
        s += math.log(fac)
    return s - math.log(23) - (m - 1) / 2 * math.log(3)

# check 0.0229 claim
bad0229 = [m for m in range(401, 20001, 2) if m * (m * m - 1) // 24 - (m // 20) * (3 * ((m * m - 1) // 8 - 1) + 1) < 0.0229 * m ** 3]
print(f"m in (401, 20000] violating c3 - k(D+1) >= 0.0229 m^3: {bad0229[:10]} (count {len(bad0229)})")
print(f"stated-tail margin at m = 401: {stated_margin(401):.3f}")
# monotonicity of the stated margin
dec = [(m, stated_margin(m), stated_margin(m + 2)) for m in range(401, 2001, 2) if stated_margin(m + 2) < stated_margin(m)]
print(f"number of odd m in [401, 2000) where the stated margin DEcreases from m to m+2: {len(dec)} (first: m={dec[0][0]}: {dec[0][1]:.3f} -> {dec[0][2]:.3f})")
mins = min((stated_margin(m), m) for m in range(403, 200001, 2))
print(f"minimum of the stated-tail margin over odd m in (401, 200000]: {mins[0]:.3f} at m = {mins[1]}")
mins2 = min((exact_single_term_margin(m), m) for m in range(403, 20001, 2))
print(f"minimum of the exact single-term margin over odd m in (401, 20000]: {mins2[0]:.3f} at m = {mins2[1]}")

# rigorous replacement for the tail: worst m in each block of constant k is m = 20k + 19.  Lower bound
#   F(20k+19) >= k (log 0.0458 + 3 log 20 + 2 log k + 1 - 10 log 3) - 0.5 log(2 pi k) - 1/(12k) - log 23 - 9 log 3
# using log k! <= k log k - k + 0.5 log(2 pi k) + 1/(12k) and 3 log(20k+19) >= 3 log(20k).
def tail_lb(k):
    return k * (math.log(0.0458) + 3 * math.log(20) + 2 * math.log(k) + 1 - 10 * math.log(3)) - 0.5 * math.log(2 * math.pi * k) - 1 / (12 * k) - math.log(23) - 9 * math.log(3)
print("rigorous tail lower bound F(20k+19) >= tail_lb(k): k=20: %.3f, k=21: %.3f, k=30: %.3f, k=100: %.3f" % (tail_lb(20), tail_lb(21), tail_lb(30), tail_lb(100)))
# derivative of tail_lb in k: (log 0.0458 + 3 log 20 + 2 log k + 1 - 10 log 3) + 2 - 0.5/k + 1/(12 k^2) > 0 for k >= 20
print("d/dk tail_lb at k = 20: %.3f (increasing for k >= 20 since the bracket grows with log k)" % (math.log(0.0458) + 3 * math.log(20) + 2 * math.log(20) + 1 - 10 * math.log(3) + 2 - 0.5 / 20))
# consistency: tail_lb(k) <= stated_margin(20k+19) for k in range
print("tail_lb(k) <= stated margin at m = 20k+19 for k = 20..2000:", all(tail_lb(k) <= stated_margin(20 * k + 19) + 1e-9 for k in range(20, 2001)))
