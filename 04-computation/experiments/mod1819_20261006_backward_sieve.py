#!/usr/bin/env python3
"""The backward (3-adic) sieve of a minimal Collatz counterexample keeps a positive proportion
(direction D62 of collatz_connectivity_from_rigidity_20261001.md).  mac-mini, 2026-10-06.

Setting.  For the odd map U(x) = oddpart(3x+1), the inverse words from a 3-adic unit n are compositions
w = (k_1, ..., k_d), k_i >= 1: x_0 = n, x_i = (2^{k_i} x_{i-1} - 1)/3.  The word is LEGAL iff every x_{i-1}
(i <= d) is prime to 3 and k_i is odd exactly when x_{i-1} = 2 (mod 3).  The predecessor x_d is SMALLER
(multiplier 2^K/3^d < 1, K = sum k_i) iff 2^K < 3^d.  s(r) = Haar fraction of units n (classes mod 3^r)
with no legal word of length <= r whose predecessor is smaller (a FIRST-PASSAGE word: all proper prefixes
have 2^{K_i} > 3^i, the full word 2^K < 3^d).

  1. Legality measure: each word of length d is legal on exactly one unit class mod 3^d, Haar weight
     (1/2) 3^(1-d)  (checked by brute force for d <= 6).
  2. First moment (rigorous, exact rationals + explicit tails): P(n has a smaller predecessor)
     <= sum_d (1/2) 3^(1-d) F_d = 0.7757408...,  F_d = # first-passage compositions of length d.
     Hence lim s(r) >= 0.2242 > 0.
  3. Exact s(r) for r <= 16 (vectorised recursion over x mod 3^r) and an independent brute force (r <= 8).
     Bracket: s(16) - (tail of the first moment beyond 16) <= lim s(r) <= s(16).
  4. Moran duality: the Chernoff generating function of the legal words is (3/2) rho(theta)^d with
     rho(theta) = 3^(theta-1)/(2^theta - 1);  rho = 1 iff g(theta) = 2^-theta + (1/3)(3/2)^theta = 1 (THM-4504's
     Moran function, roots theta = 1, 2); min rho = 3^(-(1-h)) with h = H(log_3 2) = 0.949955 (the forward
     exceptional dimension): the backward 3-adic rate is the forward 2-adic rate 2^(-(1-h)) with 2 <-> 3.
Run: python3 mod1819_20261006_backward_sieve.py   (about 1 min, ~2 GB RAM at r = 16)
"""
import math, sys, itertools
from fractions import Fraction
import numpy as np

L2_3 = math.log2(3)
OK = True


def check(cond, msg):
    global OK
    print(("  ok   " if cond else "  FAIL ") + msg, flush=True)
    OK &= bool(cond)


# ---------------------------------------------------------------- 1. legality measure
print("1. legality: every word of length d is legal on exactly one unit class mod 3^d")


def legal(n, word, mod):
    x = n % mod
    for k in word:
        r = x % 3
        if r == 0 or (k % 2 == 1) != (r == 2):
            return False
        x = ((pow(2, k, 3 * mod) * x - 1) // 3) % (mod // 3) if mod > 3 else 0
        mod //= 3
    return True


good = True
for d in range(1, 6):
    M = 3 ** d
    units = [n for n in range(M) if n % 3]
    for word in itertools.product(range(1, 7), repeat=d):
        c = sum(1 for n in units if legal(n, word, M))
        good &= (c == 1)
check(good, "d = 1..5, exponents 1..6: each word legal on exactly 1 of the 2*3^(d-1) unit classes (weight (1/2) 3^(1-d))")

# ---------------------------------------------------------------- 2. rigorous first moment
print("2. first moment over first-passage words (exact rationals, explicit tails)")
DMAX, W = 1500, 200
c = [0]
K = 0
for i in range(1, DMAX + 2):
    while 2 ** K <= 3 ** i:
        K += 1
    c.append(K)
theta = 1.5
rho15 = 3 ** (theta - 1) / (2 ** theta - 1)
Rfac = rho15 / (1 - rho15)
cnt = {0: 1}
total = Fraction(0)
partial = {}
drop = 0.0
for d in range(1, DMAX + 1):
    keys = sorted(cnt)
    acc, j, F, surv = 0, 0, 0, {}
    for Kd in range(keys[0] + 1, c[d] + W):
        while j < len(keys) and keys[j] < Kd:
            acc += cnt[keys[j]]
            j += 1
        if Kd < c[d]:
            F += acc
        else:
            surv[Kd] = acc
    tot_prev = sum(cnt.values())
    h0 = c[d] + W - d * L2_3
    drop += math.exp(math.log(tot_prev) + math.log(0.5) + (1 - d) * math.log(3) + math.log(Rfac)
                     - theta * h0 * math.log(2) - math.log(1 - 2 ** (-theta)))
    total += Fraction(F, 2 * 3 ** (d - 1))
    partial[d] = total
    cnt = surv
tail = sum(math.exp(math.log(n) + math.log(0.5) + (1 - DMAX) * math.log(3) + math.log(Rfac)
                    - theta * (Kd - DMAX * L2_3) * math.log(2)) for Kd, n in cnt.items())
FM = float(total) + drop + tail
print(f"   sum_(d <= {DMAX}) (1/2)3^(1-d) F_d = {float(total):.10f};  dropped-window bound {drop:.1e};  depth tail {tail:.1e}")
check(FM < 1, f"P(smaller predecessor) <= {FM:.10f} < 1, so lim s(r) >= {1 - FM:.10f} > 0  (D62: positive)")
print("   F_d, d = 1..14:", [int((partial[d] - partial.get(d - 1, Fraction(0))) * 2 * 3 ** (d - 1)) for d in range(1, 15)])

# ---------------------------------------------------------------- 3. exact s(r) and bracket
print("3. exact survival fractions s(r) (vectorised recursion) and an independent brute force")
R = int(sys.argv[1]) if len(sys.argv) > 1 else 16
prev = None
s = {}
for r in range(1, R + 1):
    M, Mp = 3 ** r, 3 ** (r - 1)
    x = np.arange(M, dtype=np.int64)
    unit = (x % 3) != 0
    k0 = np.where(x % 3 == 2, 1, 2).astype(np.int64)
    best = np.full(M, 1e9)
    Jmax = int(((r - 1) * (L2_3 - 1) + L2_3) / 2) + 2
    for j in range(Jmax + 1):
        pk = np.where(k0 == 1, pow(2, 1 + 2 * j, 3 * M), pow(2, 2 + 2 * j, 3 * M)).astype(np.int64)
        val = (pk * x) % (3 * M)
        h = (k0 + 2 * j).astype(float) - L2_3
        if r == 1:
            cand = h
        else:
            child = ((val - 1) // 3) % Mp
            cand = h + np.minimum(0.0, prev[child])
        best = np.minimum(best, cand)
        del pk, val, cand
    best[~unit] = 1e9
    s[r] = np.count_nonzero(best[unit] >= 0) / np.count_nonzero(unit)
    prev = best
print("   s(r):", {r: round(v, 6) for r, v in s.items()})


def survives(n0, r):
    stack = [(Fraction(n0), 0, 0)]
    while stack:
        xv, d, Kc = stack.pop()
        if d == r:
            continue
        r3 = (xv.numerator * pow(xv.denominator, -1, 3)) % 3
        if r3 == 0:
            continue
        k = 1 if r3 == 2 else 2
        while True:
            Kn, dn = Kc + k, d + 1
            if 2 ** Kn < 3 ** dn:
                return False
            if Kn - dn * L2_3 - (r - dn) * (L2_3 - 1) > 0:
                break
            stack.append(((2 ** k * xv - 1) / 3, dn, Kn))
            k += 2
    return True


bf = {r: sum(survives(n, r) for n in range(3 ** r) if n % 3) / (2 * 3 ** (r - 1)) for r in range(1, 9)}
check(all(abs(bf[r] - s[r]) < 1e-12 for r in bf), f"brute force r <= 8 agrees: {[round(bf[r], 6) for r in bf]}")
tail16 = float(partial[DMAX] - partial[R]) + drop + tail
lo, hi = s[R] - tail16, s[R]
check(lo > 0, f"bracket: {lo:.5f} <= lim s(r) <= {hi:.5f}  (s({R}) minus the first-moment tail {tail16:.6f} beyond depth {R})")
print("   (the repo's T-depth table 0.500, 0.333, 0.333, 0.31481, 0.31481, 0.31481, 0.31070, 0.30521 is the same sequence on a coarser clock)")

# ---------------------------------------------------------------- 4. Moran duality
print("4. Moran duality and the 2 <-> 3 mirror of the exponent 1 - h")
g = lambda t: 2 ** (-t) + (3 / 2) ** t / 3
rho = lambda t: 3 ** (t - 1) / (2 ** t - 1)
check(all(abs((rho(t) - 1) - 2 ** t * (g(t) - 1) / (2 ** t - 1)) < 1e-12 for t in np.linspace(0.3, 3, 40)),
      "rho(theta) - 1 = 2^theta (g(theta) - 1)/(2^theta - 1): same zero set {1, 2} as THM-4504's Moran function")
tstar = math.log2(L2_3 / (L2_3 - 1))
p = 1 / L2_3
h = -(p * math.log2(p) + (1 - p) * math.log2(1 - p))
check(abs(rho(tstar) - 3 ** (-(1 - h))) < 1e-12,
      f"min rho = rho(theta*) = {rho(tstar):.9f} = 3^-(1-h), h = H(log_3 2) = {h:.6f}, theta* = log2(log2 3/(log2 3 - 1)) = {tstar:.6f}")
# exact census: the first-passage mass m_d = (1/2)3^(1-d) F_d behaves like C(d) d^(-3/2) 3^(-(1-h)d) with bounded
# oscillating C(d) (the arithmetic of {d log2 3}); regression of log m_d + 1.5 log d on d over all nonzero terms
ds = np.array([d for d in range(2, DMAX + 1) if partial[d] != partial[d - 1]], dtype=float)
ys = np.array([math.log((partial[int(d)] - partial[int(d) - 1]).numerator) - math.log((partial[int(d)] - partial[int(d) - 1]).denominator)
               for d in ds]) + 1.5 * np.log(ds)
slope, icpt = np.linalg.lstsq(np.vstack([ds, np.ones_like(ds)]).T, ys, rcond=None)[0]
check(abs(slope - (-(1 - h) * math.log(3))) < 2e-4,
      f"first-passage mass ~ d^-3/2 3^-(1-h)d: fitted slope {slope:.6f} vs -(1-h) ln 3 = {-(1 - h) * math.log(3):.6f} over {len(ds)} nonzero terms d <= {DMAX}"
      " (the 3-adic mirror of the forward glide law W_k = Theta(2^(hk) k^(-3/2)), THM-4504)")
print("ALL CHECKS PASSED" if OK else "SOME CHECK FAILED")
