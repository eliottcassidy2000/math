#!/usr/bin/env python3
"""Auditor proof-5128-CF, check A2: THM-5128 (D1) against a direct count of triples (own code).

usage: python3 a2_mass.py K          (X = 2^k for k <= K; default 24)

Checks
  (a) the psi table against the plain recursion, 0 <= q <= h/2, |rho - (1 - phi/2)| < n/(2h+1)   (inputs of (D1))
  (b) T_s(X) for each padding s against (l_s - 2 l_(s+1) + l_(s+2)) X ln X, error normalised by (s + 1) X 2^-s
  (c) sheets: T_0(X), T_1(X) against kappa X ln X: (T - kappa X ln X)/X must stay bounded (claim O(X)), and the
      difference quotient [T(2X)/(2X) - T(X)/X]/ln 2 must tend to kappa
  (d) the word-pair weight s + 1: total constant pi^2/12.
"""
import sys
import time
import numpy as np
import mpmath as mp
from a_lib import psi_table, psi_slow, count_by_padding

K = int(sys.argv[1]) if len(sys.argv) > 1 else 24
t0 = time.time()
H = (1 << (K - 1)) + 2
psi = psi_table(H)
x = (2 * psi + 1).astype(np.int32)

# (a)
ok = all(psi[h] == psi_slow(h) for h in list(range(0, 5000)) + [H - 1, H - 2, 123456 % H, 7654321 % H])
print("(a) psi table equals the plain recursion (5000 + 4 values):", ok)
h = np.arange(H, dtype=np.int64)
q = h - psi
print("    0 <= q(h) <= h/2 for all h <", H, ":", bool(np.all(q >= 0) and np.all(2 * q <= h)))
# radical inverse phi(h) = sum_i h_i 2^(-i-1)
HH = min(H, 1 << 22)
hh = np.arange(1, HH, dtype=np.int64)
phi = np.zeros(HH - 1)
nb = np.zeros(HH - 1, dtype=np.int64)
t = hh.copy()
i = 0
while t.any():
    phi += (t & 1) * 2.0 ** (-i - 1)
    nb += (t > 0)
    t >>= 1
    i += 1
rho = (2 * psi[1:HH] + 1) / (2 * hh + 1)
dev = np.abs(rho - (1 - phi / 2)) * (2 * hh + 1) / nb
print("    max over 1 <= h < 2^22 of |rho - (1 - phi/2)| (2h+1)/n =", float(dev.max()), " (claim < 1)")

mp.mp.dps = 30


def l(j):
    return -mp.polylog(2, -mp.mpf(2) ** (-j))


SM = 12
T, diag = count_by_padding(K, x, SM)
print("census to X = 2^%d done in %.1f s" % (K, time.time() - t0))
ks = [k for k in range(10, K + 1, 2)]
if K not in ks:
    ks.append(K)

print("\n(b) per padding s:  T_s(X),  c_s = l_s - 2 l_(s+1) + l_(s+2),  E_s = (T_s - c_s X ln X) / ((s + 1) X 2^-s)")
print("  s      c_s     | " + " ".join("E_s(2^%d)" % k for k in ks))
worst = 0.0
for s in range(0, SM + 1):
    cs = float(l(s) - 2 * l(s + 1) + l(s + 2))
    row = []
    for k in ks:
        X = float(2 ** k)
        Ts = T[s, k]
        if s == 0:
            pass
        E = (Ts - cs * X * np.log(X)) / ((s + 1) * X * 2.0 ** (-s))
        worst = max(worst, abs(E))
        row.append("%+9.4f" % E)
    print(" %2d  %.8f | " % (s, cs) + " ".join(row))
print("  largest |E_s| in the table: %.4f  (the claim is that E_s is bounded uniformly in s and X)" % worst)

print("\n(c) sheets:  R_sigma(X) = (T_sigma - kappa_sigma X ln X)/X   and   slope = [T(2X)/(2X) - T(X)/X]/ln 2")
k0 = float(mp.pi ** 2 / 12 + 2 * mp.polylog(2, -mp.mpf(1) / 2) - mp.polylog(2, -mp.mpf(1) / 4))
k1 = float(mp.polylog(2, -mp.mpf(1) / 4) - mp.polylog(2, -mp.mpf(1) / 2))
T0 = T[0]
T1 = T[1:].sum(axis=0)
print("  kappa_0 = %.10f   kappa_1 = %.10f" % (k0, k1))
print("   k |      T_0(2^k)      R_0     slope_0 |      T_1(2^k)      R_1     slope_1 | share_0 (limit 0.43186)")
for k in range(10, K + 1):
    X = float(2 ** k)
    R0 = (T0[k] - k0 * X * np.log(X)) / X
    R1 = (T1[k] - k1 * X * np.log(X)) / X
    if k > 10:
        s0 = (T0[k] / X - T0[k - 1] / (X / 2)) / np.log(2)
        s1 = (T1[k] / X - T1[k - 1] / (X / 2)) / np.log(2)
        ss = "%.6f" % s0, "%.6f" % s1
    else:
        ss = "   -    ", "   -    "
    print("  %2d | %14d %+9.5f %s | %14d %+9.5f %s | %.5f" % (k, T0[k], R0, ss[0], T1[k], R1, ss[1],
                                                           T0[k] / (T0[k] + T1[k])))

print("\n(d) weight s + 1 (word pairs):  W(X) = sum_s (s + 1) T_s(X);  slope = [W(2X)/(2X) - W(X)/X]/ln 2 -> pi^2/12 = %.6f"
      % float(mp.pi ** 2 / 12))
# rows 0..SM individually, row SM+1 pools s > SM: recount those with exact weights
if K > SM:
    print("  (paddings s > %d are pooled; their weight is taken as s + 1 = %d, an underestimate of relative size < 2^-%d)"
          % (SM, SM + 2, SM))
W = np.zeros(K + 1)
for s in range(0, SM + 1):
    W += (s + 1) * T[s]
W += (SM + 2) * T[SM + 1]
for k in range(12, K + 1, 2):
    X = float(2 ** k)
    sl = (W[k] / X - W[k - 1] / (X / 2)) / np.log(2)
    print("  k = %2d  W = %14d  slope = %.6f   sheet-0 share of W = %.5f (limit 0.19641)" % (k, W[k], sl, T0[k] / W[k]))

print("\n(e) mean of r over odd N <= X against g_bar ln X:  2 (T_0 + T_1)/X - g_bar ln X   (claim: bounded)")
gbar = 2 * (k0 + k1)
for k in range(12, K + 1, 2):
    X = float(2 ** k)
    print("  k = %2d   %.5f" % (k, 2 * (T0[k] + T1[k]) / X - gbar * np.log(X)))
print("diagonal triples (h1 = h2) up to 2^%d: sheet 0: %d, sheet 1: %d" % (K, diag[0, K], diag[1:, K].sum()))
print("total time %.1f s" % (time.time() - t0))
