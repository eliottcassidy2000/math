#!/usr/bin/env python3
"""Verifier check 5 (THM-5128 (B4)(iv), added by the corrections): the regions lemma for odd moduli that are not
prime (m = 9, 15, 21, 27), on random regions Omega = {a x + b y <= X} cut with intervals:
    max over regions and residue vectors w of | #{h in B_n, I, Omega : (x, y) = w mod m} - m^-2 #(B_n, I, Omega) |
divided by 2^(n (1 - theta_m)), theta_m at the rate of the modulus plus 0.005 (rates_exact/rates_by_modulus.txt).
The ratio must stay bounded (and not grow) in n.  Also (iii) for chi_S, S = {gcd(g, m) = d}, with W = w_0(h)."""
import math
import random

import numpy as np

random.seed(51284)
LOGH = 20
H = 1 << LOGH
psi = np.zeros(H, dtype=np.int64)
for n in range(1, LOGH + 1):
    lo, hi = 1 << (n - 1), 1 << n
    h = np.arange(lo, hi, dtype=np.int64)
    par = h >> 1
    psi[lo:hi] = np.where(h & 1 == 0, psi[par] + par, psi[par] + (par >> 1) + 1)
hh = np.arange(H, dtype=np.int64)
x = 2 * psi + 1
y = 2 * hh + 1
rates = {}
for line in open("C:/tmp/wt-s91/04-computation/experiments/s91_f1_local_densities_20261110/rates_exact/rates_by_modulus.txt"):
    if line.startswith("#") or not line.strip():
        continue
    p_ = line.split()
    rates[int(p_[0])] = float(p_[1])

ok = True
for m in (9, 15, 21, 27):
    rho = rates[m] + 0.005
    theta = math.log(1 / rho) / math.log(2 / rho ** 2)
    row = []
    for n in (12, 14, 16, 18, 20):
        lo, hi = 1 << (n - 1), 1 << n
        xs, ys = x[lo:hi], y[lo:hi]
        code = (xs % m) * m + (ys % m)
        worst = 0.0
        for trial in range(30):
            a = random.randint(1, 50)
            b = random.randint(a, 4 * a + 50)
            hcut = random.uniform(lo, hi)
            X = (a * random.uniform(0.5, 1.0) + b) * (2 * hcut + 1)
            i0 = random.randint(lo, hi - 1)
            i1 = random.randint(i0, hi)
            sel = (a * xs + b * ys <= X)
            sel[: i0 - lo] = False
            sel[i1 - lo:] = False
            tot = int(sel.sum())
            cnt = np.bincount(code[sel], minlength=m * m)
            worst = max(worst, float(np.abs(cnt - tot / (m * m)).max()) / 2 ** (n * (1 - theta)))
        row.append(worst)
    grow = row[-1] <= 1.5 * max(row[:2])
    ok &= max(row) < 1.0 and grow
    print("%s m = %2d (rate %.4f, theta_m = %.4f): max |dev| / 2^(n(1 - theta_m)) for n = 12, 14, 16, 18, 20: %s"
          % ("PASS" if (max(row) < 1.0 and grow) else "FAIL", m, rates[m], theta, " ".join("%.3f" % v for v in row)))

# (iii) with chi_S for S = {gcd(g, m) = d}: normalised block sums of chi_d w_0 over B_n (W_0 = 2^-n)
g = np.gcd(x, y)
xf, yf = x.astype(float), y.astype(float)
w0 = np.log(2 * (xf + yf) / (xf + 2 * yf)) / xf
for m in (15, 9):
    rho = rates[m] + 0.005
    theta = math.log(1 / rho) / math.log(2 / rho ** 2)
    divs = [d for d in range(1, m + 1) if m % d == 0]
    worst_all = 0.0
    for d in divs:
        ind = (np.gcd(g, m) == d)
        # number of residue vectors (x, y) mod m with gcd(gcd(x, y), m) = d, by enumeration
        nS = sum(1 for a_ in range(m) for b_ in range(m) if math.gcd(math.gcd(a_, b_), m) == d)
        chi = ind - nS / (m * m)
        vals = []
        for n in range(10, LOGH + 1):
            lo, hi = 1 << (n - 1), 1 << n
            s = float((chi[lo:hi] * w0[lo:hi]).sum())
            vals.append(abs(s) * 2 ** (theta * n))          # W_0 2^(n(1-theta)) = 2^(-theta n)
        worst_all = max(worst_all, max(vals))
    ok &= worst_all < 1.0
    print("%s m = %2d: max over d | m and n = 10..20 of |sum over B_n of chi_d w_0| 2^(theta_m n) = %.4f"
          % ("PASS" if worst_all < 1.0 else "FAIL", m, worst_all))
print("ALL PASS" if ok else "SOME FAIL")
