#!/usr/bin/env python3
"""Verifier check 1 (THM-5128, corrected text of (B), 'Numerical rates').

Reads the receipt rates_exact/rates_by_modulus.txt and checks every number the corrected paragraph quotes:
  - 262 moduli: every odd m <= 521, and 1023, 1143;
  - the twelve named rates;
  - the band of the other primes up to 89 (least 0.7691 at 79, largest 0.8060 at 37);
  - the maximum 0.932039 attained exactly at the multiples of 7;
  - theta_p/2 at the bare rates and at rate + 0.005.
Also recomputes a few rates from the definition of T_k in (B) (dense eigenvalues) as an independent check.
"""
import math
import os
import sys

import numpy as np

ROOT = "C:/tmp/wt-s91/04-computation/experiments/s91_f1_local_densities_20261110"
path = os.path.join(ROOT, "rates_exact", "rates_by_modulus.txt")
rates = {}
with open(path) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split()
        rates[int(parts[0])] = float(parts[1])

ok = True


def check(name, cond):
    global ok
    print(("PASS " if cond else "FAIL ") + name)
    ok = ok and cond


print("rows:", len(rates))
check("262 moduli", len(rates) == 262)
want = set(range(3, 522, 2)) | {1023, 1143}
check("moduli = every odd m <= 521 plus 1023, 1143", set(rates) == want)
check("1039 not in the table (cited from bridges/b3_f1_data.out)", 1039 not in rates)

named = {3: 0.8660, 5: 0.7418, 7: 0.9320, 11: 0.8587, 13: 0.8272, 17: 0.8289, 19: 0.8053, 23: 0.7931,
         31: 0.9002, 43: 0.8342, 73: 0.8516, 127: 0.8754}
for m, v in named.items():
    check("named rate m=%d: table %.6f quoted %.4f" % (m, rates[m], v), abs(round(rates[m], 4) - v) < 1e-9)


def is_prime(n):
    return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))


others = [p for p in range(3, 90, 2) if is_prime(p) and p not in named]
print("other primes up to 89:", others)
vals = {p: rates[p] for p in others}
lo = min(vals, key=vals.get)
hi = max(vals, key=vals.get)
print("least %.6f at %d, largest %.6f at %d" % (vals[lo], lo, vals[hi], hi))
check("least 0.7691 at 79", lo == 79 and round(vals[lo], 4) == 0.7691)
check("largest 0.8060 at 37", hi == 37 and round(vals[hi], 4) == 0.8060)
check("band [0.769, 0.806] holds for the other primes", all(0.769 <= v <= 0.806 for v in vals.values()))
check("71, 79, 83 present in the table", all(p in rates for p in (71, 79, 83)))

mx = max(rates.values())
argmax = sorted(m for m, v in rates.items() if abs(v - mx) < 5e-7)
check("maximum 0.932039", abs(mx - 0.932039) < 5e-7)
check("attained exactly at the multiples of 7 in the table",
      argmax == sorted(m for m in rates if m % 7 == 0))
second = max(v for m, v in rates.items() if m % 7)
print("largest rate at a modulus prime to 7: %.6f at %s" % (second, sorted(m for m, v in rates.items()
                                                                           if abs(v - second) < 5e-7)))
check("no rate exceeds 0.932039 over the whole table", all(v <= 0.932039 + 5e-7 for v in rates.values()))

# theta_p/2 at the bare rates and with the margin of regions_check.py
for p, q in ((3, 0.0733), (5, 0.1157), (7, 0.0422)):
    r = rates[p]
    L = math.log(1 / r)
    th = L / (math.log(2) + 2 * L)
    r2 = r + 0.005
    L2 = math.log(1 / r2)
    th2 = L2 / (math.log(2) + 2 * L2)
    print("p=%d: theta/2 at the rate %.6f (quoted %.4f); theta at rate+0.005 = %.4f" % (p, th / 2, q, th2))
    check("theta_%d/2 supremum" % p, abs(round(th / 2, 4) - q) < 1e-9)
    check("theta_%d < 1/3" % p, th < 1 / 3)
quoted_margin = {3: 0.1425, 5: 0.2286, 7: 0.0790}
for p, q in quoted_margin.items():
    r2 = rates[p] + 0.005
    L2 = math.log(1 / r2)
    th2 = L2 / (math.log(2) + 2 * L2)
    check("theta_%d at rate + 0.005 = %.4f (quoted %.4f)" % (p, th2, q), abs(round(th2, 4) - q) < 1e-9)


# independent recomputation from the definition of T_k in (B)
def rate(m):
    inv2 = pow(2, -1, m)
    best = 0.0
    idx = lambda eta, beta: 2 * eta + beta
    for k in range(1, m):
        T = np.zeros((2 * m, 2 * m), dtype=complex)
        for eta in range(m):
            for beta in range(2):
                T[idx(eta, beta), idx((2 * eta) % m, 0)] += 1
                e = (k * (eta + beta) * inv2) % m
                T[idx(eta, beta), idx((2 * eta + 1) % m, 1)] += np.exp(2j * np.pi * e / m)
        best = max(best, max(abs(np.linalg.eigvals(T))) / 2)
    return best


for m in (3, 5, 7, 37, 71, 79, 83, 89):
    r = rate(m)
    check("recomputed rate m=%d: %.6f (table %.6f)" % (m, r, rates[m]), abs(r - rates[m]) < 2e-6)

print("ALL PASS" if ok else "SOME FAIL")
sys.exit(0 if ok else 1)
