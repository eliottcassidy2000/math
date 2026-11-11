#!/usr/bin/env python3
"""Auditor numbers-vs-outputs: read rates_exact/scan_003_101.npz (committed) and print the per-modulus maxima,
to check HYP-9660 (B), THM-5128 (B) ("other primes up to 89") and README section 3 ("followed by ... 31, 15, 63, 127, 73").
Read-only on the committed file; path passed as argument."""
import sys
from math import gcd

import numpy as np

z = np.load(sys.argv[1])
m, k, sign, est, rate = z["m"], z["k"], z["sign"], z["est"], z["rate"]
print("arrays:", {n: z[n].shape for n in z.files})
mods = sorted(set(m.tolist()))
print("moduli:", len(mods), "from", mods[0], "to", mods[-1], "all odd:", all(x % 2 == 1 for x in mods),
      "complete range:", mods == list(range(3, 102, 2)))
print("nan rates:", int(np.isnan(rate).sum()))
best = {}
bestprim = {}
for mm in mods:
    sel = m == mm
    r = rate[sel]
    kk = k[sel]
    best[mm] = float(np.nanmax(r))
    prim = np.array([gcd(int(a), mm) == 1 for a in kk])
    bestprim[mm] = float(np.nanmax(r[prim])) if prim.any() else float("nan")
print("per modulus: m, max rate over all k != 0, max over k coprime to m")
for mm in mods:
    print(f"  {mm:4d}  {best[mm]:.6f}  {bestprim[mm]:.6f}")
gmax = max(best.values())
print("global max: %.6f at m = %s" % (gmax, [mm for mm in mods if abs(best[mm] - gmax) < 5e-7]))
print("multiples of 7 in range:", [mm for mm in mods if mm % 7 == 0])
np7 = {mm: v for mm, v in best.items() if mm % 7}
top = sorted(np7.items(), key=lambda t: -t[1])[:12]
print("largest rates at moduli prime to 7:", [(mm, round(v, 6)) for mm, v in top])
ranking = sorted(best.items(), key=lambda t: -t[1])
print("all moduli ranked by rate (top 25):", [(mm, round(v, 5)) for mm, v in ranking[:25]])
primes = [p for p in mods if all(p % d for d in range(3, int(p ** 0.5) + 1, 2))]
print("prime moduli and rates:", [(p, round(best[p], 6)) for p in primes])
named = {3, 5, 7, 11, 13, 17, 19, 23, 31, 43, 73}
others = [(p, best[p]) for p in primes if p <= 89 and p not in named]
print("other primes up to 89:", [(p, round(v, 6)) for p, v in others])
print("  min %.6f  max %.6f" % (min(v for _, v in others), max(v for _, v in others)))
# composite rule: rate(m) vs max over prime-power divisors (prime powers p^a exactly dividing or dividing m)
def prime_power_divisors(n):
    out = []
    for d in range(3, n + 1, 2):
        if n % d == 0:
            # d is a prime power?
            x = d
            p = next(q for q in range(3, d + 1, 2) if d % q == 0)
            while x % p == 0:
                x //= p
            if x == 1:
                out.append(d)
    return out
print("composite moduli whose rate differs from the max over their prime-power divisors (tolerance 5e-5):")
viol = []
for mm in mods:
    pp = prime_power_divisors(mm)
    if pp == [mm] or mm in pp:
        continue
    ref = max(best[d] for d in pp)
    if abs(best[mm] - ref) > 5e-5:
        viol.append(mm)
        print(f"  m = {mm}: rate {best[mm]:.6f}; prime-power divisors {[(d, round(best[d], 4)) for d in pp]} (max {ref:.4f})")
print("list:", viol)
# proper-divisor rule: rate(m) vs max over ALL proper divisors
viol2 = []
for mm in mods:
    divs = [d for d in mods if d < mm and mm % d == 0]
    if not divs:
        continue
    ref = max(best[d] for d in divs)
    if abs(best[mm] - ref) > 5e-5:
        viol2.append((mm, round(best[mm], 4), round(ref, 4)))
print("moduli whose rate differs from the max over all proper divisors:", viol2)
for mm in (15, 45, 75, 51, 65, 85, 93, 63, 21, 35, 49, 77, 91, 31, 73, 9, 27, 81, 25):
    print(f"  m = {mm}: {best[mm]:.6f} (primitive characters: {bestprim[mm]:.6f})")
