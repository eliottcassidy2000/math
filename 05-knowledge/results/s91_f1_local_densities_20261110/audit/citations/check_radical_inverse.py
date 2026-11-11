#!/usr/bin/env python3
"""Auditor 'citations' (S91): the radical-inverse / Hammersley statements of THM-5128 (B3) and README section 6.
   rho(h) = x(h)/y(h),  phi(h) = sum_i h_i 2^(-i-1)  (van der Corput radical inverse in base 2),
   claim: |rho(h) - (1 - phi(h)/2)| < n/(2h + 1) for h with exactly n binary digits;
   the points (h / 2^n, phi(h)), h < 2^n, form the Hammersley set: every dyadic box of area 2^-n holds one point.
Exact rational arithmetic; no code of the audited session is imported."""
import sys
from fractions import Fraction

OK = True


def check(name, cond, detail=""):
    global OK
    OK = OK and bool(cond)
    print(("PASS  " if cond else "FAIL  ") + name + ("   [" + str(detail) + "]" if detail != "" else ""))


NMAX = 16
psi = [0] * (1 << NMAX)
for h in range(1, 1 << NMAX):
    k = h >> 1
    psi[h] = psi[k] + ((k >> 1) + 1 if h & 1 else k)


def phi(h):
    num = 0
    n = h.bit_length()
    for i in range(n):
        if h >> i & 1:
            num += 1 << (n - i - 1)
    return Fraction(num, 1 << n)


worst = Fraction(0)
bad = 0
lo, hi = Fraction(2), Fraction(0)
for h in range(1, 1 << NMAX):
    n = h.bit_length()
    x, y = 2 * psi[h] + 1, 2 * h + 1
    rho = Fraction(x, y)
    lo, hi = min(lo, rho), max(hi, rho)
    d = abs(rho - (1 - phi(h) / 2))
    bound = Fraction(n, y)
    if d >= bound:
        bad += 1
    worst = max(worst, d / bound)
check("(B3): |rho(h) - (1 - phi(h)/2)| < n/(2h+1) for 1 <= h < 2^16", bad == 0, "largest ratio to the bound %.4f" % float(worst))
check("(B3): y/2 < x <= y, i.e. 1/2 < rho <= 1, for 1 <= h < 2^16", Fraction(1, 2) < lo and hi <= 1, (float(lo), float(hi)))
# Hammersley net property: for every n and j, each box [a/2^j, (a+1)/2^j) x [b/2^(n-j), (b+1)/2^(n-j)) holds one point
def net_ok(n):
    pts = []
    for h in range(1 << n):
        num = 0
        for i in range(n):
            if h >> i & 1:
                num += 1 << (n - i - 1)
        pts.append((h, num))  # (h/2^n, num/2^n)
    for j in range(n + 1):
        seen = set()
        for (u, v) in pts:
            key = (u >> (n - j), v >> j)
            if key in seen:
                return False
            seen.add(key)
        if len(seen) != 1 << n:
            return False
    return True


check("the set {(h/2^n, phi_n(h))} is a (0, n, 2)-net in base 2 (Hammersley / van der Corput), n = 1..12", all(net_ok(n) for n in range(1, 13)))
# uniformity of rho on [1/2, 1]: empirical distribution function at n = 16
import bisect
vals = sorted(float(Fraction(2 * psi[h] + 1, 2 * h + 1)) for h in range(1 << 15, 1 << 16))
dev = max(abs(bisect.bisect_right(vals, 0.5 + t / 200.0) / len(vals) - t / 100.0) for t in range(1, 100))
check("rho is close to uniform on [1/2, 1] on B_16 (largest deviation of the distribution function < 0.01)", dev < 0.01, "%.5f" % dev)
print("ALL CHECKS PASSED" if OK else "SOME CHECKS FAILED")
sys.exit(0 if OK else 1)
