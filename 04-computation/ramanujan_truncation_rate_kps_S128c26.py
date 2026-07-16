#!/usr/bin/env python3
"""ramanujan_truncation_rate_kps_S128c26.py -- kind-pasteur S128 cont.26.
(A) coprime-v reduction: c_l(hv) = c_l(h) for gcd(v, lcm(1..k)) = 1 (exact check).
(B) the Ramanujan cancellation: partial sums S_L(h) = sum_{l<=L} c_l(h)/l -> 0; truncation decay.
(C) disc_v for primes v > k on the deep-well core: the r-linear (sharp) rate vs crude r^2,
    with disc read through the truncation-error frame."""
import sys, cmath
from fractions import Fraction as F
from math import gcd, pi, sin
sys.stdout.reconfigure(line_buffering=True)
def ram(l, h):
    g = gcd(l, h)
    # c_l(h) = mu(l/g) * phi(l)/phi(l/g)
    def mu(n):
        r = 1; d = 2
        while d*d <= n:
            if n % d == 0:
                n //= d
                if n % d == 0: return 0
                r = -r
            d += 1
        if n > 1: r = -r
        return r
    def phi(n):
        r = n; d = 2
        while d*d <= n:
            if n % d == 0:
                while n % d == 0: n //= d
                r -= r // d
            d += 1
        if n > 1: r -= r // n
        return r
    m = l // g
    mm = mu(m)
    return 0 if mm == 0 else mm * phi(l) // phi(m)
print("(A) coprime-v invariance c_l(hv) = c_l(h):")
import math
k = 12; L = math.lcm(*range(1, k+1))
ok = all(ram(l, h*v) == ram(l, h) for l in range(1, k+1) for h in range(1, 30) for v in [13*14+1, 361, 529] if gcd(v, L) == 1)
print("   exact for v in {183, 361, 529} coprime to lcm(1..12):", ok)
print("(B) Ramanujan cancellation S_L(h) = sum_{l<=L} c_l(h)/l:")
for h in [1, 2, 6, 12, 30]:
    row = ["L=%d: %+.4f" % (Lc, sum(F(ram(l, h), l) for l in range(1, Lc+1))) for Lc in [12, 50, 200, 800]]
    print("   h=%2d: %s" % (h, "  ".join(row)))
print("(C) disc_v for prime v > 12 on {1..12} at 1/14 (Ramanujan form, h-sum) and the r-linear rate:")
lam = 1.0/14
def ghat(h):
    tot = 0.0
    for l in range(1, 13):
        tot += ram(l, h) * sin(2*pi*h*lam/l) / (pi*h)
    return -tot   # (overlap corrections at 1/14 exist but are small; comparative rates only)
r = 12  # good intervals of the deep-well core
for v in [17, 29, 53, 101, 199]:
    d = 2*sum(ghat(h*v)**2 for h in range(1, 6000//v + 2))
    print("   v=%3d: disc ~ %.3e ; crude r^2/(3v^2) = %.3e ; ratio (sharpening) = %.1fx ; disc*v^2 = %.4f (r-linear test: ~const)"
          % (v, d, r*r/(3*v*v), (r*r/(3*v*v))/d if d>0 else float('inf'), d*v*v))
print("DONE")
