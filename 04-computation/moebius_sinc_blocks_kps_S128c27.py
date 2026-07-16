#!/usr/bin/env python3
"""moebius_sinc_blocks_kps_S128c27.py -- kind-pasteur S128 cont.27.
THE DIVISOR-BLOCK DECOMPOSITION of the Ramanujan mean square:
  A(h) = sum_{l<=k} c_l(h) sin(theta/l) = sum_{d | h, d <= k} d * M_d(theta),
  M_d(theta) = sum_{m <= k/d} mu(m) sin(theta/(d m)),  theta = 2 pi h v lambda.
(A) verify the identity exactly (floats to 1e-12); (B) Moebius partial sums |sum mu(m)/m| <= 1
    (the small-block engine) sample check; (C) block-mass table: which d-blocks carry disc_v."""
import sys
from math import gcd, pi, sin
sys.stdout.reconfigure(line_buffering=True)
def mu(n):
    r = 1; d = 2; nn = n
    while d*d <= nn:
        if nn % d == 0:
            nn //= d
            if nn % d == 0: return 0
            r = -r
        d += 1
    if nn > 1: r = -r
    return r
def phi(n):
    r = n; d = 2; nn = n
    while d*d <= nn:
        if nn % d == 0:
            while nn % d == 0: nn //= d
            r -= r // d
        d += 1
    if nn > 1: r -= r // nn
    return r
def ram(l, h):
    g = gcd(l, h); m = l // g
    mm = mu(m)
    return 0 if mm == 0 else mm * phi(l) // phi(m)
k = 12; lam = 1.0/14
print("(A) identity A(h) = sum_{d|h} d M_d: max diff over h=1..60, v in {17, 84, 182}:")
mx = 0.0
for v in [17, 84, 182]:
    for h in range(1, 61):
        theta = 2*pi*h*v*lam
        A = sum(ram(l, h) * sin(theta/l) for l in range(1, k+1))
        B = 0.0
        for d in range(1, k+1):
            if h % d == 0:
                B += d * sum(mu(m) * sin(theta/(d*m)) for m in range(1, k//d + 1))
        mx = max(mx, abs(A - B))
print("   max |A - sum d M_d| = %.2e %s" % (mx, "OK" if mx < 1e-9 else "** FAIL"))
print("(B) Moebius partial sums |sum_{m<=M} mu(m)/m| (classical <= 1):",
      ["%.4f" % abs(sum(mu(m)/m for m in range(1, M+1))) for M in [3, 6, 12, 50, 200]])
print("(C) block mass: disc_v contribution by divisor-block d (deep well, primes + resonant v):")
for v in [17, 199, 84, 182]:
    tot = 0.0
    blocks = {}
    for h in range(1, 6000//v + 2):
        theta = 2*pi*h*v*lam
        gh = -sum(ram(l, h*v) * sin(2*pi*h*v*lam/l) for l in range(1, k+1)) / (pi*h*v)
        tot += 2*gh*gh
    # block attribution via the coprime/resonant split: report d-blocks of A(hv)
    bl = {}
    for h in range(1, 6000//v + 2):
        H = h*v
        theta = 2*pi*H*lam
        for d in range(1, k+1):
            if H % d == 0:
                Md = sum(mu(m) * sin(theta/(d*m)) for m in range(1, k//d + 1))
                bl[d] = bl.get(d, 0.0) + abs(d*Md)/(pi*H)
    top = sorted(bl.items(), key=lambda x: -x[1])[:4]
    print("   v=%3d: disc ~ %.3e ; top |block| L1-mass: %s" % (v, tot, [(d, "%.3f" % s) for d, s in top]))
print("DONE")
