#!/usr/bin/env python3
"""
fold_edge_identity_boxeph_S186.py  (HYP-8580, THM-1765 sections 1-2)

THE FOLD-EDGE IDENTITY. At a fold (critical point u_c of Lambda_s, critical
value v = Lambda(u_c)) whose degeneration (s -> 0 or s -> infinity) is
governed by a TWO-TERM Newton edge  Lambda ~ c1 s^a u^{d1} + c2 s^b u^{d2}:
    Lambda''(u_c) u_c^2  =  -d1 d2 · v · (1 + o(1)).
PROOF (2 lines, in THM-1765): u_c from d1 c1 s^a u^{d1} = -d2 c2 s^b u^{d2};
then v = c1 s^a u_c^{d1}(d2-d1)/d2 and Lambda'' u^2 = c1 s^a u_c^{d1} d1(d1-d2),
ratio = -d1 d2. CONSEQUENCE ((L1) threat defused): the far-end pair-sum of
the tracked resolvent is 2/(t Lambda'' u_c^2) ~ 2v/(Lambda'' u_c^2) ->
-2/(d1 d2): O(1) UNIVERSALLY on the infinite tube; the 1/t = v weight
cancels the u_c^2 Lambda'' degeneration exactly.

TESTS (geometry-fed):
(F1) pair model Lambda = q0 + w(p1 u + p-1/u), w = sqrt(2s): edge (d1,d2) =
     (1,-1) at s->0: ratio -> +1. Exact small-s asymptotics measured.
(F2) span-2 edge: Lambda = w^2 p2 u^2 + w p-1 /u: (d1,d2) = (2,-1):
     ratio -> +2. Also the OTHER fold family of the same Lambda.
(F3) the S183 witness Lambda = h^2 + s^2 at s -> 0: fold near sqrt(b):
     identify the edge and check the ratio against -d1 d2.
(F4) a THREE-term edge (rescaling to a fixed Laurent polynomial lambda):
     ratio = lambda''(u)u^2/lambda(u) at the rescaled critical point —
     finite nonzero; measured constant vs s.
(F5) (L2) integrand bounds: |z Lambda'(z)| for zeros z of Lambda at small s
     (the 1/sqrt(s) law: z -> finite, Lambda'(z) ~ w) and at a double-zero
     s* (the |s-s*|^{-1/2} law) — both integrable weights, measured.

boxeph-2026-07-20-S186. Pure python + cmath.
"""
import cmath
import math


def newton(f, df, x0, iters=80):
    x = x0
    for _ in range(iters):
        fx = f(x)
        d = df(x)
        if d == 0:
            break
        x = x - fx / d
    return x


print("=" * 78)
print("THE FOLD-EDGE IDENTITY  Lambda'' u_c^2 = -d1 d2 v (1+o(1))")
print("=" * 78)

# ---- (F1) pair model: d1=1, d2=-1 -> ratio +1 --------------------------------
q0c, p1, pm1 = 0.3 + 0.2j, 1.1 - 0.4j, 0.7 + 0.5j
print("\n(F1) pair edge (1,-1): Lambda = q0 s + w(p1 u + pm1/u), w = sqrt(2s)")
print("     predicted ratio -d1 d2 = +1:")
for s in (1e-2, 1e-4, 1e-6):
    w = math.sqrt(2 * s)
    lam = lambda u: q0c * s + w * (p1 * u + pm1 / u)
    d1f = lambda u: w * (p1 - pm1 / (u * u))
    d2f = lambda u: w * (2 * pm1 / u ** 3)
    uc = newton(d1f, d2f, cmath.sqrt(pm1 / p1))
    v = lam(uc)
    ratio = d2f(uc) * uc * uc / v
    print("     s=%.0e  ratio = %.6f%+.6fj" % (s, ratio.real, ratio.imag))

# ---- (F2) span-2 edge: d1=2, d2=-1 -> ratio +2 -------------------------------
p2, pm1b = 0.9 + 0.3j, 1.3 - 0.2j
print("\n(F2) edge (2,-1): Lambda = w^2 p2 u^2 + w pm1 / u")
print("     predicted ratio -d1 d2 = +2:")
for s in (1e-2, 1e-4, 1e-6):
    w = math.sqrt(2 * s)
    lam = lambda u: w * w * p2 * u * u + w * pm1b / u
    d1f = lambda u: 2 * w * w * p2 * u - w * pm1b / (u * u)
    d2f = lambda u: 2 * w * w * p2 + 2 * w * pm1b / u ** 3
    # u_c^3 = pm1b/(2 w p2): scale s^{-1/6}
    uc = newton(d1f, d2f, (pm1b / (2 * w * p2)) ** (1.0 / 3))
    v = lam(uc)
    ratio = d2f(uc) * uc * uc / v
    print("     s=%.0e  ratio = %.6f%+.6fj" % (s, ratio.real, ratio.imag))

# ---- (F3) the S183 witness at small s ----------------------------------------
A_, B_ = 0.7 + 0.3j, 0.25 + 0.10j
print("\n(F3) S183 witness Lambda = (w u - a + b w/u)^2 + s^2, fold near sqrt(b):")
print("     small-s edge of Lambda': h -> -a dominates, h' edge (1,-1) in w-terms;")
print("     Lambda'' u^2 / v at the h'-fold as s -> 0 (v -> a^2: NON-degenerate,")
print("     so the identity's regime is the DEGENERATE end only — here ratio is")
print("     the finite non-edge value; measured for the record):")
for s in (1e-2, 1e-4):
    w = math.sqrt(2 * s)
    h = lambda u: w * u - A_ + B_ * w / u
    hp = lambda u: w - B_ * w / (u * u)
    hpp = lambda u: 2 * B_ * w / u ** 3
    lam = lambda u: h(u) ** 2 + s * s
    d1f = lambda u: 2 * h(u) * hp(u)
    d2f = lambda u: 2 * (hp(u) ** 2 + h(u) * hpp(u))
    uc = newton(d1f, d2f, cmath.sqrt(B_))
    v = lam(uc)
    ratio = d2f(uc) * uc * uc / v
    print("     s=%.0e  u_c=%.4f%+.4fj  |v|=%.4f  ratio = %.4f%+.4fj" %
          (s, uc.real, uc.imag, abs(v), ratio.real, ratio.imag))

# ---- (F4) three-term edge ----------------------------------------------------
c1, c2, c3 = 1.0 + 0.2j, -0.8 + 0.5j, 0.6 - 0.3j
print("\n(F4) three-term edge (all same w-scale): Lambda = w(c1 u^2 + c2 u + c3/u):")
print("     rescaled lambda(u) = c1 u^2 + c2 u + c3/u (s-free); predicted ratio =")
print("     lambda''(u*)u*^2/lambda(u*) at lambda'(u*) = 0 — constant in s:")
lamr = lambda u: c1 * u * u + c2 * u + c3 / u
d1r = lambda u: 2 * c1 * u + c2 - c3 / (u * u)
d2r = lambda u: 2 * c1 + 2 * c3 / u ** 3
ustar = newton(d1r, d2r, 0.8 - 0.2j)
pred = d2r(ustar) * ustar * ustar / lamr(ustar)
print("     rescaled u* = %.4f%+.4fj  predicted ratio = %.4f%+.4fj" %
      (ustar.real, ustar.imag, pred.real, pred.imag))
for s in (1e-2, 1e-4, 1e-6):
    w = math.sqrt(2 * s)
    lam = lambda u: w * (c1 * u * u + c2 * u + c3 / u)
    d1f = lambda u: w * (2 * c1 * u + c2 - c3 / (u * u))
    d2f = lambda u: w * (2 * c1 + 2 * c3 / u ** 3)
    uc = newton(d1f, d2f, ustar)
    ratio = d2f(uc) * uc * uc / lam(uc)
    print("     s=%.0e  ratio = %.4f%+.4fj" % (s, ratio.real, ratio.imag))

# ---- (F5) (L2) integrand weights ---------------------------------------------
print("\n(F5) (L2) weights:")
print("  small-s zeros z of the pair model: z -> sqrt(-pm1/p1), Lambda'(z) ~ 2 p1 w:")
for s in (1e-2, 1e-4, 1e-6):
    w = math.sqrt(2 * s)
    lam = lambda u: q0c * s + w * (p1 * u + pm1 / u)
    d1f = lambda u: w * (p1 - pm1 / (u * u))
    z = newton(lam, d1f, cmath.sqrt(-pm1 / p1))
    val = abs(z * d1f(z)) / w
    print("     s=%.0e  |z Lambda'(z)|/w = %.6f  (finite => residue ~ 1/(t w): integrable 1/sqrt(s))"
          % (s, val))
print("  double-zero in u at tuned s*: Lambda has zero discriminant; residues ~ |s-s*|^{-1/2}:")
# pair model zero discriminant: q0^2 s^2 = 4 w^2 p1 pm1 -> s* solves it (real tuning)
q0r, p1r, pm1r = 1.0, 1.0, 1.0
sstar = 8.0   # q0^2 s^2 = 8 s p1 pm1 with q0=1: s* = 8
for ds in (1e-1, 1e-2, 1e-3):
    s = sstar + ds
    w = math.sqrt(2 * s)
    lam = lambda u: q0r * s / 1 + w * (p1r * u + pm1r / u)
    # zeros: w p1 u^2 + q0 s u + w pm1 = 0
    disc = cmath.sqrt((q0r * s) ** 2 - 4 * w * w * p1r * pm1r)
    z1 = (-q0r * s + disc) / (2 * w * p1r)
    d1f = lambda u: w * (p1r - pm1r / (u * u))
    val = abs(z1 * d1f(z1))
    print("     s-s*=%.0e  |z Lambda'(z)| = %.6f  (ratio to sqrt(s-s*): %.4f)"
          % (ds, val, val / math.sqrt(ds)))

print("\nDONE.")
