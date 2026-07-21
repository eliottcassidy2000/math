#!/usr/bin/env python3
"""
thm1765_referee_hijack_S186r.py  (HOSTILE REFEREE of THM-1765, boxeph-S186)

Two independent breaks of THM-1765 section 1's Consequence, plus (L2) breaks.

(R1) THE PAIR-SUM FORMULA IS WRONG.  THM-1765 claims the far-end pair-sum of
     the tracked resolvent is 2/(t Lambda'' u_c^2), hence -> -2/(d1 d2) on a
     two-term edge.  The correct collision limit of the two merging residues
     of -(1/t)/(u Lambda'(u)) includes the MIDPOINT-SHIFT term (the two roots
     are NOT symmetric about u_c at second order):
        PS = (1/t) [ 2/(u_c^2 L'') + (2/3) L'''/(u_c L''^2) ],
     and the L''' term is the SAME order as the L'' term under edge scaling.
     On a two-term edge (d1, d2) this evaluates (with t = 1/v and the
     fold-edge identity, which IS correct) to
        PS_arc = -(2/3)(d1 + d2)/(d1 d2)   [NOT -2/(d1 d2)].
     Cross-check 1 (exact): for K+ = K- = 1 (charges {+1,-1} + charge-0) the
     residue sum over BOTH roots is IDENTICALLY zero -- 1/(u(1-tLambda)) has
     no pole at u = 0 or infinity, so the two residues cancel exactly.  The
     corrected constant gives 0 for (d1,d2) = (1,-1); THM-1765's gives +2.
     Cross-check 2/3: (2,-1) -> 1/3 [claimed +1]; (3,-1) -> 4/9 [claimed 2/3].
     Measured below by direct root-finding at t = (1-eps)/v.

(R2) THE IDENTITY DOES NOT COVER VALUE-HIJACKED DEGENERATING ENDS, AND THERE
     THE PAIR-SUM BLOWS UP.  The charge-0 coefficient p_0(s) of Lambda never
     enters Lambda' or Lambda'' (u-independent) but does enter v.  Polar
     bridge (THM-1645/1680): Lambda_s(u) = P(wu, w/u), w = sqrt(2s); monomial
     Z^aW^b -> w^{a+b} u^{a-b}; p_d(0) != 0 iff P contains the PURE power
     Z^d / W^{-d}.  Support P4 = ZW + Z^9 W^7 + W:
        Lambda = 2s + a_s u^2 + b_s/u,  a_s = w^16, b_s = w.
     Fold: u_c = (b_s/2a_s)^{1/3} ~ w^{-5} (drifts to infinity); v = 2s +
     3 a_s u_c^2 ~ 2s (charge-0 HIJACKS the value; edge value ~ w^6 = 8s^3).
     On the arc 1 - tLambda = t(v - Lambda): p_0 CANCELS, the root geometry
     is pure-edge, but the weight 1/t = v ~ 2s is huge relative to it:
        exact on-arc pair-sum = v/(9 a_s u_c^2) = 2^{5/3}/72 * s^{-2} + O(1)
     -> infinity, NOT integrable against e^{-s}.  u_c^2 L'' = 6 a_s u_c^2 ~
     s^3 -> 0 (the S183r threat at rate p = 3 >= 1): the 1/t = v weight
     supplies only s^1.  The identity's prediction -d1 d2 = +2 for the ratio
     L'' u_c^2 / v fails (ratio -> 0).  (L1) is NOT closed by the identity.

(R3) The colliding pair of R2 is a MIXED pinch (one 0-group + one inf-group
     root, tracked from t ~ 0): the arc is genuine, the far tube exists.

(R4) (L2) breaks.  (a) THIRD singularity mechanism: for P1 = ZW + 2 Z^3W^2 +
     3 Z^2W^3 the zeros of Lambda_s DRIFT (z1 -> 0) and |z1 L'(z1)| ~ 2s:
     per-root weight ~ s^{-1}, not the claimed 1/sqrt(s); not integrable.
     (b) the claimed small-s law |zL'(z)| ~ w assumes p_{+-1}(0) != 0; for
     P2 = 2 Z^3W^2 + 3 Z^2W^3 it is ~ w^5 (weight s^{-5/2}).  (c) the decay
     headline is false: numerically |A(iT)| ~ T^{-2/5}, slower than the
     claimed O(T^{-1/2}) (still -> 0: Liouville's need survives).

(R5) The stacked hijack P3 = (Z^3W^2 + Z^2W^3)^2 + ZW (folds at u = +-i
     sharing v = 2s exactly) has even-part pair-sum EXACTLY zero -- by the
     u -> -u symmetry plus the residue identity, NOT by anything in THM-1765
     (whose formula would predict 2 x (+2)).  Where the file's mechanism is
     "safe" it is safe for reasons the file does not state.

boxeph-referee-2026-07-20-S186r.  Pure python + cmath.  No canon edited.
"""
import cmath
import math


def newton_poly(coeffs, z, iters=200):
    """Newton for polynomial given coeff list (highest degree first)."""
    n = len(coeffs) - 1
    for _ in range(iters):
        f = 0j
        for c in coeffs:
            f = f * z + c
        df = 0j
        for k, c in enumerate(coeffs[:-1]):
            df = df * z + (n - k) * c
        if df == 0:
            break
        zn = z - f / df
        if abs(zn - z) <= 1e-15 * max(1.0, abs(z)):
            return zn
        z = zn
    return z


def distinct_roots_near(coeffs, center, radius, want=2):
    """Newton from a ring of seeds; return `want` distinct roots near center."""
    roots = []
    for k in range(16):
        seed = center * (1 + radius * cmath.exp(2j * math.pi * k / 16))
        r = newton_poly(coeffs, seed)
        if abs(r - center) > 3 * abs(center) * radius + abs(center) * 0.5:
            continue
        if all(abs(r - r0) > 1e-7 * max(1.0, abs(center)) for r0 in roots):
            roots.append(r)
        if len(roots) == want:
            break
    return roots


print("=" * 78)
print("REFEREE S186r: pair-sum law corrected; value-hijacked ends break (L1)/(L2)")
print("=" * 78)

# ------------------------------------------------------------------ R1 ----
print("\n(R1) corrected pair-sum law  PS_arc = -(2/3)(d1+d2)/(d1 d2)")
print("     [THM-1765 claims -2/(d1 d2)].  Direct measurement, t = (1-eps)/v:")

q0c, p1c, pm1c = 0.3 + 0.2j, 1.1 - 0.4j, 0.7 + 0.5j
print("  edge (1,-1): Lambda = q0 s + w(p1 u + pm1/u):"
      "  claimed +2, corrected 0, exact-residue argument 0:")
s = 1e-4
w = math.sqrt(2 * s)
uc = cmath.sqrt(pm1c / p1c)
v = q0c * s + 2 * w * cmath.sqrt(p1c * pm1c)
for eps in (1e-4, 1e-6, 1e-8):
    t = (1 - eps) / v
    A2, B2, C2 = t * w * p1c, t * q0c * s - 1, t * w * pm1c
    D = cmath.sqrt(B2 * B2 - 4 * A2 * C2)
    r1, r2 = (-B2 + D) / (2 * A2), (-B2 - D) / (2 * A2)
    tot = 0
    for u in (r1, r2):
        lp = w * (p1c - pm1c / (u * u))
        tot += -(1 / t) / (u * lp)
    print("     eps=%.0e   PS = %.3e%+.3ej" % (eps, tot.real, tot.imag))

p2c, pm1d = 0.9 + 0.3j, 1.3 - 0.2j
print("  edge (2,-1): Lambda = w^2 p2 u^2 + w pm1/u:  claimed +1, corrected 1/3:")
s = 1e-4
w = math.sqrt(2 * s)
a_s, b_s = w * w * p2c, w * pm1d
uc = (b_s / (2 * a_s)) ** (1.0 / 3)
v = a_s * uc * uc + b_s / uc
for eps in (1e-4, 1e-6, 1e-8):
    t = (1 - eps) / v
    coeffs = [t * a_s, 0, -1, t * b_s]          # t a u^3 - u + t b = 0
    pair = distinct_roots_near(coeffs, uc, 0.05)
    tot = 0
    for u in pair:
        lp = 2 * a_s * u - b_s / (u * u)
        tot += -(1 / t) / (u * lp)
    print("     eps=%.0e   PS = %.6f%+.6fj   (pair found: %d)"
          % (eps, tot.real, tot.imag, len(pair)))

p3c = 0.8 - 0.5j
print("  edge (3,-1): Lambda = w^3 p3 u^3 + w pm1/u:  claimed +2/3, corrected 4/9:")
s = 1e-4
w = math.sqrt(2 * s)
a_s, b_s = w ** 3 * p3c, w * pm1d
uc = (b_s / (3 * a_s)) ** 0.25
v = a_s * uc ** 3 + b_s / uc
for eps in (1e-4, 1e-6, 1e-8):
    t = (1 - eps) / v
    coeffs = [t * a_s, 0, 0, -1, t * b_s]       # t a u^4 - u + t b = 0
    pair = distinct_roots_near(coeffs, uc, 0.04)
    tot = 0
    for u in pair:
        lp = 3 * a_s * u * u - b_s / (u * u)
        tot += -(1 / t) / (u * lp)
    print("     eps=%.0e   PS = %.6f%+.6fj   (pair found: %d)"
          % (eps, tot.real, tot.imag, len(pair)))
print("     => the file's far-end constant -2/(d1 d2) is FALSE; the corrected")
print("        law -(2/3)(d1+d2)/(d1d2) fits all three (0, 1/3, 4/9). O(1)-ness")
print("        for value-EDGE-GOVERNED ends survives with the corrected constant")
print("        (which VANISHES identically on the (1,-1) pair edge).")

# ------------------------------------------------------------------ R2 ----
print("\n(R2) value-hijacked end with spectator root: P4 = ZW + Z^9W^7 + W")
print("     Lambda = 2s + w^16 u^2 + w/u;  identity would predict ratio")
print("     L''u_c^2 / v -> -d1 d2 = +2;  exact pair-sum = v/(9 a_s u_c^2):")
print("     s        ratio L''u^2/v    pair-sum(exact)   *s^2      2^{5/3}/72")
for s in (1e-1, 1e-2, 1e-3, 1e-4):
    w = math.sqrt(2 * s)
    a_s, b_s = w ** 16, w
    uc = (b_s / (2 * a_s)) ** (1.0 / 3)
    v = 2 * s + 3 * a_s * uc * uc
    L2u2 = 6 * a_s * uc * uc
    ps = v / (9 * a_s * uc * uc)
    print("     %.0e   %.3e         %.6e     %.6f   %.6f"
          % (s, (L2u2 / v), ps, ps * s * s, 2 ** (5.0 / 3) / 72))
print("     numeric cross-check at s = 1e-2, t = (1-eps)/v, pair near u_c:")
s = 1e-2
w = math.sqrt(2 * s)
a_s, b_s = w ** 16, w
uc = (b_s / (2 * a_s)) ** (1.0 / 3)
v = 2 * s + 3 * a_s * uc * uc
for eps in (1e-5, 1e-7):
    t = (1 - eps) / v
    coeffs = [t * a_s, 0, t * 2 * s - 1, t * b_s]
    pair = distinct_roots_near(coeffs, uc, 0.03)
    tot = 0
    for u in pair:
        lp = 2 * a_s * u - b_s / (u * u)
        tot += -(1 / t) / (u * lp)
    print("     eps=%.0e   PS_measured = %.4f   exact v/(9 a u_c^2) = %.4f  (pair: %d)"
          % (eps, tot.real, (v / (9 * a_s * uc * uc)).real, len(pair)))
print("     => pair-sum ~ 0.0441 s^{-2} -> infinity on the far tube;")
print("        int_0 e^{-s} s^{-2} ds DIVERGES: the (L1) dominated-convergence")
print("        bound fails; u_c^2 L'' ~ s^3 -> 0 is EXACTLY the S183r threat,")
print("        and the 1/t = v ~ s weight does NOT cancel it (residual s^{-2}).")

# ------------------------------------------------------------------ R3 ----
print("\n(R3) mixedness of the R2 collision (s = 0.01): track 3 roots in tau:")
s = 1e-2
w = math.sqrt(2 * s)
a_s, b_s = w ** 16, w
uc = (b_s / (2 * a_s)) ** (1.0 / 3)
v = 2 * s + 3 * a_s * uc * uc
tau0 = 1e-4
t = tau0 / v
roots = [newton_poly([t * a_s, 0, t * 2 * s - 1, t * b_s], seed)
         for seed in (b_s * t, 1 / cmath.sqrt(a_s * t), -1 / cmath.sqrt(a_s * t))]
lab = ["0-group", "inf-group", "inf-group"]
order0 = sorted(range(3), key=lambda i: abs(roots[i]))
print("     tau=1e-4: |roots| = %s  (smallest = 0-group)"
      % ", ".join("%.3e" % abs(roots[i]) for i in order0))
tau = tau0
while tau < 1 - 1e-7:
    tau = min(tau * 1.05, 1 - 1e-7)
    t = tau / v
    coeffs = [t * a_s, 0, t * 2 * s - 1, t * b_s]
    roots = [newton_poly(coeffs, r) for r in roots]
near = sorted(range(3), key=lambda i: abs(roots[i] - uc))
print("     tau~1:  roots/u_c = %s"
      % ", ".join("%.4f%+.4fj" % ((roots[i] / uc).real, (roots[i] / uc).imag)
                  for i in range(3)))
print("     colliding pair (nearest u_c): labels {%s, %s};  spectator: %s at %.3f u_c"
      % (lab[near[0]], lab[near[1]], lab[near[2]], abs(roots[near[2]] / uc)))
mixed = {lab[near[0]], lab[near[1]]} == {"0-group", "inf-group"}
print("     => MIXED pinch: %s (genuine arc; the far tube of R2 exists)." % mixed)

# ------------------------------------------------------------------ R4 ----
gam, alp, bet = 1.0, 2.0, 3.0
print("\n(R4a) (L2) third mechanism, P1 = ZW + 2Z^3W^2 + 3Z^2W^3: zero-drift:")
for s in (1e-2, 1e-3, 1e-4):
    w = math.sqrt(2 * s)
    w5 = w ** 5
    D = cmath.sqrt(4 * gam * gam * s * s - 4 * w5 * w5 * alp * bet)
    z1 = (-2 * gam * s + D) / (2 * w5 * alp)
    z1b = (w5 * bet) / ((-2 * gam * s - D) / 2) * (-1)   # product = bet/alp
    z1 = min((z1, (w5 * bet / alp) / z1 if z1 != 0 else z1), key=abs)
    # stable small zero: z_small = (b_s)/( -p0 - sqrt(...) ) * ... use product:
    zbig = (-2 * gam * s - D) / (2 * w5 * alp)
    zsm = (w5 * bet / (w5 * alp)) / zbig * (bet / bet)   # z1 z2 = bet/alp
    zsm = (bet / alp) / zbig
    lp = w5 * (alp - bet / (zsm * zsm))
    val = abs(zsm * lp)
    print("     s=%.0e   |z_small| = %.3e   |z L'(z)|/(2s) = %.6f   weight ~ %.2e"
          % (s, abs(zsm), val / (2 * s), 1 / val))
print("     => |z L'(z)| ~ 2 gamma s (zero DRIFTS to 0): weight s^{-1}: a third")
print("        mechanism, missing from the file's two-mechanism classification;")
print("        naive per-root bound NOT integrable at s -> 0.")

print("\n(R4b) coefficient-degenerate small-s law, P2 = 2Z^3W^2 + 3Z^2W^3:")
for s in (1e-2, 1e-4, 1e-6):
    w = math.sqrt(2 * s)
    w5 = w ** 5
    z = 1j * math.sqrt(bet / alp)
    lp = w5 * (alp - bet / (z * z))
    print("     s=%.0e   |z L'(z)|/w = %.3e   /w^5 = %.4f"
          % (s, abs(z * lp) / w, abs(z * lp) / w5))
print("     => the claimed law |zL'(z)| ~ 2 p1 w fails when p_{+-1}(0) = 0:")
print("        actual ~ w^5: weight s^{-5/2}, NOT integrable.")


def A_of_t(t, alp, bet):
    N1, N2 = 600, 400
    xs = []
    for i in range(N1 + 1):
        xs.append(10 ** (-12.0 + 12.0 * i / N1))
    for i in range(1, N2 + 1):
        xs.append(1.0 + 39.0 * i / N2)
    vals = []
    for s in xs:
        w5 = (2 * s) ** 2.5
        q = t * w5
        disc = cmath.sqrt(1 - 4 * q * q * alp * bet)
        u1 = 2 * q * bet / (1 + disc)
        lp = w5 * (alp - bet / (u1 * u1))
        vals.append(cmath.exp(-s) * (-(1 / t) / (u1 * lp)))
    tot = 0
    for i in range(len(xs) - 1):
        tot += (vals[i] + vals[i + 1]) / 2 * (xs[i + 1] - xs[i])
    tot += xs[0] * 1.0
    return tot


print("\n(R4c) ray-decay rate for P2 along t = iT (arcs sit on the real axis):")
print("     T          |A|           |A| T^{2/5}    |A| T^{1/2}")
for T in (1e2, 1e4, 1e6, 1e8):
    a = A_of_t(1j * T, alp, bet)
    print("     %.0e   %.4e    %10.5f    %10.5f"
          % (T, abs(a), abs(a) * T ** 0.4, abs(a) * T ** 0.5))
print("     => |A| = Theta(T^{-2/5}): the section-2 headline O(|t|^{-1/2}) is")
print("        FALSE as stated (decay survives, slower; Liouville unaffected).")

# ------------------------------------------------------------------ R5 ----
print("\n(R5) stacked hijack P3 = (Z^3W^2+Z^2W^3)^2 + ZW: Lambda = 2s + w^10(u+1/u)^2")
print("     folds u = +-i share v = 2s exactly; measure the 4-root residue sum")
print("     and each fold's 2-root sum at t = (1-eps)/v (s = 0.01):")
s = 1e-2
w10 = (2 * s) ** 5
v = 2 * s
for eps in (1e-5, 1e-7):
    t = (1 - eps) / v
    # roots: w10 t (u^2+1/u^2+2) + 2st = 1 -> quartic t w10 u^4 + (2s t - 1 + 2 t w10) u^2 + t w10 = 0
    Aq, Bq, Cq = t * w10, 2 * s * t - 1 + 2 * t * w10, t * w10
    Dq = cmath.sqrt(Bq * Bq - 4 * Aq * Cq)
    y1, y2 = (-Bq + Dq) / (2 * Aq), (-Bq - Dq) / (2 * Aq)
    us = [cmath.sqrt(y1), -cmath.sqrt(y1), cmath.sqrt(y2), -cmath.sqrt(y2)]
    tots = {"+i": 0, "-i": 0}
    for u in us:
        lp = w10 * 2 * (u + 1 / u) * (1 - 1 / (u * u))
        r = -(1 / t) / (u * lp)
        if abs(u - 1j) < 0.5:
            tots["+i"] += r
        elif abs(u + 1j) < 0.5:
            tots["-i"] += r
    print("     eps=%.0e  fold(+i) 2-root sum = %.3e%+.3ej ; fold(-i) = %.3e%+.3ej"
          % (eps, tots["+i"].real, tots["+i"].imag,
             tots["-i"].real, tots["-i"].imag))
print("     => each fold's pair-sum is ~ 0 (u -> -u symmetry + residue identity),")
print("        while THM-1765's formula would predict +2 per fold.  The safety of")
print("        this stacked hijack is NOT the fold-edge identity's doing.")

print("\nNET: section 1's Consequence is broken twice over -- wrong constant on")
print("edge-governed ends (R1), and unbounded (non-integrable) pair-sums on")
print("value-hijacked ends with spectator roots (R2, mixed per R3).  Section 2's")
print("two-mechanism classification and O(t^{-1/2}) are false as stated (R4).")
print("DONE.")
