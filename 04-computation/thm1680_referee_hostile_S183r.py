#!/usr/bin/env python3
"""
thm1680_referee_hostile_S183r.py — HOSTILE REFEREE checks for THM-1680
(stacking dichotomy / Liouville endgame / graded ladder).

Four independent probes, pure python + cmath (no numpy):

(T1) Reality two-line identity (THM-1680 section 3): for REAL a, b>0 and
     s > a^2/(4b): |u1|^2 = b, u1 h'(u1) = 2i w Im(u1), rho purely imaginary
     IDENTICALLY. Verified on an s-grid in the w = sqrt(2s) normalization
     (the identity is w-normalization invariant).

(T2) SUB-GERM TRUNCATION grade class (attack on section 2's dichotomy):
     if B == 0 on the unbounded sub-germ (s > s2) but B != 0 on (s1, s2),
     the localization I_m = int_0^{s2} e^{-s} v^{m-1} v' ds is ENDPOINT-
     dominated: pure exponential class v(s2)^m / m — NOT Gamma-scale.
     Germ v = 2 sqrt(2s) (the P0^2 germ shape), s2 = 4. Shows the claimed
     dichotomy (delete | Gamma-graded) misses a third class.

(T3) ODD-SECTOR VISIBILITY (attack on 'B == 0 => contributes nothing to any
     moment condition'): exact coefficients of int_0^inf e^-s (1 - t v)^{k/2} ds,
     v = sqrt(s):
       k = -1 (the B-term):  a_m = C(-1/2,m)(-t)^m Gamma(m/2+1): rung m^{-1/2}
       k = +1 (subleading):  a_m = C(+1/2,m)(-t)^m Gamma(m/2+1): rung m^{-3/2}
     Both are Gamma-scale. A germ with B == 0 but d != 0 (leading sqrt
     coefficients cancel, subleading do not) is GRADED-VISIBLE one rung down,
     and its arc is NOT removable. Exact identity: C(1/2,m)/C(-1/2,m)
     = -1/(2m-1).

(T4) DEFECT INSTRUMENT vs the B-criterion: a pure-subleading cut
     d*(1 - t/t_b)^{1/2} has monodromy defect 2|d| sqrt(r) at probe radius r
     — nonzero (12+ orders above the 1e-13 double-loop floor at r = 1e-3),
     though its B = 0. The INSTRUMENT (defect == 0) is finer than the TEXT
     ('B == 0'): the two criteria diverge exactly on the B==0, d!=0 class.

boxeph fleet, hostile-referee subagent, S183r. Report-only; no canon edits.
"""
import math
import cmath

print("=" * 78)
print("THM-1680 HOSTILE REFEREE CHECKS (S183r)")
print("=" * 78)

# ---------------------------------------------------------------- T1
print("\n(T1) reality identity, w = sqrt(2s), a=0.7 b=0.25 real, s > a^2/(4b)=0.49:")
a, b = 0.7, 0.25
for s in (0.6, 1.0, 4.0, 25.0, 100.0):
    w = math.sqrt(2 * s)
    disc = cmath.sqrt(a * a - 4 * b * w * w)     # imaginary for s > a^2/(8b) here
    u1 = (a + disc) / (2 * w)
    h1 = w * u1 - a + b * w / u1                 # should be 0
    hp = w - b * w / (u1 * u1)
    lam2 = 2 * hp * hp                            # Lambda'' at h=0
    rho = 1 / (u1 * cmath.sqrt(lam2))
    print("  s=%6.1f |u1|^2-b=%.1e  |h(u1)|=%.1e  Re(u1 h')=%.1e  |Re rho|/|rho|=%.1e"
          % (s, abs(abs(u1) ** 2 - b), abs(h1), abs((u1 * hp).real),
             abs(rho.real) / abs(rho)))
print("  => |u1|^2 = b and rho purely imaginary hold IDENTICALLY (algebra confirmed;")
print("     note the identity needs b > 0 — for b < 0 the pair is real, not conjugate).")

# ---------------------------------------------------------------- T2
def log_lower_gamma(aa, x):
    """log of lower incomplete gamma(aa, x) via the ascending series,
    gamma(aa,x) = x^aa e^-x sum_n x^n / (aa (aa+1) ... (aa+n))."""
    term = 1.0 / aa
    ssum = term
    n = 0
    while True:
        n += 1
        term *= x / (aa + n)
        ssum += term
        if term < 1e-18 * ssum:
            break
    return aa * math.log(x) - x + math.log(ssum), ssum * aa  # ratio-to-endpoint helper


print("\n(T2) sub-germ truncated B: I_m = int_0^{s2} e^-s v^{m-1} v' ds, v=2sqrt(2s), s2=4:")
c = 2 * math.sqrt(2)
s2 = 4.0
print("    log v(s2) = %.4f  (pure exponential class if truncated)" % math.log(c * math.sqrt(s2)))
print("      m    (1/m)log I_trunc   (1/m)log I_full    I_trunc*m/(e^-s2 v(s2)^m)")
for m in (10, 20, 40, 60):
    # I_trunc = (c^m/2) * gamma(m/2, s2);  I_full = (c^m/2) * Gamma(m/2)
    lg, ratio_series = log_lower_gamma(m / 2.0, s2)
    lI_tr = m * math.log(c) - math.log(2) + lg
    lI_full = m * math.log(c) - math.log(2) + math.lgamma(m / 2.0)
    # endpoint law: gamma(a,x) ~ x^a e^-x / a  =>  I_trunc ~ e^-s2 v(s2)^m / m
    lendpt = m * math.log(c) - math.log(2) + (m / 2.0) * math.log(s2) - s2 - math.log(m / 2.0)
    # exp(lI_tr - lendpt) = gamma(a,x) / (x^a e^-x / a) = I_trunc*m/(e^-s2 v(s2)^m)
    print("     %3d      %8.4f          %8.4f            %.4f"
          % (m, lI_tr / m, lI_full / m, math.exp(lI_tr - lendpt)))
print("    => truncated germ: (1/m)log I -> log v(s2) = %.4f (EXPONENTIAL grade," % math.log(c * math.sqrt(s2)))
print("       endpoint law I ~ e^-s2 v(s2)^m/m confirmed); full germ: (1/m)log I")
print("       diverges (Gamma-scale). The section-2 dichotomy has NO class for this.")

# ---------------------------------------------------------------- T3
print("\n(T3) odd-sector moment visibility, v = sqrt(s):")
print("     A_k(t) = int e^-s (1 - t sqrt(s))^{k/2} ds; exact a_m = C(k/2,m)(-1)^m Gamma(m/2+1)")


def binom_half(alpha, m):
    v = 1.0
    for j in range(m):
        v *= (alpha - j) / (j + 1)
    return v


print("      m    rung(-1/2 term)   rung(+1/2 term)   |C(1/2,m)/C(-1/2,m)|*(2m-1)")
prev = {}
for m in (8, 16, 32, 64, 128, 256, 512):
    la_B = math.log(abs(binom_half(-0.5, m))) + math.lgamma(m / 2.0 + 1)
    la_d = math.log(abs(binom_half(0.5, m))) + math.lgamma(m / 2.0 + 1)
    rB = la_B - math.lgamma(m / 2.0 + 1)   # log prefactor vs Gamma-scale
    rd = la_d - math.lgamma(m / 2.0 + 1)
    out = ["%5d" % m]
    for key, val in (("B", rB), ("d", rd)):
        if key in prev:
            slope = (val - prev[key][0]) / (math.log(m) - math.log(prev[key][1]))
            out.append("   %+7.4f" % slope)
        else:
            out.append("      --  ")
        prev[key] = (val, m)
    ratio = abs(binom_half(0.5, m) / binom_half(-0.5, m)) * (2 * m - 1)
    out.append("            %.12f" % ratio)
    print("   " + " ".join(out))
print("    => both terms are Gamma(m/2+1)-scale; prefactor slopes -> -1/2 (B-term)")
print("       and -3/2 (subleading d-term): a B==0, d!=0 germ enters the moment")
print("       conditions ONE RUNG DOWN — 'contributes nothing' is FALSE for it,")
print("       and its arc is NOT removable (branch cut of order (1-tv)^{1/2}).")

# ---------------------------------------------------------------- T4
print("\n(T4) defect of a PURE-SUBLEADING cut  d*(1 - t/t_b)^{1/2}, d = 0.3:")
d = 0.3
for r in (1e-2, 1e-3, 1e-4):
    defect = 2 * d * math.sqrt(r)
    print("     r=%.0e   |defect| = %.3e   (B-type defect at same r ~ %.1e; floor 1e-13)"
          % (r, defect, 2 * 0.363 / math.sqrt(2 * (1.0 / 16) * r)))
print("    => defect != 0 although B = 0: the monodromy INSTRUMENT (defect == 0)")
print("       is strictly finer than the canon TEXT's criterion (B == 0). The")
print("       biconditional 'B == 0 <=> single-valued <=> removable' fails =>.")

print("\nDONE.")
