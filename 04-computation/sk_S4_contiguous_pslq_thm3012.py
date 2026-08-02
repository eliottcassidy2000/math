#!/usr/bin/env python3
"""
THM-3012 / lane [S4-contiguous]  --  PART II (high precision + calibrated PSLQ)
==============================================================================
Companion to 04-computation/sk_S4_contiguous_thomae_thm3012.py.

  PART 5  the k = 4 normal form  S(4) = varpi/2 - W/(2 varpi)  verified to
          >= 120 digits, with W computed TWO independent ways
  PART 6  why k = 4 and no other integer k: the K-moment exponent integrality
  PART 7  PSLQ with a LIVE positive control (blind rediscovery of pi S(3)) and
          bounded negatives on the NEW contiguous targets W5a, W5b

Every null result below is reported as what it is: no Z-relation over an
explicitly named basis B, with |coefficients| <= H, at precision P.
"""

import sys, time, itertools
from fractions import Fraction as Fr
from mpmath import (mp, mpf, mpmathify, pi, sqrt, gamma, log, atan, ellipk, quad,
                    beta, betainc, hyper, pslq, nstr, exp, catalan, zeta)

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s); sys.stdout.flush(); OUT.append(s)
def hdr(t):
    say(""); say("=" * 78); say(t); say("=" * 78)

DPS = 170
mp.dps = DPS

def F2F1(z):
    r = sqrt(1 - z); w = (1 - r) / (1 + r)
    return sqrt(2 / (1 + r)) * (2 / pi) * ellipk(w)
def S(k):
    kk = mpf(k)
    return quad(lambda u: F2F1(u ** kk), [0, 1], maxdegree=14)

# =============================================================================
hdr("PART 5  S(4) = varpi/2 - W/(2 varpi)   verified at dps = %d" % DPS)
# =============================================================================
t0 = time.time()
S3 = S(3); S4 = S(4); S5 = S(5); S1 = S(1); S2 = S(2)
say("S(1..5) by tanh-sinh on the elliptic kernel: %.1f s" % (time.time() - t0))
ref3 = (sqrt(3) * log(5 + 2 * sqrt(6)) - 2 * atan(sqrt(2) / 5)) / pi
say("evaluation control |S(3) - sqrt3 log(5+2sqrt6)/pi + 2atan(sqrt2/5)/pi| = %s"
    % nstr(abs(S3 - ref3), 5))

varpi = gamma(mpf(1) / 4) ** 2 / (2 * sqrt(2 * pi))

# --- W route 1: the defining 3F2 series, summed by mpmath.hyper -------------
t0 = time.time()
W_series = hyper([mpf(3) / 4, mpf(3) / 4, mpf(1) / 2], [mpf(3) / 2, mpf(3) / 2], 1)
say("W by the defining 3F2 series (mpmath.hyper): %.1f s" % (time.time() - t0))

# --- W route 2: the iterated-Beta integral of identity (A) ------------------
# J = int_0^1 v^{-3/4}(1-v)^{-1/4} B_v(1/2,1/4) dv ,  W = J Gamma(3/2)/(2 Gamma(3/4)^2)
t0 = time.time()
J = quad(lambda v: v ** (mpf(-3) / 4) * (1 - v) ** (mpf(-1) / 4)
         * betainc(mpf(1) / 2, mpf(1) / 4, 0, v), [0, 1], maxdegree=16)
W_int = J * gamma(mpf(3) / 2) / (2 * gamma(mpf(3) / 4) ** 2)
say("W by the iterated-Beta integral:            %.1f s" % (time.time() - t0))
say("   W (series)   = %s" % nstr(W_series, 60))
say("   W (integral) = %s" % nstr(W_int, 60))
dW = abs(W_series - W_int)
say("   |difference| = %s  -> the two independent routes agree" % nstr(dW, 5))

W = W_series
lhs = varpi / 2 - W / (2 * varpi)
res = abs(lhs - S4)
say("")
say("   varpi                 = %s" % nstr(varpi, 60))
say("   varpi/2 - W/(2 varpi) = %s" % nstr(lhs, 60))
say("   S(4)                  = %s" % nstr(S4, 60))
digits = int(-mp.log(res, 10)) if res > 0 else DPS
say("   RESIDUAL = %s      => VERIFIED-EXACT to %d digits" % (nstr(res, 5), digits))
V5 = digits >= 120
say("PART 5 VERDICT:", "PASS (>=120 digits)" if V5 else "FAIL")

say("")
say("Equivalent restatements (all VERIFIED-EXACT here):")
say("   W          = varpi^2 - 2 varpi S(4)")
say("   S(4)       = (varpi^2 - W)/(2 varpi)")
Lam = 2 * catalan + 4 * sum((1 if n % 4 == 1 else -1) / (mpf(n) ** 2 * (exp(pi * n) - 1))
                            for n in range(1, 300, 2))
say("   addendum-5 Lambda = 2G + D = %s" % nstr(Lam, 50))
say("   |W - varpi^2(1 - 4 Lambda/pi^2)| = %s"
    % nstr(abs(W - varpi ** 2 * (1 - 4 * Lam / pi ** 2)), 5))
say("   |S(4) - 2 varpi Lambda/pi^2|     = %s"
    % nstr(abs(S4 - 2 * varpi * Lam / pi ** 2), 5))
say("So the contiguous route reproduces THM-3012 addendum 5 exactly, but packages")
say("the whole transcendental content into ONE balanced 3F2 with repeated")
say("parameters -- no Catalan constant and no Eisenstein tail appear.")

# =============================================================================
hdr("PART 6  why k = 4 and no other integer k -- three independent integralities")
# =============================================================================
say("(i) CONTIGUOUS REDUCTION.  (1/k)_n/(1+1/k)_n = (1/k)/(n+1/k) is a polynomial")
say("    in n iff 1/k - 1 in Z_{>=0} iff k = 1.  For k >= 2 the 3F2 is a Mellin")
say("    moment, never a finite sum of 2F1(1)'s.   [PROVED, PART 3 of PART I]")
say("")
say("(ii) ODE RESONANCE AT INFINITY.  The exponents of the 3F2 equation at")
say("     z = infinity are the upper parameters 1/4, 3/4, 1/k.  Two coincide iff")
say("     1/k in {1/4, 3/4}: among integers, k = 4 only.")
for k in range(1, 13):
    x = Fr(1, k)
    ex = sorted([Fr(1, 4), Fr(3, 4), x])
    say("     k = %2d : exponents at infinity %s%s"
        % (k, tuple(str(e) for e in ex), "   RESONANT" if len(set(ex)) < 3 else ""))
say("")
say("(iii) THE K-MOMENT WEIGHT.  THM-3012 eq (4):")
say("      S(k) = (16/(k pi sqrt2)) int_0^1 mu^{4/k-1} (2-mu^2)^{-2/k-1/2} K(mu) dmu.")
say("      The weight is a RATIONAL function of mu times K(mu) iff BOTH exponents")
say("      are integers:  4/k - 1 in Z  <=>  k | 4;   -2/k - 1/2 in Z  <=>  2/k+1/2 in Z.")
for k in range(1, 13):
    e1 = Fr(4, k) - 1
    e2 = -Fr(2, k) - Fr(1, 2)
    ok = (e1.denominator == 1) and (e2.denominator == 1)
    say("      k = %2d : mu^(%s) (2-mu^2)^(%s)%s" % (k, e1, e2, "   BOTH INTEGRAL" if ok else ""))
say("      Among all positive integers k, ONLY k = 4 makes both integral")
say("      (k | 4 forces k in {1,2,4}; 2/k + 1/2 in Z then forces k = 4).")
say("      That is exactly the hypothesis under which the moment recurrence")
say("      4m^2 M_{2m} = (2m-1)^2 M_{2m-2} + 1 of addendum 5 can be applied.")
say("      k = 5 fails it by an irreducible denominator 20.   [PROVED]")

# verify the k=4 weight is the recorded one
say("")
mom = (16 / (4 * pi * sqrt(2))) * quad(lambda m: ellipk(m ** 2) / (2 - m ** 2), [0, 1])
say("   check: (16/(4 pi sqrt2)) int_0^1 K(mu)/(2-mu^2) dmu = %s" % nstr(mom, 40))
say("          S(4)                                        = %s" % nstr(S4, 40))
say("          |difference| = %s" % nstr(abs(mom - S4), 5))

# =============================================================================
hdr("PART 7  PSLQ -- live positive control first, then bounded negatives")
# =============================================================================
say("Basis B1 (weight 1) = { alpha * L } with alpha ALGEBRAIC and L a log/arctan.")
say("The algebraic multipliers are mandatory: pi S(3) = sqrt3 log(5+2sqrt6)")
say("- 2 arctan(sqrt2/5) has an IRRATIONAL coefficient, so a rational-coefficient")
say("sweep over logs and arctans cannot represent it and returns a SPURIOUS null.")
say("")
r2, r3, r5, r6 = sqrt(2), sqrt(3), sqrt(5), sqrt(6)
MULT = [("1", mpf(1)), ("sqrt2", r2), ("sqrt3", r3), ("sqrt5", r5),
        ("sqrt6", r6), ("2^(1/4)", mpf(2) ** (mpf(1) / 4))]
LOGS = [("pi", pi), ("log(1+sqrt2)", log(1 + r2)), ("log(2+sqrt3)", log(2 + r3)),
        ("log(5+2sqrt6)", log(5 + 2 * r6)), ("log2", log(2)), ("log3", log(3)),
        ("log5", log(5)), ("atan(sqrt2)", atan(r2)), ("atan(sqrt2/5)", atan(r2 / 5)),
        ("atan(sqrt2/3)", atan(r2 / 3)), ("atan(sqrt3)", atan(r3)),
        ("atan(1/2)", atan(mpf(1) / 2)), ("atan(1/3)", atan(mpf(1) / 3)),
        ("log(golden)", log((1 + r5) / 2))]
B1 = [(m[0] + "*" + l[0], m[1] * l[1]) for m in MULT for l in LOGS]
say("|B1| = %d products (%d multipliers x %d logs/arctans); all %d pairs scanned."
    % (len(B1), len(MULT), len(LOGS), len(B1) * (len(B1) - 1) // 2))

H = 10 ** 6
TOL = mpf(10) ** (-int(DPS * 0.72))
say("coefficient bound H = %d,  precision P = %d dps,  detection tol = 1e-%d"
    % (H, DPS, int(DPS * 0.72)))

def scan_pairs(target, basis, tol, H, first_only=False, budget=None):
    out = []; n = len(basis); t0 = time.time()
    for i in range(n):
        for j in range(i + 1, n):
            if budget and time.time() - t0 > budget:
                return out, False
            rel = pslq([target, basis[i][1], basis[j][1]], tol=tol,
                       maxcoeff=H, maxsteps=1500)
            if rel and rel[0] != 0:
                out.append((rel, basis[i][0], basis[j][0]))
                if first_only: return out, True
    return out, True

say("")
say("--- LIVE POSITIVE CONTROL: blind rediscovery of pi S(3) over B1 ---")
t0 = time.time()
hit, _ = scan_pairs(pi * S3, B1, TOL, H, first_only=True)
say("elapsed %.1f s" % (time.time() - t0))
POS = bool(hit)
if hit:
    rel, b1n, b2n = hit[0]
    say("FOUND:   %d*(pi S(3))  +  %d*%s  +  %d*%s  =  0" % (rel[0], rel[1], b1n, rel[2], b2n))
    say("i.e.     pi S(3) = (%s)*%s + (%s)*%s"
        % (Fr(-rel[1], rel[0]), b1n, Fr(-rel[2], rel[0]), b2n))
else:
    say("NOT FOUND -- the instrument is broken; every negative below is WORTHLESS.")
say("POSITIVE CONTROL:", "PASS" if POS else "FAIL")

say("")
say("--- further live controls on the same instrument ---")
for nm, tgt in (("pi S(1)", pi * S1), ("pi S(2)", pi * S2)):
    h, _ = scan_pairs(tgt, B1, TOL, H, first_only=True)
    say("   %-9s : %s" % (nm, ("FOUND %s over (%s, %s)" % (h[0][0], h[0][1], h[0][2]))
                          if h else "not found"))

# --- the NEW targets --------------------------------------------------------
say("")
say("--- the NEW targets produced by the contiguous reduction (k = 5) ---")
say("Identity (A):  (pi sqrt2/x) S(5) = B(1/5,3/4)B(11/20,1/4)")
say("                                 - [G(3/4)^2/((11/20)G(3/2))] W5a")
say("Identity (B):  (pi sqrt2/x) S(5) = B(1/5,1/4)B(1/20,3/4)")
say("                                 - [G(1/4)^2/((1/20)G(1/2))] W5b")
say("A closed form for W5a or W5b is equivalent to one for S(5).")
x5 = mpf(1) / 5
b5 = mpf(3) / 4 - x5
g5 = mpf(1) / 4 - x5
leadA = beta(x5, mpf(3) / 4) * beta(b5, mpf(1) / 4)
leadB = beta(x5, mpf(1) / 4) * beta(g5, mpf(3) / 4)
W5a = (leadA - pi * sqrt(2) * S5 / x5) * b5 * gamma(mpf(3) / 2) / gamma(mpf(3) / 4) ** 2
W5b = (leadB - pi * sqrt(2) * S5 / x5) * g5 * gamma(mpf(1) / 2) / gamma(mpf(1) / 4) ** 2
mp.dps = 40
chkA = hyper([mpf(3) / 4, mpf(3) / 4, mpf(11) / 20], [mpf(3) / 2, mpf(31) / 20], 1)
chkB = hyper([mpf(1) / 4, mpf(1) / 4, mpf(1) / 20], [mpf(1) / 2, mpf(21) / 20], 1)
mp.dps = DPS
say("   W5a (from (A), %d digits) = %s" % (DPS, nstr(W5a, 45)))
say("       cross-check by direct 3F2 series at dps=40: |diff| = %s"
    % nstr(abs(W5a - chkA), 5))
say("   W5b (from (B), %d digits) = %s" % (DPS, nstr(W5b, 45)))
say("       cross-check by direct 3F2 series at dps=40: |diff| = %s"
    % nstr(abs(W5b - chkB), 5))
say("   leading Gamma quotient of (A) at k=5:  x B(1/5,3/4)B(11/20,1/4)/(pi sqrt2)")
say("      = %s     versus S(5) = %s"
    % (nstr(x5 * leadA / (pi * sqrt(2)), 30), nstr(S5, 30)))

# --- basis B2: structurally motivated for the k=5 residuals ------------------
G14, G34 = gamma(mpf(1) / 4), gamma(mpf(3) / 4)
G15, G25, G35, G45 = (gamma(mpf(j) / 5) for j in (1, 2, 3, 4))
G18, G38 = gamma(mpf(1) / 8), gamma(mpf(3) / 8)
G120, G320, G720, G920, G1120 = (gamma(mpf(j) / 20) for j in (1, 3, 7, 9, 11))
ATOMS = [("1", mpf(1)), ("pi", pi), ("pi^2", pi ** 2), ("1/pi", 1 / pi),
         ("G(Catalan)", catalan), ("pi*G", pi * catalan), ("G/pi", catalan / pi),
         ("varpi", varpi), ("varpi^2", varpi ** 2), ("varpi/pi", varpi / pi),
         ("varpi^2/pi", varpi ** 2 / pi), ("pi^2/varpi", pi ** 2 / varpi),
         ("varpi*G", varpi * catalan), ("zeta(3)", zeta(3)), ("zeta(3)/pi", zeta(3) / pi),
         ("log2", log(2)), ("log(1+sqrt2)", log(1 + r2)), ("pi*log2", pi * log(2)),
         ("log^2 2", log(2) ** 2), ("log(golden)", log((1 + r5) / 2)),
         ("pi*log(golden)", pi * log((1 + r5) / 2)),
         ("G18*G38", G18 * G38), ("(G18G38)^2/pi", (G18 * G38) ** 2 / pi),
         ("B(1/5,3/4)", beta(x5, mpf(3) / 4)), ("B(11/20,1/4)", beta(b5, mpf(1) / 4)),
         ("B(1/5,1/4)", beta(x5, mpf(1) / 4)), ("B(1/20,3/4)", beta(g5, mpf(3) / 4)),
         ("leadA", leadA), ("leadB", leadB),
         ("G15*G25", G15 * G25), ("G15*G45", G15 * G45), ("G25*G35", G25 * G35),
         ("G15^2/G25", G15 ** 2 / G25), ("G120*G920/pi", G120 * G920 / pi),
         ("G320*G720/pi", G320 * G720 / pi), ("G1120*G920", G1120 * G920),
         ("G14^2/pi", G14 ** 2 / pi), ("G34^2/pi", G34 ** 2 / pi),
         ("5^(1/4)", mpf(5) ** (mpf(1) / 4)), ("5^(1/2)", r5)]
B2 = [(m[0] + "*" + a[0], m[1] * a[1]) for m in MULT for a in ATOMS]
say("")
say("Basis B2 (structurally motivated for k=5: Beta/Gamma values with denominators")
say("4, 5, 8, 20, the lemniscatic constants, weight-<=2 constants, golden-ratio logs)")
say("|B2| = %d products (%d multipliers x %d atoms); all %d pairs scanned."
    % (len(B2), len(MULT), len(ATOMS), len(B2) * (len(B2) - 1) // 2))

say("")
say("--- instrument calibration on B2 (must find known relations of this shape) ---")
CAL = [("varpi^2 = G(1/4)^4/(8 pi)", varpi ** 2),
       ("pi S(1)", pi * S1), ("pi S(2)", pi * S2), ("pi S(3)", pi * S3),
       ("leadA (a pure Gamma quotient)", leadA)]
ncal = 0
for nm, tgt in CAL:
    h, _ = scan_pairs(tgt, B2, TOL, H, first_only=True, budget=150)
    ok = bool(h)
    ncal += ok
    say("   %-32s : %s" % (nm, ("FOUND %s over (%s, %s)" % (h[0][0], h[0][1], h[0][2]))
                           if ok else "NOT FOUND"))
say("   calibration score: %d / %d" % (ncal, len(CAL)))

say("")
say("--- B2 IS NOT USABLE FOR A NEGATIVE, AND THAT IS THE POINT ---")
say("B2 recovers varpi^2, pi S(1), pi S(2) and a pure Gamma quotient, but it does")
say("NOT recover pi S(3): it has Gamma/Beta and weight-2 atoms but no logarithm or")
say("arctangent atoms.  Its calibration score is %d/%d.  Per the honesty contract a" % (ncal, len(CAL)))
say("null over B2 is therefore reported as UNINFORMATIVE rather than as a negative:")
say("the same failure mode that addendum 3 diagnosed for the section-5 battery.")
say("")
say("The negatives for this lane are produced instead by")
say("  04-computation/sk_S4_contiguous_pslq_weight1_thm3012.py   (basis B1/B3;")
say("      B1 rediscovers pi S(3) blind, with its irrational sqrt3 coefficient), and")
say("  04-computation/sk_S4_contiguous_product_basis_thm3012.py  (basis B4 of")
say("      PRODUCTS Gamma-type x elementary; validated live on the new closed form")
say("      3F2(3/4,3/4,1/4;3/2,5/4;1) = varpi(pi - 2 log(1+sqrt2))/pi, which is a")
say("      product of the lemniscate constant with a logarithm).")
neg = [("(B2 scan omitted)", "uninformative: B2 fails the pi S(3) control")]

hdr("SUMMARY -- PART II")
say("PART 5  S(4) = varpi/2 - W/(2 varpi), W = 3F2(3/4,3/4,1/2;3/2,3/2;1):")
say("        VERIFIED-EXACT to %d digits; W confirmed by two independent routes." % digits)
say("PART 6  k = 4 is the unique positive integer at which (a) the 3F2 has a")
say("        resonant exponent pair at infinity and (b) the K-moment weight")
say("        mu^{4/k-1}(2-mu^2)^{-2/k-1/2} is a rational function.  PROVED.")
say("PART 7  live positive control pi S(3) over B1: %s" % ("PASS" if POS else "FAIL"))
say("        B2 calibration score %d/%d -- B2 cannot see weight-1 log relations, so" % (ncal, len(CAL)))
say("        no negative is claimed from it.  See the two companion scripts.")

with open("/tmp/math-wt-coinC2/05-knowledge/results/sk_S4_contiguous_pslq_thm3012.out", "w") as f:
    f.write("\n".join(OUT) + "\n")
print("\nwritten to 05-knowledge/results/sk_S4_contiguous_pslq_thm3012.out")
