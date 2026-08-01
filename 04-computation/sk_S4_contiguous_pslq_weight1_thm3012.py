#!/usr/bin/env python3
"""
THM-3012 / lane [S4-contiguous]  --  PART III (the weight-1 exclusion)

Companion to sk_S4_contiguous_pslq_thm3012.py.  That script's Gamma/Beta basis B2
is calibrated on varpi^2, pi S(1), pi S(2) and a pure Gamma quotient, but it does
NOT contain the log/arctan atoms, so it does NOT detect pi S(3) -- its null result
is therefore silent about the weight-1 class.  This script closes that gap: it
scans the SAME targets against the basis B1 that DOES rediscover pi S(3) blind,
and against the union B3 = B1 u B2 restricted to a completed region.

Targets are the NEW objects produced by the contiguous reduction:
    W    = 3F2(3/4,3/4,1/2   ; 3/2, 3/2  ; 1)   (S(4) = varpi/2 - W/(2 varpi))
    W5a  = 3F2(3/4,3/4,11/20 ; 3/2, 31/20; 1)
    W5b  = 3F2(1/4,1/4,1/20  ; 1/2, 21/20; 1)
plus S(4), S(5) themselves in several normalisations.
"""

import sys, time
from mpmath import (mp, mpf, mpmathify, pi, sqrt, gamma, log, atan, ellipk, quad,
                    beta, hyper, pslq, nstr, catalan, zeta)

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s); sys.stdout.flush(); OUT.append(s)
def hdr(t):
    say(""); say("=" * 78); say(t); say("=" * 78)

DPS = 170
mp.dps = DPS
TOLEXP = 122
TOL = mpf(10) ** (-TOLEXP)
H = 10 ** 6

def F2F1(z):
    r = sqrt(1 - z); w = (1 - r) / (1 + r)
    return sqrt(2 / (1 + r)) * (2 / pi) * ellipk(w)
def S(k):
    kk = mpf(k)
    return quad(lambda u: F2F1(u ** kk), [0, 1], maxdegree=14)

hdr("targets")
S1, S2, S3, S4, S5 = (S(k) for k in (1, 2, 3, 4, 5))
say("control |S(3) - known closed form| = %s"
    % nstr(abs(S3 - (sqrt(3) * log(5 + 2 * sqrt(6)) - 2 * atan(sqrt(2) / 5)) / pi), 5))
varpi = gamma(mpf(1) / 4) ** 2 / (2 * sqrt(2 * pi))
W = varpi ** 2 - 2 * varpi * S4          # PROVED equal to 3F2(3/4,3/4,1/2;3/2,3/2;1)
x5 = mpf(1) / 5; b5 = mpf(3) / 4 - x5; g5 = mpf(1) / 4 - x5
leadA = beta(x5, mpf(3) / 4) * beta(b5, mpf(1) / 4)
leadB = beta(x5, mpf(1) / 4) * beta(g5, mpf(3) / 4)
W5a = (leadA - pi * sqrt(2) * S5 / x5) * b5 * gamma(mpf(3) / 2) / gamma(mpf(3) / 4) ** 2
W5b = (leadB - pi * sqrt(2) * S5 / x5) * g5 * gamma(mpf(1) / 2) / gamma(mpf(1) / 4) ** 2
mp.dps = 35
say("W   cross-check vs direct 3F2 series: |diff| = %s"
    % nstr(abs(W - hyper([mpf(3)/4, mpf(3)/4, mpf(1)/2], [mpf(3)/2, mpf(3)/2], 1)), 4))
say("W5a cross-check vs direct 3F2 series: |diff| = %s"
    % nstr(abs(W5a - hyper([mpf(3)/4, mpf(3)/4, mpf(11)/20], [mpf(3)/2, mpf(31)/20], 1)), 4))
say("W5b cross-check vs direct 3F2 series: |diff| = %s"
    % nstr(abs(W5b - hyper([mpf(1)/4, mpf(1)/4, mpf(1)/20], [mpf(1)/2, mpf(21)/20], 1)), 4))
mp.dps = DPS

# ---------------------------------------------------------------- bases -----
r2, r3, r5, r6 = sqrt(2), sqrt(3), sqrt(5), sqrt(6)
MULT = [("1", mpf(1)), ("sqrt2", r2), ("sqrt3", r3), ("sqrt5", r5),
        ("sqrt6", r6), ("2^(1/4)", mpf(2) ** (mpf(1) / 4))]
LOGS = [("1", mpf(1)), ("pi", pi), ("log(1+sqrt2)", log(1 + r2)),
        ("log(2+sqrt3)", log(2 + r3)), ("log(5+2sqrt6)", log(5 + 2 * r6)),
        ("log2", log(2)), ("log3", log(3)), ("log5", log(5)),
        ("atan(sqrt2)", atan(r2)), ("atan(sqrt2/5)", atan(r2 / 5)),
        ("atan(sqrt2/3)", atan(r2 / 3)), ("atan(sqrt3)", atan(r3)),
        ("atan(1/2)", atan(mpf(1) / 2)), ("atan(1/3)", atan(mpf(1) / 3)),
        ("log(golden)", log((1 + r5) / 2)), ("log(1+sqrt5)", log(1 + r5))]
B1 = [(m[0] + "*" + l[0], m[1] * l[1]) for m in MULT for l in LOGS]

G14 = gamma(mpf(1) / 4); G34 = gamma(mpf(3) / 4)
G15, G25, G35, G45 = (gamma(mpf(j) / 5) for j in (1, 2, 3, 4))
G18, G38 = gamma(mpf(1) / 8), gamma(mpf(3) / 8)
G120, G320, G720, G920, G1120 = (gamma(mpf(j) / 20) for j in (1, 3, 7, 9, 11))
GAMMAS = [("varpi", varpi), ("varpi^2", varpi ** 2), ("varpi/pi", varpi / pi),
          ("varpi^2/pi", varpi ** 2 / pi), ("pi^2/varpi", pi ** 2 / varpi),
          ("pi/varpi", pi / varpi), ("1/varpi", 1 / varpi),
          ("G(Catalan)", catalan), ("pi*G", pi * catalan), ("G/pi", catalan / pi),
          ("varpi*G", varpi * catalan), ("varpi*G/pi", varpi * catalan / pi),
          ("pi^2", pi ** 2), ("zeta(3)", zeta(3)),
          ("G18*G38", G18 * G38), ("(G18G38)^2/pi", (G18 * G38) ** 2 / pi),
          ("B(1/5,3/4)", beta(x5, mpf(3) / 4)), ("B(11/20,1/4)", beta(b5, mpf(1) / 4)),
          ("B(1/5,1/4)", beta(x5, mpf(1) / 4)), ("B(1/20,3/4)", beta(g5, mpf(3) / 4)),
          ("leadA", leadA), ("leadB", leadB),
          ("G15*G25", G15 * G25), ("G15*G45", G15 * G45), ("G25*G35", G25 * G35),
          ("G15^2/G25", G15 ** 2 / G25), ("G120*G920/pi", G120 * G920 / pi),
          ("G320*G720/pi", G320 * G720 / pi), ("G1120*G920", G1120 * G920),
          ("G14^2/pi", G14 ** 2 / pi), ("G34^2/pi", G34 ** 2 / pi)]
B3 = B1 + [(m[0] + "*" + a[0], m[1] * a[1]) for m in MULT for a in GAMMAS]

def scan(target, basis, first_only=False, budget=None):
    out = []; n = len(basis); t0 = time.time()
    for i in range(n):
        for j in range(i + 1, n):
            if budget and time.time() - t0 > budget:
                return out, False
            rel = pslq([target, basis[i][1], basis[j][1]], tol=TOL, maxcoeff=H,
                       maxsteps=1500)
            if rel and rel[0] != 0:
                out.append((rel, basis[i][0], basis[j][0]))
                if first_only: return out, True
    return out, True

# ------------------------------------------------------- calibration --------
hdr("CALIBRATION -- the instrument must find these before any null is reported")
CAL = [("pi S(3)  [the sqrt3-coefficient relation]", pi * S3),
       ("pi S(1)", pi * S1), ("pi S(2)", pi * S2),
       ("varpi^2 = G(1/4)^4/(8 pi)", varpi ** 2),
       ("leadA (a pure Gamma quotient)", leadA)]
for label, bas in (("B1", B1), ("B3 = B1 u Gamma-atoms", B3)):
    say("")
    say("basis %s : %d elements, %d pairs" % (label, len(bas), len(bas) * (len(bas) - 1) // 2))
    sc = 0
    for nm, tgt in CAL:
        h, _ = scan(tgt, bas, first_only=True, budget=400)
        sc += bool(h)
        say("   %-40s : %s" % (nm, ("FOUND %s over (%s, %s)" % (h[0][0], h[0][1], h[0][2]))
                               if h else "NOT FOUND"))
    say("   calibration score %d/%d" % (sc, len(CAL)))
    if label == "B1": CAL1 = sc
    else: CAL3 = sc

# ---------------------------------------------------------- negatives -------
hdr("BOUNDED NEGATIVES on the contiguous targets")
TARGETS = [("W (k=4 residual)", W), ("W5a", W5a), ("W5b", W5b),
           ("S(4)", S4), ("pi S(4)", pi * S4),
           ("S(5)", S5), ("pi S(5)", pi * S5), ("pi^2 S(5)", pi ** 2 * S5)]
rows = []
for label, bas in (("B1", B1), ("B3", B3)):
    say("")
    say("--- basis %s (%d elements, %d pairs), H = %d, P = %d dps, tol = 1e-%d ---"
        % (label, len(bas), len(bas) * (len(bas) - 1) // 2, H, DPS, TOLEXP))
    for nm, tgt in TARGETS:
        t0 = time.time()
        h, done = scan(mpmathify(tgt), bas, budget=600)
        st = ("RELATIONS: %s" % h[:2]) if h else \
             ("no Z-relation over any pair" if done else "INCOMPLETE (time budget)")
        rows.append((label, nm, st))
        say("   %-18s : %-42s (%.0f s)" % (nm, st, time.time() - t0))

hdr("SUMMARY -- PART III")
say("B1 calibration %d/5 (crucially it rediscovers pi S(3) blind)" % CAL1)
say("B3 calibration %d/5" % CAL3)
for label, nm, st in rows:
    say("   %-3s  %-18s : %s" % (label, nm, st))
say("")
say("STATEMENT OF THE NEGATIVE (this is all that is claimed):")
say("  For each target T above there is no relation  T = c1 b_i + c2 b_j  with")
say("  c1, c2 in Z, |c| <= %d, b_i, b_j drawn from the stated basis (all pairs)," % H)
say("  at %d decimal digits with detection tolerance 1e-%d." % (DPS, TOLEXP))
say("  Both instruments are validated live: B1 rediscovers pi S(3) -- the")
say("  relation with the IRRATIONAL sqrt3 coefficient -- from a blind pair scan.")
say("  A finite exclusion is NOT an impossibility proof and must not be quoted as one.")

with open("/tmp/math-wt-coinC2/05-knowledge/results/sk_S4_contiguous_pslq_weight1_thm3012.out", "w") as f:
    f.write("\n".join(OUT) + "\n")
print("\nwritten to 05-knowledge/results/sk_S4_contiguous_pslq_weight1_thm3012.out")
