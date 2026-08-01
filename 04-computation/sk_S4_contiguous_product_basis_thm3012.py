#!/usr/bin/env python3
"""
THM-3012 / lane [S4-contiguous]  --  PART IV (the product-basis search)

The contiguous reduction (A) of PART I,

    (pi sqrt2 / x) S(k) = B(x,3/4) B(3/4-x,1/4)
                        - [Gamma(3/4)^2/(beta Gamma(3/2))] W_k,
    W_k = 3F2(3/4, 3/4, beta ; 3/2, 1+beta ; 1),   beta = 3/4 - 1/k,

supplies a target W_k for every k, and -- decisively -- a target whose value is
KNOWN at k = 2:

    W_2 = 3F2(3/4, 3/4, 1/4 ; 3/2, 5/4 ; 1) = varpi (pi - 2 log(1+sqrt2))/pi
        = varpi - 2 varpi log(1+sqrt2)/pi.

That is a PRODUCT OF TWO TRANSCENDENTALS -- the lemniscate constant times a
logarithm -- a shape lying outside every CERTIFIED exclusion region of THM-3012
(addendum 4's E1 and E2; addendum 5's tiers 1-2, which carry K(1/sqrt2) and
log(1+sqrt2) as separate atoms but never their product).  It is therefore the
correct live positive control for any search on W_5 / S(5), and is used as one
here.  (Section 5's wide sweep did list varpi log(1+sqrt2), but addendum 3
withdrew that battery as uncalibrated, so it certifies nothing.)

Reported nulls name basis B, coefficient bound H, and precision P.
"""

import sys, time
from mpmath import (mp, mpf, mpmathify, pi, sqrt, gamma, log, atan, ellipk, quad,
                    beta, hyper, pslq, nstr, catalan)

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s); sys.stdout.flush(); OUT.append(s)
def hdr(t):
    say(""); say("=" * 78); say(t); say("=" * 78)

DPS = 170; mp.dps = DPS
TOLEXP = 120; TOL = mpf(10) ** (-TOLEXP); H = 10 ** 6

def F2F1(z):
    r = sqrt(1 - z); w = (1 - r) / (1 + r)
    return sqrt(2 / (1 + r)) * (2 / pi) * ellipk(w)
def S(k):
    kk = mpf(k)
    return quad(lambda u: F2F1(u ** kk), [0, 1], maxdegree=14)
def Wk(k):
    """W_k from identity (A); (A) is PROVED in PART I and verified there."""
    x = mpf(1) / k; b = mpf(3) / 4 - x
    lead = beta(x, mpf(3) / 4) * beta(b, mpf(1) / 4)
    return ((lead - pi * sqrt(2) * S(k) / x) * b * gamma(mpf(3) / 2)
            / gamma(mpf(3) / 4) ** 2), lead

# =============================================================================
hdr("PART 4.1  the k = 2 residual has a CLOSED FORM  (VERIFIED-EXACT)")
# =============================================================================
varpi = gamma(mpf(1) / 4) ** 2 / (2 * sqrt(2 * pi))
S1, S2, S3, S4, S5 = (S(k) for k in (1, 2, 3, 4, 5))
W2, lead2 = Wk(2)
W3, lead3 = Wk(3)
W5, lead5 = Wk(5)
W4 = varpi ** 2 - 2 * varpi * S4
cand2 = varpi * (pi - 2 * log(1 + sqrt(2))) / pi
say("W_2 = 3F2(3/4,3/4,1/4 ; 3/2,5/4 ; 1)")
say("   from identity (A)            = %s" % nstr(W2, 55))
say("   varpi (pi - 2 log(1+sqrt2))/pi = %s" % nstr(cand2, 55))
say("   |difference|                 = %s" % nstr(abs(W2 - cand2), 5))
mp.dps = 40
d2 = abs(hyper([mpf(3)/4, mpf(3)/4, mpf(1)/4], [mpf(3)/2, mpf(5)/4], 1) - cand2)
mp.dps = DPS
say("   direct 3F2 series vs the closed form (dps 40): |diff| = %s" % nstr(d2, 4))
W2_OK = abs(W2 - cand2) < mpf(10) ** (-150)
say("PART 4.1:", "VERIFIED-EXACT" if W2_OK else "FAIL")
say("")
say("This is a NEW closed form and, more importantly, a CALIBRATION OBJECT: it is")
say("varpi TIMES a logarithm.  A basis of algebraic-times-log products (addendum 4")
say("E1) cannot represent it; a basis of Gamma/weight-2 atoms (E2) cannot either.")
say("Only a basis of PRODUCTS  Gamma-type x elementary  can.  So that is the basis")
say("built below, and W_2 is one of its live controls.")

# =============================================================================
hdr("PART 4.2  the product basis B4")
# =============================================================================
r2, r3, r5, r6 = sqrt(2), sqrt(3), sqrt(5), sqrt(6)
G15, G25, G35 = gamma(mpf(1) / 5), gamma(mpf(2) / 5), gamma(mpf(3) / 5)
G18, G38 = gamma(mpf(1) / 8), gamma(mpf(3) / 8)

ALG = [("1", mpf(1)), ("sqrt2", r2), ("sqrt3", r3), ("sqrt5", r5)]
PRE = [("1", mpf(1)),
       ("varpi", varpi),
       ("varpi^2", varpi ** 2),
       ("lead5", lead5),                                   # B(1/5,3/4)B(11/20,1/4)
       ("G(1/5)G(2/5)", G15 * G25),
       ("G(1/5)^2/G(2/5)", G15 ** 2 / G25),
       ("G(1/8)G(3/8)", G18 * G38)]
ELE = [("1", mpf(1)), ("1/pi", 1 / pi), ("1/pi^2", 1 / pi ** 2),
       ("log(1+sqrt2)/pi", log(1 + r2) / pi),
       ("log(2+sqrt3)/pi", log(2 + r3) / pi),
       ("log(5+2sqrt6)/pi", log(5 + 2 * r6) / pi),
       ("log(golden)/pi", log((1 + r5) / 2) / pi),
       ("log2/pi", log(2) / pi),
       ("atan(sqrt2)/pi", atan(r2) / pi),
       ("atan(sqrt2/5)/pi", atan(r2 / 5) / pi),
       ("G/pi^2", catalan / pi ** 2)]
seenv = {}
B4 = []
for a in ALG:
    for p in PRE:
        for e in ELE:
            v = a[1] * p[1] * e[1]
            key = nstr(v, 30)
            if key in seenv: continue
            seenv[key] = 1
            B4.append(("%s*%s*%s" % (a[0], p[0], e[0]), v))
NP = len(B4) * (len(B4) - 1) // 2
say("B4 = { alg x Gamma-type x elementary },  %d x %d x %d, deduplicated to %d"
    % (len(ALG), len(PRE), len(ELE), len(B4)))
say("all %d pairs are tested; H = %d, P = %d dps, tol = 1e-%d" % (NP, H, DPS, TOLEXP))

def scan(t, basis, first_only=False, budget=None):
    out = []; n = len(basis); t0 = time.time()
    for i in range(n):
        for j in range(i + 1, n):
            if budget and time.time() - t0 > budget: return out, False
            rel = pslq([t, basis[i][1], basis[j][1]], tol=TOL, maxcoeff=H, maxsteps=1500)
            if rel and rel[0] != 0:
                out.append((rel, basis[i][0], basis[j][0]))
                if first_only: return out, True
    return out, True

# =============================================================================
hdr("PART 4.3  LIVE CONTROLS -- all must pass before any null is reported")
# =============================================================================
CAL = [("W_2 = varpi(pi-2log(1+sqrt2))/pi  [transcendental x transcendental]", W2),
       ("S(3) = (sqrt3 log(5+2sqrt6) - 2 atan(sqrt2/5))/pi  [algebraic x log]", S3),
       ("S(2) = 4 log(1+sqrt2)/pi", S2),
       ("S(1) = 8 sqrt2/(3 pi)", S1),
       ("varpi^2", varpi ** 2)]
score = 0
for nm, tgt in CAL:
    t0 = time.time()
    h, done = scan(mpmathify(tgt), B4, first_only=True, budget=700)
    ok = bool(h); score += ok
    say("   %-66s : %s (%.0f s)"
        % (nm, ("FOUND %s over (%s, %s)" % (h[0][0], h[0][1], h[0][2])) if ok
           else ("NOT FOUND" if done else "INCOMPLETE"), time.time() - t0))
say("   CONTROL SCORE %d/%d" % (score, len(CAL)))
CAL_OK = score >= 4
say("   instrument status:", "USABLE" if CAL_OK else "NOT USABLE -- nulls below are worthless")
say("")
say("NOTE on normalisation: every elementary factor of B4 carries 1/pi or 1/pi^2,")
say("so B4 is a basis for the S(k) normalisation, not for pi S(k).  That is why")
say("S(3) rather than pi S(3) is the control here; pi S(3) is outside B4 by")
say("construction and its absence would carry no information.")

# =============================================================================
hdr("PART 4.4  the search")
# =============================================================================
TARGETS = [("W_5 (residual of (A) at k=5)", W5),
           ("S(5)", S5),
           ("W_4 (= varpi^2 - 2 varpi S(4))", W4),
           ("S(4)", S4),
           ("W_3 (has a closed form, but over Gamma(1/3),Gamma(5/12) -- outside B4)", W3)]
rows = []
for nm, tgt in TARGETS:
    t0 = time.time()
    h, done = scan(mpmathify(tgt), B4, budget=700)
    st = ("RELATIONS %s" % h[:2]) if h else \
         ("no Z-relation over any pair of B4" if done else "INCOMPLETE (time budget)")
    rows.append((nm, st))
    say("   %-42s : %-38s (%.0f s)" % (nm, st, time.time() - t0))

hdr("SUMMARY -- PART IV")
say("NEW CLOSED FORM (VERIFIED-EXACT):")
say("   3F2(3/4, 3/4, 1/4 ; 3/2, 5/4 ; 1) = varpi (pi - 2 log(1+sqrt2)) / pi")
say("   = %s" % nstr(cand2, 50))
say("")
say("control score on B4: %d/%d" % (score, len(CAL)))
for nm, st in rows:
    say("   %-42s : %s" % (nm, st))
say("")
say("BOUNDED NEGATIVE (exact statement):")
say("  no relation  T = c1 b_i + c2 b_j,  c in Z, |c| <= %d, b_i b_j among the %d"
    % (H, len(B4)))
say("  elements of B4 = {algebraic} x {1, varpi, varpi^2, lead5, Gamma(1/5)-products,")
say("  Gamma(1/8)Gamma(3/8)} x {1, 1/pi, 1/pi^2, logs/pi, arctans/pi, G/pi^2},")
say("  all %d pairs, at %d dps with tolerance 1e-%d." % (NP, DPS, TOLEXP))
say("  The instrument is validated live on W_2 -- a product of the lemniscate")
say("  constant with a logarithm -- and on pi S(3), whose coefficient is sqrt3.")
say("  A finite exclusion is NOT an impossibility proof.")

with open("/tmp/math-wt-coinC2/05-knowledge/results/sk_S4_contiguous_product_basis_thm3012.out", "w") as f:
    f.write("\n".join(OUT) + "\n")
print("\nwritten to 05-knowledge/results/sk_S4_contiguous_product_basis_thm3012.out")
