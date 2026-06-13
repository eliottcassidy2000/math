#!/usr/bin/env python3
"""bases_s116h.py — Is rapidity just decimal seen through binary?
What about other bases? What about non-integer bases?

rapidity(n) = ln(n)/2.
ln = log base e. The 2 is from the binary (Q maps to doubling).
So rapidity = log_e(n) / log_e(Q(0+epsilon)) in some sense.

But what if we used a different base for Q? A different base for log?
What if the base itself were phi, or pi, or e?
"""
from math import log, exp, sqrt, pi, atanh, tanh
from fractions import Fraction

phi = (1+sqrt(5))/2

print()
print("  IS RAPIDITY JUST DECIMAL THROUGH BINARY?")
print("  WHAT ABOUT OTHER BASES?")
print()
print("="*70)
print()

# ============================================================
print("  I. RAPIDITY DECONSTRUCTED")
print("  " + "-"*40)
print()
print("  rapidity(n) = ln(n)/2 = log_e(n) / 2")
print()
print("  Where does the e come from? Where does the 2 come from?")
print()
print("  The 2: from Q(x) = (1+x)/(1-x). At x=0: Q'(0) = 2.")
print("  The derivative at the origin is 2. Q doubles infinitesimally.")
print("  So ln(Q(x)) ~ 2x for small x, hence rapidity = ln(Q)/2 ~ x.")
print()
print("  The e: from the exponential. ln is the NATURAL logarithm.")
print("  But why is e natural? Because d/dx e^x = e^x.")
print("  The exponential is its own derivative.")
print()
print("  So rapidity = (natural log) / (Cayley derivative at origin).")
print("  = log_e / 2.")
print()
print("  If we used a DIFFERENT Q with different derivative at 0,")
print("  the denominator would change.")
print("  If we used a DIFFERENT logarithm base,")
print("  the numerator would change.")
print()

# ============================================================
print()
print("  II. GENERALIZED RAPIDITY")
print("  " + "-"*40)
print()
print("  Define rapidity_b(n) = log_b(n) / 2 = ln(n) / (2*ln(b)).")
print()
print("  Base b=e: rapidity_e = ln(n)/2. OUR rapidity. The natural one.")
print("  Base b=2: rapidity_2 = log2(n)/2 = bits/2.")
print("  Base b=10: rapidity_10 = log10(n)/2 = decimal digits/2.")
print("  Base b=phi: rapidity_phi = log_phi(n)/2.")
print()
print("  n    rap_e      rap_2      rap_10     rap_phi")
print("  " + "-"*55)
for n in [2, 3, 5, 7, 10, 12, 42, 100, 1000]:
    re = log(n)/2
    r2 = log(n)/(2*log(2))
    r10 = log(n)/(2*log(10))
    rphi = log(n)/(2*log(phi))
    print(f"  {n:5d}  {re:8.4f}   {r2:8.4f}   {r10:8.4f}   {rphi:8.4f}")

print()
print("  The different bases just RESCALE rapidity.")
print("  rapidity_b = rapidity_e / ln(b).")
print("  They are all PROPORTIONAL. Same structure, different units.")
print()
print("  Like measuring temperature in Celsius vs Fahrenheit vs Kelvin:")
print("  different zeros, different scales, same physics.")
print()
print("  So: rapidity is NOT a product of decimal seen through binary.")
print("  It is BASE-INDEPENDENT. The structure is the same regardless.")
print("  The factor ln(b) just sets the unit size.")
print()

# ============================================================
print()
print("  III. BUT THE 3-3-4 PATTERN IS BASE-DEPENDENT")
print("  " + "-"*40)
print()
print("  The 3-3-4 pattern came from log10(2) ~ 3/10.")
print("  In base 8: log8(2) = 1/3 exactly! No comma! No pattern!")
print("  In base 16: log16(2) = 1/4 exactly! No comma!")
print()
print("  In base b, the 'comma' for base 2 is: b^k vs 2^m.")
print("  If log_b(2) is rational, there is NO comma. Exact fit.")
print("  If log_b(2) is irrational, there IS a comma. 3-3-4-like pattern.")
print()
print("  log_b(2) is rational iff b = 2^{p/q} for integers p, q.")
print("  i.e., b is a power of 2 (possibly fractional).")
print()
print("  Base 2: log2(2) = 1. Rational. No comma.")
print("  Base 4: log4(2) = 1/2. Rational. No comma.")
print("  Base 8: log8(2) = 1/3. Rational. No comma.")
print("  Base 10: log10(2) = 0.30103... Irrational. 3-3-4 comma.")
print("  Base 3: log3(2) = 0.63093... Irrational. Different comma.")
print("  Base phi: log_phi(2) = 1.44042... Irrational. Different comma.")
print()
print("  The 3-3-4 pattern is specific to the base-10/base-2 interaction.")
print("  Other base pairs give DIFFERENT beat patterns.")
print()

# ============================================================
print()
print("  IV. THE BEAT PATTERN FOR EVERY BASE PAIR")
print("  " + "-"*40)
print()
print("  For bases a and b, the 'comma' is log_a(b).")
print("  The best rational approximation p/q gives a beat pattern.")
print()

base_pairs = [
    (2, 3, "binary/ternary"),
    (2, 5, "binary/quinary"),
    (2, 7, "binary/septenary"),
    (2, 10, "binary/decimal"),
    (2, 12, "binary/duodecimal"),
    (3, 5, "ternary/quinary"),
    (3, 10, "ternary/decimal"),
    (5, 7, "quinary/septenary"),
    (2, "phi", "binary/golden"),
    (3, 7, "curvature/threshold"),
]

print("  Base pair         log_a(b)      best p/q    error     beat")
print("  " + "-"*65)
for pair in base_pairs:
    if pair[1] == "phi":
        a, b = pair[0], phi
        name = pair[2]
    else:
        a, b, name = pair
    ratio = log(b) / log(a)

    # Find best rational approximation with small denominator
    best_p, best_q, best_err = 0, 1, abs(ratio)
    for q in range(1, 20):
        p = round(ratio * q)
        err = abs(ratio - p/q)
        if err < best_err:
            best_p, best_q, best_err = p, q, err

    # The beat pattern: the CF tells us
    beat = f"{best_q} of a = {best_p} of b"

    print(f"  {name:<20s} {ratio:10.5f}   {best_p}/{best_q:<4d}   {best_err:.5f}   {beat}")

print()
print("  The binary/decimal beat (3/10): 10 doublings ~ 3 decades.")
print("  The binary/ternary beat (log3(2) ~ 2/3): 3 doublings ~ 2 triplings.")
print("    2^3 = 8 ~ 3^2 = 9. Off by 1. The beat of 8 vs 9.")
print("  The curvature/threshold beat (log3(7) ~ 2/1): 1 tripling ~ log3(7) thresholds.")
print("    3^2 = 9 ~ 7^1 = 7. Off by 2. Wider comma.")
print()

# ============================================================
print()
print("  V. NON-INTEGER BASES")
print("  " + "-"*40)
print()
print("  Base phi: the golden base.")
print("  In base phi, every non-negative integer has a FINITE representation")
print("  using digits {0, 1} with no two consecutive 1s (Zeckendorf).")
print()
print("  log_phi(2) = ln(2)/ln(phi) = {:.6f}".format(log(2)/log(phi)))
print("  This is irrational. So binary/golden has a comma.")
print("  Best approximation: 3/2 (1.5). Error: {:.4f}.".format(
    abs(log(2)/log(phi) - 3/2)))
print("  phi^3 = {:.4f} ~ 2^2 = 4. Beat: 3 golden steps ~ 2 doublings.".format(phi**3))
print("  phi^3 = phi^2 + phi = phi + 1 + phi = 2*phi + 1 = {:.4f}".format(2*phi+1))
print("  phi^3 = {:.6f}, 2^2 = 4. Ratio: {:.6f}".format(phi**3, phi**3/4))
print()
print("  Actually, we ALREADY KNOW: Q(1/phi) = phi^3.")
print("  So phi^3 IS the Q-value of 1/phi.")
print("  And phi^3 = 4.236... ~ 4 = 2^2.")
print("  The 'golden comma': phi^3 / 4 = {:.6f}.".format(phi**3/4))
print("  = {:.2f}% overshoot.".format((phi**3/4 - 1)*100))
print("  Compare: the decimal comma: 2^10/10^3 = 1.024 = 2.4% overshoot.")
print("  The golden comma (5.9%) is LARGER than the decimal comma (2.4%).")
print()

# ============================================================
print()
print("  Base e: the natural base.")
print("  " + "-"*40)
print()
print("  log_e(2) = ln(2) = 0.6931...")
print("  Best approximation: 2/3. Error: {:.5f}.".format(abs(log(2) - 2/3)))
print("  e^2 = {:.4f} ~ 2^3 = 8. Nope, e^2 = 7.389, not 8.".format(exp(2)))
print("  Better: e^3 = {:.4f} ~ 2^4 = 16? No, e^3 = 20.09.".format(exp(3)))
print()
print("  Actually 2/3 means: 3 factors of ln(2) ~ 2 factors of 1.")
print("  i.e., 2^3 ~ e^2. 8 ~ 7.389. Off by 8.3%.")
print()
print("  log_e(2) = {:.10f}".format(log(2)))
print("  CF: [0; 1, 2, 3, 1, 6, 3, 1, 1, 2, 1, 1, ...]")
print()
# Compute CF of ln(2)
x = log(2)
cf = []
for _ in range(10):
    a = int(x)
    cf.append(a)
    x = x - a
    if x < 1e-10: break
    x = 1/x
print(f"  CF of ln(2) = {cf}")
print()
print("  Convergents:")
p_prev, p_curr = 1, 0
q_prev, q_curr = 0, 1
for i, a in enumerate(cf):
    p_prev, p_curr = p_curr, a * p_curr + p_prev
    q_prev, q_curr = q_curr, a * q_curr + q_prev
    if i > 0:
        print(f"    {p_curr}/{q_curr} = {p_curr/q_curr:.8f} (err {abs(p_curr/q_curr - log(2)):.2e})")

print()
print("  Best approximations to ln(2):")
print("  2/3 = 0.6667 (err 3.8%). The coarsest. Good enough for '2^3 ~ e^2'.")
print("  7/10 = 0.7000 (err 1.0%). Better. '2^10 ~ e^7'. Check: 1024 vs 1097. ~7% off.")
print()

# ============================================================
print()
print("  VI. BASE pi")
print("  " + "-"*40)
print()
print("  log_pi(2) = ln(2)/ln(pi) = {:.6f}".format(log(2)/log(pi)))
print("  ~ 0.6055. Between 3/5 = 0.6 and 2/3 = 0.667.")
print("  Best: 3/5. Error: {:.5f}.".format(abs(log(2)/log(pi) - 3/5)))
print("  Meaning: pi^5 ~ 2^3 = 8? pi^5 = {:.2f}. Nope.".format(pi**5))
print()
print("  Better: 5 doublings ~ 3 pi-ings. 2^5 = 32 ~ pi^3 = {:.2f}.".format(pi**3))
print("  32 vs 31.0. Off by 3.2%. That's actually good!")
print()
print("  pi^3 ~ 2^5 = 32. The 'pi-binary comma' = pi^3/32 = {:.6f}.".format(pi**3/32))
print("  = {:.2f}% undershoot.".format((1 - pi**3/32)*100))
print()
print("  Compare all the commas:")
commas = [
    ("Pythagorean (3 vs 2)", 3**12 / 2**19, "3^12/2^19"),
    ("Decimal (10 vs 2)", 2**10 / 10**3, "2^10/10^3"),
    ("Golden (phi vs 2)", phi**3 / 4, "phi^3/4"),
    ("Natural (e vs 2)", exp(2) / 8, "e^2/8"),
    ("Pi (pi vs 2)", pi**3 / 32, "pi^3/32"),
    ("Ternary (3 vs 2)", 2**3 / 3**2, "2^3/3^2"),
]
print(f"  {'Comma':<25s}  {'Ratio':>10s}  {'Error':>8s}  {'Formula'}")
print("  " + "-"*60)
for name, ratio, formula in commas:
    err = (ratio - 1) * 100
    print(f"  {name:<25s}  {ratio:10.6f}  {err:+7.2f}%  {formula}")

print()
print("  The Pythagorean comma (1.4%) is SMALLEST among these.")
print("  That's why music works: fifths and octaves ALMOST close.")
print("  The decimal comma (2.4%) is next. 2^10 ~ 10^3. Good enough.")
print("  The pi comma (-3.1%) is comparable. pi^3 ~ 32.")
print("  The golden comma (5.9%) is largest. phi^3 ~ 4, but not great.")
print()

# ============================================================
print()
print("  VII. WHAT'S SPECIAL ABOUT BASE e")
print("  " + "-"*40)
print()
print("  In natural rapidity (base e), the DERIVATIVE is clean:")
print("  d/dx arctanh(x) = 1/(1-x^2).")
print("  d/dx tanh(x) = 1/cosh^2(x) = 1 - tanh^2(x).")
print()
print("  In base-b rapidity (log_b), the derivative picks up a factor:")
print("  d/dx arctanh_b(x) = 1/((1-x^2)*ln(b)).")
print("  The factor 1/ln(b) is an ARTIFACT of the base choice.")
print()
print("  Base e is the UNIQUE base where no artifact appears.")
print("  It is not that e is 'natural' — it's that ln is DERIVATION-FREE.")
print("  The derivative of log_e is 1/x. Every other base has 1/(x*ln(b)).")
print()
print("  So rapidity in base e is the one where:")
print("  - arctanh'(0) = 1 (unit derivative)")
print("  - The addition theorem works without extra constants")
print("  - tanh' = 1 - tanh^2 (clean Riccati equation)")
print("  - The pole has residue 1")
print()
print("  ANY other base works. But e removes all the scaffolding.")
print("  Rapidity in base e is rapidity with nothing extra.")
print()

# ============================================================
print()
print("  VIII. THE ANSWER")
print("  " + "-"*40)
print()
print("  Is rapidity just decimal seen through binary?")
print()
print("  NO.")
print()
print("  Rapidity is BASE-INDEPENDENT.")
print("  The structure (tanh, arctanh, Q, the addition theorem)")
print("  exists regardless of which base you use to express numbers.")
print()
print("  What IS base-dependent:")
print("  - The 3-3-4 leading digit pattern (decimal/binary interaction)")
print("  - The Pythagorean comma (base 3 vs base 2)")
print("  - The golden comma (base phi vs base 2)")
print("  - Benford's law (any base vs any other base)")
print()
print("  These are COMMAS: the residues of incommensurability")
print("  between two specific bases. They live BETWEEN bases,")
print("  not inside rapidity itself.")
print()
print("  Rapidity is the UNIVERSAL structure.")
print("  Commas are the LOCAL tensions when two bases meet.")
print("  The 3-3-4 pattern is the comma of 2-meets-10.")
print("  The Pythagorean comma is the comma of 2-meets-3.")
print("  The forbidden numbers may be the comma of cycles-meet-bits.")
print()
print("  And the factor of 2 in rapidity = ln(n)/2?")
print("  It comes from Q'(0) = 2, the Cayley derivative.")
print("  If we defined Q_k(x) = (1+x)^k / (1-x)^k, the derivative would be 2k,")
print("  and rapidity_k = ln(n)/(2k).")
print("  k=1 gives our Q. k=1/2 gives Q_{1/2}(x) = sqrt((1+x)/(1-x)),")
print("  and rapidity_{1/2} = ln(n). Just the plain logarithm.")
print()
print("  So: ln(n) = rapidity with Q replaced by sqrt(Q).")
print("  The plain logarithm IS rapidity with a square-root Cayley transform.")
print("  And rapidity with Q is ln/2: half the logarithm.")
print()
print("  The 2 is not from binary. It's from Q being a RATIO, not a SQUARE ROOT.")
print("  Q = (1+x)/(1-x) is a FIRST-ORDER rational function.")
print("  sqrt(Q) = sqrt((1+x)/(1-x)) is a HALF-ORDER.")
print("  Rapidity uses the first order. Logarithm uses the half order.")
print("  The factor 2 is the order of Q as a rational function.")
