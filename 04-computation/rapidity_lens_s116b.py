#!/usr/bin/env python3
"""rapidity_lens_s116b.py — See everything through rapidity

Rapidity phi = arctanh(v). The fundamental quantity.
Q(v) = e^{2*phi}. Multiplication of Q = addition of rapidity.
Our entire theory lives in rapidity space.
"""
from math import sqrt, log, pi, cos, sin, exp, atanh, tanh, acosh, asinh, comb, factorial
from fractions import Fraction

phi_gr = (1+sqrt(5))/2  # golden ratio
tau = 1.8392867552141612  # tribonacci

def is_prime(n):
    if n < 2: return False
    for d in range(2, int(sqrt(n))+1):
        if n % d == 0: return False
    return True

primes = [p for p in range(2, 200) if is_prime(p)]

print("RAPIDITY: THE FUNDAMENTAL QUANTITY")
print("="*60)
print()
print("  Rapidity: phi = arctanh(v) = (1/2)*ln((1+v)/(1-v)) = (1/2)*ln(Q(v))")
print("  Inverse:  v = tanh(phi)")
print("  Cayley:   Q(v) = e^{2*phi}")
print()
print("  ADDITION of rapidities = MULTIPLICATION of Q-values")
print("  = COMPOSITION of velocities (relativistic addition)")
print("  Rapidity is the LOG of our theory. It linearizes everything.")
print()

# ============================================================
print("="*60)
print("1. RAPIDITY OF EVERY NATURAL NUMBER")
print("="*60)
print()
print("  The Cayley address of n is x_n = (n-1)/(n+1).")
print("  Its rapidity is arctanh(x_n) = arctanh((n-1)/(n+1)).")
print("  = (1/2)*ln(n^2/(2n/(n+1))... let's compute.")
print("  Actually: (1+x)/(1-x) = (1+(n-1)/(n+1))/(1-(n-1)/(n+1))")
print("  = (2n/(n+1))/(2/(n+1)) = n. So Q(x_n) = n.")
print("  Therefore: rapidity of n = (1/2)*ln(n).")
print()
print("  THE RAPIDITY OF n IS HALF ITS LOGARITHM.")
print()
print("  n    address x_n    rapidity = ln(n)/2     e^(2*rapidity)")
print("  " + "-"*60)
for n in range(1, 16):
    x = (n-1)/(n+1)
    rap = log(n)/2
    print(f"  {n:3d}   {x:10.6f}       {rap:10.6f}            {exp(2*rap):10.4f}")

print()
print("  Rapidity is LOGARITHMIC in n. The natural numbers are")
print("  EQUALLY SPACED on the log scale, which means their rapidities")
print("  form an ARITHMETIC PROGRESSION with common difference ln(r)/2")
print("  where r is the ratio between consecutive naturals.")
print()
print("  Between n and n+1: delta_rapidity = ln(1+1/n)/2 ~ 1/(2n).")
print("  This is the HARMONIC SERIES! The rapidity gaps between")
print("  consecutive integers decrease harmonically.")
print()

# ============================================================
print("="*60)
print("2. RAPIDITY OF MUSICAL INTERVALS")
print("="*60)
print()
print("  Musical interval (n+1)/n has Cayley pre-image 1/(2n+1).")
print("  Rapidity of the interval = arctanh(1/(2n+1))")
print("  = (1/2)*ln(Q(1/(2n+1))) = (1/2)*ln((n+1)/n)")
print()
print("  Interval      Ratio   Rapidity    Cents")
print("  " + "-"*50)

intervals = [
    ("octave",         2, 1, 1200),
    ("perfect fifth",  3, 2, 702),
    ("perfect fourth", 4, 3, 498),
    ("major third",    5, 4, 386),
    ("minor third",    6, 5, 316),
    ("major second",   9, 8, 204),
    ("minor second",  16, 15, 112),
]

for name, p, q, cents in intervals:
    rap = log(p/q)/2
    # Cents = 1200 * log2(p/q)
    actual_cents = 1200 * log(p/q) / log(2)
    print(f"  {name:16s}  {p}/{q:2d}   {rap:10.6f}   {actual_cents:7.1f}")

print()
print("  The rapidity of a musical interval = (1/2)*ln(frequency ratio).")
print("  STACKING intervals = ADDING rapidities.")
print("  Fifth + fourth = octave: ln(3/2)/2 + ln(4/3)/2 = ln(2)/2. CHECK!")
print(f"  Verify: {log(3/2)/2:.6f} + {log(4/3)/2:.6f} = {log(3/2)/2 + log(4/3)/2:.6f} = ln(2)/2 = {log(2)/2:.6f}")
print()
print("  Major third + minor third = fifth:")
print(f"  {log(5/4)/2:.6f} + {log(6/5)/2:.6f} = {log(5/4)/2+log(6/5)/2:.6f} = ln(3/2)/2 = {log(3/2)/2:.6f}")
print()
print("  CONSONANCE = LOW RAPIDITY.")
print("  The most consonant intervals have the smallest rapidities.")
print("  Rapidity IS the measure of musical tension.")
print()

# ============================================================
print("="*60)
print("3. RAPIDITY OF THE CONSTANTS")
print("="*60)
print()

constants = [
    ("1/phi (golden)",   1/phi_gr,   "arctanh(1/phi)"),
    ("1/tau (tribonacci)", 1/tau,    "arctanh(1/tau)"),
    ("1/e",              1/exp(1),   "arctanh(1/e)"),
    ("1/pi",             1/pi,       "arctanh(1/pi)"),
    ("1/sqrt(2)",        1/sqrt(2),  "arctanh(1/sqrt(2))"),
    ("1/sqrt(5)",        1/sqrt(5),  "arctanh(1/sqrt(5))"),
    ("1/3 (octave)",     1/3,        "arctanh(1/3)"),
    ("1/5 (fifth)",      1/5,        "arctanh(1/5)"),
    ("1/7 (fourth)",     1/7,        "arctanh(1/7)"),
    ("3/5 (ln2 rapidity)", 3/5,      "arctanh(3/5)"),
    ("ln(2)",            log(2),     "arctanh(ln2)"),
    ("pi/4",             pi/4,       "arctanh(pi/4)"),
]

print("  Input x          Rapidity arctanh(x)     Q(x) = e^(2*rapidity)")
print("  " + "-"*60)
for name, x, formula in constants:
    rap = atanh(x)
    qval = exp(2*rap)
    print(f"  {name:20s}   {rap:12.8f}        {qval:12.6f}")

print()
# Golden ratio rapidity
rap_phi = atanh(1/phi_gr)
print(f"  arctanh(1/phi) = {rap_phi:.10f}")
print(f"  = (3/2)*ln(phi) = {1.5*log(phi_gr):.10f}")
print(f"  Match: {abs(rap_phi - 1.5*log(phi_gr)) < 1e-10}")
print(f"  Because Q(1/phi) = phi^3, so rapidity = ln(phi^3)/2 = (3/2)*ln(phi).")
print()

# Tribonacci rapidity
rap_tau = atanh(1/tau)
print(f"  arctanh(1/tau) = {rap_tau:.10f}")
print(f"  = ln(tau) = {log(tau):.10f}")
print(f"  Match: {abs(rap_tau - log(tau)) < 1e-10}")
print(f"  Because Q(1/tau) = tau^2, so rapidity = ln(tau^2)/2 = ln(tau).")
print()

print("  THE RAPIDITY OF 1/phi IS (3/2)*ln(phi).")
print("  THE RAPIDITY OF 1/tau IS ln(tau).")
print("  The golden ratio gets a factor of 3/2 (the PERFECT FIFTH!).")
print("  The tribonacci gets a factor of 1 (the UNISON).")
print()
print("  In musical terms:")
print("  phi-rapidity = ln(phi) transposed up a fifth")
print("  tau-rapidity = ln(tau) at the fundamental")
print()

# ============================================================
print("="*60)
print("4. RAPIDITY ARITHMETIC: THE ADDITIVE STRUCTURE")
print("="*60)
print()
print("  Rapidity adds. This means we can do ARITHMETIC in rapidity space.")
print()
print("  rapidity(2) = ln(2)/2 = 0.346574")
print("  rapidity(3) = ln(3)/2 = 0.549306")
print("  rapidity(6) = ln(6)/2 = 0.895880 = rapidity(2) + rapidity(3)")
print(f"  Verify: {log(2)/2:.6f} + {log(3)/2:.6f} = {log(2)/2+log(3)/2:.6f} = {log(6)/2:.6f}")
print()
print("  MULTIPLICATION in Z maps to ADDITION in rapidity space.")
print("  This is the FUNDAMENTAL THEOREM OF ARITHMETIC restated:")
print("  rapidity(n) = sum of rapidity(p) over prime factors p of n.")
print()

# Prime factorization in rapidity space
print("  n    rapidity(n)    = sum of prime rapidities")
print("  " + "-"*55)
for n in [6, 10, 12, 15, 21, 30, 42, 105, 189, 95095]:
    rap_n = log(n)/2
    # Factor n
    temp = n
    factors = []
    d = 2
    while d*d <= temp:
        while temp % d == 0:
            factors.append(d)
            temp //= d
        d += 1
    if temp > 1:
        factors.append(temp)
    sum_raps = sum(log(p)/2 for p in factors)
    factor_str = "*".join(str(f) for f in factors)
    rap_str = " + ".join(f"ln({f})/2" for f in factors)
    print(f"  {n:6d}   {rap_n:10.6f}   = {rap_str}")

print()
print("  RAPIDITY OF 21 = ln(3)/2 + ln(7)/2 = ln(21)/2")
print(f"  = {log(21)/2:.6f}")
print("  = rapidity(curvature quantum) + rapidity(hyperbolic threshold)")
print()
print("  RAPIDITY OF 189 = ln(3^3*7)/2 = 3*ln(3)/2 + ln(7)/2")
print(f"  = {log(189)/2:.6f}")
print("  = 3*(rapidity of curvature) + rapidity of threshold")
print("  = THREE OCTAVES above the threshold!")
print()
print("  RAPIDITY OF 95095 = ln(5*7*11*13*19)/2")
print(f"  = {log(95095)/2:.6f}")
print("  = sum of rapidities of five Paley-adjacent primes")
print()

# ============================================================
print("="*60)
print("5. THE RAPIDITY SPECTRUM OF PRIMES")
print("="*60)
print()
print("  Each prime p has rapidity ln(p)/2.")
print("  The GAPS between prime rapidities:")
print()
print("  p    rapidity    gap to next    gap * 2p")
print("  " + "-"*50)
for i in range(15):
    p = primes[i]
    p_next = primes[i+1]
    rap = log(p)/2
    gap = log(p_next)/2 - log(p)/2
    # gap = ln(p_next/p)/2
    scaled = gap * 2 * p
    print(f"  {p:3d}   {rap:8.5f}   {gap:10.6f}   {scaled:8.4f}")

print()
print("  gap * 2p approaches 2 as p -> inf (by PNT: p_{n+1}/p_n -> 1).")
print("  ln(p_{n+1}/p_n) ~ (p_{n+1}-p_n)/p_n ~ 2/p_n (average gap ~ ln(p)).")
print()
print("  TWIN PRIMES in rapidity space:")
print("  Twin primes (p, p+2) have rapidity gap = ln(1+2/p)/2 ~ 1/p.")
for p in [3, 5, 11, 29, 41, 101]:
    if is_prime(p) and is_prime(p+2):
        gap = log(1+2/p)/2
        print(f"    ({p},{p+2}): rapidity gap = {gap:.8f} ~ 1/{p} = {1/p:.8f}")

print()
print("  Twin primes are points that are VANISHINGLY CLOSE in rapidity space.")
print("  The twin prime conjecture: there are infinitely many rapidity-near pairs.")
print()

# ============================================================
print("="*60)
print("6. RAPIDITY AND INFORMATION THEORY")
print("="*60)
print()
print("  The Fisher-Rao metric on Bernoulli distributions is ds^2 = dp^2/(p(1-p)).")
print("  The geodesic distance between p and q is |2*arctanh(2p-1) - 2*arctanh(2q-1)|.")
print("  Reparametrize: theta = 2*arctanh(2p-1), so p = (1+tanh(theta/2))/2.")
print("  Then ds^2 = d(theta)^2. The rapidity IS the natural parameter!")
print()
print("  A coin with bias p has Fisher-Rao position theta = 2*arctanh(2p-1).")
print("  The fair coin p=1/2 has theta = 0 (origin).")
print()
print("  Coin bias vs rapidity position:")
for p_val in [0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]:
    theta = 2*atanh(2*p_val-1)
    print(f"    p = {p_val:.2f}: theta = {theta:+8.4f}")

print()
print("  The rapidity of a biased coin is its distance from fairness")
print("  in information geometry. Maximum rapidity = certainty (p=0 or p=1).")
print()
print("  In tournaments: each arc T[i][j] is a biased coin (0 or 1).")
print("  The average rapidity of a tournament's arcs measures its")
print("  DISTANCE FROM THE RANDOM TOURNAMENT in information space.")
print()

# ============================================================
print("="*60)
print("7. THE RAPIDITY OF THE TRANSFER MATRIX EIGENVALUES")
print("="*60)
print()
print("  Transfer matrix M(x) char poly: lambda^3 - lambda^2 - x*lambda - x = 0")
print("  At x = tanh(phi), the eigenvalues live in rapidity space.")
print()
print("  The EP eigenvalue d = 1/phi has Q(d) = Q(1/phi) = phi^3.")
print(f"  EP rapidity = arctanh(1/phi) = {atanh(1/phi_gr):.10f}")
print(f"  = (3/2)*ln(phi) = {1.5*log(phi_gr):.10f}")
print()
print("  The tribonacci eigenvalue at x=1: lambda = tau.")
print("  But tau > 1, so arctanh(tau) is complex!")
print(f"  arctanh(tau) = arctanh({tau:.6f})")
print(f"  = (1/2)*ln((1+tau)/(tau-1))")
val = (1+tau)/(tau-1)
print(f"  = (1/2)*ln({val:.6f})")
print(f"  = {log(val)/2:.6f}")
print(f"  + i*pi/2  (because tau > 1, we're past the pole)")
print()
print("  The tribonacci lives PAST the rapidity singularity at v=1 (speed of light).")
print("  Its rapidity has an imaginary part = pi/2: it's a TACHYON!")
print()
print("  In our theory: x=1 is the pole of Q. For x < 1 (subluminal),")
print("  rapidity is real. For x > 1 (superluminal), rapidity gains i*pi/2.")
print("  The tribonacci constant tau ~ 1.839 is SUPERLUMINAL.")
print("  Its rapidity = ln(tau^2)/2 + i*pi/2 = ln(tau) + i*pi/2.")
print()

# ============================================================
print("="*60)
print("8. RAPIDITY LEVELS AND THE H-SPECTRUM")
print("="*60)
print()
print("  H(T) = I(Omega(T), 2) = sum over independent sets of 2^|S|.")
print("  The rapidity of H: phi_H = ln(H)/2.")
print()
print("  For the transitive tournament: H = 1, phi_H = 0 (at rest).")
print("  For T_3 (3-cycle): H = 3, phi_H = ln(3)/2 = 0.549")
print("  For T_7 (Paley): H = 189, phi_H = ln(189)/2 = 2.620")
print("  For T_11 (Paley): H = 95095, phi_H = ln(95095)/2 = 5.728")
print()

# Compute rapidity for known H values and forbidden values
print("  H-value    rapidity     notes")
print("  " + "-"*50)
special_H = [
    (1, "transitive (at rest)"),
    (3, "T_3 (3-cycle)"),
    (5, "n=4 max"),
    (7, "FORBIDDEN"),
    (9, "n=5 Paley adjacency"),
    (11, "n=5"),
    (13, "n=5"),
    (15, "n=5 cyclic"),
    (21, "FORBIDDEN"),
    (45, "n=6 max"),
    (189, "T_7 Paley max"),
    (661, "n=8 max"),
    (95095, "T_11 Paley max"),
]
for h, note in special_H:
    rap = log(h)/2
    print(f"  {h:6d}     {rap:8.5f}     {note}")

print()
print("  The forbidden rapidities:")
print(f"  rapidity(7)  = ln(7)/2  = {log(7)/2:.6f}")
print(f"  rapidity(21) = ln(21)/2 = {log(21)/2:.6f}")
print(f"  rapidity(21) - rapidity(7) = ln(3)/2 = {log(3)/2:.6f} = rapidity(3)")
print()
print("  The GAP between the two forbidden rapidities is EXACTLY")
print("  the rapidity of 3 = the curvature quantum = one OCTAVE!")
print("  H=7 and H=21 are separated by exactly one octave in rapidity space.")
print()

# ============================================================
print("="*60)
print("9. THE HARMONIC SERIES IN RAPIDITY SPACE")
print("="*60)
print()
print("  arctanh(x) = x + x^3/3 + x^5/5 + x^7/7 + ...")
print("  Each term x^k/k contributes rapidity x^k/k to the total.")
print()
print("  At x = 1/(2n+1) (the musical inputs):")
print("  arctanh(1/(2n+1)) = sum_{j=0}^inf 1/((2j+1)*(2n+1)^{2j+1})")
print("  The leading term 1/(2n+1) dominates. Higher terms decay as 1/(2n+1)^3.")
print()
print("  But for the FULL arctanh, the terms form a series of rapidities:")
print("  x + x^3/3 + x^5/5 + ...")
print("  = rapidity(fundamental) + rapidity(3rd harmonic)/3 + rapidity(5th)/5 + ...")
print()
print("  This is a FOURIER-LIKE decomposition of rapidity into harmonics!")
print("  arctanh = sum of ODD HARMONICS of the rapidity.")
print()
print("  In music: the odd harmonics 1, 3, 5, 7, ... are the overtone series")
print("  of a CLOSED PIPE (clarinet, organ stopped pipe).")
print("  Open pipes have ALL harmonics. Closed pipes have only ODD ones.")
print()
print("  arctanh IS the closed-pipe harmonic series.")
print("  artanh(x) = clarinet playing note x.")
print("  tanh(x) = the inverse: which note produces this timbre?")
print()

# ============================================================
print("="*60)
print("10. RAPIDITY COMPOSITION = VELOCITY ADDITION")
print("="*60)
print()
print("  v1 (+) v2 = (v1+v2)/(1+v1*v2) (Einstein velocity addition)")
print("  In rapidity: phi1 + phi2 (just add!)")
print()
print("  What velocities compose to give special values?")
print()

# Table: v1 (+) v2 for various musical/special velocities
special_v = [
    ("1/3", 1/3),
    ("1/5", 1/5),
    ("1/7", 1/7),
    ("3/5", 3/5),
    ("1/phi", 1/phi_gr),
]

print("  v1 (+) v2 table (relativistic addition):")
print(f"  {'':12s}", end="")
for name, _ in special_v:
    print(f"  {name:10s}", end="")
print()
for n1, v1 in special_v:
    print(f"  {n1:10s}  ", end="")
    for n2, v2 in special_v:
        composed = (v1+v2)/(1+v1*v2)
        print(f"  {composed:10.6f}", end="")
    print()

print()
print("  Key compositions:")
v_sum = (1/3 + 1/5)/(1+1/15)
print(f"  1/3 (+) 1/5 = {v_sum:.10f} = {Fraction(1,3+5).limit_denominator(100)} ???")
# (1/3+1/5)/(1+1/15) = (8/15)/(16/15) = 8/16 = 1/2
print(f"  = (8/15)/(16/15) = 1/2 EXACTLY!")
print()
print("  THE OCTAVE VELOCITY COMPOSED WITH THE FIFTH VELOCITY")
print("  GIVES v = 1/2 (rapidity = arctanh(1/2) = ln(3)/2 = rapidity of 3).")
print(f"  Check: arctanh(1/3) + arctanh(1/5) = {atanh(1/3)+atanh(1/5):.10f}")
print(f"  arctanh(1/2) = {atanh(1/2):.10f}")
print(f"  Match: {abs(atanh(1/3)+atanh(1/5)-atanh(1/2)) < 1e-10}")
print()

# More compositions
v_sum2 = (1/3 + 1/7)/(1+1/21)
print(f"  1/3 (+) 1/7 = {v_sum2:.10f}")
print(f"  = (10/21)/(22/21) = 10/22 = 5/11 EXACTLY!")
v_sum3 = (1/5 + 1/7)/(1+1/35)
print(f"  1/5 (+) 1/7 = {v_sum3:.10f}")
print(f"  = (12/35)/(36/35) = 12/36 = 1/3 EXACTLY!")
print()
print("  1/5 (+) 1/7 = 1/3!")
print("  THE FIFTH COMPOSED WITH THE FOURTH GIVES THE OCTAVE.")
print("  This is the musical identity: fifth + fourth = octave.")
print("  But now it's a RELATIVISTIC velocity addition identity!")
print()
print("  1/3 (+) 1/7 = 5/11.")
print("  The octave composed with the fourth gives 5/11.")
print("  5 and 11 are PALEY PRIMES.")
print()

# General: 1/(2a+1) (+) 1/(2b+1)
print("  GENERAL FORMULA:")
print("  1/(2a+1) (+) 1/(2b+1) = (2a+2b+2) / (2(2ab+a+b+1))")
print("  = (a+b+1) / (2ab+a+b+1)")
print()
# Verify
for a, b in [(1,2), (1,3), (2,3), (1,4), (2,4), (3,4)]:
    v1 = 1/(2*a+1)
    v2 = 1/(2*b+1)
    composed = (v1+v2)/(1+v1*v2)
    formula = (a+b+1)/(2*a*b+a+b+1)
    f = Fraction(a+b+1, 2*a*b+a+b+1)
    print(f"    1/{2*a+1} (+) 1/{2*b+1} = {f} = {float(f):.6f}   (verify: {abs(composed-float(f))<1e-10})")

print()
print("  When does 1/(2a+1) (+) 1/(2b+1) = 1/(2c+1) for integer c?")
print("  Need: (a+b+1)/(2ab+a+b+1) = 1/(2c+1)")
print("  => 2c+1 = (2ab+a+b+1)/(a+b+1)")
print("  => (2ab+a+b+1) = (2c+1)(a+b+1)")
print("  => 2ab+a+b+1 = 2ac+2bc+a+b+2c+1")
print("  => 2ab = 2c(a+b+1)")
print("  => c = ab/(a+b+1)")
print()
print("  c must be a positive integer. So (a+b+1) | ab.")
print("  Check: a=2, b=3: c = 6/6 = 1. YES! 1/5 (+) 1/7 = 1/3.")
print("  Check: a=1, b=2: c = 2/4 = 1/2. NO (not integer).")
print("  Check: a=1, b=3: c = 3/5. NO.")
print("  Check: a=3, b=5: c = 15/9 = 5/3. NO.")
print("  Check: a=5, b=9: c = 45/15 = 3. YES! 1/11 (+) 1/19 = 1/7.")
print()
print("  1/11 (+) 1/19 = 1/7!")
print("  The MINOR THIRD composed with the UNDECIMAL gives the FOURTH.")
print("  And 11, 19 are Paley primes; 7 is the threshold.")
print()
v_check = (1/11 + 1/19)/(1 + 1/(11*19))
print(f"  Verify: 1/11 (+) 1/19 = {v_check:.10f} = 1/7 = {1/7:.10f}")
print(f"  Match: {abs(v_check - 1/7) < 1e-10}")
print()

# Find ALL pairs (a,b) where c = ab/(a+b+1) is integer, with a <= b <= 20
print("  All rapidity-closed compositions 1/(2a+1) (+) 1/(2b+1) = 1/(2c+1):")
print("  (with a <= b <= 20)")
for a in range(1, 21):
    for b in range(a, 21):
        num = a*b
        den = a+b+1
        if num % den == 0:
            c = num // den
            print(f"    1/{2*a+1} (+) 1/{2*b+1} = 1/{2*c+1}  (a={a}, b={b}, c={c})")

print()

# ============================================================
print("="*60)
print("11. THE LORENTZ GROUP AND TOURNAMENT SYMMETRY")
print("="*60)
print()
print("  A Lorentz boost with rapidity phi acts as:")
print("  t' = t*cosh(phi) + x*sinh(phi)")
print("  x' = t*sinh(phi) + x*cosh(phi)")
print()
print("  For rapidity phi = ln(n)/2 (the rapidity of integer n):")
print(f"  cosh(ln(n)/2) = (n + 1/n)/2 = (n^2+1)/(2n)")
print(f"  sinh(ln(n)/2) = (n - 1/n)/2 = (n^2-1)/(2n)")
print()
print("  For n=2:")
print(f"  cosh = 5/4 = {(4+1)/4:.4f}, sinh = 3/4 = {(4-1)/4:.4f}")
print(f"  tanh = 3/5 = {3/5:.4f} ... wait, tanh(ln(2)/2) = (2-1)/(2+1) = 1/3")
print(f"  Verify: tanh(ln(2)/2) = {tanh(log(2)/2):.10f} = 1/3")
print()
print("  The velocity corresponding to 'being boosted by n=2' is v = 1/3.")
print("  Q(1/3) = 2 (the binary). The Lorentz boost by the binary")
print("  corresponds to velocity 1/3 = one third the speed of light.")
print()

# The boost matrix for n
print("  Boost matrices for small n (in rapidity ln(n)/2):")
for n in [2, 3, 5, 7]:
    c = (n + 1/n)/2
    s = (n - 1/n)/2
    v = s/c
    print(f"    n={n}: cosh={c:.4f}, sinh={s:.4f}, v=tanh={(n-1)/(n+1):.6f} = (n-1)/(n+1)")

print()
print("  The velocity for boost-n is (n-1)/(n+1) = the Cayley ADDRESS of n.")
print("  This is NOT a coincidence. It's the DEFINITION:")
print("  x_n = (n-1)/(n+1) = tanh(ln(n)/2).")
print("  The Cayley address IS the velocity. The natural number IS the boost.")
print()

# ============================================================
print("="*60)
print("12. RAPIDITY AND THE FORBIDDEN NUMBERS")
print("="*60)
print()
print("  H = 7: rapidity = ln(7)/2 = {:.6f}".format(log(7)/2))
print("  H = 21: rapidity = ln(21)/2 = {:.6f}".format(log(21)/2))
print()
print("  The gap: ln(21)/2 - ln(7)/2 = ln(3)/2 = {:.6f}".format(log(3)/2))
print("  = rapidity of 3 = one OCTAVE of musical rapidity")
print()
print("  The SUM: ln(7)/2 + ln(21)/2 = ln(147)/2 = {:.6f}".format(log(147)/2))
print("  147 = 3 * 7^2. rapidity(7) + rapidity(21) = rapidity(147).")
print()
print("  The average: (ln(7)+ln(21))/4 = ln(sqrt(147))/2 = {:.6f}".format(log(sqrt(147))/2))
print(f"  sqrt(147) = 7*sqrt(3) = {7*sqrt(3):.6f}")
print()
print("  The forbidden numbers sit at rapidities ln(7)/2 and ln(7)/2 + ln(3)/2.")
print("  Their MIDPOINT in rapidity space is at n = sqrt(147) = 7*sqrt(3).")
print("  And sqrt(3) = 2*sin(60 deg) = the flat-plane scale factor.")
print()

# Velocity of forbidden numbers
v7 = 6/8  # address of 7 = (7-1)/(7+1) = 6/8 = 3/4
v21 = 20/22  # address of 21 = 20/22 = 10/11
print(f"  Velocity of H=7:  v = (7-1)/(7+1) = 3/4 = {3/4:.6f}")
print(f"  Velocity of H=21: v = (21-1)/(21+1) = 10/11 = {10/11:.6f}")
print()
print(f"  3/4 (+) 10/11 = ?")
v_comp = (3/4 + 10/11)/(1 + 3*10/(4*11))
print(f"  = {v_comp:.10f}")
# (3/4+10/11)/(1+30/44) = (33+40)/(44) / (74/44) = 73/74
print(f"  = 73/74")
print(f"  Q(73/74) = (1+73/74)/(1-73/74) = 147 = 3*7^2")
print(f"  Verify: {(1+73/74)/(1-73/74):.1f}")
print()
print("  THE FORBIDDEN VELOCITIES COMPOSED RELATIVISTICALLY GIVE 73/74,")
print("  WHICH HAS Q-VALUE 147 = 3*7^2 = product of forbidden factors.")
print()

# ============================================================
print("="*60)
print("13. RAPIDITY AS HYPERBOLIC DISTANCE")
print("="*60)
print()
print("  On the Poincare disk model, the hyperbolic distance from 0 to r is:")
print("  d(0, r) = 2*arctanh(r) = 2*rapidity(r).")
print()
print("  The tournament lives on the interval [0,1) which IS the Poincare radius.")
print("  Each tournament parameter x has hyperbolic distance 2*arctanh(x) from 0.")
print()
print("  Hyperbolic distances of special points:")
special_points = [
    ("1/3 (octave)", 1/3),
    ("1/2", 1/2),
    ("3/5 (ln2)", 3/5),
    ("1/phi", 1/phi_gr),
    ("3/4 (H=7 addr)", 3/4),
    ("4/5", 4/5),
    ("10/11 (H=21 addr)", 10/11),
    ("EP x_c = 8-5phi", 8-5*phi_gr),
]
for name, r in special_points:
    if 0 < r < 1:
        d = 2*atanh(r)
        print(f"    {name:25s}: d = {d:.6f} = {d/log(2):.4f} * ln(2)")

print()
print("  The EP (exceptional point) distance:")
x_ep = 8 - 5*phi_gr
print(f"  x_EP = 8 - 5*phi = {x_ep:.10f}")
print(f"  d(0, x_EP) = 2*arctanh(x_EP) = {2*atanh(x_ep):.10f}")
print(f"  = ln(Q(x_EP)) = ln(phi^3 * something)...")
q_ep = (1+x_ep)/(1-x_ep)
print(f"  Q(x_EP) = {q_ep:.10f}")
print(f"  ln(Q(x_EP)) = {log(q_ep):.10f}")
print()

# ============================================================
print("="*60)
print("GRAND SUMMARY: RAPIDITY IS EVERYTHING")
print("="*60)
print()
print("  DOMAIN          RAPIDITY IS                    FORMULA")
print("  " + "-"*60)
print("  Relativity      boost parameter                arctanh(v)")
print("  Music           consonance measure             ln(interval)/2")
print("  Number theory   half-logarithm                 ln(n)/2")
print("  Information     Fisher-Rao distance            2*arctanh(2p-1)")
print("  Hyperbolic      Poincare disk distance         2*arctanh(r)")
print("  Tournament      coupling strength              arctanh(x)")
print("  Transfer mat    eigenvalue log                 ln(lambda)/2")
print("  Fibonacci       skip-3 parameter               ln(F_{n+2}/F_{n-1})/2")
print()
print("  All of these are the SAME arctanh, seen through different lenses.")
print("  Rapidity is the universal coordinate of the Cayley transform.")
print("  It linearizes multiplication, composition, and coupling.")
print("  It measures distance, consonance, and information simultaneously.")
print()
print("  And the forbidden numbers {7, 21} are separated by exactly")
print("  one octave of rapidity: ln(3)/2 = the rapidity of the")
print("  curvature quantum. The theory's impossibilities are musical.")
