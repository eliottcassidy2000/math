#!/usr/bin/env python3
"""leading_digits_s116h.py — Powers of 2, leading digits, and the geometry of 'good enough'.

1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024

Leading digits: 1, 2, 4, 8, 1, 3, 6, 1, 2, 5, 1

The leading digit CYCLES through {1,2,4,8,1,3,6,1,2,5,1,...}
with a quasi-period of ~10/3 = 3.32 steps per decade.

Because log10(2) = 0.30103... is irrational, the leading digit
NEVER exactly repeats. But it ALMOST repeats every 10 steps
(2^10 = 1024 ~ 1000). The error: 1024/1000 = 1.024.

This is Benford's law. This is the Pythagorean comma.
This is the rapidity line wrapping around the decimal circle.
"""
from math import log, log10, log2, floor, pi, exp, atanh
from fractions import Fraction

print()
print("  POWERS OF 2, LEADING DIGITS, AND THE GEOMETRY OF GOOD ENOUGH")
print()
print("="*70)
print()

# ============================================================
print("  THE TABLE")
print("  " + "-"*40)
print()
print("  n    2^n          leading   fractional     position on")
print("                    digit     part of log    decimal circle")
print("  " + "-"*65)
for n in range(0, 35):
    val = 2**n
    ld = int(str(val)[0])
    frac = (n * log10(2)) % 1.0
    print(f"  {n:2d}   {val:>12d}     {ld}        {frac:.5f}        {'|' + '#'*int(frac*40)}")

print()

# ============================================================
print()
print("  THE OBSERVATION")
print("  " + "-"*40)
print()
print("  2^0  = 1     leading digit 1")
print("  2^1  = 2                   2")
print("  2^2  = 4                   4")
print("  2^3  = 8                   8")
print("  2^4  = 16                  1  <-- back to 1, but 16 not 10")
print("  2^5  = 32                  3")
print("  2^6  = 64                  6")
print("  2^7  = 128                 1  <-- back to 1 again")
print("  2^8  = 256                 2")
print("  2^9  = 512                 5")
print("  2^10 = 1024                1  <-- back to 1 again")
print()
print("  Gaps between leading-digit-1 appearances: 4, 3, 3, 4, 3, 3, 4, 3, 3, ...")
print()
gaps = []
last = 0
for n in range(1, 100):
    if int(str(2**n)[0]) == 1:
        gaps.append(n - last)
        last = n
print(f"  Gaps: {gaps[:20]}")
print()
print("  The pattern is {4, 3, 3} repeating!")
print("  4 + 3 + 3 = 10 steps per full decade.")
print("  10 steps for 2^10 = 1024 ~ 10^3.")
print()
print("  WHY {4, 3, 3}?")
print("  Because log10(2) = 0.30103...")
print("  To get from leading digit 1 back to 1, need the fractional part")
print("  of n*log10(2) to pass through 1 (mod 1).")
print("  1/0.30103 = 3.322... steps per revolution.")
print("  So approximately every 3.32 steps, the leading digit returns to 1.")
print("  3.32 ~ 10/3. The actual pattern: 3, 3, 4, 3, 3, 4, ...")
print("  = two shorts (3) and one long (4), repeating.")
print()

# ============================================================
print()
print("  THE GEOMETRY: A CIRCLE WITH A STEP")
print("  " + "-"*40)
print()
print("  The fractional part {n * log10(2)} traces a path on the")
print("  unit circle (identified with [0, 1)).")
print("  Each step advances by log10(2) = 0.30103.")
print()
print("  After 3 steps: 3 * 0.30103 = 0.90309. Almost 1. Short by 0.097.")
print("  After 4 steps: 4 * 0.30103 = 1.20412. Past 1. Overshoot by 0.204.")
print()
print("  The DECISION at each '1': do we take 3 steps or 4?")
print("  If we've accumulated DEFICIT (we're behind), take 4 (catch up).")
print("  If we've accumulated SURPLUS (we're ahead), take 3 (coast).")
print()
print("  The deficit/surplus after each cycle:")
deficit = 0
step_log = log10(2)
print("  Cycle  Steps  End position   Deficit from integer")
for cycle in range(10):
    # Try 3 steps
    pos_3 = deficit + 3 * step_log
    pos_4 = deficit + 4 * step_log
    err_3 = 1.0 - pos_3  # how far short of 1
    err_4 = pos_4 - 1.0  # how far past 1

    if abs(err_3) < abs(err_4):
        # 3 steps is closer to integer
        steps = 3
        deficit = pos_3 - 1.0  # will be negative (short)
    else:
        steps = 4
        deficit = pos_4 - 1.0  # will be positive (past)

    print(f"  {cycle:3d}     {steps}      {deficit+1:.5f}       {deficit:+.5f}")

print()

# ============================================================
print()
print("  THE CONTINUED FRACTION OF log10(2)")
print("  " + "-"*40)
print()
print("  log10(2) = 0.30103... = [0; 3, 3, 9, 2, 1, 1, 2, 1, 3, 1, 18, ...]")
print()
# Compute CF of log10(2)
x = log10(2)
cf = []
for _ in range(12):
    a = int(x)
    cf.append(a)
    x = x - a
    if x < 1e-10:
        break
    x = 1/x
print(f"  CF = {cf}")
print()
print("  The first partial quotients: 0, 3, 3, 9, 2, ...")
print("  The convergents:")
# Compute convergents
p_prev, p_curr = 1, 0
q_prev, q_curr = 0, 1
for i, a in enumerate(cf):
    p_prev, p_curr = p_curr, a * p_curr + p_prev
    q_prev, q_curr = q_curr, a * q_curr + q_prev
    if i > 0:
        approx = p_curr / q_curr
        err = abs(approx - log10(2))
        print(f"  [{', '.join(str(c) for c in cf[:i+1])}] = {p_curr}/{q_curr} = {approx:.10f}  err = {err:.2e}")

print()
print("  The BEST rational approximations to log10(2):")
print("  3/10 = 0.3     (10 steps ~ 3 decades. 2^10 = 1024 ~ 10^3)")
print("  0.30103... exact")
print()
print("  The approximation 3/10 says: 10 doublings ~ 3 decades.")
print("  The error: 10*log10(2) - 3 = 0.0103.")
print("  2^10 = 1024 = 1000 * 1.024. The 2.4% error.")
print()
print("  The NEXT approximation: log10(2) ~ 59/196.")
print("  But the USEFUL one is 3/10: good enough for engineering.")
print("  '2^10 ~ 10^3' is the GOOD ENOUGH approximation.")
print()

# ============================================================
print()
print("  GOOD ENOUGH: THE 3-3-4 PATTERN AS A BEAT FREQUENCY")
print("  " + "-"*40)
print()
print("  Two frequencies: the DOUBLING (period 1 in log-space)")
print("  and the DECADE (period 1/log10(2) = 3.322 in log-space).")
print()
print("  They are incommensurable (log10(2) is irrational).")
print("  The BEAT between them has a quasi-period.")
print()
print("  The beat pattern 3-3-4 is the SIMPLEST approximation:")
print("  3+3+4 = 10 doublings per 3 decades.")
print("  Average: 10/3 = 3.333... doublings per decade.")
print("  Actual: 1/log10(2) = 3.322... doublings per decade.")
print("  Error: 0.011/3.322 = 0.3%.")
print()
print("  The 3-3-4 pattern is 'good enough' because:")
print("  - It uses only TWO step sizes (3 and 4)")
print("  - It repeats every 10 steps")
print("  - The error is < 1%")
print()
print("  This is EXACTLY the same as equal temperament in music!")
print("  12 semitones approximate the octave (frequency ratio 2).")
print("  7 semitones approximate the fifth (3/2).")
print("  The error: the Pythagorean comma, ~1.4%.")
print()
print("  The 3-3-4 pattern IS a temperament for log10(2).")
print("  It tempers the decade into {3,3,4} doublings.")
print("  Just as 12-TET tempers the octave into 12 equal semitones.")
print()

# ============================================================
print()
print("  THE GEOMETRY THAT MOVES WITH THIS PERIODICITY")
print("  " + "-"*40)
print()
print("  A point on a circle, advancing by angle 2*pi*log10(2) per step.")
print("  = 2*pi * 0.30103 = 1.8912 radians = 108.4 degrees per step.")
print()
print("  108 degrees is close to 360/3.333 = 108 degrees exactly.")
print("  And 108 degrees = interior angle of a regular PENTAGON!")
print()
angle = 360 * log10(2)
print(f"  Angle per step: {angle:.2f} degrees")
print(f"  Pentagon interior angle: 108.00 degrees")
print(f"  Difference: {angle - 108:.2f} degrees")
print()
print("  The step angle is 0.4 degrees MORE than the pentagon angle.")
print("  So the geometry is a SLIGHTLY OVER-ROTATING pentagon.")
print("  After 5 steps: 5 * 108.37 = 541.85 degrees = 360 + 181.85.")
print("  That's 1.5 full turns + 1.85 degrees.")
print()
print("  After 10 steps: 10 * 108.37 = 1083.7 degrees = 3 * 360 + 3.7.")
print("  Almost exactly 3 full turns. Off by 3.7 degrees.")
print()
print("  After 100 steps: 100 * 108.37 = 10837 degrees = 30 * 360 + 37.")
print("  30.1 turns. Still drifting.")
print()
print("  The SHAPE: a star polygon {10/3}.")
print("  Connect every 3rd vertex of a regular 10-gon.")
print("  This traces a star that visits 10 vertices in order")
print("  0, 3, 6, 9, 2, 5, 8, 1, 4, 7 (mod 10).")
print("  The step is 3 out of 10 = 0.3 of the circle = log10(2) approximately!")
print()
print("  The {10/3} star polygon IS the geometry of powers of 2")
print("  on the decimal circle.")
print()

# Draw it
print("  The star {10/3}:")
print("  Step 0: position 0 (digit 1)")
print("  Step 1: position 3 (digit 2... no, 8)")
# Actually the positions aren't digit positions, they're on the circle
# Let me think about this differently.
# The fractional part of n*log10(2) visits positions on [0,1)
# After 10 steps it returns to near 0 (within 0.01)
# The 10 positions visited are approximately 0/10, 3/10, 6/10, 9/10, 2/10, 5/10, 8/10, 1/10, 4/10, 7/10
# i.e., {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}/10 in the order 0,3,6,9,2,5,8,1,4,7
# This is multiplication by 3 mod 10!
print()
print("  The 10 positions (approximately):")
for n in range(10):
    pos = (n * log10(2)) % 1.0
    approx_tenth = round(pos * 10) % 10
    digit = int(str(2**n)[0])
    print(f"  n={n}: frac = {pos:.4f} ~ {approx_tenth}/10, leading digit = {digit}")

print()
print("  The positions follow the pattern 0,3,6,9,2,5,8,1,4,7 (mod 10).")
print("  This IS multiplication by 3 modulo 10.")
print("  The multiplicative group Z/10Z* = {1,3,7,9} has order 4.")
print("  But 3 generates all of Z/10Z under addition (gcd(3,10)=1).")
print()
print("  So the leading digit pattern of 2^n is controlled by the")
print("  additive structure of 3 in Z/10Z.")
print("  And 3 in Z/10Z is the CURVATURE QUANTUM modulo the decimal base.")
print()

# ============================================================
print()
print("  THE RAPIDITY CONNECTION")
print("  " + "-"*40)
print()
print("  log10(2) = ln(2)/ln(10) = (2*rapidity(2))/(2*rapidity(10))")
print("  = rapidity(2)/rapidity(10)")
print()
r2 = log(2)/2
r10 = log(10)/2
print(f"  rapidity(2) = {r2:.6f}")
print(f"  rapidity(10) = {r10:.6f}")
print(f"  ratio = {r2/r10:.6f} = log10(2)")
print()
print("  log10(2) IS the rapidity of 2 measured in units of rapidity(10).")
print("  The decimal system uses rapidity(10) as its unit.")
print("  The binary system uses rapidity(2) as its unit.")
print()
print("  Their ratio is irrational: decimal and binary are INCOMMENSURABLE.")
print("  The 3-3-4 pattern is the best 'temperament' reconciling them.")
print()
print("  In the Cayley framework:")
print("  Q(address of 2) = 2. address = 1/3.")
print("  Q(address of 10) = 10. address = 9/11.")
print()
print("  1/3 and 9/11: these are the Cayley addresses of 2 and 10.")
print("  Their rapidity ratio log10(2) determines the leading digit pattern.")
print("  The ENTIRE structure of 'how binary meets decimal'")
print("  is encoded in the ratio of two Cayley addresses.")
print()

# ============================================================
print()
print("  GOOD ENOUGH AND THE PHILOSOPHY OF APPROXIMATION")
print("  " + "-"*40)
print()
print("  The exact value: log10(2) = 0.3010299957...")
print("  The 'good enough' value: 3/10 = 0.3.")
print("  Error: 0.34%.")
print()
print("  This 0.34% error is WHY:")
print("  - 2^10 = 1024, not 1000 (2.4% overshoot)")
print("  - A kilobyte is 1024 bytes, not 1000")
print("  - The SI/IEC prefix confusion (kB vs KiB)")
print("  - The leading digit returns to 1 every 3 or 4 steps, not exactly 3.32")
print()
print("  The 'good enough' approximation 3/10 is the FIRST convergent")
print("  of the continued fraction of log10(2).")
print("  It is the BEST rational with denominator <= 10.")
print()
print("  The NEXT 'good enough': 59/196 ~ 0.30102.")
print("  Error: 0.003%. But denominator 196 is too big to be useful.")
print()
print("  In practice, 3/10 is WHERE EVERYONE STOPS.")
print("  'A kilobyte is about a thousand bytes.'")
print("  'Ten doublings is about three orders of magnitude.'")
print("  'The leading digit returns to 1 every 3-4 steps.'")
print()
print("  GOOD ENOUGH is the first convergent.")
print("  It is the point where the continued fraction FIRST gives")
print("  an approximation within human tolerance (~1%).")
print("  Beyond that, improvement is academic.")
print()

# ============================================================
print()
print("  THE STAR {10/3} AND THE PENTAGON {5/2}")
print("  " + "-"*40)
print()
print("  The star polygon {10/3} visits every 3rd vertex of a 10-gon.")
print("  The pentagram {5/2} visits every 2nd vertex of a 5-gon.")
print()
print("  {10/3} = double cover of {5/3}? No.")
print("  {10/3}: step 3 in 10-gon. Visits all 10 vertices. gcd(3,10)=1.")
print("  {5/2}: step 2 in 5-gon. Visits all 5 vertices. gcd(2,5)=1.")
print()
print("  Their ratio: {10/3} has angle step 108 degrees.")
print("               {5/2} has angle step 144 degrees.")
print("  108/144 = 3/4 = the Cayley address of 7 = the forbidden number!")
print()
print("  THE RATIO OF THE TWO STAR POLYGON ANGLES IS 3/4.")
print("  THE CAYLEY ADDRESS OF THE FORBIDDEN NUMBER.")
print()
print("  And 108 degrees = pi - pi/5 = the supplement of the golden angle's cousin.")
print("  144 degrees = pi - pi/5 = ... no, 144 = 4*36 = 4*pi/5.")
print("  108 = 3*36 = 3*pi/5.")
print("  108/144 = 3/4. YES.")
print()
print("  The decimal-binary beat frequency (108 degrees)")
print("  and the pentagram angle (144 degrees)")
print("  are in the ratio 3:4 = the Cayley address of 7.")
print()
print("  Three things meet:")
print("  1. The leading digit pattern of 2^n (log10(2) ~ 3/10)")
print("  2. The pentagram / golden ratio (108 and 144 degrees)")
print("  3. The forbidden number 7 (Cayley address 3/4)")
print()
print("  They meet because 3/4 = 3*36/4*36 = 108/144,")
print("  and 36 = 360/10 = the angular step of the decimal circle,")
print("  and log10(2) ~ 3/10 = three such steps.")
print()
print("  The forbidden number lives in the ratio of the binary step")
print("  to the pentagonal step on the decimal circle.")
