#!/usr/bin/env python3
"""tanh_understood_s116h.py — Understanding tanh and arctanh through
the geometry of 'good enough'.

tanh is the function that ACCUMULATES partial bits
and NEVER REACHES 1.

arctanh is the function that COUNTS how many partial bits
you've accumulated.

The leading-digit analysis showed: log10(2) ~ 3/10 is 'good enough'.
The 3-3-4 pattern is a temperament. The pentagon angle appears.
Now: tanh and arctanh ARE this same geometry, internalized.
"""
from math import log, exp, tanh, atanh, cosh, sinh, pi, sqrt

print()
print("  UNDERSTANDING tanh AND arctanh")
print()
print("="*70)
print()

# ============================================================
print("  I. tanh AS ACCUMULATION THAT SATURATES")
print("  " + "-"*40)
print()
print("  tanh(x) for x = 0, 0.5, 1, 1.5, 2, 3, 5, 10:")
print()
print("  x       tanh(x)    1-tanh(x)    how close to 1")
print("  " + "-"*55)
for x in [0, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
    t = tanh(x)
    gap = 1 - t
    bar = "#" * int(t * 40)
    print(f"  {x:5.1f}    {t:.8f}    {gap:.2e}     |{bar}{'.' * (40 - len(bar))}|")

print()
print("  tanh APPROACHES 1 but never reaches it.")
print("  Each unit of x gets you CLOSER, but the returns DIMINISH.")
print("  The gap to 1 shrinks exponentially: 1 - tanh(x) ~ 2*e^{-2x}.")
print()
print("  This is the same as the Mersenne ladder:")
print("  1/2, 3/4, 7/8, 15/16, 31/32, ...")
print("  = tanh(ln(3)/2), tanh(ln(7)/2), tanh(ln(15)/2), ...")
print()
print("  Verify:")
for k in range(1, 7):
    x_k = 1 - 2**(-k)
    rap = atanh(x_k)
    print(f"  1 - 2^{{-{k}}} = {x_k:.5f} = tanh({rap:.5f}) = tanh(ln({2**(k+1)-1})/2)")

print()
print("  The Mersenne ladder IS tanh evaluated at half-logs of Mersenne numbers.")
print("  Each rung: tanh(ln(M_k)/2) = 1 - 2^{-(k-1)}.")
print("  The ladder climbs toward 1 but each step is HALF the previous.")
print()

# ============================================================
print()
print("  II. tanh AS 'GOOD ENOUGH' AT EVERY SCALE")
print("  " + "-"*40)
print()
print("  At x = 1: tanh(1) = 0.7616.")
print("  'Good enough' for what?")
print("  For being MOST of the way to 1. We're 76% there.")
print()
print("  At x = 2: tanh(2) = 0.9640.")
print("  Now we're 96% there. Good enough for most purposes.")
print()
print("  At x = 3: tanh(3) = 0.9951.")
print("  99.5%. Good enough for engineering.")
print()
print("  At x = 5: tanh(5) = 0.9999.")
print("  99.99%. Good enough for everyone.")
print()
print("  The 'good enough' thresholds in rapidity:")
print()
print("  Threshold    rapidity needed    Cayley Q     Musical")
print("  " + "-"*55)
thresholds = [0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 0.9999]
for t in thresholds:
    r = atanh(t)
    q = (1+t)/(1-t)
    musical = ""
    if abs(q - 3) < 0.05: musical = "~ twelfth"
    elif abs(q - 7) < 0.1: musical = "~ forbidden!"
    elif abs(q - 19) < 1: musical = "~ 19th harmonic"
    elif abs(q - 39) < 2: musical = "~ 39th harmonic"
    elif abs(q - 199) < 10: musical = "~ 199th harmonic"
    elif abs(q - 1999) < 100: musical = "~ 1999th harmonic"
    print(f"  {t:.4f}       {r:8.4f}          {q:10.1f}    {musical}")

print()
print("  To reach 75% of the way to 1: rapidity ~ 0.55 (one twelfth).")
print("  To reach 99%: rapidity ~ 2.65.")
print("  To reach 99.99%: rapidity ~ 4.95.")
print()
print("  Each additional '9' in the precision costs ~1.15 rapidity units.")
print("  = ln(10)/2 = one DECADE of rapidity.")
print("  = 3.32 octaves.")
print()
dec = log(10)/2
print(f"  One decade of rapidity = ln(10)/2 = {dec:.4f}")
print(f"  = {dec / (log(2)/2):.2f} octaves")
print()
print("  Adding a decimal digit of precision = adding 3.32 octaves of rapidity.")
print("  The 3-3-4 pattern again! 3.32 doublings per decade.")
print("  tanh embodies the 3-3-4 rhythm at EVERY scale.")
print()

# ============================================================
print()
print("  III. arctanh AS THE INVERSE: COUNTING THE DOUBLINGS")
print("  " + "-"*40)
print()
print("  arctanh(x) = number of 'doublings toward 1' contained in x.")
print()
print("  arctanh(0) = 0. No doublings yet. Fair coin.")
print("  arctanh(0.5) = 0.549. About 1.6 octaves. (Q = 3.)")
print("  arctanh(0.75) = 0.973. About 2.8 octaves. (Q = 7. Forbidden!)")
print("  arctanh(0.9) = 1.472. About 4.2 octaves. (Q = 19.)")
print("  arctanh(0.99) = 2.647. About 7.6 octaves. (Q = 199.)")
print("  arctanh(0.999) = 3.800. About 11.0 octaves. (Q = 1999.)")
print()
print("  arctanh COUNTS how many times you've halved the gap to 1.")
print("  0.5: halved once (gap went from 1 to 0.5).")
print("  0.75: halved twice-ish (gap went from 1 to 0.25).")
print("  0.875: halved three times (gap = 0.125 = 1/8).")
print()
print("  But arctanh doesn't count EXACT halvings. It counts")
print("  CONTINUOUS halvings. The function is smooth, not stepwise.")
print()
print("  arctanh(1 - 2^{-k}) = ln(2^{k+1} - 1)/2 ~ (k+1)*ln(2)/2 ~ k * octave.")
print("  At the Mersenne points, arctanh ~ k octaves. Between them, it interpolates.")
print()

# ============================================================
print()
print("  IV. THE TAYLOR SERIES AS FRACTIONAL HALVINGS")
print("  " + "-"*40)
print()
print("  arctanh(x) = x + x^3/3 + x^5/5 + x^7/7 + ...")
print()
print("  Each term x^{2k+1}/(2k+1) is a FRACTIONAL HALVING.")
print("  The first term x gives you MOST of the way.")
print("  The correction x^3/3 gives you most of the REST.")
print("  The correction x^5/5 gives you most of what remains.")
print("  And so on.")
print()
print("  For x = 0.5:")
x = 0.5
exact = atanh(x)
cumul = 0
print(f"  Term     Value      Cumulative  Fraction of exact ({exact:.6f})")
print("  " + "-"*60)
for k in range(8):
    term = x**(2*k+1) / (2*k+1)
    cumul += term
    frac = cumul / exact
    remaining = 1 - frac
    print(f"  x^{2*k+1}/{2*k+1:<2d}   {term:10.6f}   {cumul:10.6f}    {frac:.4%}  (gap: {remaining:.4%})")

print()
print("  The first term gets you 91% of the way.")
print("  Adding the second term gets you to 99.2%.")
print("  Adding the third: 99.92%.")
print("  Each term CLOSES ~90% of the remaining gap.")
print()
print("  For x = 0.9:")
x = 0.9
exact = atanh(x)
cumul = 0
print(f"  Term     Value      Cumulative  Fraction of exact ({exact:.6f})")
print("  " + "-"*60)
for k in range(10):
    term = x**(2*k+1) / (2*k+1)
    cumul += term
    frac = cumul / exact
    print(f"  x^{2*k+1}/{2*k+1:<2d}   {term:10.6f}   {cumul:10.6f}    {frac:.4%}")

print()
print("  At x = 0.9, the convergence is SLOWER.")
print("  The first term gives only 61% of the total.")
print("  You need 7+ terms to reach 99%.")
print("  Near the pole, the partial bits are LESS effective individually")
print("  but there are MORE of them and they all matter.")
print()

# ============================================================
print()
print("  V. THE 3-3-4 RHYTHM IN tanh")
print("  " + "-"*40)
print()
print("  tanh advances by 'almost thirds' toward 1:")
print()
print("  From 0 to tanh(1) = 0.762: covered 76.2% of [0,1].")
print("  From tanh(1) to tanh(2) = 0.964: covered 84.9% of remaining gap.")
print("  From tanh(2) to tanh(3) = 0.995: covered 86.4% of remaining gap.")
print()
prev = 0
for n in range(1, 8):
    t = tanh(n)
    gap_before = 1 - prev
    gap_after = 1 - t
    covered = (gap_before - gap_after) / gap_before if gap_before > 0 else 0
    print(f"  Step {n}: tanh({n}) = {t:.6f}, "
          f"covered {covered:.1%} of remaining gap")
    prev = t

print()
print("  Each step covers ~85-97% of the remaining gap.")
print("  This is the 'good enough' principle at each level:")
print("  one more step gets you ~90% of what's left.")
print("  You asymptotically approach 1 but NEVER arrive.")
print()
print("  The NUMBER of steps to reach a given precision:")
print("  To reach 1 - epsilon: need ~ln(1/epsilon)/2 = ln(1/epsilon)/2 steps.")
print("  This is HALF the number of decimal digits of 1/epsilon.")
print()
print("  Precision 0.1: ~1.2 steps. (tanh(1.2) = 0.834)")
print("  Precision 0.01: ~2.6 steps. (tanh(2.6) = 0.989)")
print("  Precision 0.001: ~3.8 steps. (tanh(3.8) = 0.9993)")
print("  Precision 10^{-k}: ~1.15*k steps. (= k * ln(10)/2)")
print()
print("  Each DIGIT of precision costs ln(10)/2 = 1.15 units of rapidity.")
print("  = 3.32 octaves.")
print("  THE 3-3-4 PATTERN AGAIN.")
print()

# ============================================================
print()
print("  VI. tanh AND arctanh AS THE SAME GEOMETRY SEEN TWO WAYS")
print("  " + "-"*40)
print()
print("  tanh: input = rapidity (unbounded). output = velocity (bounded).")
print("  arctanh: input = velocity (bounded). output = rapidity (unbounded).")
print()
print("  tanh says: 'given this many doublings, how close to 1 am I?'")
print("  arctanh says: 'given this closeness to 1, how many doublings was it?'")
print()
print("  They are the SAME question asked in opposite directions.")
print("  tanh goes FORWARD (accumulate doublings, approach 1).")
print("  arctanh goes BACKWARD (given position, count doublings).")
print()
print("  The 'good enough' principle:")
print("  tanh: 'how much rapidity do I NEED to be good enough?'")
print("  arctanh: 'how good IS my current position?'")
print()

# ============================================================
print()
print("  VII. THE PENTAGON HIDDEN IN tanh")
print("  " + "-"*40)
print()
print("  We found: the leading digit of 2^n traces a pentagon (108 degrees).")
print("  And tanh IS the function that converts doublings to position.")
print()
print("  So tanh CONTAINS the pentagon.")
print()
print("  Where? In its relationship to the exponential:")
print("  tanh(x) = (e^{2x} - 1)/(e^{2x} + 1)")
print()
print("  At x = ln(phi)/2 (the rapidity of the golden ratio):")
x_phi = log((1+sqrt(5))/2) / 2
t_phi = tanh(x_phi)
phi = (1+sqrt(5))/2
print(f"  tanh(ln(phi)/2) = {t_phi:.10f}")
print(f"  Compare: 1/sqrt(5) = {1/sqrt(5):.10f}")
print()
print(f"  tanh(ln(phi)/2) = 1/sqrt(5)!")
print(f"  Verify: {abs(t_phi - 1/sqrt(5)) < 1e-10}")
print()
print("  PROOF: tanh(ln(phi)/2) = (phi - 1)/(phi + 1)")
print("  = (phi - 1)/(phi + 1)")
print(f"  phi - 1 = 1/phi = {1/phi:.6f}")
print(f"  phi + 1 = phi^2 = {phi**2:.6f}")
print(f"  (phi-1)/(phi+1) = (1/phi)/phi^2 = 1/phi^3 = {1/phi**3:.6f}")
print()
print("  Hmm, that gives 1/phi^3 = 0.2360. But tanh gives 0.4472 = 1/sqrt(5).")
print("  Let me recompute.")
print(f"  tanh(ln(phi)/2):")
print(f"  e^{{ln(phi)}} = phi = {phi:.6f}")
print(f"  tanh(ln(phi)/2) = (e^{{ln(phi)}} - 1)/(e^{{ln(phi)}} + 1)")
print(f"  = (phi - 1)/(phi + 1) = {(phi-1)/(phi+1):.10f}")
print(f"  = 1/(phi+1) * (phi-1) = 1/phi^2 * (phi-1) = (phi-1)/phi^2")
print(f"  = (1/phi)/phi^2 = 1/phi^3 = {1/phi**3:.10f}")
print()
print("  So tanh(ln(phi)/2) = 1/phi^3 = 0.2360..., NOT 1/sqrt(5).")
print("  My earlier computation was wrong. Let me check what gives 1/sqrt(5).")
print()
r_check = atanh(1/sqrt(5))
print(f"  arctanh(1/sqrt(5)) = {r_check:.10f}")
print(f"  ln(phi) = {log(phi):.10f}")
print(f"  ln(phi)/2 = {log(phi)/2:.10f}")
print(f"  These are NOT equal. arctanh(1/sqrt(5)) = {r_check:.6f} != ln(phi)/2 = {log(phi)/2:.6f}")
print()
# What IS arctanh(1/sqrt(5))?
print(f"  Q(1/sqrt(5)) = {(1+1/sqrt(5))/(1-1/sqrt(5)):.10f}")
print(f"  = (sqrt(5)+1)^2/(5-1) = (6+2*sqrt(5))/4 = (3+sqrt(5))/2 = {(3+sqrt(5))/2:.10f}")
print(f"  = phi + 1 = phi^2 = {phi**2:.10f}")
print()
print("  Q(1/sqrt(5)) = phi^2. So arctanh(1/sqrt(5)) = ln(phi^2)/2 = ln(phi).")
print(f"  Verify: ln(phi) = {log(phi):.10f}, arctanh(1/sqrt(5)) = {r_check:.10f}")
print(f"  Match: {abs(log(phi) - r_check) < 1e-10}")
print()
print("  So tanh(ln(phi)) = 1/sqrt(5).")
print("  And the Cayley address of phi^2 is 1/sqrt(5).")
print("  The pentagon IS in tanh: tanh(ln(phi)) = 1/sqrt(5) = cos(2*pi/5)/... ")
print()
# cos(2*pi/5) = (sqrt(5)-1)/4 = 0.309
# sin(2*pi/5) = sqrt(10+2*sqrt(5))/4 = 0.951
# 1/sqrt(5) = 0.4472
# 2*cos(2*pi/5) = (sqrt(5)-1)/2 = 1/phi = 0.618
print(f"  1/sqrt(5) = {1/sqrt(5):.6f}")
print(f"  sin(pi/5) = sin(36 deg) = {(sqrt(10-2*sqrt(5))/4):.6f}")
print(f"  Hmm, 1/sqrt(5) is not directly a pentagon trig value.")
print(f"  But sqrt(5) IS the diagonal/side of the unit square's face diagonal")
print(f"  in the golden rectangle. And phi = (1+sqrt(5))/2.")
print()
print("  The key identity: tanh(ln(phi)) = 1/sqrt(5).")
print("  Equivalently: the Cayley velocity of phi^2 = 1/sqrt(5).")
print("  phi, phi^2, and sqrt(5) are all tied together by tanh.")
print()

# ============================================================
print()
print("  VIII. WHAT tanh AND arctanh ARE")
print("  " + "-"*40)
print()
print("  tanh(x) answers: if I've taken x steps of rapidity,")
print("  what fraction of the journey to 1 have I completed?")
print()
print("  arctanh(v) answers: if I'm at position v on the journey,")
print("  how many steps of rapidity have I taken?")
print()
print("  They convert between EFFORT (rapidity) and PROGRESS (velocity).")
print("  Effort is unbounded. Progress is bounded.")
print("  The relationship is LOGARITHMIC: each unit of effort")
print("  closes a FIXED FRACTION of the remaining gap.")
print()
print("  This fraction is not 1/2 exactly. It's:")
frac = 1 - tanh(1)/1  # fraction of [0,1] closed per unit rapidity... not quite
# Actually: at rapidity x, gap = 1-tanh(x) = 2*e^{-2x}/(1+e^{-2x}) ~ 2*e^{-2x}
# From rapidity x to x+1: gap shrinks by factor e^{-2} = 0.135
# So each unit of rapidity closes (1 - e^{-2}) = 86.5% of the remaining gap
frac_close = 1 - exp(-2)
print(f"  Each unit of rapidity closes {frac_close:.1%} of the remaining gap.")
print(f"  (Because 1-tanh(x+1) = (1-tanh(x)) * e^{{-2}} / (1+small correction))")
print(f"  e^{{-2}} = {exp(-2):.4f}. So the remaining gap shrinks by factor 0.135 per step.")
print()
print("  86.5% closed per step. This is the 'good enough' per unit rapidity.")
print("  After 3 steps: closed 99.75%. After 5: 99.993%.")
print("  The approach to 1 is GEOMETRIC in rapidity, with ratio e^{{-2}} = 0.135.")
print()

# ============================================================
print()
print("  IX. THE FINAL UNDERSTANDING")
print("  " + "-"*40)
print()
print("  tanh is the function that says:")
print("  'No matter how much effort you put in, you can't reach 1.'")
print("  'But each unit of effort gets you 86.5% of what remains.'")
print("  'So practically, you're good enough after a few steps.'")
print()
print("  arctanh is the function that says:")
print("  'However close to 1 you are, I can tell you how many")
print("   steps it took to get there.'")
print("  'And the answer grows without bound as you approach 1.'")
print("  'So the journey to 1 is infinite, even though each step is finite.'")
print()
print("  Together they encode the FUNDAMENTAL ASYMMETRY of the real line:")
print("  The interval [0,1) is FINITE in length but INFINITE in rapidity.")
print("  The half-line [0, infinity) is INFINITE in length but")
print("  corresponds to exactly [0,1) in velocity.")
print()
print("  tanh and arctanh are not just functions.")
print("  They are the conversion between FINITE and INFINITE.")
print("  Between the bounded world (where we live)")
print("  and the unbounded world (where the structure lives).")
print()
print("  And the 'good enough' principle is built in:")
print("  You never need infinite rapidity.")
print("  A few octaves will do.")
print("  3 octaves gets you to 99.5%.")
print("  5 octaves gets you to 99.99%.")
print("  10 octaves gets you to 99.999999%.")
print()
print("  The rhythm: 3.32 octaves per decimal digit of precision.")
print("  The 3-3-4 pattern of powers of 2.")
print("  The pentagon angle of 108 degrees.")
print("  All the same thing: the geometry of 'close enough to 1")
print("  that the difference doesn't matter'.")
print()
print("  That is what tanh knows.")
print("  That is what arctanh counts.")
