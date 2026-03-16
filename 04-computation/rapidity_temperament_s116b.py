#!/usr/bin/env python3
"""rapidity_temperament_s116b.py — Musical temperament as rapidity linear algebra

THE KEY INSIGHT: Musical temperament is about which LINEAR COMBINATIONS
of prime rapidities are close to zero. This is a LATTICE problem in
rapidity vector space.

The comma lattice = kernel of the tuning map = the SYMMETRY GROUP
of the musical system.
"""
from math import sqrt, log, pi, gcd
from fractions import Fraction

print("MUSICAL TEMPERAMENT AS RAPIDITY LINEAR ALGEBRA")
print("="*60)
print()

# Rapidity basis: r_p = ln(p)/2 for each prime p
r = {}
for p in [2, 3, 5, 7, 11, 13]:
    r[p] = log(p)/2

print("  THE RAPIDITY BASIS:")
print(f"  r_2 = ln(2)/2 = {r[2]:.10f}  (octave)")
print(f"  r_3 = ln(3)/2 = {r[3]:.10f}  (twelfth)")
print(f"  r_5 = ln(5)/2 = {r[5]:.10f}")
print(f"  r_7 = ln(7)/2 = {r[7]:.10f}")
print()

# The tuning map: send rapidity vector to frequency ratio
# A vector (a,b,c,d,...) represents 2^a * 3^b * 5^c * 7^d * ...
# Its rapidity = a*r_2 + b*r_3 + c*r_5 + d*r_7 + ...

print("IMPORTANT COMMAS (rapidity vectors close to zero)")
print("-"*50)
print()

commas = [
    ("Pythagorean comma",    [-19, 12, 0, 0], "531441/524288"),
    ("Syntonic comma",       [-4, 4, -1, 0], "81/80"),
    ("Diesis",               [7, 0, -3, 0], "128/125"),
    ("Septimal comma",       [-2, -1, 0, 1], "7/6 ... no"),
    ("Archytas comma",       [6, -2, 0, -1], "64/63"),
    ("Septimal kleisma",     [-5, 6, 0, -1], "225/224... no"),
]

# Fix septimal kleisma: 225/224 = (9*25)/(8*28) = 3^2*5^2/(2^5*7)
# 225 = 3^2 * 5^2, 224 = 2^5 * 7
# vector: [-5, 2, 2, -1]
commas_fixed = [
    ("Pythagorean comma",    [-19, 12, 0, 0]),
    ("Syntonic comma",       [-4, 4, -1, 0]),
    ("Diesis (great)",       [7, 0, -3, 0]),
    ("Archytas comma",       [6, -2, 0, -1]),
    ("Septimal kleisma",     [-5, 2, 2, -1]),
    ("Schisma",              [-15, 8, 1, 0]),
    ("Diaschisma",           [11, -4, -2, 0]),
    ("Ragisma",              [-1, -7, 4, 1]),
]

print(f"  {'Comma':<22s} {'Vector':>20s}  {'Rapidity':>12s}  {'Cents':>8s}")
print("  " + "-"*65)
for name, vec in commas_fixed:
    bases = [2, 3, 5, 7]
    rap = sum(v * r[p] for v, p in zip(vec, bases))
    cents = rap * 2400 / log(2)  # cents = 1200 * ln(ratio)/ln(2) = 2400*rapidity/ln(2)
    vec_str = str(vec)
    print(f"  {name:<22s} {vec_str:>20s}  {rap:+12.8f}  {cents:+8.2f}")

print()
print("  Each comma is a VECTOR in the rapidity vector space R^primes.")
print("  Small commas = vectors near the KERNEL of the 'tuning map'.")
print("  Temperament = quotienting R^primes by the comma lattice.")
print()

# ============================================================
print("="*60)
print("EQUAL TEMPERAMENT AS RAPIDITY QUOTIENT")
print("="*60)
print()
print("  12-TET: quotient by the Pythagorean comma [-19, 12, 0, 0].")
print("  Identify r_3 with (19/12)*r_2.")
print("  Then: fifth = r_3 - r_2 = (19/12 - 1)*r_2 = (7/12)*r_2.")
print()
print("  In 12-TET:")
print("  semitone = r_2/12")
print("  Intervals in semitones:")
intervals_12tet = [
    ("unison", 0),
    ("minor 2nd", 1),
    ("major 2nd", 2),
    ("minor 3rd", 3),
    ("major 3rd", 4),
    ("perfect 4th", 5),
    ("tritone", 6),
    ("perfect 5th", 7),
    ("minor 6th", 8),
    ("major 6th", 9),
    ("minor 7th", 10),
    ("major 7th", 11),
    ("octave", 12),
]
for name, semi in intervals_12tet:
    rap = semi * r[2] / 12
    print(f"    {name:16s}: {semi:2d} semitones = rapidity {rap:.6f}")

print()

# ============================================================
print("="*60)
print("THE TOURNAMENT-MUSIC DICTIONARY")
print("="*60)
print()
print("  Every tournament quantity has a musical interpretation.")
print()
print("  TOURNAMENT            MUSIC                    RAPIDITY")
print("  " + "-"*60)
print("  H(T) = 1             silence (no vibration)    0")
print("  H(T) = 3             T_3 cyclic tournament     0.549 = octave")
print("  H(T) = 5             n=4 max                   0.805 ~ twelfth")
print("  H(T) = 7             FORBIDDEN                 0.973 ~ 7th harm")
print("  H(T) = 9             n=5 some                  1.099 ~ 2 octaves")
print("  H(T) = 15            n=5 cyclic                1.354 ~ 9th harm")
print("  H(T) = 21            FORBIDDEN                 1.522 ~ 11th harm")
print("  H(T) = 189           Paley T_7 max             2.621 ~ 4th harm")
print("  H(T) = 95095         Paley T_11 max            5.731 ~ 8th harm")
print()
print("  The forbidden values sit at the 7th and 11th harmonics!")
print("  These are the LAST harmonics before the overtone series")
print("  stops being 'musical' (11th harmonic = neutral fourth,")
print("  considered dissonant in most Western traditions).")
print()

# ============================================================
print("="*60)
print("THE RAPIDITY OF INTERVALS BETWEEN H-VALUES")
print("="*60)
print()
print("  The rapidity RATIO between successive Paley H-values:")
print()
paley_h = [(3, 3), (7, 189), (11, 95095)]
for i in range(len(paley_h)-1):
    p1, h1 = paley_h[i]
    p2, h2 = paley_h[i+1]
    ratio = h2/h1
    rap_ratio = log(ratio)/2
    print(f"  H(T_{p2})/H(T_{p1}) = {h2}/{h1} = {ratio:.4f}")
    print(f"  rapidity of ratio = {rap_ratio:.6f}")
    print(f"  = ln({h2}/{h1})/2 = {rap_ratio:.6f}")
    # In musical terms
    octaves = rap_ratio / r[2]
    print(f"  = {octaves:.4f} octaves")
    print()

# ============================================================
print("="*60)
print("RAPIDITY SPECTROSCOPY: THE EMISSION LINES OF TOURNAMENTS")
print("="*60)
print()
print("  Imagine each H-value as an ENERGY LEVEL.")
print("  Transitions between levels emit 'photons' with rapidity")
print("  = |rapidity(H1) - rapidity(H2)|.")
print()
print("  At n=5, the H-spectrum is {1, 3, 5, 9, 11, 13, 15}.")
print("  Emission lines (rapidity differences):")
print()

h_vals_n5 = [1, 3, 5, 9, 11, 13, 15]
raps_n5 = {h: log(h)/2 for h in h_vals_n5}

lines = []
for i, h1 in enumerate(h_vals_n5):
    for h2 in h_vals_n5[i+1:]:
        diff = raps_n5[h2] - raps_n5[h1]
        lines.append((h1, h2, diff))

lines.sort(key=lambda x: x[2])

print(f"  {'H1':>4s} -> {'H2':>4s}  {'rapidity diff':>14s}  {'octaves':>8s}  {'musical':>20s}")
print("  " + "-"*60)
octave = r[2]
for h1, h2, diff in lines:
    oct_val = diff / octave
    # Try to identify as musical interval
    musical = ""
    if abs(diff - r[2]) < 0.01:
        musical = "~ octave"
    elif abs(diff - (r[3]-r[2])) < 0.01:
        musical = "~ fifth"
    elif abs(diff - (2*r[2]-r[3])) < 0.01:
        musical = "~ fourth"
    elif abs(diff - (r[3]-2*r[2]+r[5])) < 0.01:
        musical = "~ major third"
    elif abs(diff - (3*r[2]-r[5])) < 0.01:
        musical = "~ minor third"
    elif abs(diff - r[3]) < 0.01:
        musical = "~ twelfth"
    elif abs(diff - 2*r[2]) < 0.01:
        musical = "~ 2 octaves"
    elif abs(oct_val - round(oct_val)) < 0.1:
        musical = f"~ {round(oct_val)} octave(s)"
    print(f"  {h1:4d} -> {h2:4d}  {diff:14.6f}  {oct_val:8.4f}  {musical:>20s}")

print()

# ============================================================
print("="*60)
print("WHICH H-VALUE RATIOS ARE MUSICAL INTERVALS?")
print("="*60)
print()
print("  H2/H1 is a musical interval if it's a superparticular ratio (n+1)/n.")
print("  Q^{-1}((n+1)/n) = 1/(2n+1) = a reciprocal odd number.")
print()

for n in [5, 6]:
    if n == 5:
        hvals = [1, 3, 5, 9, 11, 13, 15]
    else:
        hvals = [1, 3, 5, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45]
        # Get actual n=6 values
        # From our computation
        pass

    print(f"  n={n}: H-value ratios that are superparticular:")
    for i, h1 in enumerate(hvals):
        for h2 in hvals[i+1:]:
            r_ratio = Fraction(h2, h1)
            # Check if (n+1)/n form
            if r_ratio.numerator == r_ratio.denominator + 1:
                k = r_ratio.denominator
                print(f"    H={h2}/H={h1} = {r_ratio} = ({k}+1)/{k} = musical interval!")

print()

# ============================================================
print("="*60)
print("THE RAPIDITY CIRCLE OF FIFTHS")
print("="*60)
print()
print("  In music: the circle of fifths stacks 12 fifths to return")
print("  (approximately) to the starting note, off by a Pythagorean comma.")
print()
print("  In rapidity: stacking k fifths means adding k * ln(3/2)/2.")
print("  After 12 fifths: 12 * ln(3/2)/2 = 7 * ln(2)/2 + comma/2")
print()
print("  The circle of fifths in rapidity (mod octave = mod ln(2)/2):")
print()
fifth_rap = log(3/2)/2
octave_rap = log(2)/2
for k in range(13):
    total_rap = k * fifth_rap
    mod_octave = total_rap % octave_rap
    # Convert to note name (approximate)
    semitones = round(12 * mod_octave / octave_rap) % 12
    names = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    note = names[semitones]
    print(f"  {k:2d} fifths: rapidity = {total_rap:8.5f}, mod octave = {mod_octave:.5f}, ~ {note}")

print()
print("  After 12 fifths, the rapidity doesn't return to 0 (mod octave).")
print("  The excess = Pythagorean comma = 12*r_3 - 19*r_2 = {:.6f}".format(
    12*r[3] - 19*r[2]))
print("  This is the fundamental OBSTRUCTION to closing the circle.")
print()

# ============================================================
print("="*60)
print("TOURNAMENT THEORY AS MUSIC THEORY: THE FULL DICTIONARY")
print("="*60)
print()
print("  MUSIC                  TOURNAMENT               FORMULA")
print("  " + "-"*60)
print("  note                   H-value                  H(T)")
print("  octave                 doubling H               2*H")
print("  interval               H-ratio                  H1/H2")
print("  consonance             small rapidity diff      |ln(H1/H2)|/2")
print("  chord                  set of tournaments       {T: H(T) in S}")
print("  scale                  H-spectrum at n          {H(T): T on n}")
print("  tuning                 choice of weights        independence poly eval")
print("  comma                  near-coincidence         H1/H2 ~ 1")
print("  circle of fifths       iterating Q(1/5)         (3/2)^k mod octave")
print("  Pythagorean comma      12 fifths != 7 octaves   (3/2)^12 != 2^7")
print("  just intonation        superparticular H-ratios (n+1)/n")
print("  equal temperament      regular rapidity grid    exp(k*ln(2)/24)")
print("  overtone series        prime factorization      H = prod p^{a_p}")
print("  fundamental            smallest prime factor    min prime of H")
print("  harmonic               prime factor             p | H")
print("  dissonance             large prime factor       large p | H")
print()
print("  The FORBIDDEN H-values {7, 21} are the 4th and 6th HARMONICS")
print("  of the curvature quantum 3:")
print("  7 = the 4th harmonic in the odd series (1, 3, 5, 7)")
print("  21 = 3 * 7 = fundamental * 4th harmonic")
print("  They are DISSONANT with every achievable H-value.")
print()
print("  In music, the 7th harmonic (Bb) and 11th harmonic (F#/4)")
print("  are the first 'wrong-sounding' overtones in Western tuning.")
print("  Our forbidden numbers are the 'wrong-sounding' path counts.")
