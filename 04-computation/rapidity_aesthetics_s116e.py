#!/usr/bin/env python3
"""rapidity_aesthetics_s116e.py — Art, architecture, and aesthetics through rapidity.

The golden ratio lives at rapidity (3/2)*ln(phi).
The musical intervals ARE the alkane series ARE the dihedral groups.
What does BEAUTY look like in rapidity space?
"""
from math import sqrt, log, pi, exp, atanh

phi = (1+sqrt(5))/2

print("ART, ARCHITECTURE, AND AESTHETICS THROUGH RAPIDITY")
print("="*70)
print()

# ============================================================
print("="*70)
print("1. THE GOLDEN RATIO IN RAPIDITY")
print("="*70)
print()
print("  phi = (1+sqrt(5))/2 = {:.10f}".format(phi))
print("  rapidity(phi) = ln(phi)/2 = {:.10f}".format(log(phi)/2))
print("  = {:.4f} octaves".format(log(phi)/(2*log(2)/2)))
print()
print("  Cayley address of phi: (phi-1)/(phi+1) = 1/phi/(1+1/phi)")
addr = (phi-1)/(phi+1)
print(f"  = {addr:.10f}")
print(f"  Q({addr:.6f}) = {(1+addr)/(1-addr):.10f} = phi")
print()
print("  The golden rectangle: sides 1 and phi.")
print("  Aspect ratio = phi = 1.618...")
print("  rapidity(aspect) = ln(phi)/2 = 0.2406")
print("  In musical terms: rapidity 0.2406 = a MAJOR THIRD (5/4) minus a tiny comma.")
print(f"  rapidity(5/4) = {log(5/4)/2:.6f}")
print(f"  rapidity(phi) = {log(phi)/2:.6f}")
print(f"  difference = {log(phi)/2 - log(5/4)/2:.6f} = {(log(phi)/2 - log(5/4)/2)/(log(2)/2)*1200:.1f} cents")
print()
print("  The golden ratio is a major third plus 28 cents!")
print("  It sits BETWEEN the major third (5/4) and the minor third (6/5)")
print("  in rapidity space. The golden rectangle is musically ambiguous.")
print()

# ============================================================
print("="*70)
print("2. PAPER SIZES IN RAPIDITY")
print("="*70)
print()
print("  ISO A-series paper: aspect ratio sqrt(2):1.")
print("  rapidity(sqrt(2)) = ln(sqrt(2))/2 = ln(2)/4 = {:.6f}".format(log(2)/4))
print("  = half an octave = exactly one tritone!")
print()
print("  US Letter: 8.5 x 11, aspect = 11/8.5 = 1.2941")
print(f"  rapidity = {log(11/8.5)/2:.6f}")
print(f"  The US Letter is {log(11/8.5)/2 / (log(2)/4):.4f} times the A-series rapidity.")
print()
print("  Other paper ratios:")
papers = [
    ("A4 (ISO)", sqrt(2)),
    ("US Letter", 11/8.5),
    ("US Legal", 14/8.5),
    ("Golden rectangle", phi),
    ("Square", 1.0),
    ("16:9 screen", 16/9),
    ("4:3 screen", 4/3),
    ("2.39:1 cinema", 2.39),
    ("IMAX", 1.43),
]
print(f"  {'Format':<20s} {'Ratio':>8s} {'rapidity':>10s} {'octaves':>10s}")
print("  " + "-"*55)
for name, ratio in papers:
    if ratio > 1:
        rap = log(ratio)/2
    else:
        rap = 0
    oct_val = rap / (log(2)/2) if rap > 0 else 0
    print(f"  {name:<20s} {ratio:>8.4f} {rap:>10.6f} {oct_val:>10.4f}")

print()
print("  A4 paper = half an octave = one tritone.")
print("  4:3 screen = perfect fourth (Q(1/7) = 4/3).")
print("  16:9 screen = slightly less than a major second.")
print("  Cinema 2.39:1 = slightly more than one octave.")
print()

# ============================================================
print("="*70)
print("3. ARCHITECTURAL PROPORTIONS")
print("="*70)
print()
print("  Classical architecture uses specific proportional systems.")
print("  The Parthenon: width/height ~ phi (debated but close).")
print("  Gothic cathedrals: height/width ~ 2 (one octave).")
print()
print("  Key architectural ratios and their rapidity:")
arch_ratios = [
    ("1:1 (cube)", 1.0, "Platonic, symmetry"),
    ("1:sqrt(2)", sqrt(2), "A4 paper, ISO"),
    ("1:phi", phi, "Golden rectangle"),
    ("1:2 (double square)", 2.0, "Gothic, doubling"),
    ("1:3 (triple)", 3.0, "Classical columns"),
    ("3:4", 4/3, "Perfect fourth"),
    ("2:3", 3/2, "Perfect fifth"),
    ("5:8 ~ phi", 8/5, "Renaissance ideal"),
]
print(f"  {'Ratio':<22s} {'Value':>8s} {'rapidity':>10s} {'musical':>15s}")
print("  " + "-"*60)
for name, ratio, note in arch_ratios:
    rap = log(ratio)/2 if ratio > 1 else 0
    musical = ""
    for int_name, int_val in [("unison", 1), ("minor 3rd", 6/5), ("major 3rd", 5/4),
                                ("fourth", 4/3), ("tritone", sqrt(2)),
                                ("fifth", 3/2), ("octave", 2), ("twelfth", 3)]:
        if abs(ratio - int_val) < 0.02:
            musical = int_name
    print(f"  {name:<22s} {ratio:>8.4f} {rap:>10.6f} {musical:>15s}")

print()

# ============================================================
print("="*70)
print("4. COLOR THEORY AND RAPIDITY")
print("="*70)
print()
print("  Light wavelengths span 400-700 nm (0.81 octaves).")
print("  RGB color model: each channel 0-255 (8 octaves).")
print()
print("  The rapidity of a color channel value v (0-255):")
print("  rapidity(v) = ln(v)/2 (for v > 0)")
print()
print("  A color (R,G,B) has a rapidity VECTOR:")
print("  (rapidity(R), rapidity(G), rapidity(B))")
print()
print("  Key colors and their rapidity vectors:")
colors = [
    ("White", 255, 255, 255),
    ("Black", 1, 1, 1),
    ("Red", 255, 0, 0),
    ("Green", 0, 255, 0),
    ("Blue", 0, 0, 255),
    ("Yellow", 255, 255, 0),
    ("Cyan", 0, 255, 255),
    ("Gold", 255, 215, 0),
    ("Royal Blue", 65, 105, 225),
]
print(f"  {'Color':<14s} {'R':>4s} {'G':>4s} {'B':>4s}  {'rap(R)':>8s} {'rap(G)':>8s} {'rap(B)':>8s}")
print("  " + "-"*55)
for name, r, g, b in colors:
    rr = log(max(r,1))/2
    rg = log(max(g,1))/2
    rb = log(max(b,1))/2
    print(f"  {name:<14s} {r:4d} {g:4d} {b:4d}  {rr:8.4f} {rg:8.4f} {rb:8.4f}")

print()
print("  White = maximum rapidity on all channels = (2.77, 2.77, 2.77).")
print("  Black = minimum rapidity = (0, 0, 0).")
print("  The BRIGHTNESS of a color is the SUM of channel rapidities.")
print("  This is related to the luminance formula Y = 0.299R + 0.587G + 0.114B")
print("  but in LOG (rapidity) space rather than linear space.")
print()

# ============================================================
print("="*70)
print("5. MUSICAL TEMPERAMENT AND ARCHITECTURAL PROPORTION")
print("="*70)
print()
print("  The Pythagorean comma (12 fifths != 7 octaves) is")
print("  the SAME KIND of incommensurability as the golden ratio")
print("  being irrational (phi != p/q for any integers).")
print()
print("  In architecture: you can't tile a rectangle with golden ratio")
print("  using integer-sided squares (because phi is irrational).")
print("  This is the architectural COMMA.")
print()
print("  In music: you can't stack fifths to get exact octaves")
print("  (because ln(3)/ln(2) is irrational).")
print("  This is the musical COMMA.")
print()
print("  In both cases: the rapidity of the ratio is IRRATIONAL,")
print("  meaning it never lands on the integer rapidity lattice.")
print("  The beauty of both music and architecture comes from")
print("  APPROACHING but never reaching the integer lattice.")
print("  The near-miss IS the beauty.")
print()

# ============================================================
print("="*70)
print("6. PHOTOGRAPHY AND THE EXPOSURE TRIANGLE")
print("="*70)
print()
print("  Photography exposure: E = (f/N)^2 * t * ISO")
print("  where N = f-number, t = shutter speed, ISO = sensitivity.")
print()
print("  The STOP is the unit of exposure: 1 stop = doubling = 1 octave.")
print("  The f-stop series: f/1, f/1.4, f/2, f/2.8, f/4, f/5.6, f/8, ...")
print("  Each step is sqrt(2) = half an octave of f-number.")
print("  (Because light intensity ~ 1/N^2, so 1 stop = N*sqrt(2).)")
print()
print("  In rapidity:")
print("  1 stop of exposure = 1 octave of rapidity.")
print("  The f-stop number has rapidity that increases by ln(2)/4 per stop.")
print("  The ISO has rapidity that increases by ln(2)/2 per stop.")
print("  The shutter speed has rapidity that decreases by ln(2)/2 per stop.")
print()
print("  The exposure triangle IS a rapidity triangle:")
print("  rapidity(aperture) + rapidity(shutter) + rapidity(ISO) = const.")
print()

# ============================================================
print("="*70)
print("7. FIBONACCI SPIRALS AND THE GOLDEN RAPIDITY")
print("="*70)
print()
print("  The Fibonacci spiral: each quarter-turn expands by phi.")
print("  After n quarter-turns: expansion = phi^n.")
print("  rapidity(phi^n) = n * ln(phi)/2 = n * 0.2406.")
print()
print("  After a FULL TURN (4 quarter-turns):")
print("  expansion = phi^4 = phi^3 + phi^2 = 3*phi + 2 = {:.4f}".format(phi**4))
print("  rapidity = 4 * ln(phi)/2 = 2*ln(phi) = {:.6f}".format(2*log(phi)))
print("  = {:.4f} octaves".format(2*log(phi) / (log(2)/2)))
print()
print("  A full turn of the Fibonacci spiral expands by ~2.77 octaves.")
print("  Compare: an octave of music spans a frequency ratio of 2.")
print("  The Fibonacci spiral is MORE than a doubling per turn.")
print()

# ============================================================
print("="*70)
print("8. PERCEPTION AND RAPIDITY: STEVENS' POWER LAW")
print("="*70)
print()
print("  Stevens' power law: perceived intensity P = k * S^n")
print("  where S = stimulus, n depends on the modality.")
print()
print("  In rapidity: rapidity(P) = n * rapidity(S) + rapidity(k)")
print("  Stevens' law is LINEAR in rapidity space!")
print()
print("  Exponents for different modalities:")
modalities = [
    ("Brightness", 0.33),
    ("Loudness", 0.67),
    ("Vibration (finger)", 0.95),
    ("Length", 1.0),
    ("Weight", 1.45),
    ("Electric shock", 3.5),
    ("Temperature (warmth)", 1.6),
    ("Taste (sweetness)", 1.3),
]
print(f"  {'Modality':<25s} {'Exponent n':>10s} {'rapidity scale':>15s}")
print("  " + "-"*55)
for name, n in modalities:
    interp = ""
    if n < 0.5:
        interp = "compressive"
    elif n < 1.0:
        interp = "sub-linear"
    elif abs(n - 1.0) < 0.1:
        interp = "veridical"
    elif n < 2.0:
        interp = "expansive"
    else:
        interp = "highly expansive"
    print(f"  {name:<25s} {n:>10.2f} {interp:>15s}")

print()
print("  Brightness (n=0.33): 3 octaves of light = 1 octave of perceived brightness.")
print("  Loudness (n=0.67): 3 octaves of sound = 2 octaves of perceived loudness.")
print("  Electric shock (n=3.5): 1 octave of current = 3.5 octaves of perceived pain!")
print()
print("  Stevens' law in rapidity: P_rapidity = n * S_rapidity + const.")
print("  The exponent n IS the rapidity conversion factor between")
print("  physical stimulus and perceived intensity.")
print()

# ============================================================
print("="*70)
print("GRAND SUMMARY: BEAUTY IS RAPIDITY NEAR-COMMENSURABILITY")
print("="*70)
print()
print("  Beautiful proportions are ratios with SMALL rapidity.")
print("  Perfect symmetry = rapidity 0 (unison).")
print("  Golden ratio = rapidity 0.24 (between major and minor third).")
print("  Octave = rapidity 0.35 (the fundamental doubling).")
print()
print("  BEAUTY arises when rapidities are NEARLY but not exactly")
print("  commensurable. The golden ratio is beautiful because")
print("  it is maximally INCOMMENSURABLE (worst rational approximation).")
print("  It sits at the point of maximum distance from all simple fractions.")
print()
print("  Architecture uses rational approximations to phi (5:8, 3:5).")
print("  Music uses rational approximations to log-ratios (12-TET ~ just).")
print("  Both are TEMPERED versions of irrational rapidity ideals.")
print("  The tempering IS the beauty.")
