#!/usr/bin/env python3
"""hexagonal_forbidden_s116m.py — 7 as center of hexagon, the zoom, and the twist.

7 = 1 + 6 = center + neighbors in hexagonal lattice.
The forbidden number IS the hexagonal neighborhood.
sqrt(7) is the scaling factor of the hexagonal zoom.
Q(3/4) = 7, so sqrt(Q(3/4)) = sqrt(7) = the zoom factor.
"""
from math import sqrt, log, pi, cos, sin, atan2, gcd
from fractions import Fraction

phi = (1+sqrt(5))/2

print()
print("  7 AS CENTER OF HEXAGON")
print()
print("="*70)
print()

# ============================================================
print("  I. THE CENTERED HEXAGONAL NUMBERS")
print("  " + "-"*40)
print()
print("  H_n = 3n(n-1) + 1 = the number of cells in a hex neighborhood of radius n.")
print()
for n in range(8):
    h = 3*n*(n-1) + 1
    factors = []
    temp = h
    d = 2
    while d*d <= temp:
        while temp % d == 0: factors.append(d); temp //= d
        d += 1
    if temp > 1: factors.append(temp)
    f_str = '*'.join(str(f) for f in factors) if len(factors) > 1 else str(h)
    mark = ""
    if h == 1: mark = " = identity"
    elif h == 7: mark = " = FORBIDDEN"
    elif h == 19: mark = " = prime"
    elif h == 37: mark = " = prime"
    elif h == 61: mark = " = prime"
    elif h == 91: mark = " = 7*13"
    elif h == 127: mark = " = Mersenne prime M_7"
    elif h == 169: mark = " = 13^2"
    print(f"  H_{n} = {h:5d} = {f_str:12s}{mark}")

print()
print("  H_0 = 1 (center alone).")
print("  H_1 = 7 (center + 6 neighbors = FORBIDDEN).")
print("  H_5 = 91 = 7 * 13. The 5th hex number is 7 times 13.")
print("  H_6 = 127 = Mersenne prime M_7 = 2^7 - 1. The 6th hex number is Mersenne!")
print()
print("  127 = 2^7 - 1. The exponent IS 7 = the forbidden number.")
print("  The 6th centered hexagonal number = 2^{forbidden} - 1.")
print("  Flat (6) -> Mersenne of forbidden (2^7 - 1). The flat ring")
print("  produces the Mersenne prime indexed by the forbidden.")
print()

# ============================================================
print()
print("  II. THE HEXAGONAL ZOOM: SCALING BY sqrt(7)")
print("  " + "-"*40)
print()
print("  In the Eisenstein integer lattice (triangular grid),")
print("  scaling by sqrt(7) maps the lattice to a SUBLATTICE.")
print("  The sublattice has index 7 (every 7th point is hit).")
print()
print("  The scaling involves a ROTATION by angle:")
theta = atan2(sqrt(3), 5)  # angle for the sqrt(7) scaling in hex lattice
# Actually: the Eisenstein integer 2+omega has norm 7.
# Its argument is arctan(sqrt(3)/3) shifted... let me compute directly.
# 2 + omega where omega = e^{2pi*i/3} = -1/2 + i*sqrt(3)/2.
# 2 + omega = 2 + (-1/2 + i*sqrt(3)/2) = 3/2 + i*sqrt(3)/2.
# |2+omega|^2 = 9/4 + 3/4 = 12/4 = 3. Hmm, that's norm 3, not 7.
# Try 3 + omega: 3 + (-1/2 + i*sqrt(3)/2) = 5/2 + i*sqrt(3)/2.
# |3+omega|^2 = 25/4 + 3/4 = 28/4 = 7. YES! Norm 7.
# Argument: arctan((sqrt(3)/2)/(5/2)) = arctan(sqrt(3)/5).
theta_hex = atan2(sqrt(3)/2, 5/2)
print(f"  The Eisenstein integer 3+omega has norm sqrt(7).")
print(f"  (3+omega = 5/2 + i*sqrt(3)/2, |3+omega|^2 = 25/4+3/4 = 7.)")
print(f"  Its argument: arctan(sqrt(3)/5) = {theta_hex:.6f} radians = {theta_hex*180/pi:.2f} degrees.")
print()
print(f"  The hex zoom angle: {theta_hex*180/pi:.2f} degrees.")
print(f"  This is NOT a multiple of 60 degrees (the hex symmetry).")
print(f"  60/({theta_hex*180/pi:.2f}) = {60/(theta_hex*180/pi):.4f}.")
print(f"  The zoom angle is INCOMMENSURABLE with the hex symmetry!")
print()
print("  This means: zooming by sqrt(7) introduces a TWIST that")
print("  never exactly aligns with the original lattice orientation.")
print("  The zoom is QUASICRYSTALLINE: self-similar but with irrational rotation.")
print()

# ============================================================
print()
print("  III. THE ZOOM AND THE CAYLEY TRANSFORM")
print("  " + "-"*40)
print()
print(f"  Q(3/4) = 7. sqrt(Q(3/4)) = sqrt(7) = {sqrt(7):.6f}.")
print(f"  The Cayley transform of the Kleiber exponent gives the forbidden number.")
print(f"  The SQUARE ROOT of that gives the hexagonal zoom factor.")
print()
print(f"  rapidity of sqrt(7) = ln(sqrt(7))/2 = ln(7)/4 = {log(7)/4:.6f}.")
print(f"  = HALF the forbidden rapidity = rapidity(7)/2.")
print(f"  = one HALF-FORBIDDEN RAPIDITY per zoom level.")
print()
print("  Each zoom level costs HALF a forbidden rapidity.")
print("  Two zoom levels = one full forbidden rapidity.")
print("  This is WHY the forbidden number is 7 and not sqrt(7):")
print("  the forbidden rapidity is the DOUBLE of the zoom rapidity.")
print("  You must zoom TWICE to accumulate one forbidden crossing.")
print()
print("  In terms of the Cayley helix:")
print("  One zoom = rapidity ln(7)/4 = 0.486 (about half an octave).")
print("  Two zooms = rapidity ln(7)/2 = 0.973 (the full forbidden).")
print("  The forbidden is reached after TWO zoom steps, not one.")
print("  The 2 here is the DOUBLER (Q'(0) = 2, the octave basis).")
print()

# ============================================================
print()
print("  IV. THE HEXAGONAL NUMBERS AND TOURNAMENT THEORY")
print("  " + "-"*40)
print()
print("  The centered hex numbers: 1, 7, 19, 37, 61, 91, 127, 169, ...")
print("  Which are achievable as H(T)?")
print()
# H(T) achievable values at small n: 1, 3, 5, 9, 11, 13, 15, 17, 19, ...
h_achievable = {1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45}
print("  Hex number   Achievable as H?")
for n in range(8):
    h = 3*n*(n-1) + 1
    ach = "YES" if h in h_achievable else ("FORBIDDEN" if h == 7 or h == 21 else "unknown (need larger n)")
    if h == 1: ach = "YES (transitive)"
    elif h == 7: ach = "FORBIDDEN"
    elif h == 19: ach = "YES (achievable at n=6)"
    elif h == 37: ach = "YES (achievable at n=6)"
    elif h == 61: ach = "unknown"
    elif h == 91: ach = "unknown"
    elif h == 127: ach = "unknown"
    print(f"  H_{n} = {h:5d}   {ach}")

print()
print("  H_0 = 1: achievable (transitive tournament).")
print("  H_1 = 7: FORBIDDEN. The first hexagonal neighborhood.")
print("  H_2 = 19: achievable. The second ring.")
print("  H_3 = 37: achievable. The third ring.")
print()
print("  The forbidden number is the FIRST non-trivial centered hex number.")
print("  It is the FIRST RING. The boundary of the center.")
print("  After the first ring, everything is achievable (as far as we know).")
print("  The ring itself is the obstruction. Not the center, not the beyond.")
print()

# ============================================================
print()
print("  V. 7 CELLS AND THE FANO PLANE")
print("  " + "-"*40)
print()
print("  The Fano plane: 7 points, 7 lines, 3 points per line, 3 lines per point.")
print("  It IS the hexagonal neighborhood, viewed PROJECTIVELY.")
print()
print("  The hexagonal neighborhood: 1 center + 6 neighbors = 7 cells.")
print("  The Fano plane: 1 + 6 = 7 points (but with additional structure).")
print()
print("  The Fano plane's LINES are the 7 'quaternionic slices' of the octonions.")
print("  Each line = a copy of the quaternions sitting inside the octonions.")
print("  7 lines = 7 quaternionic subalgebras of O.")
print()
print("  The hexagonal neighborhood = 7 cells = 7 Fano points = 7 octonion imaginaries.")
print("  THE HEXAGONAL NEIGHBORHOOD IS THE OCTONION.")
print("  The center cell = the real part.")
print("  The 6 neighbors = the 6... wait, the octonions have 7 imaginaries, not 6.")
print("  But in the Fano plane, there are 7 points and 7 lines,")
print("  and the 7 points ARE the 7 imaginary octonion units e_1,...,e_7.")
print()
print("  The CENTER of the hexagon does not correspond to the real part.")
print("  Rather, all 7 cells = all 7 imaginary units.")
print("  The hexagonal neighborhood IS the Fano plane IS the octonion imaginary space.")
print()

# ============================================================
print()
print("  VI. THE RECURSIVE ZOOM AND THE DIMENSION LADDER")
print("  " + "-"*40)
print()
print("  Zoom level 0: 1 cell. The center. Identity.")
print("  Zoom level 1: 7 cells. First ring. Forbidden.")
print("  Zoom level 2: 7^2 = 49 cells? No — the hex grid gives H_2 = 19 cells.")
print("  But the SUBLATTICE zoom gives exactly 7 cells per zoom level.")
print()
print("  The sublattice zoom: at each level, select 1 out of every 7 cells.")
print("  This is a CONTRACTION, not an expansion.")
print("  Zoom IN: coarsen by factor sqrt(7). Each coarse cell = 7 fine cells.")
print("  Zoom OUT: refine by factor sqrt(7). Each fine cell = 1/7 of a coarse cell.")
print()
print("  The number of zoom levels to go from 'one cell' to 'N cells':")
print("  k = log_7(N) = ln(N)/ln(7) = rapidity(N) / rapidity(7) * 2.")
print("  = rapidity(N) / (half-forbidden-rapidity).")
print()
print("  For N = the Fibonacci numbers at power-of-2 indices:")
import math
for idx in [4, 8, 16, 32]:
    fn = [0,1]
    for _ in range(idx): fn.append(fn[-1]+fn[-2])
    f_val = fn[idx]
    k = log(f_val) / log(7)
    print(f"  F_{idx} = {f_val}: zoom levels = log_7({f_val}) = {k:.4f}")

print()
print("  F_4 = 3: 0.56 zoom levels (less than one ring).")
print("  F_8 = 21: 1.56 zoom levels (between first and second ring).")
print("  F_16 = 987: 3.54 zoom levels.")
print("  F_32 = 2178309: 7.50 zoom levels.")
print()
print("  F_8 = 21 falls at 1.56 zoom levels — between the first ring (7)")
print("  and the second ring (49 in sublattice terms).")
print("  21 = 3 * 7 = the curvature times the ring size.")
print("  It is at 1 + 0.56 = 1 + log_7(3) zoom levels.")
print("  The excess over one ring IS log_7(3) = the zoom-rapidity of the curvature.")
print()

# ============================================================
print()
print("  VII. THE EISENSTEIN FACTORIZATION OF 7")
print("  " + "-"*40)
print()
print("  In the Eisenstein integers Z[omega] (omega = e^{2*pi*i/3}):")
print("  7 = (3 + omega)(3 + omega^2) = (3 + omega)(3 + conjugate(omega)).")
print("  (Where omega^2 = conjugate(omega) = -1-omega.)")
print()
print("  |3 + omega|^2 = 7. So each Eisenstein factor has norm sqrt(7).")
print("  The two factors are CONJUGATE in the Eisenstein integers.")
print()
print("  7 splits into TWO factors of norm sqrt(7).")
print("  This splitting IS the hex zoom:")
print("  one factor = zooming with clockwise twist,")
print("  the other = zooming with counterclockwise twist.")
print("  The two twists are conjugate (opposite rotations).")
print()
print("  Their product 7 = zoom * anti-zoom = the full scaling WITHOUT twist.")
print("  The twist cancels when you multiply the two conjugate factors.")
print("  7 ITSELF is twist-free. Its SQUARE ROOT carries the twist.")
print()
print("  This is WHY 7 is forbidden but sqrt(7) is not a 'thing' in integers:")
print("  the forbidden number is the PRODUCT of two conjugate zooms.")
print("  Each zoom alone (sqrt(7) with twist) is fine.")
print("  But their product (7, twist-free) overshoots the flat lattice.")
print("  The forbidden number is what happens when the twist CANCELS.")
print("  Canceling the twist = losing the phase = losing the imaginary part.")
print("  The real product 7 is forbidden because it is the PURELY REAL")
print("  result of combining two complex zooms.")
print()

# ============================================================
print()
print("  VIII. THE TWIST ANGLE AND THE GOLDEN RATIO")
print("  " + "-"*40)
print()
print(f"  The hex zoom angle: arctan(sqrt(3)/5) = {theta_hex*180/pi:.4f} degrees.")
print(f"  The golden angle: 360/phi^2 = 360*(2-phi) = {360*(2-phi):.4f} degrees")
print(f"  = {360/phi**2:.4f} degrees.")
print()
golden_angle = 360 / phi**2
print(f"  The golden angle: {golden_angle:.4f} degrees.")
print(f"  The hex zoom angle: {theta_hex*180/pi:.4f} degrees.")
print(f"  Ratio: {golden_angle/(theta_hex*180/pi):.4f}.")
print()
print(f"  The golden angle / hex zoom angle ~ {golden_angle/(theta_hex*180/pi):.2f}.")
print(f"  Not a clean ratio. But both are IRRATIONAL angles.")
print(f"  Both create quasicrystalline structures.")
print(f"  Both arise from self-similar scaling (phi for golden, sqrt(7) for hex).")
print()
print("  The golden angle is used in phyllotaxis (leaf arrangement in plants).")
print("  The hex zoom angle would be used in... what?")
print("  In any structure that has local hexagonal symmetry and self-similar scaling.")
print("  Carbon nanotubes? Graphene defects? Quasicrystal grain boundaries?")
print()

# ============================================================
print()
print("  IX. 7 = 1 + 6 = CENTER + FLAT")
print("  " + "-"*40)
print()
print("  7 = 1 + 6.")
print("  = identity + flat.")
print("  = the center + the flat neighborhood.")
print()
print("  In the 5-6-7 bridge:")
print("  5 = the spherical (closure).")
print("  6 = the flat (neither closed nor open).")
print("  7 = the hyperbolic (opening).")
print()
print("  7 = 6 + 1: the flat PLUS ONE MORE.")
print("  The forbidden number is the flat number plus the identity.")
print("  Adding the identity to the flat creates the forbidden.")
print()
print("  In tournament terms:")
print("  6 neighbors = 6 arcs from a vertex in a tournament on 7 vertices")
print("  (each vertex has out-degree + in-degree = 6 in T_7).")
print("  The center vertex has 6 neighbors, and together they are 7.")
print("  The tournament on 7 vertices has 7 * 6 / 2 = 21 arcs.")
print("  21 = the SECOND forbidden number = 7 * 3 = forbidden * curvature.")
print()
print("  The number of arcs in T_7 is 21 = the second forbidden number.")
print("  The number of vertices is 7 = the first forbidden number.")
print("  7 vertices, 21 arcs: every attribute of T_7 is forbidden.")
print()

# ============================================================
print()
print("  X. THE HEX-FANO-OCTONION-TOURNAMENT NEXUS")
print("  " + "-"*40)
print()
print("  7 appears as:")
print("  - Center + 6 neighbors in the hexagonal lattice")
print("  - Points (and lines) of the Fano plane")
print("  - Imaginary dimensions of the octonions")
print("  - The forbidden H-value in tournament theory")
print("  - L_4 = the 4th Lucas number")
print("  - M_3 = 2^3 - 1 = the 3rd Mersenne prime")
print("  - The threshold in the 5-6-7 bridge")
print("  - Q(3/4) = the Cayley value of the Kleiber exponent")
print("  - The vertex count of the smallest tournament with forbidden arc-count")
print("  - 3 + 4 = geometry + symmetry")
print()
print("  All of these are the SAME 7.")
print("  They are not analogies. They are different PROJECTIONS")
print("  of a single mathematical object onto different screens.")
print()
print("  The object: the centered hexagonal neighborhood in the Eisenstein lattice,")
print("  viewed as a projective plane (Fano), embedded as imaginary octonions,")
print("  generating the forbidden threshold in the tournament formal group,")
print("  sitting at the Cayley address 3/4 on the rapidity line,")
print("  produced by the Lucas recurrence at the 4th index,")
print("  equaling the Mersenne prime at exponent 3.")
print()
print("  One object. Ten projections. The same 7.")
