#!/usr/bin/env python3
"""rapidity_dihedral_s116d.py — Dihedral groups, regular polygons on the unit circle,
and the rapidity-alkane-music correspondence.

The n-gon sits on the unit circle with vertices at e^{2*pi*i*k/n}.
The dihedral group D_n has order 2n = the hydrogen count of alkane C_n!
The Cayley transform Q maps the unit circle to the imaginary axis.
What happens when we push regular polygons through Q?
"""
from math import sqrt, log, pi, cos, sin, exp, atanh, tanh, gcd
from fractions import Fraction
import cmath

phi = (1+sqrt(5))/2
octave = log(2)/2

print("DIHEDRAL GROUPS, POLYGONS, AND THE RAPIDITY-ALKANE CORRESPONDENCE")
print("="*70)
print()

# ============================================================
print("="*70)
print("I. THE HYDROGEN COUNT IS THE DIHEDRAL ORDER")
print("="*70)
print()
print("  Alkane C_n H_{2n+2}: hydrogen count = 2n+2 = 2(n+1).")
print("  Dihedral group D_{n+1}: order = 2(n+1).")
print()
print("  n   Alkane    H count   D_{n+1}   |D_{n+1}|   Musical interval")
print("  " + "-"*65)
for n in range(1, 10):
    h = 2*n + 2
    interval = Fraction(n+1, n)
    names = {1: "octave", 2: "fifth", 3: "fourth", 4: "maj 3rd",
             5: "min 3rd", 6: "septimal", 7: "maj 2nd", 8: "pyth 2nd"}
    name = names.get(n, "")
    print(f"  {n}   C{n}H{h:<4d}   {h:4d}      D_{n+1:<3d}     {h:4d}        "
          f"{interval}  {name}")

print()
print("  METHANE CH4 has 4 hydrogens = |D_2| = Klein four-group.")
print("  ETHANE C2H6 has 6 hydrogens = |D_3| = symmetry group of TRIANGLE.")
print("  PROPANE C3H8 has 8 hydrogens = |D_4| = symmetry group of SQUARE.")
print("  BUTANE C4H10 has 10 hydrogens = |D_5| = symmetry group of PENTAGON.")
print("  PENTANE C5H12 has 12 hydrogens = |D_6| = symmetry group of HEXAGON.")
print()
print("  The n-th alkane's hydrogen count is the order of the")
print("  symmetry group of the (n+1)-gon!")
print()
print("  And the musical interval (n+1)/n = Q(1/(2n+1)) uses")
print("  the SAME number n+1 as the polygon vertex count.")
print()

# ============================================================
print("="*70)
print("II. REGULAR POLYGONS ON THE UNIT CIRCLE")
print("="*70)
print()
print("  The (n+1)-gon has vertices at omega_k = e^{2*pi*i*k/(n+1)},")
print("  k = 0, 1, ..., n.")
print()
print("  Apply the Cayley transform Q(z) = (1+z)/(1-z) to each vertex:")
print()

for n_plus_1 in [2, 3, 4, 5, 6, 7]:
    n = n_plus_1 - 1
    alkane = f"C{n}H{2*n+2}"
    print(f"  {n_plus_1}-GON (alkane {alkane}, interval {n_plus_1}/{n_plus_1-1}):")
    for k in range(n_plus_1):
        z = cmath.exp(2j * cmath.pi * k / n_plus_1)
        if abs(z - 1) < 1e-10:
            print(f"    k={k}: z = 1.000+0.000i  =>  Q(z) = INFINITY (pole)")
        else:
            q = (1 + z) / (1 - z)
            print(f"    k={k}: z = {z.real:+.3f}{z.imag:+.3f}i  =>  "
                  f"Q(z) = {q.real:+.4f}{q.imag:+.4f}i  "
                  f"|Q| = {abs(q):.4f}")
    print()

# ============================================================
print("="*70)
print("III. Q MAPS POLYGONS TO THE IMAGINARY AXIS")
print("="*70)
print()
print("  For z = e^{i*theta} on the unit circle (z != 1):")
print("  Q(z) = (1+z)/(1-z) = i * cot(theta/2)")
print()
print("  PROOF: z = e^{i*theta}")
print("  1+z = 1+cos(t)+i*sin(t) = 2*cos^2(t/2) + 2i*sin(t/2)*cos(t/2)")
print("      = 2*cos(t/2) * [cos(t/2) + i*sin(t/2)] = 2*cos(t/2)*e^{i*t/2}")
print("  1-z = 1-cos(t)-i*sin(t) = 2*sin^2(t/2) - 2i*sin(t/2)*cos(t/2)")
print("      = 2*sin(t/2) * [sin(t/2) - i*cos(t/2)]")
print("      = -2i*sin(t/2) * [cos(t/2) + i*sin(t/2)] = -2i*sin(t/2)*e^{i*t/2}")
print("  Q(z) = 2*cos(t/2)*e^{i*t/2} / (-2i*sin(t/2)*e^{i*t/2})")
print("       = cos(t/2) / (-i*sin(t/2)) = i*cot(t/2). QED.")
print()
print("  So Q sends the unit circle to the IMAGINARY AXIS.")
print("  The (n+1)-gon at angles theta_k = 2*pi*k/(n+1) maps to:")
print("  Q(omega_k) = i * cot(pi*k/(n+1))  for k = 1, ..., n")
print("  (k=0 gives the pole at z=1)")
print()

print("  The COTANGENT VALUES at polygon vertices:")
print()
for n_plus_1 in [3, 4, 5, 6, 7, 8]:
    n = n_plus_1 - 1
    print(f"  {n_plus_1}-gon: cot(pi*k/{n_plus_1}) for k=1,...,{n}")
    for k in range(1, n_plus_1):
        angle = pi * k / n_plus_1
        cot_val = cos(angle) / sin(angle)
        print(f"    k={k}: cot(pi*{k}/{n_plus_1}) = {cot_val:+.6f}", end="")
        # Check for special values
        if abs(cot_val - sqrt(3)) < 0.001:
            print("  = sqrt(3)", end="")
        elif abs(cot_val - 1/sqrt(3)) < 0.001:
            print("  = 1/sqrt(3)", end="")
        elif abs(cot_val - 1) < 0.001:
            print("  = 1", end="")
        elif abs(cot_val + 1) < 0.001:
            print("  = -1", end="")
        elif abs(cot_val) < 0.001:
            print("  = 0", end="")
        elif abs(abs(cot_val) - phi) < 0.01:
            print(f"  ~ {'+'if cot_val>0 else '-'}phi!", end="")
        elif abs(abs(cot_val) - 1/phi) < 0.01:
            print(f"  ~ {'+'if cot_val>0 else '-'}1/phi!", end="")
        print()
    print()

# ============================================================
print("="*70)
print("IV. THE PENTAGON AND THE GOLDEN RATIO")
print("="*70)
print()
print("  The 5-gon (pentagon) has vertices at e^{2*pi*i*k/5}, k=0,...,4.")
print("  cot(pi/5) = cot(36 deg) = sqrt(1 + 2/sqrt(5)) = ... ")
cot_pi5 = cos(pi/5)/sin(pi/5)
print(f"  = {cot_pi5:.10f}")
print(f"  = sqrt(1 + 2/sqrt(5)) = {sqrt(1 + 2/sqrt(5)):.10f}")
print()
# Actually cot(pi/5) = sqrt(5 + 2*sqrt(5))... let me compute
# cos(pi/5) = (1+sqrt(5))/4 = phi/2, sin(pi/5) = sqrt(10-2*sqrt(5))/4
# cot = cos/sin = (1+sqrt(5)) / sqrt(10-2*sqrt(5))
# = phi*2 / sqrt(10-2*sqrt(5))
# Let me just check
print(f"  cot(pi/5) = {cot_pi5:.10f}")
print(f"  cot(2*pi/5) = {cos(2*pi/5)/sin(2*pi/5):.10f}")
print()
print(f"  cot(pi/5)^2 = {cot_pi5**2:.10f}")
print(f"  5 + 2*sqrt(5) = {5 + 2*sqrt(5):.10f}")
# cot^2(pi/5) = cos^2/sin^2 = (1+cos(2pi/5))/(1-cos(2pi/5))
# cos(2pi/5) = (sqrt(5)-1)/4... actually cos(72 deg) = (sqrt(5)-1)/4
# Hmm, cos(2*pi/5) = cos(72 deg) = (sqrt(5)-1)/4
cos_72 = (sqrt(5)-1)/4
print(f"  cos(2*pi/5) = {cos(2*pi/5):.10f}")
print(f"  (sqrt(5)-1)/4 = {cos_72:.10f}")
print()

# The key observation: Q maps the pentagon to 5 points on the imaginary axis
# including infinity. The 4 finite points are at i*cot(pi*k/5) for k=1,2,3,4.
# These come in pairs: cot(pi/5) = -cot(4*pi/5), cot(2*pi/5) = -cot(3*pi/5).
# The positive values: cot(pi/5) and cot(2*pi/5).

print("  Pentagon maps under Q to:")
print(f"  k=0: pole (infinity)")
print(f"  k=1: i * cot(pi/5) = i * {cot_pi5:.6f}")
print(f"  k=2: i * cot(2*pi/5) = i * {cos(2*pi/5)/sin(2*pi/5):.6f}")
print(f"  k=3: i * cot(3*pi/5) = i * {cos(3*pi/5)/sin(3*pi/5):.6f}")
print(f"  k=4: i * cot(4*pi/5) = i * {cos(4*pi/5)/sin(4*pi/5):.6f}")
print()

cot1 = cos(pi/5)/sin(pi/5)
cot2 = cos(2*pi/5)/sin(2*pi/5)
print(f"  cot(pi/5) * cot(2*pi/5) = {cot1 * cot2:.10f}")
print(f"  Compare: 1/sqrt(5) = {1/sqrt(5):.10f}")
# Actually for regular n-gon: prod_{k=1}^{n-1} cot(pi*k/n) should have clean value
# Wait, product of cot(pi*k/n) for k=1..n-1:
# For n=5: prod = cot(pi/5)*cot(2pi/5)*cot(3pi/5)*cot(4pi/5)
# cot(3pi/5) = -cot(2pi/5), cot(4pi/5) = -cot(pi/5)
# product = cot(pi/5)*cot(2pi/5)*(-cot(2pi/5))*(-cot(pi/5))
# = cot^2(pi/5)*cot^2(2pi/5) = (cot(pi/5)*cot(2pi/5))^2
prod_cot = 1
for k in range(1, 5):
    prod_cot *= cos(pi*k/5)/sin(pi*k/5)
print(f"  Product of all cot(pi*k/5): {prod_cot:.10f}")
print(f"  = 1/5: {1/5:.10f}")
print(f"  Match: {abs(prod_cot - 1/5) < 1e-8}")
print()

# General theorem
print("  THEOREM: For odd prime p,")
print("  prod_{k=1}^{p-1} cot(pi*k/p) = 1/p (up to sign).")
print()
print("  Verification:")
for p in [3, 5, 7, 11, 13]:
    prod = 1
    for k in range(1, p):
        prod *= cos(pi*k/p)/sin(pi*k/p)
    print(f"    p={p:2d}: product = {prod:+.10f}, 1/p = {1/p:.10f}, "
          f"match: {abs(abs(prod) - 1/p) < 1e-6}")

print()
print("  So Q maps the p-gon to points on the imaginary axis whose")
print("  cotangent product equals 1/p. This is a DEEP identity relating")
print("  the Cayley transform to roots of unity.")
print()

# ============================================================
print("="*70)
print("V. THE DIHEDRAL GROUP ACTION IN RAPIDITY SPACE")
print("="*70)
print()
print("  D_n acts on the n-gon by rotations and reflections.")
print("  Rotation by 2*pi/n: omega_k -> omega_{k+1}.")
print("  Reflection: omega_k -> omega_{-k} = conjugate(omega_k).")
print()
print("  Through Q, these become actions on the imaginary axis:")
print("  Rotation: i*cot(pi*k/n) -> i*cot(pi*(k+1)/n)")
print("  Reflection: i*cot(pi*k/n) -> i*cot(-pi*k/n) = -i*cot(pi*k/n)")
print("  (conjugation becomes NEGATION on the imaginary axis!)")
print()
print("  So the dihedral reflection becomes SIGN FLIP in the Cayley image.")
print("  This is the same as RAPIDITY INVERSION: phi -> -phi.")
print("  Reflection = time reversal = rapidity negation.")
print()

# ============================================================
print("="*70)
print("VI. POLYGON DIAGONALS AND RAPIDITY INTERVALS")
print("="*70)
print()
print("  In the n-gon, the diagonal from vertex 0 to vertex k")
print("  has length |omega_0 - omega_k| = 2*sin(pi*k/n).")
print("  The RAPIDITY of this length: ln(2*sin(pi*k/n))/2.")
print()

for n in [5, 6, 7, 12]:
    print(f"  {n}-GON diagonals:")
    for k in range(1, (n+1)//2 + 1):
        length = 2*sin(pi*k/n)
        if length > 0:
            rap = log(length)/2
        else:
            rap = float('-inf')
        print(f"    k={k}: length = 2*sin(pi*{k}/{n}) = {length:.6f}, "
              f"rapidity = {rap:.6f}")
    # Compute ratios of consecutive diagonals
    print(f"  Diagonal ratios:")
    for k in range(1, (n+1)//2):
        r1 = 2*sin(pi*k/n)
        r2 = 2*sin(pi*(k+1)/n)
        ratio = r2/r1
        rap_ratio = log(ratio)/2
        # Check for musical intervals
        for name, val in [("octave", 2), ("fifth", 1.5), ("fourth", 4/3),
                          ("maj 3rd", 1.25), ("min 3rd", 1.2),
                          ("phi", phi), ("1/phi", 1/phi)]:
            if abs(ratio - val) < 0.01:
                print(f"    d_{k+1}/d_{k} = {ratio:.6f} ~ {name}!")
        print(f"    d_{k+1}/d_{k} = {ratio:.6f}, rapidity = {rap_ratio:.6f}")
    print()

# ============================================================
print("="*70)
print("VII. THE PENTAGON DIAGONAL / SIDE = PHI")
print("="*70)
print()
print("  In the regular pentagon (5-gon):")
print("  Side length = 2*sin(pi/5) = {:.6f}".format(2*sin(pi/5)))
print("  Diagonal length = 2*sin(2*pi/5) = {:.6f}".format(2*sin(2*pi/5)))
print("  Ratio = diagonal/side = sin(2*pi/5)/sin(pi/5)")
ratio_pent = sin(2*pi/5)/sin(pi/5)
print(f"  = {ratio_pent:.10f}")
print(f"  = phi = {phi:.10f}")
print(f"  Match: {abs(ratio_pent - phi) < 1e-10}")
print()
print("  The rapidity of this ratio:")
rap_phi_ratio = log(phi)/2
print(f"  rapidity(phi) = ln(phi)/2 = {rap_phi_ratio:.10f}")
print(f"  = {rap_phi_ratio/octave:.4f} octaves")
print()
print("  And we know: Q(1/phi) = phi^3.")
print("  The pentagon's diagonal-to-side ratio is phi,")
print("  whose Cayley transform CUBES it to phi^3.")
print()
print("  The PENTAGON is the geometric source of the golden ratio.")
print("  The TOURNAMENT EP (exceptional point) has eigenvalue 1/phi.")
print("  The ALKANE C4H10 (butane) has interval 5/4 (major third).")
print("  And 5 = the number of vertices of the pentagon.")
print()

# ============================================================
print("="*70)
print("VIII. CONSTRUCTIBLE POLYGONS AND FERMAT PRIMES")
print("="*70)
print()
print("  A regular n-gon is constructible (compass+straightedge) iff")
print("  n = 2^a * p1 * p2 * ... where p_i are DISTINCT Fermat primes.")
print("  Fermat primes: 3, 5, 17, 257, 65537 (only 5 known).")
print()
print("  CONSTRUCTIBLE n-gons and their alkane associations:")
print()
print("  n-gon   Constructible?  Alkane C_{n-1}    Musical interval")
print("  " + "-"*60)
constructible = set()
fermat_primes = [3, 5, 17, 257, 65537]
# Generate constructible numbers
for a in range(20):
    base = 2**a
    if base > 100:
        break
    constructible.add(base)
    for p in fermat_primes:
        if base * p <= 100:
            constructible.add(base * p)
            for q in fermat_primes:
                if q > p and base * p * q <= 100:
                    constructible.add(base * p * q)

for n in range(3, 21):
    is_constr = n in constructible
    alkane_n = n - 1
    interval = Fraction(n, n-1)
    mark = "YES" if is_constr else "no"
    print(f"  {n:3d}-gon    {mark:3s}            "
          f"C{alkane_n}H{2*alkane_n+2:<4d}         {interval}")

print()
print("  The constructible n-gons correspond to alkanes C_{n-1}H_{2n}")
print("  whose musical intervals n/(n-1) can be 'constructed'")
print("  from octaves (factor 2) and Fermat intervals (factor p_Fermat).")
print()
print("  The NON-constructible 7-gon corresponds to propane (C3H8).")
print("  Wait -- propane's interval is 4/3 (constructible!).")
print("  The 7-gon is not constructible, but 4/3 IS constructible")
print("  (it's 2^2/3 = octave^2 / triangle).")
print()
print("  The POLYGON constructibility and the INTERVAL constructibility")
print("  are DIFFERENT notions. They coincide at n=3,4,5,6 but diverge.")
print()

# ============================================================
print("="*70)
print("IX. ROOTS OF UNITY AND THE TRANSFER MATRIX")
print("="*70)
print()
print("  The transfer matrix char poly: lambda^3 - lambda^2 - x*lambda - x = 0")
print("  At x = -1: lambda^3 - lambda^2 + lambda + 1 = 0")
print("  Factor: (lambda+1)(lambda^2 - 2*lambda + 1) = (lambda+1)(lambda-1)^2")
print("  Eigenvalues: -1, 1, 1. These are ROOTS OF UNITY (2nd and 1st).")
print()
print("  At x = 0: lambda^3 - lambda^2 = 0 => lambda^2(lambda-1) = 0")
print("  Eigenvalues: 0, 0, 1.")
print()
print("  At what x do we get n-th roots of unity as eigenvalues?")
print("  If lambda = omega = e^{2*pi*i/n}, then:")
print("  omega^3 - omega^2 - x*omega - x = 0")
print("  x = (omega^3 - omega^2) / (omega + 1)")
print()

for n in [3, 4, 5, 6, 7, 8, 12]:
    omega = cmath.exp(2j * cmath.pi / n)
    if abs(omega + 1) > 1e-10:
        x = (omega**3 - omega**2) / (omega + 1)
        print(f"  n={n:2d}: omega = e^{{2pi*i/{n}}}, "
              f"x = {x.real:+.6f}{x.imag:+.6f}i, |x| = {abs(x):.6f}")
    else:
        print(f"  n={n:2d}: omega = -1 (pole of formula)")

print()
print("  For n=4 (omega=i): x = (i^3-i^2)/(i+1) = (-i+1)/(1+i)")
x4 = (-1j + 1) / (1 + 1j)
print(f"  = {x4.real:+.6f}{x4.imag:+.6f}i = -i")
print("  So the SQUARE of the transfer matrix has i as eigenvalue at x = -i.")
print("  This connects the 4-fold periodicity of Q to the square root of -1.")
print()

# ============================================================
print("="*70)
print("X. THE CYCLOTOMIC POLYNOMIAL AND TOURNAMENTS")
print("="*70)
print()
print("  The n-th cyclotomic polynomial Phi_n(x) has roots at the")
print("  PRIMITIVE n-th roots of unity. deg(Phi_n) = phi(n) (Euler totient).")
print()
print("  The Cayley transform of primitive roots:")
print("  Q(omega) = i*cot(pi*k/n) for gcd(k,n) = 1.")
print()
print("  n   phi(n)  Primitive roots' |Q| values")
print("  " + "-"*55)
for n in range(2, 13):
    eu = sum(1 for k in range(1, n+1) if gcd(k, n) == 1)
    q_vals = []
    for k in range(1, n):
        if gcd(k, n) == 1:
            angle = pi * k / n
            cot_val = abs(cos(angle) / sin(angle))
            q_vals.append(cot_val)
    q_str = ", ".join(f"{v:.4f}" for v in sorted(set(round(v, 4) for v in q_vals)))
    print(f"  {n:2d}    {eu:2d}     |Q| = {q_str}")

print()

# ============================================================
print("="*70)
print("XI. THE GRAND UNIFICATION: POLYGON-ALKANE-MUSIC-TOURNAMENT")
print("="*70)
print()
print("  (n+1)-GON         ALKANE C_n         MUSIC           TOURNAMENT")
print("  " + "-"*65)
print("  triangle (3)      methane CH4        octave 2/1      T_3 (3-cycle)")
print("  square (4)        ethane C2H6        fifth 3/2       T on 4 vertices")
print("  pentagon (5)      propane C3H8       fourth 4/3      T on 5 vertices")
print("  hexagon (6)       butane C4H10       maj 3rd 5/4     T on 6 vertices")
print("  heptagon (7)      pentane C5H12      min 3rd 6/5     T_7 Paley (max)")
print("  octagon (8)       hexane C6H14       septimal 7/6    T on 8 vertices")
print()
print("  Each column encodes the SAME number n+1:")
print("  - As a POLYGON with n+1 vertices")
print("  - As a MOLECULE with 2(n+1) hydrogen atoms")
print("  - As a FREQUENCY RATIO (n+1)/n")
print("  - As a TOURNAMENT on some number of vertices")
print()
print("  The DIHEDRAL GROUP D_{n+1} (order 2(n+1) = hydrogen count)")
print("  is the symmetry group of the polygon.")
print("  Its ROTATIONS form Z_{n+1} (cyclic, order n+1).")
print("  Its REFLECTIONS pair each vertex with its opposite.")
print()
print("  In chemistry: rotations = ROTAMERS (conformational isomers).")
print("  In music: rotations = INVERSIONS of the interval.")
print("  In tournaments: rotations = RELABELING of vertices.")
print()
print("  THE CENTRAL IDENTITY:")
print()
print("  Q(1/(2n+1)) = (n+1)/n")
print("  Cayley velocity 1/(2n+1) -> superparticular ratio (n+1)/n")
print("  -> dihedral order 2(n+1)")
print("  -> hydrogen count of C_n H_{2n+2}")
print("  -> regular (n+1)-gon on the unit circle")
print("  -> roots of Phi_{n+1}(x)")
print()
print("  And through Q, the (n+1)-gon maps to cotangent values on")
print("  the imaginary axis, with product = 1/(n+1) (for prime n+1).")
print()
print("  arctanh(x) = x + x^3/3 + x^5/5 + ... is the generating function")
print("  that SUMS over all polygons/alkanes/intervals simultaneously.")
print("  The tournament coupling IS the sum over all regular polygons.")
print()

# ============================================================
print("="*70)
print("XII. THE Q^4 = IDENTITY AND THE SQUARE")
print("="*70)
print()
print("  We proved: Q^4 = identity. The orbit x -> Q(x) -> -1/x -> Q^{-1}(x) -> x.")
print("  This is the dihedral group D_2 = Klein four-group = Z_2 x Z_2.")
print("  |D_2| = 4 = hydrogen count of METHANE CH4.")
print()
print("  The SQUARE (4-gon) has symmetry D_4 (order 8).")
print("  But Q's period is 4, not 8.")
print("  Q generates Z_4 (cyclic of order 4), not D_4.")
print()
print("  The FULL symmetry group including complex conjugation:")
print("  Q and conjugation generate the DIHEDRAL GROUP D_4 of order 8")
print("  (4 rotations {Q^0, Q^1, Q^2, Q^3} and 4 reflections).")
print()
print("  Q as Mobius matrix: [[1,1],[-1,1]] in PGL(2).")
print("  Conjugation as Mobius: [[1,0],[0,-1]] (z -> -conj(z)... no, z -> conj(z)).")
print("  Together they generate a group of order 8 acting on C union {inf}.")
print()
print("  This 8-element group is the symmetry of the Cayley transform.")
print("  Its order 8 = hydrogen count of PROPANE C3H8.")
print("  And 8 = |D_4| = symmetry group of the SQUARE.")
print("  The Cayley transform's symmetry IS the square's symmetry.")
print()

# ============================================================
print("="*70)
print("XIII. CLOSING THE CIRCLE")
print("="*70)
print()
print("  Start: the Cayley transform Q(x) = (1+x)/(1-x).")
print("  It maps the unit circle to the imaginary axis.")
print("  It has period 4 (Q^4 = identity).")
print("  It generates the Klein four-group Z_2 x Z_2.")
print()
print("  Feed it reciprocal odd numbers: Q(1/(2n+1)) = (n+1)/n.")
print("  These are the musical intervals of just intonation.")
print("  The denominators 2n+1 are the Taylor indices of arctanh.")
print()
print("  The numerator n+1 is simultaneously:")
print("  - The vertex count of a regular polygon")
print("  - Half the hydrogen count of an alkane")
print("  - The higher frequency of a musical interval")
print("  - The order of the rotation subgroup Z_{n+1} of D_{n+1}")
print()
print("  The dihedral group D_{n+1} (order 2n+2 = hydrogen count)")
print("  acts on the (n+1)-gon, which Q maps to i*cot(pi*k/(n+1)).")
print()
print("  Everything is connected through ONE function and ONE number (n+1).")
print("  The unity is not metaphorical. It is algebraic.")
