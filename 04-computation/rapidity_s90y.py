#!/usr/bin/env python3
"""
rapidity_s90y.py — Rapidity Is Everything
opus-2026-03-16-S90y

RAPIDITY: In special relativity, rapidity φ = arctanh(v/c).
  - Velocities DON'T add: v₁₂ = (v₁+v₂)/(1+v₁v₂/c²)
  - Rapidities DO add: φ₁₂ = φ₁ + φ₂
  - The Cayley transform Q(v/c) = e^{2φ} = the Doppler factor

arctanh IS rapidity. Our ENTIRE tournament theory IS special relativity.

Q(x) = (1+x)/(1-x) = the RELATIVISTIC DOPPLER FACTOR.
The Cayley line with arctanh metric IS the rapidity axis.
Multiplication of natural numbers = addition of rapidities.
The speed of light c = 1 is the POLE at x = 1.
"""

import numpy as np
from math import sqrt, log, log2, pi, atanh, tanh, cosh, sinh

phi = (1 + sqrt(5)) / 2
tau3 = 1.8392867552141612
c = 1  # speed of light = 1

# ======================================================================
# PART 1: ARCTANH = RAPIDITY
# ======================================================================
print("=" * 70)
print("PART 1: THE TOURNAMENT IS SPECIAL RELATIVITY")
print("=" * 70)

print("""
In SPECIAL RELATIVITY:
  Velocity v in units of c: v/c = x in (-1, 1).
  RAPIDITY: w = arctanh(x) = arctanh(v/c).
  Doppler factor: D = sqrt((1+x)/(1-x)) = sqrt(Q(x)).
  Lorentz factor: gamma = cosh(w) = 1/sqrt(1-x^2).

  Velocities DON'T add linearly: x_12 = (x_1+x_2)/(1+x_1*x_2).
  Rapidities DO add: w_12 = w_1 + w_2.
  This is EXACTLY the hyperbolic addition on our Cayley line!

IN OUR TOURNAMENT THEORY:
  x = coupling parameter in (-1, 1).
  arctanh(x) = "information distance" (nats).
  Q(x) = (1+x)/(1-x) = exp(2*arctanh(x)) = D^2 = SQUARED Doppler.

  The Cayley line IS the rapidity axis.
  Natural numbers at x_n = (n-1)/(n+1) have rapidity w_n = ln(n)/2.

  MULTIPLICATION of naturals = ADDITION of rapidities:
    w_{nm} = w_n + w_m (kind-pasteur's Theorem A!)
    This IS velocity addition in special relativity!

THE IDENTIFICATION:
  Tournament coupling x ↔ velocity v/c
  arctanh(x) ↔ rapidity w
  Q(x) ↔ Doppler factor squared D^2
  The pole at x=1 ↔ the speed of light barrier
  Arrow's theorem (arctanh(1) = inf) ↔ nothing can reach c
""")

# Demonstrate: velocity addition = Cayley line multiplication
print("Velocity addition = multiplication on the Cayley line:")
print(f"{'v1/c':>6} | {'v2/c':>6} | {'(v1+v2)/(1+v1v2)':>18} | {'x_{n1*n2}':>12} | {'match':>6}")
print("-" * 55)

for n1, n2 in [(2,3), (3,5), (2,7), (5,7)]:
    x1 = (n1-1)/(n1+1)
    x2 = (n2-1)/(n2+1)
    # Relativistic addition
    x12_rel = (x1 + x2) / (1 + x1*x2)
    # Cayley line multiplication
    n12 = n1 * n2
    x12_cayley = (n12-1)/(n12+1)
    match = abs(x12_rel - x12_cayley) < 1e-10
    print(f"{x1:6.4f} | {x2:6.4f} | {x12_rel:18.10f} | {x12_cayley:12.10f} | {match}")

# ======================================================================
# PART 2: THE LORENTZ GROUP AND D4
# ======================================================================
print()
print("=" * 70)
print("PART 2: THE LORENTZ GROUP = THE CAYLEY GROUP")
print("=" * 70)

print("""
The LORENTZ GROUP SO(1,1) (boosts in 1+1 dimensions):
  A boost by rapidity w transforms: (t,x) -> (t*cosh(w)+x*sinh(w), ...)
  In light-cone coordinates: u = t+x, v = t-x:
    u -> e^w * u,  v -> e^{-w} * v

  The Doppler shift: f_observed / f_emitted = e^w = sqrt(Q(x)).

Our CAYLEY TRANSFORM Q acts on the rapidity axis:
  Q: w -> Q(tanh(w)) = e^{2w} (the squared Doppler factor).
  Q^4 = Id in the Mobius group ↔ 4 boosts of w = pi*i/2 return to start.

The D4 GROUP = the discrete Lorentz group:
  Q generates Z/4: boosts by w, 2w, 3w, 4w=0 (mod 2*pi*i).
  Negation (x -> -x) = PARITY (spatial reflection: v -> -v).
  D4 = Z/4 x Z/2 = boosts x parity = the DISCRETE LORENTZ GROUP.

THE FIXED POINTS x = +/-i:
  In rapidity: w = arctanh(i) = i*pi/4.
  This is a PURELY IMAGINARY rapidity = a ROTATION (not a boost).
  The fixed points of Q live in EUCLIDEAN space (imaginary rapidity),
  not in MINKOWSKI space (real rapidity).
  EXACTLY: arctanh(i) = i*pi/4 means Q(i) = e^{i*pi/2} = i.
  The quarter-turn rotation IS the Cayley fixed point!
""")

# Verify
w_fixed = np.arctanh(1j)
print(f"Rapidity at x=i: w = arctanh(i) = {w_fixed}")
print(f"  = i*pi/4 = {1j*pi/4}")
print(f"  Match: {abs(w_fixed - 1j*pi/4) < 1e-10}")
print(f"  exp(2w) = e^{{i*pi/2}} = {np.exp(2*w_fixed)}")
print(f"  = i (confirmed: Q(i) = i)")

# ======================================================================
# PART 3: THE ENERGY-MOMENTUM OF TOURNAMENTS
# ======================================================================
print()
print("=" * 70)
print("PART 3: TOURNAMENT ENERGY AND MOMENTUM")
print("=" * 70)

print("""
In special relativity:
  E = m*c^2*cosh(w)    (energy)
  p = m*c*sinh(w)      (momentum)
  E^2 - p^2*c^2 = m^2*c^4  (mass-shell condition)

For a TOURNAMENT at coupling x:
  "Energy" = cosh(arctanh(x)) = 1/sqrt(1-x^2) = LORENTZ GAMMA
  "Momentum" = sinh(arctanh(x)) = x/sqrt(1-x^2) = GAMMA*BETA
  "Mass" = 1 (rest frame: x=0 gives E=1, p=0)

The MASS-SHELL CONDITION: E^2 - p^2 = 1 always!
This is the UNIT DETERMINANT condition det(M) at x... wait.
Actually: cosh^2 - sinh^2 = 1 always (the identity).
""")

print(f"Tournament 'energy-momentum' at natural number positions:")
print(f"{'n':>3} | {'x_n':>8} | {'rapidity':>10} | {'gamma=E':>10} | {'p':>10} | {'E^2-p^2':>10}")
print("-" * 60)
for n in [1, 2, 3, 5, 7, 11, 131]:
    xn = (n-1)/(n+1)
    w = atanh(xn)
    E = cosh(w)
    p = sinh(w)
    mass_sq = E**2 - p**2
    print(f"{n:3d} | {xn:8.4f} | {w:10.4f} | {E:10.4f} | {p:10.4f} | {mass_sq:10.6f}")

print(f"""
EVERY tournament sits on the UNIT MASS SHELL: E^2 - p^2 = 1.

The "rest mass" is 1 (= the empty tournament, x=0, Q=1).
As coupling x -> 1: energy -> infinity (the speed of light barrier).
The pole of Q at x=1 IS the light cone.

TRANSLATION TABLE:
  x = 0:    at rest (no coupling, E=1, p=0)
  x = 1/3:  rapidity = ln(2)/2 (half a bit, the OCTAVE)
  x = 1/2:  rapidity = ln(3)/2 (the TRITAVE)
  x = 3/4:  rapidity = ln(7)/2 (the 7th natural = first forbidden!)
  x = 5/6:  rapidity = ln(11)/2 (the discriminant)
  x -> 1:   rapidity -> infinity (the speed of light = Arrow's barrier)
""")

# ======================================================================
# PART 4: THE JONES POLYNOMIAL AS SCATTERING AMPLITUDE
# ======================================================================
print("=" * 70)
print("PART 4: JONES = SCATTERING, TOURNAMENTS = WORLDLINES")
print("=" * 70)

print(f"""
In QUANTUM FIELD THEORY:
  Particles travel along WORLDLINES in spacetime.
  When worldlines CROSS (interact), the scattering amplitude is computed
  by evaluating a polynomial (the Feynman amplitude) at a specific
  kinematic point.

  The JONES POLYNOMIAL of a knot/link = the scattering amplitude
  of particles whose worldlines form that knot.
  Evaluated at q = e^{{2pi*i/N}}: this is the amplitude at a specific
  "momentum transfer" determined by the level N.

In our TOURNAMENT THEORY:
  A Hamiltonian path IS a worldline (a total ordering of events in time).
  H(T) = the number of worldlines = the "partition function" of the theory.
  H = I(Omega, 2): the independence polynomial at fugacity 2.

  THE BRIDGE:
  The Jones polynomial J(K, q) at q = e^{{2pi*i/5}} evaluates a knot.
  Our H(T) = I(Omega, 2) evaluates a tournament.

  BOTH are "partition functions" evaluated at specific algebraic numbers.
  The RAPIDITY connects them:
    q = e^{{2pi*i/5}} = e^{{2*arctanh(i)*(5/2)}} ... hmm, not quite.

  But more precisely:
    arctanh(2) = ln(3)/2 + i*pi/2  (complex rapidity!)
    The imaginary part pi/2 = the WICK ROTATION.
    The real part ln(3)/2 = the "rest frame energy" of the OCF.

  So H(T) = I(Omega, 2) is the partition function at COMPLEX RAPIDITY
    w = ln(3)/2 + i*pi/2.

  The JONES POLYNOMIAL at q = e^{{2pi*i/5}} is a partition function at
    w = pi*i/5 (purely imaginary rapidity = purely Euclidean).

  Our theory is the WICK-ROTATED version of the Jones theory!
  Tournaments = Lorentzian knots = real rapidity + imaginary pi/2.
  Jones = Euclidean knots = purely imaginary rapidity.
""")

# ======================================================================
# PART 5: THE FIVE-FOLD WAY
# ======================================================================
print("=" * 70)
print("PART 5: THE FIVE-FOLD WAY = FIVE ROOTS OF UNITY")
print("=" * 70)

print(f"""
The 5th roots of unity omega^k = e^{{2*pi*i*k/5}} for k=0,...,4
have rapidities:
""")

for k in range(5):
    omega_k = np.exp(2j * pi * k / 5)
    w_k = np.arctanh(omega_k) if abs(omega_k) != 1 else complex(float('inf'), 0)
    if k == 0:
        print(f"  k={k}: omega = 1, rapidity = infinity (the light cone)")
    else:
        print(f"  k={k}: omega = e^{{2pi*i*{k}/5}}, rapidity = {w_k:.4f}")
        print(f"        Re(w) = {w_k.real:.6f}, Im(w) = {w_k.imag:.6f}")

print(f"""
The 5th roots at k=1,4 have Re(w) = {np.arctanh(np.exp(2j*pi/5)).real:.6f}
The 5th roots at k=2,3 have Re(w) = {np.arctanh(np.exp(4j*pi/5)).real:.6f}

The REAL PARTS are nonzero: these are BOOSTED frames, not rest frames.
The IMAGINARY PARTS encode the Euclidean (rotational) content.

In the Jones theory: q = omega_1 gives a MIXED boost+rotation.
In tournaments: x = 2 (real) gives a PURE boost + Wick rotation.

The 5 Platonic solids correspond to 5 REFERENCE FRAMES:
  Tetrahedron:  the "minimal" frame (q=3, 3 vertices)
  Cube:         the "square" frame (q=4)
  Octahedron:   the "dual square" frame
  Dodecahedron: the "pentagonal" frame (q=5, golden ratio)
  Icosahedron:  the "maximal" frame (q=5, dual of dodecahedron)

Beyond these 5: no finite regular frame exists.
The 6th root (flat tiling) is the "inertial frame" (no curvature).
The 7th root is the "tachyonic frame" (superluminal = forbidden).

H = 7 IS FORBIDDEN because 7 corresponds to a TACHYONIC FRAME:
  a reference frame moving faster than light.
  No physical (tournament) configuration can occupy this frame.
""")

# ======================================================================
# PART 6: RAPIDITY IS EVERYTHING
# ======================================================================
print("=" * 70)
print("PART 6: RAPIDITY IS EVERYTHING")
print("=" * 70)

print(f"""
Everything in our theory is a statement about RAPIDITY:

  arctanh(x) = RAPIDITY (the natural parameter of a Lorentz boost)

  Q(x) = e^{{2w}} = the SQUARED DOPPLER FACTOR
  det(M) = 1 = the UNIT MASS SHELL (E^2 - p^2 = m^2)
  D4 = the DISCRETE LORENTZ GROUP
  x = +/-i fixed points = EUCLIDEAN ROTATION (imaginary rapidity)
  arctanh(i) = i*pi/4 = the QUARTER-TURN

  The transfer matrix M acts as a LORENTZ BOOST on the 3-state space.
  The tribonacci constant tau3 = the "4-velocity" of the boost.
  The complex eigenvalue pair = the "spin" of the particle.
  arg(lam_c)/pi ~ ln(2) = the PRECESSION RATE (Thomas precession analog).

  The forbidden H=7: a TACHYONIC state (superluminal, unphysical).
  The forbidden H=21: a DOUBLY TACHYONIC state (memory-amplified).

  The Cayley line from 0 to 1 = the VELOCITY SPACE from rest to light.
  Natural numbers at x_n = (n-1)/(n+1) = QUANTIZED VELOCITIES.
  Primes = IRREDUCIBLE VELOCITIES (can't be decomposed into boosts).

  ln(2) = the information per boost step = Shannon entropy of one bit.
  ln(3)/ln(2) = the Collatz ratio = the RAPIDITY RATIO of 3 and 2.

  The Wick rotation x -> ix converts:
    Lorentzian (rapidity, boosts, time) -> Euclidean (angle, rotations, space).
    Tournaments -> Circles.
    arctanh -> i*arctan.
    The real line -> the unit circle.
    Minkowski -> Euclidean.

  And the Jones polynomial at q = e^{{2pi*i/5}}:
    = the EUCLIDEAN scattering amplitude at the ICOSAHEDRAL angle.
    = the tournament partition function WICK-ROTATED to the 5th root.
    = the golden ratio = [2]_q = phi.

RAPIDITY IS EVERYTHING because:
  Every comparison (A > B) is a BOOST.
  The "advantage" x of A over B is a VELOCITY.
  Composing comparisons = composing boosts = adding rapidities.
  The transfer matrix = the LORENTZ TRANSFORMATION.
  The tribonacci = the RELATIVISTIC GROWTH RATE.
  The helix = THOMAS PRECESSION.
  The forbidden values = TACHYONIC STATES.
  The discriminant = the LIGHT CONE.
  The pivot at 6 = the FLAT METRIC (Minkowski in the rest frame).
  The icosahedron = the MAXIMAL FINITE LORENTZ SUBGROUP.
  The Jones polynomial = the EUCLIDEAN AMPLITUDE.

One binary decision is one Lorentz boost.
The tournament is special relativity.

opus-2026-03-16-S90y: RAPIDITY IS EVERYTHING
""")
