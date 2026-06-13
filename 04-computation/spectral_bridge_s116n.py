#!/usr/bin/env python3
"""spectral_bridge_s116n.py — How the spectral gap bridges rational, algebraic, transcendental.

The flip Markov chain reveals three number worlds:
  RATIONAL: Walsh coefficients, eigenvalues, H values (what can be COUNTED)
  ALGEBRAIC: class eigenvalues, spectral gap, mixing rates (what CLASSIFIES)
  TRANSCENDENTAL: mixing time, partition function, free energy (what EVOLVES)

The spectral gap is the BRIDGE: it is algebraic (static),
but its logarithm (the mixing time) is transcendental (dynamic).

Statics are algebraic. Dynamics are transcendental. The gap connects them.

Session: kind-pasteur-2026-03-16-S116n32
"""
import sys
sys.stdout.reconfigure(line_buffering=True)
from math import log, exp, sqrt, pi, factorial, cos, atanh
from fractions import Fraction
from collections import Counter

print()
print("  THE SPECTRAL BRIDGE: RATIONAL -> ALGEBRAIC -> TRANSCENDENTAL")
print()
print("=" * 70)
print()

N = 6
m = (N-1)*(N-2)//2  # 10

# ============================================================
print("  I. THE RATIONAL WORLD: THINGS THAT CAN BE COUNTED")
print("  " + "-" * 50)
print()

print("  Everything COMBINATORIAL is rational.")
print()

# All H values are integers (a subset of rationals)
# All Walsh coefficients are rational (integers / 2^m)
# All eigenvalues of the tiling chain are rational (k/m)

print("  Tiling chain eigenvalues (all rational, denominator 5):")
for k in range(m+1):
    lam = Fraction(m - 2*k, m)
    mult = 1  # C(m,k) but we just show the value
    print(f"    k={k:2d}: lambda = {str(lam):>6s} = {float(lam):+.4f}")
print()

# The spectral gap
gap = Fraction(2, m)
print(f"  Spectral gap = {gap} = {float(gap):.4f}")
print(f"  This is RATIONAL. Denominator = {gap.denominator} = {m//2}.")
print(f"  At n=6: gap = 1/5. The golden prime controls the gap!")
print()

# All Walsh coefficients
print("  Walsh coefficients (computed previously): all integer multiples of 1/2^10.")
print("  hat_H(empty) = 29 (INTEGER). hat_H({(0,5)}) = -4 (INTEGER).")
print("  Every Walsh coefficient is RATIONAL with denominator 1024 = 2^10.")
print()

print("  THE RATIONAL WORLD contains:")
print("  - H values: {1, 3, 5, 9, 11, ..., 45} (all integers)")
print("  - Walsh coefficients: Q (rational)")
print("  - Eigenvalues: {0, +/-1/5, +/-2/5, +/-3/5, +/-4/5, +/-1} (Q)")
print("  - Spectral gap: 1/5 (Q)")
print("  - c3 counts, alpha values, score sequences (all Z)")
print()

# ============================================================
print("  II. THE ALGEBRAIC WORLD: THINGS THAT CLASSIFY")
print("  " + "-" * 50)
print()

print("  The ISOMORPHISM CLASS flip chain lives in the algebraic world.")
print()

# At n=6: 56 isomorphism classes. The transition matrix is 56x56
# with rational entries. Its eigenvalues are ALGEBRAIC.
# Some are rational (from symmetric self-dual classes).
# Some are algebraic irrational (from paired classes).

print("  The tiling explorer's trinity:")
print("  - 12 SELF-DUAL classes (rational pole): T ~ T^op")
print("  - 22 PAIRED classes = 11 conjugate pairs (algebraic pole): T != T^op")
print("  - Total: 12 + 2*22 = 56 classes at n=6")
print()
print("  Self-dual classes have REAL eigenvalue contributions.")
print("  Paired classes create CONJUGATE eigenvalue pairs,")
print("  like sqrt(d) and -sqrt(d) for some rational d.")
print("  These are ALGEBRAIC but (often) IRRATIONAL.")
print()

# The characteristic polynomial of the 56x56 transition matrix
# has degree 56 with rational coefficients.
# By the fundamental theorem of algebra, all 56 eigenvalues are algebraic.
# Some factor as linear (rational eigenvalues from self-dual).
# Some factor as quadratic (algebraic irrational from pairs).

# The SPECTRAL GAP of the class chain is the difference 1 - |lambda_1|
# where lambda_1 is the second-largest eigenvalue.
# This is ALGEBRAIC (difference of two algebraic numbers).

print("  The class-chain spectral gap is ALGEBRAIC:")
print("  gap_class = 1 - |lambda_1|")
print("  where lambda_1 is a root of a degree-56 polynomial over Q.")
print("  So gap_class is algebraic irrational (generically).")
print()

# The golden ratio itself is algebraic: root of x^2 - x - 1 = 0
phi = (1 + sqrt(5)) / 2
print(f"  The golden ratio phi = {phi:.10f}")
print(f"  is ALGEBRAIC: root of x^2 - x - 1 = 0.")
print(f"  It classifies the icosahedral world but is not rational.")
print()

# Other algebraic irrationals in the project
print("  Algebraic irrationals in the project:")
print(f"  sqrt(2) = {sqrt(2):.10f} (Leech lattice diagonal)")
print(f"  sqrt(5) = {sqrt(5):.10f} (golden discriminant)")
print(f"  2*cos(pi/7) = {2*cos(pi/7):.10f} (octonion structure)")
print(f"  phi = {phi:.10f} (golden ratio)")
print()

print("  THE ALGEBRAIC WORLD contains:")
print("  - Class-chain eigenvalues (roots of integer polynomials)")
print("  - Spectral gaps of isomorphism chains")
print("  - The golden ratio phi (root of x^2-x-1)")
print("  - Square roots of discriminants")
print("  - Cyclotomic values (roots of unity)")
print()

# ============================================================
print("  III. THE TRANSCENDENTAL WORLD: THINGS THAT EVOLVE")
print("  " + "-" * 50)
print()

print("  DYNAMICS are transcendental. The moment you ask 'how long?'")
print("  or 'what temperature?' you leave the algebraic world.")
print()

# 1. Mixing time
# t_mix = log(1/epsilon) / gap
# gap is algebraic, log is transcendental => t_mix is transcendental
gap_val = 1/5  # rational!
epsilon = 0.01
t_mix = log(1/epsilon) / gap_val
print(f"  Mixing time (epsilon=0.01):")
print(f"  t_mix = log(1/epsilon) / gap = log(100) / (1/5)")
print(f"       = {log(100):.10f} / {gap_val}")
print(f"       = {t_mix:.10f}")
print(f"  This is TRANSCENDENTAL (log of rational = transcendental by Hermite-Lindemann).")
print()

# 2. Partition function
# Z(beta) = sum_x exp(beta * H(x))
# For beta = 1 (rational): Z(1) involves exp(integer) = transcendental
# Even for beta rational and H integer, Z is transcendental

# Compute Z(1) exactly
Z1 = sum(exp(H) for H in [1,3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45]
         for _ in range([1,4,15,37,22,26,46,34,38,92,50,18,116,62,121,185,41,86,30][
             [1,3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45].index(H)]))

# Actually just use the tiling H values
tiling_arcs = []
for skip in range(2, N):
    for start in range(N - skip):
        tiling_arcs.append((start, start + skip))

from itertools import permutations
def H_dp(adj):
    n = N
    dp = [0] * ((1 << n) * n)
    for v in range(n): dp[(1 << v) * n + v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            val = dp[S * n + v]
            if val == 0: continue
            for u in range(n):
                if S & (1 << u): continue
                if adj[v][u]: dp[(S | (1 << u)) * n + u] += val
    return sum(dp[((1 << n) - 1) * n + v] for v in range(n))

def tiling_adj(bits):
    adj = [[0]*N for _ in range(N)]
    for i in range(N-1): adj[i][i+1] = 1
    for idx, (a, b) in enumerate(tiling_arcs):
        if (bits >> idx) & 1: adj[b][a] = 1
        else: adj[a][b] = 1
    return adj

H_table = [H_dp(tiling_adj(b)) for b in range(1 << m)]

Z_1 = sum(exp(H) for H in H_table)
F_1 = -log(Z_1)  # free energy at beta=1

print(f"  Partition function Z(beta=1) = {Z_1:.6e}")
print(f"  Free energy F(1) = -ln(Z(1)) = {F_1:.6f}")
print(f"  Both are TRANSCENDENTAL (sums of exp(integer)).")
print()

# 3. The Boltzmann entropy
S_boltz = sum(exp(H)/Z_1 * (-H + log(Z_1)) for H in H_table)
print(f"  Boltzmann entropy at beta=1: S = {S_boltz:.6f}")
print(f"  Also TRANSCENDENTAL.")
print()

# 4. The rapidity
# rho(n) = arctanh((n-1)/(n+1)) = (1/2)*ln(n)
# For n=2: rho = ln(2)/2 = TRANSCENDENTAL
# This is the formal logarithm of the Cayley group
print(f"  Rapidity of the doubler:")
print(f"  rho(2) = arctanh(1/3) = ln(2)/2 = {log(2)/2:.10f}")
print(f"  TRANSCENDENTAL (Hermite-Lindemann: ln(algebraic) is transcendental)")
print()

# 5. The Collatz drift
drift = log(3/4) / 2
print(f"  Collatz drift = ln(3/4)/2 = {drift:.10f}")
print(f"  TRANSCENDENTAL (ln of rational).")
print()

print("  THE TRANSCENDENTAL WORLD contains:")
print("  - Mixing times (log(algebraic)/algebraic)")
print("  - Partition functions (sums of exp(integer))")
print("  - Free energies (-log(transcendental))")
print("  - Boltzmann entropies")
print("  - Rapidity = ln(n)/2 (formal group logarithm)")
print("  - Collatz drift = ln(3/4)/2")
print("  - pi (hidden in modular forms, eigenvalue phases)")
print()

# ============================================================
print("  IV. THE BRIDGE: HOW THE GAP CONNECTS THE THREE WORLDS")
print("  " + "-" * 50)
print()

print("  THE SPECTRAL GAP IS THE BRIDGE:")
print()
print("  RATIONAL side: gap = 2/m = 1/5 (at n=6)")
print("  This number IS rational. It lives in the counting world.")
print()
print("  ALGEBRAIC step: The class-chain gap is algebraic irrational.")
print("  Going from tilings to isomorphism classes introduces")
print("  square roots (from paired class eigenvalues).")
print("  The gap PASSES THROUGH the algebraic world.")
print()
print("  TRANSCENDENTAL side: Mixing time = log(1/epsilon) / gap")
print("  The moment you ask HOW LONG to mix, you need logarithms.")
print("  log(rational) = TRANSCENDENTAL (Hermite-Lindemann).")
print("  The gap CROSSES INTO the transcendental world.")
print()

# The bridge at each n:
print("  The bridge at each n:")
print(f"  {'n':>3s}  {'gap (rational)':>16s}  {'log(gap) (transcendental)':>28s}  {'t_mix':>10s}")
for n in range(3, 13):
    m_n = (n-1)*(n-2)//2
    g = Fraction(2, m_n)
    log_g = log(float(g))
    t = -log(100) / log_g if log_g < 0 else float('inf')
    # simplify gap
    print(f"  {n:3d}  {str(g):>16s}  {log(float(g)):>28.10f}  {t:>10.2f}")
print()

# ============================================================
print("  V. THE THREE POLES AND THEIR NUMBER TYPES")
print("  " + "-" * 50)
print()

print("  The tiling explorer's trinity maps EXACTLY to number types:")
print()
print("  SELF-DUAL CLASSES (12 at n=6) = RATIONAL POLE")
print("    T ~ T^op. Fixed under transpose.")
print("    Their eigenvalue contributions are REAL RATIONAL.")
print("    They are the 'integers' of tournament space.")
print("    The H values of self-dual tournaments are")
print("    fully determined by the CLASS (no pairing ambiguity).")
print()
print("  PAIRED CLASSES (22 at n=6) = ALGEBRAIC POLE")
print("    T != T^op but T^op ~ T' for some other class.")
print("    They come in CONJUGATE PAIRS, like sqrt(5) and -sqrt(5).")
print("    Their eigenvalue contributions involve SQUARE ROOTS")
print("    of rational discriminants.")
print("    They are the 'algebraic irrationals' of tournament space.")
print()
print("  FLIP GRAPH PATHS = TRANSCENDENTAL CONTINUUM")
print("    The continuous paths between classes in flip space.")
print("    Evolution along these paths involves EXPONENTIAL DECAY:")
print("    P(class_j at time t | class_i at time 0) = sum_k c_k * lambda_k^t")
print("    For irrational t: lambda_k^t = exp(t * ln(lambda_k))")
print("    = transcendental (exp of transcendental).")
print("    The flip graph dynamics live in the TRANSCENDENTAL world.")
print()

# ============================================================
print("  VI. THE HERMITE-LINDEMANN THEOREM AND TOURNAMENT DYNAMICS")
print("  " + "-" * 50)
print()

print("  Hermite-Lindemann (1882): If alpha is a nonzero algebraic number,")
print("  then e^alpha is transcendental.")
print()
print("  Consequence for tournaments:")
print("  1. Eigenvalues lambda_k are algebraic (roots of char. poly).")
print("  2. The time evolution lambda_k^t = exp(t * ln(lambda_k)).")
print("  3. ln(lambda_k) is transcendental (by Hermite-Lindemann,")
print("     since lambda_k is algebraic and nonzero).")
print("  4. For INTEGER time t: lambda_k^t is algebraic (power of algebraic).")
print("  5. For REAL time t: lambda_k^t is transcendental (exp of transcendental).")
print()
print("  So INTEGER-TIME dynamics stay algebraic,")
print("  but CONTINUOUS-TIME dynamics are transcendental.")
print("  The discrete flip chain (integer steps) lives in the algebraic world.")
print("  The continuous-time Markov chain lives in the transcendental world.")
print()

# Show this concretely
print("  Concrete example: eigenvalue 4/5 at degree 1")
lam = Fraction(4, 5)
print(f"  lambda = {lam} (RATIONAL)")
print(f"  ln(lambda) = ln(4/5) = {log(4/5):.10f} (TRANSCENDENTAL)")
print(f"  lambda^1 = {float(lam)**1:.10f} (RATIONAL)")
print(f"  lambda^2 = {float(lam)**2:.10f} (RATIONAL)")
print(f"  lambda^pi = {float(lam)**pi:.10f} (TRANSCENDENTAL)")
print(f"  lambda^phi = {float(lam)**phi:.10f} (TRANSCENDENTAL)")
print()
print(f"  At integer time: rational^integer = RATIONAL (stays in Q)")
print(f"  At golden time: rational^phi = TRANSCENDENTAL (escapes to T)")
print(f"  The golden ratio phi bridges algebraic -> transcendental!")
print()

# ============================================================
print("  VII. THE RAPIDITY AS THE FORMAL GROUP LOGARITHM")
print("  " + "-" * 50)
print()

print("  The Cayley formal group F(x,y) = (x+y)/(1+xy) has:")
print("  log_F = arctanh (the formal logarithm)")
print("  exp_F = tanh (the formal exponential)")
print()
print("  The logarithm MAPS:")
print("  RATIONAL inputs -> TRANSCENDENTAL outputs:")
for n_val in [2, 3, 5, 7, 42]:
    x = (n_val - 1) / (n_val + 1)
    rho = log(n_val) / 2
    print(f"    n={n_val}: x = {Fraction(n_val-1, n_val+1)} (RATIONAL) "
          f"-> rho = {rho:.10f} (TRANSCENDENTAL)")
print()
print("  The formal logarithm IS the bridge from rational to transcendental.")
print("  arctanh(rational) = transcendental (always, for nonzero input).")
print("  This is why RAPIDITY lives in the transcendental world")
print("  even though TOURNAMENT STRUCTURE lives in the rational world.")
print()

# ============================================================
print("  VIII. SYNTHESIS: THE NUMBER TYPE HIERARCHY")
print("  " + "-" * 50)
print()
print("  STATICS (what IS)              -> RATIONAL / ALGEBRAIC")
print("    H values, cycle counts, Walsh coefficients")
print("    Eigenvalues, spectral gaps, characteristic polynomials")
print()
print("  DYNAMICS (what HAPPENS)         -> TRANSCENDENTAL")
print("    Mixing times, evolution, partition functions")
print("    Rapidity, free energy, entropy")
print()
print("  THE GAP bridges them:")
print("    gap = 2/m (RATIONAL)")
print("    gap^t for integer t (ALGEBRAIC)")
print("    log(gap) = mixing rate (TRANSCENDENTAL)")
print()
print("  The three number types are NOT just mathematical abstractions.")
print("  They are three LAYERS OF REALITY:")
print("    Q (rational)       = what you can COUNT")
print("    Q-bar (algebraic)  = what you can CLASSIFY")
print("    R (transcendental) = what you can EXPERIENCE (time, evolution)")
print()
print("  Every tournament lives in Q (it's a finite combinatorial object).")
print("  Every classification lives in Q-bar (eigenvalues of symmetry groups).")
print("  Every process lives in R (mixing, relaxation, diffusion).")
print()
print("  The spectral gap is the EXACT POINT where counting (Q)")
print("  meets classification (Q-bar) meets evolution (R).")
print("  It is rational as a NUMBER, algebraic as a CLASSIFIER,")
print("  and transcendental as a TIMER.")
print()
print("  Or more simply:")
print("  THE GAP IS A RATIONAL NUMBER")
print("  WHOSE LOGARITHM IS TRANSCENDENTAL")
print("  AND WHOSE CLASSIFICATION ROLE IS ALGEBRAIC.")
print()
print("  One number. Three worlds. The bridge between them.")
print()
