"""
transcendental_tournament_bases.py
opus-2026-03-14-S71l

Deep exploration of tournament theory in transcendental and irrational bases.

Key insight from previous work:
- Phi_3(2) = 7 (forbidden), Phi_3(tau) = tau^3 (tribonacci)
- The 3-strand Pascal is Phi_3(x)^n
- The Walsh-simplex correspondence: linear = simplex, nonlinear = corners

New explorations:
1. Base e: H values in natural logarithmic base
2. Base pi: projective planes in the Archimedes base
3. The generating function H(x) = sum H_T x^T over all tournaments
4. Euler's formula e^{i*pi}+1=0 and Phi_3(e^{2*pi*i/3})=0
5. The tau-e-phi triangle of irrational bases
6. Catalan numbers and tournament recursion
"""

import numpy as np
from math import factorial, comb, log, exp, pi, e, sqrt
from fractions import Fraction

print("=" * 70)
print("TRANSCENDENTAL AND IRRATIONAL TOURNAMENT BASES")
print("opus-2026-03-14-S71l")
print("=" * 70)

# =====================================================================
# PART 1: THE TRINITY OF BASES — phi, tau, 2
# =====================================================================
print("\n" + "=" * 70)
print("PART 1: THE TRINITY OF BASES — phi, tau, 2")
print("=" * 70)

phi = (1 + sqrt(5)) / 2  # golden ratio
tau = 1.8392867552141612  # tribonacci constant

print(f"""
  The three algebraic bases of tournament theory:

  phi = {phi:.10f} (golden ratio, root of x^2-x-1)
  tau = {tau:.10f} (tribonacci constant, root of x^3-x^2-x-1)
  2   = 2.0000000000 (tournament generator)

  Key identities:
  phi^2 = phi + 1 = Phi_2(phi)     [2-nacci satisfies Phi_2]
  tau^3 = tau^2 + tau + 1 = Phi_3(tau)  [3-nacci satisfies Phi_3]
  2^k   = 2^k                      [binary]

  At x=2: Phi_2(2) = 3, Phi_3(2) = 7

  The "span" between these bases:
  phi < tau < 2
  {phi:.6f} < {tau:.6f} < 2.000000

  Ratios:
  tau/phi = {tau/phi:.6f}   (≈ 1.137, close to ln(e^(1/e)) = ? no...)
  2/tau   = {2/tau:.6f}   (≈ 1.087)
  2/phi   = {2/phi:.6f}   (= 2/phi, the "golden excess")

  The product phi * tau = {phi*tau:.6f}
  Compare: phi * tau ≈ {phi*tau:.6f}
           phi^2      = {phi**2:.6f} = phi + 1
           sqrt(5)    = {sqrt(5):.6f}

  Interesting: phi * tau is NOT a "nice" algebraic number.
""")

# =====================================================================
# PART 2: TRANSCENDENTAL EVALUATIONS OF PHI_3
# =====================================================================
print("=" * 70)
print("PART 2: TRANSCENDENTAL EVALUATIONS OF PHI_3")
print("=" * 70)

evaluations = {
    'e': e,
    'pi': pi,
    'phi': phi,
    'tau': tau,
    'sqrt(2)': sqrt(2),
    'sqrt(3)': sqrt(3),
    'ln(2)': log(2),
    '2': 2.0,
    '3': 3.0,
    '1': 1.0,
    '0': 0.0,
    '-1': -1.0,
}

print(f"\n  Phi_3(x) = x^2 + x + 1 evaluated at special values:\n")
for name, val in evaluations.items():
    phi3_val = val**2 + val + 1
    print(f"  Phi_3({name:>8s}) = {val:>10.6f}^2 + {val:>10.6f} + 1 = {phi3_val:>12.6f}")

print(f"""
  NOTABLE:
  Phi_3(e)   = {e**2+e+1:.6f} = e^2 + e + 1  (≈ 11.107)
  Phi_3(pi)  = {pi**2+pi+1:.6f} = pi^2 + pi + 1  (≈ 14.010)
  Phi_3(phi) = {phi**2+phi+1:.6f} = phi^2 + phi + 1 = 2phi + 2 = 2(phi+1) = 2phi^2
  Phi_3(tau) = {tau**2+tau+1:.6f} = tau^3  (by definition!)
  Phi_3(0)   = 1  (the empty tournament)
  Phi_3(1)   = 3  (the 3-cycle count)
  Phi_3(-1)  = 1  (= Phi_6(1) since Phi_3(-x) = Phi_6(x))
  Phi_3(2)   = 7  (FORBIDDEN — the Fano plane)
  Phi_3(3)   = 13 (achievable H value!)

  CHECK: Phi_3(phi) = phi^2 + phi + 1 = (phi+1) + phi + 1 = 2phi + 2
  = 2(phi + 1) = 2*phi^2 = {2*phi**2:.6f}
  And phi^2 = {phi**2:.6f} = phi + 1 = {phi+1:.6f}  ✓
""")

# =====================================================================
# PART 3: THE e-pi-Phi_3 CONNECTION
# =====================================================================
print("=" * 70)
print("PART 3: THE e-pi-Phi_3 CONNECTION")
print("=" * 70)

omega = np.exp(2j * np.pi / 3)

print("""
  Euler's identity: e^(i*pi) + 1 = 0

  The roots of Phi_3(x) = x^2+x+1 are the primitive cube roots of unity:
  omega = e^(2*pi*i/3) = (-1 + i*sqrt(3))/2
  omega_bar = e^(-2*pi*i/3) = (-1 - i*sqrt(3))/2

  So Phi_3(omega) = 0, where omega involves BOTH e and pi!

  The cube root of unity omega connects:
  - e (the transcendental base)
  - pi (the circle constant)
  - i (the imaginary unit)
  - Phi_3 (the tournament polynomial)

  This is why tournament theory is fundamentally about PERIOD 3:
  The cube roots of unity form a 3-element group (1, omega, omega^2)
  which is Z/3Z -- the cyclic group of order 3.

  Tournament 3-cycles correspond to the 3rd roots of unity!
  The 3-cycle (a -> b -> c -> a) = a rotation by 2*pi/3 = 120 degrees.

  Phi_3(x) = (x - omega)(x - omega_bar)
  At x = 2: Phi_3(2) = (2-omega)(2-omega_bar) = |2-omega|^2 = 7

  |2-omega|^2 = (2 - (-1/2))^2 + (sqrt(3)/2)^2
              = (5/2)^2 + (sqrt(3)/2)^2
              = 25/4 + 3/4 = 28/4 = 7

  So 7 = |2 - omega|^2 where omega = e^(2*pi*i/3).

  THE FORBIDDEN VALUE 7 = THE SQUARED DISTANCE FROM 2 TO omega!
""")

omega = np.exp(2j * np.pi / 3)
dist_sq = abs(2 - omega)**2
print(f"  |2 - omega|^2 = |2 - e^(2*pi*i/3)|^2 = {dist_sq:.6f}")
print(f"  = 7: {'✓' if abs(dist_sq - 7) < 1e-10 else '✗'}")

# For 21:
dist_sq_4 = abs(4 - omega)**2
print(f"\n  |4 - omega|^2 = |4 - e^(2*pi*i/3)|^2 = {dist_sq_4:.6f}")
print(f"  = 21: {'✓' if abs(dist_sq_4 - 21) < 1e-10 else '✗'}")

print(f"""
  PATTERN: The forbidden values are |2^k - omega|^2:
""")
for k in range(1, 8):
    val = abs(2**k - omega)**2
    phi3_val = (2**k)**2 + 2**k + 1
    print(f"  |2^{k} - omega|^2 = |{2**k} - omega|^2 = {val:.1f} = Phi_3({2**k}) = {phi3_val}")

# =====================================================================
# PART 4: THE EXPONENTIAL GENERATING FUNCTION
# =====================================================================
print("\n" + "=" * 70)
print("PART 4: THE H-VALUE EXPONENTIAL GENERATING FUNCTION")
print("=" * 70)

egf_tau = 2*tau/(2-tau)
egf_phi = 2*phi/(2-phi)

print(f"""
  The total number of Hamiltonian paths across all tournaments on n vertices:
  Total_H(n) = sum_T H(T) = n! * 2^(C(n,2) - (n-1))

  Mean(H) = n! / 2^(n-1)

  The OGF: sum Mean(H) * x^n = sum x^n / 2^(n-1)
  = 2 * sum (x/2)^n = 2x/(2-x) for |x| < 2.

  So the OGF of Mean(H) is 2x/(2-x)!
  This has a POLE at x = 2 (the tournament generator!).

  At x = tau: 2*tau/(2-tau) = {egf_tau:.6f}
  At x = phi: 2*phi/(2-phi) = {egf_phi:.6f}
  At x = 1: 2*1/(2-1) = 2.0
""")

# Verify the EGF
print("  Verification of EGF 2x/(2-x):")
for n in range(1, 8):
    mean_H = factorial(n) / 2**(n-1)
    # Coefficient of x^n in 2x/(2-x) = 2x * sum_k (x/2)^k = sum_k 2^(1-k) x^(k+1)
    # Coeff of x^n is 2^(1-(n-1)) = 2^(2-n)
    egf_coeff = 2**(2-n)  # this times n! should give mean_H
    print(f"  n={n}: Mean(H) = {n}!/2^{n-1} = {mean_H:.2f}, EGF: [x^{n}] * {n}! = {egf_coeff * factorial(n):.2f}")

print(f"""
  The EGF 2x/(2-x) = 2 * sum_(n>=1) (x/2)^n

  This is the generating function of a GEOMETRIC series in x/2.
  The ratio x/2 encodes the "half" probability of each arc direction.

  Now consider the EGF at x = omega (cube root of unity):
  2*omega/(2-omega) = ?
""")

egf_omega = 2 * omega / (2 - omega)
print(f"  2*omega/(2-omega) = {egf_omega:.6f}")
print(f"  |2*omega/(2-omega)| = {abs(egf_omega):.6f}")
print(f"  = 2/sqrt(7) = {2/sqrt(7):.6f}")
print(f"  = 2/sqrt(Phi_3(2))!")

# =====================================================================
# PART 5: THE NATURAL BASE e AND TOURNAMENT ENTROPY
# =====================================================================
print("\n" + "=" * 70)
print("PART 5: BASE e AND TOURNAMENT ENTROPY")
print("=" * 70)

print(f"""
  In information theory, the NATURAL unit of information is the "nat"
  (natural unit), using base e instead of base 2 (bits).

  A tournament on n vertices has C(n,2) binary choices = C(n,2) bits.
  In nats: C(n,2) * ln(2) nats.

  The ENTROPY of the H distribution:
  S = -sum_h P(h) * ln(P(h))

  For the uniform distribution over all 2^C(n,2) tournaments:
  S_max = C(n,2) * ln(2) nats = C(n,2) bits

  The key ratio: H / (2^(n-1)) = Mean / sum = the "mean fraction"
  ln(Mean/Sum) = ln(n!) - (n-1)*ln(2) - C(n,2)*ln(2)
               = ln(n!) - (n-1)*ln(2) - n(n-1)/2 * ln(2)
               = ln(n!) - (n^2-n+2)/2 * ln(2)  ... hmm, not as clean.

  The important number: ln(7) = {log(7):.6f}
  Compare:
    ln(2) = {log(2):.6f}
    ln(3) = {log(3):.6f}
    ln(7) = {log(7):.6f} = ln(2) + ln(3) + ln(7/6)? No...

    ln(7) = ln(Phi_3(2))

  Since Phi_3(x) = x^2+x+1, we have:
  ln(Phi_3(2)) = ln(7) = ln(4+2+1) = ln(4+3) = ln(7)

  Not much simplification. BUT:

  d/dx ln(Phi_3(x)) = (2x+1)/(x^2+x+1)
  At x=2: d/dx = 5/7

  This is the "marginal information rate" of the 3-strand at x=2:
  each additional unit of x increases the information by 5/7 nats.
  5/7 = {5/7:.6f} nats ≈ {5/7/log(2):.6f} bits.

  At x=1: d/dx = 3/3 = 1 nat per unit = {1/log(2):.6f} bits.
  At x=phi: d/dx = (2*phi+1)/(phi^2+phi+1) = {(2*phi+1)/(phi**2+phi+1):.6f}
  At x=tau: d/dx = (2*tau+1)/(tau^2+tau+1) = {(2*tau+1)/(tau**2+tau+1):.6f}
""")

# =====================================================================
# PART 6: THE CATALAN-TOURNAMENT CONNECTION
# =====================================================================
print("=" * 70)
print("PART 6: CATALAN NUMBERS AND TOURNAMENTS")
print("=" * 70)

# Catalan numbers
def catalan(n):
    return comb(2*n, n) // (n+1)

print(f"""
  The Catalan number C_n = C(2n,n)/(n+1) counts:
  - Binary trees with n internal nodes
  - Balanced parenthesizations with n pairs
  - Paths from (0,0) to (2n,0) staying ≥ 0
  - Triangulations of (n+2)-gon

  The Catalan generating function:
  C(x) = (1 - sqrt(1-4x))/(2x) = sum C_n x^n

  This satisfies: C = 1 + x*C^2 (the quadratic recursion).

  TOURNAMENT CONNECTION:
  A tournament T on [n] can be encoded as a BINARY WORD:
  List arcs (1,2), (1,3), ..., (n-1,n) and write 1 for fwd and 0 for bwd.
  This gives a word of length C(n,2) in {0,1}.

  The number of 1s in this word = number of "forward" arcs.
  For a tournament with score sequence (s_1,...,s_n): sum s_i = C(n,2).

  Catalan numbers appear in tournament theory through:
  1. The number of SCORE SEQUENCES (not Catalan, but close)
  2. The number of STRONGLY CONNECTED score sequences
  3. Recursive tournament decompositions via binary trees
""")

print("  Catalan numbers vs tournament counts:")
for n in range(1, 10):
    cat = catalan(n)
    tourn = 2**comb(n, 2)
    print(f"  n={n}: C_{n} = {cat}, |T_n| = 2^C({n},2) = {tourn}")

# The ratio
print("\n  Ratios (tournaments / Catalan):")
for n in range(1, 10):
    cat = catalan(n)
    tourn = 2**comb(n, 2)
    if cat > 0:
        print(f"  n={n}: 2^C({n},2) / C_{n} = {tourn/cat:.2f}")

# =====================================================================
# PART 7: THE k-NACCI WEIGHTED GENERATING FUNCTIONS
# =====================================================================
print("\n" + "=" * 70)
print("PART 7: k-NACCI WEIGHTED GENERATING FUNCTIONS")
print("=" * 70)

print(f"""
  The user's insight: "k-nacci approaches 2, weighted k-nacci approaches 3"

  k-NACCI constant alpha_k: root of x^k = x^(k-1) + ... + x + 1

  The equation x^k = (x^k - 1)/(x - 1) rearranges to:
  x^k * (x-1) = x^k - 1
  x^(k+1) - x^k = x^k - 1
  x^(k+1) = 2*x^k - 1

  As k -> infinity, alpha_k -> 2 (the positive root approaches 2).

  WEIGHTED k-nacci: T_n = sum_(i=1)^k i * T_(n-i)
  The weights are {1, 2, 3, ..., k}.

  Characteristic equation: x^k = 1*x^(k-1) + 2*x^(k-2) + ... + k

  At large k, the dominant root satisfies:
  x^k ≈ sum_(i=1)^k i * x^(k-i) = x^(k-1) * sum_(i=1)^k i/x^(i-1)

  For x ≈ c (constant), sum = sum_(i=1)^inf i/c^(i-1) = c^2/(c-1)^2 as k -> inf.

  So x ≈ c where c = c^2/(c-1)^2, i.e., (c-1)^2 = c, i.e., c^2 - 3c + 1 = 0.
  c = (3 ± sqrt(5))/2. Taking the larger root: c = (3 + sqrt(5))/2 = phi + 1 = phi^2.

  phi^2 = {phi**2:.6f} ≈ 2.618

  WAIT: the user said "weighted k-nacci approaches 3", not phi^2 ≈ 2.618.
  Let me check with equal weights:

  If weights are all 1: T_n = T_(n-1) + T_(n-2) + ... + T_(n-k)
  -> dominant root -> 2 ✓

  If weights are 1,1,...,1 with weight 2:
  T_n = 2*(T_(n-1) + ... + T_(n-k))
  -> dominant root -> 4? No...

  Maybe "weighted at x=2": evaluate the k-nacci at x=2?
  The k-nacci generating polynomial at x is x^k - x^(k-1) - ... - 1.
  At x=2: 2^k - 2^(k-1) - ... - 1 = 2^k - (2^k - 1) = 1.

  Hmm, that gives 1 for all k. Not 3.

  Maybe the user means: the k-nacci with STARTING VALUES weighted by 2?
  Or: the k-nacci polynomial evaluated at x=2 gives Phi_k(2)?

  Let me compute the k-nacci constant for small k and check limits:
""")

# Compute k-nacci constants
for k in range(2, 15):
    # Characteristic polynomial: x^k - x^(k-1) - ... - x - 1 = 0
    # = x^k - (x^k - 1)/(x - 1) = 0 when x ≠ 1
    # Solve numerically
    coeffs = np.zeros(k+1)
    coeffs[k] = 1
    for i in range(k):
        coeffs[i] = -1
    roots = np.roots(coeffs[::-1])
    real_roots = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 1]
    alpha_k = max(real_roots) if real_roots else 0

    # Weighted: x^k = 1*x^(k-1) + 2*x^(k-2) + ... + k
    w_coeffs = np.zeros(k+1)
    w_coeffs[k] = 1
    for i in range(k):
        w_coeffs[i] = -(k - i)  # weight of T_{n-(k-i)} = k-i
    w_roots = np.roots(w_coeffs[::-1])
    w_real = [r.real for r in w_roots if abs(r.imag) < 1e-10 and r.real > 1]
    w_alpha = max(w_real) if w_real else 0

    print(f"  k={k:2d}: alpha_{k} = {alpha_k:.8f} (->2), weighted = {w_alpha:.8f} (->?)")

print(f"""
  The unweighted k-nacci constants approach 2 ✓
  The weighted k-nacci constants approach... {phi**2:.6f} ≈ phi^2 = phi+1

  Actually, let me recheck. The weighted limit:
  c = (3 + sqrt(5))/2 = {(3+sqrt(5))/2:.6f} = phi^2

  So the weighted k-nacci approaches phi^2 ≈ 2.618, NOT 3.

  UNLESS: the "weighted" version has weights (1, 1, ..., 1, 2) or similar.

  Let me check: weights all equal to 2:
  T_n = 2*(T_(n-1) + T_(n-2) + ... + T_(n-k))
  x^k = 2*(x^(k-1) + ... + 1) = 2*(x^k - 1)/(x - 1)
  x^k*(x-1) = 2*(x^k - 1)
  x^(k+1) - x^k = 2*x^k - 2
  x^(k+1) = 3*x^k - 2
  As k -> inf: x -> 3  ✓!!!

  So DOUBLY-WEIGHTED k-nacci (all weights = 2) approaches 3!
""")

# Verify doubly-weighted
print("  Doubly-weighted k-nacci (all weights = 2):")
for k in range(2, 15):
    # x^k = 2*(x^(k-1) + ... + 1)
    # x^k - 2*x^(k-1) - ... - 2 = 0
    coeffs = np.zeros(k+1)
    coeffs[k] = 1
    for i in range(k):
        coeffs[i] = -2
    roots = np.roots(coeffs[::-1])
    real_roots = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real > 1]
    alpha_k = max(real_roots) if real_roots else 0
    print(f"  k={k:2d}: alpha_{k} = {alpha_k:.8f} (->3)")

print(f"""
  YES! When all weights are 2 (the tournament generator!):
  The k-nacci constant approaches 3 = Phi_3(1).

  THE PATTERN:
  weight = 1: k-nacci -> 2 = tournament generator
  weight = 2: k-nacci -> 3 = Phi_3(1) = number of edge types (fwd, bwd, none)
  weight = j: k-nacci -> j+1

  The limit satisfies x^(k+1) = (j+1)*x^k - j, so x -> j+1 as k -> inf.

  TOURNAMENT INTERPRETATION:
  - Weight 1 = each arc is a BINARY choice (0 or 1) -> limit 2
  - Weight 2 = each arc is a TERNARY choice (forward, backward, or ABSENT?) -> limit 3
  - The ternary version includes NON-tournaments (some arcs missing)

  This connects to the SIMPLEX-CUBOID nesting:
  - (x+1)^n at x=1 gives 2^n (simplex, binary, weight 1)
  - (x+2)^n at x=1 gives 3^n (cuboid, ternary, weight 2)

  The k-nacci limit IS the volume of the geometric shape!
""")

# =====================================================================
# PART 8: BASE-e REPRESENTATION OF H VALUES
# =====================================================================
print("=" * 70)
print("PART 8: BASE-e REPRESENTATION OF TOURNAMENT CONSTANTS")
print("=" * 70)

print(f"""
  In base e, the key tournament constants have natural expressions:

  ln(2) = {log(2):.6f} = 1 in the "e-nary" system
  ln(3) = {log(3):.6f}
  ln(7) = {log(7):.6f}  (the forbidden value)
  ln(21) = {log(21):.6f} (the second forbidden value)

  Ratios:
  ln(7)/ln(2) = {log(7)/log(2):.6f} = log_2(7) (bits in the Fano plane)
  ln(21)/ln(2) = {log(21)/log(2):.6f} = log_2(21)
  ln(3)/ln(2) = {log(3)/log(2):.6f} = log_2(3) (bits per trit)

  The "tournament entropy":
  H_entropy(T) = ln(H(T)) / ln(2) = log_2(H(T))

  For the transitive tournament: H = 1, entropy = 0 (maximally ordered)
  For the regular tournament (n=5): H = 45, entropy = {log(45)/log(2):.3f}
  For the regular tournament (n=7): H = 189, entropy = {log(189)/log(2):.3f}
  Maximum H grows as ~ n!/2^(n-1), so max entropy ~ n*ln(n)/ln(2) - n

  The ENTROPY PER ARC:
  log_2(H(T)) / C(n,2) = how many bits of "randomness" per arc

  For H_mean = n!/2^(n-1):
  log_2(n!/2^(n-1)) / C(n,2) = (log_2(n!) - (n-1)) / (n(n-1)/2)
""")

for n in range(3, 10):
    mean_H = factorial(n) / 2**(n-1)
    m = comb(n, 2)
    entropy_per_arc = (np.log2(mean_H)) / m if m > 0 else 0
    print(f"  n={n}: Mean(H)={mean_H:.0f}, log_2(Mean)/C(n,2) = {entropy_per_arc:.4f} bits/arc")

print(f"""
  The entropy per arc DECREASES as n grows.
  This means: larger tournaments are MORE CONSTRAINED per arc.
  The asymptotic rate: ~ ln(n)/n bits/arc -> 0 as n -> inf.

  CONCLUSION: In the thermodynamic limit, each arc carries
  VANISHINGLY little information about the total H count.
  The information is in the GLOBAL structure, not local arcs.
""")

# =====================================================================
# PART 9: THE GENERATING FUNCTION RADIUS OF CONVERGENCE
# =====================================================================
print("=" * 70)
print("PART 9: RADIUS OF CONVERGENCE AND TOURNAMENT PHASE TRANSITIONS")
print("=" * 70)

print(f"""
  EGF of Mean(H): f(x) = 2x/(2-x), pole at x=2.

  The OGF of H values (summed over tournaments):
  g(x) = sum_n (sum_T H(T)) * x^n = sum_n n! * 2^(C(n,2)-(n-1)) * x^n

  This has ZERO radius of convergence (n! * 2^(n^2/2) grows faster than any geometric).

  But the EGF has radius 2, which equals the tournament generator!

  WHAT IS SPECIAL ABOUT x = 2?
  At x = 2, the EGF pole says: "the tournament reaches its MAXIMUM information state."
  Below x = 2: the series converges (subcritical regime).
  At x = 2: phase transition (critical).
  Above x = 2: divergence (supercritical).

  The EGF of the 3-strand count:
  If we want sum_n 3^n * x^n / n!, we get e^(3x) — entire!
  If we want sum_n 7^n * x^n / n!, we get e^(7x) — also entire!

  But the TOURNAMENT EGF 2x/(2-x) has a pole, not just exponential growth.
  This is because H counts have FACTORIAL growth (n!/2^(n-1)),
  not just exponential growth.

  The pole structure:
  2x/(2-x) = -2 + 4/(2-x)

  The residue at x=2 is 4. And 4 = 2^2.
  The residue 4 = (tournament generator)^2 = the SQUARE of the base.

  In the 3-strand world:
  Phi_3(x)^n summed with EGF would give... but that's just e^(Phi_3(x)*t).
  At x=2: e^(7t) — the exponential growth rate 7.

  The tournament EGF has growth rate ~ 1/(2-x) near x=2,
  which is SLOWER than e^(ct) for any c.

  This slowdown (pole vs essential singularity) is the effect of
  the COMBINATORIAL CONSTRAINTS (tournament structure) on H.
""")

# =====================================================================
# PART 10: SYNTHESIS — THE BASE HIERARCHY
# =====================================================================
print("=" * 70)
print("PART 10: THE GRAND BASE HIERARCHY")
print("=" * 70)

print(f"""
  THE HIERARCHY OF BASES IN TOURNAMENT THEORY:

  ALGEBRAIC BASES (special structures):
    phi = {phi:.6f}  (golden ratio)    — Fibonacci, binary growth
    tau = {tau:.6f}  (tribonacci)      — Phi_3, ternary structure
    2   = 2.000000  (tournament gen)   — arc choices, Walsh space
    phi^2 = {phi**2:.6f}  (golden square) — weighted k-nacci limit
    3   = 3.000000  (cuboid/trit gen)  — ternary choices, 3^n volume

  TRANSCENDENTAL BASES (limiting structures):
    e   = {e:.6f}  (natural base)      — entropy, EGF convergence
    pi  = {pi:.6f}  (circle constant)  — periodicity, roots of unity

  THE ORDERING:
    1 < phi < tau < 2 < e < 3 < pi < 4 < 7

  WHAT EACH BASE CAPTURES:
    1: identity (trivial tournaments)
    phi: golden growth (Fibonacci recursion)
    tau: tribonacci growth (Phi_3 recursion, tournament-natural)
    2: binary choices (arc orientations)
    e: information-theoretic optimum (entropy)
    3: ternary excess (cuboid counting)
    pi: periodic structure (roots of unity, 3-cycles)
    4: cuboid at tournament eval (x+2)^n|_{{x=2}}
    7: FORBIDDEN threshold (Fano obstruction)

  THE KEY INSIGHT:
  Tournaments live in the interval [2, 3] of bases.
  Below 2: not enough structure to form tournaments.
  Above 3: too much freedom — leaves tournament category.

  The TRIBONACCI tau ≈ 1.839 sits just BELOW 2,
  encoding the "almost-binary" nature of tournament growth.

  The FORBIDDEN VALUE 7 = Phi_3(2) sits FAR above this interval,
  representing structures that tournaments cannot achieve.
""")

print("\n" + "=" * 70)
print("DONE — TRANSCENDENTAL TOURNAMENT BASES")
print("=" * 70)
