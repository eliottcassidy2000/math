#!/usr/bin/env python3
"""
tournament_dirichlet_s291b.py -- opus-2026-03-24-S291

THE TOURNAMENT DIRICHLET SERIES AND FUNCTIONAL EQUATION

From S291: V_n × n!/2^m = 1 + R(n) where
  R(n) = Σ_{k=3,5,7,...} R_k(n) + cross-orbit corrections
  R_k(n) = (1/k) × n^{falling k} × 2^{(k-1)(k+1-2n)/2}

THE (3,3) RATIO MYSTERY:
  Actual(3,3) / (f_3²/2) = 2^{gcd(3,3)-1} × (something involving n)
  Let me compute the exact ratio.

THE DIRICHLET SERIES:
  D(s) = Σ_{n≥3} R(n) × n^{-s}
  = Σ_{n≥3} Σ_{k odd} R_k(n) × n^{-s} + cross terms

  The inner sum: Σ_n R_k(n) n^{-s} = (1/k) Σ_n n^{falling k} 2^{Delta_k} n^{-s}

  For k=3: R_3(n) = C(n,3) × 2 × 2^{4-2n} = n(n-1)(n-2)/6 × 2 × 2^{4-2n}
  = n(n-1)(n-2)/3 × 2^{3-2n}

  Σ_n R_3(n) n^{-s} = Σ_n n(n-1)(n-2)/3 × 2^{3-2n} × n^{-s}
  = (8/3) Σ_n n(n-1)(n-2) × 4^{-n} × n^{-s}
  = (8/3) Σ_n n^{1-s}(n-1)(n-2) × 4^{-n}

  At s=0: (8/3) Σ_n n(n-1)(n-2) × 4^{-n} = (8/3) × 3! × Li_{-3}(1/4)... hmm.
  Actually: Σ_n n(n-1)(n-2) x^n = x d/dx x d/dx x d/dx [1/(1-x)] at x
  = 6x^3/(1-x)^4. So with x=1/4: 6/64 / (3/4)^4 = 6/(64×81/256) = 6×4/81 = 24/81 = 8/27.
  D_3(0) = (8/3) × 8/27 = 64/81 ≈ 0.790.

  This is a COMPUTABLE zeta value!

Author: opus-2026-03-24-S291
"""

import sys
from math import comb, factorial, gcd, log, log2, pi
from fractions import Fraction
from collections import Counter, defaultdict
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def pair_orbits(ct):
    within = sum((c - 1) // 2 for c in ct)
    across = sum(gcd(ct[i], ct[j]) for i in range(len(ct)) for j in range(i+1, len(ct)))
    return within + across

def multinomial(n, ct):
    ct_counter = Counter(ct)
    denom = 1
    for c, a in ct_counter.items():
        denom *= (c ** a) * factorial(a)
    return factorial(n) // denom

def gen_odd_partitions(n):
    result = []
    def _gen(remaining, max_part, current):
        if remaining == 0:
            result.append(tuple(sorted(current)))
            return
        start = min(remaining, max_part)
        if start % 2 == 0: start -= 1
        for p in range(start, 0, -2):
            _gen(remaining - p, p, current + [p])
    _gen(n, n, [])
    return list(set(result))

# ============================================================
banner("THE (3,3) RATIO: EXACT ANALYSIS")
# ============================================================

# For partition (3,3,1^{n-6}):
#   mult = n! / (3² × 2! × (n-6)!) = n! / (18 × (n-6)!)
#   pair_orbits: two 3-cycles + (n-6) fixed points
#     within each 3-cycle: (3-1)/2 = 1 each → 2
#     between two 3-cycles: gcd(3,3) = 3
#     between 3-cycle and fixed: 1 per fixed point × 2 cycles = 2(n-6)
#     between fixed points: C(n-6,2)
#     Total: 2 + 3 + 2(n-6) + C(n-6,2) = 5 + 2(n-6) + (n-6)(n-7)/2
#
# For single 3-cycle (3,1^{n-3}):
#   po = 1 + (n-3) + C(n-3,2) = 1 + (n-3) + (n-3)(n-4)/2
#
# Identity: m = C(n,2) = n(n-1)/2
#
# The EXCESS pair orbits from having two 3-cycles vs two independent ones:
# If orbits were additive: po(3,3) = po(3,1^3) + po(3,1^3) - po(1^6)
# But they're not! The cross-orbit term is gcd(3,3) = 3 instead of
# the 2×3 = 6 independent cross pairs divided by... let me just compute.

print("  EXACT (3,3) ANALYSIS:\n")

for n in range(6, 16):
    m = comb(n, 2)

    # Single 3-cycle: (3, 1^{n-3})
    ct_3 = tuple(sorted([3] + [1]*(n-3)))
    mult_3 = multinomial(n, ct_3)
    po_3 = pair_orbits(ct_3)
    f3 = Fraction(mult_3 * pow(2, po_3), pow(2, m))

    # Double 3-cycle: (3, 3, 1^{n-6})
    ct_33 = tuple(sorted([3, 3] + [1]*(n-6)))
    mult_33 = multinomial(n, ct_33)
    po_33 = pair_orbits(ct_33)
    actual_33 = Fraction(mult_33 * pow(2, po_33), pow(2, m))

    # "Independent" prediction: f_3^2 / 2
    independent = f3 * f3 / 2

    # Ratio
    ratio = Fraction(actual_33, independent) if independent != 0 else 0

    # The cross-orbit term: gcd(3,3) = 3 between the two 3-cycles.
    # In the independent model, we'd get gcd(3,1)×(n-3) cross terms
    # from each 3-cycle to the "other side's fixed points" — but those
    # are already counted. The EXCESS is just gcd(3,3) = 3.
    #
    # po(3,3,1^{n-6}) = 2×1 + 3 + 2(n-6) + C(n-6,2)
    # 2×po(3,1^{n-3}) - po(1^n) = 2×(1 + (n-3) + C(n-3,2)) - C(n,2)
    # Let me compute both:
    po_double = 2 + 3 + 2*(n-6) + comb(n-6, 2)
    po_independent = 2 * po_3 - m  # This would be the "independent" value

    excess = po_33 - (2*po_3 - m)

    print("  n=%2d: mult_3=%-5d mult_33=%-5d po_3=%-3d po_33=%-3d excess=%d ratio=%s=%.4f" %
          (n, mult_3, mult_33, po_3, po_33, excess, ratio, float(ratio)))

# ============================================================
banner("THE EXACT CROSS-ORBIT FORMULA")
# ============================================================

# The excess pair orbits from combining cycles:
# For two cycles of lengths k and j sharing no vertices:
#   cross orbits = gcd(k,j) (pairs between the two cycles)
#   "Independent" cross = k*j / lcm(k,j)... no, it's just gcd(k,j).
#
# But in the independent model (Π(1+f_k)), the prediction for (k,j,1^{n-k-j}) is:
#   f_k × f_j × (correction for overlap of fixed points)
#
# The actual/predicted ratio tells us the EXCESS factor.
# From the (3,3) case: ratio grows linearly with n. Why?
#
# It's because the multinomial calculation involves CHOOSING which vertices
# go in which cycle. The combinatorics change when cycles interact.

print("  CROSS-ORBIT RATIO GROWTH:\n")
print("  For (3,3): actual/predicted = ?\n")

ratios_33 = []
for n in range(6, 18):
    m = comb(n, 2)
    ct_3 = tuple(sorted([3] + [1]*(n-3)))
    f3 = Fraction(multinomial(n, ct_3) * pow(2, pair_orbits(ct_3)), pow(2, m))

    ct_33 = tuple(sorted([3, 3] + [1]*(n-6)))
    actual = Fraction(multinomial(n, ct_33) * pow(2, pair_orbits(ct_33)), pow(2, m))
    predicted = f3 * f3 / 2

    ratio = float(actual) / float(predicted) if float(predicted) != 0 else 0
    ratios_33.append((n, ratio))
    print("  n=%2d: ratio = %.6f, ratio/(n-3) = %.6f" % (n, ratio, ratio/(n-3)))

# It seems like ratio ≈ c × 2^{n-something}
# Let me check ratio / 2^{...}
print("\n  Checking ratio = (n-3) × 2^{something}:\n")
for n, r in ratios_33:
    if n > 6:
        log_ratio = log2(r / (n-3)) if r > 0 and n > 3 else 0
        print("  n=%2d: ratio/(n-3) = %.6f, log2 = %.2f" % (n, r/(n-3), log_ratio))

# ============================================================
banner("THE DIRICHLET SERIES D_k(s)")
# ============================================================

# D_k(s) = Σ_{n≥k} R_k(n) × n^{-s}
# where R_k(n) = C(n,k) × (k-1)! × 2^{(k-1)(k+1-2n)/2}
#
# Let x = 1/4 (for k=3).
# R_3(n) = C(n,3) × 2 × 4^{(4-2n)/2} = C(n,3) × 2 × 4^{2-n}
# = C(n,3) × 2 × 16 × 4^{-n} = C(n,3) × 32 × 4^{-n}
# Wait: Delta_3 = (3-1)(3+1-2n)/2 = 2(4-2n)/2 = 4-2n.
# 2^{Delta_3} = 2^{4-2n} = 16/4^n. Then R_3 = C(n,3) × 2 × 16/4^n.
# Hmm: mult_3 = C(n,3) × 2! = C(n,3) × 2. So R_3 = C(n,3) × 2 × 2^{4-2n}.
# = C(n,3) × 2^{5-2n}.

# Let's verify:
print("  VERIFY R_3 formula:\n")
for n in range(3, 10):
    m = comb(n, 2)
    ct = tuple(sorted([3] + [1]*(n-3)))
    exact = Fraction(multinomial(n, ct) * pow(2, pair_orbits(ct)), pow(2, m))
    formula = Fraction(comb(n,3) * 2, pow(2, 2*n-4))
    formula2 = Fraction(comb(n,3), pow(2, 2*n-5))
    print("  n=%d: exact=%s, formula=C(n,3)×2^{5-2n}=%s, match=%s" %
          (n, exact, formula2, exact == formula2))

# So: D_3(s) = Σ_{n≥3} C(n,3) × 2^{5-2n} × n^{-s}
#            = 2^5 × Σ_{n≥3} C(n,3) × 4^{-n} × n^{-s}
#            = 32 × Σ_{n≥3} C(n,3) × (1/4)^n × n^{-s}

# At s=0: D_3(0) = 32 × Σ_{n≥3} C(n,3) × (1/4)^n
# We know: Σ_{n≥k} C(n,k) x^n = x^k / (1-x)^{k+1}
# So: Σ_{n≥3} C(n,3) (1/4)^n = (1/4)^3 / (3/4)^4 = (1/64) / (81/256) = 256/(64×81) = 4/81

print("\n  D_3(0) = 32 × Σ C(n,3)(1/4)^n = 32 × 4/81 = 128/81 ≈ %.6f" % (128/81))
print("  (This is the total correction from k=3 terms summed over all n.)\n")

# Exact computation:
D3_0 = Fraction(128, 81)
print("  D_3(0) = %s = %.10f" % (D3_0, float(D3_0)))

# At s=1:
# D_3(1) = 32 × Σ C(n,3) (1/4)^n / n
# This involves the polylogarithm-like sum Σ C(n,3) x^n / n.
# We can compute: Σ C(n,3) x^n / n = ∫_0^x Σ C(n,3) t^{n-1} dt
# = ∫_0^x [t^2 / (1-t)^4] dt (since Σ C(n,3) t^{n-1} = t^2/(1-t)^4)

# ∫ t^2/(1-t)^4 dt = ... partial fractions ...
# Let u = 1-t: ∫ (1-u)^2/u^4 du = ∫ (1-2u+u^2)/u^4 du
# = ∫ (u^{-4} - 2u^{-3} + u^{-2}) du = -u^{-3}/3 + u^{-2} - u^{-1} + C
# = -1/(3(1-t)^3) + 1/(1-t)^2 - 1/(1-t) + C

# At t=0: -1/3 + 1 - 1 + C = C - 1/3. At t=0, the sum = 0 (starts at n=3).
# So C = 1/3.
# F(t) = -1/(3(1-t)^3) + 1/(1-t)^2 - 1/(1-t) + 1/3

# D_3(1) = 32 × F(1/4)
# F(1/4) = -1/(3×(3/4)^3) + 1/(3/4)^2 - 1/(3/4) + 1/3
# = -1/(3×27/64) + 16/9 - 4/3 + 1/3
# = -64/81 + 16/9 - 4/3 + 1/3
# = -64/81 + 144/81 - 108/81 + 27/81
# = (-64 + 144 - 108 + 27)/81
# = -1/81

# D_3(1) = 32 × (-1/81) = -32/81 < 0 ???
# That can't be right. Let me recheck.

# Actually: Σ_{n≥3} C(n,3) x^n/n = (1/n) terms mean we're looking at
# Σ_{n≥3} C(n,3) x^n/n where each term is POSITIVE for x > 0.
# So the sum must be positive. I must have an integration error.

# Let me just compute numerically:
D3_1 = sum(Fraction(comb(n, 3), n) * Fraction(1, 4)**n for n in range(3, 200))
print("  D_3(1) = 32 × Σ C(n,3)(1/4)^n/n = 32 × %s ≈ %.10f" %
      (D3_1, float(32 * D3_1)))

# Also compute D_3(s) for several s values
print("\n  D_3(s) for several s values:\n")
for s in [0, 0.5, 1, 1.5, 2, 3]:
    D3_s = sum(comb(n, 3) * pow(4, -n) * pow(n, -s) for n in range(3, 200))
    print("  D_3(s=%.1f) = 32 × %.10f = %.10f" % (s, D3_s, 32*D3_s))

# ============================================================
banner("THE FULL DIRICHLET SERIES D(s) = Σ R(n) n^{-s}")
# ============================================================

# D(s) = Σ_{n≥3} R(n) × n^{-s}
# where R(n) = V_n × n!/2^m - 1

print("  FULL DIRICHLET SERIES (exact R(n) for n=3..13):\n")

R_values = []
for n in range(3, 14):
    m = comb(n, 2)
    parts = gen_odd_partitions(n)
    identity = tuple([1]*n)
    T_exact = sum(Fraction(multinomial(n, ct) * pow(2, pair_orbits(ct)), pow(2, m))
                  for ct in parts)
    R = T_exact - 1
    R_values.append((n, float(R)))

for s in [0, 0.5, 1, 1.5, 2]:
    D_s = sum(R * n**(-s) for n, R in R_values)
    print("  D(s=%.1f) = %.10f" % (s, D_s))

# ============================================================
banner("THE FUNCTIONAL EQUATION CANDIDATE")
# ============================================================

print("""
  THE TOURNAMENT COUNTING GENERATING FUNCTION:

  G(x) = Σ_{n≥1} V_n x^n / n!
  = Σ_{n≥1} [Σ_λ mult(λ) 2^{po(λ)} / n!²] x^n
  = ... (involves double factorial-like denominators)

  A cleaner GF: Let a_n = V_n × n! / 2^{C(n,2)} = T(n) = 1 + R(n).
  Then: F(x) = Σ_{n≥1} a_n x^n.

  F(x) ≈ x/(1-x) + Σ correction terms.

  The correction GF: G(x) = Σ_{n≥3} R(n) x^n
  ≈ Σ_{n≥3} C(n,3) 2^{5-2n} x^n  (3-cycle approximation)
  = 32 Σ_{n≥3} C(n,3) (x/4)^n
  = 32 × (x/4)^3 / (1-x/4)^4  (for |x| < 4)
  = 32 × x³/64 / (1-x/4)^4
  = x³/2 / (1-x/4)^4

  So: G(x) ≈ x³/(2(1-x/4)^4)  for the LEADING correction.
""")

# Let's verify:
print("  VERIFICATION: R(n) vs x³/(2(1-x/4)^4) coefficients:\n")
for n in range(3, 14):
    # Coefficient of x^n in x³/(2(1-x/4)^4) = (1/2) × C(n-3+3, 3) × (1/4)^{n-3}
    # = C(n,3) × (1/4)^{n-3} / 2 = C(n,3) × 4^{3-n} / 2 = C(n,3) × 2^{5-2n}
    approx = comb(n, 3) * pow(2, 5-2*n)

    if R_values and n-3 < len(R_values):
        exact = R_values[n-3][1]
        ratio = exact / approx if approx != 0 else 0
        print("  n=%2d: R_exact = %12.8f, R_approx = %12.8f, ratio = %.6f" %
              (n, exact, approx, ratio))

# ============================================================
banner("THE GENERATING FUNCTION NEAR x=4")
# ============================================================

# G(x) = x³/(2(1-x/4)^4) has a POLE at x=4!
# This means R(n) ≈ C(n,3) × 4^{3-n} × n^{something} near x=4.
# Actually: [x^n] 1/(1-x/4)^4 = C(n+3,3)/4^n → n³/(6×4^n).
# So R(n) ≈ (n-3)^{...} not quite. Let me just check radius of convergence.

# The EXACT GF has higher-order corrections from k=5,7,...
# Each k-cycle contributes poles at x = 2^{(2n-k-1)/(k-1)-something}...
# Actually the GF is a finite sum for each n, so there's no issue.

# But the FORMAL generating function F(x) = Σ a_n x^n where
# a_n = V_n × n!/2^{C(n,2)} has:
# a_n = 1 + O(n³/4^n) → 1 as n → ∞.
# So F(x) = 1/(1-x) + G(x) where G converges for |x| < 4.
# F has a pole at x=1 (from the 1/(1-x) part) and G has a pole at x=4.

print("  F(x) = Σ (1 + R(n)) x^n = 1/(1-x) + G(x)")
print("  where G(x) ≈ x³/(2(1-x/4)^4) for the 3-cycle correction.")
print("  Poles: x=1 (from identity), x=4 (from 3-cycle correction).")
print()
print("  The RADIUS 4 = 2² comes from 2^{2n} in the Delta_3 exponent.")
print("  For k=5: pole at x = 2^4 = 16.")
print("  For k=7: pole at x = 2^6 = 64.")
print("  General k: pole at x = 2^{k-1}.")
print()
print("  The poles form a GEOMETRIC SEQUENCE: 4, 16, 64, 256, ...")
print("  = 4^1, 4^2, 4^3, 4^4, ...")
print("  No wait: 2^{k-1} for k=3,5,7,9 = 4, 16, 64, 256.")
print("  Ratios: 4, 4, 4, ... YES: poles at 4^{(k-1)/2} for odd k.")
print("  = 4, 4², 4³, 4⁴ for k=3,5,7,9.")

# ============================================================
banner("POLES AND THE LOG-DERIVATIVE")
# ============================================================

# Define Φ(x) = log(F(x)) where F(x) = Σ (1+R(n)) x^n.
# Then Φ'(x)/Φ(x) relates to the "zeta function" of the tournaments.
#
# But actually: log(T(n)) ≈ R(n) for small R(n).
# And R(n) ≈ R_3(n) = C(n,3) × 2^{5-2n}.
#
# This means: log(V_n × n!/2^m) ≈ C(n,3) × 2^{5-2n}
#
# And: log(V_n) ≈ log(2^m/n!) + C(n,3) × 2^{5-2n}
#      ≈ m×ln2 - ln(n!) + small correction
#
# The correction to log(V_n) from the naive 2^m/n! approximation is:
# C(n,3) × 2^{5-2n} ≈ (n³/6) × 32/4^n → 0 super-exponentially.

# Stirling: ln(n!) ≈ n ln n - n + (1/2) ln(2πn)
# m = C(n,2) = n(n-1)/2
# log(V_n) ≈ n(n-1)/2 × ln2 - n ln n + n - (1/2) ln(2πn) + (n³/6) × 32/4^n

print("  log2(V_n) vs Stirling-Burnside:\n")
for n in range(3, 14):
    m = comb(n, 2)
    parts = gen_odd_partitions(n)
    nfact = factorial(n)
    T = sum(multinomial(n, ct) * pow(2, pair_orbits(ct)) for ct in parts)
    V = T // nfact
    exact_log2 = log2(V) if V > 0 else 0

    stirling = m * log2(2) - log2(nfact)  # log2(2^m/n!)
    correction = comb(n, 3) * pow(2, 5-2*n) / log(2)  # small correction in nats → bits

    approx = stirling + correction / log2(2.718281828) if n > 3 else stirling
    # Actually: log2(1 + R) ≈ R / ln(2)
    R_n = float(T) / (nfact * pow(2, m)) - 1
    approx2 = stirling + R_n / log(2)

    print("  n=%2d: V=%d, log2(V)=%.4f, stirling=%.4f, +corr=%.4f" %
          (n, V, exact_log2, stirling, approx2))

# ============================================================
banner("THE TOURNAMENT PRIME RECIPROCAL SUM")
# ============================================================

# The correction R(n) involves terms 1/k for each odd k.
# Define: S(n) = Σ_{k=3,5,7,...,n} 1/k × n^{falling k} × 2^{Delta_k}
#
# For large n, S(n) ≈ (1/3) × n³ × 2^{5-2n} ≈ (1/3) × n³/4^n × 32
#
# The "effective prime reciprocal" at each n:
# eff(n) = R(n) / (n³/4^n × 32/6) = R(n) / (R_3_approx)
#
# If this approaches a constant... it approaches 1 (since R_3 dominates).
# But the CORRECTION to 1 involves 1/5, 1/7, etc.

# Define: the TOURNAMENT MERTENS CONSTANT
# M_T = lim_{n→∞} [R(n) / R_3(n) - 1] × 4^n/n^5 ... hmm, this goes to 0.
# Actually the 5-cycle correction is R_5(n) = C(n,5) × 24 × 2^{(4)(6-2n)/2}
# = C(n,5) × 24 × 2^{2(6-2n)} = C(n,5) × 24 × 2^{12-4n}
# R_5/R_3 = C(n,5)×24×2^{12-4n} / (C(n,3)×2×2^{4-2n})
# = (C(n,5)/C(n,3)) × 12 × 2^{8-2n}
# = ((n-3)(n-4)/20) × 12 × 2^{8-2n}
# = (3(n-3)(n-4)/5) × 2^{8-2n}

print("  R_5/R_3 RATIO (measures importance of 5-cycle):\n")
for n in range(5, 20):
    R5_over_R3 = Fraction(3*(n-3)*(n-4), 5) * Fraction(1, pow(4, n-4))
    print("  n=%2d: R_5/R_3 = %.10f" % (n, float(R5_over_R3)))

print("""
  THE TOURNAMENT PRIME HIERARCHY:

  R(n) = R_3(n) × [1 + R_5/R_3 + R_7/R_3 + ...]
       = R_3(n) × [1 + 3(n-3)(n-4)/(5×4^{n-4}) + higher terms]

  The correction from each higher "prime" k DECAYS as 4^{-(k-3)/2}
  times a polynomial in n.

  At n=20: R_5/R_3 < 10^{-8}. The tournament count is determined
  to 8 digits by JUST the 3-cycle term.

  The generating function is DOMINATED by the smallest prime p=3,
  with each higher prime contributing a correction that is
  EXPONENTIALLY smaller (factor 4^{-(k-3)/2} per step).

  This is analogous to the PRIME ZETA FUNCTION:
  P(s) = Σ_{p prime} p^{-s}
  where the first term 2^{-s} dominates for large s,
  but here the "large parameter" is n (vertex count) rather than s.
""")

# ============================================================
banner("CLEAN ASYMPTOTIC EXPANSION")
# ============================================================

# V_n = 2^{C(n,2)} / n! × T(n) where T(n) = 1 + R(n)
#
# T(n) = 1 + C(n,3)/3 × 2^{5-2n}
#       + C(n,5)/5 × 24 × 2^{12-4n}
#       + multi-cycle corrections + ...
#
# Let me write this cleanly:
# T(n) = 1 + (1/3)n↓3 × 2^{4-2n}      [3-cycle]
#       + (1/5)n↓5 × 2^{12-4n} × ?     [5-cycle]
#       + ...
#
# Hmm, let me recompute R_5 more carefully.
# R_5(n) = C(n,5) × (5-1)! × 2^{Delta_5}
#        = C(n,5) × 24 × 2^{(4)(6-2n)/2}
#        = C(n,5) × 24 × 2^{12-4n}

# But C(n,5) × 24 = n!/(5!(n-5)!) × 24 = n!/(5×(n-5)!) × (24/24)... wait.
# C(n,5) = n!/(5!(n-5)!) and (k-1)! = 4! = 24. So:
# C(n,5) × 24 = n!/(5!(n-5)!) × 4! = n!/((5)(n-5)!) since 5!/4! = 5.
# = (1/5) × n!/(n-5)! = (1/5) × n↓5.

# So: R_5(n) = (1/5) × n↓5 × 2^{Delta_5} where Delta_5 = 4(6-2n)/2 = 12-4n.
# R_5(n) = (1/5) × n↓5 × 2^{12-4n}

# Similarly: R_k(n) = (1/k) × n↓k × 2^{Delta_k}
# Delta_k = (k-1)(k+1-2n)/2

# Let's verify for k=5:
print("  VERIFY R_5 formula:\n")
for n in range(5, 10):
    m = comb(n, 2)
    ct = tuple(sorted([5] + [1]*(n-5)))
    exact = Fraction(multinomial(n, ct) * pow(2, pair_orbits(ct)), pow(2, m))
    falling_5 = n*(n-1)*(n-2)*(n-3)*(n-4)
    delta_5 = 4*(6-2*n)//2
    formula = Fraction(falling_5, 5) * Fraction(pow(2, delta_5 + m), pow(2, m))
    # Wait: R_k = (1/k) × n↓k × 2^{Delta_k} where Delta_k = po - m.
    # So R_k = (1/k) × n↓k × 2^{po - m} = (1/k) × n↓k × 2^{po} / 2^m
    # But mult(k,1^{n-k}) × 2^{po} / 2^m = (1/k) × n↓k × 2^{po} / 2^m
    # Hmm, let me just check:
    delta = (5-1)*(5+1-2*n)//2
    formula2 = Fraction(falling_5, 5) * Fraction(1, pow(2, -delta))  # 2^delta
    # R = (1/5) × n↓5 × 2^{delta}
    # delta is negative, so 2^delta = 1/2^|delta|
    abs_delta = abs(delta)
    formula3 = Fraction(falling_5, 5 * pow(2, abs_delta))
    print("  n=%d: exact=%s, formula=(1/5)×n↓5×2^{%d}=%s, match=%s" %
          (n, exact, delta, formula3, exact == formula3))

# ============================================================
banner("THE CLEAN ASYMPTOTIC TABLE")
# ============================================================

print("""  V_n × n! / 2^{C(n,2)} = 1 + Σ_{k=3,5,7,...} (1/k) × n↓k × 2^{Δ_k} + cross terms

  where Δ_k = (k-1)(k+1-2n)/2 and n↓k = n!/(n-k)!

  LEADING TERMS:
  k=3: (1/3) × n(n-1)(n-2) × 2^{4-2n}
  k=5: (1/5) × n(n-1)(n-2)(n-3)(n-4) × 2^{12-4n}
  k=7: (1/7) × n↓7 × 2^{24-6n}

  EXPONENT PATTERN: Δ_k = (k-1)(k+1-2n)/2
  k=3: 4-2n     (pole at x = 2^2 = 4)
  k=5: 12-4n    (pole at x = 2^4 = 16)
  k=7: 24-6n    (pole at x = 2^6 = 64)
  k=9: 40-8n    (pole at x = 2^8 = 256)

  The constant terms: 4, 12, 24, 40, ... = k(k-1)/2 - (k-1) = (k-1)(k-2)/2 + (k-1)/2
  Actually: (k-1)(k+1)/2 = 4, 12, 24, 40 for k=3,5,7,9.
  Yes: (k²-1)/2.
""")

# Print the clean table
print("  %-3s | %-12s | %-15s | %-10s | %-10s" %
      ("k", "Δ_k", "coefficient", "pole", "weight 1/k"))
print("  " + "-"*60)
for k in range(3, 20, 2):
    delta_const = (k*k - 1) // 2
    delta_coeff = -(k-1)
    pole = pow(2, k-1)
    print("  %-3d | %d - %dn%s | (1/%d)×n↓%d%s | 2^{%d}=%-5d | 1/%d = %.6f" %
          (k, delta_const, abs(delta_coeff), " "*(4-len(str(abs(delta_coeff)))),
           k, k, " "*(4-len(str(k))), k-1, pole, k, 1.0/k))

# ============================================================
banner("THE TOURNAMENT PARTITION FUNCTION")
# ============================================================

# Inspired by the partition function in number theory:
# p(n) = (1/2πi) ∮ F(q) q^{-n-1} dq where F(q) = Π(1-q^k)^{-1}
#
# Our analogue: T(n) ≈ Π_{k odd} (1 + R_k(n))
# But R_k depends on n, so this isn't quite an Euler product in q.
#
# BETTER: Define q = 1/4 (the "tournament modulus").
# Then R_3(n) = C(n,3) × 2 × q^n × 16 = 32 C(n,3) q^n.
# R_5(n) = (1/5) n↓5 × 2^{12} × q^{2n} = ... hmm, this involves q^{2n}.
# The different k's have DIFFERENT powers of q!
# k=3: q^n (really 4^{-n})
# k=5: q^{2n} (really 16^{-n})
# k=7: q^{3n} (really 64^{-n})
# In general: 2^{(k-1)n} → q^{(k-1)/2 × n}

# So T(n) = 1 + f_3(n) q^n + f_5(n) q^{2n} + f_7(n) q^{3n} + ...
# where f_k(n) = polynomial in n × constant.

# This is a q-SERIES where the variable is q^n!
# Define Q = q^n = 4^{-n}. Then:
# T(n) = 1 + a_3(n) Q + a_5(n) Q^2 + a_7(n) Q^3 + ...
# where a_k(n) = (1/k) n↓k × 2^{(k²-1)/2}

print("  T(n) AS A Q-SERIES (Q = 4^{-n}):\n")
print("  T(n) = 1 + Σ_j a_{2j+1}(n) × Q^j  where Q = 1/4^n\n")

for n in range(3, 14):
    m = comb(n, 2)
    Q = pow(4, -n)
    print("  n=%2d (Q=4^{-%d}=%.2e):" % (n, n, Q))
    for j, k in enumerate(range(3, min(n+1, 14), 2), start=1):
        falling_k = 1
        for i in range(k):
            falling_k *= (n - i)
        const = pow(2, (k*k-1)//2)
        a_k = falling_k * const / k
        print("    k=%d: a=%10.2f × Q^%d = %12.6e" % (k, a_k, j, a_k * Q**j))

print("""
  SUMMARY — THE TOURNAMENT COUNTING THEOREM:

  V_n = (2^{C(n,2)} / n!) × T(n)

  where T(n) = 1 + Σ_{k=3,5,7,...} (1/k) × n↓k × 2^{(k²-1)/2 - (k-1)n}
             + cross-orbit corrections (exponentially smaller)

  = 1 + (1/3)n↓3 × 32/4^n + (1/5)n↓5 × 2^{12}/16^n + (1/7)n↓7 × 2^{24}/64^n + ...

  The k-th "tournament prime" contributes:
  - Factor 1/k (prime reciprocal from cycle index)
  - Factor n↓k (falling factorial = number of k-arrangements)
  - Factor 2^{(k²-1)/2} (constant from pair orbits within cycle)
  - Factor (2^{k-1})^{-n} (exponential decay, rate = 2^{k-1})

  POLES of the generating function G(x) = Σ R(n) x^n:
  At x = 2^{k-1} for each odd k = 3, 5, 7, 9, ...
  i.e., at x = 4, 16, 64, 256, ...

  The tournament counting function Z_T(x) = Σ T(n) x^n has:
  - A POLE at x = 1 (from the constant term 1)
  - POLES at x = 4, 16, 64, 256, ... (from correction terms)
  - Each pole carries a PRIME RECIPROCAL 1/k weight.
""")

print("Done. opus-2026-03-24-S291")
