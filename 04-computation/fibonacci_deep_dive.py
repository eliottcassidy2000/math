#!/usr/bin/env python3
"""
fibonacci_deep_dive.py — opus-2026-03-12-S67c

DEEP EXPLORATION OF THE FIBONACCI CONNECTION

Key known facts:
- prod(1+Q_k) = F_p for Interval tournament
- F_p = B_m(1) where B_m is Morgan-Voyce polynomial
- Transfer matrix [[3,-1],[1,0]] has eigenvalues φ², ψ²
- Interval MINIMIZES prod(1+Q_k) yet MAXIMIZES H for p≥13

NEW DIRECTIONS TO EXPLORE:

1. FIBONACCI NUMBERS IN OTHER TOURNAMENT INVARIANTS
   - Is there a Fibonacci pattern in alpha_k (independent set counts)?
   - Does H(Int)/p relate to Fibonacci in a deeper way?

2. LUCAS NUMBERS AND THE COMPLEMENT
   - L_p = F_{p-1} + F_{p+1} = trace of transfer matrix ^m
   - Does prod(Q_k) = 1 (already known) relate to L_p?
   - prod(2+Q_k) = ?

3. FIBONACCI IN THE WALSH SPECTRUM
   - Do Walsh coefficients satisfy Fibonacci-type recurrences?
   - The degree-0 term h_hat[∅] = mean H — is it Fibonacci-related?

4. ZECKENDORF REPRESENTATION
   - Every positive integer has a unique Zeckendorf (Fibonacci base) representation
   - Can H be decomposed in the Fibonacci basis meaningfully?

5. THE FIBONACCI LATTICE
   - Fibonacci numbers appear in phyllotaxis as optimal packing
   - The 1D Fibonacci lattice has quasicrystalline order
   - Connection to the tournament's cycle packing structure?

6. FIBONACCI POLYNOMIALS
   - Generalize F_n to F_n(x) via recurrence F_{n+1}(x) = x F_n(x) + F_{n-1}(x)
   - prod(1+Q_k) at general x = ?
   - Does F_p(x) appear as a generating function for something?

7. BINET'S FORMULA AND GALOIS CONJUGATES
   - F_p = (φ^p - ψ^p)/√5 where ψ = (1-√5)/2
   - Since Q_k are algebraic integers in Q(cos 2π/p),
     how does the Galois action interact with Fibonacci?

8. PISANO PERIODS
   - F_n mod m has period π(m) (Pisano period)
   - π(p) divides p² - 1 (and equals p-1 or 2(p+1) depending on p mod 5)
   - Does this constrain H mod F_p?
"""

import numpy as np
from math import comb, gcd
from collections import Counter

def fib(n):
    a, b = 0, 1
    for _ in range(n):
        a, b = b, a+b
    return a

def lucas(n):
    a, b = 2, 1
    for _ in range(n):
        a, b = b, a+b
    return a

phi = (1 + np.sqrt(5)) / 2
psi = (1 - np.sqrt(5)) / 2

def legendre(a, p):
    a = a % p
    if a == 0: return 0
    ls = pow(a, (p-1)//2, p)
    return -1 if ls == p-1 else ls

def interval_Q(p):
    m = (p-1)//2
    return np.array([
        (np.sin(m*np.pi*k/p)/np.sin(np.pi*k/p))**2
        for k in range(1, m+1)
    ])

def morgan_voyce_B(m, x):
    """B_m(x) = sum_{j=0}^m C(m+j, 2j) x^j"""
    return sum(comb(m+j, 2*j) * x**j for j in range(m+1))

def morgan_voyce_b(m, x):
    """b_m(x) = sum_{j=0}^m C(m+j+1, 2j+1) x^j (the "small" MV polynomial)"""
    return sum(comb(m+j+1, 2*j+1) * x**j for j in range(m+1))

print("=" * 70)
print("FIBONACCI DEEP DIVE — opus-2026-03-12-S67c")
print("=" * 70)
print()

# ============================================================
# 1. FIBONACCI IN prod(a + Q_k) for various a
# ============================================================

print("PART 1: prod(a + Q_k) FOR VARIOUS a")
print("=" * 70)
print()
print("prod(1+Q_k) = F_p. What about prod(a+Q_k) for other a?")
print()

for p in [5, 7, 11, 13, 17]:
    m = (p-1)//2
    Qs = interval_Q(p)

    Fp = fib(p)
    Lp = lucas(p)

    print(f"p={p}, m={m}:")

    for a in [-2, -1, 0, 1, 2, 3, 4, 5]:
        prod_val = np.prod(a + Qs)
        # Check if this is a known sequence
        # prod(0+Q) = prod(Q) = 1
        # prod(1+Q) = F_p
        # prod(-1+Q) = prod(Q-1) = P_m(1) = period-6 pattern
        # prod(2+Q) = ?

        # Morgan-Voyce: B_m(x) = prod(Q_k + 1) when x=1 is wrong...
        # Actually P_m(t) = prod(t - Q_k) = sum (-1)^j e_j t^{m-j}
        # prod(a+Q_k) = (-1)^m P_m(-a) (since prod(a+Q_k) = prod(-(−a−Q_k)) = (-1)^m prod(-a-Q_k))
        # Wait: prod(a+Q_k) = prod(a+Q_k). And prod(t-Q_k) = P_m(t).
        # So prod(a+Q_k) = prod(-(Q_k-a)) = (-1)^m prod(Q_k - a) = (-1)^m P_m(a)...
        # No: P_m(a) = prod(a - Q_k), so prod(a + Q_k) ≠ P_m(a).
        # prod(a + Q_k) = prod(-(Q_k - (-a))) = ... hmm.
        # Let's just compute P_m(-a) = prod(-a - Q_k) = (-1)^m prod(a+Q_k)
        # So prod(a+Q_k) = (-1)^m P_m(-a)

        Pm_neg_a = sum((-1)**j * comb(m+j, 2*j) * (-a)**(m-j) for j in range(m+1))
        expected = (-1)**m * Pm_neg_a

        print(f"  prod({a:+d}+Q) = {prod_val:>14.1f}  P_m({-a:+d}) = {Pm_neg_a:>14}  (-1)^m·P = {expected:>14}")

    # Check specific values
    prod_2 = np.prod(2 + Qs)
    prod_3 = np.prod(3 + Qs)

    # Is prod(2+Q) related to Lucas?
    Bm_2 = morgan_voyce_B(m, 2)  # B_m(2)
    print(f"  B_m(1) = {morgan_voyce_B(m, 1)} = F_p = {Fp}")
    print(f"  B_m(2) = {Bm_2}")
    print(f"  prod(2+Q) = {prod_2:.1f}")
    print(f"  prod(3+Q) = {prod_3:.1f}")

    # B_m(x) = prod(1 + x*Q_k)? No... let me think.
    # Actually P_m(t) = prod(t - Q_k) and the e_j = C(m+j,2j)
    # prod(1+Q_k) = (-1)^m P_m(-1)
    # P_m(-1) = sum (-1)^j C(m+j,2j) (-1)^{m-j} = (-1)^m sum C(m+j,2j) = (-1)^m B_m(1)
    # So prod(1+Q_k) = B_m(1) = F_{2m+1} ✓

    # prod(2+Q_k) = (-1)^m P_m(-2)
    # P_m(-2) = sum (-1)^j C(m+j,2j) (-2)^{m-j}
    #         = (-1)^m sum (-1)^{m-j} (-1)^j C(m+j,2j) 2^{m-j} * (-1)^{m-j}
    # Hmm let me just compute directly
    # P_m(-2) = sum_{j=0}^m (-1)^j C(m+j,2j) (-2)^{m-j}
    #         = sum_{j=0}^m (-1)^j C(m+j,2j) (-1)^{m-j} 2^{m-j}
    #         = (-1)^m sum_{j=0}^m C(m+j,2j) 2^{m-j}   [(-1)^j (-1)^{m-j} = (-1)^m]
    #         = (-1)^m B_m^*(2) where B_m^*(x) = sum C(m+j,2j) x^{m-j}
    # So prod(2+Q_k) = (-1)^m * (-1)^m * B_m^*(2) = B_m^*(2) = sum C(m+j,2j) 2^{m-j}

    # This is a different polynomial! It's x^m B_m(1/x) evaluated at x=2
    # = 2^m B_m(1/2)
    prod_2_formula = 2**m * morgan_voyce_B(m, 0.5)
    print(f"  2^m B_m(1/2) = {prod_2_formula:.1f}")

    # More generally: prod(a+Q_k) = a^m B_m(1/a)
    # Because P_m(-a) = sum (-1)^j C(m+j,2j) (-a)^{m-j}
    # = (-1)^m sum C(m+j,2j) a^{m-j} = (-1)^m a^m sum C(m+j,2j) (1/a)^j
    # = (-1)^m a^m B_m(1/a)
    # And prod(a+Q_k) = (-1)^m P_m(-a) = a^m B_m(1/a)

    print(f"\n  IDENTITY: prod(a+Q_k) = a^m · B_m(1/a)")
    for a in [1, 2, 3, 4, 5]:
        lhs = np.prod(a + Qs)
        rhs = a**m * morgan_voyce_B(m, 1/a)
        print(f"    a={a}: prod = {lhs:.4f}, a^m B_m(1/a) = {rhs:.4f}, match = {abs(lhs-rhs) < 0.01}")

    print()

# ============================================================
# 2. THE GENERATING FUNCTION B_m(x) AND ITS ROOTS
# ============================================================

print("=" * 70)
print("PART 2: B_m(x) AS A GENERATING FUNCTION")
print("=" * 70)
print()
print("prod(a + Q_k) = a^m B_m(1/a)")
print("Setting x = 1/a: prod(1/x + Q_k) = (1/x)^m B_m(x)")
print("So B_m(x) = x^m prod(Q_k + 1/x)")
print()
print("The ZEROS of B_m(x) are at x_k = -1/Q_k!")
print("Since all Q_k > 0, all zeros of B_m are REAL and NEGATIVE.")
print()

for p in [7, 11, 13]:
    m = (p-1)//2
    Qs = interval_Q(p)

    # Zeros of B_m(x)
    zeros = [-1/Q for Q in Qs]
    print(f"  p={p}: B_m zeros = {np.sort(zeros)}")
    print(f"         Q_k = {np.sort(Qs)[::-1]}")

    # B_m(x) at special values
    print(f"  B_m(-1) = {morgan_voyce_B(m, -1)}")
    print(f"  B_m(0)  = {morgan_voyce_B(m, 0)} (should be 1)")
    print(f"  B_m(1)  = {morgan_voyce_B(m, 1)} = F_p = {fib(p)}")
    print(f"  B_m(2)  = {morgan_voyce_B(m, 2)}")
    print(f"  B_m(3)  = {morgan_voyce_B(m, 3)}")
    print(f"  B_m(4)  = {morgan_voyce_B(m, 4)}")

    # Pell numbers? B_m(2) might relate to Pell
    # Pell: P_0=0, P_1=1, P_{n+1} = 2P_n + P_{n-1}
    # Companion Pell (half-companion): H_n: H_0=1, H_1=1, H_{n+1}=2H_n+H_{n-1}
    def pell(n):
        a, b = 0, 1
        for _ in range(n): a, b = b, 2*b+a
        return a

    def pell_companion(n):
        a, b = 1, 1
        for _ in range(n): a, b = b, 2*b+a
        return a

    print(f"  Pell_{{2m+1}} = {pell(2*m+1)}")
    print(f"  B_m(2) = {morgan_voyce_B(m, 2)}")

    # Check B_m(2) = P_{2m+1} (Pell)?
    # B_m(x) satisfies B_{m+1}(x) = (2+x)B_m(x) - B_{m-1}(x)
    # At x=2: B_{m+1}(2) = 4·B_m(2) - B_{m-1}(2)
    # Pell: P_{n+1} = 2P_n + P_{n-1}... different recurrence!

    # Actually at x=2: recurrence is a_{n+1} = 4a_n - a_{n-1}
    # This gives the bisection of Pell: P_{2n}
    # Let me check
    print(f"  Check: is B_m(2) related to Pell?")
    Bm2_seq = [morgan_voyce_B(k, 2) for k in range(8)]
    print(f"    B_0..B_7 at x=2: {Bm2_seq}")
    pell_seq = [pell(k) for k in range(16)]
    print(f"    Pell_0..15: {pell_seq}")

    # a_n: 1, 3, 11, 41, 153, 571, ...
    # These are denominators of convergents of sqrt(2)
    # a_n = ((1+sqrt(2))^{2n+1} + (1-sqrt(2))^{2n+1}) / 2
    # = half-integer Chebyshev T at x=3/2?

    print()

# ============================================================
# 3. B_m(x) GENERATING FUNCTION FOR GENERAL x
# ============================================================

print("=" * 70)
print("PART 3: PROD(a+Q_k) = a^m B_m(1/a) — THE FULL PICTURE")
print("=" * 70)
print()
print("This gives us a PARAMETERIZED FAMILY of identities:")
print()

for p in [7, 11, 13]:
    m = (p-1)//2

    print(f"p={p}, m={m}:")
    print(f"  a=1: prod(1+Q) = F_{{2m+1}} = F_{p} = {fib(p)}")

    # a → ∞: prod(a+Q) ≈ a^m (since Q_k << a)
    # a → 0: prod(Q_k) = 1

    # Derivative at a=1: d/da [a^m B_m(1/a)] at a=1
    # = m a^{m-1} B_m(1/a) + a^m B_m'(1/a) (-1/a²)
    # At a=1: m B_m(1) - B_m'(1)
    # = m F_p - B_m'(1)

    # B_m'(x) = sum_{j=1}^m j C(m+j, 2j) x^{j-1}
    Bm_prime_1 = sum(j * comb(m+j, 2*j) for j in range(1, m+1))

    # Also: d/da prod(a+Q_k) at a=1 = prod(1+Q_k) * sum 1/(1+Q_k)
    # = F_p * sum 1/(1+Q_k)
    sum_recip_1pQ = sum(1/(1+Q) for Q in interval_Q(p))
    deriv_direct = fib(p) * sum_recip_1pQ
    deriv_formula = m * fib(p) - Bm_prime_1

    print(f"  d/da prod(a+Q)|_{{a=1}} = {deriv_direct:.4f}")
    print(f"  = m·F_p - B'_m(1) = {m}·{fib(p)} - {Bm_prime_1} = {deriv_formula}")
    print(f"  F_p · sum 1/(1+Q_k) = {fib(p)} · {sum_recip_1pQ:.6f} = {deriv_direct:.4f}")
    print(f"  sum 1/(1+Q_k) = {sum_recip_1pQ:.6f}")

    # Is sum 1/(1+Q_k) a nice number?
    # We know sum 1/Q_k = p-2 and sum Q_k = m(m+1)/2
    # sum 1/(1+Q_k) = sum Q_k/(Q_k(1+Q_k))... not obviously nice

    # But from B_m: sum 1/(1+Q_k) = (m F_p - B'_m(1)) / F_p
    print(f"  = m - B'_m(1)/F_p = {m} - {Bm_prime_1/fib(p):.6f}")
    print(f"  B'_m(1)/F_p = {Bm_prime_1/fib(p):.6f}")

    # Check: B'_m(1)/B_m(1) is the log derivative
    # This relates to the digamma function of the Fibonacci
    print(f"  log-deriv B'_m(1)/B_m(1) = {Bm_prime_1/morgan_voyce_B(m, 1):.6f}")
    print()

# ============================================================
# 4. FIBONACCI WORD AND TOURNAMENT STRUCTURE
# ============================================================

print("=" * 70)
print("PART 4: FIBONACCI WORD AND TOURNAMENT MORPHISM")
print("=" * 70)
print()
print("The Fibonacci word is the infinite binary word obtained from")
print("the substitution S: 0→01, 1→0.")
print("It encodes the 1D quasicrystal (Fibonacci chain).")
print()
print("Connection: the Interval tournament has generating set")
print("S = {1,...,m} = the 'interval' in Z_p. The complement is")
print("{m+1,...,p-1}. This partition of Z_p* into two halves")
print("is the simplest possible — analogous to the Fibonacci")
print("word's partition into 0s and 1s.")
print()

# Generate Fibonacci words
def fib_word(n):
    """Generate n-th Fibonacci word by substitution."""
    w = "0"
    for _ in range(n):
        w = w.replace("0", "01x").replace("1", "0x").replace("x", "")
    return w

for k in range(2, 8):
    w = fib_word(k)
    print(f"  F-word {k}: {w[:50]}... (length {len(w)})")

print()

# The Fibonacci word has letter frequencies ratio φ:1
# For the Interval tournament at prime p, the generating set has
# |S| = m = (p-1)/2, and the complement has m elements too.
# So the frequencies are equal (50-50), unlike the Fibonacci word.

# But the INTERNAL structure of S = {1,...,m} IS Fibonacci-like:
# The gaps in the circulant are all 1 (within S) and then one big gap.
# The Dirichlet kernel concentrates power at low frequencies,
# like the Fibonacci lattice concentrates diffraction peaks.

print("Key analogy:")
print("  Fibonacci chain: 01001010010010100101...")
print("  has diffraction peaks at positions n/φ mod 1")
print()
print("  Interval tournament: eigenvalues Q_k peak at k=1")
print("  with Q_1/Q_2 → ∞ as p → ∞")
print()
print("  Both exhibit 'quasicrystalline order':")
print("  - Not periodic (the Q_k are all different)")
print("  - Not random (they follow a precise Chebyshev formula)")
print("  - Self-similar (Morgan-Voyce recurrence = self-similar structure)")
print()

# ============================================================
# 5. FIBONACCI IN H: DECOMPOSITION
# ============================================================

print("=" * 70)
print("PART 5: H IN THE FIBONACCI BASIS")
print("=" * 70)
print()

# Every integer has a unique Zeckendorf representation as sum of
# non-consecutive Fibonacci numbers.
def zeckendorf(n):
    """Return Zeckendorf (Fibonacci base) representation."""
    if n == 0: return []
    fibs = []
    a, b = 1, 2
    while b <= n:
        a, b = b, a+b
    # Now a is the largest Fibonacci ≤ n
    result = []
    while n > 0:
        if a <= n:
            result.append(a)
            n -= a
        a_new = (a + round(a / phi)) // 2  # Previous Fibonacci... hacky
        # Better: recompute
        while a > n and a > 1:
            a, b = round(a / phi), a
        if a <= 0: break
    return result

# Actually let me do this properly
def zeckendorf_repr(n):
    """Return Zeckendorf representation as list of Fibonacci indices."""
    if n <= 0: return []
    # Build Fibonacci list
    fibs = [1, 2]
    while fibs[-1] < n:
        fibs.append(fibs[-1] + fibs[-2])

    result = []
    for f in reversed(fibs):
        if f <= n:
            result.append(f)
            n -= f
    return result

# H values for Interval
known_H_int = {5: 15, 7: 175, 11: 93027, 13: 3711175}

for p, H in known_H_int.items():
    m = (p-1)//2
    zeck = zeckendorf_repr(H)
    Fp = fib(p)

    print(f"  p={p}: H = {H}")
    print(f"    Zeckendorf: {' + '.join(str(z) for z in zeck)}")
    print(f"    H / F_p = {H / Fp:.6f}")
    print(f"    H mod F_p = {H % Fp}")

    # Fibonacci quotient and remainder
    q, r = divmod(H, Fp)
    print(f"    H = {q} · F_{p} + {r}")

    # Is H/p a Fibonacci or Lucas number?
    h = H // p
    print(f"    H/p = {h}")

    # Check nearby Fibonacci and Lucas numbers
    k = 1
    while fib(k) < h: k += 1
    print(f"    Nearest F: F_{k-1}={fib(k-1)}, F_{k}={fib(k)}")

    k = 1
    while lucas(k) < h: k += 1
    print(f"    Nearest L: L_{k-1}={lucas(k-1)}, L_{k}={lucas(k)}")
    print()

# ============================================================
# 6. THE GOLDEN RATIO IN THE Q-SPECTRUM
# ============================================================

print("=" * 70)
print("PART 6: GOLDEN RATIO IN THE EIGENVALUE SPECTRUM")
print("=" * 70)
print()

for p in [7, 11, 13, 17, 23, 29, 43]:
    m = (p-1)//2
    Qs = interval_Q(p)
    Qs_sorted = np.sort(Qs)[::-1]

    # Ratios of consecutive Q values
    ratios = [Qs_sorted[i] / Qs_sorted[i+1] for i in range(min(5, len(Qs_sorted)-1))]

    # Check for golden ratio patterns
    # Q_1 / sum(Q_k) = Q_1 / (m(m+1)/2) → 2/(m+1) * Q_1/m
    frac = Qs_sorted[0] / np.sum(Qs)

    # The key ratio: prod(1+Q_k)^{1/m} = F_p^{1/m} → φ²
    geom_mean = fib(p) ** (1/m)

    print(f"  p={p:2d}: F_p^(1/m) = {geom_mean:.6f} (→ φ² = {phi**2:.6f})")
    print(f"    Q_1/Q_2 = {ratios[0]:.4f}")
    if len(ratios) > 1:
        print(f"    Q_2/Q_3 = {ratios[1]:.4f}")

print()
print(f"  ASYMPTOTIC: F_p^(1/m) → φ² = {phi**2:.6f}")
print(f"  because F_p ≈ φ^p/√5 and p = 2m+1, so")
print(f"  F_p^(1/m) ≈ (φ^(2m+1)/√5)^(1/m) → φ^2 = {phi**2:.6f}")
print()

# ============================================================
# 7. FIBONACCI DIVISIBILITY AND TOURNAMENT STRUCTURE
# ============================================================

print("=" * 70)
print("PART 7: FIBONACCI DIVISIBILITY LATTICE")
print("=" * 70)
print()
print("Key property: F_m | F_n iff m | n.")
print("So F_p is always coprime to F_q for distinct primes p, q > 5.")
print()
print("The Pisano period π(p) = period of F_n mod p:")
print()

for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]:
    # Compute Pisano period
    a, b = 0, 1
    period = 0
    for i in range(1, 6*p):
        a, b = b, (a+b) % p
        if a == 0 and b == 1:
            period = i
            break

    Fp_mod_p = fib(p) % p
    # Legendre symbol (5/p) determines whether p | F_{p-1} or p | F_{p+1}
    leg5 = legendre(5, p)

    print(f"  p={p:2d}: π(p) = {period:4d}, F_p mod p = {Fp_mod_p}, (5/p) = {leg5:+d}")
    print(f"    p divides F_{{p-{1 if leg5==1 else '(p+1)'}}} = F_{p-1 if leg5==1 else p+1}")

print()
print("  F_p ≡ (5/p) mod p (by Fibonacci reciprocity law)")
print("  So prod(1+Q_k) ≡ (5/p) mod p")
print()
print("  This constrains the Q-spectrum modularly!")
print()

# ============================================================
# 8. CREATIVE SYNTHESIS: THE FIBONACCI-TOURNAMENT DICTIONARY
# ============================================================

print("=" * 70)
print("PART 8: FIBONACCI-TOURNAMENT DICTIONARY")
print("=" * 70)
print()
print("| Fibonacci/Golden Ratio     | Tournament                          |")
print("|----------------------------|-------------------------------------|")
print("| F_p = prod(1+Q_k)          | Eigenvalue product of Interval      |")
print("| φ² = lim F_p^{1/m}         | Geometric mean of 1+Q_k             |")
print("| B_m(x) = Morgan-Voyce      | Generating function of e_j           |")
print("| [[3,-1],[1,0]]^m           | Transfer matrix of ladder network   |")
print("| Zeckendorf representation  | (Open) Decomposition of H?          |")
print("| F_m | F_n iff m | n        | (Open) Divisibility in H?           |")
print("| Pisano period π(p)         | Periodicity of F_p mod p            |")
print("| Fibonacci word (quasicryst)| Dirichlet kernel (aperiodic order)  |")
print("| φ = (1+√5)/2               | Eigenvalue of transfer matrix       |")
print("| F_{p-1} + F_{p+1} = L_p    | Lucas = trace of T^m                |")
print("| 1/φ + 1/φ² = 1             | Eigenvalue reciprocal sum = 1       |")
print("| F_n² + F_{n+1}² = F_{2n+1} | Parseval-type identity for Q_k?     |")
print()
print("DEEPEST CONNECTION:")
print("  The Interval tournament IS a 1D Fibonacci quasicrystal")
print("  in frequency space. Its eigenvalues (Q_k) are the")
print("  'phonon frequencies' of the Fibonacci chain. The identity")
print("  prod(1+Q_k) = F_p says: the thermal partition function")
print("  of this quasicrystal at temperature T = 1/ln(2) equals")
print("  the p-th Fibonacci number.")
print()
print("UNEXPLORED DIRECTIONS:")
print("  1. Does the Fibonacci LATTICE (2D, Penrose tiling)")
print("     have an analogue for 2D tournaments?")
print("  2. The golden ratio appears in KAM theory (stability")
print("     of Hamiltonian systems). Does tournament stability")
print("     (under perturbation of S) follow golden-ratio scaling?")
print("  3. Fibonacci numbers count tilings of 1×n rectangles.")
print("     Does F_p count tilings related to the tournament?")
print("  4. The Binet formula F_p = (φ^p - ψ^p)/√5 means")
print("     prod(1+Q_k) splits into TWO Galois-conjugate parts.")
print("     This might connect to the Galois structure of Q(cos 2π/p).")
print()
