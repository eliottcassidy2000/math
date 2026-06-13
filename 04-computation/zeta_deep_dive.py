"""Deep dive: k-nacci traces, forbidden values, and L-functions."""
import numpy as np
from math import factorial, gcd
from fractions import Fraction

# ═══════════════════════════════════════════════════════════
# 1. COMPLETE k-NACCI TRACE TABLE — WHERE DO 7,21,42 APPEAR?
# ═══════════════════════════════════════════════════════════

print("="*70)
print("k-NACCI TRACES: COMPLETE TABLE (k=2..8, n=1..15)")
print("="*70)

forbidden = {7, 21}
special = {42, 31, 63, 127, 255}

for k in range(2, 9):
    M = np.zeros((k, k), dtype=np.int64)
    for i in range(k-1):
        M[i][i+1] = 1
    for j in range(k):
        M[k-1][j] = 1
    
    Mpow = np.eye(k, dtype=np.int64)
    traces = []
    for n in range(1, 16):
        Mpow = Mpow @ M
        tr = int(np.trace(Mpow))
        traces.append(tr)
    
    print(f"k={k}: {traces}")
    for n, tr in enumerate(traces, 1):
        if tr in forbidden:
            print(f"  *** Tr(M_{k}^{n}) = {tr} (FORBIDDEN H) ***")
        if tr == 42:
            print(f"  *** Tr(M_{k}^{n}) = 42 (BASE-42!) ***")

# Key finding from above: Tr(M_3^5) = 21
print("\n--- CRITICAL: Tr(M_3^n) sequence ---")
M3 = np.array([[0,1,0],[0,0,1],[1,1,1]], dtype=np.int64)
Mpow = np.eye(3, dtype=np.int64)
print("Tribonacci traces: ", end="")
trib_traces = []
for n in range(1, 20):
    Mpow = Mpow @ M3
    tr = int(np.trace(Mpow))
    trib_traces.append(tr)
print(trib_traces)

# The tribonacci numbers themselves
print("\nTribonacci numbers: ", end="")
trib = [0, 0, 1]
for i in range(20):
    trib.append(trib[-1] + trib[-2] + trib[-3])
print(trib[3:23])

# Relationship: Tr(M_3^n) vs tribonacci(n)
print("\nTr(M_3^n) vs T(n+2):")
for n in range(1, 15):
    print(f"  n={n:2d}: Tr={trib_traces[n-1]:>6d}, T(n+2)={trib[n+2]:>6d}, "
          f"ratio={trib_traces[n-1]/trib[n+2] if trib[n+2] else 'n/a':.4f}")

# ═══════════════════════════════════════════════════════════
# 2. FORBIDDEN VALUES AS k-NACCI EVALUATIONS
# ═══════════════════════════════════════════════════════════

print("\n" + "="*70)
print("FORBIDDEN VALUES IN k-NACCI STRUCTURE")
print("="*70)

# 7 = Tr(M_3^3) = 2^3 - 1 (THM-227)
# 21 = Tr(M_3^5) = ???
# What's the pattern?

print("\nTr(M_3^n) mod small numbers:")
for mod in [2, 3, 6, 7, 42]:
    residues = [t % mod for t in trib_traces[:20]]
    print(f"  mod {mod:2d}: {residues}")

# Factor 21 = 3 × 7
# Tr(M_3^3) = 7, Tr(M_3^5) = 21 = 3 × 7
# Is there a multiplicative relationship?
print(f"\nTr(M_3^3) = 7")
print(f"Tr(M_3^5) = 21 = 3 × Tr(M_3^3)")
print(f"Tr(M_3^7) = {trib_traces[6]}")
print(f"Tr(M_3^9) = {trib_traces[8]}")

# Check if Tr(M_3^{2k+1}) has pattern
print("\nTr(M_3^n) for ODD n:")
for n in range(1, 20, 2):
    if n <= len(trib_traces):
        print(f"  Tr(M_3^{n:2d}) = {trib_traces[n-1]:>6d}", end="")
        # Factor it
        val = trib_traces[n-1]
        if val > 1:
            factors = []
            v = val
            for p in range(2, val+1):
                while v % p == 0:
                    factors.append(p)
                    v //= p
                if v == 1:
                    break
            print(f"  = {'×'.join(map(str, factors))}", end="")
        print()

# ═══════════════════════════════════════════════════════════
# 3. NEWTON'S IDENTITIES AND POWER SUM SYMMETRIC FUNCTIONS
# ═══════════════════════════════════════════════════════════

print("\n" + "="*70)
print("NEWTON'S IDENTITIES: TRACES AS POWER SUM SYMMETRIC FUNCTIONS")
print("="*70)

# M_k has characteristic polynomial x^k - x^{k-1} - ... - x - 1
# Eigenvalues λ_1, ..., λ_k satisfy:
#   e_1 = λ_1+...+λ_k = 1
#   e_2 = Σ λ_iλ_j = -1
#   ...
#   e_k = λ_1·...·λ_k = (-1)^{k+1}
# (all elementary symmetric functions are ±1)
#
# Newton: p_n = Σ λ_i^n (power sums = traces)
# p_1 = e_1 = 1
# p_2 = e_1·p_1 - 2e_2 = 1 + 2 = 3
# p_3 = e_1·p_2 - e_2·p_1 + 3e_3 = 3+1+3 = 7 (for k≥3)
# p_n = p_{n-1} + p_{n-2} + ... + p_{n-k} (for n > k)

# For k=3: p_n = p_{n-1} + p_{n-2} + p_{n-3}
# p_1=1, p_2=3, p_3=7
# p_4 = 7+3+1 = 11
# p_5 = 11+7+3 = 21
# p_6 = 21+11+7 = 39

print("\nFor k=3 (tribonacci companion matrix):")
print("Eigenvalues satisfy x³ - x² - x - 1 = 0")
print("Power sums p_n = Tr(M_3^n) satisfy: p_n = p_{n-1} + p_{n-2} + p_{n-3}")
print()
print("p_1 = 1")
print("p_2 = 3  (via Newton: e_1·p_1 - 2e_2 = 1·1 - 2·(-1) = 3)")
print("p_3 = 7  (via Newton: e_1·p_2 - e_2·p_1 + 3e_3 = 3 + 1 + 3 = 7)")
print("p_4 = 11 (tribonacci recurrence: 7+3+1)")
print("p_5 = 21 (tribonacci recurrence: 11+7+3)")
print("p_6 = 39 (tribonacci recurrence: 21+11+7)")
print("p_7 = 71 (tribonacci recurrence: 39+21+11)")
print()
print("KEY INSIGHT: 7 and 21 are CONSECUTIVE TRIBONACCI POWER SUMS!")
print("  7 = p_3 = 2³ - 1 (Mersenne, THM-227)")
print("  21 = p_5 = 7 + 11 + 3 = p_3 + p_4 + p_2")
print()
print("The forbidden pair {7, 21} corresponds to {p_3, p_5} in the")
print("tribonacci trace sequence. Both are at ODD indices.")

# ═══════════════════════════════════════════════════════════
# 4. WHY p_3 AND p_5 ARE SPECIAL
# ═══════════════════════════════════════════════════════════

print("\n" + "="*70)
print("WHY p_3=7 AND p_5=21 ARE FORBIDDEN BUT p_7=71 IS NOT")
print("="*70)

# At n = 2k+1 vertices, the achievable H values need at least k+1 
# In general, tournament on n vertices has H ≥ 1 and H odd
# Forbidden means: no tournament on ANY n has H = that value

# p_1 = 1: achievable (transitive tournament on any n)
# p_2 = 3: achievable (C_3, or transitive on n≥4)
# p_3 = 7: FORBIDDEN (proved via cycle argument)
# p_4 = 11: achievable (various tournaments at n=5)
# p_5 = 21: FORBIDDEN (proved via poisoning DAG)
# p_6 = 39: achievable (at n=6)
# p_7 = 71: achievable (at n=7)

print("Tribonacci trace p_n and H-achievability:")
achievable_check = {1: True, 3: True, 7: False, 11: True, 21: False,
                    39: True, 71: True}
for n in range(1, 8):
    tr = trib_traces[n-1]
    ach = achievable_check.get(tr, "?")
    print(f"  p_{n} = {tr:>4d}  achievable: {ach}")

# Compute: what's special about indices 3 and 5?
# THM-227: p_n = 2^n - 1 for n ≤ k=3
# So p_3 = 7 = 2^3-1 is "maximal Mersenne" in the Newton range
# For n > k: p_n follows the tribonacci recurrence instead
#
# p_5 = p_4 + p_3 + p_2 = 11 + 7 + 3 = 21 = 3·7
# The key: 21 = p_2 · p_3 = 3 × 7!
print(f"\n21 = p_2 × p_3 = 3 × 7")
print(f"Is this a general pattern? p_{2}·p_{3} = {trib_traces[1]*trib_traces[2]}")
print(f"p_{3}·p_{4} = {trib_traces[2]*trib_traces[3]} = 77")

# Check: is p_m · p_n ever = p_k for other indices?
print("\nProduct relationships p_a · p_b = p_c:")
for a in range(1, 10):
    for b in range(a, 10):
        prod = trib_traces[a-1] * trib_traces[b-1]
        for c in range(1, 15):
            if c-1 < len(trib_traces) and trib_traces[c-1] == prod:
                print(f"  p_{a} · p_{b} = {trib_traces[a-1]} × {trib_traces[b-1]} = {prod} = p_{c}")

# ═══════════════════════════════════════════════════════════
# 5. RIEMANN-VON MANGOLDT ANALOGY
# ═══════════════════════════════════════════════════════════

print("\n" + "="*70)
print("RIEMANN-VON MANGOLDT ANALOGY")
print("="*70)

print("""
Von Mangoldt: log ζ(s) = Σ_p Σ_k p^{-ks}/k  (Euler product log)
              -ζ'/ζ(s) = Σ_n Λ(n)/n^s  (Dirichlet series for Λ)

Tournament analogue:
  log H(T) = log I(Ω(T), 2) 
            = Σ_{cliques C in Ω} log(1 + 2^{|C|})
            ≈ Σ_{cycles c} 2^{|independent set containing c|}/H(T)

The "primes" of tournament structure are ODD CYCLES.
The "von Mangoldt function" Λ_T counts their contribution.

Explicit analogy table:
  ζ(s)         ↔  H(T)         (the central object)
  primes p     ↔  odd cycles C  (the building blocks)
  p^{-s}       ↔  2^{|C|}       (the weight)
  log ζ        ↔  log H          
  Λ(n)         ↔  α_k           (cycle independence counts)
  trivial zeros ↔  even cycles   (contribute nothing)
  RH           ↔  THM-H         (cancellation bound)
  critical strip ↔ Walsh degrees (where information lives)
""")

# ═══════════════════════════════════════════════════════════
# 6. DIRICHLET L-FUNCTION AND PALEY BETTI FORMULA
# ═══════════════════════════════════════════════════════════

print("="*70)
print("DIRICHLET L-FUNCTIONS AND PALEY BETTI NUMBERS")
print("="*70)

def legendre(a, p):
    if a % p == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1

# For T_p, the Betti numbers are β_m = m(m-3)/2, β_{m+1} = C(m+1,2)
# where m = (p-1)/2
# 
# The eigenvalues of the adjacency matrix are:
#   λ_0 = (p-1)/2 = m (the "trivial" eigenvalue)
#   λ_j = (-1 + χ(j)√p)/2 for j=1,...,p-1
#   where χ(j) = (j/p) is the Legendre symbol
#
# So there are exactly 2 distinct non-trivial eigenvalues:
#   λ_+ = (-1+√p)/2  with multiplicity m (for QRs)
#   λ_- = (-1-√p)/2  with multiplicity m (for NQRs)
#
# The Gauss sum g = Σ χ(a) ω^a satisfies g² = χ(-1)p = -p (for p≡3 mod 4)
# So g = i√p (up to sign)

for p in [7, 11, 23, 43]:
    if not all(p % i != 0 for i in range(2, int(p**0.5)+1)):
        continue
    m = (p-1)//2
    
    # L-function values
    L_vals = {}
    for s in [1, 2, 3]:
        L = sum(legendre(n, p) / n**s for n in range(1, 100000))
        L_vals[s] = L
    
    # Class number (p ≡ 3 mod 4)
    if p % 4 == 3:
        h = round(abs(L_vals[1]) * p**0.5 / np.pi)
    else:
        h = "N/A"
    
    beta_m = m*(m-3)//2 if m >= 3 else 0
    beta_mp1 = m*(m+1)//2
    chi = 1 - beta_m + beta_mp1
    
    print(f"\nPaley T_{p} (m={m}):")
    print(f"  β_{m} = {beta_m} = m(m-3)/2")
    print(f"  β_{m+1} = {beta_mp1} = C(m+1,2)")
    print(f"  χ = {chi} (should = p = {p})")
    print(f"  h(-{p}) = {h}")
    print(f"  L(1,χ) = {L_vals[1]:.6f}")
    print(f"  L(2,χ) = {L_vals[2]:.6f}")
    
    # KEY QUESTION: Is β_m related to L-function?
    # β_m = m(m-3)/2 = (m²-3m)/2
    # h(-p) = w·√p·L(1,χ)/(2π) where w = number of roots of unity
    # For p>3: w=2, so h = √p·L(1,χ)/π
    if p % 4 == 3 and p > 3:
        h_exact = p**0.5 * abs(L_vals[1]) / np.pi
        print(f"  h(-{p}) exact = {h_exact:.4f}")
        print(f"  β_{m} / h = {beta_m / max(h_exact, 0.01):.4f}")
        print(f"  β_{m+1} / h = {beta_mp1 / max(h_exact, 0.01):.4f}")
        # Ratio β_{m+1}/β_m
        if beta_m > 0:
            print(f"  β_{m+1}/β_{m} = {beta_mp1/beta_m:.4f} = (m+1)/(m-3) = {(m+1)/(m-3):.4f}")

# ═══════════════════════════════════════════════════════════
# 7. THE CRITICAL STRIP ANALOGY
# ═══════════════════════════════════════════════════════════

print("\n" + "="*70)
print("THE CRITICAL STRIP ANALOGY")
print("="*70)

print("""
In Riemann zeta:
  - Critical strip: 0 < Re(s) < 1
  - All nontrivial zeros lie in 0 ≤ Re(s) ≤ 1
  - RH: all nontrivial zeros have Re(s) = 1/2
  - Trivial zeros at s = -2, -4, -6, ...

In tournament Walsh spectrum:
  - Walsh degrees d = 0, 2, 4, ..., n-1
  - Nonzero amplitudes only at EVEN degrees (for odd n)
  - "Trivial zeros" at odd Walsh degrees (vanish identically)
  - The "critical strip" is Walsh degrees 2 ≤ d ≤ n-3
  - THM-J says universality ↔ the degree-1 amplitude isn't too suppressed
  
The analogy deepens:
  - ζ(s) has trivial zeros at negative even integers
  - H_hat[S] = 0 whenever S has an odd-length component
  - Both: "half" of all degrees are trivially zero
  
  - RH controls the error in π(x) ≈ Li(x)  
  - THM-J controls the error in S(T) ≈ 0 (universality)
  - Both: a single spectral condition determines precision
""")

# Compute the "critical line" for Walsh spectrum
print("Walsh amplitude decay (n=7):")
for d in range(0, 7, 2):
    # |H_hat| at degree d with r components, k = d/2
    k = d // 2
    amp = factorial(7 - 2*k) / 2**(7-1)
    # v_2 of the amplitude
    s2 = bin(7-2-d).count('1') if d < 5 else 0
    v2 = -d - s2
    print(f"  degree {d}: |amp| = {amp:.6f}, v_2 = {v2}")

