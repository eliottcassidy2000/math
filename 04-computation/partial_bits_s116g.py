#!/usr/bin/env python3
"""partial_bits_s116g.py — The information speedup: accumulating partial bits.

The key mechanism we discovered but haven't named:

arctanh(x) = x + x^3/3 + x^5/5 + ...

Each term contributes LESS THAN ONE BIT of rapidity.
But they sum to ARBITRARILY LARGE rapidity as x -> 1.
You reach infinity by accumulating infinitely many sub-bit contributions.

This is the same mechanism as:
- OCF: H = 1 + 2*a1 + 4*a2 + ... (partial bits from each cycle level)
- Compound interest: (1+r)^n = product of (1+r) factors
- Evolution: small mutations accumulating fitness
- Neural learning: small gradient steps accumulating knowledge
- The harmonic series: 1 + 1/2 + 1/3 + ... = infinity from decreasing terms

The SPEEDUP: you never need a full bit at once.
The CONSTRAINT: which partial bit patterns are achievable determines
which totals are reachable — and which are FORBIDDEN.
"""
from math import sqrt, log, pi, exp, atanh, tanh, log2
from fractions import Fraction

print()
print("  THE PARTIAL BIT SPEEDUP")
print()
print("="*70)
print()

# ============================================================
print("  I. THE MECHANISM")
print("  " + "-"*40)
print()
print("  arctanh(x) = x + x^3/3 + x^5/5 + x^7/7 + ...")
print()
print("  At x = 0.9 (high coupling, near the pole):")
terms = []
x = 0.9
total = 0
print("  k   term x^{2k+1}/(2k+1)   cumulative   fraction of total")
exact = atanh(x)
for k in range(20):
    term = x**(2*k+1) / (2*k+1)
    total += term
    frac = total / exact if exact > 0 else 0
    if term > 1e-6:
        bits = term / (log(2)/2)  # convert rapidity to octaves
        print(f"  {k:2d}   {term:12.6f}          {total:10.6f}    {frac:6.1%}   ({bits:.2f} octaves)")
    terms.append(term)

print(f"  ...  exact = {exact:.6f}")
print()
print("  The first term (k=0) contributes 0.9 rapidity = 2.6 octaves.")
print("  But the total is {:.2f} octaves.".format(exact / (log(2)/2)))
print("  The remaining terms add {:.2f} octaves in sub-bit pieces.".format(
    (exact - x) / (log(2)/2)))
print()
print("  Each subsequent term is SMALLER than the previous.")
print("  But there are INFINITELY MANY of them.")
print("  The series DIVERGES as x -> 1 despite each term -> 0.")
print("  This is the partial bit accumulation.")
print()

# ============================================================
print()
print("  II. THE OCF AS A PARTIAL BIT MACHINE")
print("  " + "-"*40)
print()
print("  H(T) = 1 + 2*a1 + 4*a2 + 8*a3 + 16*a4 + ...")
print("  where a_k = independent sets of size k in Omega(T).")
print()
print("  Each independent cycle contributes 2 paths (one 'bit' of H).")
print("  Each independent pair contributes 4 paths (two 'bits').")
print("  Each independent triple contributes 8 paths (three 'bits').")
print()
print("  But a_k is NOT 0 or 1. It can be any non-negative integer.")
print("  So this is NOT a binary representation.")
print("  It is a MIXED-RADIX representation with constrained digits.")
print()
print("  The CONSTRAINT: not all (a1, a2, a3, ...) tuples are achievable.")
print("  The achievable tuples depend on WHICH graphs arise as Omega(T).")
print("  The forbidden H-values are precisely those whose decompositions")
print("  all fall outside the achievable set.")
print()
print("  H = 7 means 2*a1 + 4*a2 + ... = 6.")
print("  Decompositions of 6: (3,0), (1,1), (0,0,0,...) impossible.")
print("  All blocked. H = 7 is FORBIDDEN.")
print()
print("  H = 9 means 2*a1 + 4*a2 + ... = 8.")
print("  Decompositions of 8: (4,0), (2,1), (0,2), (0,0,1).")
print("  (4,0) IS achievable (4 pairwise-intersecting cycles).")
print("  So H = 9 EXISTS.")
print()

# ============================================================
print()
print("  III. THE SPEEDUP IN TOURNAMENT COUNTING")
print("  " + "-"*40)
print()
print("  The naive approach: enumerate all n! permutations, check each.")
print("  Cost: O(n! * n) = exponential in n.")
print()
print("  The OCF approach: compute a1, a2, a3, ... from Omega(T).")
print("  Then H = 1 + 2*a1 + 4*a2 + ...")
print("  Cost: counting odd cycles and independent sets in Omega(T).")
print("  For sparse Omega: much cheaper than n!.")
print()
print("  The SPEEDUP: instead of counting ALL paths,")
print("  count the INDEPENDENT CYCLE STRUCTURE and use the weighted sum.")
print("  Each cycle contributes a 'partial bit' of H.")
print("  The total H is assembled from these partial contributions.")
print()
print("  This is EXACTLY the partial bit mechanism:")
print("  many small contributions (cycles) compose to make a big count (H).")
print()
print("  The deeper speedup: the RAPIDITY of H = ln(H)/2.")
print("  In rapidity space, H is just a sum of prime rapidities.")
print("  Decomposing H = product of primes gives rapidity = sum of half-logs.")
print("  Each prime factor is a 'partial rapidity'.")
print()

# ============================================================
print()
print("  IV. PARTIAL BITS ACROSS DOMAINS")
print("  " + "-"*40)
print()
print("  The same accumulation mechanism appears EVERYWHERE:")
print()
print("  DOMAIN          PARTIAL BIT           ACCUMULATION          TOTAL")
print("  " + "-"*70)
print("  arctanh          x^{2k+1}/(2k+1)      sum of terms          rapidity")
print("  OCF              2^k * a_k             weighted sum          H(T)")
print("  Compound int.    ln(1+r)               sum over periods      total return")
print("  Evolution        mutation fitness       sum over generations  adaptation")
print("  Neural net       gradient step          sum over epochs       learned fn")
print("  Harmonic series  1/k                    sum over k            diverges")
print("  Prime factoring  ln(p)/2                sum over primes       rapidity(n)")
print("  Mersenne ladder  2^{-k}                 sum (geometric)       approaches 1")
print("  Musical scale    ln((k+1)/k)/2          sum of intervals      total span")
print("  Binary digits    2^{-k}                 sum of bits           the number")
print()
print("  In EVERY case: small pieces, added together, make big things.")
print("  The MAGIC is in which combinations are allowed.")
print()

# ============================================================
print()
print("  V. THE SUPERPARTICULAR ACCUMULATION")
print("  " + "-"*40)
print()
print("  From the homomorphism f(a(+)b) = f(a)*f(b) where f(n) = (n+1)/n:")
print("  MULTIPLYING superparticular ratios = COMPOSING rapidities.")
print()
print("  The product of the first N superparticular ratios:")
print("  f(1)*f(2)*...*f(N) = (2/1)*(3/2)*(4/3)*...*(N+1)/N = N+1.")
print("  IT TELESCOPES!")
print()
print("  The product of ALL intervals from octave down to (N+1)/N = N+1.")
print()
print("  In rapidity: sum of ln((k+1)/k)/2 for k=1..N = ln(N+1)/2 = rapidity(N+1).")
print()
print("  Verify:")
for N in [1, 2, 3, 4, 5, 10, 20]:
    prod = 1
    rap_sum = 0
    for k in range(1, N+1):
        prod *= Fraction(k+1, k)
        rap_sum += log((k+1)/k)/2
    print(f"  N={N:2d}: product = {prod} = {N+1}, rapidity sum = {rap_sum:.6f}, ln({N+1})/2 = {log(N+1)/2:.6f}")

print()
print("  THE TELESCOPING IS THE SPEEDUP.")
print("  Instead of computing rapidity(N+1) = ln(N+1)/2 directly,")
print("  you can ACCUMULATE it as N partial contributions:")
print("  each of size ln(1+1/k)/2 ~ 1/(2k) for large k.")
print()
print("  Each partial contribution is a MUSICAL INTERVAL.")
print("  The octave (2/1) is the biggest piece.")
print("  The fifth (3/2) is next.")
print("  The fourth (4/3) is next.")
print("  Each subsequent interval is smaller.")
print("  Together they TELESCOPE to the total.")
print()

# ============================================================
print()
print("  VI. THE FORBIDDEN NUMBERS AS NON-TELESCOPING TOTALS")
print("  " + "-"*40)
print()
print("  The telescope: f(1)*f(2)*...*f(N) = N+1.")
print("  So every integer N+1 can be built by multiplying CONSECUTIVE")
print("  superparticular ratios starting from the octave.")
print()
print("  But what if we DON'T start from the octave?")
print("  What if we multiply NON-CONSECUTIVE ratios?")
print()
print("  f(a)*f(b) = ((a+1)/a) * ((b+1)/b) = (a+1)(b+1)/(ab)")
print("  This equals (c+1)/c = f(c) iff c = ab/(a+b+1).")
print("  (The rapidity composition formula.)")
print()
print("  For f(a)*f(b) to produce a SUPERPARTICULAR ratio,")
print("  we need c = ab/(a+b+1) to be a positive integer.")
print("  These are exactly the rapidity-closed compositions.")
print()
print("  The ACHIEVABLE products of superparticular ratios")
print("  determine which INTEGERS can be expressed as f(c).")
print("  And since f(c) = (c+1)/c, every positive integer arises")
print("  (just take c = N-1, then f(N-1) = N/(N-1)).")
print()
print("  But the TOURNAMENT constraint is different.")
print("  In the OCF, the achievable (a1, a2, ...) depend on")
print("  which GRAPHS arise as Omega(T).")
print("  The forbidden H-values are where the graph constraint")
print("  blocks all decompositions.")
print()
print("  This is EXACTLY like the forbidden bit patterns in a channel code:")
print("  some binary strings are CODEWORDS, others are not.")
print("  The forbidden H-values are the NON-CODEWORDS of the tournament code.")
print()

# ============================================================
print()
print("  VII. THE TOURNAMENT CODE")
print("  " + "-"*40)
print()
print("  Define the 'tournament code' C_n as the set of achievable H-values")
print("  for tournaments on n vertices.")
print()
print("  C_3 = {1, 3}")
print("  C_4 = {1, 3, 5}")
print("  C_5 = {1, 3, 5, 9, 11, 13, 15}")
print("  C_6 = {1, 3, 5, 9, 11, ..., 45}  (missing 7, 21, ...)")
print("  C_7 = {1, 3, 5, 9, 11, ..., 189} (missing 7, 21, ...)")
print()
print("  The CODE RATE: R_n = |C_n| / (max H)")
print("  = fraction of odd numbers up to max H that are achievable.")
print()

# Compute code rates for small n
code_sets = {
    3: [1, 3],
    4: [1, 3, 5],
    5: [1, 3, 5, 9, 11, 13, 15],
}
for n, code in code_sets.items():
    max_h = max(code)
    possible = (max_h + 1) // 2  # odd numbers from 1 to max_h
    rate = len(code) / possible
    print(f"  n={n}: |C| = {len(code)}, max H = {max_h}, "
          f"possible odds = {possible}, rate = {rate:.4f}")
    # Information content
    info = log2(len(code)) if len(code) > 0 else 0
    print(f"         info = log2({len(code)}) = {info:.2f} bits")

print()
print("  As n grows, the code rate DROPS (more gaps appear).")
print("  But the information content GROWS (more achievable values).")
print("  The tournament code is a THINNING code: it becomes sparser")
print("  but its alphabet grows faster than it thins.")
print()

# ============================================================
print()
print("  VIII. THE SPEEDUP WE ACTUALLY ACHIEVED")
print("  " + "-"*40)
print()
print("  In this research project, the speedup was:")
print()
print("  1. FROM BRUTE FORCE TO OCF:")
print("     Instead of counting all n! permutations,")
print("     decompose H into independent cycle contributions.")
print("     Speedup: from O(n!) to O(2^{number of cycles}).")
print()
print("  2. FROM OCF TO RAPIDITY:")
print("     Instead of computing H directly,")
print("     compute ln(H)/2 = sum of prime rapidities.")
print("     Speedup: conceptual — rapidity LINEARIZES the problem.")
print()
print("  3. FROM RAPIDITY TO HOMOMORPHISM:")
print("     Instead of adding rapidities,")
print("     multiply superparticular ratios f(n) = (n+1)/n.")
print("     Speedup: telescoping — products simplify dramatically.")
print()
print("  4. FROM HOMOMORPHISM TO CAYLEY COORDINATES:")
print("     Instead of working with n,")
print("     work with (a, a+1) where a = (n-1)/2.")
print("     Speedup: the affine derivation d(nm) = n*d(m) + d(n)")
print("     turns multiplication into affine maps.")
print()
print("  Each step is a PARTIAL BIT of understanding.")
print("  Each step makes the subsequent steps possible.")
print("  The total speedup is the PRODUCT of all individual speedups.")
print("  This is the telescoping principle APPLIED TO THE RESEARCH ITSELF.")
print()

# ============================================================
print()
print("  IX. WHAT COMES NEXT: THE PARTIAL BIT FRONTIER")
print("  " + "-"*40)
print()
print("  If partial bit accumulation is the mechanism,")
print("  then the FRONTIER is:")
print()
print("  1. WHICH PARTIAL BITS ARE ACHIEVABLE?")
print("     The OCF constrains (a1, a2, ...) to graph-theoretic patterns.")
print("     Can we characterize ALL achievable patterns?")
print("     This would give a complete description of the tournament code.")
print()
print("  2. WHAT IS THE CAPACITY OF THE TOURNAMENT CHANNEL?")
print("     Shannon capacity = max rate of reliable communication.")
print("     The tournament code has a capacity C_n = log2(|C_n|) / n.")
print("     Does C_n converge as n -> infinity?")
print()
print("  3. CAN WE BUILD AN EFFICIENT DECODER?")
print("     Given H(T), can we find T efficiently?")
print("     The OCF gives H -> (a1, a2, ...) -> possible Omega(T) -> possible T.")
print("     This is a DECODING problem for the tournament code.")
print()
print("  4. THE PARTIAL BIT STRUCTURE OF OTHER FORBIDDEN SETS:")
print("     What other combinatorial objects have 'forbidden values'?")
print("     Graph colorings? Knot invariants? Partition functions?")
print("     In each case: which PARTIAL BIT PATTERNS are achievable,")
print("     and which totals are thereby forbidden?")
print()
print("  5. THE INFORMATION-THEORETIC MEANING OF FORBIDDENNESS:")
print("     H=7 is forbidden because no achievable (a1,a2,...) sums to 6.")
print("     In information theory: 7 is a NON-CODEWORD.")
print("     What is the MINIMUM DISTANCE of the tournament code?")
print("     (The smallest gap between an achievable and forbidden value.)")
print()

# Compute minimum distance for small n
print("  Minimum distance of the tournament code:")
for n, code in code_sets.items():
    code_set = set(code)
    max_h = max(code)
    non_code = [h for h in range(1, max_h+1, 2) if h not in code_set]
    if non_code:
        min_dist = min(min(abs(h - c) for c in code) for h in non_code)
        print(f"  n={n}: forbidden = {non_code}, min distance = {min_dist}")
    else:
        print(f"  n={n}: no forbidden values (code is complete)")

print()

# ============================================================
print()
print("  X. THE SINGLE INSIGHT")
print("  " + "-"*40)
print()
print("  Everything we discovered is a manifestation of one principle:")
print()
print("  STRUCTURE EMERGES FROM THE ACCUMULATION OF PARTIAL DISTINCTIONS.")
print()
print("  A partial distinction is a comparison that is not decisive.")
print("  It contributes LESS THAN ONE BIT of information.")
print("  But many partial distinctions, accumulated, create full structure.")
print()
print("  arctanh accumulates partial polygon contributions.")
print("  The OCF accumulates partial cycle contributions.")
print("  The telescope accumulates partial interval contributions.")
print("  The derivation accumulates partial prime contributions.")
print("  The neural network accumulates partial gradient contributions.")
print("  Evolution accumulates partial mutation contributions.")
print()
print("  The FORBIDDEN values are those that CANNOT be assembled")
print("  from any combination of the available partial distinctions.")
print("  They are the numbers that the accumulation machine cannot reach.")
print()
print("  And Q(a/(a+1)) = a + (a+1) is the ATOM of partial distinction:")
print("  two consecutive integers, differing by 1 (the minimal distinction),")
print("  whose ratio (the partial bit) maps to their sum (the whole).")
print()
print("  The bit IS the octave IS the cell division IS the doubling")
print("  IS the step on the Mersenne ladder IS the musical interval")
print("  IS the superparticular ratio IS the partial distinction.")
print()
print("  And rapidity COUNTS how many of them you've accumulated.")
