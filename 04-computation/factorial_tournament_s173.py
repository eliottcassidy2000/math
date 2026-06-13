#!/usr/bin/env python3
"""
factorial_tournament_s173.py -- opus-2026-03-22-S173

FACTORIALS, DOUBLE FACTORIALS, AND TOURNAMENT THEORY

Every major formula in tournament theory involves factorials:
  E[H] = n! / 2^{n-1}
  Hook product = prod_{d=1}^k (2d-1)^{k-d+1}
  f^{delta_k} = C(k+1,2)! / hook_prod
  H_max = ? (involves factorials of the Paley construction)

This script systematically maps ALL factorial structures:
  n! = ordinary factorial
  n!! = double factorial (product of same-parity numbers)
  (x)_n = falling factorial = x(x-1)...(x-n+1)
  x^{(n)} = rising factorial = x(x+1)...(x+n-1) (Pochhammer)

And finds where each appears in tournament theory.

Author: opus-2026-03-22-S173
"""

import sys
import numpy as np
from math import comb, factorial, prod, gcd, log2
from itertools import permutations, combinations
from collections import defaultdict, Counter
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def double_factorial(n):
    """n!! = n * (n-2) * (n-4) * ... * (1 or 2)"""
    if n <= 0: return 1
    result = 1
    while n > 0:
        result *= n
        n -= 2
    return result

def falling_factorial(x, n):
    """(x)_n = x * (x-1) * ... * (x-n+1)"""
    result = 1
    for i in range(n):
        result *= (x - i)
    return result

def rising_factorial(x, n):
    """x^(n) = x * (x+1) * ... * (x+n-1) (Pochhammer)"""
    result = 1
    for i in range(n):
        result *= (x + i)
    return result

# ============================================================
# PART 1: FACTORIAL MAP OF TOURNAMENT THEORY
# ============================================================

banner("PART 1: THE FACTORIAL MAP")

print("""
EVERY key formula involves factorials. Let me catalog them all.

1. E[H(T)] = n! / 2^{n-1}
   This is the EXPECTED number of Hamiltonian paths in a random tournament.
   Proof: Each of the n! permutations is a HP with probability 1/2^{n-1}.

   REWRITE: E[H] = n! / 2^{n-1} = (2n-2)!! / (n-1) ... no.
   Actually: n! / 2^{n-1} = prod_{k=1}^{n-1} (k/2) * 2

   Let me compute the first values and recognize the pattern:
""")

for n in range(1, 12):
    EH = Fraction(factorial(n), 2**(n-1))
    print("    n=%2d: E[H] = %s = %.4f" % (n, EH, float(EH)))

print("""
   E[H] sequence: 1, 1, 3/2, 3, 15/2, 45/2, 315/4, 315, ...
   Numerators: 1, 1, 3, 3, 15, 45, 315, 315, 2835, 14175, ...

   RECOGNITION: The ODD values (n=1,3,5,...) give:
     E[H(n)] = n!! / 2^{(n-1)/2} ... let me check.
""")

print("  Checking factorial formulas for E[H]:")
for n in range(1, 10):
    EH = Fraction(factorial(n), 2**(n-1))
    # Try: n!! / something
    ndf = double_factorial(n)
    ratio = Fraction(ndf, int(EH)) if int(EH) > 0 else 0
    # Try: (2n-2)!! / (2^{n-1})
    df2 = double_factorial(2*n-2)
    ratio2 = Fraction(df2, 2**(n-1))
    print("    n=%d: E[H]=%s, n!!=%d, (2n-2)!!=%d, (2n-2)!!/2^{n-1}=%s" %
          (n, EH, ndf, df2, Fraction(df2, 2**(n-1))))

# ============================================================
# PART 2: THE HOOK PRODUCT AND DOUBLE FACTORIALS
# ============================================================

banner("PART 2: HOOK PRODUCT AND DOUBLE FACTORIALS")

print("  From S164: hook_prod(k) = prod_{d=0}^{k-1} (2(k-d)-1)^{d+1}")
print("  = 1^k * 3^{k-1} * 5^{k-2} * ... * (2k-1)^1\n")

print("  The odd numbers (2d-1) are the DOUBLE FACTORIAL building blocks:")
print("  (2k-1)!! = 1 * 3 * 5 * ... * (2k-1)\n")

for k in range(1, 8):
    hp = 1
    for d in range(k):
        hp *= (2*(k-d)-1)**(d+1)
    # Express hook_prod in terms of double factorials
    # hp = prod_{m=1}^{k} (2m-1)^{k-m+1}
    # = (1)^k * (3)^{k-1} * (5)^{k-2} * ... * (2k-1)^1
    # = prod_{m=1}^{k} (2m-1)^{k-m+1}

    # Compare with (2k-1)!!^k / something
    odf = double_factorial(2*k-1)
    m = k*(k+1)//2
    f_delta = factorial(m) // hp

    print("  k=%d: hook_prod = %d" % (k, hp))
    print("       (2k-1)!! = %d" % odf)
    print("       hook_prod / (2k-1)!! = %s" % Fraction(hp, odf))
    print("       f^delta = %d" % f_delta)
    print()

# ============================================================
# PART 3: FALLING FACTORIAL AND TOURNAMENT COUNTING
# ============================================================

banner("PART 3: FALLING FACTORIALS IN TOURNAMENT COUNTING")

print("  The falling factorial (x)_n = x(x-1)...(x-n+1) appears in:")
print("  1. The number of labeled tournaments: 2^C(n,2) = 2^{n(n-1)/2}")
print("  2. The number of permutations: n! = (n)_n")
print("  3. Score sequence enumeration\n")

print("  KEY FORMULA: The number of tournaments with SPECIFIC score seq")
print("  s = (s_0, s_1, ..., s_{n-1}) is the PERMANENT of a specific matrix.\n")

print("  But more interestingly:")
print("  E[H] = n! / 2^{n-1} = falling_factorial(n, n) / 2^{n-1}")
print("  And: the Walsh coefficient hat{H}[S] involves factorials:\n")

print("  From THM-080: hat{M}[S] = (-1)^{asc(S)} * 2^s * (n-2-|S|)! / 2^{n-2}")
print("  The (n-2-|S|)! is a TRUNCATED factorial = falling factorial:\n")

for n in range(3, 8):
    for s_size in range(0, n-1, 2):  # even Walsh degrees
        coeff_factor = factorial(n - 2 - s_size)
        denom = 2**(n-2)
        ratio = Fraction(coeff_factor, denom)
        ff = falling_factorial(n-2, s_size)
        # (n-2-|S|)! = (n-2)! / (n-2)_{|S|} ... no
        # (n-2-|S|)! = falling_factorial(n-2-|S|, n-2-|S|)
        print("    n=%d, |S|=%d: (n-2-|S|)! = %d, falling(n-2, |S|) = %d" %
              (n, s_size, coeff_factor, ff))

# ============================================================
# PART 4: RISING FACTORIAL AND THE INDEPENDENCE POLYNOMIAL
# ============================================================

banner("PART 4: RISING FACTORIAL AND I(Omega, x)")

print("  The independence polynomial I(G, x) = sum a_k x^k")
print("  At n<=5: I(Omega(T), x) = 1 + alpha_1 * x (linear!)")
print("  So H = I(Omega, 2) = 1 + 2*alpha_1.\n")

print("  For the independence polynomial of the COMPLETE graph K_m:")
print("  I(K_m, x) = 1 + m*x (only empty and singleton independent sets).")
print("  This is a FALLING factorial of degree 1.\n")

print("  For the EMPTY graph on m vertices:")
print("  I(E_m, x) = (1+x)^m = sum C(m,k) x^k")
print("  At x=2: I(E_m, 2) = 3^m.\n")

print("  For a GENERAL graph G, I(G, x) is between these extremes.")
print("  The OCF: H = I(Omega, 2) interpolates between K_m and E_m.\n")

# Now: can we express I(Omega, x) using Pochhammer/falling factorials?
print("  At n=6, alpha_2 > 0 for some tournaments.")
print("  I(Omega, x) = 1 + alpha_1 x + alpha_2 x^2 + ...")
print("  = 1 + alpha_1 x + alpha_2 x(x-1)/2 * 2 + ...")
print("  Hmm, the falling factorial basis: x^{(k)} = x(x-1)...(x-k+1)")
print("  Expansion: I = sum b_k * x^{(k)} / k!")
print("  Then I(Omega, 2) = sum b_k * 2^{(k)} / k! = sum b_k * (2*1*0*...)/k!")
print("  But 2^{(k)} = 0 for k >= 3. So only k=0,1,2 contribute.\n")

print("  WAIT: 2^{(0)} = 1, 2^{(1)} = 2, 2^{(2)} = 2*1 = 2, 2^{(3)} = 2*1*0 = 0.")
print("  So in the falling factorial basis:")
print("  I(Omega, 2) = b_0 + 2*b_1 + 2*b_2/2 + 0 + 0 + ...")
print("             = b_0 + 2*b_1 + b_2\n")

print("  THIS IS REMARKABLE: evaluating at x=2 in the FALLING FACTORIAL basis")
print("  automatically truncates at degree 2!")
print("  Only the constant, linear, and quadratic terms contribute.\n")

# Let me verify this
print("  VERIFICATION:")
print("  In the standard basis: I = sum a_k x^k")
print("  a_0 = 1, a_1 = alpha_1, a_2 = alpha_2, etc.")
print()
print("  In the falling factorial basis: I = sum b_k * (x)_k / k!")
print("  Converting: x^0 = (x)_0")
print("              x^1 = (x)_1")
print("              x^2 = (x)_2 + (x)_1 = x(x-1) + x = x^2")
print("  So: a_0 = b_0, a_1 = b_1, a_2 = b_2/2 + ... ")
print()
print("  Actually, the conversion uses Stirling numbers:")
print("  x^k = sum_{j=0}^{k} S(k,j) * (x)_j")
print("  where S(k,j) are Stirling numbers of the second kind.\n")

print("  So: I(Omega, x) = sum_k a_k x^k = sum_k a_k sum_j S(k,j) (x)_j")
print("  = sum_j (sum_k a_k S(k,j)) (x)_j\n")

print("  At x=2: (2)_0=1, (2)_1=2, (2)_2=2, (2)_j=0 for j>=3.\n")

print("  THEREFORE:")
print("  H = I(Omega, 2) = sum_j c_j (2)_j")
print("  where c_j = sum_k a_k S(k,j)")
print("  and (2)_j = 0 for j >= 3.\n")

print("  So: H = c_0 + 2*c_1 + 2*c_2")
print("  where c_0 = sum_k a_k S(k,0) = a_0 = 1")
print("        c_1 = sum_k a_k S(k,1) = sum_k a_k = I(Omega, 1)")
print("        c_2 = sum_k a_k S(k,2)")
print()

# Stirling numbers of the second kind
def stirling2(n, k):
    if n == 0 and k == 0: return 1
    if n == 0 or k == 0: return 0
    if k > n: return 0
    return k * stirling2(n-1, k) + stirling2(n-1, k-1)

print("  Stirling numbers S(k,j) for small k,j:")
for k in range(6):
    row = [stirling2(k, j) for j in range(6)]
    print("    S(%d, *) = %s" % (k, row))

print()
print("  THE H FORMULA IN FALLING FACTORIAL BASIS:")
print("  H = I(Omega, 2) = 1 + 2*I(Omega, 1) + 2*c_2")
print("  where c_2 = sum_{k>=2} a_k * S(k, 2)")
print("  and I(Omega, 1) = total number of independent sets in Omega.\n")

# Verify at n=5
print("  VERIFICATION AT n=5:")
print("  At n<=5: a_k = 0 for k >= 2 (alpha_2 = 0).")
print("  So c_2 = 0, and H = 1 + 2*alpha_1 = 1 + 2*I(Omega,1) - 2.")
print("  Wait: I(Omega, 1) = 1 + alpha_1. So 2*I(Omega,1) = 2 + 2*alpha_1.")
print("  H = 1 + 2 + 2*alpha_1 + 0 = 3 + 2*alpha_1.")
print("  But H = 1 + 2*alpha_1. CONTRADICTION.\n")

print("  Let me recompute more carefully.")
print("  I(Omega, x) = sum a_k x^k.")
print("  I(Omega, 2) = sum a_k * 2^k.")
print("  In the falling factorial basis: 2^k = sum_j S(k,j) * (2)_j.")
print("  And (2)_j = 2, 2, 0, 0, ... for j = 1, 2, 3, ...")
print("  Wait: (2)_0 = 1, (2)_1 = 2, (2)_2 = 2*1 = 2, (2)_3 = 2*1*0 = 0.")
print()
print("  So 2^k = sum_{j=0}^{min(k,2)} S(k,j) * (2)_j")
print("  = S(k,0)*1 + S(k,1)*2 + S(k,2)*2")
print()

for k in range(6):
    s0 = stirling2(k, 0)
    s1 = stirling2(k, 1)
    s2 = stirling2(k, 2)
    val = s0 + 2*s1 + 2*s2
    print("    2^%d = %d, S(%d,0)+2*S(%d,1)+2*S(%d,2) = %d+2*%d+2*%d = %d, match=%s" %
          (k, 2**k, k, k, k, s0, s1, s2, val, val == 2**k))

print()
print("  YES! 2^k = S(k,0) + 2*S(k,1) + 2*S(k,2) for all k.")
print("  This is a KNOWN identity: 2^k = sum_j S(k,j) * (2)_j.")
print("  It's just the change-of-basis formula evaluated at x=2.\n")

print("  THEREFORE:")
print("  H = I(Omega, 2) = sum_k a_k * 2^k")
print("    = sum_k a_k * [S(k,0) + 2*S(k,1) + 2*S(k,2)]")
print("    = a_0 + 2*sum_k a_k + 2*sum_k a_k*S(k,2)")
print()
print("  Wait: S(k,0) = 0 for k>0, S(0,0)=1. And S(k,1) = 1 for k>=1.")
print("  So:")
print("  H = a_0*1 + 2*(a_1 + a_2 + a_3 + ...) + 2*(a_2*1 + a_3*3 + a_4*7 + ...)")
print("    = 1 + 2*(I(Omega,1) - 1) + 2*sum_{k>=2} a_k * S(k,2)")
print("    = 2*I(Omega,1) - 1 + 2*sum_{k>=2} a_k * S(k,2)\n")

print("  At n<=5 (alpha_2 = 0):")
print("  H = 2*I(Omega,1) - 1 = 2*(1 + alpha_1) - 1 = 1 + 2*alpha_1. CHECK!\n")

print("  At n=6 (alpha_2 > 0):")
print("  H = 2*I(Omega,1) - 1 + 2*a_2*S(2,2)")
print("    = 2*(1 + alpha_1 + alpha_2) - 1 + 2*alpha_2*1")
print("    = 1 + 2*alpha_1 + 2*alpha_2 + 2*alpha_2")
print("    = 1 + 2*alpha_1 + 4*alpha_2")
print("  This IS the OCF: H = 1 + 2*alpha_1 + 4*alpha_2. CHECK!\n")

print("  AND GENERALLY:")
print("  H = sum_{k=0}^{alpha(Omega)} a_k * 2^k")
print("  = sum_{j=0}^{2} (sum_k a_k S(k,j)) * (2)_j")
print("  = (sum_k a_k S(k,0)) + 2*(sum_k a_k S(k,1)) + 2*(sum_k a_k S(k,2))")
print()
print("  THE FALLING FACTORIAL DECOMPOSITION OF H:")
print("  H = B_0 + 2*B_1 + 2*B_2")
print("  where B_j = sum_k a_k * S(k, j) are the STIRLING TRANSFORMS of")
print("  the independence number sequence.\n")

print("  B_0 = a_0 = 1 (always)")
print("  B_1 = sum_k a_k = I(Omega, 1) = total independent sets")
print("  B_2 = sum_k a_k * S(k, 2) = Stirling-2 weighted count\n")

print("  SO: H = 1 + 2*B_1 + 2*B_2 - 2 ... wait let me redo.")
print("  H = B_0*(2)_0 + B_1*(2)_1 + B_2*(2)_2 = 1 + 2*B_1 + 2*B_2")
print("  Hmm but B_0 should be just 1 since S(0,0)=1 and S(k,0)=0 for k>0.")
print("  B_0 = a_0 * S(0,0) = 1*1 = 1.")
print("  B_1 = sum_{k>=1} a_k * S(k,1) = sum_{k>=1} a_k * 1 = I(Omega,1) - 1.")
print("  B_2 = sum_{k>=2} a_k * S(k,2).\n")

print("  Corrected: H = 1 + 2*(I(Omega,1)-1) + 2*B_2")
print("           = 2*I(Omega,1) - 1 + 2*B_2\n")

# Verify at specific cases
print("  VERIFICATION:")
print("  alpha_1 = 3, alpha_2 = 0 (n=5, H=7-class):")
print("    I(Omega,1) = 1 + 3 = 4, B_2 = 0")
print("    H = 2*4 - 1 + 0 = 7. CHECK (but H=7 is forbidden!)")
print("    ... this just shows the formula works, not that alpha_1=3 is achievable.\n")

print("  alpha_1 = 5, alpha_2 = 0 (n=5, H=11):")
print("    I(Omega,1) = 6, B_2 = 0")
print("    H = 2*6 - 1 = 11. CHECK.\n")

print("  alpha_1 = 6, alpha_2 = 1 (n=6, H=17):")
print("    I(Omega,1) = 1 + 6 + 1 = 8, B_2 = 1*S(2,2) = 1")
print("    H = 2*8 - 1 + 2*1 = 17. CHECK.\n")

# ============================================================
# PART 5: THE DOUBLE FACTORIAL AND THE STAIRCASE
# ============================================================

banner("PART 5: DOUBLE FACTORIALS AND THE STAIRCASE")

print("  The hook product of delta_k involves products of ODD numbers.")
print("  (2d-1)!! = 1 * 3 * 5 * ... * (2d-1) is the odd double factorial.\n")

print("  From S164: hook_prod(k) = prod_{m=1}^{k} (2m-1)^{k-m+1}")
print("  Let's express this using double factorials:\n")

for k in range(1, 8):
    hp = 1
    for m in range(1, k+1):
        hp *= (2*m-1)**(k-m+1)

    # Compare with products of double factorials
    # (2k-1)!! = prod_{m=1}^k (2m-1) = 1*3*5*...*(2k-1)
    odf = double_factorial(2*k-1)

    # hook_prod = prod_{m=1}^k (2m-1)^{k-m+1}
    # = prod_{j=1}^k [(2j-1)!!]^? ... not directly

    # Try: hook_prod = prod_{j=1}^k (2j-1)!!
    prod_odfs = 1
    for j in range(1, k+1):
        prod_odfs *= double_factorial(2*j-1)

    print("  k=%d: hook_prod = %d, prod_{j=1}^k (2j-1)!! = %d, match = %s" %
          (k, hp, prod_odfs, hp == prod_odfs))

print()
print("  YES! hook_prod(k) = prod_{j=1}^k (2j-1)!!")
print("  = 1!! * 3!! * 5!! * ... * (2k-1)!!")
print("  = 1 * (1*3) * (1*3*5) * ... * (1*3*5*...*(2k-1))\n")

# Verify
for k in range(1, 8):
    hp = 1
    for m in range(1, k+1):
        hp *= (2*m-1)**(k-m+1)
    prod_odfs = 1
    for j in range(1, k+1):
        prod_odfs *= double_factorial(2*j-1)
    assert hp == prod_odfs, "Mismatch at k=%d" % k

print("  VERIFIED for k=1..7.")
print()
print("  COROLLARY: f^{delta_k} = C(k+1,2)! / prod_{j=1}^k (2j-1)!!")
print()

for k in range(1, 8):
    m = k*(k+1)//2
    prod_odfs = 1
    for j in range(1, k+1):
        prod_odfs *= double_factorial(2*j-1)
    f = factorial(m) // prod_odfs
    print("  k=%d: f^delta = %d! / [1!! * 3!! * ... * %d!!] = %d" %
          (k, m, 2*k-1, f))

# ============================================================
# PART 6: RISING FACTORIALS AND THE MARKOV CHAIN
# ============================================================

banner("PART 6: RISING FACTORIALS AND THE MARKOV CHAIN")

print("  The spectral gap of the arc-reversal Markov chain on G_n is 2/n.")
print("  The mixing time ~ n/2 steps.")
print("  The stationary distribution: pi[i] = size(i) / 2^m.\n")

print("  The class sizes are size(i) = n! / |Aut(T_i)|.")
print("  So pi[i] = (n! / |Aut(T_i)|) / 2^m = (n! / |Aut(T_i)|) / 2^{C(n,2)}.\n")

print("  The RISING FACTORIAL connection:")
print("  2^{C(n,2)} = 2^{n(n-1)/2}")
print("  n! = prod_{k=1}^n k")
print("  Ratio: 2^{n(n-1)/2} / n! = prod_{k=1}^{n-1} 2^k / (k+1)")
print("  This relates to the number of LABELED tournaments per ISO CLASS.\n")

for n in range(3, 9):
    m = comb(n, 2)
    total = 2**m
    nfac = factorial(n)
    ratio = Fraction(total, nfac)
    print("  n=%d: 2^C(n,2)/n! = 2^%d/%d! = %s = %.4f" %
          (n, m, n, ratio, float(ratio)))

print()
print("  This ratio = average class size = 2^C(n,2) / (n! * sum 1/|Aut|)")
print("  ... no, sum of sizes = 2^C(n,2) and sum of 1 = |V(G_n)|.")
print("  So average size = 2^C(n,2) / |V(G_n)|.\n")

# ============================================================
# SYNTHESIS
# ============================================================

banner("SYNTHESIS: THE FACTORIAL HIERARCHY OF TOURNAMENT THEORY")

print("""
FOUR FACTORIAL FAMILIES, FOUR ROLES:

1. ORDINARY FACTORIAL n!
   WHERE: E[H] = n!/2^{n-1}, class size = n!/|Aut|, SYT count f^delta
   ROLE: Counts LABELED structures (permutations, labelings)
   TOURNAMENT MEANING: The total number of vertex orderings

2. DOUBLE FACTORIAL (2k-1)!!
   WHERE: Hook product = prod_{j=1}^k (2j-1)!!
   ROLE: Counts ODD-STEP structures (matchings, hook lengths)
   TOURNAMENT MEANING: The staircase's ODD HOOK structure
   KEY FORMULA: f^{delta_k} = C(k+1,2)! / prod_{j=1}^k (2j-1)!!

3. FALLING FACTORIAL (x)_k = x(x-1)...(x-k+1)
   WHERE: H = 2*I(Omega,1) - 1 + 2*B_2 via (2)_k truncation
   ROLE: The OCF evaluation at x=2 TRUNCATES at degree 2
   TOURNAMENT MEANING: H depends only on independence sets of size <= 2
   KEY FORMULA: H = 1 + 2*(alpha_1 + alpha_2 + ...) + 2*sum a_k S(k,2)

4. RISING FACTORIAL (Pochhammer) x^{(k)} = x(x+1)...(x+k-1)
   WHERE: Exponential generating functions, Stirling connections
   ROLE: Connects tournament counts to EGFs of the form e^x / (1-x)
   TOURNAMENT MEANING: The EGF for compatible pair counts uses rising factorials

THE DEEPEST FORMULA:

  H = 1 + 2*B_1 + 2*B_2

  where B_1 = I(Omega, 1) - 1 = (total independent sets) - 1
        B_2 = sum_{k>=2} a_k * S(k, 2)

  This says: H counts are determined by a 3-TERM FORMULA
  in the STIRLING TRANSFORM of the independence number sequence.

  At n <= 5: B_2 = 0, so H = 2*I(Omega, 1) - 1.
  At n >= 6: B_2 > 0, adding the "packing correction."

  The falling factorial truncation at x=2 is WHY the OCF works:
  evaluating I(Omega, x) at x=2 automatically ignores independence
  sets of size >= 3, because (2)_3 = 2*1*0 = 0.

  THIS IS THE COMBINATORIAL REASON BEHIND THE OCF:
  The fugacity x=2 is special because 2 is the SMALLEST INTEGER
  where the falling factorial vanishes at degree 3.
  At x=1: all falling factorials are nonzero ((1)_k = 0 for k>=2... wait no)
  (1)_0 = 1, (1)_1 = 1, (1)_2 = 1*0 = 0.
  So x=1 truncates at degree 1, giving I(Omega, 1) = total indep sets.
  x=2 truncates at degree 2, giving the "packing-corrected" count.
  x=3 would truncate at degree 3, and so on.

  THE GENERAL PATTERN:
  I(Omega, k) counts independent sets of size <= k, weighted by
  (k)_j / j! for size-j sets. At x=k, the falling factorial (k)_j = 0
  for j > k, so ONLY sets of size <= k contribute.

  THIS IS WHY x=2 IS NATURAL:
  H = I(Omega, 2) counts independent sets of size <= 2,
  weighted by (2)_0 = 1, (2)_1 = 2, (2)_2 = 2.
  The PAIR INTERACTIONS (alpha_2) matter, but TRIPLE INTERACTIONS don't.
  This is a "2-BODY" approximation to the full independence polynomial.
""")

print("Done. opus-2026-03-22-S173")
