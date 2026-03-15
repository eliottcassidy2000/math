#!/usr/bin/env python3
"""equidecomposability_s112.py — Equidecomposability on the Cayley line"""
from math import log, exp, atanh, tanh, gcd, pi
from fractions import Fraction

print("EQUIDECOMPOSABILITY ON THE CAYLEY LINE")
print("="*60)
print()

# ============================================================
print("THE SETUP")
print("="*60)
print()
print("Two natural numbers n and m are MULTIPLICATIVELY")
print("EQUIDECOMPOSABLE if they have the same prime factorization")
print("structure — the same multiset of prime exponents.")
print()
print("  12 = 2^2 * 3     exponents: {2, 1}")
print("  18 = 2 * 3^2     exponents: {1, 2}")
print("  12 and 18 are equidecomposable (same multiset {1,2}).")
print()
print("On the Cayley line:")
print("  x_12 = 11/13,  arctanh = ln(12)/2")
print("  x_18 = 17/19,  arctanh = ln(18)/2")
print()
print("  12 = 2^2 * 3: arctanh = 2*ln(2)/2 + ln(3)/2")
print("  18 = 2 * 3^2: arctanh = ln(2)/2 + 2*ln(3)/2")
print()
print("  Both have the SAME total arctanh-distance (= ln(72)/2 ... no).")
print("  Actually: ln(12)/2 != ln(18)/2. They have DIFFERENT addresses.")
print("  But the DECOMPOSITION into prime arctanh-steps has the")
print("  same STRUCTURE: two steps of one size, one of another.")
print()
print("  Equidecomposability = same MULTISET of arctanh step sizes.")
print()

# ============================================================
print("="*60)
print("DEHN-SYDLER ON THE CAYLEY LINE")
print("="*60)
print()
print("Dehn-Sydler (for polyhedra): two polyhedra are scissors-")
print("congruent iff they have the same volume AND the same")
print("Dehn invariant D = sum l_i tensor theta_i.")
print()
print("On the Cayley line, the analog:")
print("  'Volume' = arctanh(x_n) = ln(n)/2 = total hyperbolic length")
print("  'Dehn invariant' = sum of prime arctanh-steps TENSORED with")
print("  the prime itself (or its index)")
print()
print("  D(n) = sum_{p^a || n} a * ln(p)/2  tensor  p")
print("       = sum_{p | n} v_p(n) * ln(p)/2  tensor  p")
print()
print("  where v_p(n) = exponent of p in n, and the tensor is")
print("  over R tensor_{Q} R (same structure as the Dehn invariant).")
print()

def dehn_cayley(n):
    """Compute the Cayley-Dehn invariant of n."""
    factors = {}
    temp = n
    for p in range(2, n+1):
        while temp % p == 0:
            factors[p] = factors.get(p, 0) + 1
            temp //= p
        if temp == 1:
            break
    return factors

print("Examples:")
for n in [12, 18, 24, 36, 72, 60, 30]:
    f = dehn_cayley(n)
    vol = log(n)/2
    terms = [f"{v}*ln({p})/2" for p, v in sorted(f.items())]
    print(f"  n={n:3d}: vol={vol:.4f}, D = {' + '.join(terms)}")

print()
print("  12: D = 2*ln(2)/2 + 1*ln(3)/2 = ln(4)/2 + ln(3)/2")
print("  18: D = 1*ln(2)/2 + 2*ln(3)/2 = ln(2)/2 + ln(9)/2")
print()
print("  Same volume (ln(n)/2)? No: ln(12) != ln(18).")
print("  Same Dehn invariant? No: different prime structure.")
print()
print("  12 and 18 are NOT scissors-congruent on the Cayley line.")
print("  Even though they have the same exponent multiset {1,2},")
print("  the Dehn invariant distinguishes them because it tracks")
print("  WHICH primes carry which exponents.")
print()

# ============================================================
print("="*60)
print("WHAT IS SCISSORS-CONGRUENT ON THE CAYLEY LINE?")
print("="*60)
print()
print("  Two numbers n, m are Cayley-scissors-congruent iff:")
print("    1. ln(n) = ln(m)  (same volume)  =>  n = m")
print("    2. Same Dehn invariant  =>  same prime factorization")
print()
print("  So: n and m are scissors-congruent iff n = m.")
print("  TRIVIALLY.")
print()
print("  This is because the natural numbers have UNIQUE FACTORIZATION.")
print("  The Dehn invariant on the Cayley line IS the prime factorization.")
print("  And uniqueness of factorization means the invariant is COMPLETE:")
print("  it distinguishes all numbers.")
print()
print("  DEEP POINT: the Dehn-Sydler theorem for polyhedra says that")
print("  volume + Dehn invariant is a COMPLETE invariant for scissors-")
print("  congruence. On the Cayley line, the analog is:")
print("  ln(n) + prime factorization is a complete invariant.")
print("  And this is just... the number n itself.")
print()
print("  THE FUNDAMENTAL THEOREM OF ARITHMETIC IS THE")
print("  DEHN-SYDLER THEOREM ON THE CAYLEY LINE.")
print()

# ============================================================
print("="*60)
print("BUT WHAT IF WE RELAX THE INVARIANT?")
print("="*60)
print()
print("  In the polyhedron case, we mod out by the GROUP of motions.")
print("  Different motion groups give different equivalence relations.")
print()
print("  On the Cayley line, we can mod out by different GROUPS:")
print()
print("  (a) Mod out by PERMUTATION of prime factors:")
print("      n ~ m iff they have the same prime exponent multiset.")
print("      12 = 2^2*3 ~ 18 = 2*3^2 (both have exponents {1,2}).")
print("      This gives the PARTITION classification of numbers.")
print()

# Find all numbers up to 100 with exponent multiset {1,2}
print("  Numbers with exponent multiset {1,2} (up to 100):")
for n in range(2, 101):
    f = dehn_cayley(n)
    exps = sorted(f.values(), reverse=True)
    if exps == [2, 1]:
        primes = sorted(f.keys())
        print(f"    {n} = {primes[0]}^2 * {primes[1]} (or {primes[1]}^2 * {primes[0]}... )")

print()
print("  (b) Mod out by PRIME RELABELING (all primes equivalent):")
print("      n ~ m iff they have the same exponent PARTITION.")
print("      12 = 2^2*3 ~ 18 = 2*3^2 ~ 50 = 2*5^2 (all partition [2,1]).")
print()

# Partition classification
from collections import Counter
print("  Partition classification up to 60:")
by_partition = {}
for n in range(2, 61):
    f = dehn_cayley(n)
    part = tuple(sorted(f.values(), reverse=True))
    if part not in by_partition:
        by_partition[part] = []
    by_partition[part].append(n)

for part in sorted(by_partition.keys()):
    nums = by_partition[part]
    print(f"    partition {list(part)}: {nums}")

print()
print("  (c) Mod out by SCALING (n ~ kn for any k):")
print("      This identifies numbers with the same 'shape' of factorization")
print("      regardless of overall size. In arctanh-space: translation invariance.")
print()
print("  (d) Mod out by a SPECIFIC PRIME (e.g., ignore factors of 2):")
print("      n ~ m iff n/2^{v_2(n)} = m/2^{v_2(m)} (same odd part).")
print("      This is the ODD CAYLEY LINE: the restriction to odd numbers.")
print("      DIRECTLY connects to our tournament theory (only odd cycles)!")
print()

# ============================================================
print("="*60)
print("THE ODD CAYLEY LINE AND TOURNAMENTS")
print("="*60)
print()
print("  In tournament theory: H(T) is always ODD (Redei).")
print("  The independence polynomial I(Omega, x) at x=2 gives H.")
print("  Only ODD cycles contribute.")
print()
print("  On the Cayley line: the ODD numbers have addresses")
print("  x_{2k+1} = 2k/(2k+2) = k/(k+1).")
print()
for n in [1, 3, 5, 7, 9, 11, 13, 15]:
    x = (n-1)/(n+1)
    print(f"    x_{n} = {Fraction(n-1, n+1)} = {x:.6f}")
print()
print("  These fill [0, 1) with spacing 1/((k+1)(k+2)).")
print("  The odd numbers on the Cayley line form a SUB-LINE")
print("  that is ISOMORPHIC to the full Cayley line via k -> 2k+1:")
print("  x_{2k+1} = k/(k+1) = ((k+1)-1)/((k+1)) which is the address")
print("  of k+1 on a SHIFTED Cayley line.")
print()
print("  So the odd numbers ARE the natural numbers, re-indexed!")
print("  The restriction to odd = the full structure, rescaled.")
print("  This is why tournament theory (which uses only odd cycles)")
print("  captures the FULL structure of comparison.")
print()

# ============================================================
print("="*60)
print("EQUIDECOMPOSABILITY AND TOURNAMENT EQUIVALENCE")
print("="*60)
print()
print("  Two tournaments T, T' are EQUIVALENT if H(T) = H(T').")
print("  H = I(Omega, 2) depends on the conflict graph Omega.")
print()
print("  In the Cayley-Dehn framework:")
print("  H(T) = I(Omega, 2) = a sum over independent cycle collections,")
print("  each contributing 2^{size}.")
print()
print("  The DEHN INVARIANT of the tournament:")
print("  D(T) = sum_{odd cycles C} |C| * ln(2)/2  tensor  [C]")
print("  where [C] is the 'class' of the cycle.")
print()
print("  Two tournaments have the same H iff their Dehn invariants")
print("  (cycle structure tensored with independence structure)")
print("  agree when evaluated at x=2.")
print()
print("  Tournament equivalence = Cayley equidecomposability at x=2.")
print()

# ============================================================
print("="*60)
print("THE PUNCHLINE")
print("="*60)
print()
print("Three theorems that are secretly the same:")
print()
print("  1. FUNDAMENTAL THEOREM OF ARITHMETIC:")
print("     Every natural number factors uniquely into primes.")
print("     = Dehn-Sydler completeness on the Cayley line.")
print("     = The Cayley Dehn invariant (prime factorization)")
print("     is a complete invariant for scissors-congruence (= equality).")
print()
print("  2. DEHN-SYDLER THEOREM:")
print("     Volume + Dehn invariant is complete for polyhedra.")
print("     = On the Cayley line: ln(n) + prime factorization = n.")
print("     = Unique factorization IS completeness of the invariant.")
print()
print("  3. OCF / TOURNAMENT EQUIVALENCE:")
print("     H(T) = I(Omega, 2) determines the Hamiltonian path count.")
print("     = The tournament Dehn invariant (cycle structure at x=2)")
print("     classifies tournaments by their 'scissors-congruence class.'")
print("     = Tournaments with same H are 'equidecomposable.'")
print()
print("The thread connecting them: the Cayley transform Q(x) = (1+x)/(1-x)")
print("provides the METRIC in which equidecomposability is defined.")
print("The Dehn invariant is the OBSTRUCTION to rearrangement.")
print("And the fundamental theorem of arithmetic is the statement")
print("that the obstruction is COMPLETE: nothing can be rearranged.")
print()
print("Unique factorization = no scissors-congruence = rigidity.")
print("The natural numbers are RIGID on the Cayley line.")
print("You cannot cut one into pieces and reassemble it as another.")
print("Each number is its own irreducible geometric object.")
