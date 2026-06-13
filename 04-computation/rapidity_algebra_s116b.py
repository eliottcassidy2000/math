#!/usr/bin/env python3
"""rapidity_algebra_s116b.py — The algebra of rapidity-closed compositions

When does 1/(2a+1) (+) 1/(2b+1) = 1/(2c+1)?
Answer: c = ab/(a+b+1), so (a+b+1) | ab.

This is a NUMBER-THEORETIC condition with deep structure.
"""
from math import sqrt, log, gcd
from fractions import Fraction

print("THE RAPIDITY ALGEBRA: WHEN MUSICAL INTERVALS CLOSE")
print("="*60)
print()
print("  Condition: c = ab/(a+b+1) is a positive integer.")
print("  Equivalently: (a+b+1) | ab.")
print()

# ============================================================
print("STRUCTURE THEOREM")
print("-"*40)
print()
print("  Let d = gcd(a, a+b+1). Then d | a and d | (a+b+1).")
print("  Since d | a, d | (a+b+1-a) = b+1. So d | gcd(a, b+1).")
print("  Similarly let e = gcd(b, a+b+1). Then e | gcd(b, a+1).")
print()
print("  (a+b+1) | ab iff (a+b+1)/gcd(a,a+b+1) | b")
print("  iff (a+b+1)/gcd(a,b+1) | b")
print()
print("  So the condition depends on gcd(a, b+1) and gcd(b, a+1).")
print()

# Enumerate ALL solutions with a <= b <= 50
solutions = []
for a in range(1, 51):
    for b in range(a, 200):
        if (a*b) % (a+b+1) == 0:
            c = (a*b) // (a+b+1)
            solutions.append((a, b, c))

print(f"  Found {len(solutions)} solutions with a <= 50, b <= 200:")
print()
print("  a    b    c    1/(2a+1) (+) 1/(2b+1) = 1/(2c+1)    gcd(a,b+1) gcd(b,a+1)")
print("  " + "-"*75)
for a, b, c in solutions[:30]:
    g1 = gcd(a, b+1)
    g2 = gcd(b, a+1)
    print(f"  {a:3d}  {b:3d}  {c:3d}    1/{2*a+1:3d} (+) 1/{2*b+1:3d} = 1/{2*c+1:3d}           {g1:3d}       {g2:3d}")

print(f"  ... ({len(solutions)} total)")
print()

# ============================================================
print("PATTERN ANALYSIS")
print("-"*40)
print()

# Group by c
from collections import defaultdict
by_c = defaultdict(list)
for a, b, c in solutions:
    by_c[c].append((a, b))

print("  Grouped by output c:")
for c in sorted(by_c.keys())[:15]:
    pairs = by_c[c]
    print(f"  c={c:3d} (1/{2*c+1:3d}): {len(pairs)} decompositions")
    for a, b in pairs[:5]:
        print(f"         {a:3d} (+) {b:3d}  i.e. 1/{2*a+1} (+) 1/{2*b+1}")
    if len(pairs) > 5:
        print(f"         ... {len(pairs)-5} more")
print()

# ============================================================
print("FAMILIES OF SOLUTIONS")
print("-"*40)
print()

# Family 1: a = n(n+1)-1, b = n(n+2), c = n^2-1
# Or: consecutive pairs
print("  FAMILY 1: CONSECUTIVE (a, b) = (k, k+1) -> c = k(k+1)/(2k+2) = k/2")
print("  Only works when k is even: (k, k+1, k/2)")
print()
for k in range(2, 20, 2):
    a, b = k, k+1
    if (a*b) % (a+b+1) == 0:
        c = (a*b) // (a+b+1)
        print(f"    ({a}, {b}) -> c = {c}")

print()
print("  FAMILY 2: (a, b) where b = a^2 + a - 1")
print("  c = a(a^2+a-1)/(a^2+a) = a - 1 + something... checking")
for a in range(2, 15):
    b = a*a + a - 1
    if (a*b) % (a+b+1) == 0:
        c = (a*b) // (a+b+1)
        print(f"    ({a}, {b}) -> c = {c}")

print()
# Check: when b = (a+1)(2c+1) - a - 1 = 2c(a+1) + a
# c = ab/(a+b+1). If we fix a, then b = c(a+b+1)/a, so b(a-c) = c(a+1), b = c(a+1)/(a-c).
# For b > 0: need a > c.
print("  PARAMETRIC FAMILY: Fix a, then b = c(a+1)/(a-c) for c < a.")
print("  b is integer iff (a-c) | c(a+1).")
print()

for a in range(2, 8):
    print(f"  a = {a}:")
    for c in range(1, a):
        if (c*(a+1)) % (a-c) == 0:
            b = (c*(a+1)) // (a-c)
            if b >= a:
                print(f"    c={c}: b = {c}*{a+1}/{a-c} = {b}  -> 1/{2*a+1} (+) 1/{2*b+1} = 1/{2*c+1}")

print()

# ============================================================
print("THE PALEY PRIME COMPOSITIONS")
print("-"*40)
print()
print("  Which compositions involve Paley primes {3, 7, 11}?")
print("  In our notation: 3 -> a=1 (1/3), 7 -> a=3 (1/7), 11 -> a=5 (1/11)")
print()
paley_a = {1: 3, 3: 7, 5: 11}  # a -> 2a+1
for a, b, c in solutions:
    if (2*a+1) in [3, 5, 7, 11, 13, 19, 23] or (2*b+1) in [3, 5, 7, 11, 13, 19, 23] or (2*c+1) in [3, 5, 7, 11, 13, 19, 23]:
        names = []
        for val in [2*a+1, 2*b+1, 2*c+1]:
            if val in [3, 7, 11]:
                names.append(f"PALEY({val})")
            elif val in [5, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
                names.append(f"prime({val})")
            else:
                names.append(str(val))
        print(f"  1/{2*a+1} (+) 1/{2*b+1} = 1/{2*c+1}   [{', '.join(names)}]")

print()

# ============================================================
print("="*60)
print("THE RAPIDITY MONOID")
print("="*60)
print()
print("  The set {1/(2n+1) : n >= 1} under (+) is NOT closed.")
print("  But the SUBSET of compositions that ARE closed forms a structure.")
print()
print("  The closure condition c = ab/(a+b+1) being integer means")
print("  the rapidity-closed triples form a PARTIAL BINARY OPERATION.")
print()
print("  Is it associative? Check: (1/5 (+) 1/7) (+) 1/9 = 1/3 (+) 1/9 = 3/7.")
print("  And: 1/5 (+) (1/7 (+) 1/9) = 1/5 (+) 1/4 = ?")
v1 = (1/7 + 1/9)/(1 + 1/63)
print(f"  1/7 (+) 1/9 = {v1:.10f} = 1/4 = {1/4:.10f}? {abs(v1-1/4)<1e-10}")
v2 = (1/5 + 1/4)/(1 + 1/20)
print(f"  1/5 (+) 1/4 = {v2:.10f}")
f = Fraction(1,5)
g = Fraction(1,4)
h = (f+g)/(1+f*g)
print(f"  = {h} = {float(h):.6f}")
print()
print("  1/5 (+) 1/4 = 9/19. Is 19 odd? YES. Is (19-1)/2 = 9 an integer? YES.")
print("  So 1/5 (+) 1/4 = 9/19 = 1/(19/9)... not of form 1/(2c+1).")
print("  BUT 1/4 is NOT of the form 1/(2n+1). It's 1/4 with 4 even.")
print()
print("  The rapidity composition IS associative (it's velocity addition).")
print("  But the SET {1/(2n+1)} is not closed under it.")
print("  The closed triples are scattered through the space.")
print()

# ============================================================
print("="*60)
print("GRAPH STRUCTURE OF RAPIDITY COMPOSITIONS")
print("="*60)
print()
print("  Draw a directed graph: edge a -> c if there exists b with")
print("  1/(2a+1) (+) 1/(2b+1) = 1/(2c+1) and a <= b.")
print("  (So a is always the smaller component.)")
print()

# Build adjacency
adj = defaultdict(set)
for a, b, c in solutions:
    adj[a].add(c)
    adj[b].add(c)

print("  Reachability from each node (up to a=15):")
for a in range(1, 16):
    if adj[a]:
        targets = sorted(adj[a])
        target_str = ", ".join(f"{c}(1/{2*c+1})" for c in targets[:8])
        if len(targets) > 8:
            target_str += f" ...+{len(targets)-8}"
        print(f"    a={a:2d} (1/{2*a+1:3d}) -> {target_str}")

print()

# ============================================================
print("="*60)
print("FORBIDDEN COMPOSITION: CAN WE REACH 7 AND 21?")
print("="*60)
print()
# Can we compose to get H=7? That means 1/(2c+1) = 1/7, c=3.
# Need a,b with ab/(a+b+1)=3, so ab=3(a+b+1)=3a+3b+3.
# ab-3a-3b=3, (a-3)(b-3)=12.
# Factor 12: 1*12, 2*6, 3*4.
# (a-3, b-3) in {(1,12),(2,6),(3,4)} -> (a,b) in {(4,15),(5,9),(6,7)}
print("  To compose to c=3 (i.e., 1/7 = the fourth):")
print("  Need (a-3)(b-3) = 12.")
print("  12 = 1*12 = 2*6 = 3*4")
print("  (a,b) = (4,15), (5,9), (6,7)")
for a, b in [(4,15), (5,9), (6,7)]:
    c = (a*b)//(a+b+1)
    print(f"    1/{2*a+1} (+) 1/{2*b+1} = 1/{2*c+1}: c = {c}")
print()

# c=10 means 1/21 (the forbidden H value scaled)
# (a-10)(b-10) = 10*11 = 110
# 110 = 1*110, 2*55, 5*22, 10*11
print("  To compose to c=10 (i.e., 1/21):")
print("  Need (a-10)(b-10) = 110.")
print("  110 = 1*110 = 2*55 = 5*22 = 10*11")
for d1, d2 in [(1,110),(2,55),(5,22),(10,11)]:
    a, b = 10+d1, 10+d2
    c = (a*b)//(a+b+1)
    print(f"    1/{2*a+1} (+) 1/{2*b+1} = 1/{2*c+1}: c = {c}")
print()

# The DEEP question: (a-c)(b-c) = c(c+1) = c^2+c
print("  GENERAL: c = ab/(a+b+1) <=> (a-c)(b-c) = c(c+1) = c^2+c")
print("  PROOF: ab = c(a+b+1). Then ab-ca-cb = c.")
print("  a(b-c) - c(b-c) = c + c^2 - c^2 = c ... hmm let me redo.")
print("  ab - c(a+b+1) = 0 => ab - ca - cb - c = 0")
print("  => a(b-c) = c(b+1) => (a-c)(b-c) = c(b+1) - c(b-c)")
print("  Actually: ab = ca + cb + c => ab - ca - cb + c^2 = c^2 + c")
print("  => (a-c)(b-c) = c^2 + c = c(c+1). QED!")
print()
print("  So the solutions are EXACTLY the factorizations of c(c+1):")
print("  c(c+1) = d1 * d2 with d1 <= d2, a = c+d1, b = c+d2.")
print()

# Display for small c
print("  c     c(c+1)    factorizations       (a,b) pairs")
print("  " + "-"*65)
for c in range(1, 12):
    cc = c*(c+1)
    facts = []
    for d1 in range(1, int(sqrt(cc))+1):
        if cc % d1 == 0:
            d2 = cc // d1
            a, b = c+d1, c+d2
            facts.append((d1, d2, a, b))
    fact_str = "; ".join(f"{d1}*{d2}->{a},{b}" for d1,d2,a,b in facts)
    print(f"  {c:3d}   {cc:5d}     {fact_str}")

print()
print("  THE NUMBER OF COMPOSITIONS TO c IS tau(c(c+1))/2")
print("  (half the number of divisors of c(c+1), since d1 <= d2).")
print()
print("  tau(c(c+1)) for c = 1,...,20:")
for c in range(1, 21):
    cc = c*(c+1)
    ndiv = sum(1 for d in range(1, cc+1) if cc % d == 0)
    ncomp = sum(1 for d in range(1, int(sqrt(cc))+1) if cc % d == 0)
    print(f"    c={c:2d}: c(c+1) = {cc:4d}, tau = {ndiv:3d}, compositions = {ncomp}")

print()
print("  Highly composite c(c+1) give many compositions.")
print("  c=6: 6*7=42=2*3*7, tau=8, 4 compositions.")
print("  c=12: 12*13=156=4*3*13, tau=12, 6 compositions.")
print()

# ============================================================
print("="*60)
print("RAPIDITY CHAINS: ITERATED COMPOSITIONS")
print("="*60)
print()
print("  Start with the octave 1/3 (c=1) and iterate:")
print("  Can we reach all 1/(2k+1) by composing with themselves?")
print()

# 1/3 (+) 1/3 = (2/3)/(1+1/9) = (2/3)/(10/9) = 6/10 = 3/5
# Q(3/5) = 4. rapidity doubles.
print("  1/3 (+) 1/3 = 3/5  (NOT of form 1/(2n+1), since 5 is odd but 3/5 != 1/k)")
print("  So the octave cannot even compose with itself to stay in the set!")
print()
print("  1/5 (+) 1/5 = (2/5)/(1+1/25) = (2/5)/(26/25) = 50/130 = 5/13")
print("  5/13: is this 1/(2c+1)? 13/5 = 2.6. Not integer.")
print()
print("  SELF-COMPOSITION 1/(2a+1) (+) 1/(2a+1):")
print("  = 2/(2a+1) / (1 + 1/(2a+1)^2)")
print("  = 2(2a+1) / ((2a+1)^2 + 1)")
print("  = (4a+2) / (4a^2+4a+2)")
print("  = (2a+1) / (2a^2+2a+1)")
print()
print("  For this to be 1/(2c+1), need (2a^2+2a+1)/(2a+1) = 2c+1")
print("  => 2a^2+2a+1 = (2a+1)(2c+1) = 4ac+2a+2c+1")
print("  => 2a^2 = 4ac+2c = 2c(2a+1)")
print("  => c = a^2/(2a+1)")
print("  Integer iff (2a+1) | a^2.")
print()
print("  gcd(a, 2a+1) = gcd(a, 1) = 1. So gcd(a^2, 2a+1) = 1.")
print("  Therefore (2a+1) | a^2 only if 2a+1 = 1, i.e., a = 0.")
print()
print("  THEOREM: No musical interval can compose with itself to")
print("  produce another musical interval. Self-composition ALWAYS")
print("  leaves the set {1/(2n+1)}.")
print()
print("  This is because gcd(a^2, 2a+1) = 1: consecutive-ish")
print("  numbers are always coprime. The PARITY obstruction")
print("  prevents self-doubling in rapidity space.")
print()

# ============================================================
print("="*60)
print("THE c(c+1) STRUCTURE AND PRONIC NUMBERS")
print("="*60)
print()
print("  The number of ways to compose to interval c is tau(c(c+1))/2")
print("  where c(c+1) is a PRONIC NUMBER (product of consecutive integers).")
print()
print("  Pronic numbers: 2, 6, 12, 20, 30, 42, 56, 72, 90, 110, 132, ...")
print("  = n(n+1) = 2*T_n where T_n is the n-th triangular number.")
print()
print("  c(c+1) is ALWAYS EVEN (one of c, c+1 is even).")
print("  c(c+1) is NEVER a perfect square (consecutive integers can't be).")
print()
print("  The rapidity algebra is governed by PRONIC FACTORIZATION.")
print("  Each factorization c(c+1) = d1*d2 gives a pair (c+d1, c+d2)")
print("  of musical interval indices whose relativistic composition")
print("  equals interval c.")
print()
print("  For the FORBIDDEN c=3 (interval 1/7):")
print("  c(c+1) = 12 = 1*12 = 2*6 = 3*4")
print("  Pairs: (4,15), (5,9), (6,7)")
print("  i.e., 1/9 (+) 1/31, 1/11 (+) 1/19, 1/13 (+) 1/15 all give 1/7.")
print("  THREE ways to reach the threshold.")
print()
print("  For the FORBIDDEN c=10 (interval 1/21):")
print("  c(c+1) = 110 = 1*110 = 2*55 = 5*22 = 10*11")
print("  Pairs: (11,120), (12,65), (15,32), (20,21)")
print("  FOUR ways to reach the forbidden number 21.")
print("  And 110 = 2*5*11 = the Paley primes {5, 11} times 2!")
print()
