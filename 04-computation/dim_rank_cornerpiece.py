#!/usr/bin/env python3
"""
dim_rank_cornerpiece.py — opus-2026-03-14-S79
=============================================

THE dim/rank RATIOS OF EXCEPTIONAL LIE GROUPS
AND THEIR CONNECTION TO CORNER PIECE NUMBERS 3ⁿ-2ⁿ

Discovery: dim(E₇)/rank(E₇) = 133/7 = 19 = 3³-2³

Is this part of a pattern? Let's check ALL exceptional Lie groups
and the classical series.

Parts:
1. dim/rank for all exceptionals
2. Connection to 3ⁿ-2ⁿ corner piece sequence
3. Classical Lie algebras dim/rank
4. The dim formula in terms of (2,3)
5. Grand pattern: dim = rank × (some recurrence value)
"""

from math import gcd
from fractions import Fraction

KEY1, KEY2 = 2, 3
f = lambda z: (z - KEY1) * (z - KEY2)

print("=" * 70)
print("PART 1: dim/rank FOR ALL EXCEPTIONAL LIE GROUPS")
print("=" * 70)
print()

exceptionals = [
    ("G₂", 2, 14, 6),
    ("F₄", 4, 52, 12),
    ("E₆", 6, 78, 12),
    ("E₇", 7, 133, 18),
    ("E₈", 8, 248, 30),
]

print(f"  {'Type':<5s} {'rank':>5s} {'dim':>5s} {'h':>5s} {'dim/rank':>10s} "
      f"{'as fraction':>12s} {'corner?':>10s}")
print(f"  {'-'*65}")

corner_pieces = {3**n - 2**n: n for n in range(10)}

for typ, r, d, h in exceptionals:
    ratio = Fraction(d, r)
    cp = corner_pieces.get(int(ratio), "—")
    print(f"  {typ:<5s} {r:>5d} {d:>5d} {h:>5d} {float(ratio):>10.4f} "
          f"{str(ratio):>12s} {'3^'+str(cp)+'-2^'+str(cp) if cp != '—' else '—':>10s}")

print()
print("  Corner piece sequence 3ⁿ-2ⁿ: ", end="")
for n in range(8):
    print(f"{3**n-2**n}", end=", ")
print("...")
print()

# Check which ratios match
print("  MATCHES:")
print(f"    dim(G₂)/rank(G₂) = 14/2 = 7 = 2²+3 = KEY₁²+KEY₂")
print(f"    dim(F₄)/rank(F₄) = 52/4 = 13 = KEY₁²+KEY₂²")
print(f"    dim(E₆)/rank(E₆) = 78/6 = 13 = KEY₁²+KEY₂²")
print(f"    dim(E₇)/rank(E₇) = 133/7 = 19 = 3³-2³ = KEY₂³-KEY₁³ ✓")
print(f"    dim(E₈)/rank(E₈) = 248/8 = 31 = 2⁵-1 = KEY₁^5-1")
print()

# Now 7, 13, 13, 19, 31 — are these ALL in the (2,3) world?
print("  ALL dim/rank ratios in terms of (2,3):")
print(f"     7 = 2²+3¹   = KEY₁² + KEY₂")
print(f"    13 = 2²+3²   = KEY₁² + KEY₂² (appears TWICE)")
print(f"    19 = 3³-2³   = KEY₂³ - KEY₁³ (corner piece n=3)")
print(f"    31 = 2⁵-1    = KEY₁^(KEY₁+KEY₂) - 1 (Mersenne)")
print()

# Alternative expressions
print("  DEEPER PATTERNS:")
print(f"     7 = 2³-1    (Mersenne prime, M₃)")
print(f"    13 = 2⁴-3    (power of 2 minus KEY₂)")
print(f"    19 = 2⁵-13   ... no. 19 = 3³-2³ better")
print(f"    31 = 2⁵-1    (Mersenne prime, M₅)")
print()

# As sums/differences of powers of 2 and 3
print("  As a^i ± b^j where {a,b} ⊆ {2,3}:")
for d_over_r in [7, 13, 19, 31]:
    reps = []
    for a in [2, 3]:
        for b in [2, 3]:
            for i in range(1, 8):
                for j in range(0, 8):
                    if a**i + b**j == d_over_r:
                        reps.append(f"{a}^{i}+{b}^{j}")
                    if a**i - b**j == d_over_r:
                        reps.append(f"{a}^{i}-{b}^{j}")
    print(f"    {d_over_r:>3d} = {', '.join(reps[:5])}")

print()
print("=" * 70)
print("PART 2: THE (2,3) RECURRENCE AT dim/rank")
print("=" * 70)
print()

# General solution: a(n) = A·2ⁿ + B·3ⁿ
# What (A,B) gives 7, 13, 19, 31?

print("  Can we find A,B such that A·2ⁿ + B·3ⁿ = dim/rank for n=1,...?")
print()

# If dim/rank follows A·2ⁿ + B·3ⁿ for some sequence of n-values...
# But the groups are discrete, not a continuous sequence

# Instead: which sequences A·2ⁿ + B·3ⁿ contain all four values?
targets = [7, 13, 19, 31]

# For each pair of targets, solve for (A,B)
print("  For each pair of consecutive values, solve A·2ⁿ + B·3ⁿ:")
for n1, n2 in [(1,2), (2,3), (3,4), (1,3), (2,4), (3,5)]:
    for i, v1 in enumerate(targets):
        for j, v2 in enumerate(targets):
            if i >= j:
                continue
            # v1 = A·2^n1 + B·3^n1
            # v2 = A·2^n2 + B·3^n2
            det = 2**n1 * 3**n2 - 2**n2 * 3**n1
            if det != 0:
                A = Fraction(v1 * 3**n2 - v2 * 3**n1, det)
                B = Fraction(v2 * 2**n1 - v1 * 2**n2, det)
                # Check if A, B are nice
                if A.denominator <= 6 and B.denominator <= 6:
                    seq = [int(A * 2**n + B * 3**n) for n in range(1, 7)]
                    hits = sum(1 for t in targets if t in seq)
                    if hits >= 2:
                        print(f"    n=({n1},{n2}) for ({v1},{v2}): "
                              f"A={A}, B={B}, seq={seq}, hits={hits}")

print()
# Direct check: are 7,13,19,31 part of 3ⁿ-2ⁿ?
print("  Direct check against 3ⁿ-2ⁿ:")
for n in range(1, 10):
    val = 3**n - 2**n
    marker = " ← dim/rank!" if val in targets else ""
    print(f"    n={n}: 3^{n}-2^{n} = {val}{marker}")

print()
print(f"  Only 19 = 3³-2³ matches.")
print()

# Check 2ⁿ-1 (Mersenne)
print("  Direct check against 2ⁿ-1 (Mersenne):")
for n in range(1, 10):
    val = 2**n - 1
    marker = " ← dim/rank!" if val in targets else ""
    print(f"    n={n}: 2^{n}-1 = {val}{marker}")

print()
print(f"  7 = 2³-1 (M₃) and 31 = 2⁵-1 (M₅) match.")
print(f"  Both are Mersenne PRIMES!")
print()

# Check 2ⁿ+3ⁿ (sum of key powers)
print("  Direct check against KEY₁ⁿ+KEY₂ⁿ = 2ⁿ+3ⁿ:")
for n in range(0, 8):
    val = 2**n + 3**n
    marker = " ← dim/rank!" if val in targets else ""
    print(f"    n={n}: 2^{n}+3^{n} = {val}{marker}")

print()
print(f"  13 = 2²+3² matches! (for both F₄ and E₆)")

print()
print("=" * 70)
print("PART 3: UNIFIED FORMULA — dim/rank AND THE KEYS")
print("=" * 70)
print()

# Summary:
# G₂: dim/rank = 7 = 2³-1
# F₄: dim/rank = 13 = 2²+3²
# E₆: dim/rank = 13 = 2²+3²
# E₇: dim/rank = 19 = 3³-2³
# E₈: dim/rank = 31 = 2⁵-1

print("  FORMULA CANDIDATES:")
print()

# Pattern 1: dim/rank = 2^a + (-1)^b · 3^c
print("  dim/rank = 2^a ± 3^c form:")
print(f"    G₂:  7 = 2³ - 3⁰ = 8-1")
print(f"    F₄: 13 = 2² + 3² = 4+9")
print(f"    E₆: 13 = 2² + 3² = 4+9")
print(f"    E₇: 19 = -2³ + 3³ = -8+27")
print(f"    E₈: 31 = 2⁵ - 3⁰ = 32-1")
print()

# The exponents (a,c) are:
print("  Exponents in 2^a ± 3^c:")
print(f"    G₂: (3,0) with sign -")
print(f"    F₄: (2,2) with sign +")
print(f"    E₆: (2,2) with sign +")
print(f"    E₇: (3,3) with sign +  (reading as 3³-2³)")
print(f"    E₈: (5,0) with sign -")
print()
print(f"  The exponents: (3,0), (2,2), (2,2), (3,3), (5,0)")
print(f"  For G₂ and E₈: second exponent is 0 (pure Mersenne)")
print(f"  For F₄ and E₆: equal exponents (sum of equal powers)")
print(f"  For E₇: equal exponents 3 (difference of equal powers)")
print()
print(f"  The first exponents: 3, 2, 2, 3, 5")
print(f"  = KEY₂, KEY₁, KEY₁, KEY₂, KEY₁+KEY₂ !")
print()

# Pattern 2: look at dim - rank (= number of roots)
print("  dim - rank = number of positive + negative roots:")
for typ, r, d, h in exceptionals:
    roots = d - r
    print(f"    {typ}: dim-rank = {d}-{r} = {roots} = "
          f"2·{roots//2} positive roots")

print()
print("  Number of POSITIVE roots:")
for typ, r, d, h in exceptionals:
    pos = (d - r) // 2
    print(f"    {typ}: {pos} positive roots = "
          f"rank·(dim/rank - 1)/2 = {r}·{(d//r - 1)//2 if d % r == 0 and (d//r-1) % 2 == 0 else Fraction(d-r, 2*r)}")

print()

# The h · rank connection
print("  Checking: dim = rank · h + something?")
for typ, r, d, h in exceptionals:
    rh = r * h
    diff = d - rh
    print(f"    {typ}: rank·h = {r}·{h} = {rh}, dim = {d}, diff = {diff}")

print()
print(f"  None of them has dim = rank·h exactly.")
print(f"  But: dim = rank·(h-1) + rank + something?")
for typ, r, d, h in exceptionals:
    print(f"    {typ}: rank·(h-1) = {r}·{h-1} = {r*(h-1)}, "
          f"dim-rank·(h-1) = {d - r*(h-1)}")

print()
# dim = rank + 2·(#positive roots) = rank + (number of roots)
# For a simple Lie algebra, dim = rank + |Φ| where Φ is the root system
# |Φ| = 2·(#positive roots)
# For E₈: |Φ| = 240, rank = 8, dim = 248 ✓

# The number of positive roots:
print("  Number of positive roots = dim/rank formula:")
print(f"  For type X_r: #Φ⁺ = (dim-rank)/2")
for typ, r, d, h in exceptionals:
    phi_plus = (d - r) // 2
    # This should equal r·h/2 (well-known formula: |Φ⁺| = rh/2)
    rh2 = r * h // 2
    check = "✓" if phi_plus == rh2 else "✗"
    print(f"    {typ}: #Φ⁺ = {phi_plus}, r·h/2 = {rh2} {check}")

print()
print(f"  YES! #Φ⁺ = r·h/2 always.")
print(f"  Therefore: dim = rank + 2·(rank·h/2) = rank·(h+1)")
print(f"  Wait: dim = rank + rank·h = rank·(1+h)?")
for typ, r, d, h in exceptionals:
    check = "✓" if d == r * (h + 1) else "✗"
    print(f"    {typ}: rank·(h+1) = {r}·{h+1} = {r*(h+1)}, dim = {d} {check}")

print()
print(f"  WRONG. dim ≠ rank·(h+1) in general.")
print(f"  The correct formula: dim = rank + r·h = rank·(1+h)... no wait")
print(f"  dim = rank + |Φ| = rank + r·h")
print(f"  But |Φ| = 2·|Φ⁺| = 2·(r·h/2) = r·h")
for typ, r, d, h in exceptionals:
    check = "✓" if d == r + r * h else "✗"
    rh = r * h
    print(f"    {typ}: rank + rank·h = {r} + {rh} = {r + rh}, dim = {d} {check}")

print()
print(f"  So dim = rank·(1 + h) FOR ALL SIMPLE LIE ALGEBRAS!")
print(f"  This means: dim/rank = h + 1")
print()

for typ, r, d, h in exceptionals:
    print(f"    {typ}: dim/rank = {d//r} = h+1 = {h}+1 = {h+1} "
          f"{'✓' if d//r == h+1 else '✗'}")

print()
print(f"  WAIT: dim/rank should be h+1?")
print(f"    G₂: 14/2 = 7, h+1 = 7 ✓")
print(f"    F₄: 52/4 = 13, h+1 = 13 ✓")
print(f"    E₆: 78/6 = 13, h+1 = 13 ✓")
print(f"    E₇: 133/7 = 19, h+1 = 19 ✓")
print(f"    E₈: 248/8 = 31, h+1 = 31 ✓")
print()
print(f"  YES!!! dim/rank = h+1 FOR ALL EXCEPTIONAL LIE GROUPS!")
print(f"  This is because dim = rank + |Φ| = rank + rank·h = rank·(h+1)")
print()

# So the "corner piece" connection is really about h+1!
print("=" * 70)
print("PART 4: THE h+1 VALUES AND (2,3) EXPRESSIONS")
print("=" * 70)
print()

print("  The Coxeter numbers plus 1:")
for typ, r, d, h in exceptionals:
    print(f"    {typ}: h+1 = {h+1}", end="")
    val = h + 1
    reps = []
    for a in [2, 3]:
        for i in range(0, 8):
            for b in [2, 3]:
                for j in range(0, 8):
                    if a != b or i != j:
                        if a**i + b**j == val:
                            reps.append(f"{a}^{i}+{b}^{j}")
                        if a**i - b**j == val and a**i > b**j:
                            reps.append(f"{a}^{i}-{b}^{j}")
    print(f" = {', '.join(reps[:3])}")

print()
print(f"  So the real question is: what are the h-values?")
print(f"  h(G₂) = 6 = KEY₁·KEY₂")
print(f"  h(F₄) = 12 = KEY₁²·KEY₂")
print(f"  h(E₆) = 12 = KEY₁²·KEY₂")
print(f"  h(E₇) = 18 = KEY₁·KEY₂²")
print(f"  h(E₈) = 30 = KEY₁·KEY₂·(KEY₁+KEY₂)")
print()

# The h-values factored by 2 and 3
print("  h-values in KEY₁^a · KEY₂^b · (KEY₁+KEY₂)^c form:")
print(f"    h(G₂) = 6  = 2¹·3¹ → (a,b,c) = (1,1,0)")
print(f"    h(F₄) = 12 = 2²·3¹ → (a,b,c) = (2,1,0)")
print(f"    h(E₆) = 12 = 2²·3¹ → (a,b,c) = (2,1,0)")
print(f"    h(E₇) = 18 = 2¹·3² → (a,b,c) = (1,2,0)")
print(f"    h(E₈) = 30 = 2¹·3¹·5¹ → (a,b,c) = (1,1,1)")
print()
print(f"  BEAUTIFUL: E₈ is the ONLY one that needs KEY₁+KEY₂ = 5")
print(f"  G₂, F₄, E₆ use only KEY₁ and KEY₂")
print(f"  E₇ uses KEY₁¹·KEY₂² (more KEY₂ than KEY₁)")
print(f"  E₈ uses all three: KEY₁·KEY₂·(KEY₁+KEY₂)")

print()
print("  dim = rank · (h+1):")
for typ, r, d, h in exceptionals:
    h1 = h + 1
    print(f"    dim({typ}) = {r}·{h1} = {r*h1}")
    # Factor r and h+1
    print(f"      rank = {r}, h+1 = {h1}", end="")
    if h1 in corner_pieces:
        print(f" = 3^{corner_pieces[h1]}-2^{corner_pieces[h1]}", end="")
    if (2**r - 1) == h1:
        print(f" = 2^rank-1 (Mersenne!)", end="")
    print()

print()
print(f"  FOR G₂: h+1 = 7 = 2³-1 = 2^rank+1 - 1 (rank=2, so 2^3-1) ✓")
print(f"  FOR E₈: h+1 = 31 = 2⁵-1 = 2^rank-3 - 1 (rank=8, so 2^5-1) ✓")
print(f"  Both G₂ and E₈ have h+1 = Mersenne prime!")
print()

# The connection to dim(E₇)/rank(E₇) = 19 = 3³-2³
print("  THE FINAL CONNECTION:")
print(f"  dim(E₇)/rank(E₇) = 19 = 3³-2³ = h(E₇)+1")
print(f"  This is TRIVIALLY TRUE because dim/rank = h+1 always!")
print(f"  But the INTERESTING fact is that h(E₇)+1 = 19 = 3³-2³")
print(f"  i.e., h(E₇) = 18 = 3³-2³-1 = 27-8-1")
print(f"  Or equivalently: h(E₇) = KEY₂³ - KEY₁³ - 1")
print()

# Check all h values against corner pieces
print("  Which h values are related to 3ⁿ-2ⁿ?")
for typ, r, d, h in exceptionals:
    for n in range(1, 8):
        cp = 3**n - 2**n
        if h == cp:
            print(f"    h({typ}) = {h} = 3^{n}-2^{n} ✓")
        elif h + 1 == cp:
            print(f"    h({typ})+1 = {h+1} = 3^{n}-2^{n}")
        elif h - 1 == cp:
            print(f"    h({typ})-1 = {h-1} = 3^{n}-2^{n}")

print()
print("=" * 70)
print("PART 5: THE CLASSICAL LIE ALGEBRAS — SAME FORMULA")
print("=" * 70)
print()

# Classical: A_n (sl(n+1)), B_n (so(2n+1)), C_n (sp(2n)), D_n (so(2n))
classicals = [
    ("A₁=B₁=C₁", 1, 3, 2),
    ("A₂", 2, 8, 3),
    ("B₂=C₂", 2, 10, 4),
    ("A₃=D₃", 3, 15, 4),
    ("B₃", 3, 21, 6),
    ("C₃", 3, 21, 6),
    ("A₄", 4, 24, 5),
    ("B₄", 4, 36, 8),
    ("C₄", 4, 36, 8),
    ("D₄", 4, 28, 6),
    ("A₅", 5, 35, 6),
    ("D₅", 5, 45, 8),
]

print(f"  {'Type':<12s} {'rank':>4s} {'dim':>4s} {'h':>3s} {'dim/rank':>8s} {'h+1':>4s} {'match':>5s}")
for typ, r, d, h in classicals:
    ratio = Fraction(d, r)
    h1 = h + 1
    match = "✓" if int(ratio) == h1 else "✗"
    print(f"  {typ:<12s} {r:>4d} {d:>4d} {h:>3d} {str(ratio):>8s} {h1:>4d} {match:>5s}")

print()
print(f"  dim/rank = h+1 holds for ALL simple Lie algebras!")
print(f"  This is the well-known formula: dim(g) = rank(g)·(h(g)+1)")
print(f"  derived from dim(g) = rank(g) + |Φ| = rank(g) + rank(g)·h(g)")
print()
print(f"  Where |Φ| = rank·h comes from the fact that the")
print(f"  average squared root length times rank gives |Φ|·h... ")
print(f"  Actually: |Φ| = rank·h is well-known for simply-laced types")
print(f"  and extends to all types with appropriate normalization.")

print()
print("=" * 70)
print("GRAND SYNTHESIS")
print("=" * 70)
print()

print("The dim/rank = h+1 identity means the 'corner piece' connection")
print("is really about Coxeter numbers:")
print()
print(f"  h(G₂) + 1 = 7 = 2³-1 (Mersenne prime M₃)")
print(f"  h(F₄) + 1 = 13 = 2²+3² (sum of key squares)")
print(f"  h(E₆) + 1 = 13 = 2²+3² (same!)")
print(f"  h(E₇) + 1 = 19 = 3³-2³ (corner piece at n=3)")
print(f"  h(E₈) + 1 = 31 = 2⁵-1 (Mersenne prime M₅)")
print()
print("The Coxeter numbers themselves:")
print(f"  h(G₂) = 6 = 2·3 = KEY₁·KEY₂")
print(f"  h(F₄) = 12 = 4·3 = KEY₁²·KEY₂")
print(f"  h(E₆) = 12 = 4·3 = KEY₁²·KEY₂")
print(f"  h(E₇) = 18 = 2·9 = KEY₁·KEY₂²")
print(f"  h(E₈) = 30 = 2·3·5 = KEY₁·KEY₂·(KEY₁+KEY₂)")
print()
print("The 2^a·3^b decomposition of h:")
print(f"  G₂: 2¹3¹ → exponent vector (1,1)")
print(f"  F₄: 2²3¹ → exponent vector (2,1)")
print(f"  E₆: 2²3¹ → exponent vector (2,1)")
print(f"  E₇: 2¹3² → exponent vector (1,2)")
print(f"  E₈: 2¹3¹5¹ → needs the THIRD key")
print()
print("Pattern: exponent vectors (a,b) where h = 2^a · 3^b:")
print(f"  G₂: (1,1) — balanced")
print(f"  F₄ = E₆: (2,1) — KEY₁-heavy")
print(f"  E₇: (1,2) — KEY₂-heavy")
print(f"  E₈: needs 5 = KEY₁+KEY₂ (transcends the two-key system)")
print()
print("  G₂ and E₈ are the BOOKENDS:")
print(f"    G₂ starts with h = KEY₁·KEY₂ (the product)")
print(f"    E₈ ends with h = KEY₁·KEY₂·(KEY₁+KEY₂) (the full triple)")
print(f"    Ratio: h(E₈)/h(G₂) = 30/6 = 5 = KEY₁+KEY₂")
print()
print("  F₄/E₆ and E₇ are DUAL:")
print(f"    h(F₄) = h(E₆) = 12 = KEY₁²·KEY₂")
print(f"    h(E₇) = 18 = KEY₁·KEY₂²")
print(f"    Ratio: h(E₇)/h(E₆) = 18/12 = 3/2 = KEY₂/KEY₁")
print(f"    This is the GOLDEN RATIO of the keys!")
