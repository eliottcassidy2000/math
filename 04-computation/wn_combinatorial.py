#!/usr/bin/env python3
"""
Investigate combinatorial interpretations of W(n).

W(n) = Σ_{σ ∈ NUD(n)} 2^{adj1(σ)}

Since 2^{adj1(σ)} = Σ_{S ⊆ asc_adj(σ)} 1, we can interpret W(n) as:
W(n) = #{(σ, S) : σ ∈ NUD(n), S ⊆ asc_positions(σ)}

where asc_positions(σ) = {i : σ(i+1) = σ(i) + 1}.

Equivalently: count pairs (σ, S) where σ has no unit descents and S marks
a subset of unit ascent positions.

This can also be interpreted as: count permutations σ with each unit ascent
position colored either red or blue, and no unit descents.

Another interpretation: W(n) = #{words in a modified permutation structure}

Let me also check: does W(n) count anything simpler?

W(1)=1, W(2)=2, W(3)=8, W(4)=32, W(5)=158, W(6)=928

Is W(n)/2^n integer?
W(3)/8=1, W(4)/16=2, W(5)/32=4.9375. No.

Is there a connection to derangements? D(n) = n!(1 - 1 + 1/2! - ... + (-1)^n/n!)
D(3)=2, D(4)=9, D(5)=44, D(6)=265. Not matching.

Let me check OEIS for W(n) sequence: 1, 2, 8, 32, 158, 928, 6350, 49752
"""
from fractions import Fraction as F
import math

wn = [0, 1, 2, 8, 32, 158, 928, 6350, 49752, 439670, 4327904]

# Check various normalizations
print("=== W(n) / various normalizations ===")
for n in range(1, 11):
    w = wn[n]
    nfact = math.factorial(n)
    print(f"n={n}: W={w}, W/n!={F(w,nfact)}, W/(n+1)!={F(w,math.factorial(n+1))}")

# Check: W(n) = number of (something) on [n+1] or [n-1]?
# W(3)=8 = 2^3. W(4)=32 = 2^5. W(5) = 158. Hmm.

# Let me check if W(n) counts labeled structures on n elements.
# W(n)/n!: 1, 1, 4/3, 4/3, 79/60, ...
print("\n=== W(n)/n! exact ===")
for n in range(1, 11):
    print(f"  W({n})/n! = {F(wn[n], math.factorial(n))}")

# W(3)/3! = 8/6 = 4/3
# W(4)/4! = 32/24 = 4/3
# W(5)/5! = 158/120 = 79/60
# Not obviously clean.

# Let me check W(n) - n! (the "excess"):
print("\n=== W(n) - n! ===")
for n in range(1, 11):
    excess = wn[n] - math.factorial(n)
    print(f"  n={n}: {excess}")
# W(n) - n!: 0, 0, 2, 8, 38, 208, 1330, 9432, 76910, ...
# This sequence: let me check what it matches
excess = [wn[n] - math.factorial(n) for n in range(3, 11)]
print(f"\nExcess sequence (n>=3): {excess}")

# Check: is this 2*(n-1)! + something?
print("\n=== (W(n) - n!) / (n-1)! ===")
for n in range(3, 11):
    e = wn[n] - math.factorial(n)
    print(f"  n={n}: {F(e, math.factorial(n-1))}")

# Check if excess satisfies subfactorial-like recurrence
print("\n=== Excess recurrence check ===")
exc = {n: wn[n] - math.factorial(n) for n in range(1, 11)}
for n in range(3, 11):
    # (n-1) * exc(n-1) + (n-2) * exc(n-2)?
    pred = (n-1)*exc[n-1] + (n-2)*exc[n-2]
    actual = exc[n]
    print(f"  n={n}: (n-1)*exc(n-1)+(n-2)*exc(n-2) = {pred}, actual = {actual}, diff = {actual - pred}")

# The diff is the correction W(n) - (n-1)W(n-1) - (n-2)W(n-2) shifted by factorials
# Since W satisfies NUD recurrence PLUS correction, and n! satisfies n*n!/(n+1)... hmm

# Let me try: does W(n) = NUD(n) * something nice?
nud = {1: 1, 2: 1, 3: 3, 4: 11, 5: 53, 6: 309, 7: 2119, 8: 16687, 9: 148329, 10: 1468457}
print("\n=== W(n)/NUD(n) ===")
for n in range(1, 11):
    print(f"  n={n}: W/NUD = {F(wn[n], nud[n])} = {float(F(wn[n], nud[n])):.10f}")

# W/NUD: 1, 2, 8/3, 32/11, 158/53, 928/309, 6350/2119, ...
# → e ≈ 2.718...

# Check: W(n) = NUD(n) * (1 + 1 + 1/2! + 1/3! + ... + 1/(n-1)!) exactly?
# No, because W(n)/NUD(n) approaches e from above for small n
# Actually from below: 1, 2, 2.667, 2.909, 2.981, 3.003(? no, 928/309=3.003?)
print(f"\n  928/309 = {928/309:.6f}")  # Should be ~3.003 — but e=2.718. Let me recheck.
# Wait, 928/309 = 3.003, but e = 2.718. So W/NUD does NOT approach e from above.
# Let me check the data table from THM-219.
# THM-219 says W(5)/NUD(5) = 158/53 ≈ 2.981, W(10)/NUD(10) ≈ 2.947
# So it's decreasing toward e. At n=5: 2.981 > e. OK.

# But W(6)/NUD(6) = 928/309 ≈ 3.003?? That's HIGHER than n=5.
print(f"  158/53 = {158/53:.6f}")
print(f"  928/309 = {928/309:.6f}")
# Hmm, 158/53 = 2.981, 928/309 = 3.003. That's going UP. But THM-219 data shows decrease.
# Let me check NUD values carefully.
print(f"\n  NUD(6) = 309? Recurrence: (5)*NUD(5) + (4)*NUD(4) = 5*53 + 4*11 = 265+44 = 309. Yes.")
print(f"  W(6) = 928, W/NUD = {928/309:.6f}")
# 3.003. But the THM-219 table shows W(5)/NUD(5)=2.981 and W(10)/NUD(10)=2.947
# So the ratio goes 2.981, 3.003, then decreases? Non-monotone convergence to e.
print("\n=== W(n)/NUD(n) detailed ===")
for n in range(1, 11):
    ratio = F(wn[n], nud[n])
    print(f"  n={n}: {float(ratio):.10f}")
