#!/usr/bin/env python3
"""
THM-J SPECTRAL CONNECTION
===========================
We proved: v_2(M_amp(n,d,s=0)) = -d - s_2(n-2-d)

At d=1: v_2(M_amp(n,1,0)) = -1 - s_2(n-3)

THM-J says: S(T) mod 2^{n-1} is universal iff s_2(n-3) <= 1.

So the universality criterion = "degree-1 M amplitude has v_2 >= -2".

Is this a coincidence or does the Walsh spectrum CONTROL the 2-adic structure
of S(T) in a precise way?

Also explore: what patterns exist across ALL degrees?

opus-2026-03-16-S73
"""
from math import factorial, log2
from fractions import Fraction

def s2(n):
    """Binary digit sum of n."""
    return bin(n).count('1')

def v2(n):
    """2-adic valuation of n."""
    if n == 0: return float('inf')
    n = abs(n)
    v = 0
    while n % 2 == 0:
        v += 1
        n //= 2
    return v

print("=" * 70)
print("THM-J SPECTRAL CONNECTION")
print("=" * 70)

print("\n--- v_2(M_amp(n,d,s=0)) = -d - s_2(n-2-d) ---")
print(f"\n{'n':>4} {'d':>4} {'f=n-2-d':>8} {'s_2(f)':>7} {'v_2(amp)':>9} {'amp=f!/2^(n-2)':>16}")

for n in range(3, 22, 2):  # odd n
    for d in range(1, n-1, 2):  # odd d (M degrees for odd n)
        f = n - 2 - d
        if f < 0: continue
        s2f = s2(f)
        v2_amp = -d - s2f
        amp = Fraction(factorial(f), 2**(n-2))
        print(f"{n:>4} {d:>4} {f:>8} {s2f:>7} {v2_amp:>9} {str(amp):>16}")
    print()

print("=" * 70)
print("KEY: v_2 at degree 1 vs THM-J criterion")
print("=" * 70)

print(f"\n{'n':>4} {'n-3':>5} {'s_2(n-3)':>9} {'v_2(d=1)':>9} {'Universal?':>11} {'Interpretation':>30}")
for n in range(3, 40, 2):
    s2_nm3 = s2(n - 3)
    v2_d1 = -1 - s2_nm3
    univ = "YES" if s2_nm3 <= 1 else "NO"
    # Interpretation
    if s2_nm3 == 0:
        interp = f"n-3=0: trivial (n=3)"
    elif s2_nm3 == 1:
        interp = f"n-3=2^{int(log2(n-3))}: power of 2"
    else:
        interp = f"s_2={s2_nm3}: complex binary"
    print(f"{n:>4} {n-3:>5} {s2_nm3:>9} {v2_d1:>9} {univ:>11} {interp:>30}")

print("\n" + "=" * 70)
print("THE PATTERN: v_2 at each degree d across n")
print("=" * 70)

print("\n    ", end="")
for d in range(1, 18, 2):
    print(f"{'d='+str(d):>6}", end="")
print()
for n in range(3, 24, 2):
    print(f"n={n:>2}:", end="")
    for d in range(1, min(n-1, 18), 2):
        f = n - 2 - d
        if f >= 0:
            val = -d - s2(f)
            print(f"{val:>6}", end="")
        else:
            print(f"{'—':>6}", end="")
    print()

print("\n" + "=" * 70)
print("DISCOVERY: The v_2 staircase")
print("=" * 70)

print("""
At fixed d, v_2(M_amp(n,d,0)) = -d - s_2(n-2-d) depends on n through s_2.

s_2(m) has a fractal structure:
  m = 0: s_2 = 0
  m = 1: s_2 = 1
  m = 2: s_2 = 1
  m = 3: s_2 = 2
  m = 4: s_2 = 1
  m = 5: s_2 = 2
  m = 6: s_2 = 2
  m = 7: s_2 = 3
  m = 8: s_2 = 1

At d=1: v_2 = -1 - s_2(n-3). The MINIMUM v_2 (most 2-adic precision)
occurs when s_2(n-3) is maximized, i.e., n-3 = 2^k - 1 (all 1s in binary).
  n = 4: v_2 = -2 (n-3=1, s_2=1)
  n = 6: v_2 = -3 (n-3=3, s_2=2)  <-- n=6 AGAIN!
  n = 10: v_2 = -4 (n-3=7, s_2=3)
  n = 18: v_2 = -5 (n-3=15, s_2=4)

The n where v_2(d=1) is MOST negative is n = 2^k + 3:
  n = 4, 6, 10, 18, 34, 66, ...

These are the LEAST universal n for the signed HP permanent!
""")

# Verify the "least universal" claim
print("Least-universal n (s_2(n-3) maximized for each bit-length):")
for k in range(1, 10):
    n = (1 << k) + 2  # n-3 = 2^k - 1
    s2_val = s2(n - 3)
    print(f"  n={n}: n-3={n-3}={bin(n-3)}, s_2={s2_val}, v_2(d=1)={-1-s2_val}")

print("\n" + "=" * 70)
print("PRODUCT STRUCTURE: Total 2-adic weight of Walsh spectrum")
print("=" * 70)

for n in [5, 7, 9, 11]:
    print(f"\n  n={n}:")
    total_v2 = 0
    count = 0
    for d in range(1, n-1, 2):
        f = n - 2 - d
        if f < 0: continue
        v2_base = -d - s2(f)
        # Count valid monomials at this degree (approximation: just s=0)
        print(f"    d={d}: v_2={v2_base}, amp={Fraction(factorial(f), 2**(n-2))}")
        total_v2 += v2_base
        count += 1
    print(f"    Sum of v_2 values = {total_v2}, count = {count}, mean = {total_v2/count:.3f}")

print("\n" + "=" * 70)
print("EULER CHARACTERISTIC CONNECTION")
print("=" * 70)
print("""
The Euler characteristic of the chain complex Omega at dimension d
involves alternating sums of dim(Omega_p). The v_2 staircase at each
degree encodes how "binary-complex" the free vertex count is.

Key observation: the v_2 values at d=max (the top degree) are always
-(n-2) for odd n. This is because f=0, s_2(0)=0, so v_2 = -(n-2).

The v_2 values at d=1 (bottom degree) are -1 - s_2(n-3).

The RANGE of v_2 across degrees is:
  max(v_2) - min(v_2) = (-1 - s_2(n-3)) - (-(n-2)) = n - 3 - s_2(n-3)
  = v_2((n-3)!)  [by Legendre's formula!]
""")

for n in range(3, 22, 2):
    d_min, d_max = 1, n - 2
    v2_top = -1 - s2(n-3)  # at d=1
    v2_bot = -(n-2)          # at d=n-2 (wait, d must be odd...)
    # Actually at d=n-2: f=0, v_2=-d-0=-(n-2)
    v2_range = v2_top - v2_bot
    leg = (n-3) - s2(n-3)  # v_2((n-3)!)
    print(f"  n={n}: v_2(d=1)={v2_top}, v_2(d={n-2})={v2_bot}, range={v2_range}, v_2((n-3)!)={leg}, match={v2_range==leg}")

print("\n" + "=" * 70)
print("THE SPECTRAL LEGENDRE IDENTITY")
print("=" * 70)
print("""
v_2(M_amp(n,1,0)) - v_2(M_amp(n,n-2,0)) = v_2((n-3)!)

This says: the 2-adic "spread" of the M Walsh spectrum equals
the 2-adic valuation of (n-3)! — exactly Legendre's formula applied
to the free vertex count at the bottom of the spectrum.

More precisely: the v_2 range across the M spectrum ENCODES the
factorial structure of the free vertices. The binary digit sum s_2
acts as a "complexity penalty" that reduces the effective spread.
""")

print("\nDone.")
