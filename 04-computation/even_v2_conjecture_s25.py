#!/usr/bin/env python3
"""
opus-2026-04-05-S25: Even-n v_2 conjecture.

DISCOVERY: For even n, v_2(T(n)) = n/2 + v_2(c(n))
where c(n) = #{odd k : 1 <= k < n/2, gcd(k,n) = 1}
         = number of coprime 2-part all-odd partitions of n.

Key insight: ALL individual odd coefficients at minimum v_2 level are ODD
(since n_odd/(ab) with a,b odd coprime is always odd).
So the sum of c(n) odd numbers has v_2 = v_2(c(n)).

Formula for c(n): with n = 2^a * m (m odd):
  a >= 2: c = 2^{a-2} * phi(m)
  a = 1:  c = phi(m)/2
"""

from math import factorial, gcd
from collections import Counter

def v2(n):
    if n == 0: return float('inf')
    v = 0
    while n % 2 == 0: v += 1; n //= 2
    return v

def euler_phi(n):
    result = n
    p = 2
    temp = n
    while p * p <= temp:
        if temp % p == 0:
            while temp % p == 0:
                temp //= p
            result -= result // p
        p += 1
    if temp > 1:
        result -= result // temp
    return result

T_val = {
    1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880,
    9:191536, 10:9733056, 11:903753248, 12:154108311168,
    13:48542114686912, 14:28401423719122304,
    15:31021002160355166848,
    16:63530415842308265100288,
    17:244912778438520759443245824,
    18:1783398846284777975419600287232,
    19:24605641171260376770598003978281472,
    20:645022068557873570931850526424042500096,
}

def c_n(n):
    """Number of coprime 2-part all-odd partitions of even n."""
    return sum(1 for k in range(1, n//2) if k % 2 == 1 and gcd(k, n) == 1)

def c_n_formula(n):
    """Closed form: n = 2^a * m, c = 2^{a-2}*phi(m) if a>=2, phi(m)/2 if a=1."""
    a = v2(n)
    m = n >> a
    if m == 1:
        return 2**(a-2) if a >= 2 else 0
    if a >= 2:
        return 2**(a-2) * euler_phi(m)
    else:  # a = 1
        return euler_phi(m) // 2

print("=" * 75)
print("CONJECTURE: v_2(T(n)) = n/2 + v_2(c(n)) for even n >= 4")
print("  c(n) = #{odd k : 1 <= k < n/2, gcd(k,n) = 1}")
print("=" * 75)
print()

print(f"{'n':>4} {'v_2(T)':>7} {'n/2':>4} {'c(n)':>5} {'v_2(c)':>6} "
      f"{'predicted':>9} {'match':>6} {'c formula':>10}")
print("-" * 60)

all_match = True
for n in range(4, 21, 2):
    actual = v2(T_val[n])
    c = c_n(n)
    c_form = c_n_formula(n)
    v2c = v2(c) if c > 0 else 0
    predicted = n // 2 + v2c
    match = "✓" if predicted == actual else "✗"
    c_match = "✓" if c == c_form else "✗"
    if predicted != actual: all_match = False
    print(f"{n:>4} {actual:>7} {n//2:>4} {c:>5} {v2c:>6} {predicted:>9} {match:>6} "
          f"{c_form:>8} {c_match:>2}")

print()
if all_match:
    print("✓ CONJECTURE VERIFIED for all even n from 4 to 20.")
else:
    print("✗ CONJECTURE FAILS somewhere.")

print()
print("=" * 75)
print("UNIFIED THEOREM (proved for odd, conjectured for even):")
print("=" * 75)
print()
print("For n >= 3:")
print()
print("  v_2(T(n)) = { (n-1)/2                     if n is odd    [PROVED]")
print("              { n/2 + v_2(c(n))              if n is even   [CONJECTURED]")
print()
print("where c(n) = #{odd k : 1 ≤ k < n/2, gcd(k,n) = 1}.")
print()
print("For even n = 2^a × m (m odd):")
print("  c(n) = 2^{a-2} × φ(m)  if a ≥ 2")
print("  c(n) = φ(m)/2           if a = 1")
print()
print("v_2(c(n)) counts the 2-adic content of the coprime partition count.")
print()
print("DEEP REASON: The excess for even n comes from CANCELLATION at the")
print("minimum v_2 level. There are c(n) terms with coprime 2-part partitions,")
print("each contributing an ODD coefficient (since n_odd/(ab) is odd for odd a,b).")
print("Their sum has v_2 = v_2(c(n)) by the 'sum of c odd numbers' principle.")
