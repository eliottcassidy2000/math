#!/usr/bin/env python3
"""
opus-2026-04-05-S25: Even-n v_2 analysis and odd-part factorization.

For ODD n: v_2(T(n)) = (n-1)/2 (PROVED).
For EVEN n: v_2(T(n)) ≥ n/2, with excess 0, 0, 1, 1, 1, 0, 2, 0, 2 for n=4..20.

Questions:
1. What determines the even-n excess?
2. The odd part T(n)/2^{(n-1)/2} for odd n: are 11971 and 28242289 really prime?
3. Is the n-cycle still the v_2 minimizer for even n? (No — n-cycle doesn't contribute!)
"""

from math import factorial, gcd, comb
from collections import Counter

def v2(n):
    if n == 0: return float('inf')
    v = 0
    while n % 2 == 0: v += 1; n //= 2
    return v

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    i = 5
    while i*i <= n:
        if n%i==0 or n%(i+2)==0: return False
        i += 6
    return True

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

def odd_partitions(n, max_part=None):
    if max_part is None: max_part = n
    if n == 0:
        yield ()
        return
    start = min(n, max_part)
    if start % 2 == 0: start -= 1
    for k in range(start, 0, -2):
        for rest in odd_partitions(n - k, k):
            yield (k,) + rest

def z_lambda(lam):
    counts = Counter(lam)
    z = 1
    for c in lam:
        z *= c
    for m in counts.values():
        z *= factorial(m)
    return z

def e_lambda(lam):
    k = len(lam)
    cross = sum(gcd(lam[i], lam[j]) for i in range(k) for j in range(i+1, k))
    internal = sum((c - 1) // 2 for c in lam)
    return cross + internal

# ═══════════════════════════════════════════════════════════════
# Part 1: Even-n v_2 minimizer analysis
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("EVEN n: Which partition minimizes v_2?")
print("  (For even n, the n-cycle doesn't contribute since n is even)")
print("=" * 70)
print()

for n in [4, 6, 8, 10, 12]:
    nfact = factorial(n)
    entries = []
    total = 0
    for lam in odd_partitions(n):
        z = z_lambda(lam)
        e = e_lambda(lam)
        coeff = nfact // z
        contrib = coeff * (2 ** e)
        total += contrib
        v = v2(contrib)
        entries.append((v, lam, contrib))

    entries.sort()
    T = total // nfact

    print(f"n = {n}: T({n}) = {T}, v_2 = {v2(T)}, n/2 = {n//2}, excess = {v2(T) - n//2}")
    # Show top 5 minimum v_2 terms
    for i, (v, lam, contrib) in enumerate(entries[:5]):
        oc = contrib >> v
        print(f"  v_2={v}: {lam}, odd coeff = {oc}")
    print()

# ═══════════════════════════════════════════════════════════════
# Part 2: For even n, the (n-1, 1) partition
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("EVEN n: The (n-1, 1) partition contribution")
print("=" * 70)
print()

for n in range(4, 22, 2):
    # (n-1, 1): z = (n-1)*1*1!*1! = n-1. n!/z = n!/(n-1) = n*(n-2)!
    # e = gcd(n-1,1) + (n-2)/2 + 0 = 1 + (n-2)/2 = n/2
    z = n - 1
    coeff = factorial(n) // z
    e = 1 + (n - 2) // 2  # = n/2
    contrib = coeff * (2 ** e)
    v = v2(contrib)
    v_in_T = v - v2(factorial(n))
    print(f"n={n:>3}: (n-1,1) contrib v_2 = {v}, in T: v_2 = {v_in_T} = n/2 = {n//2}")

print()

# ═══════════════════════════════════════════════════════════════
# Part 3: Even-n cancellation analysis
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("EVEN n: Terms at minimum v_2 level — cancellation check")
print("=" * 70)
print()

for n in [4, 6, 8, 10, 12, 14]:
    nfact = factorial(n)
    v2_nfact = v2(nfact)

    # Find minimum v_2 among all terms
    min_v2 = float('inf')
    term_list = []
    for lam in odd_partitions(n):
        z = z_lambda(lam)
        e = e_lambda(lam)
        contrib = (nfact // z) * (2 ** e)
        v = v2(contrib)
        term_list.append((v, lam, contrib))
        if v < min_v2:
            min_v2 = v

    # Collect all terms at minimum v_2
    min_terms = [(lam, contrib >> min_v2) for v, lam, contrib in term_list if v == min_v2]

    # Sum of odd coefficients
    odd_sum = sum(oc for _, oc in min_terms)
    v2_sum = v2(odd_sum)

    actual_v2T = v2(T_val[n])
    predicted = min_v2 + v2_sum - v2_nfact

    print(f"n={n}: minimum v_2 in sum = {min_v2}")
    for lam, oc in min_terms:
        print(f"  {lam}: odd coeff = {oc}")
    print(f"  Sum of odd coeffs = {odd_sum}, v_2 of sum = {v2_sum}")
    print(f"  Predicted v_2(T) = {min_v2} + {v2_sum} - {v2_nfact} = {predicted}")
    print(f"  Actual v_2(T) = {actual_v2T}, match = {'✓' if predicted == actual_v2T else '✗'}")
    print()

# ═══════════════════════════════════════════════════════════════
# Part 4: Odd parts T(n)/2^{(n-1)/2} for odd n
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("ODD PARTS: T(n) / 2^{(n-1)/2} for odd n")
print("=" * 70)
print()

for n in range(3, 21, 2):
    t = T_val[n]
    half = (n - 1) // 2
    odd = t >> half
    # Factor
    remaining = odd
    factors = []
    for p in range(3, min(int(remaining**0.5) + 2, 10**6)):
        while remaining % p == 0:
            factors.append(p)
            remaining //= p
    if remaining > 1:
        factors.append(remaining)
    if not factors:
        factors = [1]

    if len(factors) == 1 and factors[0] == odd and odd > 1:
        tag = "PRIME"
    else:
        tag = " × ".join(str(f) for f in factors)

    print(f"  n={n:>2}: T/{n-1}2^{half} = {odd:>20}  =  {tag}")

# Double-check the odd part sequence against OEIS: 1, 3, 57, 11971, 28242289, ...
print()
print("Odd part sequence: ", [T_val[n] >> ((n-1)//2) for n in range(3, 21, 2)])
print()

# ═══════════════════════════════════════════════════════════════
# Part 5: v_2(T(n)) for even n — attempt to find formula
# ═══════════════════════════════════════════════════════════════

print("=" * 70)
print("EVEN n: Searching for a v_2 formula")
print("=" * 70)
print()

print("n   v_2(T)  n/2  excess  v_2(n)  v_2(n-1)  v_2(n!)  Legendre bits")
print("-" * 70)
for n in range(4, 21, 2):
    t = T_val[n]
    val = v2(t)
    excess = val - n//2
    v2n = v2(n)
    v2nm1 = v2(n-1)  # always 0 since n even => n-1 odd
    v2nfact = v2(factorial(n))

    # The (n-1) cycle would give (n-2)/2 ... but it doesn't contribute
    # Instead (n-1, 1) gives n/2

    # What about looking at floor(log_2(n-1))?
    from math import log2
    lb = int(log2(n)) if n > 0 else 0

    print(f"{n:>3}  {val:>5}  {n//2:>3}  {excess:>6}  {v2n:>5}  {v2nm1:>7}  {v2nfact:>6}  log2(n)={lb}")
