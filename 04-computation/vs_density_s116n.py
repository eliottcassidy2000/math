#!/usr/bin/env python3
"""vs_density_s116n.py — Density of the von Staudt convergent basin.

Which n converge to 1806 under the VS map n -> prod{p : (p-1)|n}?
Key insight: ALL odd primes p > 2 satisfy vs_prod(p) = 2 (since only (p-1)|p is d=1,
and d+1=2 is prime, plus d=p gives d+1=p+1 which is even hence composite for p>2).
Wait — that's wrong. Let me think more carefully.

For prime p: divisors of p are {1, p}. So d+1 values are {2, p+1}.
2 is always prime. p+1 is prime iff p is a twin prime's smaller member.
So vs_prod(p) = 2 if p+1 is NOT prime, and vs_prod(p) = 2*(p+1) if p+1 IS prime.
For p=2: vs_prod(2) = {d+1 : d|2 and d+1 prime} = {2, 3}. Product = 6.
For p=3: divisors {1,3}, d+1={2,4}. 4 not prime. vs_prod(3) = 2.
For p=5: divisors {1,5}, d+1={2,6}. 6 not prime. vs_prod(5) = 2.
For p=7: divisors {1,7}, d+1={2,8}. 8 not prime. vs_prod(7) = 2.
For p=11: divisors {1,11}, d+1={2,12}. 12 not prime. vs_prod(11) = 2.
For p=13: divisors {1,13}, d+1={2,14}. 14 not prime. vs_prod(13) = 2.

So ALL primes p >= 3 have vs_prod(p) = 2, hence converge via 2->6->42->1806.
p=2 has vs_prod(2) = 6, converges via 6->42->1806.

This means the basin includes ALL primes! And the density of primes is ~1/ln(n).
But many composites also converge. What's the total density?

Session: kind-pasteur-2026-03-16-S116n32
"""
import sys
sys.stdout.reconfigure(line_buffering=True)

from math import isqrt, log

def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0 or n % 3 == 0: return False
    d = 5
    while d * d <= n:
        if n % d == 0 or n % (d+2) == 0: return False
        d += 6
    return True

def vs_prod(n):
    result = 1
    for d in range(1, n+1):
        if n % d == 0 and is_prime(d+1):
            result *= (d+1)
    return result

print()
print("  DENSITY OF THE VON STAUDT CONVERGENT BASIN")
print()
print("=" * 70)
print()

# ============================================================
print("  I. WHY ALL PRIMES CONVERGE")
print("  " + "-" * 50)
print()
print("  For prime p >= 3: divisors of p = {1, p}.")
print("  d+1 values: {2, p+1}.")
print("  2 is always prime.")
print("  p+1 is EVEN (since p odd), hence composite for p >= 3.")
print("  So vs_prod(p) = 2 for ALL primes p >= 3.")
print("  Then: 2 -> 6 -> 42 -> 1806.")
print()
print("  For p=2: divisors = {1,2}. d+1 = {2,3}. vs_prod(2) = 6.")
print("  Then: 6 -> 42 -> 1806.")
print()
print("  ALL PRIMES converge to 1806 in at most 5 steps.")
print()

# ============================================================
print("  II. WHICH COMPOSITES CONVERGE?")
print("  " + "-" * 50)
print()

# A composite n converges if vs_prod(n) either:
# (a) equals a value that converges (recursively), or
# (b) equals 1806 directly

# The key question: when does vs_prod(n) stay small enough to converge?
# vs_prod(n) = product of ALL primes p with (p-1)|n.
# If n has many divisors, vs_prod(n) can explode.

# Count converging vs diverging in ranges
limit = 10**8  # Upper bound for detecting divergence
max_check = 1000

converge_count = [0] * (max_check + 1)  # 1 if converges, 0 if not

for n in range(1, max_check + 1):
    current = n
    converged = False
    for _ in range(12):  # max iterations
        if current == 1806:
            converged = True
            break
        if current > limit:
            break
        nxt = vs_prod(current)
        if nxt == current and nxt != 1806:  # stuck at non-1806 fixed point
            break
        current = nxt
    converge_count[n] = 1 if converged else 0

# Density by range
print(f"  Convergence density up to {max_check}:")
for end in [50, 100, 200, 500, 1000]:
    total = sum(converge_count[1:end+1])
    print(f"    [1, {end:4d}]: {total}/{end} = {total/end:.4f}")
print()

# Identify the NON-converging values
non_converging = [n for n in range(1, max_check+1) if converge_count[n] == 0]
print(f"  Non-converging values up to {max_check}: {len(non_converging)}")
print(f"  First 40: {non_converging[:40]}")
print()

# Factor the non-converging values
print("  Factoring first 20 non-converging values:")
def factorize(n):
    factors = {}
    temp = n
    for p in range(2, isqrt(temp)+2):
        while temp % p == 0:
            factors[p] = factors.get(p, 0) + 1
            temp //= p
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

for n in non_converging[:20]:
    f = factorize(n)
    fstr = ' * '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(f.items()))
    # Count divisors
    tau = 1
    for e in f.values():
        tau *= (e + 1)
    print(f"    {n:4d} = {fstr:20s}  tau={tau:3d}")
print()

# ============================================================
print("  III. WHAT DETERMINES DIVERGENCE?")
print("  " + "-" * 50)
print()

# Hypothesis: n diverges if its number of divisors tau(n) is large enough
# that vs_prod(n) picks up too many primes.

# The critical observation: n diverges if vs_prod(n) has more prime factors
# than n itself, causing unbounded growth.

# Check: do all non-converging have even tau?
even_tau = sum(1 for n in non_converging if len([d for d in range(1, n+1) if n % d == 0]) % 2 == 0)
print(f"  Non-converging with even tau: {even_tau}/{len(non_converging)}")

# Check: do non-converging have 4 | n?
div_by_4 = sum(1 for n in non_converging if n % 4 == 0)
div_by_6 = sum(1 for n in non_converging if n % 6 == 0)
div_by_8 = sum(1 for n in non_converging if n % 8 == 0)
div_by_10 = sum(1 for n in non_converging if n % 10 == 0)
div_by_12 = sum(1 for n in non_converging if n % 12 == 0)
print(f"  Non-converging with 4|n: {div_by_4}/{len(non_converging)} ({div_by_4/len(non_converging):.3f})")
print(f"  Non-converging with 6|n: {div_by_6}/{len(non_converging)} ({div_by_6/len(non_converging):.3f})")
print(f"  Non-converging with 8|n: {div_by_8}/{len(non_converging)} ({div_by_8/len(non_converging):.3f})")
print(f"  Non-converging with 10|n: {div_by_10}/{len(non_converging)} ({div_by_10/len(non_converging):.3f})")
print(f"  Non-converging with 12|n: {div_by_12}/{len(non_converging)} ({div_by_12/len(non_converging):.3f})")
print()

# What's special about 4? vs_prod(4) = vs_prod(2^2)
# Divisors of 4: {1,2,4}. d+1 = {2,3,5}. All prime!
# vs_prod(4) = 2*3*5 = 30. Now 30 has many divisors...
# vs_prod(30): divisors of 30 = {1,2,3,5,6,10,15,30}
# d+1 = {2,3,4,6,7,11,16,31}. Primes: {2,3,7,11,31}.
# vs_prod(30) = 2*3*7*11*31 = 14322.
# This explodes because 30 attracted 5 primes at once.

print("  KEY: n=4 attracts prime 5 (since 4+1=5 is prime).")
print("  vs_prod(4) = 2*3*5 = 30.")
print("  30 has 8 divisors, attracting 5 primes: {2,3,7,11,31}.")
print("  vs_prod(30) = 14322. This explodes.")
print()
print("  The trigger: if ANY divisor d of n gives d+1 = 5, i.e., 4|n,")
print("  then 5 enters the product, leading to much larger next values.")
print("  Since 5 is NOT in the Hurwitz set {2,3,7}, its presence")
print("  opens up new branches that never converge to 1806.")
print()

# Verify: is 4|n necessary for divergence?
non_conv_not_div4 = [n for n in non_converging if n % 4 != 0]
print(f"  Non-converging with 4 NOT dividing n: {len(non_conv_not_div4)}")
if non_conv_not_div4:
    print(f"    Values: {non_conv_not_div4[:20]}")
print()

# ============================================================
print("  IV. THE COMPLETE CHARACTERIZATION")
print("  " + "-" * 50)
print()

# Check if n converges iff 4 does not divide n AND vs_prod(n) converges
# Actually: odd numbers always converge (either prime->2->... or composite
# with only odd divisors, whose d+1 are even hence mostly composite).

# Check: do ALL odd n converge?
odd_non_conv = [n for n in non_converging if n % 2 == 1]
print(f"  Odd non-converging values up to {max_check}: {odd_non_conv[:20]}")
print(f"  Total: {len(odd_non_conv)}")
print()

# For odd n: all divisors d are odd, so d+1 is even.
# d+1 prime iff d+1=2 iff d=1.
# So vs_prod(odd n) = 2 for ALL odd n > 1!
# Then 2 -> 6 -> 42 -> 1806.
print("  For ODD n > 1: every divisor d is odd, so d+1 is even.")
print("  The only even prime is 2, from d=1.")
print("  So vs_prod(odd n) = 2 for ALL odd n > 1.")
print("  => ALL odd integers converge to 1806.")
print()

# For even n = 2m: n has divisor 2, so d+1=3 enters.
# n has divisor 1, so d+1=2 enters.
# What else? Depends on other divisors.
# If 4|n, then d=4 gives d+1=5 (prime!), which enters.
# If 4 does not divide n (n = 2*odd), then d=2 gives d+1=3.
# Other even divisors: d=2k, d+1=2k+1. This could be prime.

# For n = 2*(odd prime q): divisors = {1, 2, q, 2q}
# d+1 = {2, 3, q+1, 2q+1}
# q+1 is even (q odd), so composite unless q+1=2 (q=1, not prime).
# 2q+1: prime? Depends on q.
# If 2q+1 is NOT prime: vs_prod(2q) = 2*3 = 6 -> converges!
# If 2q+1 IS prime (Sophie Germain): vs_prod(2q) = 2*3*(2q+1) -> ???

print("  For n = 2*p (p odd prime):")
print("  Divisors: {1, 2, p, 2p}. d+1 = {2, 3, p+1, 2p+1}.")
print("  p+1 even -> composite. 2p+1 may be prime (Sophie Germain).")
print()

# Check Sophie Germain cases
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]:
    if not is_prime(p): continue
    n = 2*p
    vsp = vs_prod(n)
    sg = is_prime(2*p+1)
    conv = converge_count[n] if n <= max_check else "?"
    print(f"    p={p:3d}, n=2p={n:4d}: 2p+1={2*p+1:4d} "
          f"({'prime (SG)' if sg else 'composite':12s}), "
          f"vs_prod={vsp:8d}, converges={'yes' if conv==1 else 'no' if conv==0 else '?'}")
print()

# n=2*p where 2p+1 is prime:
# vs_prod(2p) = 2*3*(2p+1) = 6*(2p+1)
# Then we need vs_prod(6*(2p+1)) to converge.
# 6*(2p+1) has factors including 6 and (2p+1).
# For p=29: vs_prod(58) = 2*3*59 = 354. vs_prod(354) = ?
print("  Checking: n=58=2*29, 2*29+1=59 (prime).")
print(f"  vs_prod(58) = {vs_prod(58)}")
print(f"  vs_prod(354) = {vs_prod(354)}")
print(f"  354 -> {vs_prod(354)} -> converges? {converge_count[354] if 354 <= max_check else '?'}")
print()

# ============================================================
print("  V. DENSITY SUMMARY")
print("  " + "-" * 50)
print()

# The converging set includes:
# - ALL odd integers (density 1/2)
# - Even integers 2m where vs_prod(2m) converges
# - In particular: 2p for all primes p with 2p+1 composite (density ~1/2 of primes)

total_conv = sum(converge_count[1:max_check+1])
total_odd_conv = sum(1 for n in range(1, max_check+1, 2) if converge_count[n] == 1)
total_even_conv = sum(1 for n in range(2, max_check+1, 2) if converge_count[n] == 1)
total_odd = (max_check + 1) // 2
total_even = max_check // 2

print(f"  Up to {max_check}:")
print(f"    Odd converging: {total_odd_conv}/{total_odd} = {total_odd_conv/total_odd:.4f}")
print(f"    Even converging: {total_even_conv}/{total_even} = {total_even_conv/total_even:.4f}")
print(f"    Total: {total_conv}/{max_check} = {total_conv/max_check:.4f}")
print()

# Non-converging even numbers up to 200
non_conv_even = [n for n in non_converging if n % 2 == 0 and n <= 200]
print(f"  Non-converging even numbers up to 200: {non_conv_even}")
print()

# Pattern: these all have 4|n
has_4 = [n for n in non_conv_even if n % 4 == 0]
no_4 = [n for n in non_conv_even if n % 4 != 0]
print(f"  With 4|n: {has_4[:20]}")
print(f"  Without 4|n: {no_4[:20]}")
print()

# So the conjecture is: n converges iff n is odd OR (n=2m and vs_prod(2m) converges)
# And the main obstruction is 4|n (which brings in the non-Hurwitz prime 5).

print("  CONJECTURE: n converges to 1806 iff the VS orbit of n")
print("  never encounters a value divisible by 4.")
print("  Equivalently: the prime 5 never enters the VS prime set along the orbit.")
print()
print("  The density of convergent n is approximately:")
print(f"  {total_conv/max_check:.3f} (empirically, up to {max_check})")
print()

# Check if ALL multiples of 4 diverge
mult4_conv = sum(1 for n in range(4, max_check+1, 4) if converge_count[n] == 1)
print(f"  Multiples of 4 that converge: {mult4_conv}/{max_check//4}")
if mult4_conv > 0:
    conv4 = [n for n in range(4, max_check+1, 4) if converge_count[n] == 1]
    print(f"    Examples: {conv4[:10]}")
print()
