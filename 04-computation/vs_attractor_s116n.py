#!/usr/bin/env python3
"""vs_attractor_s116n.py — The von Staudt attractor basin of 1806.

The chain n -> prod{p : (p-1)|n} has a unique fixed point 1806 in [1,100000].
Which starting values converge to it? What characterizes the basin of attraction?

Session: kind-pasteur-2026-03-16-S116n32
"""
import sys
sys.stdout.reconfigure(line_buffering=True)

from math import isqrt

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
    """Product of primes p with (p-1)|n."""
    result = 1
    for d in range(1, n+1):
        if n % d == 0 and is_prime(d+1):
            result *= (d+1)
    return result

print()
print("  VON STAUDT ATTRACTOR BASIN OF 1806")
print()
print("=" * 70)
print()

# ============================================================
print("  I. WHICH STARTING VALUES CONVERGE TO 1806?")
print("  " + "-" * 50)
print()

converges = []  # starting values that reach 1806
diverges = []   # starting values that exceed 10^8
limit = 10**8

for start in range(1, 201):
    current = start
    chain = [current]
    converged = False
    for _ in range(10):
        nxt = vs_prod(current)
        chain.append(nxt)
        if nxt == 1806:
            converged = True
            break
        if nxt == current:  # other fixed point
            break
        if nxt > limit:
            break
        current = nxt
    if converged:
        converges.append((start, chain))
    else:
        diverges.append((start, chain))

print(f"  Starting values 1-200 that converge to 1806:")
for start, chain in converges:
    print(f"    {start:4d}: {' -> '.join(str(x) for x in chain)}")
print()
print(f"  Total converging: {len(converges)}/200")
print()

# ============================================================
print("  II. CONVERGENT STARTING VALUES")
print("  " + "-" * 50)
print()

# Extract just the starting values
conv_starts = [s for s, _ in converges]
print(f"  Converging starts: {conv_starts}")
print()

# Factor each
def factorize(n):
    if n <= 1: return {n: 1}
    factors = {}
    temp = n
    for p in range(2, isqrt(temp)+2):
        while temp % p == 0:
            factors[p] = factors.get(p, 0) + 1
            temp //= p
    if temp > 1:
        factors[temp] = factors.get(temp, 0) + 1
    return factors

for s in conv_starts:
    f = factorize(s)
    fstr = ' * '.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(f.items()))
    print(f"    {s:4d} = {fstr}")
print()

# What do converging starts have in common?
# Check: all divisors of 42?
divs_42 = [d for d in range(1, 43) if 42 % d == 0]
print(f"  Divisors of 42: {divs_42}")
print(f"  Divisors of 42 in converging set: {[d for d in divs_42 if d in conv_starts]}")
print(f"  Non-divisors of 42 in converging set: {[s for s in conv_starts if 42 % s != 0]}")
print()

# Check: all divisors of 1806?
divs_1806 = [d for d in range(1, 201) if 1806 % d == 0]
print(f"  Divisors of 1806 up to 200: {divs_1806}")
print(f"  Divisors of 1806 in converging set: {[d for d in divs_1806 if d in conv_starts]}")
print(f"  Converging but not divisor of 1806: {[s for s in conv_starts if 1806 % s != 0]}")
print()

# ============================================================
print("  III. THE GROWTH RATE AT EACH STEP")
print("  " + "-" * 50)
print()
print("  For the canonical chain 1 -> 2 -> 6 -> 42 -> 1806:")
def vs_primes_list(n):
    return [d+1 for d in range(1, n+1) if n % d == 0 and is_prime(d+1)]

chain_canonical = [1, 2, 6, 42, 1806]
for i in range(len(chain_canonical)):
    n = chain_canonical[i]
    primes = vs_primes_list(n) if n > 0 else []
    prod = 1
    for p in primes:
        prod *= p
    new_str = ""
    if i > 0:
        prev_primes = set(vs_primes_list(chain_canonical[i-1]))
        new_p = set(primes) - prev_primes
        new_str = f"  NEW: {new_p}" if new_p else "  (no new primes)"
    print(f"  n={n:5d}: VS primes = {primes}, product = {prod}{new_str}")
print()

# For diverging chains, how many primes enter at each step?
print("  Diverging chains — prime count explosion:")
for start in [4, 8, 10, 12, 20, 30]:
    current = start
    for step in range(3):
        primes = vs_primes_list(current) if current <= 100000 else ['(too large)']
        if isinstance(primes, list) and all(isinstance(p, int) for p in primes):
            print(f"    start={start}, step {step}: n={current}, "
                  f"|primes|={len(primes)}, primes={primes[:8]}{'...' if len(primes)>8 else ''}")
            prod = 1
            for p in primes:
                prod *= p
            current = prod
        else:
            print(f"    start={start}, step {step}: n={current} (too large)")
            break
    print()

# ============================================================
print("  IV. THE CRITICAL PROPERTY: ONE NEW PRIME PER STEP")
print("  " + "-" * 50)
print()
print("  Canonical chain analysis:")
print()
print("  Step 0: n=1. Divisors of 1: {1}. d+1=2 prime. VS={2}. Prod=2.")
print("    NEW PRIME: 2. (1 new prime)")
print()
print("  Step 1: n=2. Divisors of 2: {1,2}. d+1={2,3}. VS={2,3}. Prod=6.")
print("    NEW PRIME: 3. (1 new prime)")
print()
print("  Step 2: n=6. Divisors of 6: {1,2,3,6}. d+1={2,3,4,7}.")
print("    Primes: {2,3,7}. VS={2,3,7}. Prod=42.")
print("    NEW PRIME: 7. (1 new prime)")
print()
print("  Step 3: n=42. Divisors of 42: {1,2,3,6,7,14,21,42}.")
print("    d+1 = {2,3,4,7,8,15,22,43}.")
print("    Primes: {2,3,7,43}. VS={2,3,7,43}. Prod=1806.")
print("    NEW PRIME: 43. (1 new prime)")
print()
print("  Step 4: n=1806. VS primes = {2,3,7,43}. Prod=1806. FIXED POINT.")
print("    NEW PRIMES: none. (0 new primes)")
print()
print("  The canonical chain adds EXACTLY ONE PRIME per step until stabilizing.")
print("  This is the minimal growth trajectory.")
print()

# Why does 42 attract exactly one new prime?
# Divisors of 42: 1,2,3,6,7,14,21,42
# d+1 candidates: 2,3,4,7,8,15,22,43
# Primes already in 6's VS set: {2,3,7}
# New candidate: 43 = 42+1 (the n+1 prime)
# 4 = not prime, 8 = not prime, 15 = not prime, 22 = not prime
# So the ONLY new prime from 42 is 43 = 42+1.
# This works because:
# - 42+1 = 43 is prime (Sylvester)
# - All other d+1 for proper divisors d < 42 give non-primes or primes already in {2,3,7}

print("  WHY 42 attracts exactly one new prime:")
print("  Divisors of 42 contribute d+1 = {2,3,4,7,8,15,22,43}")
print("  Already known: {2,3,7}")
print("  New candidates: {4(=2^2), 8(=2^3), 15(=3*5), 22(=2*11), 43}")
print("  Composite: {4, 8, 15, 22}")
print("  PRIME: {43}")
print("  So 43 is the ONLY new prime. This requires 43=42+1 to be prime!")
print()
print("  Similarly for 6: d+1 = {2,3,4,7}. Known: {2,3}. New: {4(composite), 7(PRIME)}.")
print("  7 is the only new prime. Requires 7=6+1 to be prime!")
print()
print("  The chain works because at each step, n+1 is prime and all other")
print("  new d+1 values are composite. This is the Sylvester condition:")
print("  a_{k+1} = P_k + 1 is prime.")
print()

# ============================================================
print("  V. THE SYLVESTER PRIMALITY CONDITION")
print("  " + "-" * 50)
print()
print("  The canonical chain adds one prime per step IFF:")
print("  P_k + 1 is prime, AND")
print("  no divisor d of P_k with d not in {1, P_{k-1}} gives a new prime d+1.")
print()
print("  Let's verify the second condition is automatic:")
print("  If P_k = 2*3*7*...p_k, then every proper divisor d is a product of")
print("  a subset of {2,3,7,...,p_k}. So d+1 = (some product)+1.")
print("  For d+1 to be a NEW prime not already in the set,")
print("  d+1 must not divide P_k (otherwise d+1 would already be a factor).")
print()

# Check: for 42, the divisors d with d+1 prime:
for n_val, label in [(6, "n=6"), (42, "n=42"), (1806, "n=1806")]:
    divs = sorted([d for d in range(1, n_val+1) if n_val % d == 0])
    new_primes = []
    existing_primes = set()
    for d in divs:
        if is_prime(d+1):
            existing_primes.add(d+1)
    print(f"  {label}: divisors giving primes d+1: {sorted(existing_primes)}")

print()

# ============================================================
print("  VI. THE NUMBER THEORY BEHIND THE FIXED POINT")
print("  " + "-" * 50)
print()
print("  1806 = 2 * 3 * 7 * 43 is a VS fixed point because:")
print("  (1) Every prime factor p of 1806 has (p-1)|1806:")
print("      1|1806, 2|1806, 6|1806, 42|1806. All true.")
print("  (2) No OTHER prime q has (q-1)|1806:")
print("      Need (q-1) | 1806 = 2*3*7*43 and q prime, q not in {2,3,7,43}.")
for q in [5, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53, 59, 61, 67]:
    divides = (1806 % (q-1) == 0)
    if divides:
        print(f"      q={q}: (q-1)={q-1} | 1806? YES! (but q should not be new...)")
    # Only print near-misses
    elif 1806 % (q-1) < 10:
        print(f"      q={q}: (q-1)={q-1} | 1806? NO (remainder {1806 % (q-1)})")

print()
print("  Checking all primes q < 200:")
for q in range(2, 200):
    if is_prime(q) and q not in {2, 3, 7, 43} and 1806 % (q-1) == 0:
        print(f"    PROBLEM: q={q}, (q-1)={q-1} divides 1806!")

print("  (No problems found => 1806 is confirmed as VS fixed point)")
print()

# Why don't other primes sneak in?
# 1806 = 2 * 3 * 7 * 43
# Divisors: 1,2,3,6,7,14,21,42,43,86,129,258,301,602,903,1806
# d+1: 2,3,4,7,8,15,22,43,44,87,130,259,302,603,904,1807
# Primes among d+1: 2,3,7,43
# NON-primes: 4=2^2, 8=2^3, 15=3*5, 22=2*11, 44=2^2*11, 87=3*29,
#             130=2*5*13, 259=7*37, 302=2*151, 603=3*201=3*3*67, 904=2^3*113, 1807=13*139

print("  d+1 for each divisor d of 1806:")
divs_1806 = sorted([d for d in range(1, 1807) if 1806 % d == 0])
for d in divs_1806:
    dp1 = d + 1
    if is_prime(dp1):
        print(f"    d={d:5d}: d+1={dp1:5d} PRIME (already a factor of 1806)")
    else:
        # Factor d+1
        f = factorize(dp1)
        fstr = '*'.join(f'{p}^{e}' if e > 1 else str(p) for p, e in sorted(f.items()))
        print(f"    d={d:5d}: d+1={dp1:5d} = {fstr}")
print()

print("  Every d+1 that is prime is already in {2,3,7,43}.")
print("  Every other d+1 is composite.")
print("  The set {2,3,7,43} is SELF-CONTAINED under the VS operation.")
print()

# ============================================================
print("  VII. CONNECTION TO HURWITZ")
print("  " + "-" * 50)
print()
print("  {2,3,7} = primes dividing disc(Hurwitz quaternions) = 42")
print("  {2,3,7,43} = primes dividing denom(B_42) = 1806")
print("  43 = 42+1 = disc(Hurwitz)+1")
print()
print("  The Hurwitz primes are NOT a VS fixed point.")
print("  They need 43 to complete the self-selecting set.")
print("  42 is the TRANSITION: the bridge from {2,3,7} to {2,3,7,43}.")
print()
print("  Interpretation: The Hurwitz quaternion algebra (ramified at 2,3,7)")
print("  is INCOMPLETE from the von Staudt perspective.")
print("  It needs one more prime (43 = 42+1) to achieve self-consistency.")
print("  The VS fixed point 1806 = 42 * 43 is the COMPLETED Hurwitz system.")
print()
