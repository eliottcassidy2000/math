#!/usr/bin/env python3
"""
signed_permanent_ocf_connection.py — Explore whether S(T) can be expressed
via the independence polynomial I(Omega(T), x) at some point.

We know:
  H(T) = I(Omega(T), 2)  (OCF, proved by Grinberg-Stanley)
  S(T) = sum_P prod_{i} B[P_i][P_{i+1}]  where B[i][j] = 2A[i][j] - 1

Expanding: S(T) = sum_{k=0}^{n-1} (-1)^{n-1-k} * 2^k * D_k
where D_k counts (permutation, k-forward-edges) pairs.

Question: Is S(T)/2^{n-1} related to I(Omega(T), x) at some x?

At n=5: S/16 = H - 3*t3 = I(Omega,2) - 3*alpha_1
At n=7: S/64 has fractional part 3/4, so NOT an integer multiple of I(Omega,2).

NEW IDEA: What about I(Omega(T), -2)? Or I(Omega(T), 1)?

I(Omega(T), 1) = number of independent sets in Omega(T) (graph-theoretic invariant)
I(Omega(T), -1) = related to chromatic polynomial evaluation

Let's compute I(Omega(T), x) for several x values and compare with S(T).

Author: opus-2026-03-07-S43b
"""
from itertools import permutations, combinations
import math

def tournament_from_bits(bits, n):
    """Create adjacency matrix from bit encoding."""
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def hamiltonian_paths(A, n):
    """Count H(T) via DP."""
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            cnt = dp[(mask, v)]
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + cnt
    return sum(dp.get((full, v), 0) for v in range(n))

def signed_hp_permanent(A, n):
    """Compute S(T) = sum_P prod B[P_i][P_{i+1}] where B = 2A-1."""
    S = 0
    for P in permutations(range(n)):
        prod = 1
        for i in range(n-1):
            prod *= (2*A[P[i]][P[i+1]] - 1)
        S += prod
    return S

def find_odd_cycles(A, n):
    """Find all odd directed cycles in the tournament."""
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            # Check all cyclic orderings
            for perm in permutations(verts[1:]):
                path = (verts[0],) + perm
                is_cycle = True
                for i in range(length):
                    if not A[path[i]][path[(i+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    # Normalize: minimum vertex first
                    min_idx = path.index(min(path))
                    normalized = path[min_idx:] + path[:min_idx]
                    if normalized not in cycles:
                        cycles.append(normalized)
    return cycles

def independence_polynomial(cycles, x):
    """Compute I(Omega(T), x) where Omega(T) has cycles as vertices."""
    m = len(cycles)
    if m > 20:
        return None  # Too many cycles

    # Build adjacency: two cycles conflict if they share a vertex
    cycle_sets = [set(c) for c in cycles]

    # Enumerate independent sets
    result = 0
    for mask in range(1 << m):
        # Check independence
        bits = []
        for i in range(m):
            if mask & (1 << i):
                bits.append(i)

        independent = True
        for i in range(len(bits)):
            for j in range(i+1, len(bits)):
                if cycle_sets[bits[i]] & cycle_sets[bits[j]]:
                    independent = False
                    break
            if not independent:
                break

        if independent:
            result += x ** len(bits)

    return result

# Test at n=5
print("=== n=5: S(T) vs I(Omega, x) at various x ===")
print(f"{'bits':>6} {'H':>5} {'S':>6} {'S/16':>5} {'I(2)':>5} {'I(-2)':>6} {'I(1)':>5} {'I(-1)':>6} {'t3':>3} {'alpha1':>6}")

n = 5
m = n*(n-1)//2
for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    H = hamiltonian_paths(A, n)
    S = signed_hp_permanent(A, n)
    cycles = find_odd_cycles(A, n)
    t3 = sum(1 for c in cycles if len(c) == 3)
    alpha1 = len(cycles)

    I2 = independence_polynomial(cycles, 2)
    Im2 = independence_polynomial(cycles, -2)
    I1 = independence_polynomial(cycles, 1)
    Im1 = independence_polynomial(cycles, -1)

    # Only print one representative per isomorphism class (use H as proxy)
    # Actually print first few
    if bits < 50 or bits % 200 == 0:
        print(f"{bits:6d} {H:5d} {S:6d} {S//16:5d} {I2:5d} {Im2:6d} {I1:5d} {Im1:6d} {t3:3d} {alpha1:6d}")

# Now look for a formula
print("\n=== Checking S/16 = I(Omega, 2) - 3*alpha_1 at n=5 ===")
all_match = True
for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    H = hamiltonian_paths(A, n)
    S = signed_hp_permanent(A, n)
    cycles = find_odd_cycles(A, n)
    alpha1 = len(cycles)

    I2 = independence_polynomial(cycles, 2)
    predicted = I2 - 3*alpha1

    if S//16 != predicted:
        print(f"  MISMATCH at bits={bits}: S/16={S//16}, predicted={predicted}")
        all_match = False
        break

if all_match:
    print("  CONFIRMED: S/16 = I(Omega,2) - 3*alpha_1 for ALL n=5 tournaments")

# Check: S/16 = H - 3*t3 (from the skeleton document)
print("\n=== Checking S/16 = H - 3*t3 at n=5 ===")
all_match = True
for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    H = hamiltonian_paths(A, n)
    S = signed_hp_permanent(A, n)
    cycles = find_odd_cycles(A, n)
    t3 = sum(1 for c in cycles if len(c) == 3)

    if S//16 != H - 3*t3:
        print(f"  MISMATCH: S/16={S//16}, H-3t3={H - 3*t3}")
        all_match = False
        break

if all_match:
    print("  CONFIRMED: S/16 = H - 3*t3 for ALL n=5 tournaments")

# Now check: what is H - 3*t3 in terms of independence polynomial?
# H = I(Omega, 2) = 1 + 2*alpha1 + 4*alpha2
# So H - 3*alpha1 = 1 - alpha1 + 4*alpha2
# = I(Omega, 2) - 3*alpha1
# But alpha1 = t3 + t5 at n=5
# H - 3*t3 = H - 3*alpha1 + 3*t5
# So H - 3*t3 = 1 - alpha1 + 4*alpha2 + 3*t5

# At n=5, t5 = 0 or 1 depending on whether there's a 5-cycle
# So if t5 = 0: S/16 = 1 - alpha1 + 4*alpha2
# if t5 = 1: S/16 = 4 - alpha1 + 4*alpha2 = 1 - (alpha1-3) + 4*alpha2

print("\n=== Detailed: alpha1, t3, t5, alpha2, S/16 at n=5 ===")
seen = set()
for bits in range(1 << m):
    A = tournament_from_bits(bits, n)
    H = hamiltonian_paths(A, n)
    S = signed_hp_permanent(A, n)
    cycles = find_odd_cycles(A, n)
    t3 = sum(1 for c in cycles if len(c) == 3)
    t5 = sum(1 for c in cycles if len(c) == 5)
    alpha1 = len(cycles)

    # Compute alpha2 = number of independent pairs
    cycle_sets = [set(c) for c in cycles]
    alpha2 = sum(1 for i in range(len(cycles)) for j in range(i+1, len(cycles))
                 if not (cycle_sets[i] & cycle_sets[j]))

    key = (H, t3, t5, alpha1, alpha2, S//16)
    if key not in seen:
        seen.add(key)
        # Check: S/16 = 1 + 2*alpha1 + 4*alpha2 - 3*alpha1 = 1 - alpha1 + 4*alpha2
        formula = 1 - alpha1 + 4*alpha2
        match = "OK" if S//16 == formula else "FAIL"
        print(f"  H={H:3d} t3={t3} t5={t5} a1={alpha1:2d} a2={alpha2:2d} S/16={S//16:3d} | 1-a1+4*a2={formula:3d} {match}")

# Now the key question: Is 1 - alpha1 + 4*alpha2 = I(Omega, x) for some x?
# I(Omega, x) = 1 + alpha1*x + alpha2*x^2 + ...
# So I(Omega, -1) = 1 - alpha1 + alpha2 - alpha3 + ...
# And 4*I(Omega, -1) = 4 - 4*alpha1 + 4*alpha2 - ...
# Our formula: 1 - alpha1 + 4*alpha2
# This is NOT I(Omega, x) at any single x, since the x and x^2 coefficients
# would need different values.

# But: I(Omega, 2) - 3*I(Omega, 1) + 3*1 = (1+2a1+4a2) - 3*(1+a1+a2) + 3
#    = 1+2a1+4a2 - 3-3a1-3a2 + 3 = 1 - a1 + a2
# Not quite right either.

# Let's try: I(Omega, 2) - 3*alpha1 = 1 + 2a1 + 4a2 - 3a1 = 1 - a1 + 4a2 ✓
# So S/16 = I(Omega, 2) - 3*alpha1 (at n=5)

# But alpha1 = I'(Omega, 0) = coefficient of x in I(Omega, x)
# So S/16 = I(Omega, 2) - 3*I'(Omega, 0)

# Or: S/16 = I(Omega, 2) - 3*(I(Omega,1) - 1) + 3*(alpha2 - alpha3 + ...)
# Hmm, not clean.

# What about: define J(x) = I(Omega, x) at even n, and check S/2^{n-1} = J(2) - ...?

# Let's try n=7 to see if the pattern extends
print("\n=== n=7: S(T) analysis (sampling) ===")
import random
n = 7
m = n*(n-1)//2
random.seed(42)

print(f"{'H':>5} {'S':>8} {'S/64':>7} {'I(2)':>5} {'a1':>4} {'I2-3a1':>7} {'frac':>5}")
for trial in range(20):
    bits = random.getrandbits(m)
    A = tournament_from_bits(bits, n)
    H = hamiltonian_paths(A, n)
    S = signed_hp_permanent(A, n)
    cycles = find_odd_cycles(A, n)
    t3 = sum(1 for c in cycles if len(c) == 3)
    alpha1 = len(cycles)

    I2 = independence_polynomial(cycles, 2)

    s_over = S / 64
    formula = I2 - 3*alpha1

    print(f"{H:5d} {S:8d} {s_over:7.2f} {I2:5d} {alpha1:4d} {formula:7d} {s_over - formula:7.2f}")
