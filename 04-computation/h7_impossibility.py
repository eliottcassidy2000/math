#!/usr/bin/env python3
"""
Why is H=7 impossible? Analyze via OCF.

H(T) = I(Omega(T), 2) = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
where alpha_k = number of independent sets of size k in Omega(T).

For H=7: need 2*alpha_1 + 4*alpha_2 = 6
Only option with small alpha: alpha_1=3, alpha_2=0
(alpha_1=1, alpha_2=1 impossible since 1 vertex can't have independent pair)

So: need exactly 3 directed odd cycles, all pairwise sharing a vertex.

Question: can a tournament have exactly 3 directed odd cycles?

kind-pasteur-2026-03-06-S21
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from itertools import combinations

def count_directed_odd_cycles(T, n):
    """Count all directed odd cycles in tournament T."""
    total = 0
    for k in range(3, n+1, 2):
        for combo in combinations(range(n), k):
            verts = list(combo)
            # Count directed Hamiltonian cycles on this subset
            dp = {}
            dp[(1, 0)] = 1
            for mask in range(1, 1 << k):
                if not (mask & 1):
                    continue
                for vi in range(k):
                    if not (mask & (1 << vi)):
                        continue
                    c = dp.get((mask, vi), 0)
                    if c == 0:
                        continue
                    for ui in range(k):
                        if mask & (1 << ui):
                            continue
                        if T[verts[vi]][verts[ui]]:
                            key = (mask | (1 << ui), ui)
                            dp[key] = dp.get(key, 0) + c
            full = (1 << k) - 1
            for vi in range(1, k):
                c = dp.get((full, vi), 0)
                if c > 0 and T[verts[vi]][verts[0]]:
                    total += c
    return total

# Collect alpha_1 distribution at each n
for n in range(3, 7):
    m = n*(n-1)//2
    alpha1_counts = {}
    h_by_alpha1 = {}
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        a1 = count_directed_odd_cycles(T, n)
        alpha1_counts[a1] = alpha1_counts.get(a1, 0) + 1
        if a1 not in h_by_alpha1:
            h_by_alpha1[a1] = set()
        h_by_alpha1[a1].add(h)

    print(f"\nn={n}: alpha_1 distribution")
    for a1 in sorted(alpha1_counts.keys()):
        h_vals = sorted(h_by_alpha1[a1])
        print(f"  alpha_1={a1}: count={alpha1_counts[a1]}, H values={h_vals}")

    # Check: is alpha_1=3 achievable?
    if 3 in alpha1_counts:
        print(f"  >>> alpha_1=3 IS achievable at n={n}")
    else:
        print(f"  >>> alpha_1=3 NOT achievable at n={n}")

# Key insight: if alpha_1 is always even, then H mod 4 is constrained
print("\n\n=== KEY OBSERVATION ===")
print("If alpha_1 is always even, then 2*alpha_1 is always 0 mod 4")
print("So H = 1 + 0 + 4*alpha_2 + ... = 1 mod 4 always")
print("But H=3 = 3 mod 4 exists, so alpha_1 CAN be odd!")
print("Let's check parity of alpha_1...")

for n in range(3, 7):
    m = n*(n-1)//2
    odd_count = 0
    total = 0
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        a1 = count_directed_odd_cycles(T, n)
        total += 1
        if a1 % 2 == 1:
            odd_count += 1
    print(f"  n={n}: alpha_1 odd in {odd_count}/{total} = {100*odd_count/total:.1f}%")

print("\nDone.")
