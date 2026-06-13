#!/usr/bin/env python3
"""
Converse of Redei's theorem: which odd integers arise as H(T)?

Redei says H(T) is always odd. The converse asks: for which odd k
does there exist a tournament T with H(T) = k?

Known: H=7 and H=21 are NOT achievable (mentioned in web search).
Let's verify this computationally and find all gaps.

kind-pasteur-2026-03-06-S21
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count

# Collect all achievable H values at each n
print("=" * 60)
print("CONVERSE OF REDEI: ACHIEVABLE H VALUES")
print("=" * 60)

max_H_by_n = {}
all_H_values = set()

for n in range(2, 9):  # n=8 is 2^28, too large. Do n=2..7
    m = n * (n-1) // 2
    h_values = set()
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        h_values.add(h)
        all_H_values.add(h)
    max_H_by_n[n] = max(h_values)
    h_sorted = sorted(h_values)
    print(f"\nn={n}: {len(h_values)} distinct H values, range [{min(h_values)}, {max(h_values)}]")
    print(f"  Values: {h_sorted}")

# Find gaps: odd numbers not achieved
max_val = max(all_H_values)
achieved_odd = {h for h in all_H_values if h % 2 == 1}
all_odd = {k for k in range(1, max_val + 1) if k % 2 == 1}
gaps = sorted(all_odd - achieved_odd)

print(f"\n{'='*60}")
print(f"Odd integers 1 to {max_val} NOT achieved as H(T) for n <= 7:")
print(f"Total gaps: {len(gaps)}")
if len(gaps) <= 50:
    print(f"Gaps: {gaps}")
else:
    print(f"First 30: {gaps[:30]}")
    print(f"Last 10: {gaps[-10:]}")

# Check which small gaps are "permanent" (not achievable at ANY n)
# H(T) for n-vertex tournament can only take certain values
# At n=2: H=1 only
# At n=3: H={1,3}
# At n=4: H={1,3,5}
# At n=5: ...
# Key question: is 7 achievable? Is 21 achievable?
print(f"\nSpecific checks:")
print(f"  H=7 achievable? {7 in all_H_values}")
print(f"  H=21 achievable? {21 in all_H_values}")
print(f"  H=9 achievable? {9 in all_H_values}")
print(f"  H=11 achievable? {11 in all_H_values}")

# Maximum H at each n
print(f"\nMax H by n:")
for n in sorted(max_H_by_n):
    print(f"  n={n}: max H = {max_H_by_n[n]}")

print("\nDone.")
