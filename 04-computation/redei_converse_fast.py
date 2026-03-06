#!/usr/bin/env python3
"""
Converse of Redei's theorem: which odd integers arise as H(T)?
Fast version - up to n=6 exhaustive, n=7 uses DP.

kind-pasteur-2026-03-06-S21
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count

print("=" * 60)
print("CONVERSE OF REDEI: ACHIEVABLE H VALUES")
print("=" * 60)

all_H_values = set()

for n in range(2, 7):  # n=2 through n=6
    m = n * (n-1) // 2
    h_values = set()
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        h = hamiltonian_path_count(T)
        h_values.add(h)
    all_H_values.update(h_values)
    h_sorted = sorted(h_values)
    print(f"\nn={n}: {len(h_values)} distinct H values, range [{min(h_values)}, {max(h_values)}]")
    print(f"  Values: {h_sorted}")

# For n=7, sample rather than exhaustive
print(f"\nn=7: sampling regular + random tournaments...")
import random
random.seed(42)
n = 7
m = n*(n-1)//2
h_values_7 = set()

# All regular tournaments (score sequence 3,3,3,3,3,3,3)
# We know from previous work there are 2640 of them
# But scanning all 2^21 takes too long
# Instead, sample 50000 random tournaments
for _ in range(50000):
    bits = random.randint(0, (1 << m) - 1)
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    h_values_7.add(h)

all_H_values.update(h_values_7)
print(f"  {len(h_values_7)} distinct H values from 50000 samples")
print(f"  Range: [{min(h_values_7)}, {max(h_values_7)}]")
print(f"  Values: {sorted(h_values_7)}")

# Overall gaps
max_val = max(all_H_values)
achieved_odd = {h for h in all_H_values if h % 2 == 1}
all_odd_up_to = set(range(1, max_val + 1, 2))
gaps = sorted(all_odd_up_to - achieved_odd)

print(f"\n{'='*60}")
print(f"Odd integers 1..{max_val} NOT achieved (up to n=6 exhaustive + n=7 sample):")
print(f"Total gaps: {len(gaps)}")
if len(gaps) <= 80:
    print(f"Gaps: {gaps}")
else:
    print(f"First 40: {gaps[:40]}")

# Specific checks
for k in [7, 21, 23, 29, 31, 37, 39, 41, 43, 47]:
    status = "ACHIEVED" if k in all_H_values else "GAP"
    print(f"  H={k}: {status}")

print("\nDone.")
