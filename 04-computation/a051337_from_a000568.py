#!/usr/bin/env python3
"""
Compute A051337 (strongly connected tournaments) from A000568 via OGF inversion.

If B(x) = sum T(n) x^n is the OGF for non-iso tournaments,
then C(x) = 1 - 1/B(x) gives SC(n) = strongly connected non-iso tournaments.

Uses our computed A000568 values from the b-file (n=0 to 150+).

Author: opus-2026-03-08-S48
"""

import os
from fractions import Fraction

def load_bfile(path):
    """Load OEIS b-file format: lines of 'n value'."""
    vals = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                n, v = int(parts[0]), int(parts[1])
                vals[n] = v
    return vals

# Load A000568 values
bfile = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'b000568.txt')
T = load_bfile(bfile)
max_n = max(T.keys())

print(f"Loaded A000568 values for n=0..{max_n}")

# Compute 1/B(x) = sum b_n x^n where B(x) = sum T(n) x^n
# b[0] = 1/T[0] = 1
# For n >= 1: b[n] = -sum_{k=1}^n T[k] * b[n-k]
b = [Fraction(0)] * (max_n + 1)
b[0] = Fraction(1)
for n in range(1, max_n + 1):
    s = Fraction(0)
    for k in range(1, n + 1):
        if k in T:
            s += T[k] * b[n - k]
    b[n] = -s

# SC(n) = -b[n] for n >= 1 (from C(x) = 1 - 1/B(x))
print(f"\nA051337: Strongly connected non-isomorphic tournaments")
print(f"{'n':>4}  {'SC(n)':>20}  {'T(n)':>20}  {'SC/T %':>8}")
print("-" * 60)

sc_values = {}
for n in range(1, min(max_n + 1, 201)):
    sc_n = int(-b[n])
    t_n = T[n]
    pct = sc_n / t_n * 100 if t_n > 0 else 0
    sc_values[n] = sc_n
    if n <= 20 or n % 10 == 0:
        print(f"{n:4d}  {sc_n:>20d}  {t_n:>20d}  {pct:>7.3f}%")

# Save as b-file
outpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'b051337.txt')
with open(outpath, 'w') as f:
    for n in sorted(sc_values.keys()):
        f.write(f"{n} {sc_values[n]}\n")
print(f"\nSaved {len(sc_values)} values to {outpath}")

# Check against known OEIS values
known_sc = {1: 1, 2: 0, 3: 1, 4: 1, 5: 6, 6: 28, 7: 282,
            8: 5765, 9: 185086, 10: 9698086}
print("\nValidation against known OEIS values:")
for n, expected in sorted(known_sc.items()):
    got = sc_values.get(n, "?")
    match = "OK" if got == expected else f"FAIL (got {got})"
    print(f"  SC({n}) = {expected}  {match}")
