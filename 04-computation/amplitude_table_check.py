#!/usr/bin/env python3
"""
AMPLITUDE TABLE CONSISTENCY CHECK
===================================
THM-080 states: |hat{M[a,b]}[S]| = 2^s * (n-2-|S|)! / 2^{n-2}

The theorem file contains an amplitude table (lines 156-161).
Let's check if the table entries are consistent with the stated formula.

opus-2026-03-16-S73
"""
from math import factorial

print("=" * 70)
print("AMPLITUDE TABLE vs FORMULA CHECK")
print("=" * 70)
print()
print("Formula: amplitude = 2^s * (n-2-d)! / 2^{n-2}")
print("where d = |S| = Walsh degree, s = # unrooted even-length components")
print()

# Table from THM-080 (lines 156-161)
table = {
    (3, 1, 0): "1/2",
    (5, 1, 0): "1/4",
    (5, 3, 0): "1/8",
    (7, 1, 0): "3/4",
    (7, 3, 0): "1/16",
    (7, 3, 1): "1/8",
    (7, 5, 0): "1/32",
    (9, 1, 0): "3/2",
    (9, 3, 0): "3/8",
    (9, 3, 1): "3/4",
    (9, 5, 0): "1/16",
    (9, 7, 0): "1/128",
}

def frac_str(num, den):
    """Convert to simplified fraction string."""
    from math import gcd
    g = gcd(abs(num), abs(den))
    return f"{num//g}/{den//g}"

print(f"{'(n,d,s)':<12} {'Table':>10} {'Formula':>12} {'Match':>7}")
print("-" * 45)

mismatches = []
for (n, d, s), table_val in sorted(table.items()):
    # Compute from formula
    num = (2**s) * factorial(n - 2 - d)
    den = 2**(n - 2)
    formula_val = frac_str(num, den)

    match = table_val == formula_val
    flag = "  ✓" if match else "  ✗ !!!"
    print(f"({n},{d},{s}){'':<5} {table_val:>10} {formula_val:>12} {flag}")

    if not match:
        mismatches.append((n, d, s, table_val, formula_val, num/den))

if mismatches:
    print(f"\n{'='*70}")
    print(f"FOUND {len(mismatches)} MISMATCHES!")
    print(f"{'='*70}")
    for n, d, s, tv, fv, fval in mismatches:
        # Parse table value
        parts = tv.split('/')
        tval = int(parts[0]) / int(parts[1])
        ratio = tval / fval
        print(f"  n={n}, d={d}, s={s}: table={tv}={tval}, formula={fv}={fval}, ratio={ratio:.6f}")

    print(f"\nAll mismatches are at n=9. The formula works perfectly for n=3,5,7.")
    print(f"The n=9 row appears to contain incorrect values.")

    # Show what the CORRECT table should be
    print(f"\nCORRECT n=9 amplitudes from the formula:")
    for d in [1, 3, 5, 7]:
        for s in [0, 1, 2]:
            num = (2**s) * factorial(9 - 2 - d)
            den = 2**7
            val = num / den
            if val < 100 and (d, s) != (1, 1):  # reasonable entries
                print(f"  d={d}, s={s}: {frac_str(num, den)} = {val}")
                if val < 0.001:
                    break

print()
print("=" * 70)
print("LADDER RATIO CHECK: H_amp(d,r) / M_amp(d-1,s=r-1)")
print("Predicted: = n - d always")
print("=" * 70)

# H formula: 2^r * (n-|S|)! / 2^{n-1}  (even Walsh degrees)
# M formula: 2^s * (n-2-|S|)! / 2^{n-2}  (odd Walsh degrees)
# Ladder: H at degree d with r components / M at degree d-1 with s=r-1 components

print(f"\n{'n':>3} {'H_deg':>6} {'r':>3} {'H_amp':>12} {'M_deg':>6} {'s':>3} {'M_amp':>12} {'Ratio':>8} {'n-d':>5}")
print("-" * 70)

for n in [5, 7, 9]:
    for d_H in range(2, n, 2):  # H even degrees
        for r in range(1, d_H // 2 + 1):
            H_num = (2**r) * factorial(n - d_H)
            H_den = 2**(n - 1)
            H_amp = H_num / H_den

            d_M = d_H - 1  # M odd degree, one less
            s = r - 1       # one fewer unrooted component
            if n - 2 - d_M < 0:
                continue
            M_num = (2**s) * factorial(n - 2 - d_M)
            M_den = 2**(n - 2)
            M_amp = M_num / M_den

            if M_amp > 0:
                ratio = H_amp / M_amp
                expected = n - d_H
                match = abs(ratio - expected) < 1e-10
                flag = "✓" if match else "✗"
                print(f"{n:>3} {d_H:>6} {r:>3} {H_amp:>12.6f} {d_M:>6} {s:>3} {M_amp:>12.6f} {ratio:>8.2f} {expected:>5} {flag}")

print()
print("=" * 70)
print("BIT ACCUMULATION: 2-adic structure of the ladder")
print("=" * 70)

print(f"\nThe ladder ratio H_amp / M_amp = n - d. Let's look at its 2-adic valuation:")
for n in [5, 7, 9, 11]:
    print(f"\n  n={n}:")
    for d_H in range(2, min(n, 10), 2):
        ratio = n - d_H
        v2 = 0
        r = ratio
        while r % 2 == 0 and r > 0:
            v2 += 1
            r //= 2
        print(f"    d={d_H}: ratio = {ratio} = 2^{v2} * {ratio // 2**v2}")

print("\nDone.")
