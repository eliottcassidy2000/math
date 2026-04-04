"""
a000568_extend_and_analyze.py — opus-2026-04-04-S14

Use the existing C enumerator to push A000568 further and analyze patterns.

GOALS:
1. Compute a(n) mod small primes for n up to 500+ (fast modular arithmetic)
2. Look for patterns: a(n) mod 2, mod 3, mod 7, mod 42
3. Growth rate analysis: log(a(n)) / C(n,2) → log(2)
4. Connection to recent session insights
5. The fiber bundle predicts: a(n) ≈ Σ orbits(classes at n-1) / n
"""
import subprocess
import time
from math import comb, factorial, log, log2
from collections import Counter

def compute_a000568(n, threads=8):
    """Call the C enumerator and return the result."""
    result = subprocess.run(
        ['04-computation/burnside_enum_v2', '568', str(n), str(threads)],
        capture_output=True, text=True, timeout=300
    )
    # Parse the output to get a(n)
    for line in result.stdout.split('\n'):
        if line.startswith(f'a({n})'):
            parts = line.split('=')
            if len(parts) >= 2:
                return int(parts[1].strip().split('\n')[0])
    return None

# ==============================================================
# 1. Compute a(n) for a range and analyze mod structure
# ==============================================================
print("A000568 EXTENSION AND ANALYSIS")
print("=" * 60)

# Compute for a range of n
values = {}
for n in range(3, 101):
    t0 = time.time()
    val = compute_a000568(n)
    elapsed = time.time() - t0
    if val is not None:
        values[n] = val
        if n <= 20 or n % 10 == 0:
            print(f"  a({n:3d}) = {val} ({elapsed:.3f}s)")

# ==============================================================
# 2. MOD STRUCTURE: a(n) mod 2, 3, 7, 42
# ==============================================================
print(f"\n{'='*60}")
print(f"MOD STRUCTURE")
print(f"{'='*60}")

for p in [2, 3, 5, 7, 11, 42]:
    residues = [values[n] % p for n in sorted(values)]
    print(f"\n  a(n) mod {p} (n=3..100):")
    # Show first 30
    print(f"    {residues[:30]}")
    # Distribution
    dist = Counter(residues)
    print(f"    Distribution: {dict(sorted(dist.items()))}")

    # Period detection: is there a period?
    for period in range(1, 25):
        if all(residues[i] == residues[i + period]
               for i in range(min(30, len(residues) - period))):
            print(f"    *** PERIODIC with period {period} ***")
            print(f"    Pattern: {residues[:period]}")
            break

# ==============================================================
# 3. GROWTH RATE: log2(a(n)) vs C(n,2) * log2(2) - log2(n!)
# ==============================================================
print(f"\n{'='*60}")
print(f"GROWTH RATE ANALYSIS")
print(f"{'='*60}")

for n in [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]:
    if n in values and values[n] > 0:
        a = values[n]
        cn2 = comb(n, 2)
        nfact = factorial(n)
        ratio = log2(a) / cn2
        # a(n) ~ 2^{C(n,2)} / n!
        predicted_log = cn2 * log(2) - sum(log(k) for k in range(1, n+1))
        actual_log = log(a)
        correction = actual_log - predicted_log
        print(f"  n={n:3d}: log2(a) = {log2(a):.2f}, C(n,2) = {cn2}, "
              f"log2(a)/C(n,2) = {ratio:.6f}, "
              f"correction = {correction:.4f}")

# ==============================================================
# 4. CONSECUTIVE RATIOS: a(n)/a(n-1)
# ==============================================================
print(f"\n{'='*60}")
print(f"CONSECUTIVE RATIOS")
print(f"{'='*60}")

for n in range(4, 101):
    if n in values and n-1 in values and values[n-1] > 0:
        ratio = values[n] / values[n-1]
        # Expected: roughly 2^{(n-1)/2}
        expected = 2**((n-1)/2)
        if n <= 20 or n % 10 == 0:
            print(f"  a({n})/a({n-1}) = {ratio:.4f}, "
                  f"2^((n-1)/2) = {expected:.1f}, "
                  f"ratio/expected = {ratio/expected:.4f}")

# ==============================================================
# 5. THE FIBER BUNDLE PREDICTION
# ==============================================================
print(f"\n{'='*60}")
print(f"FIBER BUNDLE PREDICTION")
print(f"{'='*60}")

# From S12: merging_ratio ≈ n (approaches n for large n)
# So: a(n) ≈ Σ orbits / n where Σ orbits = a(n-1) * 2^{n-1} / avg_aut
# If avg_aut ≈ 1 (most tournaments have trivial Aut):
# Σ orbits ≈ a(n-1) * 2^{n-1}
# a(n) ≈ a(n-1) * 2^{n-1} / n

# Check this approximation
print(f"  Testing a(n) ≈ a(n-1) * 2^(n-1) / n:")
for n in range(4, 30):
    if n in values and n-1 in values:
        predicted = values[n-1] * 2**(n-1) / n
        actual = values[n]
        ratio = actual / predicted
        if n <= 15 or n % 5 == 0:
            print(f"    n={n}: predicted={predicted:.1f}, actual={actual}, "
                  f"ratio={ratio:.4f}")

# Better: a(n) ≈ a(n-1) * 2^{n-1} / n * correction(n)
# where correction accounts for the non-trivial Aut tournaments
# At n=5: a(4)=4, 2^4/5*4 = 12.8 vs actual 12. Correction = 12/12.8 = 0.9375
# At n=6: a(5)=12, 2^5/6*12 = 64 vs actual 56. Correction = 56/64 = 0.875

# ==============================================================
# 6. a(n) mod 7 DEEP ANALYSIS (connection to forbidden values)
# ==============================================================
print(f"\n{'='*60}")
print(f"a(n) mod 7 — CONNECTION TO FORBIDDEN VALUES")
print(f"{'='*60}")

mod7_vals = [(n, values[n] % 7) for n in sorted(values)]
print(f"  a(n) mod 7 for n=3..50:")
for i in range(0, min(48, len(mod7_vals)), 8):
    chunk = mod7_vals[i:i+8]
    print(f"    {' '.join(f'{n}:{v}' for n, v in chunk)}")

# Is a(n) ever divisible by 7?
div7 = [n for n, v in mod7_vals if v == 0]
print(f"\n  a(n) divisible by 7 at: {div7[:20]}")

# a(n) mod 42
mod42_vals = [(n, values[n] % 42) for n in sorted(values)]
print(f"\n  a(n) mod 42 for n=3..20:")
for n, v in mod42_vals[:18]:
    print(f"    a({n}) mod 42 = {v}")

# ==============================================================
# 7. DIGIT GROWTH: How many digits does a(n) have?
# ==============================================================
print(f"\n{'='*60}")
print(f"DIGIT GROWTH")
print(f"{'='*60}")

for n in range(5, 101, 5):
    if n in values:
        digits = len(str(values[n]))
        cn2 = comb(n, 2)
        predicted_digits = int(cn2 * log(2) / log(10) - sum(log(k)/log(10) for k in range(1, n+1)))
        print(f"  a({n:3d}): {digits:5d} digits (predicted ~{predicted_digits})")

if __name__ == "__main__":
    pass  # Already runs at import
