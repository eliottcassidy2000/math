#!/usr/bin/env python3
"""
THM-340 Verification: Q(d,k) = [x^d] B(x)^k
where B(x) = sum_{a>=2} SC(a+1) * x^a.

The exact d-good-cut formula is:
  exactly-d-good(n) = sum_{k=1}^{floor(d/2)} C(n-d, k) * Q(d,k)
where Q(d,k) counts labeled arrangements of k non-overlapping
strongly-connected tiling blocks of total width d.

This session proves: Q(d,k) = [x^d] B(x)^k.
"""

from math import comb
from fractions import Fraction

# SC tiling counts: SC[n] = #{path-fixed SC tilings on n vertices}
# n=2: 1 (base), n=3: 1, n=4: 5, n=5: 50, n=6: 903, n=7: 30773, ...
SC = {2: 1, 3: 1, 4: 5, 5: 50, 6: 903, 7: 30773, 8: 2032504,
      9: 264271477, 10: 68184627441, 11: 35047197032002}

def B_power_coeff(d, k, max_a=12):
    """Compute [x^d] B(x)^k where B(x) = sum_{a>=2} SC(a+1)*x^a."""
    # Iterate over all compositions of d into k parts, each >= 2
    result = 0
    def recurse(remaining, parts_left, current_product):
        nonlocal result
        if parts_left == 0:
            if remaining == 0:
                result += current_product
            return
        for a in range(2, remaining - 2*(parts_left-1) + 1):
            s = a + 1
            if s in SC:
                recurse(remaining - a, parts_left - 1, current_product * SC[s])
    recurse(d, k, 1)
    return result

# The Q(d,k) values from the previous session (verified)
Q_known = {
    (2,1): 1, (3,1): 5, (4,1): 50, (5,1): 903, (6,1): 30773,
    (7,1): 2032504, (8,1): 264271477, (9,1): 68184627441,
    (4,2): 1, (5,2): 10, (6,2): 125, (7,2): 2306, (8,2): 73076,
    (9,2): 4463038, (10,2): 552760703,
    (6,3): 1, (7,3): 15, (8,3): 225, (9,3): 4334, (10,3): 130659,
    (8,4): 1, (9,4): 20, (10,4): 350,
    (10,5): 1
}

print("=" * 65)
print("THM-340 VERIFICATION: Q(d,k) = [x^d] B(x)^k")
print("=" * 65)

all_match = True
for (d, k), expected in sorted(Q_known.items()):
    computed = B_power_coeff(d, k)
    match = (computed == expected)
    if not match:
        all_match = False
    status = "✓" if match else "✗ MISMATCH"
    print(f"  Q({d},{k}): expected={expected:>15}, computed={computed:>15}  {status}")

print()
if all_match:
    print("ALL Q(d,k) VALUES VERIFIED: Q(d,k) = [x^d]B(x)^k  ✓")
else:
    print("MISMATCH DETECTED!")

# Extend Q triangle to more terms
print()
print("=" * 65)
print("EXTENDED Q(d,k) TRIANGLE (d=2..15, k=1..floor(d/2))")
print("=" * 65)
print()
print("B(x) = SC(3)x^2 + SC(4)x^3 + SC(5)x^4 + ... = x^2 + 5x^3 + 50x^4 + ...")
print()

for d in range(2, 16):
    row = []
    for k in range(1, d//2 + 1):
        q = B_power_coeff(d, k)
        row.append(q)
    row_str = ", ".join(str(x) for x in row)
    print(f"  d={d:2d}: [{row_str}]")

# Verify the exactly-d-good formula using Q(d,k)
print()
print("=" * 65)
print("EXACTLY-d-GOOD FORMULA VERIFICATION (using Q(d,k) = [x^d]B^k)")
print("=" * 65)

# SC tiling counts via IE
def sc_count(n):
    """Count path-fixed SC tilings on n vertices via IE formula."""
    m = (n-1)*(n-2)//2
    if n <= 1:
        return 1 if n >= 0 else 0
    # f(S) = sum over nonempty T⊆S of (-1)^{|T|+1} h(T)
    # h({k}) = k(n-k)-1 (# tiles that cross cut k)
    # h(T) = min(T)*(n-max(T)) for |T|>=2
    total = 0
    n_cuts = n - 1
    # Inclusion-exclusion: SC = sum over all subsets S of cuts (sign * 2^{f(S)})
    # where f(S) = m - |{tiles crossing some cut in S}|
    from itertools import combinations
    def cuts_crossed_by_set(S):
        if not S:
            return 0
        from itertools import combinations as comb2
        # |tiles crossing ANY cut in S| by IE
        total_crossed = 0
        cuts = sorted(S)
        for r in range(1, len(cuts)+1):
            sign = (-1)**(r+1)
            for subset in comb2(cuts, r):
                lo = min(subset)
                hi = max(subset)
                if r == 1:
                    val = subset[0] * (n - subset[0]) - 1
                else:
                    val = lo * (n - hi)
                total_crossed += sign * val
        return total_crossed
    
    total = 2**m  # all tilings
    # Subtract non-SC via cuts
    for s_size in range(1, n_cuts+1):
        sign = (-1)**s_size
        for S in combinations(range(1, n_cuts+1), s_size):
            f_val = cuts_crossed_by_set(S)
            total += sign * (2**f_val)
    return total

# Get exact d-good counts from the tiling enumeration
def exact_d_good_formula(n, d, Q_cache={}):
    """Compute exactly-d-good(n) using Q(d,k)."""
    result = 0
    for k in range(1, d//2 + 1):
        key = (d, k)
        if key not in Q_cache:
            Q_cache[key] = B_power_coeff(d, k)
        result += comb(n-d, k) * Q_cache[key]
    return result

def exact_d_good_brute(n, d):
    """Brute-force count tilings with exactly d good cuts."""
    m = (n-1)*(n-2)//2
    count = 0
    for mask in range(2**m):
        # Count good cuts: cut k is good if any tile (x,y) with x>=k>y is upward
        # Tile ordering: for y=1..n-2, x=n down to y+2
        tiles = []
        for y in range(1, n-1):
            for x in range(n, y+1, -1):
                tiles.append((x, y))
        # Assign upward/downward based on mask
        tile_up = [(mask >> i) & 1 for i in range(len(tiles))]
        # Count good cuts
        good_cuts = 0
        for k in range(1, n):
            for i, (x, y) in enumerate(tiles):
                if x >= k > y and tile_up[i]:
                    good_cuts += 1
                    break
        if good_cuts == d:
            count += 1
    return count

print()
print("Spot-check formula vs brute force:")
errors = 0
for n in range(4, 10):
    for d in range(0, n):
        formula_val = exact_d_good_formula(n, d)
        brute_val = exact_d_good_brute(n, d)
        if formula_val != brute_val:
            print(f"  MISMATCH n={n}, d={d}: formula={formula_val}, brute={brute_val}")
            errors += 1

if errors == 0:
    print("  All formula values match brute force for n=4..9, d=0..n-1 ✓")

# Show the exact-d-good table
print()
print("exactly-d-good(n,d) table (formula):")
print(f"{'n':>4}", end="")
for d in range(0, 12):
    print(f"  d={d:2d}", end="")
print()
for n in range(3, 14):
    print(f"  n={n:2d}", end="")
    for d in range(0, min(n, 12)):
        val = exact_d_good_formula(n, d)
        print(f"  {val:5}", end="")
    print()

