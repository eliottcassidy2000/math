#!/usr/bin/env python3
"""
DOES THE BIBD ARRANGEMENT MAXIMIZE DISJOINT PAIRS?

For a collection of b triples from [p] (with each vertex in r triples and each
pair in lambda triples), the number of disjoint pairs is:

  D = C(b,2) - sum_v C(d(v),2) + sum_e C(lambda_e, 2)

For BIBD: D_bibd = C(b,2) - p*C(r,2) + C(p,2)*C(lambda,2)

Question: is D_bibd the maximum over ALL arrangements with b triples?

The two Jensen terms work in opposite directions:
  - min sum_v C(d(v),2) at d(v) = r (INCREASES D)
  - min sum_e C(lambda_e,2) at lambda_e = lambda (DECREASES D)

So it's not obvious. Let's test computationally at small sizes.

kind-pasteur-2026-03-06-S18h
"""
from math import comb
from itertools import combinations
from collections import defaultdict

def count_disjoint_pairs(triples, n):
    """Count pairs of triples that share no vertices."""
    b = len(triples)
    count = 0
    for i in range(b):
        si = set(triples[i])
        for j in range(i+1, b):
            if not (si & set(triples[j])):
                count += 1
    return count

def stats(triples, n):
    """Compute d(v), lambda_e stats and verify formula."""
    b = len(triples)
    d = defaultdict(int)
    lam = defaultdict(int)
    for t in triples:
        for v in t:
            d[v] += 1
        for i in range(3):
            for j in range(i+1, 3):
                pair = (min(t[i], t[j]), max(t[i], t[j]))
                lam[pair] += 1

    sum_Cd = sum(comb(d[v], 2) for v in range(n))
    sum_Clam = sum(comb(lam[e], 2) for e in lam)

    D_formula = comb(b, 2) - sum_Cd + sum_Clam
    D_actual = count_disjoint_pairs(triples, n)

    return {
        'd_vals': sorted(d[v] for v in range(n)),
        'lam_vals': sorted(lam.values()),
        'D': D_actual,
        'D_formula': D_formula,
        'sum_Cd': sum_Cd,
        'sum_Clam': sum_Clam,
    }

# ============================================================
# TEST 1: n=7, all possible sets of 14 triples
# ============================================================
# This is way too large (C(35, 14) ~ 3 billion).
# Instead, test with actual tournament 3-cycles.

# Generate all n=7 tournaments and their 3-cycles
print("=" * 70)
print("TEST: BIBD vs NON-BIBD DISJOINT PAIRS")
print("=" * 70)

# First, let's verify the formula
print("\nVerifying D = C(b,2) - sum C(d(v),2) + sum C(lam_e, 2):")

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits

def find_3cycles_T(T):
    n = len(T)
    cycles = []
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    cycles.append((i,j,k))
                elif T[i][k] and T[k][j] and T[j][i]:
                    cycles.append((i,j,k))
    return cycles

# Small test at n=5
n = 5
m = n*(n-1)//2
print(f"\nn={n}: formula verification (all {1 << m} tournaments)")
errors = 0
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    c3 = find_3cycles_T(T)
    if len(c3) < 2:
        continue
    s = stats(c3, n)
    if s['D'] != s['D_formula']:
        errors += 1
        print(f"  ERROR at bits={bits}: D={s['D']}, formula={s['D_formula']}")
print(f"  Errors: {errors}")

# ============================================================
# TEST 2: At n=6, among regular tournaments (score (2,2,2,3,3,3)
# or (3,3,3,3,3) at n=5), compare BIBD-like vs others
# ============================================================
print(f"\n{'=' * 70}")
print("n=6: DISJOINT PAIRS BY LAMBDA UNIFORMITY")
print("=" * 70)

n = 6
m = n*(n-1)//2
by_c3 = defaultdict(list)  # c3 -> list of (D, d_vals, lam_vals)

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    c3 = find_3cycles_T(T)
    if len(c3) < 2:
        continue
    s = stats(c3, n)
    by_c3[len(c3)].append(s)

for c3 in sorted(by_c3.keys()):
    data = by_c3[c3]
    D_vals = [s['D'] for s in data]
    max_D = max(D_vals)
    # Check if max D is achieved by most uniform lambda
    max_D_entries = [s for s in data if s['D'] == max_D]
    # Compute lambda variance for each
    for entry in max_D_entries[:2]:
        lam_var = 0
        lam_vals = entry['lam_vals']
        if lam_vals:
            mean_lam = sum(lam_vals) / len(lam_vals)
            lam_var = sum((l - mean_lam)**2 for l in lam_vals) / len(lam_vals)
        print(f"  c3={c3}: max D = {max_D}, D range [{min(D_vals)},{max_D}], "
              f"max_D lambda_vals={entry['lam_vals'][:10]}, lambda_var={lam_var:.2f}")
        break

# ============================================================
# TEST 3: THE KEY QUESTION AT n=7
# Among REGULAR tournaments (c3=14 forced), does the Paley (BIBD)
# arrangement maximize D?
# ============================================================
print(f"\n{'=' * 70}")
print("n=7 REGULAR: BIBD vs NON-BIBD ALPHA_2")
print("=" * 70)

n = 7
m = n*(n-1)//2

# Find regular tournaments (score all 3s)
regular_data = []
count = 0
for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    scores = sorted(sum(T[i]) for i in range(n))
    if scores != [3]*7:
        continue
    count += 1
    c3 = find_3cycles_T(T)
    s = stats(c3, n)
    regular_data.append(s)
    if count % 50 == 0:
        print(f"  Found {count} regular tournaments...", flush=True)

print(f"\nTotal regular n=7 tournaments: {len(regular_data)}")
D_vals = [s['D'] for s in regular_data]
print(f"Alpha_2 (disjoint pairs) range: [{min(D_vals)}, {max(D_vals)}]")

# Distribution
D_dist = defaultdict(int)
for d in D_vals:
    D_dist[d] += 1
for d in sorted(D_dist.keys(), reverse=True):
    lam_examples = [s for s in regular_data if s['D'] == d]
    # Check lambda variance
    variances = []
    for s in lam_examples[:5]:
        lam_vals = s['lam_vals']
        if lam_vals:
            mean_l = sum(lam_vals) / len(lam_vals)
            var_l = sum((l - mean_l)**2 for l in lam_vals) / len(lam_vals)
            variances.append(var_l)
    avg_var = sum(variances)/len(variances) if variances else 0
    print(f"  D={d}: {D_dist[d]} tournaments, avg lambda_variance={avg_var:.3f}")

# Check if ALL D=7 tournaments have uniform lambda (BIBD)
max_D = max(D_vals)
bibd_entries = [s for s in regular_data if s['D'] == max_D]
print(f"\nTournaments with max D={max_D}:")
for i, s in enumerate(bibd_entries[:5]):
    print(f"  #{i+1}: d_vals={s['d_vals']}, lambda_vals={s['lam_vals']}")

# Check: are ALL max-D regular tournaments BIBDs?
all_bibd = True
for s in bibd_entries:
    if len(set(s['lam_vals'])) != 1:
        all_bibd = False
        break
print(f"All max-D regular tournaments are BIBDs? {all_bibd}")

# Check: Paley T_7 specifically
def build_paley(p):
    qr = set()
    for x in range(1, p):
        qr.add((x * x) % p)
    T = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and (j - i) % p in qr:
                T[i][j] = 1
    return T

T7 = build_paley(7)
c3_paley = find_3cycles_T(T7)
s_paley = stats(c3_paley, 7)
print(f"\nPaley T_7: D={s_paley['D']}, d_vals={s_paley['d_vals']}, lambda_vals={s_paley['lam_vals']}")

print("\nDone.")
