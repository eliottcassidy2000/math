"""
vitali_n8_decomposition.py -- kind-pasteur-2026-03-13-S61
Verify: delta_H = 2*delta_c7 + 4*delta_{disjoint (3,5)-pairs} at n=8.

At n=8, the OCF gives H = 1 + 2*alpha_1 + 4*alpha_2 where:
- alpha_1 = total odd directed cycles (c3 + c5 + c7)
- alpha_2 = disjoint pairs of odd cycles in Omega(T)

Since c3 and c5 are preserved, delta_alpha_1 = delta_c7.
alpha_2 has contributions from:
- Disjoint (3,3)-pairs: need 6 vertices. PRESERVED (proved at n=7, likely at n=8 too).
- Disjoint (3,5)-pairs: need exactly 8 vertices. CAN CHANGE (new at n=8).
- Disjoint (5,5)-pairs: need 10 > 8. IMPOSSIBLE.
- Anything involving 7-cycles: need >= 10. IMPOSSIBLE.

So delta_alpha_2 = delta_{disjoint (3,5)-pairs} + delta_{disjoint (3,3)-pairs}
=> delta_H = 2*delta_c7 + 4*(delta_33 + delta_35)

Test whether delta_33 = 0 and the full formula holds.
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def find_3_cycles(A, n):
    """Find all directed 3-cycle vertex sets."""
    cycles = []
    for combo in combinations(range(n), 3):
        i, j, k = combo
        if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]):
            cycles.append(frozenset(combo))
    return cycles

def find_5_cycles(A, n):
    """Find all directed 5-cycle vertex sets."""
    cycles = set()
    for combo in combinations(range(n), 5):
        has_cycle = False
        for perm in permutations(combo[1:]):
            path = (combo[0],) + perm
            valid = True
            for k in range(5):
                if A[path[k]][path[(k+1) % 5]] != 1:
                    valid = False
                    break
            if valid:
                has_cycle = True
                break
        if has_cycle:
            cycles.add(frozenset(combo))
    return list(cycles)

def count_7_cycles(A, n):
    """Count directed 7-cycles on each 7-vertex subset."""
    total = 0
    for combo in combinations(range(n), 7):
        for perm in permutations(combo[1:]):
            path = (combo[0],) + perm
            valid = True
            for k in range(7):
                if A[path[k]][path[(k+1) % 7]] != 1:
                    valid = False
                    break
            if valid:
                total += 1
    return total

def count_disjoint_pairs(cycles_a, cycles_b):
    """Count disjoint pairs between two cycle lists."""
    count = 0
    for c1 in cycles_a:
        for c2 in cycles_b:
            if len(c1 & c2) == 0:
                count += 1
    return count

def count_disjoint_same(cycles):
    """Count disjoint pairs within one cycle list."""
    count = 0
    for i in range(len(cycles)):
        for j in range(i+1, len(cycles)):
            if len(cycles[i] & cycles[j]) == 0:
                count += 1
    return count

n = 8
total_bits = n * (n - 1) // 2
np.random.seed(2026)

print("=" * 70)
print("VERIFYING: delta_H = 2*delta_c7 + 4*(delta_33 + delta_35) at n=8")
print("=" * 70)

verified = 0
failed = 0
all_data = []

for trial in range(1500):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    lam_orig = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue

        B = reverse_subtournament(A, n, list(subset))
        lam_new = lambda_graph(B, n)

        if not np.array_equal(lam_orig, lam_new):
            continue

        H_A = count_ham_paths(A, n)
        H_B = count_ham_paths(B, n)
        delta_H = H_B - H_A

        # 3-cycles
        c3_A = find_3_cycles(A, n)
        c3_B = find_3_cycles(B, n)

        # 5-cycles
        c5_A = find_5_cycles(A, n)
        c5_B = find_5_cycles(B, n)

        # 7-cycle count
        c7_A = count_7_cycles(A, n)
        c7_B = count_7_cycles(B, n)
        delta_c7 = c7_B - c7_A

        # Disjoint (3,3)-pairs
        d33_A = count_disjoint_same(c3_A)
        d33_B = count_disjoint_same(c3_B)
        delta_33 = d33_B - d33_A

        # Disjoint (3,5)-pairs
        d35_A = count_disjoint_pairs(c3_A, c5_A)
        d35_B = count_disjoint_pairs(c3_B, c5_B)
        delta_35 = d35_B - d35_A

        # Disjoint (5,5)-pairs (should be 0 at n=8, need 10 vertices)
        d55_A = count_disjoint_same(c5_A)
        d55_B = count_disjoint_same(c5_B)
        delta_55 = d55_B - d55_A

        # Formula: delta_H = 2*delta_c7 + 4*(delta_33 + delta_35)
        predicted = 2 * delta_c7 + 4 * (delta_33 + delta_35)

        if predicted == delta_H:
            verified += 1
        else:
            failed += 1
            print(f"  FAIL: dH={delta_H}, dc7={delta_c7}, d33={delta_33}, d35={delta_35}, d55={delta_55}, predicted={predicted}")

        all_data.append({
            'delta_H': delta_H,
            'delta_c7': delta_c7,
            'delta_33': delta_33,
            'delta_35': delta_35,
            'delta_55': delta_55,
        })

    if trial % 300 == 0 and trial > 0:
        print(f"  ... {trial}/1500, verified={verified}, failed={failed}")

print(f"\nRESULT: verified={verified}, failed={failed}")
print(f"Formula delta_H = 2*delta_c7 + 4*(delta_33 + delta_35): {'CONFIRMED' if failed == 0 else 'FAILS'}")

# Distribution analysis
print(f"\n  delta_33 distribution: {Counter(d['delta_33'] for d in all_data)}")
print(f"  delta_35 distribution: {Counter(d['delta_35'] for d in all_data)}")
print(f"  delta_55 distribution: {Counter(d['delta_55'] for d in all_data)}")

# Show cases where delta_33 != 0
d33_nonzero = [d for d in all_data if d['delta_33'] != 0]
if d33_nonzero:
    print(f"\n  delta_33 != 0 cases: {len(d33_nonzero)}")
    for d in d33_nonzero[:5]:
        print(f"    dH={d['delta_H']}, dc7={d['delta_c7']}, d33={d['delta_33']}, d35={d['delta_35']}")
else:
    print(f"\n  delta_33 = 0 ALWAYS (disjoint 3-3 pairs preserved) - CONFIRMED at n=8")

# Show the full decomposition
print(f"\n  Full decomposition table:")
print(f"  {'dH':>4} {'dc7':>4} {'d33':>4} {'d35':>4} {'d55':>4} {'2dc7+4(d33+d35)':>16} {'match':>6}")
for d in all_data:
    pred = 2 * d['delta_c7'] + 4 * (d['delta_33'] + d['delta_35'])
    match = 'OK' if pred == d['delta_H'] else 'FAIL'
    if d['delta_H'] != 0:  # only show changing cases
        print(f"  {d['delta_H']:>4} {d['delta_c7']:>4} {d['delta_33']:>4} {d['delta_35']:>4} {d['delta_55']:>4} {pred:>16} {match:>6}")

print("\n" + "=" * 70)
print("THE HIERARCHICAL STRUCTURE")
print("=" * 70)
print("""
At each n, the Vitali atom activates different "channels":

n=5: No change (trivially gauge-invariant)
n=6: No change (c5 invariant, no 7-cycles)
n=7: delta_H = 2*delta_c7 (only c7 channel, delta=+-1 so dH=+-2)
n=8: delta_H = 2*delta_c7 + 4*delta_35 (two channels: c7 and disjoint (3,5))
     delta_c7 in {-1,0,1,2}, delta_35 in {-1,0,1}
     => delta_H in {-6,-4,-2,0,2,4,6}

Prediction for n=9:
  New channels: delta_c9 (9-cycles on all 9 vertices)
  alpha_3 = disjoint triples: three 3-cycles need 9 = n vertices!
  So delta_H = 2*delta_c7 + 2*delta_c9 + 4*delta_35 + 4*delta_37 + 8*delta_333
  where delta_333 = change in disjoint 3-3-3 triples (requires exactly n=9)

This is the {2,1,0} overlap weight hierarchy in action:
- Each new n opens a new "channel" in the OCF
- The Vitali atom pushes through ALL open channels simultaneously
- The total delta_H is a WEIGHTED SUM over channels
- Weights follow the OCF powers of 2: alpha_k gets weight 2^k
""")

print("\nDone.")
