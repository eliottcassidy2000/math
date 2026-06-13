"""
vitali_n8_extension.py -- kind-pasteur-2026-03-13-S61
Extend Vitali atom analysis to n=8.

At n=8: new phenomena expected because:
1. 8-cycles don't exist in odd-cycle counting (8 is even)
2. But 5-cycles on 5 of 8 vertices + 4-subset reversal could now
   change c5 (unlike n=6 where c5 was invariant)
3. The 4-subset leaves 4 external vertices (not 3), so the completion
   space is richer
4. Disjoint 3-cycle pairs can involve both subset and ext vertices

Questions:
- Does the lambda-preserving (1,1,2,2) reversal still only change H by +-2?
- Or can larger deltas occur?
- Does the cyclic ext rule still hold?
- What about 4-vertex external tournaments?
- Does c5 now become non-invariant?
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

def count_directed_cycles(A, n, length):
    if length > n:
        return 0
    count = 0
    for combo in combinations(range(n), length):
        for perm in permutations(combo[1:]):
            path = (combo[0],) + perm
            valid = True
            for k in range(length):
                if A[path[k]][path[(k+1) % length]] != 1:
                    valid = False
                    break
            if valid:
                count += 1
    return count

n = 8
total_bits = n * (n - 1) // 2
np.random.seed(2026)

print("=" * 70)
print("PART 1: LAMBDA-PRESERVING (1,1,2,2) REVERSALS AT n=8")
print("=" * 70)

results = []
sample_count = 0

for trial in range(1000):
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

        H_orig = count_ham_paths(A, n)
        H_new = count_ham_paths(B, n)
        delta_H = H_new - H_orig

        ext = [v for v in range(n) if v not in subset]
        sub_list = list(subset)

        # mixed_c3
        mixed_c3 = 0
        for s1, s2 in combinations(sub_list, 2):
            for e in ext:
                if (A[s1][s2] and A[s2][e] and A[e][s1]) or (A[s2][s1] and A[s1][e] and A[e][s2]):
                    mixed_c3 += 1

        # ext tournament type (4 vertices now!)
        ext_scores = tuple(sorted([sum(A[ext[i]][ext[j]] for j in range(4) if j != i) for i in range(4)]))

        # c5 counts
        c5_orig = count_directed_cycles(A, n, 5)
        c5_new = count_directed_cycles(B, n, 5)
        delta_c5 = c5_new - c5_orig

        # c7 counts
        c7_orig = count_directed_cycles(A, n, 7)
        c7_new = count_directed_cycles(B, n, 7)
        delta_c7 = c7_new - c7_orig

        results.append({
            'delta_H': delta_H,
            'delta_c5': delta_c5,
            'delta_c7': delta_c7,
            'mixed_c3': mixed_c3,
            'ext_scores': ext_scores,
            'H_orig': H_orig,
        })

        sample_count += 1

    if trial % 100 == 0 and trial > 0:
        print(f"  ... {trial}/1000, {sample_count} valid reversals so far")

print(f"\nTotal valid lambda-preserving (1,1,2,2) reversals: {sample_count}")

changing = [r for r in results if r['delta_H'] != 0]
preserving = [r for r in results if r['delta_H'] == 0]

print(f"  Changing: {len(changing)} ({len(changing)/sample_count*100:.1f}%)")
print(f"  Preserving: {len(preserving)} ({len(preserving)/sample_count*100:.1f}%)")

# Delta_H distribution
print(f"\n  delta_H distribution:")
delta_dist = Counter(r['delta_H'] for r in results)
for d in sorted(delta_dist.keys()):
    print(f"    delta_H={d}: {delta_dist[d]}")

# Delta_c5 and delta_c7
print(f"\n  delta_c5 distribution (changing only):")
c5_dist = Counter(r['delta_c5'] for r in changing)
for d in sorted(c5_dist.keys()):
    print(f"    delta_c5={d}: {c5_dist[d]}")

print(f"\n  delta_c7 distribution (changing only):")
c7_dist = Counter(r['delta_c7'] for r in changing)
for d in sorted(c7_dist.keys()):
    print(f"    delta_c7={d}: {c7_dist[d]}")

# Check: is delta_H = 2*delta_c5 + 4*delta_c7 + higher?
# By OCF: H = 1 + 2*alpha_1 + 4*alpha_2 + ... where alpha_1 = c3+c5+c7, alpha_2 = disjoint pairs
# Actually H = I(Omega,2) and these alpha terms don't decompose that simply.
# But delta_H should relate to delta of odd cycle counts.
print(f"\n  Checking delta_H = 2*delta_c5 + 4*delta_c7 + ...:")
for r in changing[:10]:
    # delta_H = 2*(delta_c5 + delta_c7) + 4*(delta_alpha_2)?
    # This formula is approximate. Let's just check the data.
    predicted = 2 * r['delta_c5'] + 4 * r['delta_c7']
    print(f"    dH={r['delta_H']}, dc5={r['delta_c5']}, dc7={r['delta_c7']}, 2*dc5+4*dc7={predicted}")

# mixed_c3 analysis
print(f"\n  mixed_c3 distribution:")
mc_dist_ch = Counter(r['mixed_c3'] for r in changing)
mc_dist_pr = Counter(r['mixed_c3'] for r in preserving)
all_mc = sorted(set(list(mc_dist_ch.keys()) + list(mc_dist_pr.keys())))
for mc in all_mc:
    c = mc_dist_ch.get(mc, 0)
    p = mc_dist_pr.get(mc, 0)
    pct = c/(c+p)*100 if c+p > 0 else 0
    print(f"    mixed_c3={mc}: changing={c}, preserving={p}, rate={pct:.1f}%")

# External tournament type
print(f"\n  External tournament type (4-vertex scores):")
ext_dist_ch = Counter(r['ext_scores'] for r in changing)
ext_dist_pr = Counter(r['ext_scores'] for r in preserving)
all_ext = sorted(set(list(ext_dist_ch.keys()) + list(ext_dist_pr.keys())))
for ext in all_ext:
    c = ext_dist_ch.get(ext, 0)
    p = ext_dist_pr.get(ext, 0)
    pct = c/(c+p)*100 if c+p > 0 else 0
    print(f"    {ext}: changing={c}, preserving={p}, rate={pct:.1f}%")

print("\n" + "=" * 70)
print("PART 2: CROSS-TAB AT n=8 — MIXED_C3 vs EXT TYPE")
print("=" * 70)

# At n=8, external tournament on 4 vertices can be:
# (0,1,2,3) = transitive
# (1,1,2,2) = non-trans, 2 directed 3-cycles
# Other score sequences

cross = Counter()
for r in results:
    key = (r['mixed_c3'], r['ext_scores'], r['delta_H'] != 0)
    cross[key] += 1

for mc in sorted(set(r['mixed_c3'] for r in results)):
    for ext in sorted(set(r['ext_scores'] for r in results)):
        c = cross.get((mc, ext, True), 0)
        p = cross.get((mc, ext, False), 0)
        if c + p > 0:
            pct = c/(c+p)*100
            print(f"  mixed_c3={mc}, ext={ext}: changing={c}, preserving={p}, rate={pct:.1f}%")

print("\n" + "=" * 70)
print("PART 3: DOES c5 NOW CHANGE? (New at n=8)")
print("=" * 70)

c5_changed = [r for r in results if r['delta_c5'] != 0]
print(f"  Reversals with delta_c5 != 0: {len(c5_changed)} out of {sample_count}")
if c5_changed:
    print(f"  delta_c5 values: {Counter(r['delta_c5'] for r in c5_changed)}")
    # Is c5 change independent of c7 change?
    c5_only = [r for r in c5_changed if r['delta_c7'] == 0]
    c7_only = [r for r in results if r['delta_c5'] == 0 and r['delta_c7'] != 0]
    both = [r for r in results if r['delta_c5'] != 0 and r['delta_c7'] != 0]
    print(f"  c5 only: {len(c5_only)}, c7 only: {len(c7_only)}, both: {len(both)}")
else:
    print(f"  c5 is STILL invariant at n=8!")
    print(f"  This means the phase transition for c5 is at n>=9 (if it exists)")

print("\n" + "=" * 70)
print("PART 4: DELTA DECOMPOSITION — WHAT DRIVES delta_H AT n=8?")
print("=" * 70)

# At n=8, H = 1 + 2*Omega_1 + 4*Omega_2 + 8*Omega_3 (OCF with independence polynomial)
# where Omega_k = independence number of size k in conflict graph
# The conflict graph Omega has one vertex per odd directed cycle
# Two cycles conflict if they share >= 2 vertices (W>=2 overlap)
# So independent sets are collections of "barely overlapping" cycles

# Check: is delta_H always even?
print(f"  delta_H values: {Counter(r['delta_H'] for r in results)}")
all_even = all(r['delta_H'] % 2 == 0 for r in results)
print(f"  All delta_H even? {all_even}")

# Check: delta_H distribution for changing cases
if changing:
    print(f"\n  |delta_H| distribution:")
    abs_dist = Counter(abs(r['delta_H']) for r in changing)
    for d in sorted(abs_dist.keys()):
        print(f"    |delta_H|={d}: {abs_dist[d]}")

print("\n" + "=" * 70)
print("PART 5: THE DIMENSIONAL HIERARCHY AT n=8")
print("=" * 70)

# At n=7: 4 subset + 3 ext = 7 vertices -> 7-cycles are the key
# At n=8: 4 subset + 4 ext = 8 vertices -> 7-cycles and 5-cycles on various subsets
# But 8-cycles don't exist in OCF (even length)
# The key question: do 7-cycles on 7-of-8 vertices now contribute?
# There are C(8,7)=8 possible 7-vertex subsets, each containing
# at least 3 subset vertices (since |subset|=4 and |7-subset intersect subset| >= 4+7-8 = 3)

# Count c7 on each 7-vertex subset
if changing:
    print("\n  Detailed cycle analysis for H-changing examples at n=8:")
    for idx in range(min(5, len(changing))):
        r = changing[idx]
        # Reconstruct the tournament
        # We don't have bits stored, so let's just show the data we have
        print(f"\n  Example {idx+1}: dH={r['delta_H']}, dc5={r['delta_c5']}, dc7={r['delta_c7']}")
        print(f"    mixed_c3={r['mixed_c3']}, ext={r['ext_scores']}")

# Summary
print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)
print(f"""
At n=8:
- Lambda-preserving (1,1,2,2) reversals: {sample_count} valid
- H-changing: {len(changing)} ({len(changing)/max(1,sample_count)*100:.1f}%)
- c5 invariance: {'HOLDS' if not c5_changed else f'FAILS ({len(c5_changed)} changes)'}
- delta_H values: {dict(Counter(r['delta_H'] for r in changing))}
- All delta_H even: {all_even}
""")

print("Done.")
