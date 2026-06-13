"""
SC Maximizer: algebraic proof attempt via Z_2 action on conflict graph.

Key insight to prove: For SC tournament T* with anti-aut sigma,
the Z_2 action on Omega(T*) creates cycle-pair structure that
boosts I(Omega, 2) beyond what any NSC tournament can achieve.

Strategy:
1. Show sigma acts as automorphism of Omega(T*)
2. Use the Z_2 action to decompose I(Omega, 2) into orbit contributions
3. Show the orbit structure maximizes I at x=2

Also: test the "T-value" = alpha_1 + 2*alpha_2 + 4*alpha_3 + ... = (H-1)/2
across SC vs NSC tournaments at n=6,7 to understand the mechanism.

kind-pasteur-2026-03-25-S1
"""
import sys
import time
from collections import defaultdict
from itertools import combinations, permutations
from math import factorial, comb

sys.stdout.reconfigure(line_buffering=True)

def tournament_from_bits(n, bits):
    """Create adjacency matrix from upper-triangle bits."""
    T = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    return T

def score_seq(T, n):
    return tuple(sorted(sum(T[i]) for i in range(n)))

def is_sc_score(s, n):
    """Check if score sequence is self-complementary."""
    return all(s[i] + s[n-1-i] == n-1 for i in range(n))

def count_ham_paths(T, n):
    """Held-Karp DP."""
    dp = {}
    for v in range(n):
        dp[(1<<v, v)] = 1
    for sz in range(2, n+1):
        for mask in range(1<<n):
            if bin(mask).count('1') != sz:
                continue
            for v in range(n):
                if not (mask & (1<<v)):
                    continue
                prev = mask ^ (1<<v)
                for u in range(n):
                    if u == v or not (prev & (1<<u)):
                        continue
                    if T[u][v]:
                        dp[(mask,v)] = dp.get((mask,v), 0) + dp.get((prev,u), 0)
    full = (1<<n)-1
    return sum(dp.get((full,v), 0) for v in range(n))

def find_odd_cycles(T, n, max_len=None):
    """Find all directed odd cycles."""
    if max_len is None:
        max_len = n
    cycles = []
    for length in range(3, max_len+1, 2):
        for verts in combinations(range(n), length):
            # Check all cyclic orderings
            vset = set(verts)
            for perm in permutations(verts):
                is_cycle = True
                for i in range(length):
                    if not T[perm[i]][perm[(i+1)%length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    # Normalize: start from min vertex
                    min_idx = perm.index(min(perm))
                    normalized = perm[min_idx:] + perm[:min_idx]
                    cycles.append(normalized)
                    break  # Only need one direction per vertex set... wait, need all
            # Actually, count ALL directed cycles on this vertex set
            for perm in permutations(verts):
                is_cycle = True
                for i in range(length):
                    if not T[perm[i]][perm[(i+1)%length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    min_idx = perm.index(min(perm))
                    normalized = perm[min_idx:] + perm[:min_idx]
                    if normalized not in cycles:
                        cycles.append(normalized)
    # Remove duplicates
    return list(set(cycles))

def find_directed_odd_cycles(T, n, max_len=None):
    """Find all directed odd cycles more efficiently."""
    if max_len is None:
        max_len = n
    all_cycles = []
    for length in range(3, max_len + 1, 2):
        for verts in combinations(range(n), length):
            vlist = list(verts)
            # Try all cyclic permutations starting from vlist[0]
            for perm in permutations(vlist[1:]):
                cycle = (vlist[0],) + perm
                is_cycle = True
                for i in range(length):
                    if not T[cycle[i]][cycle[(i+1) % length]]:
                        is_cycle = False
                        break
                if is_cycle:
                    all_cycles.append(cycle)
    return all_cycles

def compute_independence_poly(cycles, n):
    """Compute independence polynomial coefficients of conflict graph.
    Two cycles conflict if they share a vertex."""
    # Build conflict graph
    num_cycles = len(cycles)
    cycle_verts = [set(c) for c in cycles]

    # Find all independent sets
    alpha = [0] * (num_cycles + 1)
    alpha[0] = 1

    # For small number of cycles, enumerate subsets
    if num_cycles <= 25:
        for mask in range(1 << num_cycles):
            selected = [i for i in range(num_cycles) if (mask >> i) & 1]
            k = len(selected)
            # Check independence: no two share a vertex
            independent = True
            for a in range(len(selected)):
                for b in range(a+1, len(selected)):
                    if cycle_verts[selected[a]] & cycle_verts[selected[b]]:
                        independent = False
                        break
                if not independent:
                    break
            if independent:
                alpha[k] += 1
    else:
        # Use greedy/backtracking for larger cases
        # Simple recursive approach
        def count_indep(idx, used_verts, k):
            alpha[k] += 1
            for i in range(idx, num_cycles):
                if not (used_verts & cycle_verts[i]):
                    count_indep(i + 1, used_verts | cycle_verts[i], k + 1)
        count_indep(0, set(), 0)

    return alpha

def is_sc_tournament(T, n):
    """Check if tournament is self-complementary (has anti-automorphism)."""
    # T^op[i][j] = T[j][i]
    Top = [[T[j][i] for j in range(n)] for i in range(n)]
    # Check if T and T^op are isomorphic
    for perm in permutations(range(n)):
        match = True
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                if Top[perm[i]][perm[j]] != T[i][j]:
                    match = False
                    break
            if not match:
                break
        if match:
            return True
    return False

def analyze_score_class(n, target_score, sample_limit=None):
    """Analyze all tournaments with a given score sequence."""
    total_bits = comb(n, 2)
    results = []
    count = 0

    for bits in range(1 << total_bits):
        T = tournament_from_bits(n, bits)
        s = score_seq(T, n)
        if s != target_score:
            continue
        count += 1
        if sample_limit and count > sample_limit:
            break

        H = count_ham_paths(T, n)
        sc = is_sc_tournament(T, n)

        results.append({
            'bits': bits,
            'H': H,
            'sc': sc,
        })

    return results

# ==== MAIN ANALYSIS ====

print("=" * 60)
print("SC MAXIMIZER PROOF ANALYSIS")
print("kind-pasteur-2026-03-25-S1")
print("=" * 60)

# n=5: exhaustive
n = 5
print(f"\n{'#'*60}")
print(f"  n = {n}")
print(f"{'#'*60}")

total_bits = comb(n, 2)
sc_scores = set()
all_scores = defaultdict(list)

t0 = time.time()
for bits in range(1 << total_bits):
    T = tournament_from_bits(n, bits)
    s = score_seq(T, n)
    H = count_ham_paths(T, n)
    all_scores[s].append((bits, H))

print(f"  Enumerated {1 << total_bits} tournaments in {time.time()-t0:.1f}s")

# Find SC score sequences
for s in sorted(all_scores.keys()):
    if is_sc_score(s, n):
        sc_scores.add(s)

print(f"  SC score sequences: {sorted(sc_scores)}")

for s in sorted(sc_scores):
    entries = all_scores[s]
    max_H = max(H for _, H in entries)

    # Check which maximizers are SC
    maximizers = [(bits, H) for bits, H in entries if H == max_H]
    sc_maximizers = 0
    nsc_maximizers = 0

    for bits, H in maximizers:
        T = tournament_from_bits(n, bits)
        if is_sc_tournament(T, n):
            sc_maximizers += 1
        else:
            nsc_maximizers += 1

    # Also find max H among SC and NSC separately
    sc_max = 0
    nsc_max = 0
    for bits, H in entries:
        T = tournament_from_bits(n, bits)
        if is_sc_tournament(T, n):
            sc_max = max(sc_max, H)
        else:
            nsc_max = max(nsc_max, H)

    print(f"\n  Score {s}: {len(entries)} tournaments, max H = {max_H}")
    print(f"    SC maximizers: {sc_maximizers}, NSC maximizers: {nsc_maximizers}")
    print(f"    SC max H = {sc_max}, NSC max H = {nsc_max}")
    if sc_max >= nsc_max:
        print(f"    [OK] SC achieves max")
    else:
        print(f"    [FAIL] NSC beats SC!")

# n=6: exhaustive but slower
n = 6
print(f"\n{'#'*60}")
print(f"  n = {n}")
print(f"{'#'*60}")

total_bits = comb(n, 2)
all_scores_6 = defaultdict(list)

t0 = time.time()
for bits in range(1 << total_bits):
    T = tournament_from_bits(n, bits)
    s = score_seq(T, n)
    H = count_ham_paths(T, n)
    all_scores_6[s].append((bits, H))

print(f"  Enumerated {1 << total_bits} tournaments in {time.time()-t0:.1f}s")

sc_scores_6 = set()
for s in sorted(all_scores_6.keys()):
    if is_sc_score(s, n):
        sc_scores_6.add(s)

print(f"  SC score sequences ({len(sc_scores_6)} total): {sorted(sc_scores_6)}")

for s in sorted(sc_scores_6):
    entries = all_scores_6[s]
    max_H = max(H for _, H in entries)

    # Sample to check SC property (full is_sc check is expensive)
    # Just check maximizers
    maximizers = [(bits, H) for bits, H in entries if H == max_H]

    sc_max = 0
    nsc_max = 0
    sc_count = 0
    nsc_count = 0

    # For speed, only check SC status for maximizers and a sample
    for bits, H in maximizers[:50]:  # Check up to 50 maximizers
        T = tournament_from_bits(n, bits)
        if is_sc_tournament(T, n):
            sc_count += 1
            sc_max = max(sc_max, H)
        else:
            nsc_count += 1
            nsc_max = max(nsc_max, H)

    # Also check the highest-H NSC
    sorted_entries = sorted(entries, key=lambda x: -x[1])
    for bits, H in sorted_entries[:100]:
        T = tournament_from_bits(n, bits)
        if not is_sc_tournament(T, n):
            nsc_max = max(nsc_max, H)
            break

    verdict = "[OK]" if sc_max >= nsc_max and sc_count > 0 else "[CHECK]"
    print(f"  Score {s}: {len(entries)} tours, max H={max_H}, SC in max: {sc_count}/{sc_count+nsc_count} {verdict}")

# n=7: sample-based (too many tournaments for exhaustive)
n = 7
print(f"\n{'#'*60}")
print(f"  n = {n} (sampling)")
print(f"{'#'*60}")

import random
random.seed(42)

# SC score sequences at n=7
sc_scores_7 = []
for s in combinations(range(n), (n-1)//2):
    # Generate all score sequences... too complex. Just use known SC scores.
    pass

# Known: regular score (3,3,3,3,3,3,3) is SC
# Test: for regular tournaments, is max H always SC?
print("\n  Testing regular score (3,3,3,3,3,3,3):")
total_bits = comb(n, 2)  # 21 bits = 2M tournaments

# Sample and check
t0 = time.time()
best_H = 0
best_bits = 0
sample_count = 0
sc_best = 0
nsc_best = 0
H_dist = defaultdict(int)

for _ in range(50000):
    bits = random.randint(0, (1 << total_bits) - 1)
    T = tournament_from_bits(n, bits)
    s = score_seq(T, n)
    if s != (3,3,3,3,3,3,3):
        continue
    sample_count += 1
    H = count_ham_paths(T, n)
    H_dist[H] += 1
    if H > best_H:
        best_H = H
        best_bits = bits

print(f"  Sampled {sample_count} regular tournaments in {time.time()-t0:.1f}s")
print(f"  Best H = {best_H}")
print(f"  H distribution (top 5): {sorted(H_dist.items(), key=lambda x: -x[1])[:5]}")

# Check if the best is SC
if sample_count > 0:
    T = tournament_from_bits(n, best_bits)
    sc = is_sc_tournament(T, n)
    print(f"  Best tournament SC: {sc}")

# Now the key analysis: for the maximizers at n=5,6, compute the Z_2
# orbit structure on Omega and show how it creates extra independence.

print(f"\n{'='*60}")
print("  Z_2 ORBIT ANALYSIS ON OMEGA")
print(f"{'='*60}")

# At n=5, take the SC maximizer (H=15, Paley T_5 = circulant {1,2})
n = 5
# Paley T_5: connection set {1,2} on Z_5
T5 = [[0]*5 for _ in range(5)]
QR = {1, 2}  # Quadratic residues mod 5 (since 1^2=1, 2^2=4=4 mod 5... hmm)
# Actually QR mod 5 = {1, 4} (1^2=1, 2^2=4)
# For p=5 ≡ 1 mod 4, Paley tournament doesn't exist in the standard sense
# Use circulant with S = {1, 2} instead
S = {1, 2}
for i in range(5):
    for j in range(5):
        if i != j and ((j - i) % 5) in S:
            T5[i][j] = 1

H5 = count_ham_paths(T5, 5)
print(f"\n  Circulant T_5 with S={{1,2}}: H = {H5}")

cycles5 = find_directed_odd_cycles(T5, 5)
print(f"  Directed odd cycles: {len(cycles5)}")
for c in sorted(cycles5, key=lambda x: (len(x), x)):
    print(f"    {c}")

# Find the anti-automorphism sigma
print("\n  Anti-automorphisms of T_5:")
Top5 = [[T5[j][i] for j in range(5)] for i in range(5)]
for perm in permutations(range(5)):
    match = True
    for i in range(5):
        for j in range(5):
            if i == j: continue
            if Top5[perm[i]][perm[j]] != T5[i][j]:
                match = False
                break
        if not match:
            break
    if match:
        # Check if involution
        is_invol = all(perm[perm[i]] == i for i in range(5))
        print(f"    sigma = {perm}, involution: {is_invol}")

print("\nDONE.")
