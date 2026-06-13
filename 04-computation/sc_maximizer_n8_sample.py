#!/usr/bin/env python3
"""
SC Maximizer sampling test at n=8.

OPEN-Q-016: Within each self-complementary score class, is max H always
achieved by an SC tournament?

At n=8 exhaustive is infeasible (2^28 = 268M). We sample and check.

Strategy:
1. Sample 200K random tournaments
2. Group by score sequence
3. For each score class, track max H and whether it's SC
4. Focus on self-complementary score classes (palindromic scores)

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import random_tournament, hamiltonian_path_count, opposite_tournament
from itertools import permutations
from collections import defaultdict
import random

def score_sequence(T):
    n = len(T)
    return tuple(sorted([sum(T[i]) for i in range(n)], reverse=True))

def is_sc_score(ss):
    """Check if score sequence is self-complementary (palindromic)."""
    n = len(ss)
    comp = tuple(sorted([n-1-s for s in ss], reverse=True))
    return ss == comp

def is_self_converse_fast(T):
    """Quick check: is T isomorphic to its opposite?
    Uses invariant-based filtering first."""
    n = len(T)
    T_opp = [[1 - T[i][j] if i != j else 0 for j in range(n)] for i in range(n)]

    # Quick invariant check: sorted score sequence (already same for SC score class)
    # Check sorted (out-degree, in-degree) multisets
    scores_T = tuple(sorted([sum(T[i]) for i in range(n)]))
    scores_opp = tuple(sorted([sum(T_opp[i]) for i in range(n)]))
    if scores_T != scores_opp:
        return False

    # Check 3-cycle count (invariant)
    c3_T = 0
    c3_opp = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if T[i][j] and T[j][k] and T[k][i]:
                    c3_T += 1
                elif T[i][k] and T[k][j] and T[j][i]:
                    c3_T += 1
                if T_opp[i][j] and T_opp[j][k] and T_opp[k][i]:
                    c3_opp += 1
                elif T_opp[i][k] and T_opp[k][j] and T_opp[j][i]:
                    c3_opp += 1
    if c3_T != c3_opp:
        return False

    # Full isomorphism check (brute force for n=8)
    # Use sorted score profile as vertex invariant
    # Score profile: (out_degree, sorted_neighbor_scores)
    def vertex_profile(M):
        profiles = []
        for i in range(n):
            out_deg = sum(M[i])
            # Sorted out-degrees of out-neighbors
            out_nbrs = tuple(sorted([sum(M[j]) for j in range(n) if j != i and M[i][j]]))
            in_nbrs = tuple(sorted([sum(M[j]) for j in range(n) if j != i and M[j][i]]))
            profiles.append((out_deg, out_nbrs, in_nbrs))
        return profiles

    prof_T = vertex_profile(T)
    prof_opp = vertex_profile(T_opp)

    # Try to find a permutation mapping T to T_opp
    # Group vertices by profile
    from collections import Counter
    if Counter([p for p in prof_T]) != Counter([p for p in prof_opp]):
        return False

    # For each vertex in T, possible mappings to T_opp
    # Group by profile
    profile_groups_T = defaultdict(list)
    profile_groups_opp = defaultdict(list)
    for i, p in enumerate(prof_T):
        profile_groups_T[p].append(i)
    for i, p in enumerate(prof_opp):
        profile_groups_opp[p].append(i)

    # Generate candidate mappings using profile matching
    # For small groups, enumerate. For large groups, sample.
    profiles = sorted(profile_groups_T.keys())

    def check_mapping(perm):
        for i in range(n):
            for j in range(n):
                if i != j and T[i][j] != T_opp[perm[i]][perm[j]]:
                    return False
        return True

    # Simple approach: try random permutations within profile groups
    from itertools import product

    groups_T = [profile_groups_T[p] for p in profiles]
    groups_opp = [profile_groups_opp[p] for p in profiles]

    # If any group has > 6 elements, just sample
    max_group = max(len(g) for g in groups_T) if groups_T else 0

    if max_group <= 4:
        # Enumerate all mappings
        from itertools import permutations as perms
        options = []
        for gt, go in zip(groups_T, groups_opp):
            options.append(list(perms(go)))

        from itertools import product as iprod
        for combo in iprod(*options):
            perm = [0] * n
            for gt, mapping in zip(groups_T, combo):
                for i, j in zip(gt, mapping):
                    perm[i] = j
            if check_mapping(perm):
                return True
        return False
    else:
        # Sample random mappings
        for _ in range(1000):
            perm = [0] * n
            for gt, go in zip(groups_T, groups_opp):
                go_shuffled = list(go)
                random.shuffle(go_shuffled)
                for i, j in zip(gt, go_shuffled):
                    perm[i] = j
            if check_mapping(perm):
                return True
        return False  # Probably not SC, but not certain

print("=" * 70)
print("SC MAXIMIZER TEST AT n=8 (SAMPLING)")
print("=" * 70)

n = 8
num_samples = 100000

# Track per score class: max_H, whether SC achieves it
class_data = defaultdict(lambda: {'max_H': 0, 'max_H_sc': 0, 'max_H_is_sc': None,
                                    'count': 0, 'sc_count': 0})

for trial in range(num_samples):
    T = random_tournament(n)
    ss = score_sequence(T)
    H = hamiltonian_path_count(T)

    info = class_data[ss]
    info['count'] += 1

    if H > info['max_H']:
        info['max_H'] = H
        # Check if this tournament is SC
        if is_sc_score(ss):
            info['max_H_is_sc'] = is_self_converse_fast(T)

    # Also check if any SC tournament in this class has high H
    if is_sc_score(ss) and is_self_converse_fast(T):
        info['sc_count'] += 1
        if H > info['max_H_sc']:
            info['max_H_sc'] = H

    if trial > 0 and trial % 20000 == 0:
        print(f"  {trial}/{num_samples} done...", flush=True)

print(f"\nProcessed {num_samples} tournaments")
print(f"Score classes found: {len(class_data)}")

# Report on self-complementary score classes
print(f"\n{'='*60}")
print("SELF-COMPLEMENTARY SCORE CLASSES")
print(f"{'='*60}")

sc_classes = [(ss, d) for ss, d in class_data.items() if is_sc_score(ss)]
sc_classes.sort(key=lambda x: -x[1]['count'])

violations = 0
for ss, d in sc_classes:
    if d['sc_count'] == 0:
        continue  # Not enough SC samples

    gap = d['max_H'] - d['max_H_sc']
    is_violation = gap > 0
    if is_violation:
        violations += 1

    if d['count'] >= 20 or is_violation:
        status = "VIOLATION!" if is_violation else "OK"
        print(f"  {ss}: count={d['count']}, sc_count={d['sc_count']}, "
              f"max_H={d['max_H']}, max_H_sc={d['max_H_sc']}, gap={gap} {status}")

print(f"\nViolations (max_H from NSC > max_H from SC): {violations}")
if violations == 0:
    print("SC MAXIMIZER CONJECTURE HOLDS in this sample!")

# Also check: what's the global max H at n=8?
all_H = []
for ss, d in class_data.items():
    all_H.append(d['max_H'])

print(f"\nGlobal max H observed: {max(all_H)}")
print(f"(Known: a(8) = 661)")

print(f"\n{'='*70}")
print("DONE")
print("=" * 70)
