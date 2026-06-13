#!/usr/bin/env python3
"""
SC Maximizer Test at n=7 — Fast Version
=========================================
Uses hash-based pruning to speed up the is_self_converse check.

Strategy:
1. Phase 1: Enumerate all 2^21 tournaments, compute (score_seq, H)
2. Phase 2: For each score class, check if max-H achiever is SC
   - Use degree-sequence of the tournament as a quick filter
   - Use a certificate-based approach: if T and T^op have different
     invariants (score seq is same by construction), they're NOT isomorphic
   - Only do full canonical form as last resort, and only for ONE tournament per class

kind-pasteur-2026-03-06-S18
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import (tournament_from_bits, hamiltonian_path_count,
                             opposite_tournament)
from itertools import permutations
from collections import defaultdict
import time

def score_sequence(T):
    n = len(T)
    return tuple(sorted([sum(T[i]) for i in range(n)], reverse=True))

def tournament_invariant(T):
    """Compute a fast isomorphism invariant (not complete but useful for pruning).
    Returns sorted tuple of (out_degree, sorted_neighbor_out_degrees)."""
    n = len(T)
    scores = [sum(T[i]) for i in range(n)]
    # For each vertex: (its score, sorted scores of its out-neighbors)
    profiles = []
    for i in range(n):
        out_nbrs = [j for j in range(n) if j != i and T[i][j]]
        in_nbrs = [j for j in range(n) if j != i and T[j][i]]
        out_scores = tuple(sorted([scores[j] for j in out_nbrs]))
        in_scores = tuple(sorted([scores[j] for j in in_nbrs]))
        profiles.append((scores[i], out_scores, in_scores))
    return tuple(sorted(profiles))

def canonical_form_fast(T):
    """Canonical form with pruning: only check permutations consistent with score partition."""
    n = len(T)
    scores = [sum(T[i]) for i in range(n)]

    # Group vertices by score
    score_groups = defaultdict(list)
    for i in range(n):
        score_groups[scores[i]].append(i)

    # If all scores are distinct, there's only one valid permutation mapping
    sorted_scores = sorted(score_groups.keys(), reverse=True)

    # Generate only permutations that map vertices to same-score positions
    # This is MUCH smaller than n! when scores are spread out
    from itertools import product

    # For each score value, we need to try all permutations of vertices with that score
    groups = [score_groups[s] for s in sorted_scores]
    group_sizes = [len(g) for g in groups]

    # Target positions: first group gets positions 0..len-1, etc.
    target_positions = []
    pos = 0
    for g in groups:
        target_positions.append(list(range(pos, pos + len(g))))
        pos += len(g)

    # Sort vertices within each group for the "target" side
    # We'll fix target positions and permute source vertices
    best = None

    # Generate all valid permutations
    from itertools import permutations as perms
    group_perms = [list(perms(g)) for g in groups]

    for combo in product(*group_perms):
        # Build permutation: perm[target_pos] = source_vertex
        perm = []
        for group_perm in combo:
            perm.extend(group_perm)

        form = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form

    return best

def is_self_converse_fast(T):
    """Check T ≅ T^op with fast pruning."""
    Top = opposite_tournament(T)

    # Quick check: invariants must match (they will for T vs T^op by construction
    # since score sequences are the same, but the fine invariant might differ)
    inv_T = tournament_invariant(T)
    inv_Top = tournament_invariant(Top)
    if inv_T != inv_Top:
        return False

    # Full check with score-based pruning
    return canonical_form_fast(T) == canonical_form_fast(Top)

# ══════════════════════════════════════════════════════════════
# Phase 1: Load results from the previous run or recompute
# ══════════════════════════════════════════════════════════════
n = 7
m = n * (n - 1) // 2
total = 1 << m

print(f"n={n}, m={m}, total={total}")
print("Phase 1: Computing H for all tournaments...")
t0 = time.time()

score_groups = defaultdict(lambda: {'max_h': 0, 'max_bits': [], 'count': 0,
                                     'h_dist': defaultdict(int)})

for bits in range(total):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    ss = score_sequence(T)

    g = score_groups[ss]
    g['count'] += 1
    g['h_dist'][h] += 1
    if h > g['max_h']:
        g['max_h'] = h
        g['max_bits'] = [bits]
    elif h == g['max_h']:
        g['max_bits'].append(bits)

    if bits % 200000 == 0 and bits > 0:
        elapsed = time.time() - t0
        pct = bits / total * 100
        print(f"  {bits}/{total} ({pct:.1f}%), {elapsed:.0f}s elapsed", flush=True)

t1 = time.time()
print(f"Phase 1 complete: {t1-t0:.1f}s, {len(score_groups)} score classes")

# ══════════════════════════════════════════════════════════════
# Phase 2: Check SC for max-H — just check FIRST achiever per class
# ══════════════════════════════════════════════════════════════
print(f"\nPhase 2: SC check for max-H tournaments (one per score class)...")

all_sc_max = True
results = []

for ss in sorted(score_groups.keys()):
    g = score_groups[ss]
    max_h = g['max_h']
    max_bits = g['max_bits']

    # Check first few max-H achievers
    checked = 0
    found_sc = False
    for bits in max_bits[:20]:  # check up to 20
        T = tournament_from_bits(n, bits)
        t_start = time.time()
        sc = is_self_converse_fast(T)
        t_check = time.time() - t_start
        checked += 1
        if sc:
            found_sc = True
            break

    results.append((ss, max_h, len(max_bits), g['count'], found_sc, checked))
    status = "SC" if found_sc else f"no SC in {checked} checked"
    print(f"  {ss}: max H={max_h}, {len(max_bits)} achievers, {status} (check took {t_check:.2f}s)")

    if not found_sc:
        all_sc_max = False
        # Check more if needed
        print(f"    WARNING: checking more achievers...")
        for bits in max_bits[20:200]:
            T = tournament_from_bits(n, bits)
            sc = is_self_converse_fast(T)
            checked += 1
            if sc:
                found_sc = True
                print(f"    Found SC at attempt {checked}")
                break
        if not found_sc:
            print(f"    COUNTEREXAMPLE after {checked} checks!")

t2 = time.time()
print(f"\nPhase 2 complete: {t2-t1:.1f}s")

# ══════════════════════════════════════════════════════════════
# Results
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"RESULTS: SC MAXIMIZER TEST AT n={n}")
print(f"{'='*70}")
print(f"Score classes: {len(results)}")
print(f"SC always achieves max H within score class? {all_sc_max}")

# Global max
global_max_h = max(r[1] for r in results)
print(f"\nGlobal max H = {global_max_h}")
for ss, max_h, n_max, n_total, found_sc, checked in results:
    if max_h == global_max_h:
        print(f"  Achieved by score {ss}, SC found: {found_sc}")

# H distribution summary
print(f"\nH distribution by score class:")
for ss in sorted(score_groups.keys()):
    g = score_groups[ss]
    h_vals = sorted(g['h_dist'].keys())
    print(f"  {ss}: H range [{min(h_vals)}, {max(h_vals)}], max H={g['max_h']}, "
          f"{len(h_vals)} distinct H values, {g['count']} tournaments")

print(f"\nTotal time: {time.time()-t0:.1f}s")
