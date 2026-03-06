#!/usr/bin/env python3
"""
SC Maximizer Test at n=7 (Exhaustive)
======================================
Test: Within each score sequence class, is max H always achieved by a
self-converse tournament?

Strategy (for speed):
1. Enumerate all 2^21 = 2,097,152 tournaments at n=7
2. Compute score_seq and H for each (both fast at n=7)
3. Group by score_seq, find max H per group
4. Only check is_self_converse for max-H tournaments (expensive but rare)

kind-pasteur-2026-03-06-S18
"""

import sys
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code')
from tournament_lib import (tournament_from_bits, hamiltonian_path_count,
                             opposite_tournament)
from itertools import permutations
from collections import defaultdict
import time

def score_sequence(T):
    n = len(T)
    return tuple(sorted([sum(T[i]) for i in range(n)], reverse=True))

def canonical_form(T):
    """Canonical form via brute-force permutation enumeration."""
    n = len(T)
    best = None
    for perm in permutations(range(n)):
        form = tuple(T[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

def is_self_converse(T):
    Top = opposite_tournament(T)
    return canonical_form(T) == canonical_form(Top)

# ══════════════════════════════════════════════════════════════
# Phase 1: Enumerate all tournaments, group by score sequence
# ══════════════════════════════════════════════════════════════
n = 7
m = n * (n - 1) // 2  # 21
total = 1 << m  # 2,097,152

print(f"n={n}, m={m}, total={total}")
print("Phase 1: Computing H for all tournaments...")
t0 = time.time()

# score_seq -> {max_h, max_bits_list, all_h_values}
score_groups = defaultdict(lambda: {'max_h': 0, 'max_bits': [], 'count': 0,
                                      'sc_checked': False})

for bits in range(total):
    T = tournament_from_bits(n, bits)
    h = hamiltonian_path_count(T)
    ss = score_sequence(T)

    g = score_groups[ss]
    g['count'] += 1
    if h > g['max_h']:
        g['max_h'] = h
        g['max_bits'] = [bits]
    elif h == g['max_h']:
        g['max_bits'].append(bits)

    if bits % 100000 == 0 and bits > 0:
        elapsed = time.time() - t0
        pct = bits / total * 100
        eta = elapsed / pct * (100 - pct)
        print(f"  {bits}/{total} ({pct:.1f}%), {elapsed:.0f}s elapsed, ~{eta:.0f}s remaining")

t1 = time.time()
print(f"Phase 1 complete: {t1-t0:.1f}s, {len(score_groups)} score classes")

# ══════════════════════════════════════════════════════════════
# Phase 2: Check is_self_converse for max-H tournaments only
# ══════════════════════════════════════════════════════════════
print(f"\nPhase 2: Checking self-converse for max-H tournaments...")

# How many max-H tournaments total?
total_max = sum(len(g['max_bits']) for g in score_groups.values())
print(f"  Total max-H tournaments to check: {total_max}")

all_sc_max = True
results = []

for ss in sorted(score_groups.keys()):
    g = score_groups[ss]
    max_h = g['max_h']
    max_bits = g['max_bits']

    # Check if ANY max-H tournament is SC
    any_sc = False
    any_nsc = False
    for bits in max_bits:
        T = tournament_from_bits(n, bits)
        sc = is_self_converse(T)
        if sc:
            any_sc = True
        else:
            any_nsc = True
        # Once we find one SC, we can stop checking if we just want "any SC"
        if any_sc:
            break

    results.append((ss, max_h, len(max_bits), g['count'], any_sc, any_nsc))

    if not any_sc:
        all_sc_max = False
        print(f"  COUNTEREXAMPLE: scores={ss}, max H={max_h}, {len(max_bits)} achievers, NONE SC!")

t2 = time.time()
print(f"Phase 2 complete: {t2-t1:.1f}s")

# ══════════════════════════════════════════════════════════════
# Results
# ══════════════════════════════════════════════════════════════
print(f"\n{'='*70}")
print(f"RESULTS: SC MAXIMIZER TEST AT n={n}")
print(f"{'='*70}")
print(f"Score classes: {len(results)}")
print(f"SC always achieves max H within score class? {all_sc_max}")

print(f"\nPer score class:")
for ss, max_h, n_max, n_total, any_sc, any_nsc in results:
    sc_label = "SC" if any_sc else "NSC!!!"
    print(f"  {ss}: max H={max_h:4d}, {n_max:5d} achievers / {n_total:6d} total, max is {sc_label}")

# ══════════════════════════════════════════════════════════════
# Phase 3: Also check if global max is SC
# ══════════════════════════════════════════════════════════════
global_max_h = max(g['max_h'] for g in score_groups.values())
global_max_ss = [ss for ss, g in score_groups.items() if g['max_h'] == global_max_h]
print(f"\nGlobal max H = {global_max_h}, achieved by score classes: {global_max_ss}")

# Check all global max achievers
for ss in global_max_ss:
    g = score_groups[ss]
    for bits in g['max_bits'][:5]:  # check first 5
        T = tournament_from_bits(n, bits)
        sc = is_self_converse(T)
        print(f"  bits={bits}, SC={sc}")

print(f"\nTotal time: {time.time()-t0:.1f}s")
