#!/usr/bin/env python3
"""
TOPOLOGICAL PHASES OF TOURNAMENTS

Three types emerge:
  Type P (point/contractible): β = (1,0,0,...) — most common
  Type C (circle/S^1): β = (1,1,0,...) — has a directed 1-hole
  Type S (3-sphere/S^3): β = (1,0,0,1,...) — has a directed 3-hole

Questions:
1. What invariant determines which phase a tournament is in?
2. Is there a "phase transition" as density parameters change?
3. Does β_3 at n=7 require a specific substructure?
4. Is the Ω structure (transitive triples vs cycles) the key?
"""
import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict
import random
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    path_betti_numbers, enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix, count_3cycles, ham_path_count
)

random.seed(42)

def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

def count_transitive_triples(A, n):
    """Count transitive triples: a→b→c with a→c."""
    count = 0
    for a in range(n):
        for b in range(n):
            if a == b or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[a][c]:
                    count += 1
    return count

def count_3cycles_directed(A, n):
    """Count directed 3-cycles: a→b→c→a."""
    count = 0
    for a in range(n):
        for b in range(n):
            if a == b or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[c][a]:
                    count += 1
    return count // 3

def transitivity_defect(A, n):
    """Number of 2-paths (a,b,c) with a→b→c but c→a (anti-transitive)."""
    count = 0
    for a in range(n):
        for b in range(n):
            if a == b or not A[a][b]: continue
            for c in range(n):
                if c == a or c == b: continue
                if A[b][c] and A[c][a]:
                    count += 1
    return count

# ===== n=7: Full characterization =====
print("=" * 70)
print("n=7: WHAT DETERMINES THE TOPOLOGICAL PHASE?")
print("=" * 70)

n = 7
data = defaultdict(list)
total = 0
for trial in range(500):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=6)

    t3 = count_3cycles(A, n)
    tt = count_transitive_triples(A, n)
    td = transitivity_defect(A, n)
    H = ham_path_count(A, n)
    scores = tuple(sorted([sum(A[i]) for i in range(n)]))

    # dim(Ω_2) = tt (transitive triples)
    # The ratio tt/|A_2| measures "how transitive" the tournament is
    a2_count = sum(1 for a in range(n) for b in range(n) for c in range(n)
                   if a!=b and b!=c and a!=c and A[a][b] and A[b][c])

    phase = "P" if betti[1]==0 and betti[3]==0 else ("C" if betti[1]>0 else "S")
    data[phase].append({
        't3': t3, 'tt': tt, 'td': td, 'H': H,
        'scores': scores, 'a2': a2_count,
        'trans_ratio': tt/a2_count if a2_count > 0 else 0
    })
    total += 1

    if trial % 100 == 99:
        print(f"  ... {trial+1} done", flush=True)

print(f"\nPhase distribution:")
for phase in ['P', 'C', 'S']:
    if data[phase]:
        print(f"  {phase}: {len(data[phase])} ({100*len(data[phase])/total:.1f}%)")

print(f"\nBy phase — mean values:")
for phase in ['P', 'C', 'S']:
    if not data[phase]:
        continue
    d = data[phase]
    print(f"\n  Phase {phase} ({len(d)} tournaments):")
    print(f"    t3:           mean={np.mean([x['t3'] for x in d]):.1f}, range=[{min(x['t3'] for x in d)},{max(x['t3'] for x in d)}]")
    print(f"    trans_triples: mean={np.mean([x['tt'] for x in d]):.1f}")
    print(f"    trans_defect:  mean={np.mean([x['td'] for x in d]):.1f}")
    print(f"    H:            mean={np.mean([x['H'] for x in d]):.1f}")
    print(f"    trans_ratio:  mean={np.mean([x['trans_ratio'] for x in d]):.4f}")
    print(f"    scores:       {Counter([x['scores'] for x in d]).most_common(5)}")

# ===== Key discriminator: transitivity ratio =====
print("\n\n" + "=" * 70)
print("TRANSITIVITY RATIO AS PHASE DISCRIMINATOR")
print("=" * 70)

# Bin by transitivity ratio
bins = defaultdict(lambda: Counter())
for phase in data:
    for x in data[phase]:
        bin_val = round(x['trans_ratio'], 2)
        bins[bin_val][phase] += 1

print(f"\n  trans_ratio -> phase distribution:")
for ratio in sorted(bins.keys()):
    total_bin = sum(bins[ratio].values())
    if total_bin >= 3:
        parts = ", ".join(f"{p}:{bins[ratio][p]}" for p in ['P','C','S'] if bins[ratio][p] > 0)
        print(f"    {ratio:.2f}: {parts} (total {total_bin})")

# ===== t3 vs phase =====
print("\n\n" + "=" * 70)
print("t3 (3-CYCLES) VS PHASE")
print("=" * 70)

t3_by_phase = defaultdict(lambda: Counter())
for phase in data:
    for x in data[phase]:
        t3_by_phase[x['t3']][phase] += 1

for t3 in sorted(t3_by_phase.keys()):
    total_t3 = sum(t3_by_phase[t3].values())
    if total_t3 >= 3:
        parts = ", ".join(f"{p}:{t3_by_phase[t3][p]}" for p in ['P','C','S'] if t3_by_phase[t3][p] > 0)
        print(f"  t3={t3:2d}: {parts}")

# ===== Score sequence vs phase =====
print("\n\n" + "=" * 70)
print("SCORE SEQUENCE VS PHASE")
print("=" * 70)

score_by_phase = defaultdict(lambda: Counter())
for phase in data:
    for x in data[phase]:
        score_by_phase[x['scores']][phase] += 1

for scores in sorted(score_by_phase.keys()):
    total_s = sum(score_by_phase[scores].values())
    if total_s >= 5:
        parts = ", ".join(f"{p}:{score_by_phase[scores][p]}" for p in ['P','C','S'] if score_by_phase[scores][p] > 0)
        print(f"  {scores}: {parts}")

# ===== NEW INVARIANT: count of "fully transitive" 4-subsets =====
print("\n\n" + "=" * 70)
print("NEW INVARIANT: ACYCLIC 4-SUBSETS")
print("=" * 70)

print("\nA 4-subset is 'acyclic' if the induced subtournament has no 3-cycle.")
print("This might predict β_3.")

n = 7
for trial in range(100):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=6)
    phase = "P" if betti[1]==0 and betti[3]==0 else ("C" if betti[1]>0 else "S")

    # Count acyclic 4-subsets
    acyclic_4 = 0
    total_4 = 0
    for quad in combinations(range(n), 4):
        total_4 += 1
        # Check if any triple in quad forms a 3-cycle
        has_cycle = False
        for triple in combinations(quad, 3):
            i, j, k = triple
            if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                has_cycle = True
                break
        if not has_cycle:
            acyclic_4 += 1

    if trial < 20:
        print(f"  Trial {trial}: phase={phase}, acyclic_4={acyclic_4}/{total_4}, β={betti[:4]}", flush=True)

print("\nDone.")
