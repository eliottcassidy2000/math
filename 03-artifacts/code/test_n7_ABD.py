#!/usr/bin/env python3
"""
test_n7_ABD.py -- Test the A-B-D structure at n=7 (OPEN-Q-010, OPEN-Q-011).

Definitions (per FINAL_FINDINGS.md / OPEN-Q-011):
  For a fixed (T, v):
  A = sum_{P' in Ham(T-v)} #{TypeII positions in P'}
      = sum_{P'} (inshat(v,P')-1)/2  [algebraic identity]
  B = sum_{P'} sum_{TypeII pos j} mu(v, P'[j], P'[j+1])
      = sum over 3-cycles (v,a,b) through v: mu(v,a,b) * #{P': (a,b) consecutive}
  D = sum_{C odd cycle through v} mu(C)  [Claim A's RHS/2]

Claim A: H(T) - H(T-v) = 2*D

At n=6: A - B ~ -5.88, B - D ~ +5.88 (near-cancellation, OPEN-Q-011)
At n=7: OPEN-Q-010/012 ask: do A, B, D still nearly cancel? Is A = D?

Key insight at n=7:
  mu(3-cycle) can be > 1 (4 complement vertices => possible 3-cycles)
  mu(5-cycle) = 1 always (only 2 complement vertices)
  mu(7-cycle) = 1 always (0 complement vertices)

If A = D exactly at n=7, that would be a remarkable theorem.
If A - B = -(B - D) approximately, the near-cancellation extends to n=7.
"""

import sys, random
sys.path.insert(0, '.')
from tournament_lib import (
    random_tournament, hamiltonian_path_count, delete_vertex,
    find_odd_cycles, conflict_graph, independence_poly_at, mu
)
from itertools import permutations


def ham_paths(T):
    """All Hamiltonian paths of T."""
    n = len(T)
    result = []
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[i+1]] for i in range(n-1)):
            result.append(perm)
    return result


def typeii_positions(T, v, path):
    """Count TypeII positions: v->P[j] and P[j+1]->v."""
    return sum(1 for j in range(len(path)-1)
               if T[v][path[j]] and T[path[j+1]][v])


def compute_ABD(T, v):
    """Compute A, B, D for tournament T and vertex v."""
    n = len(T)
    Tv, labels = delete_vertex(T, v)
    others = [u for u in range(n) if u != v]

    # Enumerate Ham paths of T-v
    paths_tv = []
    for perm in permutations(others):
        if all(T[perm[i]][perm[i+1]] for i in range(len(perm)-1)):
            paths_tv.append(perm)

    A = 0
    B = 0
    for path in paths_tv:
        for j in range(len(path)-1):
            a, b = path[j], path[j+1]
            if T[v][a] and T[b][v]:  # TypeII: v->a, b->v; and a->b (given by valid path)
                A += 1
                B += mu(T, v, (v, a, b))

    # D: sum of mu(C) over all odd cycles through v
    D = 0
    all_cycles = find_odd_cycles(T)
    for cyc in all_cycles:
        if v in cyc:
            D += mu(T, v, cyc)

    return A, B, D, len(paths_tv)


def analyze_n7(n_trials=150, seed=42):
    rng = random.Random(seed)
    n = 7
    print(f"=== n={n} Analysis: A, B, D values ===")
    print(f"Testing {n_trials} random tournaments x {n} vertices = up to {n_trials*n} (T,v) pairs")
    print()

    results = []
    for _ in range(n_trials):
        T = random_tournament(n, rng)
        for v in range(n):
            A, B, D, npaths = compute_ABD(T, v)
            H_T = hamiltonian_path_count(T)
            H_Tv, _ = delete_vertex(T, v)
            H_Tv_count = hamiltonian_path_count(H_Tv)
            claim_a_lhs = H_T - H_Tv_count
            results.append((A, B, D, claim_a_lhs // 2, A-B, B-D, A-D))

    import statistics
    As  = [r[0] for r in results]
    Bs  = [r[1] for r in results]
    Ds  = [r[2] for r in results]
    ADs = [r[6] for r in results]
    ABs = [r[4] for r in results]
    BDs = [r[5] for r in results]
    claim_a_checks = [r[3] for r in results]

    print(f"Pairs analyzed: {len(results)}")
    print()
    print(f"A (TypeII count):  mean={statistics.mean(As):.3f}, std={statistics.stdev(As):.3f}")
    print(f"B (mu-weighted 3): mean={statistics.mean(Bs):.3f}, std={statistics.stdev(Bs):.3f}")
    print(f"D (Claim A RHS/2): mean={statistics.mean(Ds):.3f}, std={statistics.stdev(Ds):.3f}")
    print()
    print(f"A - B: mean={statistics.mean(ABs):.4f}, min={min(ABs)}, max={max(ABs)}")
    print(f"B - D: mean={statistics.mean(BDs):.4f}, min={min(BDs)}, max={max(BDs)}")
    print(f"A - D: mean={statistics.mean(ADs):.4f}, min={min(ADs)}, max={max(ADs)}")
    print()

    # Check A = D in all cases
    a_eq_d = sum(1 for r in results if r[6] == 0)
    print(f"A == D (exact): {a_eq_d}/{len(results)} ({100*a_eq_d/len(results):.1f}%)")
    print(f"A == B (exact): {sum(1 for r in results if r[4]==0)}/{len(results)}")
    print(f"B == D (exact): {sum(1 for r in results if r[5]==0)}/{len(results)}")
    print()

    # Near-cancellation check
    near_cancel = sum(1 for r in results if abs(r[4] + r[5]) <= 1)
    print(f"Near-cancellation |A-B + B-D| <= 1: {near_cancel}/{len(results)}")
    # Equivalently, |A-D| <= 1
    print(f"|A-D| <= 1: {sum(1 for r in results if abs(r[6])<=1)}/{len(results)}")
    print()

    if a_eq_d == len(results):
        print("THEOREM CANDIDATE: A = D exactly for ALL (T,v) pairs at n=7!")
        print("This would mean: sum_{P'} TypeII_count = sum_{odd cycles C through v} mu(C)")
        print("Equivalently: the naive TypeII count formula is EXACT at n=7.")
    else:
        dist = {}
        for r in results:
            dist[r[6]] = dist.get(r[6], 0) + 1
        print(f"A - D distribution: {sorted(dist.items())}")


if __name__ == "__main__":
    analyze_n7(n_trials=150)
