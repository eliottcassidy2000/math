#!/usr/bin/env python3
"""
sigma* analysis: anti-automorphism action on directed cycle set.
Corrects earlier analysis that used frozenset (wrong for 5-cycles on same vertex set).

kind-pasteur-2026-03-06-S18e
"""
import sys
sys.path.insert(0, r'C:\Users\Eliott\Documents\GitHub\math\03-artifacts\code')
from tournament_lib import (tournament_from_bits, opposite_tournament,
                             find_odd_cycles, conflict_graph, hamiltonian_path_count)
from itertools import permutations

def canonical_form(T):
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

def find_involutory_anti_aut(T):
    n = len(T)
    for perm in permutations(range(n)):
        if any(perm[perm[i]] != i for i in range(n)):
            continue
        ok = True
        for i in range(n):
            for j in range(n):
                if i != j and T[perm[i]][perm[j]] != (1 - T[i][j]):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return perm
    return None

def apply_sigma_to_cycle(sigma, cycle):
    """Apply anti-automorphism to directed cycle: map vertices and reverse direction."""
    mapped = tuple(sigma[v] for v in cycle)
    reversed_mapped = tuple(reversed(mapped))
    # Canonicalize: rotate so minimum vertex is first
    min_val = min(reversed_mapped)
    min_idx = reversed_mapped.index(min_val)
    return reversed_mapped[min_idx:] + reversed_mapped[:min_idx]

def build_sigma_star(sigma, cycles):
    """Build sigma* permutation on cycle indices."""
    cycle_to_idx = {c: i for i, c in enumerate(cycles)}
    sigma_star = []
    for c in cycles:
        mapped = apply_sigma_to_cycle(sigma, c)
        idx = cycle_to_idx.get(mapped)
        if idx is None:
            return None  # sigma(c) not in cycle list
        sigma_star.append(idx)
    return sigma_star

for n in [5, 6]:
    m = n * (n - 1) // 2
    print(f"\n{'='*70}")
    print(f"n={n}: CORRECTED sigma* analysis (directed cycles)")
    print("=" * 70)

    seen = {}
    for bits in range(1 << m):
        T = tournament_from_bits(n, bits)
        cf = canonical_form(T)
        if cf in seen:
            continue
        if not is_self_converse(T):
            seen[cf] = None
            continue

        h = hamiltonian_path_count(T)
        sigma = find_involutory_anti_aut(T)
        if sigma is None:
            seen[cf] = None
            continue

        cycles = find_odd_cycles(T)
        if not cycles:
            seen[cf] = None
            continue

        cg = conflict_graph(cycles)
        sigma_star = build_sigma_star(sigma, cycles)

        if sigma_star is None:
            print(f"  bits={bits} H={h}: sigma* maps outside cycle set!")
            seen[cf] = None
            continue

        # Properties
        is_inv = all(sigma_star[sigma_star[i]] == i for i in range(len(cycles)))

        # Check if sigma* preserves conflict graph
        preserves_cg = True
        for i in range(len(cycles)):
            for j in range(len(cycles)):
                if cg[i][j] != cg[sigma_star[i]][sigma_star[j]]:
                    preserves_cg = False
                    break
            if not preserves_cg:
                break

        fixed = sum(1 for i in range(len(cycles)) if sigma_star[i] == i)
        paired = sum(1 for i in range(len(cycles)) if sigma_star[i] > i)

        # Check: are paired non-adjacent cycles vertex-disjoint?
        paired_disjoint = 0
        paired_adjacent = 0
        for i in range(len(cycles)):
            j = sigma_star[i]
            if j > i:
                if cg[i][j]:
                    paired_adjacent += 1
                else:
                    paired_disjoint += 1

        print(f"  bits={bits} H={h:3d}: |Omega|={len(cycles):2d}, sigma*: inv={is_inv}, "
              f"preserves_cg={preserves_cg}, fixed={fixed}, paired={paired} "
              f"(disjoint={paired_disjoint}, adjacent={paired_adjacent})")

        seen[cf] = True

print("\nDone.")
