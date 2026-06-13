#!/usr/bin/env python3
"""
INV-004c: Algebraic formula for delta_H under single-vertex cut-flip.

Key finding: delta_H(phi_v, T) is fully determined by the "link" of v
(the sub-tournament T\v together with the arc directions from v).

This means: delta_H(phi_v, T) = H(T_flipped) - H(T) can be computed
from local information. Can we find an explicit formula?

Observation: phi_v reverses all arcs through v. So the sub-tournament T\v
is unchanged, and only the in/out partition of N(v) flips.

This means: H(T) - H(phi_v(T)) = H(T; v's arcs) - H(T; v's arcs flipped)
where "v's arcs" refers to the direction of arcs incident to v.

Key structural insight: A Hamiltonian path uses v exactly once. The path
enters v from some predecessor and exits to some successor (unless v is
first or last). Flipping v's arcs swaps which arcs go in/out.

Let T' = sub-tournament on V\{v}. For each pair (a,b) of positions
where v can be inserted, v can be placed between a and b in the path
if a->v and v->b (i.e., a beats v and v beats b). After flipping,
it's v->a and b->v... wait, that's wrong. Let me think more carefully.

Author: opus-2026-03-06-S18
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '03-artifacts', 'code'))
from tournament_lib import tournament_from_bits, hamiltonian_path_count
from itertools import permutations
from collections import defaultdict

def get_ham_paths(T):
    n = len(T)
    paths = []
    for perm in permutations(range(n)):
        if all(T[perm[i]][perm[i+1]] for i in range(n-1)):
            paths.append(perm)
    return paths

def cut_flip_vertex(T, v):
    """Flip all arcs through vertex v."""
    n = len(T)
    T_new = [row[:] for row in T]
    for u in range(n):
        if u == v:
            continue
        T_new[v][u] = 1 - T[v][u]
        T_new[u][v] = 1 - T[u][v]
    return T_new

print("=" * 70)
print("CUT-FLIP FORMULA SEARCH")
print("=" * 70)

# Detailed path analysis: how does phi_v affect individual paths?
print("\n--- Path-level analysis at n=4 ---")
n = 4
m = n * (n - 1) // 2

for bits in [4, 6]:  # H=5 cases with interesting delta_H
    T = tournament_from_bits(n, bits)
    paths = get_ham_paths(T)
    H = len(paths)
    scores = [sum(T[i][j] for j in range(n) if j != i) for i in range(n)]

    print(f"\n  bits={bits}, H={H}, scores={scores}")
    print(f"  Paths: {paths}")

    for v in range(n):
        T_flip = cut_flip_vertex(T, v)
        paths_flip = get_ham_paths(T_flip)
        H_flip = len(paths_flip)

        # Which paths survived, which are new?
        survived = [p for p in paths if p in paths_flip]
        lost = [p for p in paths if p not in paths_flip]
        gained = [p for p in paths_flip if p not in paths]

        print(f"\n  phi_{v}: H {H} -> {H_flip} (delta={H_flip-H})")
        print(f"    Survived: {survived}")
        print(f"    Lost: {lost}")
        print(f"    Gained: {gained}")

        # For each lost/gained path, what's v's position?
        for p in lost:
            pos = p.index(v)
            print(f"    Lost {p}: v at position {pos}")
        for p in gained:
            pos = p.index(v)
            print(f"    Gained {p}: v at position {pos}")

# Key structural observation: paths through v can be decomposed
# v at position 0: v -> ... (v must beat its successor)
# v at position k: ...->a->v->b->... (a beats v, v beats b)
# v at position n-1: ...->v (v's predecessor beats v)

print("\n\n--- Insertion decomposition ---")
print("  H(T) = sum over positions k of: (# ways to have v at position k)")
print("  After phi_v, the in/out neighborhoods swap.")
print("  So paths with v at position k that required v->b now need b->v.")

# Compute H via insertion at each position
n = 5
m = n * (n - 1) // 2

formula_works = 0
formula_fails = 0

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)
    H = hamiltonian_path_count(T)

    for v in range(n):
        T_flip = cut_flip_vertex(T, v)
        H_flip = hamiltonian_path_count(T_flip)
        dH = H_flip - H

        # Compute candidate formula: based on the idea that
        # phi_v maps path (..a->v->b..) to (..b->v->a..) in the reverse
        # Actually let's count paths by v's position
        paths = get_ham_paths(T)
        paths_flip = get_ham_paths(T_flip)

        # Count paths by position of v
        pos_count = defaultdict(int)
        pos_count_flip = defaultdict(int)
        for p in paths:
            pos_count[p.index(v)] += 1
        for p in paths_flip:
            pos_count_flip[p.index(v)] += 1

        # After phi_v, a path with v at position k in T corresponds to
        # v at position (n-1-k) in phi_v(T) (because in/out swap)
        # So pos_count_flip[k] should equal pos_count[n-1-k]
        formula_ok = True
        for k in range(n):
            if pos_count_flip[k] != pos_count[n-1-k]:
                formula_ok = False
                break

        if formula_ok:
            formula_works += 1
        else:
            formula_fails += 1
            if formula_fails <= 5:
                print(f"  FAIL: bits={bits}, v={v}")
                print(f"    pos_count = {dict(pos_count)}")
                print(f"    pos_count_flip = {dict(pos_count_flip)}")
                print(f"    expected = {dict((k, pos_count[n-1-k]) for k in range(n))}")

print(f"\n  Position-reversal formula: {formula_works} pass, {formula_fails} fail")

if formula_fails == 0:
    print("  => phi_v maps paths with v at position k to paths with v at position n-1-k!")
    print("  => delta_H = sum_k (pos_count[n-1-k] - pos_count[k])")
    print("             = sum_k pos_count[n-1-k] - sum_k pos_count[k] = 0")
    print("  Wait, that gives delta_H = 0 always, which is wrong.")
    print("  The mapping is not 1-1; it's a count identity.")

# Let me try a different approach: what if we decompose H by
# the position of v and the arc direction at each boundary?
print("\n\n--- Boundary arc analysis ---")
n = 4
m = n * (n - 1) // 2

for bits in [4]:
    T = tournament_from_bits(n, bits)
    paths = get_ham_paths(T)
    H = len(paths)

    for v in range(n):
        # For each path, record:
        # - Position of v
        # - Left neighbor (if exists) and whether left->v
        # - Right neighbor (if exists) and whether v->right
        for p in paths:
            pos = p.index(v)
            left = p[pos-1] if pos > 0 else None
            right = p[pos+1] if pos < n-1 else None
            left_arc = f"{left}->v" if left is not None else "START"
            right_arc = f"v->{right}" if right is not None else "END"
            print(f"  bits={bits}, path={p}: v={v} at pos {pos}, {left_arc}, {right_arc}")

# New idea: H(T) via deletion-contraction on vertex v
print("\n\n--- H via position counting ---")
print("  Let h_k(T,v) = # Hamiltonian paths with v at position k (0-indexed)")
print("  Then H(T) = sum_k h_k(T,v)")
print("  After phi_v: h_k(phi_v(T), v) = h_{n-1-k}(T, v)?  (hypothesis)")

n = 5
m = n * (n - 1) // 2
reversal_holds = 0
reversal_fails = 0

for bits in range(1 << m):
    T = tournament_from_bits(n, bits)

    for v in range(n):
        paths = get_ham_paths(T)
        T_flip = cut_flip_vertex(T, v)
        paths_flip = get_ham_paths(T_flip)

        pos_count = [0] * n
        pos_count_flip = [0] * n
        for p in paths:
            pos_count[p.index(v)] += 1
        for p in paths_flip:
            pos_count_flip[p.index(v)] += 1

        ok = all(pos_count_flip[k] == pos_count[n-1-k] for k in range(n))
        if ok:
            reversal_holds += 1
        else:
            reversal_fails += 1

print(f"  Position-reversal: {reversal_holds} hold, {reversal_fails} fail")
if reversal_fails == 0:
    print("  THEOREM: h_k(phi_v(T), v) = h_{n-1-k}(T, v) for all T, v, k")
    print("  COROLLARY: delta_H(phi_v) = H(phi_v(T)) - H(T)")
    print("           = sum_k h_{n-1-k}(T,v) - sum_k h_k(T,v) = 0")
    print("  Wait — this means delta_H = 0 always? That contradicts the data!")
    print("  Let me recheck...")
elif reversal_fails > 0:
    print("  Position-reversal is FALSE.")
    print("  Need a different structural decomposition.")

# Actually, let me just verify delta_H values more carefully
print("\n\n--- Sanity check: delta_H values ---")
n = 4
T = tournament_from_bits(n, 4)
H = hamiltonian_path_count(T)
print(f"  T from bits=4, H={H}")
for v in range(n):
    T_flip = cut_flip_vertex(T, v)
    H_flip = hamiltonian_path_count(T_flip)
    print(f"  phi_{v}: H={H} -> {H_flip}, delta={H_flip - H}")

# H(T) by position of v
paths = get_ham_paths(T)
for v in range(n):
    counts = [0] * n
    for p in paths:
        counts[p.index(v)] += 1
    print(f"  v={v}: position counts = {counts}")

    T_flip = cut_flip_vertex(T, v)
    paths_flip = get_ham_paths(T_flip)
    counts_flip = [0] * n
    for p in paths_flip:
        counts_flip[p.index(v)] += 1
    print(f"  v={v}: flip position counts = {counts_flip}")
    print(f"  v={v}: reversed original = {counts[::-1]}")

print(f"\n{'='*70}")
print("DONE")
print("=" * 70)
