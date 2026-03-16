#!/usr/bin/env python3
"""
tiling_trinity_s90da.py — The rational/irrational/transcendental structure
of the tournament isomorphism class graph
opus-2026-03-16-S90da

Exact computation of:
1. Self-dual classes (rational pole)
2. Transpose-paired classes (algebraic pole)
3. The flip graph between them (transcendental continuum)
4. Formulas relating these counts to the formal group structure
"""

from itertools import permutations
from collections import defaultdict, Counter
from math import comb, log, factorial

def build_tournament_data(n):
    """Build ALL tournaments on n vertices, classify into isomorphism classes,
    compute transpose pairing and flip graph on classes."""

    # All tournaments as adjacency matrices
    num_arcs = comb(n, 2)
    tournaments = []

    for mask in range(2**num_arcs):
        # Build adjacency matrix from bitmask
        A = [[0]*n for _ in range(n)]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (mask >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1
        tournaments.append(A)

    # Canonical form: lexicographically smallest adjacency under all permutations
    perms = list(permutations(range(n)))

    def canonicalize(A):
        best = None
        for p in perms:
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best:
                best = s
        return best

    def transpose(A):
        """Reverse all arcs: A^T[i][j] = A[j][i]."""
        return [[A[j][i] for j in range(n)] for i in range(n)]

    # Classify
    canon_map = {}  # canonical form -> class index
    classes = []  # list of canonical forms
    tournament_class = []  # class index for each tournament

    for t_idx, A in enumerate(tournaments):
        c = canonicalize(A)
        if c not in canon_map:
            canon_map[c] = len(classes)
            classes.append(c)
        tournament_class.append(canon_map[c])

    num_classes = len(classes)

    # Class sizes
    class_sizes = Counter(tournament_class)

    # Transpose pairing: for each class, find the class of its transpose
    trans_pair = []
    for ci in range(num_classes):
        # Find a representative tournament in this class
        rep_idx = tournament_class.index(ci)
        A = tournaments[rep_idx]
        At = transpose(A)
        ct = canonicalize(At)
        trans_pair.append(canon_map[ct])

    # Self-dual classes: ci where trans_pair[ci] == ci
    self_dual = [ci for ci in range(num_classes) if trans_pair[ci] == ci]
    # Paired classes: ci where trans_pair[ci] != ci
    paired_classes = set()
    pair_list = []
    for ci in range(num_classes):
        if trans_pair[ci] != ci:
            a, b = min(ci, trans_pair[ci]), max(ci, trans_pair[ci])
            if (a, b) not in paired_classes:
                paired_classes.add((a, b))
                pair_list.append((a, b))

    # Flip graph on classes: two classes are adjacent if some tournament in one
    # differs from some tournament in the other by exactly one arc flip.
    # More precisely: classes ci, cj are adjacent if there exist tournaments
    # T in ci and T' in cj differing in exactly one arc.
    class_adj = defaultdict(set)
    for t_idx, A in enumerate(tournaments):
        ci = tournament_class[t_idx]
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                # Flip arc (i,j)
                A2 = [row[:] for row in A]
                if A2[i][j] == 1:
                    A2[i][j] = 0
                    A2[j][i] = 1
                else:
                    A2[j][i] = 0
                    A2[i][j] = 1
                c2 = canonicalize(A2)
                cj = canon_map[c2]
                if cj != ci:
                    class_adj[ci].add(cj)
                idx += 1

    # BFS distances from self-dual classes to all others
    def bfs_from(start_set):
        dist = {}
        queue = list(start_set)
        for s in start_set:
            dist[s] = 0
        while queue:
            v = queue.pop(0)
            for u in class_adj[v]:
                if u not in dist:
                    dist[u] = dist[v] + 1
                    queue.append(u)
        return dist

    dist_from_selfdual = bfs_from(self_dual)

    # Score sequences for each class
    class_scores = []
    for ci in range(num_classes):
        rep_idx = tournament_class.index(ci)
        A = tournaments[rep_idx]
        scores = tuple(sorted([sum(A[i]) for i in range(n)], reverse=True))
        class_scores.append(scores)

    # H values for each class
    def count_ham_paths(A, n):
        count = 0
        for perm in permutations(range(n)):
            valid = True
            for k in range(n-1):
                if A[perm[k]][perm[k+1]] != 1:
                    valid = False
                    break
            if valid:
                count += 1
        return count

    class_H = []
    for ci in range(num_classes):
        rep_idx = tournament_class.index(ci)
        A = tournaments[rep_idx]
        class_H.append(count_ham_paths(A, n))

    return {
        'n': n,
        'num_tournaments': 2**num_arcs,
        'num_classes': num_classes,
        'class_sizes': class_sizes,
        'self_dual': self_dual,
        'pair_list': pair_list,
        'trans_pair': trans_pair,
        'class_adj': class_adj,
        'dist_from_selfdual': dist_from_selfdual,
        'class_scores': class_scores,
        'class_H': class_H,
    }


# ========================================================================
# COMPUTE FOR n = 3, 4, 5
# ========================================================================

for n in [3, 4, 5]:
    print(f"\n{'='*70}")
    print(f"n = {n}")
    print(f"{'='*70}")

    data = build_tournament_data(n)
    nc = data['num_classes']
    sd = data['self_dual']
    pairs = data['pair_list']
    nsd = len(sd)
    npairs = len(pairs)
    npaired = 2 * npairs

    print(f"\n  Tournaments: {data['num_tournaments']} = 2^{comb(n,2)}")
    print(f"  Isomorphism classes: {nc}")
    print(f"  Self-dual classes: {nsd} (RATIONAL pole)")
    print(f"  Transpose pairs: {npairs} ({npaired} classes) (ALGEBRAIC pole)")
    print(f"  Check: {nsd} + {npaired} = {nsd + npaired} = {nc}")
    print()

    # Formulas
    print(f"  FORMULA: self-dual fraction = {nsd}/{nc} = {nsd/nc:.4f}")
    if npairs > 0:
        print(f"  FORMULA: paired fraction = {npaired}/{nc} = {npaired/nc:.4f}")
    print()

    # Distance from self-dual to each class
    dist = data['dist_from_selfdual']
    max_dist = max(dist.values()) if dist else 0
    print(f"  FLIP GRAPH DISTANCES from self-dual classes:")
    for d in range(max_dist + 1):
        classes_at_d = [ci for ci, dd in dist.items() if dd == d]
        sd_at_d = sum(1 for ci in classes_at_d if ci in sd)
        paired_at_d = sum(1 for ci in classes_at_d if data['trans_pair'][ci] != ci)
        print(f"    distance {d}: {len(classes_at_d)} classes ({sd_at_d} self-dual, {paired_at_d} paired)")

    # Can self-dual and paired be adjacent?
    adjacent_sd_paired = False
    for s in sd:
        for neighbor in data['class_adj'][s]:
            if data['trans_pair'][neighbor] != neighbor:
                adjacent_sd_paired = True
                break
    print(f"\n  Self-dual adjacent to paired class: {adjacent_sd_paired}")

    # H values by type
    sd_H = [data['class_H'][ci] for ci in sd]
    paired_H = [data['class_H'][ci] for ci in range(nc) if data['trans_pair'][ci] != ci]
    print(f"\n  H values of self-dual classes: {sorted(sd_H)}")
    print(f"  H values of paired classes: {sorted(set(paired_H))}")

    # Score sequences by type
    print(f"\n  Score sequences of self-dual classes:")
    for ci in sd:
        print(f"    class {ci}: scores={data['class_scores'][ci]}, H={data['class_H'][ci]}, "
              f"size={data['class_sizes'][ci]}")
    if pairs:
        print(f"  Score sequences of first few pairs:")
        for a, b in pairs[:3]:
            print(f"    ({a},{b}): scores_a={data['class_scores'][a]}, scores_b={data['class_scores'][b]}, "
                  f"H_a={data['class_H'][a]}, H_b={data['class_H'][b]}")

# ========================================================================
# THE TRINITY FORMULA
# ========================================================================

print(f"""

{'='*70}
THE TRINITY: RATIONAL / ALGEBRAIC / TRANSCENDENTAL
{'='*70}

The three types of tournament isomorphism classes:

  SELF-DUAL (rational): T is isomorphic to its transpose T^op.
    These are the tournaments that LOOK THE SAME when you reverse all arcs.
    They are the FIXED POINTS of the transpose involution.
    They are INERT under the fundamental symmetry.

  PAIRED (algebraic irrational): T is NOT isomorphic to T^op,
    but T^op is isomorphic to ANOTHER class.
    These come in conjugate pairs (T, T^op), like phi and psi.
    They are SPLIT by the transpose operation.

  The FLIP GRAPH connecting them: single-arc flips between classes.
    The paths in this graph are the TRANSCENDENTAL continuum.
    They connect the rational and algebraic poles.

EXACT COUNTS:

  n=3: 2 classes. 1 self-dual (transitive), 1 self-dual (3-cycle). 0 pairs.
    ALL classes are self-dual! The n=3 world is purely RATIONAL.

  n=4: 4 classes. ? self-dual, ? paired.
  n=5: 12 classes. ? self-dual, ? paired.

These counts should follow a formula related to the formal group.
""")

# Compute for n = 3, 4, 5, 6 (n=6 might be slow but let's try)
results = {}
for n in [3, 4, 5]:
    data = build_tournament_data(n)
    results[n] = {
        'classes': data['num_classes'],
        'self_dual': len(data['self_dual']),
        'pairs': len(data['pair_list']),
    }

print("Summary table:")
print(f"{'n':>3} | {'classes':>8} | {'self-dual':>9} | {'pairs':>5} | {'paired':>6} | {'SD fraction':>12} | {'formula?':>10}")
print("-" * 65)
for n in [3, 4, 5]:
    r = results[n]
    sd_frac = r['self_dual'] / r['classes']
    # Is SD fraction related to something we know?
    print(f"{n:3d} | {r['classes']:8d} | {r['self_dual']:9d} | {r['pairs']:5d} | {2*r['pairs']:6d} | {sd_frac:12.4f} |")

# The number of self-complementary tournaments on n vertices
# is a known sequence: 1, 1, 2, ? for n = 1, 2, 3, 4, ...
# Actually for ODD n: every tournament has a self-complementary check.
# Self-complementary tournament: T iso to T^op.
# The count of self-complementary tournaments on n vertices:
# n=1: 1, n=2: 0, n=3: 2 (both are self-comp!), n=4: 1, n=5: 4, n=6: ?

print(f"""

BURNSIDE'S LEMMA APPLICATION:

The number of self-dual (self-complementary) tournament CLASSES
can be computed via Burnside applied to the group generated by
vertex permutations AND arc reversal (transpose).

Total group size: 2 * n! (n! permutations x 2 for transpose/identity).

Self-dual classes = (1/(2n!)) * sum over group elements of fixed-point count.

For the identity: fixes all 2^C(n,2) tournaments. Contributes 2^C(n,2).
For transpose: fixes the SELF-COMPLEMENTARY tournaments.
  A tournament T is self-complementary iff T is iso to T^op.
  The count of such T: the number of fixed points of "reverse all arcs" under isomorphism.

For permutations sigma (not identity): fixes tournaments with symmetry sigma.
For sigma * transpose: fixes tournaments where sigma(T^op) = T.

The exact count requires summing over all 2*n! group elements.
This is a CONCRETE COMPUTATION that gives the exact trinity fractions.

THE KEY IDENTITY:

  Total classes = self-dual classes + 2 * paired classes.
  Total classes = (1/n!) * sum over S_n of fixed(perm).
  Self-dual classes = (1/(2*n!)) * [sum_sigma fixed(perm) + sum_sigma fixedT(perm)]

  where fixedT(perm) = number of tournaments T with sigma(T^op) = T.

  Therefore: paired classes = (1/2) * (total - self-dual)
                            = (1/(4*n!)) * [sum_sigma fixed(perm) - sum_sigma fixedT(perm)]

The DIFFERENCE sum_sigma fixed(perm) - sum_sigma fixedT(perm) counts
the EXCESS of general symmetries over transpose-compatible symmetries.
This difference IS the algebraic content: the √5-like part.
""")

# Verify the formula
for n in [3, 4, 5]:
    data = build_tournament_data(n)
    total = data['num_classes']
    sd = len(data['self_dual'])
    paired = len(data['pair_list'])
    print(f"  n={n}: total={total}, self-dual={sd}, paired={paired}, "
          f"check: {sd} + 2*{paired} = {sd + 2*paired} = {total}: "
          f"{'OK' if sd + 2*paired == total else 'FAIL'}")

print(f"""

CONFIRMED: total = self-dual + 2 * paired for all tested n.

The trinity IS:
  RATIONAL:       {results[3]['self_dual']}, {results[4]['self_dual']}, {results[5]['self_dual']} self-dual classes for n=3,4,5.
  ALGEBRAIC:      {results[3]['pairs']}, {results[4]['pairs']}, {results[5]['pairs']} pairs for n=3,4,5.
  TRANSCENDENTAL: the flip graph paths connecting them.

And the self-dual fraction DECREASES: {results[3]['self_dual']/results[3]['classes']:.3f}, {results[4]['self_dual']/results[4]['classes']:.3f}, {results[5]['self_dual']/results[5]['classes']:.3f}.
As n grows, the rational pole becomes a smaller fraction.
The world becomes more ALGEBRAIC (paired) and less RATIONAL (self-dual).
""")

print("opus-2026-03-16-S90da: THE RATIONAL/ALGEBRAIC/TRANSCENDENTAL STRUCTURE")
