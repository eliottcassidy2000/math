#!/usr/bin/env python3
"""
UPSET MONOTONICITY: DEEP STRUCTURE AND PROOF STRATEGIES
opus-2026-03-14-S89

The Upset Monotonicity Conjecture: If i->j is an arc with
score(i) < score(j), then H(T) > H(T^{ij}).

This script explores:
1. The quantitative relationship: how much does H decrease?
2. Connection to deletion-contraction
3. The "path through arc" counting
4. Connection to the Harary-Moser formula
5. Potential proof via the permanent interpretation
"""

from itertools import permutations
from math import factorial, comb
from collections import Counter, defaultdict
from fractions import Fraction

def compute_H(n, adj):
    count = 0
    for p in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if adj[(p[k], p[k+1])] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count

def tournament_from_bits(n, bits):
    adj = {}
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
            idx += 1
    return adj

print("=" * 70)
print("UPSET MONOTONICITY: DEEP STRUCTURE")
print("opus-2026-03-14-S89")
print("=" * 70)

# ======================================================================
# PART 1: QUANTITATIVE dH FOR UPSET FLIPS
# ======================================================================
print("\n" + "=" * 70)
print("PART 1: HOW MUCH DOES H DECREASE WHEN FLIPPING AN UPSET ARC?")
print("=" * 70)

for n in range(3, 7):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    # Collect dH values for upset flips, categorized by score difference
    dH_by_score_diff = defaultdict(list)

    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        scores = [sum(adj.get((i,j), 0) for j in range(n) if j != i) for i in range(n)]

        H_orig = compute_H(n, adj)

        for e_idx, (ei, ej) in enumerate(edges):
            # Check if this is an upset arc
            if adj[(ei, ej)] == 1 and scores[ei] < scores[ej]:
                winner, loser = ei, ej
            elif adj[(ej, ei)] == 1 and scores[ej] < scores[ei]:
                winner, loser = ej, ei
            else:
                continue

            score_diff = scores[loser] - scores[winner]

            nbr = bits ^ (1 << e_idx)
            H_flip = compute_H(n, tournament_from_bits(n, nbr))
            dH = H_orig - H_flip  # positive means H decreased after flip
            dH_by_score_diff[score_diff].append(dH)

    print(f"\n  n={n}: dH for upset flips (grouped by score difference):")
    for sd in sorted(dH_by_score_diff.keys()):
        vals = dH_by_score_diff[sd]
        avg = sum(vals) / len(vals)
        min_v = min(vals)
        max_v = max(vals)
        all_pos = all(v > 0 for v in vals)
        print(f"    Score diff {sd}: count={len(vals)}, "
              f"avg dH={avg:.2f}, min={min_v}, max={max_v}, "
              f"all positive? {all_pos}")

# ======================================================================
# PART 2: DELETION-CONTRACTION VIEW OF dH
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: PATH COUNTING — PATHS THROUGH THE FLIPPED ARC")
print("=" * 70)

print("""
  When we flip arc i->j to j->i:
  dH = H(T) - H(T^{ij})
     = (# paths using i->j in T) - (# paths using j->i in T^{ij})

  Because paths NOT using the arc i->j (or j->i) are the same in both.

  So: dH = #{paths through i->j in T} - #{paths through j->i in T^{ij}}

  Now, #{paths through a->b} = #{(A,B) : A ends at a, B starts at b,
  A union B = V\\{a,b}, A cap B = empty}

  Let P(a, S) = # Hamiltonian paths on S ending at a (in the induced sub-T on S cup {a}).
  Let Q(b, S) = # Hamiltonian paths on S starting at b.

  Then #{paths through a->b} = sum over partitions (A,B) of V\\{a,b}
  of P(a, A) * Q(b, B).

  This is a CONVOLUTION: the product of two partition functions.
""")

# For n=4, let's compute the path-through-arc counts
n = 4
m = 6
edges = []
for i in range(n):
    for j in range(i+1, n):
        edges.append((i, j))

for bits in [0b000010, 0b000101]:  # Two H=3 tournaments
    adj = tournament_from_bits(n, bits)
    scores = [sum(adj.get((i,j), 0) for j in range(n) if j != i) for i in range(n)]
    H_orig = compute_H(n, adj)

    print(f"\n  T={bits:06b}, H={H_orig}, scores={scores}")

    for e_idx, (ei, ej) in enumerate(edges):
        # Determine arc direction
        if adj[(ei, ej)] == 1:
            a, b = ei, ej
        else:
            a, b = ej, ei

        # Count paths through a->b
        through_count = 0
        for p in permutations(range(n)):
            ok = True
            for k in range(n-1):
                if adj[(p[k], p[k+1])] != 1:
                    ok = False
                    break
            if ok:
                # Check if path uses a->b
                uses_arc = any(p[k] == a and p[k+1] == b for k in range(n-1))
                if uses_arc:
                    through_count += 1

        is_upset = (adj[(a,b)] == 1 and scores[a] < scores[b])
        upset_str = " [UPSET]" if is_upset else ""

        print(f"    Arc {a}->{b}: {through_count}/{H_orig} paths through it, "
              f"scores ({scores[a]},{scores[b]}){upset_str}")

# ======================================================================
# PART 3: THE PERMANENT INTERPRETATION
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: H AS A PERMANENT AND dH AS A MINOR DIFFERENCE")
print("=" * 70)

print("""
  H(T) = permanent of a specific matrix derived from T.

  The adjacency matrix A has A[i,j] = 1 if i->j, 0 otherwise.
  The "path matrix" P has P[i,j] = A[i,j] for consecutive positions.

  More precisely, H = sum over permutations sigma of
  prod_{k=1}^{n-1} A[sigma(k), sigma(k+1)].

  This is NOT the standard permanent (which would be
  sum prod A[i, sigma(i)]), but a "path permanent."

  The path permanent is the permanent of the (n-1) x (n-1) matrix
  M where M[k, (i,j)] = [sigma assigns (i,j) to step k].

  Actually, H = sum_{sigma in S_n} prod_{k=0}^{n-2} A[sigma(k), sigma(k+1)]

  When we flip arc i->j, we change A[i,j] from 1 to 0 and A[j,i] from 0 to 1.

  dH = H(T) - H(T^{ij})
     = sum_sigma (prod_k A[sigma(k),sigma(k+1)] - prod_k A'[sigma(k),sigma(k+1)])
  where A' is the flipped matrix.

  A path sigma contributes to dH iff it uses the arc (i,j) in T
  or the arc (j,i) in T'. Paths using neither contribute 0.

  So dH = (# paths using i->j in T) - (# paths using j->i in T')
""")

# ======================================================================
# PART 4: SCORE-BASED BOUND ON PATH COUNTS
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: SCORE-BASED BOUND ON PATH COUNTS THROUGH AN ARC")
print("=" * 70)

# For the upset monotonicity proof, we need:
# #{paths through i->j in T} > #{paths through j->i in T^{ij}}
# when score(i) < score(j).

# Intuition: after flipping i->j to j->i:
# - j has one more outgoing edge (to i), score(j) increases by 1
# - i has one fewer outgoing edge (lost j), score(i) decreases by 1
# - The new scores satisfy: score'(i) = score(i) - 1, score'(j) = score(j) + 1

# Before: score(i) < score(j), so score(i) <= score(j) - 1
# After: score'(i) = score(i) - 1 <= score(j) - 2
#         score'(j) = score(j) + 1

# The key quantity: paths through a->b depend on:
# 1. How many paths can ARRIVE at a (enter from some predecessor)
# 2. How many paths can LEAVE b (continue to some successor)
# Higher in-degree at a => more ways to arrive
# Higher out-degree at b => more ways to leave

# For i->j in T: in-degree(i) = n-1-score(i) (many arrivals since low score)
#                 out-degree(j) = score(j) (many departures since high score)

# For j->i in T': in-degree(j) = n-1-score'(j) = n-2-score(j) (fewer arrivals)
#                  out-degree(i) = score'(i) = score(i)-1 (fewer departures)

# So: paths through i->j benefit from HIGH in-degree at i and HIGH out-degree at j.
#     paths through j->i suffer from LOW in-degree at j and LOW out-degree at i.

# This is the HANDSHAKE ARGUMENT:
# The upset arc i->j is "assisted" by its score environment.
# The corrected arc j->i is "opposed" by the new scores.

for n in range(3, 7):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    # For each upset arc, compute the "score product" heuristic
    # and compare with actual dH
    score_product_correct = 0
    score_product_wrong = 0

    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        scores = [sum(adj.get((i,j), 0) for j in range(n) if j != i) for i in range(n)]
        H_orig = compute_H(n, adj)

        for e_idx, (ei, ej) in enumerate(edges):
            if adj[(ei, ej)] == 1 and scores[ei] < scores[ej]:
                winner, loser = ei, ej
            elif adj[(ej, ei)] == 1 and scores[ej] < scores[ei]:
                winner, loser = ej, ei
            else:
                continue

            # Score product heuristic for "path through arc":
            # P(i->j) ~ in_deg(i) * out_deg(j) = (n-1-s_i) * s_j
            # P(j->i) ~ in_deg(j) * out_deg(i) (after flip)
            # After flip: s'_i = s_i - 1, s'_j = s_j + 1
            # P(j->i) ~ (n-1-s'_j) * s'_i = (n-2-s_j) * (s_i-1)

            s_w = scores[winner]  # low score, beats loser
            s_l = scores[loser]   # high score, but loses

            p_original = (n - 1 - s_w) * s_l
            p_flipped = (n - 2 - s_l) * (s_w - 1) if s_w > 0 else 0

            nbr = bits ^ (1 << e_idx)
            H_flip = compute_H(n, tournament_from_bits(n, nbr))
            actual_dH = H_orig - H_flip

            if (p_original > p_flipped) == (actual_dH > 0):
                score_product_correct += 1
            else:
                score_product_wrong += 1

    total = score_product_correct + score_product_wrong
    print(f"\n  n={n}: Score-product heuristic predicts dH sign correctly:")
    if total > 0:
        print(f"    Correct: {score_product_correct}/{total} ({score_product_correct/total:.4f})")
        print(f"    Wrong: {score_product_wrong}/{total}")
    else:
        print(f"    No upset arcs at this n (all non-transitive have equal scores)")

# ======================================================================
# PART 5: THE STEEPEST DESCENT PATH — GREEDY TO TRANSITIVE
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: GREEDY DESCENT — ALWAYS FLIP BIGGEST UPSET")
print("=" * 70)

# Since every upset flip decreases H, we can greedily flip upsets
# to reach a transitive tournament. How many steps does it take?

for n in range(3, 7):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    step_counts = Counter()
    max_steps = 0

    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        scores = [sum(adj.get((i,j), 0) for j in range(n) if j != i) for i in range(n)]

        steps = 0
        current_bits = bits

        while True:
            current_adj = tournament_from_bits(n, current_bits)
            current_scores = [sum(current_adj.get((i,j), 0) for j in range(n) if j != i) for i in range(n)]

            if sorted(current_scores) == list(range(n)):
                break  # transitive

            # Find the biggest upset (largest score difference)
            best_upset = None
            best_diff = 0

            for e_idx, (ei, ej) in enumerate(edges):
                if current_adj[(ei, ej)] == 1 and current_scores[ei] < current_scores[ej]:
                    diff = current_scores[ej] - current_scores[ei]
                    if diff > best_diff:
                        best_diff = diff
                        best_upset = e_idx
                elif current_adj[(ej, ei)] == 1 and current_scores[ej] < current_scores[ei]:
                    diff = current_scores[ei] - current_scores[ej]
                    if diff > best_diff:
                        best_diff = diff
                        best_upset = e_idx

            if best_upset is None:
                break  # shouldn't happen

            current_bits ^= (1 << best_upset)
            steps += 1

            if steps > m * 2:
                print(f"  WARNING: too many steps for {bits:0{m}b}")
                break

        step_counts[steps] += 1
        max_steps = max(max_steps, steps)

    print(f"\n  n={n}: Steps to reach transitive via greedy upset-flip:")
    print(f"    Max steps: {max_steps}")
    for s in range(max_steps + 1):
        if step_counts[s] > 0:
            print(f"    {s} steps: {step_counts[s]} tournaments ({step_counts[s]/2**m:.4f})")

# ======================================================================
# PART 6: THE BUBBLESORT CONNECTION
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: THE BUBBLESORT CONNECTION")
print("=" * 70)

print("""
  INSIGHT: Flipping upset arcs in a tournament is EXACTLY like
  performing bubble sort on the score sequence!

  In bubble sort: if element i > element j and i is before j,
  swap them. Each swap reduces the number of inversions by 1.

  In tournament: if score(i) < score(j) but i beats j,
  flip the arc. Each flip makes the tournament "more transitive."

  The number of upset arcs = number of inversions in the score sequence.
  (Not exactly, because the score sequence can have ties. But for
  tournaments with distinct scores, this is exact.)

  The MAXIMUM number of steps = maximum number of inversions
  = C(n,2) for a "reverse-sorted" tournament.
  But can a tournament have ALL arcs as upsets? Only if every arc
  goes from low-score to high-score, which is... the transitive tournament
  in REVERSE order? But the reverse of a transitive tournament IS transitive
  (on the reverse ordering). So a tournament can't have all arcs as upsets.

  Actually: the maximum number of inversions in a tournament's score
  sequence is C(n,2)/2 or so (half the arcs).
""")

# Count upset arcs for each tournament
for n in range(3, 7):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    upset_counts = Counter()
    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        scores = [sum(adj.get((i,j), 0) for j in range(n) if j != i) for i in range(n)]

        upsets = 0
        for ei, ej in edges:
            if adj[(ei, ej)] == 1 and scores[ei] < scores[ej]:
                upsets += 1
            elif adj[(ej, ei)] == 1 and scores[ej] < scores[ei]:
                upsets += 1

        upset_counts[upsets] += 1

    print(f"\n  n={n}: Distribution of upset arc count:")
    for u in sorted(upset_counts.keys()):
        print(f"    {u} upsets: {upset_counts[u]} tournaments ({upset_counts[u]/2**m:.4f})")

print("\n" + "=" * 70)
print("DONE -- UPSET MONOTONICITY DEEP ANALYSIS")
print("=" * 70)
