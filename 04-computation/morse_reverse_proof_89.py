#!/usr/bin/env python3
"""
PROVING THE REVERSE DIRECTION OF THM-206
opus-2026-03-14-S89

If T is not transitive, we need to find an arc whose flip decreases H.

Key idea: For a non-transitive tournament, there exist vertices i,j
where i beats j but score(i) < score(j). This is an "upset" arc.
Flipping upsets should decrease H.

Alternative: Use the Walsh derivative formula:
dH_e(T) = 2 * sum_{S: e in S} H_hat[S] * chi_S(T)
If we can show this sum is positive for some e in a 3-cycle, done.

Alternative 2: Direct counting argument via Hamiltonian path structure.
"""

from itertools import permutations
from math import factorial, comb
from collections import Counter, defaultdict

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
print("PROVING THE REVERSE: NON-TRANSITIVE => HAS DECREASING FLIP")
print("opus-2026-03-14-S89")
print("=" * 70)

# ======================================================================
# PART 1: THE "UPSET" ARC APPROACH
# ======================================================================
print("\n" + "=" * 70)
print("PART 1: DOES FLIPPING AN 'UPSET' ARC ALWAYS DECREASE H?")
print("=" * 70)

# An "upset" arc is i->j where score(i) < score(j).
# In a transitive tournament, there are no upset arcs.
# In a non-transitive tournament, at least one upset exists.

for n in range(3, 7):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    total_non_trans = 0
    upset_always_decreases = 0
    upset_sometimes_decreases = 0
    upset_never_decreases = 0

    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        scores = [sum(adj.get((i,j), 0) for j in range(n) if j != i) for i in range(n)]

        # Check if transitive
        if sorted(scores) == list(range(n)):
            continue
        total_non_trans += 1

        H_orig = compute_H(n, adj)

        # Find upset arcs and check if flipping them decreases H
        upset_decreases = []
        for e_idx, (i, j) in enumerate(edges):
            if adj[(i,j)] == 1 and scores[i] < scores[j]:
                # i beats j but i has lower score: upset
                nbr = bits ^ (1 << e_idx)
                H_flip = compute_H(n, tournament_from_bits(n, nbr))
                upset_decreases.append(H_flip < H_orig)
            elif adj[(j,i)] == 1 and scores[j] < scores[i]:
                # j beats i but j has lower score: upset
                nbr = bits ^ (1 << e_idx)
                H_flip = compute_H(n, tournament_from_bits(n, nbr))
                upset_decreases.append(H_flip < H_orig)

        if all(upset_decreases):
            upset_always_decreases += 1
        elif any(upset_decreases):
            upset_sometimes_decreases += 1
        else:
            upset_never_decreases += 1

    print(f"\n  n={n}: {total_non_trans} non-transitive tournaments")
    print(f"    Flipping ANY upset decreases H: {upset_always_decreases}")
    print(f"    Flipping SOME upset decreases H: {upset_sometimes_decreases}")
    print(f"    Flipping NO upset decreases H: {upset_never_decreases}")

# ======================================================================
# PART 2: THE "CYCLE ARC" APPROACH
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: DOES FLIPPING A 3-CYCLE ARC ALWAYS DECREASE H?")
print("=" * 70)

# If T has a 3-cycle (i,j,k) with i->j->k->i, does flipping one
# of the three arcs always decrease H?

for n in range(3, 7):
    m = n * (n - 1) // 2
    edges = []
    edge_idx = {}
    for i in range(n):
        for j in range(i+1, n):
            edge_idx[(i,j)] = len(edges)
            edge_idx[(j,i)] = len(edges)
            edges.append((i, j))

    total_non_trans = 0
    cycle_arc_always = 0
    cycle_arc_sometimes = 0
    cycle_arc_never = 0

    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        scores = [sum(adj.get((i,j), 0) for j in range(n) if j != i) for i in range(n)]

        if sorted(scores) == list(range(n)):
            continue
        total_non_trans += 1

        H_orig = compute_H(n, adj)

        # Find all 3-cycles
        cycle_arcs = set()
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                for k in range(n):
                    if k == i or k == j:
                        continue
                    if adj.get((i,j)) and adj.get((j,k)) and adj.get((k,i)):
                        # Found 3-cycle i->j->k->i
                        cycle_arcs.add(edge_idx[(i,j)])
                        cycle_arcs.add(edge_idx[(j,k)])
                        cycle_arcs.add(edge_idx[(k,i)])

        if not cycle_arcs:
            continue  # shouldn't happen for non-transitive

        # Check if flipping any cycle arc decreases H
        decreases = []
        for e_idx in cycle_arcs:
            nbr = bits ^ (1 << e_idx)
            H_flip = compute_H(n, tournament_from_bits(n, nbr))
            decreases.append(H_flip < H_orig)

        if all(decreases):
            cycle_arc_always += 1
        elif any(decreases):
            cycle_arc_sometimes += 1
        else:
            cycle_arc_never += 1

    print(f"\n  n={n}: {total_non_trans} non-transitive")
    print(f"    ALL cycle arcs decrease H: {cycle_arc_always}")
    print(f"    SOME cycle arcs decrease H: {cycle_arc_sometimes}")
    print(f"    NO cycle arc decreases H: {cycle_arc_never}")

# ======================================================================
# PART 3: THE "WORST EDGE" APPROACH
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: THE MAXIMUM-DECREASE EDGE -- WHAT CHARACTERIZES IT?")
print("=" * 70)

# For each non-transitive tournament, find the edge whose flip
# decreases H the most. What characterizes this edge?

for n in range(3, 6):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    max_decrease_properties = defaultdict(int)
    total_non_trans = 0

    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        scores = [sum(adj.get((i,j), 0) for j in range(n) if j != i) for i in range(n)]

        if sorted(scores) == list(range(n)):
            continue
        total_non_trans += 1

        H_orig = compute_H(n, adj)

        best_decrease = 0
        best_edge = None

        for e_idx in range(m):
            nbr = bits ^ (1 << e_idx)
            H_flip = compute_H(n, tournament_from_bits(n, nbr))
            decrease = H_orig - H_flip
            if decrease > best_decrease:
                best_decrease = decrease
                best_edge = e_idx

        if best_edge is not None:
            i, j = edges[best_edge]
            # Is this an upset arc?
            if adj[(i,j)] == 1 and scores[i] < scores[j]:
                max_decrease_properties['upset'] += 1
            elif adj[(j,i)] == 1 and scores[j] < scores[i]:
                max_decrease_properties['upset'] += 1
            else:
                max_decrease_properties['not_upset'] += 1

            # Is it part of a 3-cycle?
            in_cycle = False
            winner, loser = (i, j) if adj[(i,j)] == 1 else (j, i)
            for k in range(n):
                if k == i or k == j:
                    continue
                if adj.get((winner, loser)) and adj.get((loser, k)) and adj.get((k, winner)):
                    in_cycle = True
                    break
                if adj.get((k, winner)) and adj.get((winner, loser)) and adj.get((loser, k)):
                    in_cycle = True
                    break
            max_decrease_properties['in_cycle' if in_cycle else 'not_in_cycle'] += 1

    print(f"\n  n={n}: Properties of the maximum-decrease edge:")
    for prop, cnt in sorted(max_decrease_properties.items()):
        print(f"    {prop}: {cnt}/{total_non_trans} ({cnt/total_non_trans:.4f})")

# ======================================================================
# PART 4: DIRECT PROOF VIA PATH COUNTING
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: PATH COUNTING ARGUMENT")
print("=" * 70)

print("""
  Consider a non-transitive tournament T with a 3-cycle i->j->k->i.

  A Hamiltonian path P either:
  (a) Uses at most one arc of the 3-cycle (doesn't "traverse" the cycle)
  (b) Uses exactly two arcs of the 3-cycle (partially traverses)

  If P uses two arcs of the cycle, say i->j and j->k, then P contains
  the sub-path ...->i->j->k->... The third arc k->i is not used.

  When we flip k->i to i->k:
  - Paths using i->j and j->k still work (they don't use k->i)
  - The path REVERSAL that used k->i and j->k is disrupted
  - New paths using i->k become available

  But how does this CHANGE THE COUNT?

  Key insight from the n=3 case:
  C_3 has H=3: the three paths are 0->1->2, 1->2->0, 2->0->1.
  Flipping 2->0 to 0->2 gives the transitive T_3 with H=1.
  The path 0->1->2 SURVIVES.
  The paths 1->2->0 and 2->0->1 are KILLED (they used arc 2->0).
  The new path using 0->2 would be ...->0->2->..., but 0->1 and 0->2
  means 0 has two outgoing edges, and 1->2 is also present, so the
  only Hamiltonian path is 0->1->2. Net change: lost 2 paths.

  This suggests: flipping a cycle arc kills paths that use it,
  and the creation of new paths doesn't fully compensate.
""")

# Let's count more carefully for n=4
print("  Detailed path analysis for n=4, T with one 3-cycle:")
n = 4
m = 6
edges = []
edge_idx = {}
for i in range(n):
    for j in range(i+1, n):
        edge_idx[(i,j)] = len(edges)
        edge_idx[(j,i)] = len(edges)
        edges.append((i, j))

# Take a specific H=3 tournament
# From Part 4 output: T=000010 has one 3-cycle
bits = 0b000010
adj = tournament_from_bits(n, bits)
print(f"\n  Tournament T = {bits:06b}:")
for i in range(n):
    for j in range(n):
        if i != j:
            if adj[(i,j)]:
                print(f"    {i} -> {j}")

H_orig = compute_H(n, adj)
print(f"  H(T) = {H_orig}")

# List all Hamiltonian paths
print(f"  Hamiltonian paths of T:")
paths = []
for p in permutations(range(n)):
    ok = True
    for k in range(n-1):
        if adj[(p[k], p[k+1])] != 1:
            ok = False
            break
    if ok:
        paths.append(p)
        arcs_used = [(p[k], p[k+1]) for k in range(n-1)]
        print(f"    {'->'.join(map(str, p))}, arcs: {arcs_used}")

# For each edge flip, show which paths survive and which die
print(f"\n  Path analysis for each edge flip:")
for e_idx in range(m):
    i_e, j_e = edges[e_idx]
    # Which direction is the arc?
    if adj[(i_e, j_e)] == 1:
        arc = (i_e, j_e)
    else:
        arc = (j_e, i_e)
    rev_arc = (arc[1], arc[0])

    nbr = bits ^ (1 << e_idx)
    adj_flip = tournament_from_bits(n, nbr)
    H_flip = compute_H(n, adj_flip)

    # Which paths use this arc?
    paths_using = [p for p in paths if any((p[k], p[k+1]) == arc for k in range(n-1))]
    paths_not_using = [p for p in paths if all((p[k], p[k+1]) != arc for k in range(n-1))]

    # New paths in flipped tournament
    new_paths = []
    for p in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if adj_flip[(p[k], p[k+1])] != 1:
                ok = False
                break
        if ok and p not in paths:
            new_paths.append(p)

    print(f"    Flip {arc[0]}->{arc[1]} to {arc[1]}->{arc[0]}: "
          f"H {H_orig}->{H_flip}, "
          f"killed {len(paths_using)} paths, "
          f"kept {len(paths_not_using)}, "
          f"created {len(new_paths)} new")

# ======================================================================
# PART 5: THE "MAX-SCORE VERTEX IN CYCLE" STRATEGY
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: STRATEGY -- FLIP THE ARC INTO THE MAX-SCORE VERTEX")
print("=" * 70)

# In a 3-cycle i->j->k->i, the vertex with the highest score
# is the one that "should" beat the others in a transitive order.
# Flipping the arc that goes INTO this vertex (making it "more dominant")
# should decrease H.

# Let's test: for each 3-cycle, flip the arc that points AWAY from
# the highest-score vertex (i.e., make the highest-score vertex win).

for n in range(3, 7):
    m = n * (n - 1) // 2
    edges_list = []
    e_idx_map = {}
    for i in range(n):
        for j in range(i+1, n):
            e_idx_map[(i,j)] = len(edges_list)
            e_idx_map[(j,i)] = len(edges_list)
            edges_list.append((i, j))

    strategy_works = 0
    strategy_fails = 0
    total_non_trans = 0

    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        scores = [sum(adj.get((i,j), 0) for j in range(n) if j != i) for i in range(n)]

        if sorted(scores) == list(range(n)):
            continue
        total_non_trans += 1

        H_orig = compute_H(n, adj)

        # Find a 3-cycle
        found_decrease = False
        for i in range(n):
            if found_decrease:
                break
            for j in range(n):
                if j == i or found_decrease:
                    continue
                for k in range(n):
                    if k == i or k == j or found_decrease:
                        continue
                    if adj.get((i,j)) and adj.get((j,k)) and adj.get((k,i)):
                        # 3-cycle i->j->k->i
                        # Find vertex with highest score in {i,j,k}
                        triple = [(scores[i], i), (scores[j], j), (scores[k], k)]
                        triple.sort(reverse=True)
                        max_v = triple[0][1]

                        # The arc that max_v LOSES in the cycle
                        # If max_v = i: k beats i (arc k->i). Flip to i->k.
                        # If max_v = j: i beats j... wait, j loses to...
                        # In cycle i->j->k->i, j beats k, and i beats j, k beats i.
                        # j LOSES to i. So flip i->j to j->i? That makes j beat i.

                        # The losing arc for max_v:
                        if max_v == i:
                            # i loses to k (k->i). Flip k->i to i->k.
                            flip_edge = e_idx_map[(i, k)]
                        elif max_v == j:
                            # j loses to i (i->j). Flip i->j to j->i.
                            flip_edge = e_idx_map[(i, j)]
                        else:  # max_v == k
                            # k loses to j (j->k). Flip j->k to k->j.
                            flip_edge = e_idx_map[(j, k)]

                        nbr = bits ^ (1 << flip_edge)
                        H_flip = compute_H(n, tournament_from_bits(n, nbr))
                        if H_flip < H_orig:
                            found_decrease = True
                            strategy_works += 1
                        else:
                            strategy_fails += 1
                        break

        if not found_decrease and total_non_trans > strategy_works:
            # Try ALL 3-cycles with this strategy
            pass  # already counted as fail

    print(f"\n  n={n}: 'Flip max-score losing arc' strategy:")
    print(f"    Works (first cycle tried): {strategy_works}/{total_non_trans}")
    print(f"    Fails: {strategy_fails}")
    success_rate = strategy_works / total_non_trans if total_non_trans > 0 else 0
    print(f"    Success rate: {success_rate:.4f}")

# ======================================================================
# PART 6: THE "ANY CYCLE ARC" STRATEGY WITH SCORE ORDERING
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: ALWAYS FIND A CYCLE ARC THAT DECREASES H?")
print("=" * 70)

# For each non-transitive tournament, check ALL 3-cycles.
# For each 3-cycle, does AT LEAST ONE arc flip decrease H?

for n in range(3, 7):
    m = n * (n - 1) // 2
    edges_list = []
    e_idx_map = {}
    for i in range(n):
        for j in range(i+1, n):
            e_idx_map[(i,j)] = len(edges_list)
            e_idx_map[(j,i)] = len(edges_list)
            edges_list.append((i, j))

    total_non_trans = 0
    all_cycles_have_decrease = 0
    some_cycles_miss = 0

    for bits in range(2**m):
        adj = tournament_from_bits(n, bits)
        scores = [sum(adj.get((i,j), 0) for j in range(n) if j != i) for i in range(n)]

        if sorted(scores) == list(range(n)):
            continue
        total_non_trans += 1

        H_orig = compute_H(n, adj)

        # Find all 3-cycles
        cycles = []
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    for a, b, c in [(i,j,k), (i,k,j), (j,i,k), (j,k,i), (k,i,j), (k,j,i)]:
                        if adj.get((a,b)) and adj.get((b,c)) and adj.get((c,a)):
                            cycles.append((a,b,c))
                            break  # one orientation per triple

        # For each cycle, check if any arc flip decreases H
        all_have = True
        for (a, b, c) in cycles:
            arcs = [(a,b), (b,c), (c,a)]
            any_decrease = False
            for arc in arcs:
                e = e_idx_map[arc]
                nbr = bits ^ (1 << e)
                H_flip = compute_H(n, tournament_from_bits(n, nbr))
                if H_flip < H_orig:
                    any_decrease = True
                    break
            if not any_decrease:
                all_have = False

        if all_have:
            all_cycles_have_decrease += 1
        else:
            some_cycles_miss += 1

    print(f"\n  n={n}: {total_non_trans} non-transitive")
    print(f"    EVERY 3-cycle has a decreasing arc: {all_cycles_have_decrease}")
    print(f"    SOME 3-cycle has no decreasing arc: {some_cycles_miss}")

print("\n" + "=" * 70)
print("DONE -- REVERSE DIRECTION ANALYSIS")
print("=" * 70)
