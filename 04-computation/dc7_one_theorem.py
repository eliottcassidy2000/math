"""
dc7_one_theorem.py -- kind-pasteur-2026-03-13-S61

WHY is |dc7| <= 1 at n=7 under Vitali reversal?

This is a remarkable constraint. We know:
1. Every 7-cycle uses at least 1 S-S arc (pigeonhole)
2. c7(A) and c7(B) are disjoint (no cycle in both)
3. dc7 = c7(B) - c7(A) in {-1, 0, 1}

This means: c7(A) and c7(B) differ by at most 1!
c7(B) = c7(A) + dc7 where dc7 in {-1, 0, 1}.

Hypothesis: the 7-cycles can be PAIRED between A and B,
with at most 1 unpaired cycle. Each A-cycle maps to a B-cycle
via a "conjugation" that swaps the S-arcs' directions.

This is like a MATCHING in a bipartite graph:
  Left = {7-cycles in A}, Right = {7-cycles in B}
  Edge = cycle in A can be "transformed" to cycle in B by
         switching the direction of S-S arcs.

If the matching is almost perfect (size = min(|L|, |R|) - at most 1),
then |dc7| <= 1.

This script explores:
1. The explicit matching between A-cycles and B-cycles
2. WHY the matching is almost perfect
3. What determines whether dc7 = +1, 0, or -1
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def sub_scores(A, n, subset):
    k = len(subset)
    return tuple(sorted([sum(A[subset[i]][subset[j]] for j in range(k) if i != j) for i in range(k)]))

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print("WHY |dc7| <= 1 AT n=7")
print("=" * 60)

np.random.seed(42)
examples = []

for trial in range(5000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    L = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        if not np.array_equal(L, lambda_graph(B, n)):
            continue

        S = set(subset)
        S_sorted = sorted(S)
        outside = sorted(set(range(n)) - S)

        # Find ALL directed Hamiltonian cycles in A and B
        cycles_A = set()
        cycles_B = set()

        for perm in permutations(range(1, n)):
            cycle = (0,) + perm
            valid_A = all(A[cycle[i]][cycle[(i+1)%n]] for i in range(n))
            valid_B = all(B[cycle[i]][cycle[(i+1)%n]] for i in range(n))

            # Canonicalize: rotate to start with smallest vertex
            # and compare with reverse
            if valid_A:
                min_start = min(range(n), key=lambda i: cycle[i])
                canonical = tuple(cycle[(min_start + j) % n] for j in range(n))
                cycles_A.add(canonical)

            if valid_B:
                min_start = min(range(n), key=lambda i: cycle[i])
                canonical = tuple(cycle[(min_start + j) % n] for j in range(n))
                cycles_B.add(canonical)

        c7_A = len(cycles_A) // 7  # each cycle appears 7 times
        c7_B = len(cycles_B) // 7
        dc7 = c7_B - c7_A

        # For each A-cycle, try to find a matching B-cycle
        # A matching: swap the S-S arc directions in the cycle
        # i.e., at each point where the cycle traverses an S-S arc,
        # reverse the direction.

        # The "S-arc conjugation" of a cycle:
        # For each arc (cycle[i], cycle[i+1]) where both are in S,
        # the reverse tournament has the opposite arc.
        # So the cycle in A can be converted to a potential cycle in B
        # by reversing all S-S arcs' endpoints.

        # But this doesn't create a valid cycle in general!
        # The cycle visits vertices in order, and reversing some arcs
        # breaks the path.

        # Let me think about this differently.
        # A 7-cycle is a cyclic permutation pi of {0,...,6}.
        # It's valid in A iff A[pi(i)][pi(i+1)] = 1 for all i.
        # In B: B[pi(i)][pi(i+1)] = 1 for all i.
        # B differs from A only on S-S entries.
        # So the cycle works in B iff:
        #   For S-X and X-X arcs: same as A (unchanged)
        #   For S-S arcs: B[s1][s2] = A[s2][s1] (reversed)
        # So the cycle works in B iff:
        #   All non-S-S arcs: pi(i)->pi(i+1) valid in A
        #   All S-S arcs: pi(i+1)->pi(i) valid in A (reversed)

        # This means: for a cycle to work in B, at every S-S step,
        # the cycle must traverse AGAINST the A-direction.

        # Now: consider the "conjugate" cycle:
        # For a cycle in A, reverse the order of each S-segment.
        # But that's not well-defined as a simple transformation.

        # Instead: the number of S-S arcs in each cycle is k >= 1.
        # The S-vertices divide the cycle into "S-segments" (maximal runs of S-vertices)
        # and "X-segments" (runs of X-vertices between S-segments).

        # For each S-segment of length m, there are m-1 S-S arcs.
        # Under reversal, these arcs flip. The segment s1->s2->...->sm
        # in A becomes sm->...->s2->s1 in B (direction reversed).

        # For the FULL cycle to work in B:
        # Each S-segment must run in reverse direction.
        # But the connections between segments (S->X and X->S) are unchanged.

        # So the B-cycle, if it exists, reverses each S-segment while keeping
        # the rest of the cycle intact.

        # At n=7 with 4 S-vertices and 3 X-vertices:
        # The possible S-segment structures are:
        # (4): all 4 S-vertices consecutive => 3 S-S arcs, 1 S-segment
        # (3,1): 3 consecutive + 1 isolated => 2 S-S arcs, 2 S-segments
        # (2,2): two pairs => 2 S-S arcs, 2 S-segments
        # (2,1,1): one pair + two isolated => 1 S-S arc, 3 S-segments
        # (1,1,1,1): impossible at n=7 (pigeonhole)

        # For each A-cycle, the "S-reversed" cycle:
        # - Reverse each S-segment
        # - Keep X-segments unchanged
        # This gives a SPECIFIC permutation. Is it a valid cycle in B?

        # Let me implement this.
        matched = 0
        unmatched_A = []
        unmatched_B = []

        # Get unique directed cycles (normalize by starting from vertex 0's position)
        unique_A = set()
        unique_B = set()
        for perm in permutations(range(1, n)):
            cycle = (0,) + perm
            if all(A[cycle[i]][cycle[(i+1)%n]] for i in range(n)):
                unique_A.add(cycle)
            if all(B[cycle[i]][cycle[(i+1)%n]] for i in range(n)):
                unique_B.add(cycle)

        # Build the S-reversal mapping
        # For a cycle c = (c0, c1, ..., c6), find S-segments and reverse them.
        def s_reverse(cycle, S_set):
            """Reverse S-segments in a cycle."""
            n = len(cycle)
            result = list(cycle)
            # Find S-segments (maximal runs of S-vertices)
            i = 0
            while i < n:
                if cycle[i] in S_set:
                    # Start of S-segment
                    j = i
                    while j < n and cycle[j % n] in S_set:
                        j += 1
                    # But handle wrap-around for cycle
                    if j <= n:
                        # Segment from i to j-1
                        segment = [cycle[k] for k in range(i, j)]
                        segment.reverse()
                        for k, v in enumerate(segment):
                            result[i + k] = v
                        i = j
                    else:
                        break
                else:
                    i += 1
            return tuple(result)

        # Actually, this is tricky because cycles wrap around.
        # Let me just check directly: for each A-cycle, is its S-reversal in B?

        # Simpler approach: for each A-cycle, does the S-reversed version exist as a B-cycle?
        # And vice versa.

        # Even simpler: just count the intersection.
        # Shared cycles should be 0 (proved).
        shared = unique_A & unique_B
        only_A = unique_A - shared
        only_B = unique_B - shared

        # Each undirected cycle appears in 7 directed forms (rotations).
        # And each directed cycle has 7 rotations.
        # But we're not de-rotating. Let me just divide by 7.

        examples.append({
            'bits': bits,
            'subset': subset,
            'dc7': dc7,
            'c7_A': len(unique_A) // 7,
            'c7_B': len(unique_B) // 7,
            'n_A': len(unique_A),
            'n_B': len(unique_B),
            'shared': len(shared),
        })
        break

    if len(examples) >= 50:
        break

print(f"Analyzed {len(examples)} Vitali pairs")

# Verify no shared cycles
print(f"\nShared cycles: {all(ex['shared'] == 0 for ex in examples)}")

# dc7 distribution with c7 counts
print(f"\ndc7 distribution:")
for dc7_val in sorted(set(ex['dc7'] for ex in examples)):
    entries = [ex for ex in examples if ex['dc7'] == dc7_val]
    c7_A_vals = [ex['c7_A'] for ex in entries]
    print(f"  dc7={dc7_val:+d}: {len(entries)} examples, "
          f"mean c7_A={np.mean(c7_A_vals):.1f}, range=[{min(c7_A_vals)},{max(c7_A_vals)}]")

# WHY |dc7| <= 1?
# Let's analyze the S-reversal structure in detail.
print(f"\n{'='*60}")
print("S-SEGMENT STRUCTURE ANALYSIS")
print(f"{'='*60}")

# For each A-cycle, classify by its S-segment structure
for ex_idx in range(min(5, len(examples))):
    ex = examples[ex_idx]
    if ex['dc7'] == 0 and ex_idx < 3:
        continue  # Skip boring dc7=0 cases for now

    bits = ex['bits']
    A = bits_to_adj(bits, n)
    B = reverse_subtournament(A, n, list(ex['subset']))
    S = set(ex['subset'])

    print(f"\nExample {ex_idx}: dc7={ex['dc7']}, c7_A={ex['c7_A']}, c7_B={ex['c7_B']}")

    # Enumerate A-cycles with S-segment structure
    cycle_structures = defaultdict(list)
    for perm in permutations(range(1, n)):
        cycle = (0,) + perm
        if all(A[cycle[i]][cycle[(i+1)%n]] for i in range(n)):
            # Find S-segment structure
            in_S = [cycle[i] in S for i in range(n)]
            # Find segments (maximal runs of True/False)
            segments = []
            i = 0
            while i < n:
                if in_S[i]:
                    j = i
                    while j < n and in_S[j % n]:
                        j += 1
                    segments.append(('S', j - i))
                    i = j
                else:
                    j = i
                    while j < n and not in_S[j % n]:
                        j += 1
                    segments.append(('X', j - i))
                    i = j

            # Normalize: rotate to start with first S-segment
            seg_key = tuple(s for s in segments if s[0] == 'S')
            s_lengths = tuple(sorted(s[1] for s in segments if s[0] == 'S'))
            cycle_structures[s_lengths].append(cycle)

    for struct, cycles in sorted(cycle_structures.items()):
        print(f"  S-segments {struct}: {len(cycles)//7} cycles (x7 rotations)")

    # Same for B-cycles
    cycle_structures_B = defaultdict(list)
    for perm in permutations(range(1, n)):
        cycle = (0,) + perm
        if all(B[cycle[i]][cycle[(i+1)%n]] for i in range(n)):
            in_S = [cycle[i] in S for i in range(n)]
            segments = []
            i = 0
            while i < n:
                if in_S[i]:
                    j = i
                    while j < n and in_S[j % n]:
                        j += 1
                    segments.append(('S', j - i))
                    i = j
                else:
                    j = i
                    while j < n and not in_S[j % n]:
                        j += 1
                    segments.append(('X', j - i))
                    i = j
            s_lengths = tuple(sorted(s[1] for s in segments if s[0] == 'S'))
            cycle_structures_B[s_lengths].append(cycle)

    print(f"  B-cycles:")
    for struct, cycles in sorted(cycle_structures_B.items()):
        print(f"    S-segments {struct}: {len(cycles)//7} cycles")

    # Compare: same structure type in A and B?
    for struct in sorted(set(cycle_structures.keys()) | set(cycle_structures_B.keys())):
        a_count = len(cycle_structures.get(struct, [])) // 7
        b_count = len(cycle_structures_B.get(struct, [])) // 7
        diff = b_count - a_count
        if a_count + b_count > 0:
            print(f"    {struct}: A={a_count}, B={b_count}, diff={diff:+d}")

# Key insight: the dc7 should come from the (4,) segment type
# (all 4 S-vertices consecutive). Under reversal, reversing a 4-segment
# creates a different ordering, and there's either 0 or 1 extra cycle.

print(f"\n{'='*60}")
print("THE (4,) SEGMENT: ALL S CONSECUTIVE")
print(f"{'='*60}")

# When all 4 S-vertices are consecutive: ...S1 S2 S3 S4 X1 X2 X3 ...
# Under reversal: ...S4 S3 S2 S1 X1 X2 X3 ... (S-segment reversed)
# The reversal creates a SPECIFIC new cycle.
# This cycle exists in B iff:
# S4->S3->S2->S1 in B (which is S1->S2->S3->S4 in A, always true)
# AND the X-segments are unchanged.
# So the reversed cycle ALWAYS exists in B!

# Wait, but the X-to-S connections:
# Original: ...Xprev -> S1 and S4 -> Xnext
# Reversed: ...Xprev -> S4 and S1 -> Xnext
# But Xprev->S4 is unchanged (X-S arc), and S1->Xnext is unchanged (S-X arc).
# So we need A[Xprev][S4] = 1 and A[S1][Xnext] = 1.
# In the original: A[Xprev][S1] = 1 and A[S4][Xnext] = 1.
# But in the reversed: we need A[Xprev][S4] and A[S1][Xnext].
# These are DIFFERENT arcs! They may or may not be 1.

# So the reversed (4,) cycle does NOT always exist in B.
# This is where the +/-1 comes from!

# Let me verify: for a (4,) cycle S1->S2->S3->S4->X1->X2->X3->S1
# The reversed version is: S4->S3->S2->S1->X1->X2->X3->S4
# Changed connections: Xprev=X3, Xnext=X1
# Original: X3->S1, S4->X1
# Reversed: X3->S4, S1->X1
# These are DIFFERENT! X3->S4 and S1->X1 must hold.

# For (2,2) and (3,1) segments: similar analysis.
# The PAIRING between A-cycles and B-cycles may leave exactly 0 or 1 unmatched.

# Let me check explicitly: do (2,2) and (2,1,1) cycles ALWAYS pair up?
print("\nDetailed pairing analysis:")
for ex_idx in range(min(20, len(examples))):
    ex = examples[ex_idx]
    if ex['c7_A'] == 0 and ex['c7_B'] == 0:
        continue

    bits = ex['bits']
    A = bits_to_adj(bits, n)
    B = reverse_subtournament(A, n, list(ex['subset']))
    S = set(ex['subset'])

    # Classify all cycles by structure AND by the S-ordering within the cycle
    for perm in permutations(range(1, n)):
        cycle = (0,) + perm
        valid_A = all(A[cycle[i]][cycle[(i+1)%n]] for i in range(n))
        valid_B = all(B[cycle[i]][cycle[(i+1)%n]] for i in range(n))

        if not (valid_A or valid_B):
            continue

        # Get S-vertices in cycle order
        s_order = [v for v in cycle if v in S]
        # Get the S-segment lengths
        in_S = [cycle[i] in S for i in range(n)]
        s_lengths = []
        i = 0
        while i < n:
            if in_S[i]:
                j = i
                while j < n and in_S[j]:
                    j += 1
                s_lengths.append(j - i)
                i = j
            else:
                i += 1
        s_struct = tuple(sorted(s_lengths))

        # For the first few non-trivial examples, show the structure
        if valid_A and not valid_B:
            pass  # A-only cycle
        elif valid_B and not valid_A:
            pass  # B-only cycle

    # Just show the dc7 and structure counts
    break

# Let me try a different approach: direct proof attempt
print(f"\n{'='*60}")
print("PROOF SKETCH: WHY |dc7| <= 1 AT n=7")
print(f"{'='*60}")

# At n=7 with |S|=4, |X|=3:
# A Hamiltonian cycle visits all 7 vertices.
# The S-vertices appear in some positions in the cycle.
# Between consecutive S-vertices in the cycle, there are 0 or more X-vertices.
# 4 gaps between S-vertices sum to 7, each >= 1.
# Possible gap sequences (unordered): (1,1,1,4), (1,1,2,3), (1,2,2,2)
# (Higher gap = more X-vertices between S-vertices)
#
# S-S arcs = #{gaps of size 1}:
# (1,1,1,4): 3 S-S arcs (but 4 gaps with one having 4 X-vertices)
# Wait: 4 gaps summing to 7. The gaps represent the TOTAL number of
# positions between consecutive S-vertices (including the S-vertices themselves?
# No: gaps = number of X-vertices between consecutive S-vertices + 1.
# Actually, let me reconsider.
#
# In a cycle of 7 vertices, 4 are S-vertices.
# They divide the cycle into 4 arcs (between consecutive S-vertices going around).
# Each arc has at least 0 X-vertices.
# Total X-vertices = 3.
# So the 4 arc lengths (number of edges) sum to 7.
# An arc of length 1 = direct S-S edge (no X-vertex between).
# An arc of length 2 = one X-vertex between two S-vertices.
# Etc.
#
# 4 arc lengths >= 1 summing to 7.
# Compositions of 7 into 4 parts, each >= 1:
# After substituting a_i = b_i + 1: 4 values b_i >= 0 summing to 3.
# Solutions: C(6,3) = 20 ordered compositions.
# Unordered partitions: 3 = 3+0+0+0, 2+1+0+0, 1+1+1+0.
# So arc lengths: (4,1,1,1), (3,2,1,1), (2,2,2,1).
# Number of S-S arcs = #{arcs of length 1} = 3, 2, 1 respectively.
#
# Under reversal, ALL S-internal arcs flip.
# This means: for each cycle, the S-ordering reverses within each S-segment.
#
# The KEY observation: the (1,1,2,2) tournament on S has EXACTLY ONE
# labeling up to the automorphism group (|Aut|=1).
# So there are 24 distinct labelings, and the reversal maps between them.
#
# For the SPECIFIC S-ordering within a cycle segment:
# If the segment has k S-vertices, there are k! possible orderings.
# Under reversal, the ordering reverses: (s1,s2,...,sk) -> (sk,...,s2,s1).
# For this to be valid in B: we need B[sk][s_{k-1}] = ... = 1.
# Which means A[s_{k-1}][sk] = ... = 1 (opposite direction).
# But in A, the segment was s1->s2->...->sk.
# In B, we need sk->s_{k-1}->...->s1.
# This is exactly the REVERSE of the A-segment, which always works in B
# (by definition of reversal: B[si][sj] = A[sj][si]).
#
# BUT: the connections from X to S at the segment boundaries change!
# If the A-cycle has ...Xprev -> S1 -> S2 -> ... -> Sk -> Xnext ...
# Then the B-"conjugate" has ...Xprev -> Sk -> ... -> S2 -> S1 -> Xnext ...
# We need: A[Xprev][Sk] = 1 and A[S1][Xnext] = 1.
# In the original: A[Xprev][S1] = 1 and A[Sk][Xnext] = 1.
# So the B-conjugate exists iff BOTH:
#   (a) A[Xprev][Sk] = 1 (X-vertex connects to LAST S in segment)
#   (b) A[S1][Xnext] = 1 (FIRST S connects to next X-vertex)

print("""
The conjugation maps each A-cycle to a potential B-cycle by reversing S-segments.
The conjugated cycle exists in B iff the boundary connections work:
  A[Xprev][S_last] = 1 and A[S_first][Xnext] = 1.

At n=7, each cycle has at most 3 S-segments (one per gap of size 1).
The conjugation can fail at 2-4 boundary points.

CLAIM: the number of "boundary failures" is always even,
so the parity of |A-cycles matching to B-cycles| is preserved,
giving |dc7| <= 1 (at most one unmatched cycle).
""")

# Verify: compute the conjugation success rate
print("Conjugation analysis:")
match_stats = Counter()

for ex in examples[:30]:
    bits = ex['bits']
    A = bits_to_adj(bits, n)
    B = reverse_subtournament(A, n, list(ex['subset']))
    S = set(ex['subset'])

    for perm in permutations(range(1, n)):
        cycle = (0,) + perm
        if not all(A[cycle[i]][cycle[(i+1)%n]] for i in range(n)):
            continue

        # Find S-segments in the cycle
        # A segment is a maximal run of S-vertices in the cycle
        in_S = [cycle[i] in S for i in range(n)]

        # For each S-segment, identify boundary vertices
        # and check if conjugation works
        # Build the conjugated cycle
        conj = list(cycle)
        i = 0
        segments = []
        while i < n:
            if in_S[i]:
                j = i
                while j < n and in_S[j]:
                    j += 1
                seg = list(cycle[i:j])
                seg.reverse()
                for k, v in enumerate(seg):
                    conj[i + k] = v
                segments.append((i, j))
                i = j
            else:
                i += 1

        # Check if conjugated cycle is valid in B
        conj = tuple(conj)
        valid_in_B = all(B[conj[i]][conj[(i+1)%n]] for i in range(n))

        # Count boundary failures
        # For each S-segment [i, j), check:
        # Xprev -> S_last: A[cycle[(i-1)%n]][conj[i]] (= A[Xprev][S_last])
        # S_first -> Xnext: A[conj[j-1]][cycle[j%n]] (= A[S_first][Xnext])
        failures = 0
        for seg_start, seg_end in segments:
            xprev = cycle[(seg_start - 1) % n]
            s_first_new = conj[seg_start]
            s_last_new = conj[seg_end - 1]
            xnext = cycle[seg_end % n]

            if not in_S[(seg_start - 1) % n]:
                if not A[xprev][s_first_new]:
                    failures += 1
            if seg_end < n and not in_S[seg_end % n]:
                if not A[s_last_new][xnext]:
                    failures += 1

        match_stats[(valid_in_B, failures)] += 1

print(f"  (valid_in_B, boundary_failures): {dict(sorted(match_stats.items()))}")

print("\nDone.")
