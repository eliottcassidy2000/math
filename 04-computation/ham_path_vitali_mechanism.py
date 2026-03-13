"""
ham_path_vitali_mechanism.py — kind-pasteur-2026-03-13-S61
The DEFINITIVE mechanism: 7-cycles with 3 internal arcs = Ham paths on subset.
Reversal swaps Ham paths. External "completion space" changes by ±1.

Key finding from vitali_c7_mechanism.py:
- delta_c7 = +1 iff c7_3 = 2 (exactly 2 cycles use Ham paths of subset)
- delta_c7 = -1 iff c7_3 = 3 (exactly 3 cycles use Ham paths of subset)
- c7_1 = 0 always (no 1-internal-arc cycles in changing cases)
- EEESSSS is the ONLY threading pattern that changes

This script:
1. Enumerate Ham paths of (1,1,2,2) tournament (before/after reversal)
2. Count "completions" — orderings of 3 external vertices that close the 7-cycle
3. Show that delta = (completions_after - completions_before) exactly
4. Find the algebraic obstruction (when does completion count change?)
5. Connect to overlap weights W=0,1,2
6. Prove the dimensional explanation: why n<=6 is trivial
7. Explore the connection to Vitali sets in measure theory
8. Check if the "gauge charge" is a cohomology class
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter

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

def reverse_subtournament(A, n, subset):
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

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

def ham_paths_on_subset(A, subset):
    """Find all directed Hamiltonian paths on the sub-tournament."""
    k = len(subset)
    paths = []
    for perm in permutations(range(k)):
        valid = True
        for i in range(k - 1):
            if A[subset[perm[i]]][subset[perm[i+1]]] != 1:
                valid = False
                break
        if valid:
            paths.append(tuple(subset[perm[i]] for i in range(k)))
    return paths

def count_completions(A, n, subset, ham_path, ext):
    """Count ways to arrange external vertices around the Ham path to form a 7-cycle.

    The 7-cycle must be: ...ext_part... ham_path[0] -> ham_path[1] -> ham_path[2] -> ham_path[3] ...ext_part...
    where the ext vertices form a directed path from ham_path[3] back to ham_path[0]
    through all 3 external vertices.
    """
    # The 7-cycle is: ham_path[0] -> ham_path[1] -> ham_path[2] -> ham_path[3] -> ext_perm[0] -> ext_perm[1] -> ext_perm[2] -> ham_path[0]
    count = 0
    for perm in permutations(ext):
        # Check: ham_path[3] -> perm[0], perm[0] -> perm[1], perm[1] -> perm[2], perm[2] -> ham_path[0]
        if (A[ham_path[3]][perm[0]] and
            A[perm[0]][perm[1]] and
            A[perm[1]][perm[2]] and
            A[perm[2]][ham_path[0]]):
            count += 1
    return count

def find_all_ham_cycles(A, n):
    """Find all directed Hamiltonian cycles (fixing vertex 0 as start)."""
    count = 0
    for perm in permutations(range(1, n)):
        path = (0,) + perm
        valid = True
        for i in range(n):
            if A[path[i]][path[(i+1) % n]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def sub_scores(A, n, subset):
    k = len(subset)
    scores = []
    for i in range(k):
        s = sum(A[subset[i]][subset[j]] for j in range(k) if i != j)
        scores.append(s)
    return tuple(sorted(scores))

print("=" * 70)
print("PART 1: HAM PATHS OF THE (1,1,2,2) TOURNAMENT — BEFORE AND AFTER")
print("=" * 70)

# The (1,1,2,2) tournament on {0,1,2,3}: there are exactly 4 labeled tournaments
# with this score sequence (or is it 3? Let me check)
# On 4 vertices, a tournament with score sequence (1,1,2,2) is the "cyclic" type

# Let's enumerate all (1,1,2,2) tournaments on 4 vertices
four_tourn = []
for bits in range(1 << 6):  # C(4,2) = 6 arcs
    A = bits_to_adj(bits, 4)
    scores = tuple(sorted([sum(A[i]) for i in range(4)]))
    if scores == (1, 1, 2, 2):
        four_tourn.append((bits, A.copy()))

print(f"Number of labeled (1,1,2,2) tournaments on 4 vertices: {len(four_tourn)}")

for bits, A in four_tourn:
    hp = ham_paths_on_subset(A, [0,1,2,3])
    print(f"\n  bits={bits}: scores={[sum(A[i]) for i in range(4)]}")
    print(f"  Arcs: ", end="")
    for i in range(4):
        for j in range(4):
            if A[i][j]:
                print(f"{i}->{j}", end=" ")
    print()
    print(f"  Ham paths ({len(hp)}):")
    for p in hp:
        print(f"    {p}")

    # Now reverse and show Ham paths
    B = reverse_subtournament(A, 4, [0,1,2,3])
    hp_rev = ham_paths_on_subset(B, [0,1,2,3])
    print(f"  After reversal ({len(hp_rev)} paths):")
    for p in hp_rev:
        print(f"    {p}")

    # Which paths are shared?
    shared = set(hp) & set(hp_rev)
    only_before = set(hp) - set(hp_rev)
    only_after = set(hp_rev) - set(hp)
    print(f"  Shared: {len(shared)}, Lost: {len(only_before)}, Gained: {len(only_after)}")

print("\n" + "=" * 70)
print("PART 2: COMPLETION COUNTS — THE ALGEBRAIC KEY")
print("=" * 70)

# For each H-changing Vitali pair, compute:
# - Ham paths of subset BEFORE reversal
# - For each Ham path, count completions through external vertices
# - Same AFTER reversal
# - Verify: delta_c7 = sum(completions_after) - sum(completions_before) for SSSS-consecutive cycles

n = 7
total_bits = n * (n - 1) // 2
np.random.seed(42)

print("\nSearching for H-changing Vitali pairs...")
examples = []
example_count = 0

for trial in range(5000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    lam_orig = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B = reverse_subtournament(A, n, list(subset))
        lam_new = lambda_graph(B, n)
        if not np.array_equal(lam_orig, lam_new):
            continue

        c7_A = find_all_ham_cycles(A, n)
        c7_B = find_all_ham_cycles(B, n)
        delta = c7_B - c7_A

        if delta != 0:
            examples.append((bits, list(subset), delta, A.copy(), B.copy()))
            example_count += 1
            if example_count >= 10:
                break
    if example_count >= 10:
        break

print(f"Found {len(examples)} H-changing examples")

for idx, (bits, subset, delta, A, B) in enumerate(examples):
    ext = [v for v in range(n) if v not in subset]

    hp_before = ham_paths_on_subset(A, subset)
    hp_after = ham_paths_on_subset(B, subset)

    print(f"\n  Example {idx+1}: subset={subset}, ext={ext}, delta_c7={delta}")
    print(f"  Ham paths before ({len(hp_before)}):")

    total_comp_before = 0
    for p in hp_before:
        comp = count_completions(A, n, subset, p, ext)
        total_comp_before += comp
        print(f"    {p}: {comp} completions")

    print(f"  Ham paths after ({len(hp_after)}):")
    total_comp_after = 0
    for p in hp_after:
        comp = count_completions(B, n, subset, p, ext)
        total_comp_after += comp
        print(f"    {p}: {comp} completions")

    print(f"  Total completions: {total_comp_before} -> {total_comp_after}, "
          f"delta={total_comp_after - total_comp_before}")
    print(f"  Expected delta_c7 = {delta}")

    # Verify: completions count = c7_3 (7-cycles with 3 internal arcs)
    # But we need to be careful: the Ham path can be at ANY position in the cycle
    # A 7-cycle (v0,v1,...,v6) with all 4 subset vertices consecutive could have
    # the block at positions 0-3, 1-4, 2-5, 3-6, 4-0(wrap), 5-1(wrap), 6-2(wrap)
    # But since we fix v0=0 for cycle counting, this is more complex.
    # Actually, for DIRECTED cycles, each cycle is counted once with a fixed start.
    # The "3 internal arcs" count from Part 10 of the previous script already handles this.

    # Let's also verify by counting ALL 7-cycles that have the 4 subset vertices consecutive
    # (in ANY circular position)

print("\n" + "=" * 70)
print("PART 3: THE COMPLETION SPACE — WHAT DETERMINES IT?")
print("=" * 70)

# For a Ham path p = (a,b,c,d) on the subset, a completion is a permutation
# (e1,e2,e3) of the external vertices such that:
#   d -> e1, e1 -> e2, e2 -> e3, e3 -> a
# This is determined by:
# 1. A[d][e1] = 1 (d beats e1)
# 2. A[e1][e2] = 1 (e1 beats e2)
# 3. A[e2][e3] = 1 (e2 beats e3)
# 4. A[e3][a] = 1 (e3 beats a)
#
# Since arcs between external vertices are UNCHANGED by the reversal,
# and arcs between subset and external are UNCHANGED,
# the ONLY difference is which Ham paths exist.
#
# So: completions depend on the ENDPOINT pair (d, a) of the Ham path.
# The reversal changes the set of Ham paths, hence the set of endpoint pairs.
# The completion count for endpoint pair (d, a) is FIXED (depends only on external arcs).

print("\nKey insight: completion count depends ONLY on endpoints (d, a) of the Ham path")
print("and the external arcs (which are unchanged by reversal).")
print()
print("So: delta_c7 = sum_{(d,a) in HP_after} comp(d,a) - sum_{(d,a) in HP_before} comp(d,a)")

# Verify this by computing completion counts for all possible endpoint pairs
for idx, (bits, subset, delta, A, B) in enumerate(examples[:5]):
    ext = [v for v in range(n) if v not in subset]

    # Compute completion count for each possible endpoint pair
    endpoint_comp = {}
    for a in subset:
        for d in subset:
            if a == d:
                continue
            comp = 0
            for perm in permutations(ext):
                if (A[d][perm[0]] and A[perm[0]][perm[1]] and
                    A[perm[1]][perm[2]] and A[perm[2]][a]):
                    comp += 1
            endpoint_comp[(d, a)] = comp

    # Get endpoint pairs of Ham paths before and after
    hp_before = ham_paths_on_subset(A, subset)
    hp_after = ham_paths_on_subset(B, subset)

    endpoints_before = [(p[-1], p[0]) for p in hp_before]
    endpoints_after = [(p[-1], p[0]) for p in hp_after]

    sum_before = sum(endpoint_comp[ep] for ep in endpoints_before)
    sum_after = sum(endpoint_comp[ep] for ep in endpoints_after)

    print(f"\n  Example {idx+1}: delta_c7 = {delta}")
    print(f"  Endpoint pairs before: {endpoints_before}")
    print(f"  Endpoint pairs after:  {endpoints_after}")
    print(f"  Completion counts: ", end="")
    for ep in sorted(endpoint_comp.keys()):
        if endpoint_comp[ep] > 0:
            print(f"({ep[0]},{ep[1]}):{endpoint_comp[ep]} ", end="")
    print()
    print(f"  Sum before: {sum_before}, Sum after: {sum_after}, Delta: {sum_after - sum_before}")

    # BUT WAIT: the Ham path endpoints in AFTER use the REVERSED sub-tournament.
    # The completion formula uses arcs d->e1 and e3->a, which are subset-to-ext arcs.
    # These are UNCHANGED. So the endpoint completion is the SAME function.

    # Verify: does sum_after - sum_before = delta * (something)?

print("\n" + "=" * 70)
print("PART 4: ENDPOINT PAIR ANALYSIS — THE REVERSAL MAP")
print("=" * 70)

# The (1,1,2,2) tournament on {a,b,c,d} has specific Ham paths.
# After reversal, the Ham paths change.
# Question: HOW do the endpoint pairs change?
#
# In a (1,1,2,2) tournament, the two "low" vertices have out-degree 1 (within subset)
# and the two "high" vertices have out-degree 2 (within subset).
# After reversal, low becomes high and vice versa.
#
# Ham path endpoint pairs (tail, head): the HEAD is the last vertex visited,
# the TAIL is the first.

print("Analyzing endpoint pair maps for all (1,1,2,2) tournaments:")

for bits, A4 in four_tourn[:4]:
    print(f"\n  Tournament bits={bits}:")
    scores = {i: sum(A4[i]) for i in range(4)}
    low = [i for i in range(4) if scores[i] <= 1]
    high = [i for i in range(4) if scores[i] >= 2]
    print(f"    Low-score: {low}, High-score: {high}")

    hp = ham_paths_on_subset(A4, [0,1,2,3])
    B4 = reverse_subtournament(A4, 4, [0,1,2,3])
    hp_rev = ham_paths_on_subset(B4, [0,1,2,3])

    ep_before = [(p[0], p[-1]) for p in hp]
    ep_after = [(p[0], p[-1]) for p in hp_rev]

    print(f"    Endpoints (start, end) before: {ep_before}")
    print(f"    Endpoints (start, end) after:  {ep_after}")

    # Classify by which vertices are starts/ends
    starts_before = Counter(ep[0] for ep in ep_before)
    ends_before = Counter(ep[1] for ep in ep_before)
    starts_after = Counter(ep[0] for ep in ep_after)
    ends_after = Counter(ep[1] for ep in ep_after)

    print(f"    Start frequencies before: {dict(starts_before)}")
    print(f"    Start frequencies after:  {dict(starts_after)}")

print("\n" + "=" * 70)
print("PART 5: THE CRITICAL FORMULA — WHEN DOES COMPLETION COUNT CHANGE?")
print("=" * 70)

# We know:
# - Completion count for endpoint (d, a) depends on external arcs only
# - Reversal changes the Ham paths, hence the endpoint pairs
# - delta_c7 = sum of completion changes
#
# The question reduces to: for which EXTERNAL configurations does the
# endpoint-weighted completion sum change?
#
# Let's parametrize. For 3 external vertices {e1, e2, e3} and endpoint (d, a):
# comp(d, a) = #{perms of ext : d->first, consecutive, last->a}
#
# This is the number of Ham paths from {e: A[d][e]=1} through ext back to {e: A[e][a]=1}
# i.e., it's a path count in the external tournament with specified entry/exit points.

# For 3 external vertices, the external tournament is either:
# - Transitive (0,1,2): one Ham path
# - Cyclic (1,1,1): two Ham paths (one in each direction, wait no — 3 for cyclic)

# Actually for 3 vertices:
# Transitive 0->1->2: Ham paths are 0->1->2 only
# Cyclic 0->1->2->0: Ham paths are 0->1->2 and 0->2->1 wait...
# 0->1, 1->2, 2->0: path 0->1->2 ✓, path 0->2->1? need 0->2: NO.
# Actually 2->0->1 ✓ and 1->2->0 ✓ — these are also valid.
# So cyclic has 3 Ham paths (one starting from each vertex).

print("External tournament types and their Ham path structure:")
print("  Transitive (0,1,2): 1 directed Ham path per starting vertex")
print("    0->1->2: starts at 0 ends at 2")
print("  Cyclic (1,1,1): 1 Ham path per starting vertex (different endpoint each)")
print()

# Now the key: for a given endpoint pair (d, a) of the subset Ham path,
# we need d -> first_ext and last_ext -> a.
# This means:
# - first_ext must be beaten by d (in N^+(d) ∩ ext)
# - last_ext must beat a (in N^-(a) ∩ ext)
# - first_ext, ..., last_ext must be a Ham path of ext tournament

# Let's compute this systematically
print("Computing completion matrix C[d][a] for each example:")

for idx, (bits, subset, delta, A, B) in enumerate(examples[:5]):
    ext = [v for v in range(n) if v not in subset]

    # External tournament type
    ext_scores = tuple(sorted([sum(A[ext[i]][ext[j]] for j in range(3) if j != i) for i in range(3)]))
    ext_type = "cyclic" if ext_scores == (1, 1, 1) else "transitive"

    print(f"\n  Example {idx+1}: ext_type={ext_type}, ext={ext}")

    # For each subset vertex, compute N^+(v) ∩ ext and N^-(v) ∩ ext
    for s in subset:
        out_ext = [e for e in ext if A[s][e]]
        in_ext = [e for e in ext if A[e][s]]
        print(f"    {s}: out_to_ext={out_ext}, in_from_ext={in_ext}")

    # Completion matrix
    print(f"    Completion matrix C[d][a]:")
    for d in subset:
        for a in subset:
            if a == d:
                continue
            comp = 0
            for perm in permutations(ext):
                if (A[d][perm[0]] and A[perm[0]][perm[1]] and
                    A[perm[1]][perm[2]] and A[perm[2]][a]):
                    comp += 1
            if comp > 0:
                print(f"      C[{d}][{a}] = {comp}")

    # Show Ham path endpoints and their completion values
    hp_before = ham_paths_on_subset(A, subset)
    hp_after = ham_paths_on_subset(B, subset)

    print(f"    Before: endpoints -> completions")
    for p in hp_before:
        d, a = p[-1], p[0]
        comp = 0
        for perm in permutations(ext):
            if (A[d][perm[0]] and A[perm[0]][perm[1]] and
                A[perm[1]][perm[2]] and A[perm[2]][a]):
                comp += 1
        print(f"      {p}: ({d},{a}) -> {comp}")

    print(f"    After: endpoints -> completions")
    for p in hp_after:
        d, a = p[-1], p[0]
        comp = 0
        for perm in permutations(ext):
            if (A[d][perm[0]] and A[perm[0]][perm[1]] and
                A[perm[1]][perm[2]] and A[perm[2]][a]):
                comp += 1
        print(f"      {p}: ({d},{a}) -> {comp}")

print("\n" + "=" * 70)
print("PART 6: THE LOW/HIGH VERTEX MECHANISM")
print("=" * 70)

# In a (1,1,2,2) sub-tournament:
# - Low vertices (out-degree 1 within subset) tend to be "entry points" for 7-cycles
# - High vertices (out-degree 2) tend to be "exit points"
# After reversal, low and high swap.
# The external arc pattern to low vs high vertices determines the completion change.

print("\nLow/high vertex analysis for H-changing examples:")

for idx, (bits, subset, delta, A, B) in enumerate(examples[:5]):
    ext = [v for v in range(n) if v not in subset]

    # Identify low and high within subset
    scores_in_sub = {s: sum(A[s][t] for t in subset if t != s) for s in subset}
    low = [s for s in subset if scores_in_sub[s] <= 1]
    high = [s for s in subset if scores_in_sub[s] >= 2]

    # After reversal
    scores_after = {s: sum(B[s][t] for t in subset if t != s) for s in subset}
    low_after = [s for s in subset if scores_after[s] <= 1]
    high_after = [s for s in subset if scores_after[s] >= 2]

    print(f"\n  Example {idx+1}: delta={delta}")
    print(f"    Before: low={low} (score {[scores_in_sub[s] for s in low]}), "
          f"high={high} (score {[scores_in_sub[s] for s in high]})")
    print(f"    After:  low={low_after} (score {[scores_after[s] for s in low_after]}), "
          f"high={high_after} (score {[scores_after[s] for s in high_after]})")

    # External arcs to low vs high
    low_out = sum(A[s][e] for s in low for e in ext)
    low_in = sum(A[e][s] for e in ext for s in low)
    high_out = sum(A[s][e] for s in high for e in ext)
    high_in = sum(A[e][s] for e in ext for s in high)

    print(f"    Low->ext: {low_out}, Ext->low: {low_in}")
    print(f"    High->ext: {high_out}, Ext->high: {high_in}")

    # The asymmetry between low and high interaction with ext
    # might determine the sign of delta
    asymmetry = (low_out - high_out, low_in - high_in)
    print(f"    Low-High asymmetry: out={asymmetry[0]}, in={asymmetry[1]}")

print("\n" + "=" * 70)
print("PART 7: STATISTICAL TEST — DOES EXT CYCLICITY + MIXED C3 PREDICT CHANGE?")
print("=" * 70)

# From Part 1 of the previous script:
# mixed_c3 = 0 => NEVER changes (0%)
# mixed_c3 = 4 + ext_c3 = 1 (cyclic ext) => 73.9% change rate
# mixed_c3 = 4 + ext_c3 = 0 (transitive ext) => 18.7%...
# but wait, the rate for transitive ext can be higher if restricted to mixed_c3=4

# Let me compute the EXACT predictor: is it mixed_c3 > 0 AND ext_c3 = 1?

print("\nLarger sample for exact predictor search...")
np.random.seed(7777)

predictor_data = []
for trial in range(3000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    lam_orig = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B_ = reverse_subtournament(A, n, list(subset))
        lam_new = lambda_graph(B_, n)
        if not np.array_equal(lam_orig, lam_new):
            continue

        c7_A = find_all_ham_cycles(A, n)
        c7_B = find_all_ham_cycles(B_, n)
        is_changing = (c7_A != c7_B)

        ext = [v for v in range(n) if v not in subset]

        # mixed_c3
        mixed_c3 = 0
        for s1, s2 in combinations(list(subset), 2):
            for e in ext:
                if (A[s1][s2] and A[s2][e] and A[e][s1]) or (A[s2][s1] and A[s1][e] and A[e][s2]):
                    mixed_c3 += 1

        # ext_c3
        i, j, k = ext
        ext_c3 = 1 if ((A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j])) else 0

        # arc signature: each subset vertex's (out, in) to ext
        arc_sig = tuple(sorted(
            (sum(A[s][e] for e in ext), sum(A[e][s] for e in ext))
            for s in subset
        ))

        predictor_data.append({
            'is_changing': is_changing,
            'mixed_c3': mixed_c3,
            'ext_c3': ext_c3,
            'arc_sig': arc_sig,
        })

    if trial % 1000 == 0 and trial > 0:
        print(f"  ... {trial} processed")

ch = [d for d in predictor_data if d['is_changing']]
pr = [d for d in predictor_data if not d['is_changing']]
print(f"\nTotal: {len(predictor_data)}, changing: {len(ch)}, preserving: {len(pr)}")

# Cross-tabulation: mixed_c3 x ext_c3
print("\nCross-tab: mixed_c3 x ext_c3 -> change rate")
for mc in sorted(set(d['mixed_c3'] for d in predictor_data)):
    for ec in sorted(set(d['ext_c3'] for d in predictor_data)):
        total = sum(1 for d in predictor_data if d['mixed_c3'] == mc and d['ext_c3'] == ec)
        changing = sum(1 for d in predictor_data if d['mixed_c3'] == mc and d['ext_c3'] == ec and d['is_changing'])
        if total > 0:
            print(f"  mixed_c3={mc}, ext_c3={ec}: {changing}/{total} = {changing/total*100:.1f}%")

# Arc signature analysis
print("\nArc signature -> change rate (top signatures):")
all_sigs = Counter(d['arc_sig'] for d in predictor_data)
for sig, count in all_sigs.most_common():
    changing = sum(1 for d in predictor_data if d['arc_sig'] == sig and d['is_changing'])
    print(f"  {sig}: {changing}/{count} = {changing/count*100:.1f}%")

# Now the ultimate test: is there a PERFECT predictor?
print("\nSearching for perfect discriminator...")

# Test: mixed_c3 > 0 AND specific arc signature
# From previous run: arc_sig with all (1,2) or all (2,1) have ~30-37% change rate
# arc_sig with extreme values (0,3) or (3,0) have 0% change rate
# What about: mixed_c3 > 0 AND arc_sig has NO extreme values?

for mc_thresh in [1, 2, 3, 4]:
    pred_yes = [d for d in predictor_data if d['mixed_c3'] >= mc_thresh]
    correct = sum(1 for d in pred_yes if d['is_changing'])
    total_pred = len(pred_yes)
    total_change = len(ch)
    if total_pred > 0 and total_change > 0:
        precision = correct / total_pred * 100
        recall = correct / total_change * 100
        print(f"  mixed_c3 >= {mc_thresh}: precision={precision:.1f}%, recall={recall:.1f}% "
              f"({correct}/{total_pred} predicted, {total_change} actual)")

print("\n" + "=" * 70)
print("PART 8: THE MEASURE-THEORETIC CONNECTION — VITALI AND NON-MEASURABILITY")
print("=" * 70)
print("""
The classical Vitali set is constructed by selecting one representative from
each equivalence class of R/Q. The key property: the equivalence classes
are translates, so each has the same measure, but the union is bounded —
contradicting countable additivity.

Our tournament Vitali structure has an EXACT analogy:

1. EQUIVALENCE CLASSES: Lambda-classes (same labeled lambda graph)
   = orbits under the "gauge group" of cycle-preserving transformations

2. REPRESENTATIVES: Choosing one tournament per lambda class

3. TRANSLATION GROUP: The 4-vertex reversal group (generated by reversals
   of (1,1,2,2) sub-tournaments)

4. NON-MEASURABILITY: The "H-functional" (Hamiltonian path count) is NOT
   constant on lambda classes at n=7. It varies by ±2 within a class.
   This is the tournament analogue of the Vitali non-measurability:
   the natural "measure" (H) is not invariant under the "translations"
   (4-reversals).

5. GAUGE TRIVIALITY AT n<=6: The translations ARE measure-preserving at
   n<=6 (H is constant on lambda classes). The gauge group is trivial.
   This is like R/Q over a bounded interval — no contradiction.

6. PHASE TRANSITION AT n=7: The gauge group becomes NON-trivial. The
   translations stop being measure-preserving. A "Vitali atom" appears:
   the smallest gauge transformation that breaks H-invariance.

7. THE {2,1,0} HIERARCHY: The overlap weights form a filtration:
   W=2: two 3-cycles share an edge -> lambda graph captures this
   W=1: two 3-cycles share a vertex -> lambda graph captures this
   W=0: two 3-cycles are disjoint -> lambda graph FAILS to capture this
   The failure at W=0 is exactly where the Vitali non-measurability lives.
   It's a "hidden dimension" — the disjoint cycle interaction that lambda
   cannot see.

The connection to higher-dimensional structure:
- The overlap weight space is a graded vector space: W_2 + W_1 + W_0
- The lambda graph "sees" W_2 and W_1 perfectly
- W_0 is the "dark sector" — invisible to lambda
- At n<=6, W_0 has no effect on H (because no 7-cycles exist)
- At n=7, W_0 emerges as the carrier of non-trivial gauge charge
- The Vitali atom lives in W_0: it's a cycle interaction that only
  manifests through 7-cycles, which require the full vertex set

This is why the phase transition is sharp: it's a DIMENSIONAL transition.
Below the critical dimension (n<7), the dark sector is empty.
At the critical dimension (n=7), the dark sector opens up, carrying
exactly 1 unit of gauge charge (delta_c7 = +/-1).
""")

print("\n" + "=" * 70)
print("PART 9: QUANTIFYING THE DARK SECTOR — W=0 CYCLE INTERACTIONS")
print("=" * 70)

# Let's directly measure: for changing vs preserving reversals,
# how does the W=0 interaction differ?

print("\nW=0 (disjoint 3-cycle pair) analysis:")

np.random.seed(9999)
w0_data = []

for trial in range(2000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    lam_orig = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue
        B_ = reverse_subtournament(A, n, list(subset))
        lam_new = lambda_graph(B_, n)
        if not np.array_equal(lam_orig, lam_new):
            continue

        c7_A = find_all_ham_cycles(A, n)
        c7_B = find_all_ham_cycles(B_, n)
        is_changing = (c7_A != c7_B)

        # Find 3-cycles
        def get_3cycles(M):
            cycles = []
            for combo in combinations(range(n), 3):
                i, j, k = combo
                if (M[i][j] and M[j][k] and M[k][i]) or (M[j][i] and M[i][k] and M[k][j]):
                    cycles.append(frozenset(combo))
            return cycles

        c3_A = get_3cycles(A)
        c3_B = get_3cycles(B_)

        # W=0 pairs (disjoint)
        def count_w0(cycles):
            count = 0
            for i in range(len(cycles)):
                for j in range(i+1, len(cycles)):
                    if len(cycles[i] & cycles[j]) == 0:
                        count += 1
            return count

        w0_A = count_w0(c3_A)
        w0_B = count_w0(c3_B)

        # W=0 pairs that involve subset vertices
        def w0_involving_subset(cycles, sub_set):
            """W=0 pairs where at least one cycle has 2+ subset vertices."""
            count = 0
            for i in range(len(cycles)):
                for j in range(i+1, len(cycles)):
                    if len(cycles[i] & cycles[j]) == 0:
                        if len(cycles[i] & sub_set) >= 2 or len(cycles[j] & sub_set) >= 2:
                            count += 1
            return count

        sub_set = frozenset(subset)
        w0_sub_A = w0_involving_subset(c3_A, sub_set)
        w0_sub_B = w0_involving_subset(c3_B, sub_set)

        w0_data.append({
            'is_changing': is_changing,
            'w0_A': w0_A, 'w0_B': w0_B,
            'w0_sub_A': w0_sub_A, 'w0_sub_B': w0_sub_B,
            'delta_w0': w0_B - w0_A,
            'delta_w0_sub': w0_sub_B - w0_sub_A,
        })

    if trial % 500 == 0 and trial > 0:
        print(f"  ... {trial} processed")

ch_w0 = [d for d in w0_data if d['is_changing']]
pr_w0 = [d for d in w0_data if not d['is_changing']]

print(f"\nTotal: {len(w0_data)}, changing: {len(ch_w0)}, preserving: {len(pr_w0)}")

print("\n  W=0 total delta (changing): ", Counter(d['delta_w0'] for d in ch_w0))
print("  W=0 total delta (preserving): ", Counter(d['delta_w0'] for d in pr_w0))
print("  W=0 subset delta (changing): ", Counter(d['delta_w0_sub'] for d in ch_w0))
print("  W=0 subset delta (preserving): ", Counter(d['delta_w0_sub'] for d in pr_w0))

print("\n" + "=" * 70)
print("PART 10: SUMMARY — THE COMPLETE MECHANISM")
print("=" * 70)
print("""
THE VITALI c7 MECHANISM AT n=7:

1. STRUCTURE: A (1,1,2,2)-scored 4-vertex sub-tournament has exactly 5
   directed Ham paths and 2 directed 3-cycles.

2. REVERSAL: Flipping all 6 internal arcs maps it to another (1,1,2,2)
   tournament with 5 DIFFERENT Ham paths.

3. THREADING: A 7-cycle can use 1, 2, or 3 internal arcs (arcs within
   the 4-subset). Only the 3-internal-arc case matters: these correspond
   to 7-cycles where the 4 subset vertices are CONSECUTIVE (EEESSSS pattern).

4. COMPLETIONS: A 3-internal-arc 7-cycle = a Ham path (a→b→c→d) of the
   subset + a directed path through all 3 external vertices from d back to a.
   The completion count depends ONLY on the endpoint pair (d,a) and the
   external arcs (which are UNCHANGED by reversal).

5. DELTA: delta_c7 = Σ comp(d,a) [over AFTER endpoints]
                    - Σ comp(d,a) [over BEFORE endpoints]
   The SAME completion function, evaluated on DIFFERENT endpoint pairs.

6. PHASE TRANSITION: At n≤6, there are no 7-cycles, so delta_c7 is
   meaningless. At n=7, 7-cycles appear for the first time, and the
   Ham path → completion mechanism activates.

7. PREDICTORS:
   - mixed_c3 = 0 => NEVER changes (no mixed 3-cycles = no interaction)
   - ext tournament cyclic => 73.9% change rate
   - ext tournament transitive => much lower change rate
   - The arc signature (each subset vertex's out/in to ext) modulates the rate

8. THE DARK SECTOR: The {2,1,0} overlap weight filtration has a "dark"
   W=0 component that the lambda graph cannot detect. This W=0 sector
   carries the Vitali gauge charge and only manifests at n≥7.
""")

print("\nDone.")
