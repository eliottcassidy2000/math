#!/usr/bin/env python3
"""
pos_blueblack_formulas_s111.py — Exact Formulas for PoS, Blueself, Blackself
opus-2026-03-21-S111

THE GOAL: Find EXACT formulas for:
1. The SIZE of the PoS class (number of labeled tournaments with the PoS score)
2. The number of BLACKSELF tournaments in the PoS (odd n)
3. The number of BLUESELF tournaments in the PoS (even n)
4. The H DISTRIBUTION within the PoS
5. How these relate to the recursive structure (source-sink-middle)

THE KEY STRUCTURE (from S110):
At n=5, PoS score (1,2,2,2,3):
  Type A (sink→source, cyclic middle): 40 tours, H=15
  Type B (source→sink, cyclic middle): 120 tours, H=11
  Type B (source→sink, transitive middle): 120 tours, H=13

At n=5, the PoS has 280 labeled tournaments.
280 = 5!/2 · C(3,?) ... let me derive this.

COUNTING THE PoS:
The PoS score at n=5 is (1,2,2,2,3).
The number of labeled tournaments with this score = ?

APPROACH: Use the source-sink-middle decomposition.
A tournament with score (1,2,2,2,3) has:
  - One vertex with score 1 (the "sink")
  - One vertex with score 3 (the "source")
  - Three vertices with score 2 (the "middle")

Step 1: Choose which vertex is the sink: n = 5 ways
Step 2: Choose which vertex is the source: n-1 = 4 ways
Step 3: The middle tournament on 3 vertices: 2^3 = 8 possibilities
Step 4: The source-sink arc direction: 2 ways (sink→source or source→sink)

BUT NOT ALL COMBINATIONS GIVE SCORE (1,2,2,2,3).
The constraint is that the scores work out correctly.
"""

from math import comb, factorial
from fractions import Fraction
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

print("=" * 72)
print("  EXACT FORMULAS FOR PoS, BLUESELF, AND BLACKSELF")
print("  opus-2026-03-21-S111")
print("=" * 72)
print()

# =============================================================================
# HELPER
# =============================================================================

def H_dp(adj_flat, n):
    dp = [0] * ((1 << n) * n)
    for v in range(n): dp[(1 << v) * n + v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)): continue
            val = dp[S * n + v]
            if val == 0: continue
            for u in range(n):
                if S & (1 << u): continue
                if adj_flat[v * n + u]:
                    dp[(S | (1 << u)) * n + u] += val
    full = (1 << n) - 1
    return sum(dp[full * n + v] for v in range(n))

# =============================================================================
# PART 1: EXACT PoS CLASS SIZE FORMULA
# =============================================================================

print("=" * 72)
print("  PART 1: COUNTING THE PoS CLASS EXACTLY")
print("=" * 72)
print()

# Build all n=5 tournaments and verify counts
for n in [5, 6]:
    m = comb(n, 2)
    pos_score_candidates = []

    # Find the PoS score: the self-complementary score with minimum nonzero S₂
    all_scores = defaultdict(int)
    for mask in range(2**m):
        scores = [0]*n
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (mask >> idx) & 1: scores[i] += 1
                else: scores[j] += 1
                idx += 1
        sc = tuple(sorted(scores))
        all_scores[sc] += 1

    # Find self-complementary scores with min nonzero S₂
    best_pos = None
    best_S2 = float('inf')
    for sc, count in all_scores.items():
        comp = tuple(sorted(n-1-s for s in sc))
        if sc != comp:
            continue
        S2 = sum((s - (n-1)/2)**2 for s in sc)
        if 0 < S2 < best_S2:
            best_S2 = S2
            best_pos = sc

    print(f"  n={n}: PoS score = {best_pos}, S₂ = {best_S2}")
    print(f"  PoS class size (labeled): {all_scores.get(best_pos, 0)}")

    # Now derive the size combinatorially
    # The PoS score has a specific multiplicity structure
    score_mult = Counter(best_pos)
    print(f"  Score multiplicities: {dict(sorted(score_mult.items()))}")

    # The number of labeled tournaments = n! / (Aut of score) × (# per vertex assignment)
    # Actually, the formula involves a more careful counting...
    # Let's count by the source-sink-middle decomposition.

    # Identify source and sink scores
    min_score = min(best_pos)
    max_score = max(best_pos)
    mid_scores_with_shift = tuple(s for s in best_pos if s != min_score and s != max_score)
    if len(mid_scores_with_shift) != n - 2:
        # Handle case where min or max appears multiple times
        # At n=5: min=1 (once), max=3 (once), mid=2,2,2
        # At n=6: min=1 (once), max=4 (once), mid=2,2,3,3
        remaining = list(best_pos)
        remaining.remove(min_score)
        remaining.remove(max_score)
        mid_scores_with_shift = tuple(sorted(remaining))

    mid_internal = tuple(s - 1 for s in mid_scores_with_shift)
    print(f"  Source score: {max_score}, Sink score: {min_score}")
    print(f"  Middle scores (full): {mid_scores_with_shift}")
    print(f"  Middle internal scores: {mid_internal}")
    print()

    # Count by structure
    # Choose sink vertex: n ways
    # Choose source vertex: n-1 ways
    # Source-sink direction: 2 ways (but constrained by scores)
    # Middle tournament: #{tournaments on n-2 vertices with specific score}

    # Actually, let's count more carefully for each source-sink polarity
    if n == 5:
        # Type A: sink→source. Sink score = 1 (beats only source).
        # Source score = 3 (loses only to sink).
        # Middle: each beats sink, loses to source. Internal score = total - 1.
        # Middle internal must be (1,1,1) = regular C₃.
        # #{C₃ on 3 vertices} = 2 (two orientations of the cycle)

        # But wait: any tournament on 3 vertices with all scores 1 is regular.
        # There are only 2 regular tournaments on 3 vertices (the two C₃).

        n_sink = 5  # Choose sink vertex
        n_source = 4  # Choose source vertex (not sink)
        n_mid_reg = 2  # Regular tournaments on 3 middle vertices
        type_A = n_sink * n_source * n_mid_reg
        print(f"  Type A (sink→source, regular middle): {type_A}")

        # Type B: source→sink. Source score = 3 (beats all middle + sink).
        # Wait: if source→sink, then source beats sink, so source score includes
        # one win from beating sink. But source score is 3 = wins from 3 middle
        # + ... hmm, let me recount.

        # If source→sink: source beats sink and all middle. Score = n-1 = 4.
        # But PoS source has score 3, not 4! So source→sink gives score 4, NOT 3.

        # Wait: at n=5, source has score 3. If source beats sink, source beats
        # sink + some middle. Score 3 = 1(sink) + 2(middle). So source beats
        # 2 of 3 middle vertices.

        # And sink doesn't beat source. Sink has score 1. Sink must beat
        # exactly 1 vertex. If source→sink, sink doesn't beat source.
        # So sink beats exactly 1 middle vertex.

        # This is NOT the simple "source beats all middle" structure!
        # It's more complex.

        print()
        print("  CORRECTION: Source→sink doesn't mean source beats all middle!")
        print("  Let me recount using the ACTUAL tournament data.")
        print()

# =============================================================================
# PART 2: STRUCTURAL DECOMPOSITION — EXACT COUNTS
# =============================================================================

print("=" * 72)
print("  PART 2: EXACT STRUCTURAL DECOMPOSITION AT n=5")
print("=" * 72)
print()

n = 5
m = comb(n, 2)
pos_score = (1, 2, 2, 2, 3)

# Enumerate ALL tournaments in the PoS
pos_tours = []
for mask in range(2**m):
    adj = [[0]*n for _ in range(n)]
    scores = [0]*n
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (mask >> idx) & 1: adj[i][j] = 1; scores[i] += 1
            else: adj[j][i] = 1; scores[j] += 1
            idx += 1
    if tuple(sorted(scores)) != pos_score:
        continue

    adj_flat = [0]*(n*n)
    for i in range(n):
        for j in range(n):
            adj_flat[i*n+j] = adj[i][j]
    H = H_dp(adj_flat, n)

    # Identify sink (score 1) and source (score 3)
    sink = scores.index(1)
    source = scores.index(3)

    # Direction between source and sink
    s2s = adj[sink][source]  # Does sink beat source?

    # Middle vertices
    middle = [v for v in range(n) if v != sink and v != source]

    # How many middle vertices does the source beat?
    source_mid_wins = sum(adj[source][v] for v in middle)

    # How many middle vertices does the sink beat?
    sink_mid_wins = sum(adj[sink][v] for v in middle)

    # Middle internal scores
    mid_internal = tuple(sorted(sum(adj[v][w] for w in middle if w != v) for v in middle))

    # Check if self-flip (blackself)
    # Flip = complement non-backbone arcs
    # For now, just check if T ≅ flip(T)
    flip_adj = [[0]*n for _ in range(n)]
    for i in range(n):
        flip_adj[i][(i+1) % n if i < n-1 else 0] = adj[i][(i+1) % n if i < n-1 else 0]  # keep backbone
    # Actually, the "flip" in tiling terms flips all non-backbone arcs
    # Backbone arcs: (i, i+1) for i=0..n-2, so the path n-1→n-2→...→1→0
    # Hmm, the backbone is a specific path. Let me just compute flip simply:
    # flip(T)[i][j] = 1 - T[i][j] for non-backbone arcs
    # backbone arcs: (k, k-1) for k=1..n-1 (the path n-1→...→0)

    pos_tours.append({
        'mask': mask, 'H': H, 'sink': sink, 'source': source,
        'sink_beats_source': bool(s2s),
        'source_mid_wins': source_mid_wins,
        'sink_mid_wins': sink_mid_wins,
        'mid_internal': mid_internal,
        'scores': tuple(scores),
    })

print(f"  Total PoS tournaments: {len(pos_tours)}")
print()

# Full structural breakdown
print("  STRUCTURAL BREAKDOWN:")
print(f"  {'sink→src':>8} | {'src_mid_w':>9} | {'snk_mid_w':>9} | {'mid_int':>12} | {'H':>3} | {'count':>5}")
print("  " + "-" * 55)

structure_counter = Counter()
for t in pos_tours:
    key = (t['sink_beats_source'], t['source_mid_wins'], t['sink_mid_wins'], t['mid_internal'])
    structure_counter[key] += 1

for key in sorted(structure_counter.keys()):
    sbs, smw, skw, mi = key
    # Find H for this structure
    H_vals = [t['H'] for t in pos_tours
              if (t['sink_beats_source'], t['source_mid_wins'], t['sink_mid_wins'], t['mid_internal']) == key]
    H_val = H_vals[0]  # Should be constant within structure
    H_unique = len(set(H_vals)) == 1
    count = structure_counter[key]
    print(f"  {'T' if sbs else 'F':>8} | {smw:>9} | {skw:>9} | {str(mi):>12} | {H_val:>3} | {count:>5}")

print()

# =============================================================================
# PART 3: THE EXACT COUNTING FORMULA
# =============================================================================

print("=" * 72)
print("  PART 3: DERIVING THE EXACT COUNTING FORMULA")
print("=" * 72)
print()

print("""
  From the structural breakdown at n=5:

  CASE 1: sink→source (sink beats source)
    Source beats all middle, loses to sink → source score = 3
    Sink beats source, loses to all middle → sink score = 1
    Middle: each beats sink, loses to source
    Middle score = internal_score + 1
    → Middle internal must be regular (1,1,1) for PoS score (2,2,2)

    Count: (5 choose sink) × (4 choose source) × (#{regular tournaments on 3})
    = 5 × 4 × 2 = 40

  CASE 2: source→sink (source beats sink)
    Source beats sink + some middle. Source score = 3 = 1(sink) + 2(middle).
    So source beats exactly 2 of 3 middle vertices.
    Sink loses to source, beats exactly 1 middle vertex. Sink score = 1.

    The middle structure is MORE COMPLEX:
    - Source beats 2 middle, loses to 1 middle
    - Sink beats 1 middle, loses to 2 middle
    - The middle's internal + external arcs must give total scores (2,2,2)

    Internal score + (beat sink?) + (beat source?) = total score
    For middle vertex v:
      total = internal + (1 if v beats sink else 0) + (1 if v beats source else 0)
      = internal + (1 - adj[sink][v]) + (1 - adj[source][v])

    Wait: "beats sink" means adj[v][sink]=1, and "beats source" means adj[v][source]=1.
    Sink beats 1 middle → adj[sink][v]=1 for exactly 1 middle v.
    So for that v: adj[v][sink]=0, "doesn't beat sink"
    For other 2 middle: adj[v][sink]=1, "beats sink"

    Source beats 2 middle → adj[source][v]=1 for exactly 2 middle v.
    So for those 2: adj[v][source]=0, "doesn't beat source"
    For the remaining 1: adj[v][source]=1, "beats source"
""")

# Let me just count directly what structures exist
print("  Exact count verification:")
case1 = sum(1 for t in pos_tours if t['sink_beats_source'])
case2 = sum(1 for t in pos_tours if not t['sink_beats_source'])
print(f"    Case 1 (sink→source): {case1}")
print(f"    Case 2 (source→sink): {case2}")
print(f"    Total: {case1 + case2} = {len(pos_tours)}")
print()

# Derive Case 2 count
# Choose sink: 5 ways
# Choose source: 4 ways
# Source beats 2 of 3 middle: C(3,2) = 3 ways
# Sink beats 1 of 3 middle: CONSTRAINED — which one?
# The middle vertex that beats source also might be the one sink beats.

# Actually, let me count by the detailed structure
for sbs, smw, skw, mi in sorted(structure_counter.keys()):
    count = structure_counter[(sbs, smw, skw, mi)]
    # Derive the count from vertex choices
    n_sink_choices = 5
    n_source_choices = 4

    if sbs:  # sink→source
        # Source beats all 3 middle
        # Sink beats only source
        # Middle internal = mi
        # Number of middle tournaments with given internal scores
        if mi == (1, 1, 1):
            n_mid = 2  # Two cyclic C₃
        elif mi == (0, 1, 2):
            n_mid = 6  # Transitive on 3 vertices: 3! = 6? No, only 1 up to iso
            # Actually 3! labelings of the transitive tournament = 6/1 = 6? Hmm.
            # On 3 labeled vertices: there are 6 transitive tournaments (3!)
            n_mid = 6
        else:
            n_mid = '?'
        predicted = n_sink_choices * n_source_choices * (n_mid if isinstance(n_mid, int) else 0)
        print(f"    sink→src={sbs}, smw={smw}, skw={skw}, mi={mi}: actual={count}, "
              f"predicted=5×4×{n_mid}={predicted}")
    else:
        # More complex — need to count the number of ways to assign arcs
        print(f"    sink→src={sbs}, smw={smw}, skw={skw}, mi={mi}: actual={count}")

print()

# =============================================================================
# PART 4: THE PoS SIZE FORMULA
# =============================================================================

print("=" * 72)
print("  PART 4: THE PoS SIZE FORMULA")
print("=" * 72)
print()

print("""
  At n=5, the PoS class (1,2,2,2,3) has 280 tournaments.

  Derivation attempt:
  - Choose the sink and source: 5 × 4 = 20 ordered (sink, source) pairs
  - For each pair, the middle has 3 vertices
  - The middle tournament + source/sink arcs are constrained

  For Case 1 (sink→source):
    - Middle internal must be regular (1,1,1): 2 tournaments
    - Total: 20 × 2 = 40

  For Case 2 (source→sink):
    - Source beats 2 of 3 middle, sink beats 1 of 3 middle
    - Multiple sub-cases depending on overlap
    - Total: 240

  280 = 40 + 240

  Let's verify 240 = 20 × 12:
    For each (sink, source) pair: 12 valid configurations in Case 2.
    12 = C(3,2) × C(3,1) × (valid middle tournaments)
    C(3,2) = 3 (which 2 middle vertices the source beats)
    × something = 12... so "something" = 4.

  Let me count: for a given (sink, source) with source→sink:
""")

# Count Case 2 per (sink, source) pair
per_pair = Counter()
for t in pos_tours:
    if not t['sink_beats_source']:
        per_pair[(t['sink'], t['source'])] += 1

if per_pair:
    counts = list(per_pair.values())
    print(f"  Case 2 tournaments per (sink, source) pair:")
    print(f"    Values: {sorted(set(counts))}")
    print(f"    All same? {len(set(counts)) == 1}")
    if len(set(counts)) == 1:
        print(f"    Count per pair: {counts[0]}")
        print(f"    20 pairs × {counts[0]} = {20 * counts[0]}")

print()

# =============================================================================
# PART 5: BLACKSELF COUNT FORMULA AT n=5
# =============================================================================

print("=" * 72)
print("  PART 5: BLACKSELF COUNT IN THE PoS AT n=5")
print("=" * 72)
print()

print("""
  From S109: 36 blackself tournaments in the PoS at n=5.
  H distribution: {11: 12, 13: 20, 15: 4}

  A tournament is BLACKSELF if its flip (complement of all non-backbone arcs)
  gives an isomorphic tournament, where the isomorphism is NOT grid-symmetric.

  At odd n, ALL self-flip tournaments are blackself (no blueself possible).

  36 out of 280 PoS tournaments are blackself.
  36/280 = 9/70.

  Can we derive 36 from the structure?

  36 = 4 × 9? Or 36 = 6 × 6? Or 36 = C(9,2)?
  36 = 6² = (3!)² = 36.

  Or: the backbone path has n! / n = (n-1)! = 24 possible starting points...
  Hmm, not obvious.
""")

# Verify blackself count by checking self-flip
# The backbone path is (n-1, n-2, ..., 1, 0)
# The tiling bits are non-backbone arcs

backbone_arcs = set()
for i in range(n-1):
    backbone_arcs.add((n-1-i, n-2-i))  # Path n-1→n-2→...→0

# For each PoS tournament, check if flip gives an isomorphic tournament
blackself_count = 0
blackself_by_H = Counter()

for t in pos_tours:
    mask = t['mask']

    # Build adjacency
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (mask >> idx) & 1: adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1

    # Flip: reverse all non-backbone arcs
    flip = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if (i, j) in backbone_arcs:
                flip[i][j] = adj[i][j]  # Keep backbone
            else:
                flip[i][j] = 1 - adj[i][j]  # Flip non-backbone

    # Check isomorphism: does there exist a permutation σ with flip = σ(adj)?
    flip_tuple = tuple(tuple(row) for row in flip)

    is_iso = False
    for perm in permutations(range(n)):
        permuted = tuple(tuple(adj[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if permuted == flip_tuple:
            is_iso = True
            break

    if is_iso:
        blackself_count += 1
        blackself_by_H[t['H']] += 1

print(f"  Blackself count in PoS: {blackself_count}")
print(f"  By H: {dict(sorted(blackself_by_H.items()))}")
print()

# =============================================================================
# PART 6: EXTENDING TO n=6
# =============================================================================

print("=" * 72)
print("  PART 6: PoS AT n=6 — STRUCTURAL DECOMPOSITION")
print("=" * 72)
print()

n = 6
m = comb(n, 2)
pos_score_6 = (1, 2, 2, 3, 3, 4)

# Find all PoS tournaments at n=6
pos_tours_6 = []
t0 = time.time()
for mask in range(2**m):
    adj = [[0]*n for _ in range(n)]
    scores = [0]*n
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (mask >> idx) & 1: adj[i][j] = 1; scores[i] += 1
            else: adj[j][i] = 1; scores[j] += 1
            idx += 1
    if tuple(sorted(scores)) != pos_score_6:
        continue

    adj_flat = [0]*(n*n)
    for i in range(n):
        for j in range(n):
            adj_flat[i*n+j] = adj[i][j]
    H = H_dp(adj_flat, n)

    # Identify sink (score 1) and source (score 4)
    sink = scores.index(1)
    source = scores.index(4)
    sbs = adj[sink][source]

    pos_tours_6.append({
        'mask': mask, 'H': H, 'sink': sink, 'source': source,
        'sink_beats_source': bool(sbs),
        'scores': tuple(scores),
    })

elapsed = time.time() - t0
print(f"  n=6: Found {len(pos_tours_6)} PoS tournaments in {elapsed:.1f}s")
print(f"  H distribution: {dict(sorted(Counter(t['H'] for t in pos_tours_6).items()))}")
print()

# Case 1 vs Case 2
case1_6 = sum(1 for t in pos_tours_6 if t['sink_beats_source'])
case2_6 = sum(1 for t in pos_tours_6 if not t['sink_beats_source'])
print(f"  Case 1 (sink→source): {case1_6}")
print(f"  Case 2 (source→sink): {case2_6}")
print()

# H by case
print(f"  H by case:")
for case_name, case_filter in [("sink→source", True), ("source→sink", False)]:
    H_dist = Counter(t['H'] for t in pos_tours_6 if t['sink_beats_source'] == case_filter)
    print(f"    {case_name}: {dict(sorted(H_dist.items()))}")
print()

# Per (sink, source) pair count
per_pair_6 = Counter()
for t in pos_tours_6:
    per_pair_6[(t['sink'], t['source'], t['sink_beats_source'])] += 1
case1_per_pair = [c for (s,r,sbs), c in per_pair_6.items() if sbs]
case2_per_pair = [c for (s,r,sbs), c in per_pair_6.items() if not sbs]
print(f"  Case 1 per pair: {sorted(set(case1_per_pair)) if case1_per_pair else 'none'}")
print(f"  Case 2 per pair: {sorted(set(case2_per_pair)) if case2_per_pair else 'none'}")
print()

# =============================================================================
# PART 7: THE SIZE FORMULAS
# =============================================================================

print("=" * 72)
print("  PART 7: SIZE FORMULAS — COLLECTING THE PATTERN")
print("=" * 72)
print()

print(f"""
  PoS class sizes:
    n=5: |PoS(5)| = 280
    n=6: |PoS(6)| = {len(pos_tours_6)}

  Decomposition:
    n=5: Case1 = 40, Case2 = 240. Total = 280.
    n=6: Case1 = {case1_6}, Case2 = {case2_6}. Total = {len(pos_tours_6)}.

  PATTERNS:
    n=5: 280 = 5 × 4 × 14 = 20 × 14
    n=6: {len(pos_tours_6)} = 6 × 5 × {len(pos_tours_6)//30 if len(pos_tours_6) % 30 == 0 else '?'}

    Actually: 280 = C(5,1) × C(4,1) × 14
    And 14 = number of tournaments per (sink,source) pair

    Can we express 14 in terms of n?
    At n=5: 14 tournaments per pair = 2 (Case 1) + 12 (Case 2) = 14
    At n=6: {case1_6//30 if case1_6 % 30 == 0 else case1_6/30:.1f} (Case 1 per pair) + {case2_6//30 if case2_6 % 30 == 0 else case2_6/30:.1f} (Case 2 per pair) per pair
""")

# Verify per-pair counts
if case1_per_pair:
    c1pp = case1_per_pair[0] if len(set(case1_per_pair)) == 1 else 'varies'
else:
    c1pp = 0
if case2_per_pair:
    c2pp = case2_per_pair[0] if len(set(case2_per_pair)) == 1 else 'varies'
else:
    c2pp = 0

n_pairs = n * (n-1)  # ordered (sink, source) pairs

print(f"  At n=5: {20} pairs × {14} per pair = {280}")
print(f"  At n=6: {n_pairs} pairs × ({c1pp} + {c2pp}) per pair = {n_pairs * (c1pp + c2pp) if isinstance(c1pp, int) and isinstance(c2pp, int) else '?'}")
print()

# Case 1 per pair at n=5: 2 (two regular C₃ orientations)
# At n=6: ? (number of regular/near-regular tournaments on n-2 vertices)
print(f"  Case 1 per pair: # middle tournaments with 'right' internal scores")
print(f"    n=5: 2 (= #{'{'}regular C₃{'}'} on 3 vertices)")
print(f"    n=6: {c1pp} (= #{'{'}specific tournaments on 4 vertices{'}'})")
print()

# =============================================================================
# PART 8: BLACKSELF/BLUESELF AT n=6
# =============================================================================

print("=" * 72)
print("  PART 8: BLUESELF AND BLACKSELF AT n=6")
print("=" * 72)
print()

# At n=6 (even), BOTH blueself and blackself can exist
# Check self-flip for each PoS tournament
backbone_arcs_6 = set()
for i in range(n-1):
    backbone_arcs_6.add((n-1-i, n-2-i))

blueself_6 = 0
blackself_6 = 0
selfflip_by_H = Counter()

for t in pos_tours_6[:2000]:  # Sample for speed
    mask = t['mask']
    adj = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (mask >> idx) & 1: adj[i][j] = 1
            else: adj[j][i] = 1
            idx += 1

    # Flip
    flip = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if (i, j) in backbone_arcs_6:
                flip[i][j] = adj[i][j]
            else:
                flip[i][j] = 1 - adj[i][j]

    flip_tuple = tuple(tuple(row) for row in flip)
    is_iso = False
    for perm in permutations(range(n)):
        permuted = tuple(tuple(adj[perm[i]][perm[j]] for j in range(n)) for i in range(n))
        if permuted == flip_tuple:
            is_iso = True
            break

    if is_iso:
        selfflip_by_H[t['H']] += 1

print(f"  Self-flip tournaments in PoS at n=6 (sampled {min(2000, len(pos_tours_6))}):")
print(f"    By H: {dict(sorted(selfflip_by_H.items()))}")
print(f"    Total self-flip: {sum(selfflip_by_H.values())}")
print()

# =============================================================================
# PART 9: SYNTHESIS — THE EXACT FORMULAS
# =============================================================================

print("=" * 72)
print("  PART 9: SYNTHESIS — TOWARD EXACT FORMULAS")
print("=" * 72)
print()

print(f"""
  WHAT WE CAN STATE EXACTLY:

  1. PoS CLASS SIZE:
     |PoS(n)| = n(n-1) × T_per_pair(n)
     where T_per_pair(n) = #{'{'}valid middle + arc configurations per (sink,source) pair{'}'}

     n=5: T_per_pair = 14 = 2 + 12 (Case1 + Case2)
     n=6: T_per_pair = {(case1_6 + case2_6) // (n*(n-1)) if (case1_6 + case2_6) % (n*(n-1)) == 0 else '?'}

  2. H DISTRIBUTION within PoS:
     Determined by the source-sink polarity AND the middle tournament type.
     At n=5: 3 H values from 2 structural types × 2 middle types.
     At n=6: more H values from more complex middle structures.

  3. BLACKSELF COUNT:
     At n=5: 36 blackself in PoS (= {36}/{280} = {Fraction(36,280)} of class)
     The blackself fraction might be 1/(n-1) or similar.
     36/280 = 9/70. Hmm.

  4. THE RECURSIVE OBSERVATION:
     PoS(n) has middle = tournament on n-2 vertices
     PoS(n) score = (1, PoS(n-2) scores + 1, n-2)
     This recursion holds for the SCORE SEQUENCE.

     BUT: the H recursion is more complex (depends on polarity).

  5. THE KEY FORMULA TO DERIVE:
     Var(H|PoS(n)) = f(n, Var(H|PoS(n-2)), ...)
     This would give the OCR recursion.
""")

print("opus-2026-03-21-S111: EXACT FORMULAS FOR PoS, BLUESELF, BLACKSELF")
