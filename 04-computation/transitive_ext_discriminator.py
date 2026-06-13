"""
transitive_ext_discriminator.py -- kind-pasteur-2026-03-13-S61
The remaining mystery: when mixed_c3=4, ext_c3=0 (transitive external),
only 29.3% of lambda-preserving (1,1,2,2) reversals change H.
What distinguishes the 29.3% from the 70.7%?

Key insight from Part 5 of previous script: completion count C[d][a]
depends on the "reachability" from d to a through the 3-vertex external
path. For transitive ext, there's only ONE Ham path through the ext,
so C[d][a] = 1 iff d beats the ext-source AND ext-sink beats a.

This means the completion matrix has a specific rank-1 structure for
transitive ext, and the question reduces to an endpoint combinatorial
property of the (1,1,2,2) sub-tournament.

Also: larger sample to confirm the 100% rate for cyclic ext.
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

def sub_scores(A, n, subset):
    k = len(subset)
    scores = []
    for i in range(k):
        s = sum(A[subset[i]][subset[j]] for j in range(k) if i != j)
        scores.append(s)
    return tuple(sorted(scores))

def count_ham_paths(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask_size in range(2, n+1):
        for mask in range(1 << n):
            if bin(mask).count('1') != mask_size:
                continue
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                prev_mask = mask ^ (1 << v)
                total = 0
                for u in range(n):
                    if (prev_mask & (1 << u)) and A[u][v]:
                        total += dp.get((prev_mask, u), 0)
                if total:
                    dp[(mask, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def ham_paths_on_subset(A, subset):
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

def find_all_ham_cycles(A, n):
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

n = 7
total_bits = n * (n - 1) // 2

print("=" * 70)
print("PART 1: LARGE SAMPLE — CONFIRM THE 100% CYCLIC RATE")
print("=" * 70)

np.random.seed(31415)
results = {'cyclic_changing': 0, 'cyclic_preserving': 0,
           'trans_changing': 0, 'trans_preserving': 0}

# Also collect details for transitive cases
trans_details = []

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

        ext = [v for v in range(n) if v not in subset]
        sub_list = list(subset)

        # mixed_c3
        mixed_c3 = 0
        for s1, s2 in combinations(sub_list, 2):
            for e in ext:
                if (A[s1][s2] and A[s2][e] and A[e][s1]) or (A[s2][s1] and A[s1][e] and A[e][s2]):
                    mixed_c3 += 1
        if mixed_c3 == 0:
            continue  # these never change, skip

        # ext_c3
        i, j, k = ext
        ext_c3 = 1 if ((A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j])) else 0

        H_orig = count_ham_paths(A, n)
        H_new = count_ham_paths(B, n)
        is_changing = (H_orig != H_new)

        if ext_c3 == 1:
            if is_changing:
                results['cyclic_changing'] += 1
            else:
                results['cyclic_preserving'] += 1
        else:
            if is_changing:
                results['trans_changing'] += 1
            else:
                results['trans_preserving'] += 1

            # Collect detail for transitive case
            # What is the structure of the external transitive ordering?
            # Find source and sink of ext
            ext_scores_map = {e: sum(A[e][f] for f in ext if f != e) for e in ext}
            ext_source = [e for e in ext if ext_scores_map[e] == 2][0]  # beats both
            ext_middle = [e for e in ext if ext_scores_map[e] == 1][0]
            ext_sink = [e for e in ext if ext_scores_map[e] == 0][0]  # beaten by both

            # The unique Ham path through ext is: source -> middle -> sink
            # Completion C[d][a] = 1 iff A[d][source] and A[sink][a]

            # Who does d beat? d must beat ext_source
            d_beats_source = {s: A[s][ext_source] for s in sub_list}
            a_beaten_by_sink = {s: A[ext_sink][s] for s in sub_list}

            # Ham paths before and after
            hp_before = ham_paths_on_subset(A, sub_list)
            hp_after = ham_paths_on_subset(B, sub_list)

            # Endpoint pairs
            ep_before = [(p[-1], p[0]) for p in hp_before]
            ep_after = [(p[-1], p[0]) for p in hp_after]

            # Completions before and after
            comp_before = sum(1 for d, a in ep_before if d_beats_source[d] and a_beaten_by_sink[a])
            comp_after = sum(1 for d, a in ep_after if d_beats_source[d] and a_beaten_by_sink[a])

            # What subset vertices beat source / are beaten by sink?
            beats_source = [s for s in sub_list if d_beats_source[s]]
            beaten_by_sink = [s for s in sub_list if a_beaten_by_sink[s]]

            # Internal scores
            int_scores = {s: sum(A[s][t] for t in sub_list if t != s) for s in sub_list}
            low = [s for s in sub_list if int_scores[s] <= 1]
            high = [s for s in sub_list if int_scores[s] >= 2]

            trans_details.append({
                'is_changing': is_changing,
                'beats_source': beats_source,
                'beaten_by_sink': beaten_by_sink,
                'low': low, 'high': high,
                'comp_before': comp_before, 'comp_after': comp_after,
                'delta_comp': comp_after - comp_before,
                'n_beats_source': len(beats_source),
                'n_beaten_by_sink': len(beaten_by_sink),
                'source_from_low': sum(1 for s in beats_source if s in low),
                'source_from_high': sum(1 for s in beats_source if s in high),
                'sink_to_low': sum(1 for s in beaten_by_sink if s in low),
                'sink_to_high': sum(1 for s in beaten_by_sink if s in high),
            })

    if trial % 1000 == 0 and trial > 0:
        print(f"  ... {trial}/5000, cyclic={results['cyclic_changing']}+{results['cyclic_preserving']}, "
              f"trans={results['trans_changing']}+{results['trans_preserving']}")

print(f"\nCyclic external: changing={results['cyclic_changing']}, preserving={results['cyclic_preserving']}")
c_total = results['cyclic_changing'] + results['cyclic_preserving']
if c_total > 0:
    print(f"  Rate: {results['cyclic_changing']/c_total*100:.1f}%")
    print(f"  {'CONFIRMED 100%' if results['cyclic_preserving'] == 0 else 'NOT 100%!'}")

print(f"\nTransitive external: changing={results['trans_changing']}, preserving={results['trans_preserving']}")
t_total = results['trans_changing'] + results['trans_preserving']
if t_total > 0:
    print(f"  Rate: {results['trans_changing']/t_total*100:.1f}%")

print("\n" + "=" * 70)
print("PART 2: TRANSITIVE EXT — WHAT DISTINGUISHES CHANGING FROM PRESERVING?")
print("=" * 70)

ch_trans = [d for d in trans_details if d['is_changing']]
pr_trans = [d for d in trans_details if not d['is_changing']]

print(f"\nTotal transitive: {len(trans_details)}, changing: {len(ch_trans)}, preserving: {len(pr_trans)}")

# Feature: how many subset vertices beat ext_source?
print("\n--- n_beats_source (how many subset vertices beat ext-source) ---")
ch_counts = Counter(d['n_beats_source'] for d in ch_trans)
pr_counts = Counter(d['n_beats_source'] for d in pr_trans)
for k in sorted(set(list(ch_counts.keys()) + list(pr_counts.keys()))):
    c = ch_counts.get(k, 0)
    p = pr_counts.get(k, 0)
    pct = c/(c+p)*100 if c+p > 0 else 0
    print(f"  n={k}: changing={c}, preserving={p}, rate={pct:.1f}%")

# Feature: how many subset vertices are beaten by ext_sink?
print("\n--- n_beaten_by_sink ---")
ch_counts = Counter(d['n_beaten_by_sink'] for d in ch_trans)
pr_counts = Counter(d['n_beaten_by_sink'] for d in pr_trans)
for k in sorted(set(list(ch_counts.keys()) + list(pr_counts.keys()))):
    c = ch_counts.get(k, 0)
    p = pr_counts.get(k, 0)
    pct = c/(c+p)*100 if c+p > 0 else 0
    print(f"  n={k}: changing={c}, preserving={p}, rate={pct:.1f}%")

# Cross-tab: (n_beats_source, n_beaten_by_sink)
print("\n--- Cross-tab (n_beats_source, n_beaten_by_sink) ---")
cross = Counter()
for d in trans_details:
    key = (d['n_beats_source'], d['n_beaten_by_sink'], d['is_changing'])
    cross[key] += 1
for bs in sorted(set(d['n_beats_source'] for d in trans_details)):
    for bk in sorted(set(d['n_beaten_by_sink'] for d in trans_details)):
        c = cross.get((bs, bk, True), 0)
        p = cross.get((bs, bk, False), 0)
        if c + p > 0:
            print(f"  beats_source={bs}, beaten_by_sink={bk}: "
                  f"changing={c}, preserving={p}, rate={c/(c+p)*100:.1f}%")

# Feature: do source-beaters overlap with low or high vertices?
print("\n--- source_from_low, source_from_high ---")
for combo in [(0,0), (0,1), (0,2), (1,0), (1,1), (1,2), (2,0), (2,1), (2,2)]:
    subset_data = [d for d in trans_details if d['source_from_low'] == combo[0] and d['source_from_high'] == combo[1]]
    if subset_data:
        c = sum(1 for d in subset_data if d['is_changing'])
        p = len(subset_data) - c
        print(f"  low_beat_source={combo[0]}, high_beat_source={combo[1]}: "
              f"changing={c}, preserving={p}, rate={c/(c+p)*100:.1f}%")

# Feature: delta_comp values
print("\n--- delta_comp distribution ---")
print(f"  Changing: {Counter(d['delta_comp'] for d in ch_trans)}")
print(f"  Preserving: {Counter(d['delta_comp'] for d in pr_trans)}")

print("\n" + "=" * 70)
print("PART 3: THE RANK-1 COMPLETION MATRIX ANALYSIS")
print("=" * 70)

# For transitive ext, C[d][a] = 1 iff d beats source AND sink beats a
# So the completion matrix is rank-1: C = v_out * v_in^T
# where v_out[d] = A[d][source], v_in[a] = A[sink][a]
#
# Total completions = sum over endpoints (d,a) of v_out[d] * v_in[a]
# = #{(d,a) endpoint pair : d beats source AND sink beats a}
#
# This is a BILINEAR FORM on the endpoint set.
# The change is:
# delta = sum_{(d,a) in AFTER} v_out[d]*v_in[a] - sum_{(d,a) in BEFORE} v_out[d]*v_in[a]
#
# Let f(S) = sum_{(d,a) endpoints of ham-paths-of-S} v_out[d]*v_in[a]
# Then delta = f(S') - f(S) where S' is the reversed sub-tournament

print("\nFor transitive ext: C[d][a] is rank-1 (v_out * v_in^T)")
print("delta_comp = f(S') - f(S) where f = bilinear form on HP endpoint pairs\n")

# Can we express f in terms of the sub-tournament structure?
# The endpoints of the 5 ham paths of S depend on which vertices are
# "sources" and "sinks" of S.
# S has scores (1,1,2,2). The two score-2 vertices are "highs", score-1 are "lows".
# HP endpoints: the START tends to be high (more out-arcs -> more paths starting there)
# and END tends to be low (more in-arcs -> more paths ending there)?

# Let's check: for each (1,1,2,2) tournament, what are the start/end distributions?
print("Ham path endpoint distribution for (1,1,2,2) tournaments:")

for bits in range(1 << 6):
    A4 = bits_to_adj(bits, 4)
    scores = tuple(sorted([sum(A4[i]) for i in range(4)]))
    if scores != (1, 1, 2, 2):
        continue

    hp = ham_paths_on_subset(A4, [0,1,2,3])
    starts = [p[0] for p in hp]
    ends = [p[-1] for p in hp]

    int_scores = {i: sum(A4[i]) for i in range(4)}
    low = [i for i in range(4) if int_scores[i] <= 1]
    high = [i for i in range(4) if int_scores[i] >= 2]

    start_low = sum(1 for s in starts if s in low)
    start_high = sum(1 for s in starts if s in high)
    end_low = sum(1 for e in ends if e in low)
    end_high = sum(1 for e in ends if e in high)

    print(f"  bits={bits}: low={low}, high={high}, "
          f"starts: {start_low}L+{start_high}H, ends: {end_low}L+{end_high}H")

    # After reversal
    B4 = reverse_subtournament(A4, 4, [0,1,2,3])
    hp_rev = ham_paths_on_subset(B4, [0,1,2,3])
    starts_rev = [p[0] for p in hp_rev]
    ends_rev = [p[-1] for p in hp_rev]

    # After reversal, low and high SWAP
    new_int_scores = {i: sum(B4[i]) for i in range(4)}
    new_low = [i for i in range(4) if new_int_scores[i] <= 1]
    new_high = [i for i in range(4) if new_int_scores[i] >= 2]

    start_low_rev = sum(1 for s in starts_rev if s in low)  # using ORIGINAL low
    start_high_rev = sum(1 for s in starts_rev if s in low)

    # Just show first few
    if bits > 20:
        break

print("\n" + "=" * 70)
print("PART 4: THE EXACT CONDITION FOR TRANSITIVE EXT")
print("=" * 70)

# For transitive ext with source=e_s, middle=e_m, sink=e_k:
# Unique Ham path through ext: e_s -> e_m -> e_k
# C[d][a] = 1 iff A[d][e_s]=1 (d beats source) AND A[e_k][a]=1 (sink beats a)
#
# Let B_s = {subset vertices that beat e_s}
# Let B_k = {subset vertices beaten by e_k}
#
# Then comp = #{HP with endpoint pair (d,a) : d in B_s, a in B_k}
#
# For the (1,1,2,2) sub-tournament T:
# f(T) = #{HP (p_0,...,p_3) : p_3 in B_s, p_0 in B_k}
#
# After reversal to T':
# f(T') = #{HP' (p_0,...,p_3) : p_3 in B_s, p_0 in B_k}
#
# delta = f(T') - f(T)
#
# NOTE: B_s and B_k are the SAME before/after (external arcs unchanged).
# Only the HP set changes.

# Let's compute |B_s|, |B_k| for all transitive-ext cases
print("\nDistribution of |B_s| and |B_k| for transitive ext:")
print("(B_s = subset vertices beating ext-source, B_k = beaten by ext-sink)")

bs_bk_change = Counter()
bs_bk_preserve = Counter()

for d in trans_details:
    key = (d['n_beats_source'], d['n_beaten_by_sink'])
    if d['is_changing']:
        bs_bk_change[key] += 1
    else:
        bs_bk_preserve[key] += 1

all_keys = sorted(set(list(bs_bk_change.keys()) + list(bs_bk_preserve.keys())))
for key in all_keys:
    c = bs_bk_change.get(key, 0)
    p = bs_bk_preserve.get(key, 0)
    pct = c/(c+p)*100 if c+p > 0 else 0
    print(f"  |B_s|={key[0]}, |B_k|={key[1]}: changing={c}, preserving={p}, rate={pct:.1f}%")

# Now: which B_s, B_k subsets lead to changes?
# The key constraint: |B_s| and |B_k| must be EXACTLY 2 for H to potentially change.
# (From arc_sig analysis: all subset vertices have (2,1) or (1,2) to ext.
#  For transitive ext, if each subset vertex beats exactly 1 of 3 ext:
#  beating source is harder than beating sink. But our data shows uniform (2,1).)
#
# Actually, from the arc_sig: for CHANGING cases, each subset vertex has
# exactly 2 arcs out to ext and 1 arc in. For 3 ext vertices with
# transitive order source > middle > sink:
# A subset vertex beats source iff it beats middle and sink as well?
# No! Arcs from subset to ext are independent of the ext internal order.

# Let's check: do B_s and B_k overlap with low/high?
print("\nB_s intersection with low/high vertices:")
for d in trans_details[:20]:
    bs = set(d['beats_source'])
    bk = set(d['beaten_by_sink'])
    low_set = set(d['low'])
    high_set = set(d['high'])

    bs_low = bs & low_set
    bs_high = bs & high_set
    bk_low = bk & low_set
    bk_high = bk & high_set

    print(f"  {'CHANGE' if d['is_changing'] else 'KEEP  '}: "
          f"|B_s|={len(bs)}, B_s_low={len(bs_low)}, B_s_high={len(bs_high)}, "
          f"|B_k|={len(bk)}, B_k_low={len(bk_low)}, B_k_high={len(bk_high)}, "
          f"delta={d['delta_comp']}")

print("\n" + "=" * 70)
print("PART 5: REFINED CROSS-TAB — B_s/B_k vs LOW/HIGH")
print("=" * 70)

# For transitive ext, the question is: where do B_s and B_k land
# relative to the low/high partition?
# There are limited possibilities since |B_s| is constrained.

refined_data = Counter()
for d in trans_details:
    bs = set(d['beats_source'])
    bk = set(d['beaten_by_sink'])
    low_set = set(d['low'])
    high_set = set(d['high'])

    key = (len(bs & low_set), len(bs & high_set),
           len(bk & low_set), len(bk & high_set),
           d['is_changing'])
    refined_data[key] += 1

print(f"\n{'bs_low':>6} {'bs_hi':>5} {'bk_low':>6} {'bk_hi':>5} {'changing':>8} {'preserv':>8} {'rate':>8}")
for bs_l, bs_h, bk_l, bk_h in sorted(set((k[0],k[1],k[2],k[3]) for k in refined_data.keys())):
    c = refined_data.get((bs_l, bs_h, bk_l, bk_h, True), 0)
    p = refined_data.get((bs_l, bs_h, bk_l, bk_h, False), 0)
    if c + p > 0:
        pct = c / (c+p) * 100
        print(f"{bs_l:>6} {bs_h:>5} {bk_l:>6} {bk_h:>5} {c:>8} {p:>8} {pct:>7.1f}%")

print("\n" + "=" * 70)
print("PART 6: THE COMPLETION FORMULA AS BILINEAR PAIRING")
print("=" * 70)

# For a (1,1,2,2) tournament T on subset S = {a,b,c,d}:
# Let (l1, l2) = low vertices, (h1, h2) = high vertices
# T has 5 Ham paths. Let E(T) = set of endpoint pairs (end, start)
#
# For transitive ext: comp = |{(d,a) in E(T) : d in B_s, a in B_k}|
# = |E(T) cap (B_s x B_k)|
#
# After reversal: comp' = |E(T') cap (B_s x B_k)|
# delta = comp' - comp
#
# The question: when does this bilinear pairing change?
#
# Key structural fact: reversal swaps low <-> high.
# So if B_s = {l1, h2} (one low, one high), after reversal the
# "same vertex" is now high and low respectively, but B_s doesn't change!
#
# The Ham path ENDPOINTS change though. Before: high vertices appear
# more often as starts (2 of 5 start at each high vs 0.5 at each low).
# After reversal: the new highs (former lows) appear more as starts.

# Let me verify this start/end distribution claim
print("Start/end frequencies by low/high status:")

for bits in [4, 12, 17]:  # representative (1,1,2,2) tournaments
    A4 = bits_to_adj(bits, 4)
    scores = tuple(sorted([sum(A4[i]) for i in range(4)]))
    if scores != (1, 1, 2, 2):
        continue

    int_scores = {i: sum(A4[i]) for i in range(4)}
    low = [i for i in range(4) if int_scores[i] <= 1]
    high = [i for i in range(4) if int_scores[i] >= 2]

    hp = ham_paths_on_subset(A4, [0,1,2,3])
    ep = [(p[-1], p[0]) for p in hp]

    # Count: how often does each vertex appear as end? as start?
    end_count = Counter(e for e, s in ep)
    start_count = Counter(s for e, s in ep)

    print(f"\n  bits={bits}: low={low}, high={high}")
    for v in range(4):
        role = 'L' if v in low else 'H'
        print(f"    v={v} ({role}): starts={start_count[v]}, ends={end_count[v]}")

    # Same for reversed
    B4 = reverse_subtournament(A4, 4, [0,1,2,3])
    hp_rev = ham_paths_on_subset(B4, [0,1,2,3])
    ep_rev = [(p[-1], p[0]) for p in hp_rev]

    new_int = {i: sum(B4[i]) for i in range(4)}
    new_low = [i for i in range(4) if new_int[i] <= 1]
    new_high = [i for i in range(4) if new_int[i] >= 2]

    end_count_rev = Counter(e for e, s in ep_rev)
    start_count_rev = Counter(s for e, s in ep_rev)

    print(f"  After reversal: new_low={new_low}, new_high={new_high}")
    for v in range(4):
        role = 'L' if v in new_low else 'H'
        print(f"    v={v} ({role}): starts={start_count_rev[v]}, ends={end_count_rev[v]}")

print("\n" + "=" * 70)
print("PART 7: SEARCH FOR PERFECT TRANSITIVE PREDICTOR")
print("=" * 70)

# Let me add more features specific to the transitive case

# Feature: does the ext-source beat any low vertex?
# Feature: does the ext-sink lose to any high vertex?
# Feature: B_s = subset of low vertices only? high only? mixed?

trans_features = []
for d in trans_details:
    bs = set(d['beats_source'])
    bk = set(d['beaten_by_sink'])
    low_set = set(d['low'])
    high_set = set(d['high'])

    # B_s type: 'all_low', 'all_high', 'mixed', 'empty'
    if len(bs) == 0:
        bs_type = 'empty'
    elif bs <= low_set:
        bs_type = 'all_low'
    elif bs <= high_set:
        bs_type = 'all_high'
    else:
        bs_type = 'mixed'

    # Same for B_k
    if len(bk) == 0:
        bk_type = 'empty'
    elif bk <= low_set:
        bk_type = 'all_low'
    elif bk <= high_set:
        bk_type = 'all_high'
    else:
        bk_type = 'mixed'

    trans_features.append({
        'is_changing': d['is_changing'],
        'bs_type': bs_type,
        'bk_type': bk_type,
        'n_bs': len(bs), 'n_bk': len(bk),
        'delta': d['delta_comp'],
    })

# Cross-tab: bs_type x bk_type
print(f"\n{'bs_type':>10} {'bk_type':>10} {'changing':>10} {'preserv':>10} {'rate':>8}")
cross2 = Counter()
for f in trans_features:
    cross2[(f['bs_type'], f['bk_type'], f['is_changing'])] += 1

for bt in ['empty', 'all_low', 'all_high', 'mixed']:
    for kt in ['empty', 'all_low', 'all_high', 'mixed']:
        c = cross2.get((bt, kt, True), 0)
        p = cross2.get((bt, kt, False), 0)
        if c + p > 0:
            pct = c/(c+p)*100
            print(f"{bt:>10} {kt:>10} {c:>10} {p:>10} {pct:>7.1f}%")

# Does the bilinear pairing have a deterministic formula?
# comp = |E(T) intersect (B_s x B_k)|
# For B_s = {l1, l2} (both low) and B_k = {h1, h2} (both high):
# Need endpoint pairs (d,a) where d is low and a is high.
# How many such pairs exist?

print("\n--- Endpoint pairs by (end_type, start_type) ---")
for bits in [4, 12, 17]:
    A4 = bits_to_adj(bits, 4)
    scores = tuple(sorted([sum(A4[i]) for i in range(4)]))
    if scores != (1, 1, 2, 2):
        continue

    int_scores = {i: sum(A4[i]) for i in range(4)}
    low = set(i for i in range(4) if int_scores[i] <= 1)
    high = set(i for i in range(4) if int_scores[i] >= 2)

    hp = ham_paths_on_subset(A4, [0,1,2,3])
    ep = [(p[-1], p[0]) for p in hp]

    types = Counter()
    for e, s in ep:
        e_type = 'L' if e in low else 'H'
        s_type = 'L' if s in low else 'H'
        types[(e_type, s_type)] += 1

    print(f"  bits={bits}: {dict(types)}")

    # Reversed
    B4 = reverse_subtournament(A4, 4, [0,1,2,3])
    hp_rev = ham_paths_on_subset(B4, [0,1,2,3])
    ep_rev = [(p[-1], p[0]) for p in hp_rev]

    # Using ORIGINAL low/high labels
    types_rev = Counter()
    for e, s in ep_rev:
        e_type = 'L' if e in low else 'H'
        s_type = 'L' if s in low else 'H'
        types_rev[(e_type, s_type)] += 1

    print(f"    After: {dict(types_rev)}")
    print(f"    Delta (LL): {types_rev.get(('L','L'),0) - types.get(('L','L'),0)}, "
          f"(LH): {types_rev.get(('L','H'),0) - types.get(('L','H'),0)}, "
          f"(HL): {types_rev.get(('H','L'),0) - types.get(('H','L'),0)}, "
          f"(HH): {types_rev.get(('H','H'),0) - types.get(('H','H'),0)}")

print("\n" + "=" * 70)
print("PART 8: UNIVERSAL ENDPOINT TYPE DISTRIBUTION")
print("=" * 70)

# Check: is the (end_type, start_type) distribution UNIVERSAL for all (1,1,2,2)?
all_ep_types = []
for bits in range(1 << 6):
    A4 = bits_to_adj(bits, 4)
    scores = tuple(sorted([sum(A4[i]) for i in range(4)]))
    if scores != (1, 1, 2, 2):
        continue

    int_scores = {i: sum(A4[i]) for i in range(4)}
    low = set(i for i in range(4) if int_scores[i] <= 1)
    high = set(i for i in range(4) if int_scores[i] >= 2)

    hp = ham_paths_on_subset(A4, [0,1,2,3])
    ep = [(p[-1], p[0]) for p in hp]

    types = Counter()
    for e, s in ep:
        e_type = 'L' if e in low else 'H'
        s_type = 'L' if s in low else 'H'
        types[(e_type, s_type)] += 1

    type_sig = tuple(sorted(types.items()))
    all_ep_types.append(type_sig)

unique_sigs = Counter(tuple(s) for s in all_ep_types)
print(f"Universal endpoint type distribution across all 24 labeled (1,1,2,2) tournaments:")
for sig, count in unique_sigs.items():
    print(f"  {dict(sig)}: {count} tournaments")

print("\n  This means: for ANY (1,1,2,2) tournament:")
print("  - 1 endpoint pair has (end=L, start=L)")
print("  - 1 endpoint pair has (end=H, start=L)")
print("  - 1 endpoint pair has (end=L, start=H)")
print("  - 2 endpoint pairs have (end=H, start=H)")
print("  (Or the symmetric variant)")

print("\n  After reversal (L and H swap roles):")
print("  - 2 pairs (L,L), 1 pair (L,H), 1 pair (H,L), 1 pair (H,H)")
print("  i.e., (L,L) goes 1->2, (H,H) goes 2->1")
print()
print("  So for B_s x B_k bilinear form:")
print("  If B_s=low, B_k=low: delta = (end=L,start=L) change = +1")
print("  If B_s=high, B_k=high: delta = (end=H,start=H) change = -1")
print("  If B_s=low, B_k=high: delta = (end=L,start=H) change = 0")
print("  If B_s=high, B_k=low: delta = (end=H,start=L) change = 0")
print("  Mixed B_s/B_k: depends on specifics")

# Verify with our data
print("\n  VERIFICATION against actual data:")
for f in trans_features:
    predicted = 0
    if f['bs_type'] == 'all_low' and f['bk_type'] == 'all_low':
        predicted = 1
    elif f['bs_type'] == 'all_high' and f['bk_type'] == 'all_high':
        predicted = -1
    actual = f['delta']
    if predicted != 0:
        match = "OK" if predicted == actual else "MISMATCH"
        if match == "MISMATCH":
            print(f"  MISMATCH: bs={f['bs_type']}, bk={f['bk_type']}, predicted={predicted}, actual={actual}")

# Count how many cases fall into each category
print("\n  Category counts:")
cats = Counter()
for f in trans_features:
    cat = (f['bs_type'], f['bk_type'])
    cats[cat] += 1
for cat, count in sorted(cats.items()):
    print(f"    B_s={cat[0]}, B_k={cat[1]}: {count}")

print("\nDone.")
