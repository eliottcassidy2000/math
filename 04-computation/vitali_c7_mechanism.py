"""
vitali_c7_mechanism.py — kind-pasteur-2026-03-13-S61
WHY does the (1,1,2,2) 4-vertex reversal change c7 at n=7?
What distinguishes the 28.6% that change H from the 71.4% that don't?

The key insight to find: the 4 reversed vertices + 3 remaining vertices form
7-cycles. The reversal changes WHICH 7-cycles exist. We need to understand
the interaction pattern.

Analysis plan:
1. For each H-changing Vitali pair, identify the 4-vertex subset and 3 external vertices
2. Classify the arc pattern between the 4-subset and the 3 external vertices
3. Compare with H-preserving reversals — what's different?
4. Count exactly which 7-cycles appear/disappear
5. Look for the geometric/combinatorial obstruction at n<=6
6. Explore the connection to overlap weights W=0,1,2
7. Check if external vertex tournament type matters
8. Examine the "threading" of 7-cycles through the 4-subset
9. Look for a cohomological obstruction (boundary of what?)
10. Test whether the H-change can be predicted from LOCAL data
"""

import numpy as np
from itertools import combinations, permutations
import sys

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

def adj_to_bits(A, n):
    bits = 0
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if A[i][j] == 1:
                bits |= (1 << idx)
            idx += 1
    return bits

def count_ham_paths(A, n):
    """Count Hamiltonian paths via DP."""
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

def count_directed_cycles(A, n, length):
    """Count directed cycles of given length."""
    if length > n:
        return 0
    count = 0
    for combo in combinations(range(n), length):
        sub = [[A[combo[i]][combo[j]] for j in range(length)] for i in range(length)]
        # Count Hamiltonian cycles on this subgraph
        # Fix vertex 0 as start to avoid overcounting
        for perm in permutations(range(1, length)):
            path = (0,) + perm
            valid = True
            for k in range(length):
                if sub[path[k]][path[(k+1) % length]] != 1:
                    valid = False
                    break
            if valid:
                count += 1
    return count

def count_cycles_on_subset(A, n, subset):
    """Count directed Hamiltonian cycles on a specific vertex subset."""
    k = len(subset)
    sub = [[A[subset[i]][subset[j]] for j in range(k)] for i in range(k)]
    count = 0
    for perm in permutations(range(1, k)):
        path = (0,) + perm
        valid = True
        for i in range(k):
            if sub[path[i]][path[(i+1) % k]] != 1:
                valid = False
                break
        if valid:
            count += 1
    return count

def reverse_subtournament(A, n, subset):
    """Reverse all arcs within the subset."""
    B = A.copy()
    for i in subset:
        for j in subset:
            if i != j:
                B[i][j] = A[j][i]
    return B

def lambda_graph(A, n):
    """Compute lambda_{uv} = # 3-cycles containing both u and v."""
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(n):
            if u == v:
                continue
            for w in range(n):
                if w == u or w == v:
                    continue
                # Check if u,v,w form a 3-cycle
                sub = [A[u][v], A[v][w], A[w][u], A[v][u], A[u][w], A[w][v]]
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
    # Each cycle counted once per pair orientation, divide by... actually
    # lambda_{uv} counts cycles, each has one direction through u->v or v->u
    # Actually for undirected lambda, we count the number of 3-cycles containing edge {u,v}
    # A 3-cycle uvw appears as u->v->w->u or u->w->v->u (one is the actual direction)
    # Standard: lambda_{uv} = #{w : uvw form a 3-cycle} = #{w : (u->v->w->u) or (v->u->w->v)}
    # But actually lambda_{uv} = |N^+(u) ∩ N^+(v)| + |N^-(u) ∩ N^-(v)|... no.
    # Simpler: lambda_{uv} = #{w : {u,v,w} forms a 3-cycle (directed)}
    # A directed 3-cycle on {u,v,w}: exactly 2 of the 3 possible directed 3-cycles
    # Let's just count directly.
    L2 = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            cnt = 0
            for w in range(n):
                if w == u or w == v:
                    continue
                # Does {u,v,w} support a directed 3-cycle?
                edges = (A[u][v], A[v][w], A[w][u], A[v][u], A[u][w], A[w][v])
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    cnt += 1
            L2[u][v] = L2[v][u] = cnt
    return L2

def sub_scores(A, n, subset):
    """Get scores of vertices within subset."""
    k = len(subset)
    scores = []
    for i in range(k):
        s = sum(A[subset[i]][subset[j]] for j in range(k) if i != j)
        scores.append(s)
    return tuple(sorted(scores))

def external_arc_pattern(A, n, subset):
    """Get the arc pattern between subset and external vertices."""
    ext = [v for v in range(n) if v not in subset]
    # For each external vertex, record (in_from_subset, out_to_subset)
    pattern = []
    for e in ext:
        in_count = sum(A[s][e] for s in subset)
        out_count = sum(A[e][s] for s in subset)
        pattern.append((in_count, out_count))
    return tuple(sorted(pattern)), ext

def external_tournament_type(A, n, ext_vertices):
    """Get the tournament type on external vertices."""
    k = len(ext_vertices)
    scores = []
    for i in range(k):
        s = sum(A[ext_vertices[i]][ext_vertices[j]] for j in range(k) if i != j)
        scores.append(s)
    return tuple(sorted(scores))

print("=" * 70)
print("PART 1: ANATOMY OF H-CHANGING vs H-PRESERVING 4-REVERSALS AT n=7")
print("=" * 70)

n = 7
total_bits = n * (n - 1) // 2
# We'll sample and find Vitali pairs
np.random.seed(42)

changing_data = []  # (bits, subset, delta_H, ext_pattern, ext_tour_type, c7_before, c7_after)
preserving_data = []  # same structure

sample_size = 2000
for trial in range(sample_size):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    H_orig = count_ham_paths(A, n)
    lam_orig = lambda_graph(A, n)

    for subset in combinations(range(n), 4):
        ss = sub_scores(A, n, list(subset))
        if ss != (1, 1, 2, 2):
            continue

        B = reverse_subtournament(A, n, list(subset))
        lam_new = lambda_graph(B, n)

        if not np.array_equal(lam_orig, lam_new):
            continue

        H_new = count_ham_paths(B, n)
        c7_orig = count_cycles_on_subset(A, n, list(range(n)))  # full 7-cycle
        c7_new = count_cycles_on_subset(B, n, list(range(n)))

        ext_pat, ext_verts = external_arc_pattern(A, n, list(subset))
        ext_type = external_tournament_type(A, n, ext_verts)

        record = (bits, subset, H_new - H_orig, ext_pat, ext_type, c7_orig, c7_new)

        if H_new != H_orig:
            changing_data.append(record)
        else:
            preserving_data.append(record)

    if trial % 200 == 0 and trial > 0:
        print(f"  ... {trial}/{sample_size} tournaments processed, "
              f"{len(changing_data)} changing, {len(preserving_data)} preserving")

print(f"\nTotal H-changing reversals: {len(changing_data)}")
print(f"Total H-preserving reversals: {len(preserving_data)}")
if changing_data:
    print(f"Ratio changing: {len(changing_data)/(len(changing_data)+len(preserving_data)):.3f}")

# Analyze external arc patterns
print("\n--- External arc patterns (sorted (in,out) per external vertex) ---")
from collections import Counter

if changing_data:
    ch_patterns = Counter(r[3] for r in changing_data)
    pr_patterns = Counter(r[3] for r in preserving_data)

    all_patterns = set(list(ch_patterns.keys()) + list(pr_patterns.keys()))
    print(f"\n{'Pattern':<30} {'Changing':>10} {'Preserving':>10} {'Change%':>10}")
    for pat in sorted(all_patterns):
        ch = ch_patterns.get(pat, 0)
        pr = pr_patterns.get(pat, 0)
        pct = ch / (ch + pr) * 100 if (ch + pr) > 0 else 0
        print(f"{str(pat):<30} {ch:>10} {pr:>10} {pct:>9.1f}%")

# External tournament type
print("\n--- External tournament type (3-vertex sub-tournament scores) ---")
if changing_data:
    ch_ext = Counter(r[4] for r in changing_data)
    pr_ext = Counter(r[4] for r in preserving_data)

    all_ext = set(list(ch_ext.keys()) + list(pr_ext.keys()))
    print(f"\n{'Ext Type':<20} {'Changing':>10} {'Preserving':>10} {'Change%':>10}")
    for et in sorted(all_ext):
        ch = ch_ext.get(et, 0)
        pr = pr_ext.get(et, 0)
        pct = ch / (ch + pr) * 100 if (ch + pr) > 0 else 0
        print(f"{str(et):<20} {ch:>10} {pr:>10} {pct:>9.1f}%")

# c7 change analysis
print("\n--- 7-cycle count changes ---")
if changing_data:
    deltas_c7 = Counter((r[5], r[6], r[2]) for r in changing_data)
    print("  (c7_before, c7_after, delta_H): count")
    for key in sorted(deltas_c7.keys()):
        print(f"  {key}: {deltas_c7[key]}")

print("\n" + "=" * 70)
print("PART 2: THREADING ANALYSIS — How 7-cycles pass through the 4-subset")
print("=" * 70)

# For a small number of H-changing pairs, trace EXACTLY which 7-cycles exist
# before and after the reversal

def find_all_ham_cycles(A, n):
    """Find all directed Hamiltonian cycles (fixing vertex 0 as start)."""
    cycles = []
    for perm in permutations(range(1, n)):
        path = (0,) + perm
        valid = True
        for i in range(n):
            if A[path[i]][path[(i+1) % n]] != 1:
                valid = False
                break
        if valid:
            cycles.append(path)
    return cycles

if changing_data:
    print("\nDetailed analysis of first 5 H-changing reversals:")
    for idx, (bits, subset, dH, ext_pat, ext_type, c7b, c7a) in enumerate(changing_data[:5]):
        A = bits_to_adj(bits, n)
        B = reverse_subtournament(A, n, list(subset))

        cycles_before = find_all_ham_cycles(A, n)
        cycles_after = find_all_ham_cycles(B, n)

        ext = [v for v in range(n) if v not in subset]

        print(f"\n  Example {idx+1}: subset={subset}, ext={ext}, dH={dH}")
        print(f"    7-cycles before: {len(cycles_before)}, after: {len(cycles_after)}")

        # For each cycle, record how it passes through the subset
        # A "threading" is the sequence of subset/ext membership: e.g., SSSEESS
        def threading(cycle, subset_set):
            return tuple('S' if v in subset_set else 'E' for v in cycle)

        subset_set = set(subset)

        # Count threadings
        thr_before = Counter()
        for c in cycles_before:
            t = threading(c, subset_set)
            # Normalize: circular, so find canonical rotation
            min_t = min(t[i:] + t[:i] for i in range(len(t)))
            thr_before[min_t] += 1

        thr_after = Counter()
        for c in cycles_after:
            t = threading(c, subset_set)
            min_t = min(t[i:] + t[:i] for i in range(len(t)))
            thr_after[min_t] += 1

        all_threadings = set(list(thr_before.keys()) + list(thr_after.keys()))
        print(f"    Threading patterns (S=subset, E=external):")
        for t in sorted(all_threadings):
            b = thr_before.get(t, 0)
            a = thr_after.get(t, 0)
            change = "  CHANGED" if b != a else ""
            print(f"      {''.join(t)}: {b} -> {a}{change}")

        # Also look at the arcs BETWEEN subset and external
        print(f"    Arcs subset->ext: ", end="")
        for s in subset:
            for e in ext:
                if A[s][e]:
                    print(f"{s}->{e}", end=" ")
        print()
        print(f"    Arcs ext->subset: ", end="")
        for e in ext:
            for s in subset:
                if A[e][s]:
                    print(f"{e}->{s}", end=" ")
        print()

print("\n" + "=" * 70)
print("PART 3: WHAT DISTINGUISHES THE 28.6% — LOCAL PREDICTABILITY")
print("=" * 70)

# For each lambda-preserving (1,1,2,2) reversal, compute local features
# and see if H-change is predictable

if changing_data or preserving_data:
    # Features to test:
    # 1. Number of arcs from subset to ext vs ext to subset
    # 2. Score sequence of ext vertices (in full tournament)
    # 3. Number of 3-cycles containing exactly 2 subset + 1 ext vertex
    # 4. Number of 5-cycles containing all 4 subset vertices + 1 ext
    # 5. Number of 5-cycles containing exactly 3 subset + 2 ext
    # 6. Whether the 4-subset has a source or sink

    # Let's compute simpler features on a fresh sample
    print("\nComputing predictive features on fresh sample...")

    feature_records = []  # (features_dict, is_changing)

    np.random.seed(123)
    for trial in range(1500):
        bits = np.random.randint(0, 1 << total_bits)
        A = bits_to_adj(bits, n)
        H_orig = count_ham_paths(A, n)
        lam_orig = lambda_graph(A, n)

        for subset in combinations(range(n), 4):
            ss = sub_scores(A, n, list(subset))
            if ss != (1, 1, 2, 2):
                continue

            B = reverse_subtournament(A, n, list(subset))
            lam_new = lambda_graph(B, n)

            if not np.array_equal(lam_orig, lam_new):
                continue

            H_new = count_ham_paths(B, n)
            is_changing = (H_new != H_orig)

            ext = [v for v in range(n) if v not in subset]
            sub_list = list(subset)

            # Feature 1: arc balance (subset->ext minus ext->subset)
            out_arcs = sum(A[s][e] for s in sub_list for e in ext)
            in_arcs = sum(A[e][s] for e in ext for s in sub_list)
            arc_balance = out_arcs - in_arcs  # Always 0? (4*3=12 arcs total, out+in=12)

            # Feature 2: full scores of subset vertices
            full_scores = tuple(sorted([sum(A[s][j] for j in range(n) if j != s) for s in sub_list]))

            # Feature 3: number of 3-cycles with exactly 2 from subset + 1 from ext
            mixed_c3 = 0
            for s1, s2 in combinations(sub_list, 2):
                for e in ext:
                    if (A[s1][s2] and A[s2][e] and A[e][s1]) or (A[s2][s1] and A[s1][e] and A[e][s2]):
                        mixed_c3 += 1

            # Feature 4: c3 count within ext (0 or 1 for 3 vertices)
            ext_c3 = count_cycles_on_subset(A, n, ext)

            # Feature 5: whether any subset vertex is a source/sink in full tournament
            has_source_sink = any(
                sum(A[s][j] for j in range(n) if j != s) in (0, n-1)
                for s in sub_list
            )

            # Feature 6: arc pattern signature
            arc_sig = tuple(sorted(
                (sum(A[s][e] for e in ext), sum(A[e][s] for e in ext))
                for s in sub_list
            ))

            # Feature 7: how many ext vertices beat ALL of low-score subset vertices
            low_score_verts = [s for s in sub_list
                              if sum(A[s][j] for j in sub_list if j != s) <= 1]
            high_score_verts = [s for s in sub_list
                               if sum(A[s][j] for j in sub_list if j != s) >= 2]

            # Feature 8: ext-to-subset score variance
            ext_sub_scores = [sum(A[e][s] for s in sub_list) for e in ext]
            ext_sub_var = np.var(ext_sub_scores) if ext_sub_scores else 0

            feature_records.append({
                'arc_balance': arc_balance,
                'full_scores': full_scores,
                'mixed_c3': mixed_c3,
                'ext_c3': ext_c3,
                'has_source_sink': has_source_sink,
                'arc_sig': arc_sig,
                'ext_sub_var': ext_sub_var,
                'out_arcs': out_arcs,
                'is_changing': is_changing,
                'delta_H': H_new - H_orig
            })

        if trial % 300 == 0 and trial > 0:
            print(f"  ... {trial} processed, {len(feature_records)} valid reversals")

    ch_recs = [r for r in feature_records if r['is_changing']]
    pr_recs = [r for r in feature_records if not r['is_changing']]
    print(f"\nTotal valid: {len(feature_records)}, changing: {len(ch_recs)}, preserving: {len(pr_recs)}")

    # Analyze each feature
    print("\n--- Feature: arc_balance (out-in) ---")
    ch_bal = Counter(r['arc_balance'] for r in ch_recs)
    pr_bal = Counter(r['arc_balance'] for r in pr_recs)
    for b in sorted(set(list(ch_bal.keys()) + list(pr_bal.keys()))):
        c = ch_bal.get(b, 0)
        p = pr_bal.get(b, 0)
        pct = c/(c+p)*100 if c+p > 0 else 0
        print(f"  balance={b}: changing={c}, preserving={p}, rate={pct:.1f}%")

    print("\n--- Feature: full_scores of 4-subset ---")
    ch_fs = Counter(r['full_scores'] for r in ch_recs)
    pr_fs = Counter(r['full_scores'] for r in pr_recs)
    for fs in sorted(set(list(ch_fs.keys()) + list(pr_fs.keys()))):
        c = ch_fs.get(fs, 0)
        p = pr_fs.get(fs, 0)
        pct = c/(c+p)*100 if c+p > 0 else 0
        print(f"  scores={fs}: changing={c}, preserving={p}, rate={pct:.1f}%")

    print("\n--- Feature: mixed_c3 (3-cycles with 2 subset + 1 ext) ---")
    ch_mc = Counter(r['mixed_c3'] for r in ch_recs)
    pr_mc = Counter(r['mixed_c3'] for r in pr_recs)
    for mc in sorted(set(list(ch_mc.keys()) + list(pr_mc.keys()))):
        c = ch_mc.get(mc, 0)
        p = pr_mc.get(mc, 0)
        pct = c/(c+p)*100 if c+p > 0 else 0
        print(f"  mixed_c3={mc}: changing={c}, preserving={p}, rate={pct:.1f}%")

    print("\n--- Feature: ext_c3 (3-cycle in external vertices) ---")
    ch_ec = Counter(r['ext_c3'] for r in ch_recs)
    pr_ec = Counter(r['ext_c3'] for r in pr_recs)
    for ec in sorted(set(list(ch_ec.keys()) + list(pr_ec.keys()))):
        c = ch_ec.get(ec, 0)
        p = pr_ec.get(ec, 0)
        pct = c/(c+p)*100 if c+p > 0 else 0
        print(f"  ext_c3={ec}: changing={c}, preserving={p}, rate={pct:.1f}%")

    print("\n--- Feature: arc_sig (subset vertex (out_to_ext, in_from_ext) sorted) ---")
    ch_as = Counter(r['arc_sig'] for r in ch_recs)
    pr_as = Counter(r['arc_sig'] for r in pr_recs)
    for sig in sorted(set(list(ch_as.keys()) + list(pr_as.keys()))):
        c = ch_as.get(sig, 0)
        p = pr_as.get(sig, 0)
        pct = c/(c+p)*100 if c+p > 0 else 0
        print(f"  arc_sig={sig}: changing={c}, preserving={p}, rate={pct:.1f}%")

    print("\n--- Feature: out_arcs (total arcs from subset to ext) ---")
    ch_oa = Counter(r['out_arcs'] for r in ch_recs)
    pr_oa = Counter(r['out_arcs'] for r in pr_recs)
    for oa in sorted(set(list(ch_oa.keys()) + list(pr_oa.keys()))):
        c = ch_oa.get(oa, 0)
        p = pr_oa.get(oa, 0)
        pct = c/(c+p)*100 if c+p > 0 else 0
        print(f"  out_arcs={oa}: changing={c}, preserving={p}, rate={pct:.1f}%")

print("\n" + "=" * 70)
print("PART 4: THE HIDDEN DIMENSION — {2,1,0} OVERLAP WEIGHTS AND VITALI")
print("=" * 70)

# Connection to overlap weights: W=2 (share 2 vertices), W=1 (share 1), W=0 (disjoint)
# Under the Vitali reversal, the 3-cycle VERTEX SETS change but overlap distribution preserves
# The question: does the overlap structure between 3-cycles and 7-cycles explain the phase transition?

# For a few examples, compute:
# - The 3-cycles that exist before/after
# - Their overlap graph before/after
# - The 7-cycle that appears/disappears
# - How it "threads" through the overlap structure

if changing_data:
    print("\nOverlap weight analysis for H-changing examples:")

    for idx in range(min(3, len(changing_data))):
        bits, subset, dH, ext_pat, ext_type, c7b, c7a = changing_data[idx]
        A = bits_to_adj(bits, n)
        B = reverse_subtournament(A, n, list(subset))

        # Find all 3-cycles before and after
        def find_3_cycles(M, n):
            cycles = []
            for combo in combinations(range(n), 3):
                i, j, k = combo
                if (M[i][j] and M[j][k] and M[k][i]) or (M[j][i] and M[i][k] and M[k][j]):
                    cycles.append(frozenset(combo))
            return cycles

        c3_before = find_3_cycles(A, n)
        c3_after = find_3_cycles(B, n)

        c3_only_before = set(c3_before) - set(c3_after)
        c3_only_after = set(c3_after) - set(c3_before)
        c3_common = set(c3_before) & set(c3_after)

        print(f"\n  Example {idx+1}: subset={subset}, dH={dH}")
        print(f"    3-cycles: {len(c3_before)} -> {len(c3_after)} (net change: {len(c3_after)-len(c3_before)})")
        print(f"    Lost: {[set(c) for c in c3_only_before]}")
        print(f"    Gained: {[set(c) for c in c3_only_after]}")
        print(f"    Common: {len(c3_common)}")

        # Overlap weights between changed cycles and common cycles
        if c3_only_before or c3_only_after:
            print(f"    Overlap of lost cycles with common cycles:")
            for lost in c3_only_before:
                overlaps = Counter()
                for common in c3_common:
                    w = len(lost & common)
                    overlaps[w] += 1
                print(f"      {set(lost)}: W=0:{overlaps[0]}, W=1:{overlaps[1]}, W=2:{overlaps[2]}")

            print(f"    Overlap of gained cycles with common cycles:")
            for gained in c3_only_after:
                overlaps = Counter()
                for common in c3_common:
                    w = len(gained & common)
                    overlaps[w] += 1
                print(f"      {set(gained)}: W=0:{overlaps[0]}, W=1:{overlaps[1]}, W=2:{overlaps[2]}")

        # Now look at 5-cycles
        def find_5_cycles(M, n):
            cycles = set()
            for combo in combinations(range(n), 5):
                for perm in permutations(combo[1:]):
                    path = (combo[0],) + perm
                    valid = True
                    for k in range(5):
                        if M[path[k]][path[(k+1) % 5]] != 1:
                            valid = False
                            break
                    if valid:
                        # Canonical form
                        min_idx = path.index(min(path))
                        canonical = path[min_idx:] + path[:min_idx]
                        cycles.add(canonical)
            return cycles

        c5_before = find_5_cycles(A, n)
        c5_after = find_5_cycles(B, n)

        c5_only_before = c5_before - c5_after
        c5_only_after = c5_after - c5_before

        print(f"    5-cycles: {len(c5_before)} -> {len(c5_after)} (net change: {len(c5_after)-len(c5_before)})")
        print(f"    Lost: {len(c5_only_before)}, Gained: {len(c5_only_after)}")

        # 7-cycles
        c7_before_list = find_all_ham_cycles(A, n)
        c7_after_list = find_all_ham_cycles(B, n)

        print(f"    7-cycles: {len(c7_before_list)} -> {len(c7_after_list)} (delta={len(c7_after_list)-len(c7_before_list)})")

        # Detail of lost/gained 7-cycles
        set_before = set(c7_before_list)
        set_after = set(c7_after_list)
        lost_7 = set_before - set_after
        gained_7 = set_after - set_before

        subset_set = set(subset)
        if lost_7:
            print(f"    Lost 7-cycles: {len(lost_7)}")
            for c in list(lost_7)[:3]:
                # Show which arcs within subset are used
                internal_arcs = []
                for i in range(7):
                    if c[i] in subset_set and c[(i+1)%7] in subset_set:
                        internal_arcs.append(f"{c[i]}->{c[(i+1)%7]}")
                print(f"      {c}, internal arcs: {internal_arcs}")

        if gained_7:
            print(f"    Gained 7-cycles: {len(gained_7)}")
            for c in list(gained_7)[:3]:
                internal_arcs = []
                for i in range(7):
                    if c[i] in subset_set and c[(i+1)%7] in subset_set:
                        internal_arcs.append(f"{c[i]}->{c[(i+1)%7]}")
                print(f"      {c}, internal arcs: {internal_arcs}")

print("\n" + "=" * 70)
print("PART 5: n=6 IMPOSSIBILITY — WHY c5 and c3 ALONE DETERMINE H")
print("=" * 70)

# At n=6, c7 doesn't exist (only odd cycles up to c5 matter for H via OCF)
# H = 1 + 2*alpha_1 + 4*alpha_2 where alpha_1 = c3+c5, alpha_2 = disjoint pairs
# Under lambda-preserving reversal: c3 preserved (always), c5 preserved (proved exhaustive at n=6)
# => alpha_1 preserved => if alpha_2 also preserved => H preserved
#
# But WHY is c5 preserved at n=6? Let's understand the mechanism.
# At n=6, a 5-cycle uses 5 of 6 vertices. The 4-vertex reversal affects 4 of 6.
# So a 5-cycle always has at least 3 of its vertices in the 4-subset.

print("\nExhaustive analysis of why c5 is invariant under lambda-preserving 4-reversal at n=6:")

n6 = 6
total_bits_6 = n6 * (n6-1) // 2
c5_changed = 0
c5_preserved = 0

for bits in range(min(1 << total_bits_6, 1 << total_bits_6)):  # exhaustive at n=6
    A = bits_to_adj(bits, n6)
    lam_orig = lambda_graph(A, n6)

    for subset in combinations(range(n6), 4):
        ss = sub_scores(A, n6, list(subset))
        if ss != (1, 1, 2, 2):
            continue

        B = reverse_subtournament(A, n6, list(subset))
        lam_new = lambda_graph(B, n6)

        if not np.array_equal(lam_orig, lam_new):
            continue

        c5_orig = count_directed_cycles(A, n6, 5)
        c5_new = count_directed_cycles(B, n6, 5)

        if c5_orig != c5_new:
            c5_changed += 1
        else:
            c5_preserved += 1

print(f"  Lambda-preserving (1,1,2,2) 4-reversals at n=6: {c5_changed + c5_preserved}")
print(f"  c5 changed: {c5_changed}")
print(f"  c5 preserved: {c5_preserved}")
print(f"  => c5 invariance at n=6: {'PROVED (exhaustive)' if c5_changed == 0 else 'FAILS'}")

# Now at n=7: how many 5-vertex subsets containing 3+ members of the 4-subset exist?
# And do their c5 counts change?
print("\n  At n=7: a 5-cycle can intersect the 4-subset in {3, 4} vertices")
print("  and the 3-cycle complement in {1, 2} vertices.")
print("  There are C(4,3)*C(3,2) + C(4,4)*C(3,1) = 4*3 + 1*3 = 15 such 5-vertex subsets")
print("  (all 21 5-subsets intersect the 4-subset in at least 2 vertices)")
print("  But only those with 3+ can have their internal arcs changed.")

print("\n" + "=" * 70)
print("PART 6: THE DIMENSIONAL EXPLANATION")
print("=" * 70)
print("""
The {2,1,0} overlap weight structure connects to a hidden higher-dimensional
structure in tournaments. Think of it this way:

- At n=5: The 4-vertex subset spans "almost everything" (4 of 5 vertices).
  Only 1 external vertex. The overlap weight space is 1-dimensional.
  The reversal is gauge-trivial because there's no room for non-trivial
  interaction between the 4-subset and the exterior.

- At n=6: The 4-vertex subset leaves 2 external vertices.
  The overlap weight space is 2-dimensional.
  But c5 is still preserved because every 5-cycle must use 4 of 6 vertices,
  and the lambda graph constrains the 5-cycle structure completely.

- At n=7: The 4-vertex subset leaves 3 external vertices.
  The overlap weight space is 3-dimensional.
  Now 7-cycles exist that THREAD through both the 4-subset AND all 3 external
  vertices. The lambda graph constrains c3 and c5, but the 7-cycle threading
  has a degree of freedom that lambda doesn't capture.

This is the "hidden higher-dimensional structure buried in tournaments of n 2,1,0":
the overlap weights W=2,1,0 describe a filtration of cycle interactions,
and the Vitali atom lives in the kernel of the W=2 and W=1 constraints
but is detected by the W=0 (disjoint) structure at the 7-cycle level.
""")

print("\n" + "=" * 70)
print("PART 7: COHOMOLOGICAL INTERPRETATION — IS THE VITALI ATOM A COBOUNDARY?")
print("=" * 70)

# The key structure: we have a "gauge transformation" (the 4-reversal) that
# preserves a "field" (the lambda graph) but may or may not preserve an "action" (H).
# This is reminiscent of a cohomological obstruction.
#
# Let's test: is there a cochain complex where:
# - 0-cochains = vertex labels
# - 1-cochains = arc orientations
# - 2-cochains = 3-cycle parities
# and the Vitali atom is a 2-cocycle that's a coboundary iff H is preserved?

# First: compute H1 of the lambda graph (as a simplicial complex)
# The lambda graph is an undirected weighted graph on n vertices
# Its clique complex has H0, H1, H2...

# Actually, let's check something more concrete:
# Is the H-change related to the NUMBER OF 7-CYCLES passing through ALL 4 subset vertices?

print("\nIs delta_H determined by c7 threading through all 4 subset vertices?")

if changing_data:
    for idx in range(min(5, len(changing_data))):
        bits, subset, dH, _, _, _, _ = changing_data[idx]
        A = bits_to_adj(bits, n)
        B = reverse_subtournament(A, n, list(subset))

        # Count 7-cycles using all 4 subset vertices
        cycles_A = find_all_ham_cycles(A, n)
        cycles_B = find_all_ham_cycles(B, n)

        subset_set = set(subset)
        c7_through_all_A = sum(1 for c in cycles_A if subset_set.issubset(set(c)))
        c7_through_all_B = sum(1 for c in cycles_B if subset_set.issubset(set(c)))

        # (Every 7-cycle uses all 7 vertices, so this is always = total c7)
        print(f"  Example {idx+1}: c7_A={len(cycles_A)}, c7_B={len(cycles_B)}, "
              f"through_all_A={c7_through_all_A}, through_all_B={c7_through_all_B}")
        print(f"    (All 7-cycles pass through all 4 subset vertices trivially: n=7)")

# Since n=7, every 7-cycle uses ALL vertices. So the question is really:
# which arcs within the 4-subset does each 7-cycle use?
# The reversal flips ALL internal arcs. So we need to count:
# 7-cycles that use an even vs odd number of internal arcs

print("\n  Key insight: since every 7-cycle uses all 7 vertices,")
print("  the reversal changes a 7-cycle iff it uses an ODD number of internal arcs.")
print("  (An arc from subset vertex i to subset vertex j is reversed; if the cycle")
print("   uses an even number of such arcs, the reversal cancels out.)")

print("\n  Counting 7-cycles by number of INTERNAL arcs used:")

if changing_data:
    for idx in range(min(3, len(changing_data))):
        bits, subset, dH, _, _, _, _ = changing_data[idx]
        A = bits_to_adj(bits, n)

        cycles = find_all_ham_cycles(A, n)
        subset_set = set(subset)

        by_internal = Counter()
        for c in cycles:
            internal = 0
            for i in range(7):
                if c[i] in subset_set and c[(i+1)%7] in subset_set:
                    internal += 1
            by_internal[internal] += 1

        print(f"\n  Example {idx+1} (dH={dH}): {len(cycles)} total 7-cycles")
        for k in sorted(by_internal.keys()):
            print(f"    {k} internal arcs: {by_internal[k]} cycles")

        # Now same for the reversed tournament
        B = reverse_subtournament(A, n, list(subset))
        cycles_B = find_all_ham_cycles(B, n)

        by_internal_B = Counter()
        for c in cycles_B:
            internal = 0
            for i in range(7):
                if c[i] in subset_set and c[(i+1)%7] in subset_set:
                    internal += 1
            by_internal_B[internal] += 1

        print(f"  After reversal: {len(cycles_B)} total 7-cycles")
        for k in sorted(by_internal_B.keys()):
            print(f"    {k} internal arcs: {by_internal_B[k]} cycles")

print("\n" + "=" * 70)
print("PART 8: THE EXACT MECHANISM — WHICH INTERNAL ARC CONFIGURATIONS SURVIVE?")
print("=" * 70)

# A 7-cycle visits all 7 vertices. The 4 subset vertices are visited in some ORDER.
# This order is a permutation of the 4 subset vertices within the cycle.
# The key question: for a given visit order, how many internal arcs does the cycle use?
#
# If subset = {a,b,c,d}, and the cycle visits them in order a...b...c...d...
# then the internal arcs are those where consecutive cycle vertices are BOTH in subset.
# The 4 subset vertices partition the cycle into 4 arcs (segments), and internal arcs
# are the "zero-length" segments where two subset vertices are adjacent in the cycle.

# A 7-cycle has 7 positions. 4 are for subset, 3 for external.
# The number of internal arcs = number of adjacent subset vertex pairs in the cycle.
# This equals 4 - (number of "runs" of external vertices between consecutive subset vertices)
# If all 3 external vertices are between different subset vertex pairs: 0 internal arcs
# If 2 external in one gap, 1 in another: 1 internal arc
# etc.

# Possible distributions of 3 external vertices into 4 gaps:
# (3,0,0,0): 3 in one gap, others empty -> 3 internal arcs
# (2,1,0,0): -> 2 internal arcs
# (1,1,1,0): -> 1 internal arc
# Can't have (1,1,1,1) since only 3 external and 4 gaps -> impossible
# So possible internal arc counts: 1, 2, or 3

print("\nPossible internal arc counts in a 7-cycle with 4 subset + 3 external vertices:")
print("  Distribution of 3 ext in 4 gaps -> internal arcs:")
print("  (3,0,0,0) -> 3 internal arcs")
print("  (2,1,0,0) -> 2 internal arcs")
print("  (1,1,1,0) -> 1 internal arc")
print("  (cannot fill all 4 gaps with only 3 ext vertices)")
print()
print("  So internal arc count is in {1, 2, 3}.")
print("  Reversal flips all internal arcs.")
print("  A cycle SURVIVES reversal iff ALL its internal arcs are in BOTH directions")
print("  in the original tournament (impossible for tournament: each arc has one direction).")
print("  So NO cycle survives the reversal in its original form.")
print("  Instead, cycles map to OTHER potential cycles.")
print()
print("  The delta_c7 = c7(B) - c7(A) counts the NET change in valid 7-cycles.")

# Let's verify and deepen this
print("\nVerifying with concrete example:")
if changing_data:
    bits, subset, dH, _, _, _, _ = changing_data[0]
    A = bits_to_adj(bits, n)
    B = reverse_subtournament(A, n, list(subset))

    cycles_A = find_all_ham_cycles(A, n)
    cycles_B = find_all_ham_cycles(B, n)

    subset_set = set(subset)

    # For each cycle in A, try to apply reversal and see if result is in B
    print(f"  Tournament A: {len(cycles_A)} 7-cycles, B: {len(cycles_B)} 7-cycles")

    # A cycle c = (v0, v1, ..., v6, v0) is valid iff A[vi][v_{i+1}]=1 for all i
    # After reversal on subset S, arc (u,v) becomes:
    #   A[v][u] if u,v both in S (reversed)
    #   A[u][v] if at most one of u,v in S (unchanged)
    # So c is valid in B iff:
    #   For internal arcs: A[v_{i+1}][vi] = 1 (reversed direction works)
    #   For external arcs: A[vi][v_{i+1}] = 1 (unchanged, same as before)

    # Count cycles that go from valid to valid/invalid
    survived = 0
    lost = 0
    for c in cycles_A:
        valid_in_B = True
        for i in range(7):
            u, v = c[i], c[(i+1)%7]
            if u in subset_set and v in subset_set:
                # Internal arc: need B[u][v]=1, which is A[v][u]
                if A[v][u] != 1:
                    valid_in_B = False
                    break
            else:
                # External arc: unchanged
                if A[u][v] != 1:
                    valid_in_B = False
                    break
        if valid_in_B:
            survived += 1
        else:
            lost += 1

    print(f"  Of A's {len(cycles_A)} cycles: {survived} survive in B, {lost} lost")

    # Similarly check cycles in B that aren't from A
    gained = len(cycles_B) - survived
    print(f"  Gained in B: {gained}")
    print(f"  Net: {gained - lost} = delta_c7")

print("\n" + "=" * 70)
print("PART 9: INTERNAL ARC PARITY AND THE COBOUNDARY MAP")
print("=" * 70)

# The key structural question: when does a lambda-preserving (1,1,2,2) reversal
# change the 7-cycle count by exactly +1 or -1?
#
# Let's think about this combinatorially. The (1,1,2,2) sub-tournament has
# exactly 2 directed 3-cycles. After reversal, it STILL has 2 directed 3-cycles
# (reversal preserves c3 count within the sub-tournament).
#
# The 4-vertex sub-tournament has 6 arcs. Reversal flips all 6.
# A 7-cycle uses k of these internal arcs (k in {1,2,3}).
# After reversal, the cycle's k internal arcs are all flipped.
# The cycle survives iff the reversed internal arcs still form valid edges,
# which they do iff the ORIGINAL cycle had those arcs in the "wrong" direction.
#
# Wait, that's not quite right. Let me think more carefully.
#
# After reversal, a 7-cycle c is valid in B iff:
# - For each arc (u,v) in c with both u,v in S: A[v][u] = 1 (i.e., original arc goes v->u)
# - But if c has arc u->v and A[u][v]=1 (it's valid in A), then A[v][u]=0.
# So a cycle that uses an internal arc u->v in A CANNOT survive in B
# (because B has arc v->u there, but the cycle needs u->v).
#
# Conversely, a cycle c' that is valid in B but not A must have:
# - Some internal arcs (u,v) where B[u][v]=1 but A[u][v]=0
# - i.e., A[v][u]=1 (the arc goes v->u in A, but c' needs u->v)
#
# So the mapping is: c in A with internal arcs -> flip those arcs -> check if still valid
# This is NOT the same cycle; it's a different sequence of vertices.

# Let me just directly compute: for each 7-cycle in A, reverse its internal arcs
# and check if the resulting sequence is a valid 7-cycle in B.

print("Direct reversal mapping of 7-cycles:")

if changing_data:
    for idx in range(min(3, len(changing_data))):
        bits, subset, dH, _, _, _, _ = changing_data[idx]
        A = bits_to_adj(bits, n)
        B = reverse_subtournament(A, n, list(subset))

        cycles_A = find_all_ham_cycles(A, n)
        cycles_B = find_all_ham_cycles(B, n)

        subset_set = set(subset)

        # For each cycle in A, classify by internal arc pattern
        print(f"\n  Example {idx+1} (dH={dH}):")

        # The cycle visits subset vertices in some order. Let's track this.
        for c in cycles_A[:5]:  # first 5 cycles
            internal_arcs = []
            external_arcs = []
            for i in range(7):
                u, v = c[i], c[(i+1)%7]
                if u in subset_set and v in subset_set:
                    internal_arcs.append((u,v))
                else:
                    external_arcs.append((u,v))

            # Order of subset visits
            sub_order = [v for v in c if v in subset_set]

            print(f"    Cycle {c}: sub_order={sub_order}, "
                  f"internal={internal_arcs}, #ext={len(external_arcs)}")

print("\n" + "=" * 70)
print("PART 10: THE DEFINITIVE TEST — WHAT PREDICTS H-CHANGE?")
print("=" * 70)

# Based on all the above, the most promising predictor is the NET difference
# in 7-cycle count under reversal. But we want a PREDICTOR from local data.
#
# Hypothesis: H changes iff the number of 7-cycles with an ODD number of
# internal arcs differs from those with an EVEN number.
# Wait — all cycles have 1, 2, or 3 internal arcs. Even = {2}, Odd = {1, 3}.
#
# Let's compute: for each case, count c7_odd (1 or 3 internal) and c7_even (2 internal)
# H-change might relate to c7_odd vs c7_even balance.

print("\nTesting: is delta_c7 = c7_even - c7_odd ?")

if changing_data:
    correct = 0
    total_test = 0

    for data_list, label in [(changing_data, "CHANGING"), (preserving_data, "PRESERVING")]:
        for idx in range(min(20, len(data_list))):
            bits, subset, dH, _, _, _, _ = data_list[idx]
            A = bits_to_adj(bits, n)
            B = reverse_subtournament(A, n, list(subset))

            cycles_A = find_all_ham_cycles(A, n)
            subset_set = set(subset)

            c7_odd = 0
            c7_even = 0
            for c in cycles_A:
                internal = sum(1 for i in range(7) if c[i] in subset_set and c[(i+1)%7] in subset_set)
                if internal % 2 == 0:
                    c7_even += 1
                else:
                    c7_odd += 1

            c7_A = len(cycles_A)
            c7_B = len(find_all_ham_cycles(B, n))
            delta = c7_B - c7_A
            predicted = c7_even - c7_odd

            total_test += 1
            if delta == predicted:
                correct += 1

    print(f"  Correct predictions: {correct}/{total_test}")

# Try another formula: delta_c7 related to internal arc structure
print("\nTesting: is delta_c7 = 2*(c7_even - c7_odd)?")
if changing_data:
    correct2 = 0
    total_test2 = 0
    deltas_check = []

    for data_list, label in [(changing_data, "CHANGING"), (preserving_data, "PRESERVING")]:
        for idx in range(min(30, len(data_list))):
            bits, subset, dH, _, _, _, _ = data_list[idx]
            A = bits_to_adj(bits, n)
            B = reverse_subtournament(A, n, list(subset))

            cycles_A = find_all_ham_cycles(A, n)
            subset_set = set(subset)

            c7_1 = sum(1 for c in cycles_A
                       if sum(1 for i in range(7)
                              if c[i] in subset_set and c[(i+1)%7] in subset_set) == 1)
            c7_2 = sum(1 for c in cycles_A
                       if sum(1 for i in range(7)
                              if c[i] in subset_set and c[(i+1)%7] in subset_set) == 2)
            c7_3 = sum(1 for c in cycles_A
                       if sum(1 for i in range(7)
                              if c[i] in subset_set and c[(i+1)%7] in subset_set) == 3)

            c7_A = len(cycles_A)
            c7_B = len(find_all_ham_cycles(B, n))
            delta = c7_B - c7_A

            total_test2 += 1
            deltas_check.append({
                'label': label,
                'delta': delta,
                'c7_1': c7_1, 'c7_2': c7_2, 'c7_3': c7_3,
                'dH': dH
            })

            if delta == 2 * (c7_2 - c7_1 - c7_3):
                correct2 += 1

    print(f"  Correct: {correct2}/{total_test2}")

    # Just display the data to find the pattern
    print(f"\n  {'Label':>10} {'delta':>6} {'c7_1':>5} {'c7_2':>5} {'c7_3':>5} {'dH':>4}")
    for d in deltas_check:
        print(f"  {d['label']:>10} {d['delta']:>6} {d['c7_1']:>5} {d['c7_2']:>5} {d['c7_3']:>5} {d['dH']:>4}")

print("\n\nDone. See vitali_c7_mechanism.out for full results.")
