"""
alpha10_forcing_n7.py -- kind-pasteur-2026-03-14-S66

At n=7, WHY does alpha_1=10 force alpha_2>=2?

At n=6, the mechanism is rigid: score (1,2,2,3,3,4), c3=6, 4 directed
5-cycles, 2 complementary 3-cycle pairs. All alpha_1=10 tournaments identical.

At n=7, the picture is richer. We need to understand:
1. Which cycle size decompositions give alpha_1=10?
2. Why can't all 10 cycles be pairwise conflicting?
3. What's the minimum alpha_2, and why is it >= 2?

APPROACH: Sample n=7 tournaments, filter for alpha_1 near 10,
analyze cycle structure and disjoint pair mechanisms.
"""

import numpy as np
from itertools import combinations

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=np.int8)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def count_directed_hamcycles(A, vertices):
    k = len(vertices)
    if k < 3 or k % 2 == 0:
        return 0
    vlist = list(vertices)
    sub = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            if i != j:
                sub[i][j] = int(A[vlist[i]][vlist[j]])
    full = (1 << k) - 1
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if dp[mask][v] == 0:
                continue
            for u in range(1, k):
                if mask & (1 << u):
                    continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    total = 0
    for v in range(1, k):
        if dp[full][v] and sub[v][0]:
            total += dp[full][v]
    return total

def full_cycle_analysis(A, n):
    """Complete cycle analysis: all odd cycles, their directed counts, disjoint pairs."""
    cycles_by_size = {}  # size -> list of (vertex_set, directed_count)

    for size in range(3, n+1, 2):
        cycles_by_size[size] = []
        for subset in combinations(range(n), size):
            cnt = count_directed_hamcycles(A, list(subset))
            if cnt > 0:
                cycles_by_size[size].append((frozenset(subset), cnt))

    # Flatten all cycles
    all_cycles = []
    for size in sorted(cycles_by_size.keys()):
        for vs, cnt in cycles_by_size[size]:
            all_cycles.append((vs, cnt, size))

    alpha_1 = sum(cnt for _, cnt, _ in all_cycles)

    # Disjoint pairs (vertex-set level, weighted by directed cycle count)
    alpha_2 = 0
    disjoint_detail = []
    for i in range(len(all_cycles)):
        for j in range(i+1, len(all_cycles)):
            vs_i, cnt_i, sz_i = all_cycles[i]
            vs_j, cnt_j, sz_j = all_cycles[j]
            if len(vs_i & vs_j) == 0:
                alpha_2 += cnt_i * cnt_j
                disjoint_detail.append((sz_i, sz_j, vs_i, vs_j, cnt_i, cnt_j))

    return {
        'cycles_by_size': cycles_by_size,
        'all_cycles': all_cycles,
        'alpha_1': alpha_1,
        'alpha_2': alpha_2,
        'disjoint_detail': disjoint_detail
    }

def main():
    print("=" * 70)
    print("alpha_1=10 FORCING MECHANISM AT n=7")
    print("=" * 70)

    n = 7
    rng = np.random.default_rng(2026_03_14_10)
    num_samples = 100000

    # Collect alpha_1=10 tournaments
    a10_tours = []
    # Also collect alpha_1 near 10 for context
    a_near_10 = {a1: {'min_a2': float('inf'), 'count': 0} for a1 in range(7, 14)}

    for trial in range(num_samples):
        A = random_tournament(n, rng)

        # Fast pre-filter: count 3-cycles only first
        c3 = 0
        for a, b, c in combinations(range(n), 3):
            if (A[a][b] and A[b][c] and A[c][a]) or \
               (A[a][c] and A[c][b] and A[b][a]):
                c3 += 1

        # c3 alone is often >= 10, so alpha_1 >= c3 >= 10
        # We need alpha_1 EXACTLY 10, so c3 <= 10
        if c3 > 12:  # some slack for 5-cycle contribution
            continue

        info = full_cycle_analysis(A, n)
        a1 = info['alpha_1']
        a2 = info['alpha_2']

        if 7 <= a1 <= 13:
            if a1 in a_near_10:
                a_near_10[a1]['count'] += 1
                a_near_10[a1]['min_a2'] = min(a_near_10[a1]['min_a2'], a2)

        if a1 == 10:
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))
            size_decomp = {}
            for size, cycs in info['cycles_by_size'].items():
                dc = sum(cnt for _, cnt in cycs)
                if dc > 0:
                    size_decomp[size] = (len(cycs), dc)
            a10_tours.append({
                'scores': scores,
                'a2': a2,
                'size_decomp': size_decomp,
                'disjoint_detail': info['disjoint_detail'],
                'c3': c3
            })

        if (trial + 1) % 25000 == 0:
            print(f"  {trial+1}/{num_samples}, found {len(a10_tours)} with alpha_1=10")

    print(f"\nTotal alpha_1=10 found: {len(a10_tours)}")

    # Context: min alpha_2 for alpha_1 near 10
    print(f"\nMin alpha_2 for alpha_1 near 10:")
    for a1 in sorted(a_near_10.keys()):
        d = a_near_10[a1]
        if d['count'] > 0:
            print(f"  alpha_1={a1:2d}: min_a2={d['min_a2']}, "
                  f"count={d['count']}")

    if not a10_tours:
        print("\nNo alpha_1=10 found in sample. Increase sample size.")
        return

    # Analyze alpha_2 distribution
    a2_dist = {}
    for t in a10_tours:
        a2 = t['a2']
        a2_dist[a2] = a2_dist.get(a2, 0) + 1

    print(f"\nalpha_2 distribution for alpha_1=10:")
    for a2 in sorted(a2_dist.keys()):
        print(f"  alpha_2={a2}: {a2_dist[a2]} tournaments")

    # Analyze score sequences
    score_dist = {}
    for t in a10_tours:
        s = t['scores']
        score_dist[s] = score_dist.get(s, 0) + 1

    print(f"\nScore sequences for alpha_1=10:")
    for s in sorted(score_dist.keys(), key=lambda x: -score_dist[x])[:10]:
        c3 = 35 - sum(si*(si-1)//2 for si in s)  # C(7,3)=35
        print(f"  {s}: {score_dist[s]} tours, c3={c3}")

    # Analyze size decompositions
    decomp_dist = {}
    for t in a10_tours:
        key = tuple(sorted(t['size_decomp'].items()))
        decomp_dist[key] = decomp_dist.get(key, 0) + 1

    print(f"\nCycle size decompositions (size: (vertex_sets, directed_cycles)):")
    for key in sorted(decomp_dist.keys(), key=lambda x: -decomp_dist[x])[:10]:
        desc = ", ".join(f"sz{s}: {vs} sets/{dc} dir" for s, (vs, dc) in key)
        print(f"  {desc}: {decomp_dist[key]} tours")

    # Analyze disjoint pair types
    pair_type_dist = {}
    for t in a10_tours:
        pair_types = tuple(sorted([(s1, s2) for s1, s2, _, _, _, _ in t['disjoint_detail']]))
        pair_type_dist[pair_types] = pair_type_dist.get(pair_types, 0) + 1

    print(f"\nDisjoint pair types (cycle sizes):")
    for pt in sorted(pair_type_dist.keys(), key=lambda x: -pair_type_dist[x])[:10]:
        desc = ", ".join(f"{s1}x{s2}" for s1, s2 in pt) if pt else "(none)"
        # Get corresponding alpha_2
        print(f"  [{desc}]: {pair_type_dist[pt]} tours")

    # Deep dive: first few examples
    print(f"\n{'='*70}")
    print("DETAILED EXAMPLES")
    print(f"{'='*70}")

    for idx, t in enumerate(a10_tours[:5]):
        print(f"\nExample {idx+1}: scores={t['scores']}, c3={t['c3']}, alpha_2={t['a2']}")
        for size, (vs_count, dc) in sorted(t['size_decomp'].items()):
            print(f"  {size}-cycles: {vs_count} vertex sets, {dc} directed cycles")
        if t['disjoint_detail']:
            print(f"  Disjoint pairs:")
            for s1, s2, vs1, vs2, c1, c2 in t['disjoint_detail']:
                print(f"    {set(vs1)}({s1}-cyc,{c1}d) x {set(vs2)}({s2}-cyc,{c2}d) -> {c1*c2}")

    # Part 2: The alpha_1=3 comparison
    print(f"\n{'='*70}")
    print("COMPARISON: alpha_1=3 FORCING MECHANISM")
    print(f"{'='*70}")

    rng2 = np.random.default_rng(2026_03_14_03)
    a3_tours = []

    for trial in range(num_samples):
        A = random_tournament(n, rng2)
        c3 = 0
        for a, b, c in combinations(range(n), 3):
            if (A[a][b] and A[b][c] and A[c][a]) or \
               (A[a][c] and A[c][b] and A[b][a]):
                c3 += 1
        if c3 > 5:
            continue
        info = full_cycle_analysis(A, n)
        if info['alpha_1'] == 3:
            scores = tuple(sorted([sum(A[i]) for i in range(n)]))
            a3_tours.append({
                'scores': scores,
                'a2': info['alpha_2'],
                'c3': c3,
                'size_decomp': {sz: (len(cycs), sum(c for _, c in cycs))
                                for sz, cycs in info['cycles_by_size'].items()
                                if cycs},
                'disjoint_detail': info['disjoint_detail']
            })

        if (trial + 1) % 25000 == 0 and a3_tours:
            pass

    print(f"\nTotal alpha_1=3 found: {len(a3_tours)}")
    if a3_tours:
        a2_dist_3 = {}
        for t in a3_tours:
            a2 = t['a2']
            a2_dist_3[a2] = a2_dist_3.get(a2, 0) + 1

        print(f"alpha_2 distribution for alpha_1=3:")
        for a2 in sorted(a2_dist_3.keys()):
            print(f"  alpha_2={a2}: {a2_dist_3[a2]} tournaments")

        score_dist_3 = {}
        for t in a3_tours:
            score_dist_3[t['scores']] = score_dist_3.get(t['scores'], 0) + 1
        print(f"\nScore sequences for alpha_1=3:")
        for s in sorted(score_dist_3.keys(), key=lambda x: -score_dist_3[x])[:10]:
            c3 = 35 - sum(si*(si-1)//2 for si in s)
            print(f"  {s}: {score_dist_3[s]} tours, c3={c3}")

        print(f"\nDetailed examples:")
        for idx, t in enumerate(a3_tours[:3]):
            print(f"\nExample {idx+1}: scores={t['scores']}, c3={t['c3']}, alpha_2={t['a2']}")
            for size, (vs_count, dc) in sorted(t['size_decomp'].items()):
                print(f"  {size}-cycles: {vs_count} vertex sets, {dc} directed cycles")
            if t['disjoint_detail']:
                print(f"  Disjoint pairs:")
                for s1, s2, vs1, vs2, c1, c2 in t['disjoint_detail']:
                    print(f"    {set(vs1)}({s1}-cyc,{c1}d) x {set(vs2)}({s2}-cyc,{c2}d) -> {c1*c2}")

    # Part 3: The STRUCTURAL argument
    print(f"\n{'='*70}")
    print("STRUCTURAL ARGUMENT ATTEMPT")
    print(f"{'='*70}")
    print()
    print("For alpha_1=3:")
    print("  c3=3 at n>=7: three 3-cycles that pairwise share a vertex.")
    print("  Helly property: 3 pairwise intersecting sets of size 3 in universe of")
    print("  size >= 7 need NOT have a common element.")
    print("  Without common vertex: can arrange the 3 cycles to use 7 vertices,")
    print("  but then the tournament on those 7 vertices must have c3=3 exactly.")
    print("  Score seq constraint: c3 = C(7,3) - sum C(s_i,2) = 35 - sum C(s_i,2)")
    print("  For c3=3: sum C(s_i,2) = 32. Near-transitive (lots of score spread).")
    print("  But 3 non-common-vertex 3-cycles at n=7 almost always generate")
    print("  5-cycles and additional 3-cycles, pushing alpha_1 above 3.")
    print("  So alpha_1=3 requires very specific structure => forces alpha_2>=2.")
    print()
    print("For alpha_1=10:")
    print("  At n=6: RIGID. Only score (1,2,2,3,3,4), c3=6+4dc5=10.")
    print("  At n>=7: Need 10 directed odd cycles total (3+5+7-cycles).")
    print("  With alpha_2=0: ALL 10 cycles pairwise share a vertex.")
    print("  10 pairwise overlapping odd cycles on 7+ vertices:")
    print("  Each cycle uses >= 3 vertices. If all are 3-cycles, c3=10,")
    print("  and c3 = C(n,3) - sum C(s_i,2).")
    print("  At n=7: c3=10 means sum C(s_i,2)=25 out of C(7,3)=35.")
    print("  score (0,2,2,3,3,4,7-1=6)? -> sum = 0+1+1+3+3+6+?")
    print()
    # Compute which scores give c3=10 at n=7
    print("  Score sequences with c3=10 at n=7 (sum C(s_i,2) = 25):")
    from itertools import product as iprod
    count_valid = 0
    for scores in combinations(range(7), 7):
        pass  # This is wrong, scores are out-degrees summing to C(7,2)=21

    # Enumerate valid score sequences at n=7
    # Score s_0 <= s_1 <= ... <= s_6, sum = 21
    valid_scores_c3_10 = []
    def gen_scores(pos, remaining, lo, acc):
        if pos == 7:
            if remaining == 0:
                c3 = 35 - sum(s*(s-1)//2 for s in acc)
                if c3 == 10:
                    valid_scores_c3_10.append(tuple(acc))
            return
        hi = min(remaining, 6)
        for s in range(lo, hi+1):
            if remaining - s >= (6 - pos) * s:  # enough remaining for rest
                gen_scores(pos+1, remaining-s, s, acc + [s])

    gen_scores(0, 21, 0, [])
    for s in valid_scores_c3_10:
        print(f"    {s}: c3=10, sum C(s_i,2)={sum(si*(si-1)//2 for si in s)}")

if __name__ == "__main__":
    main()
