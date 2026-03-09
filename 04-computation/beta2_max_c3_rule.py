"""
Verify: vertex participating in MOST 3-cycles is always a good vertex
when b1(T) = 0.

This was 100% at n=6 exhaustive. Verify at n=7,8,9 sampled.
Also investigate WHY this works and attempt algebraic characterization.
"""

import numpy as np
from itertools import combinations
from collections import defaultdict

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

def find_3cycles(A, n):
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            cycles.append((a, b, c))
        if A[a][c] and A[c][b] and A[b][a]:
            cycles.append((a, c, b))
    return cycles

def shared_directed_edge(c1, c2):
    edges1 = {(c1[0],c1[1]), (c1[1],c1[2]), (c1[2],c1[0])}
    edges2 = {(c2[0],c2[1]), (c2[1],c2[2]), (c2[2],c2[0])}
    return len(edges1 & edges2) > 0

def cycle_graph_components(cycles):
    if not cycles:
        return []
    nc = len(cycles)
    adj = [[] for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if shared_directed_edge(cycles[i], cycles[j]):
                adj[i].append(j)
                adj[j].append(i)
    visited = [False]*nc
    components = []
    for start in range(nc):
        if visited[start]:
            continue
        comp = []
        stack = [start]
        while stack:
            v = stack.pop()
            if visited[v]:
                continue
            visited[v] = True
            comp.append(v)
            for u in adj[v]:
                if not visited[u]:
                    stack.append(u)
        components.append(comp)
    return components

def is_dominated(A, n, cyc):
    verts = set(cyc)
    for d in range(n):
        if d in verts:
            continue
        a, b, c = cyc
        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
            return True
    return False

def compute_b1(A, n):
    cycles = find_3cycles(A, n)
    if not cycles:
        return 0
    comps = cycle_graph_components(cycles)
    free_count = 0
    for comp in comps:
        if all(not is_dominated(A, n, cycles[ci]) for ci in comp):
            free_count += 1
    return free_count

def delete_vertex(A, n, v):
    keep = [i for i in range(n) if i != v]
    B = np.zeros((n-1, n-1), dtype=int)
    for i, ki in enumerate(keep):
        for j, kj in enumerate(keep):
            B[i][j] = A[ki][kj]
    return B

def count_3cycles_per_vertex(A, n):
    """Count directed 3-cycles involving each vertex."""
    counts = [0]*n
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            counts[a] += 1; counts[b] += 1; counts[c] += 1
        elif A[a][c] and A[c][b] and A[b][a]:
            counts[a] += 1; counts[b] += 1; counts[c] += 1
    return counts

def main():
    print("=" * 70)
    print("MAX 3-CYCLE PARTICIPATION RULE: Always good vertex?")
    print("=" * 70)

    for n in [5, 6, 7, 8, 9]:
        total = 2 ** (n*(n-1)//2)
        if n <= 6:
            sample_range = range(total)
            method = "exhaustive"
        else:
            rng = np.random.RandomState(42)
            nsamp = {7: 5000, 8: 2000, 9: 500}[n]
            sample_range = [rng.randint(0, total) for _ in range(nsamp)]
            method = f"sampled({nsamp})"

        n_b1_0 = 0
        max_c3_good = 0
        max_c3_good_strict = 0  # when there's a UNIQUE max
        any_max_c3_good = 0  # some vertex tied for max is good
        min_c3_bad_always = 0  # min c3 vertex is always bad?
        n_with_bad = 0

        # Additional rules to test
        rule_results = defaultdict(int)

        for bits in sample_range:
            A = bits_to_adj(bits, n)
            if compute_b1(A, n) != 0:
                continue
            n_b1_0 += 1

            c3 = count_3cycles_per_vertex(A, n)
            scores = [sum(A[i]) for i in range(n)]

            # Compute b1(T\v) for each v
            b1_tv = []
            for v in range(n):
                B = delete_vertex(A, n, v)
                b1_tv.append(compute_b1(B, n-1))

            bad = set(v for v in range(n) if b1_tv[v] > 0)
            if not bad:
                max_c3_good += 1
                max_c3_good_strict += 1
                any_max_c3_good += 1
                for rule in rule_results:
                    pass
                n_b1_0  # already counted
                continue

            n_with_bad += 1

            # Max c3 vertex
            max_c3 = max(c3)
            max_c3_verts = [v for v in range(n) if c3[v] == max_c3]

            if any(v not in bad for v in max_c3_verts):
                any_max_c3_good += 1
            if len(max_c3_verts) == 1 and max_c3_verts[0] not in bad:
                max_c3_good_strict += 1

            # Test: vertex with max score*(n-1-score) = max product
            products = [scores[v]*(n-1-scores[v]) for v in range(n)]
            max_prod_v = max(range(n), key=lambda v: products[v])
            if max_prod_v not in bad:
                rule_results['max_product'] += 1

            # Test: vertex with max (c3(v) - dominated_c3(v))
            # = vertex with most FREE 3-cycles through it
            free_c3 = [0]*n
            cycles = find_3cycles(A, n)
            for cyc in cycles:
                if not is_dominated(A, n, cyc):
                    for v in cyc:
                        free_c3[v] += 1
            max_free_c3_v = max(range(n), key=lambda v: free_c3[v])
            if max_free_c3_v not in bad:
                rule_results['max_free_c3'] += 1

            # Test: vertex in most DOMINATED cycles
            dom_c3 = [c3[v] - free_c3[v] for v in range(n)]
            max_dom_c3_v = max(range(n), key=lambda v: dom_c3[v])
            if max_dom_c3_v not in bad:
                rule_results['max_dom_c3'] += 1

            # Test: vertex with score closest to (n-1)/2 (most "central")
            central_dist = [abs(scores[v] - (n-1)/2) for v in range(n)]
            min_dist = min(central_dist)
            central_verts = [v for v in range(n) if central_dist[v] == min_dist]
            if any(v not in bad for v in central_verts):
                rule_results['most_central'] += 1

        print(f"\nn={n} ({method}): {n_b1_0} with b1=0, {n_with_bad} have bad vertices")
        print(f"  Any max-c3 vertex good:     {any_max_c3_good}/{n_b1_0} = {100*any_max_c3_good/n_b1_0:.1f}%")
        if n_with_bad > 0:
            for rule, count in sorted(rule_results.items()):
                print(f"  {rule} good: {count}/{n_with_bad} ({100*count/n_with_bad:.1f}%)")

    # Deep investigation: WHY does max-c3 vertex always work?
    print(f"\n{'='*70}")
    print("WHY DOES MAX-C3 WORK? Detailed analysis at n=6")
    print(f"{'='*70}")

    n = 6
    total = 2**(n*(n-1)//2)

    # For bad vertices: compare c3 count to good vertices
    bad_c3_gaps = []  # (c3_bad, c3_good_max, gap)
    bad_free_c3 = []

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        c3 = count_3cycles_per_vertex(A, n)
        b1_tv = []
        for v in range(n):
            B = delete_vertex(A, n, v)
            b1_tv.append(compute_b1(B, n-1))

        bad = [v for v in range(n) if b1_tv[v] > 0]
        good = [v for v in range(n) if b1_tv[v] == 0]

        if bad and good:
            max_bad_c3 = max(c3[v] for v in bad)
            min_good_c3 = min(c3[v] for v in good)
            gap = min_good_c3 - max_bad_c3
            bad_c3_gaps.append(gap)

    if bad_c3_gaps:
        from collections import Counter
        gap_dist = Counter(bad_c3_gaps)
        print(f"\nGap = min_good_c3 - max_bad_c3:")
        for gap in sorted(gap_dist.keys()):
            print(f"  gap={gap}: {gap_dist[gap]} tournaments")
        print(f"  Min gap: {min(bad_c3_gaps)}")
        print(f"  ALWAYS positive? {all(g > 0 for g in bad_c3_gaps)}")
        print(f"  ALWAYS non-negative? {all(g >= 0 for g in bad_c3_gaps)}")

    # Check: is c3 count for bad vertex strictly LESS than all good vertices?
    print(f"\n{'='*70}")
    print("STRICT SEPARATION: Bad vertices have strictly fewer 3-cycles?")
    print(f"{'='*70}")

    n = 6
    strict_sep = 0
    non_strict = 0
    total_with_both = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        c3 = count_3cycles_per_vertex(A, n)
        b1_tv = []
        for v in range(n):
            B = delete_vertex(A, n, v)
            b1_tv.append(compute_b1(B, n-1))

        bad = [v for v in range(n) if b1_tv[v] > 0]
        good = [v for v in range(n) if b1_tv[v] == 0]

        if bad and good:
            total_with_both += 1
            max_bad = max(c3[v] for v in bad)
            min_good = min(c3[v] for v in good)
            if max_bad < min_good:
                strict_sep += 1
            else:
                non_strict += 1
                if non_strict <= 3:
                    print(f"  NON-STRICT: bad c3={[c3[v] for v in bad]}, good c3={[c3[v] for v in good]}")
                    print(f"    scores: {[sum(A[i]) for i in range(n)]}")

    print(f"\nStrict separation: {strict_sep}/{total_with_both} ({100*strict_sep/total_with_both:.1f}%)")
    print(f"Non-strict: {non_strict}/{total_with_both}")

    # FORMULA for c3(v) in terms of score
    print(f"\n{'='*70}")
    print("C3(v) FORMULA: c3(v) = s(v)*(n-1-s(v)) - something?")
    print(f"{'='*70}")

    n = 6
    # c3(v) for a vertex v equals: sum_{a in out(v), b in in(v)} I(a->b) - sum_{a in out(v), b in in(v)} I(b->a)
    # Actually: c3(v) = #{(a,b) : v->a, a->b, b->v} = #{directed 2-paths from v back to v}
    # = sum_{a: v->a} |out(a) ∩ in(v) \ {v}|

    # But there's a nice formula: let s = s(v), then
    # The edges among out(v) contribute to NOT forming 3-cycles through v
    # Specifically: c3(v) = s*(n-1-s) - (#{edges from out(v) to out(v)} - #{edges from in(v) to in(v)}) ... hmm

    # Let's just empirically check: c3(v) vs s(v)*(n-1-s(v))
    diffs = defaultdict(list)
    for bits in range(min(total, 10000)):
        A = bits_to_adj(bits, n)
        c3 = count_3cycles_per_vertex(A, n)
        for v in range(n):
            s = sum(A[v])
            product = s*(n-1-s)
            d = product - c3[v]
            diffs[(s, c3[v])].setdefault(d, 0)
            diffs[(s, c3[v])][d] = diffs[(s, c3[v])].get(d, 0) + 1

    # Actually, I know the formula. For vertex v with out-set S:
    # c3(v) = s(n-1-s) - |{directed edges within out(v)}|_{excess}
    # Let me use: c3(v) = s(n-1-s) - sum_{a in out(v)} s_S(a) where s_S(a) = |{b in out(v): a->b}|
    # No, more precisely: c3(v) = sum_{a:v->a} |{b: a->b, b->v, b!=a}| = sum_{a:v->a} |out(a) ∩ in(v)|

    # For each v: Let W = out(v), W_bar = in(v). |W| = s, |W_bar| = n-1-s.
    # c3(v) = sum_{a in W} |out_T(a) ∩ W_bar| = sum_{a in W} (# of W_bar vertices that a beats)
    #        = |{(a,b) : a in W, b in W_bar, a->b}|
    #        = (# of edges from W to W_bar)
    # But total edges between W and W_bar = s*(n-1-s) (tournament: exactly one direction each).
    # Edges W->W_bar: some number m. Edges W_bar->W: s*(n-1-s) - m.
    # c3(v) = m = |{(a,b) : v->a, a->b, b->v}|
    # The complement: s*(n-1-s) - m = |{(a,b) : a in W, b in W_bar, b->a}|
    # = |{(b,a) : b in W_bar, a in W, b->a}| which are "transitive triples with v"
    # i.e. v->a, b->v, b->a => transitive at (b,v,a).
    # So c3(v) = s(n-1-s) - #{transitive triples through v as middle}
    # = s(n-1-s) - sum_{a in W} |{b in W_bar : b->a}|
    # = s(n-1-s) - sum_{a in W} |in_T(a) ∩ W_bar|

    # Anyway, the key point is: c3(v) = m where m depends on the cross-edges between out(v) and in(v).
    # For the MEDIAN score vertex, s(n-1-s) is maximized, giving a "budget" for c3.

    print("Formula: c3(v) = #{edges from out(v) to in(v)}")
    print("Upper bound: c3(v) <= s(v)*(n-1-s(v))")
    print("This is maximized at s = (n-1)/2 (median score)")

    # Can we prove: bad vertex v has c3(v) < c3(w) for some good w?
    # Or equivalently: max c3 vertex is good?

    # Key insight: if v is bad, T\v has a free component spanning all n-2 vertices.
    # The cycles freed by v's deletion must have had v as their only dominator.
    # v dominates cycle {a,b,c} means v->a,b,c or a,b,c->v.
    # So v is either a source or sink relative to {a,b,c}.
    # If v beats all of {a,b,c}, then {a,b,c} ⊆ out(v).
    # Each such cycle contributes to the count of 3-cycles NOT through v
    # (since {a,b,c} are a 3-cycle among themselves, in out(v)).

    # WAIT. Cycles among out(v) vertices are NOT counted in c3(v).
    # c3(v) counts 3-cycles CONTAINING v.
    # Dominated-by-v cycles are in out(v) or in(v), NOT containing v.

    # So removing v doesn't reduce c3 of other vertices much... let me think.

    print(f"\n{'='*70}")
    print("PROVING THE MAX-C3 RULE")
    print(f"{'='*70}")

    # Hypothesis: For a bad vertex v, deleting v "frees" cycles in out(v) or in(v).
    # These freed cycles exist AMONG the out-neighbors or in-neighbors of v.
    # But a vertex w with many c3(w) cycles through it is "entangled" with many
    # different subsets of vertices — it's hard for w to be the unique dominator
    # of a connected component of cycles.

    # Let me count: for bad vertex v, how many cycles does v uniquely dominate?
    n = 6
    v_uniq_dom_counts = []
    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        cycles = find_3cycles(A, n)
        b1_tv = []
        for v_test in range(n):
            B = delete_vertex(A, n, v_test)
            b1_tv.append(compute_b1(B, n-1))

        for v in range(n):
            if b1_tv[v] == 0:
                continue
            # v is bad. Count cycles uniquely dominated by v.
            uniq_count = 0
            for cyc in cycles:
                if v in cyc:
                    continue
                a, b, c = cyc
                if (A[v][a] and A[v][b] and A[v][c]) or (A[a][v] and A[b][v] and A[c][v]):
                    # v dominates cyc. Check uniqueness.
                    other_dom = False
                    for d in range(n):
                        if d == v or d in {a,b,c}:
                            continue
                        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
                            other_dom = True
                            break
                    if not other_dom:
                        uniq_count += 1
            v_uniq_dom_counts.append(uniq_count)

    if v_uniq_dom_counts:
        from collections import Counter
        dist = Counter(v_uniq_dom_counts)
        print(f"\nUniquely-v-dominated cycle count for bad vertices (n=6):")
        for k in sorted(dist.keys()):
            print(f"  {k} cycles: {dist[k]}")

    # Check: cycles in out(v) vs in(v)
    print(f"\nBad vertex: cycles in out(v) vs in(v) subtournaments:")
    out_in_counts = defaultdict(int)
    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue

        c3 = count_3cycles_per_vertex(A, n)
        b1_tv = []
        for v_test in range(n):
            B = delete_vertex(A, n, v_test)
            b1_tv.append(compute_b1(B, n-1))

        for v in range(n):
            if b1_tv[v] == 0:
                continue
            s = sum(A[v])
            out_v = [u for u in range(n) if u != v and A[v][u]]
            in_v = [u for u in range(n) if u != v and A[u][v]]
            # 3-cycles among out(v)
            c3_out = 0
            for a, b, c in combinations(out_v, 3):
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    c3_out += 1
            # 3-cycles among in(v)
            c3_in = 0
            for a, b, c in combinations(in_v, 3):
                if (A[a][b] and A[b][c] and A[c][a]) or (A[a][c] and A[c][b] and A[b][a]):
                    c3_in += 1
            out_in_counts[(s, c3_out, c3_in)] += 1

    print(f"\n  (score, c3_in_out, c3_in_in): count")
    for key in sorted(out_in_counts.keys()):
        print(f"  {key}: {out_in_counts[key]}")

if __name__ == '__main__':
    main()
