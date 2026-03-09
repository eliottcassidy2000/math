"""
Find a vertex selection rule that ALWAYS picks a good vertex.

From earlier analysis:
- max c3(v): 100% at n=5,6 but ~98.9% at n=7,8
- max/min score: 92.6% at n=6
- We need BETTER rules.

New ideas to test:
1. Vertex NOT in any uniquely-dominated cycle (as a member)
2. Vertex with max # of dominated cycles through it
3. Vertex in the most "cross-component" position (connected to most dominated cycles)
4. Vertex with highest "free cycle count" not through it (highest f(v))
5. Vertex whose deletion minimally changes component structure
6. "Anti-bridge": vertex whose cycles DON'T bridge free and dominated parts
"""
import numpy as np
from itertools import combinations, permutations
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

def random_tournament(n, rng):
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if rng.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A

def find_3cycles(A, n):
    cycles = []
    for a, b, c in combinations(range(n), 3):
        if A[a][b] and A[b][c] and A[c][a]:
            cycles.append((a, b, c))
        if A[a][c] and A[c][b] and A[b][a]:
            cycles.append((a, c, b))
    return cycles

def get_dom(A, n, cyc):
    doms = []
    a, b, c = cyc
    for d in range(n):
        if d in {a,b,c}: continue
        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
            doms.append(d)
    return doms

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
        if visited[start]: continue
        comp = []
        stack = [start]
        while stack:
            v = stack.pop()
            if visited[v]: continue
            visited[v] = True
            comp.append(v)
            for u in adj[v]:
                if not visited[u]: stack.append(u)
        components.append(comp)
    return components

def is_dominated(A, n, cyc):
    return len(get_dom(A, n, cyc)) > 0

def compute_b1(A, n):
    cycles = find_3cycles(A, n)
    if not cycles:
        return 0
    comps = cycle_graph_components(cycles)
    return sum(1 for comp in comps if all(not is_dominated(A, n, cycles[ci]) for ci in comp))

def delete_vertex(A, n, v):
    keep = [i for i in range(n) if i != v]
    B = np.zeros((n-1, n-1), dtype=int)
    for i, ki in enumerate(keep):
        for j, kj in enumerate(keep):
            B[i][j] = A[ki][kj]
    return B

def is_sc(A, n):
    for start in range(n):
        visited = set()
        stack = [start]
        while stack:
            v = stack.pop()
            if v in visited: continue
            visited.add(v)
            for u in range(n):
                if A[v][u] and u not in visited:
                    stack.append(u)
        if len(visited) < n:
            return False
    return True

def compute_vertex_features(A, n, cycles):
    """Compute multiple features for each vertex."""
    features = {}
    for v in range(n):
        # c3: number of 3-cycles through v
        c3 = sum(1 for cyc in cycles if v in set(cyc))

        # dom_thru: number of dominated cycles through v
        dom_thru = sum(1 for cyc in cycles if v in set(cyc) and is_dominated(A, n, cyc))

        # free_thru: number of free cycles through v
        free_thru = c3 - dom_thru

        # unique_dom_of: number of cycles uniquely dominated BY v
        unique_dom_of = sum(1 for cyc in cycles
                           if v not in set(cyc) and get_dom(A, n, cyc) == [v])

        # f(v): freed cycle count
        f_v = sum(1 for cyc in cycles
                  if v not in set(cyc) and len(get_dom(A, n, cyc)) <= 1
                  and (len(get_dom(A, n, cyc)) == 0 or get_dom(A, n, cyc) == [v]))

        # score
        score = int(sum(A[v]))

        # dom_above_v: cycles uniquely dominated by v from above (v beats all)
        dom_above = 0
        dom_below = 0
        for cyc in cycles:
            if v in set(cyc):
                continue
            doms = get_dom(A, n, cyc)
            if doms != [v]:
                continue
            a, b, c = cyc
            if A[v][a] and A[v][b] and A[v][c]:
                dom_above += 1
            else:
                dom_below += 1

        # multi_dom_thru: cycles through v that have multiple dominators
        multi_dom_thru = sum(1 for cyc in cycles
                            if v in set(cyc) and len(get_dom(A, n, cyc)) >= 2)

        # Bridge test: does removing v's cycles split the cycle graph?
        restricted = [cyc for cyc in cycles if v not in set(cyc)]
        if restricted:
            r_comps = cycle_graph_components(restricted)
            all_comps = cycle_graph_components(cycles)
            delta_comp = len(r_comps) - len(all_comps)
        else:
            delta_comp = 0

        features[v] = {
            'c3': c3, 'score': score, 'dom_thru': dom_thru,
            'free_thru': free_thru, 'unique_dom_of': unique_dom_of,
            'f': f_v, 'dom_above': dom_above, 'dom_below': dom_below,
            'multi_dom_thru': multi_dom_thru, 'delta_comp': delta_comp,
        }
    return features

def test_rule(rule_name, selector, tournaments, n):
    """Test a selection rule. selector(features, A, n) -> selected vertex."""
    successes = 0
    total = 0
    for A, bad_set in tournaments:
        cycles = find_3cycles(A, n)
        if not cycles:
            continue
        features = compute_vertex_features(A, n, cycles)
        selected = selector(features, A, n)
        total += 1
        if selected not in bad_set:
            successes += 1
    return successes, total

def main():
    rng = np.random.RandomState(42)

    for n in [6, 7, 8]:
        print(f"\n{'='*70}")
        print(f"TESTING SELECTION RULES AT n={n}")
        print(f"{'='*70}")

        # Collect tournaments with bad vertices
        tournaments = []

        if n <= 6:
            total = 2**(n*(n-1)//2)
            sample = range(total)
        else:
            nsamp = {7: 5000, 8: 2000}[n]
            sample = range(nsamp)

        for idx in sample:
            if n <= 6:
                A = bits_to_adj(idx, n)
            else:
                A = random_tournament(n, rng)

            if compute_b1(A, n) != 0:
                continue

            bad = set()
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad.add(v)

            if len(bad) > 0:
                tournaments.append((A, bad))

        print(f"  {len(tournaments)} tournaments with >= 1 bad vertex")

        # Define selection rules
        rules = {
            'max_c3': lambda f, A, n: max(range(n), key=lambda v: f[v]['c3']),
            'min_c3': lambda f, A, n: min(range(n), key=lambda v: f[v]['c3']),
            'max_score': lambda f, A, n: max(range(n), key=lambda v: f[v]['score']),
            'min_score': lambda f, A, n: min(range(n), key=lambda v: f[v]['score']),
            'max_dom_thru': lambda f, A, n: max(range(n), key=lambda v: f[v]['dom_thru']),
            'min_dom_thru': lambda f, A, n: min(range(n), key=lambda v: f[v]['dom_thru']),
            'max_free_thru': lambda f, A, n: max(range(n), key=lambda v: f[v]['free_thru']),
            'min_free_thru': lambda f, A, n: min(range(n), key=lambda v: f[v]['free_thru']),
            'max_unique_dom': lambda f, A, n: max(range(n), key=lambda v: f[v]['unique_dom_of']),
            'min_unique_dom': lambda f, A, n: min(range(n), key=lambda v: f[v]['unique_dom_of']),
            'max_f': lambda f, A, n: max(range(n), key=lambda v: f[v]['f']),
            'min_f': lambda f, A, n: min(range(n), key=lambda v: f[v]['f']),
            'max_multi_dom_thru': lambda f, A, n: max(range(n), key=lambda v: f[v]['multi_dom_thru']),
            'min_delta_comp': lambda f, A, n: min(range(n), key=lambda v: f[v]['delta_comp']),
            'max_dom_above': lambda f, A, n: max(range(n), key=lambda v: f[v]['dom_above']),
            'max_dom_below': lambda f, A, n: max(range(n), key=lambda v: f[v]['dom_below']),
            # Composite rules
            'max_c3_break_score': lambda f, A, n: max(range(n),
                key=lambda v: (f[v]['c3'], f[v]['score'])),
            'max_c3_break_dom': lambda f, A, n: max(range(n),
                key=lambda v: (f[v]['c3'], f[v]['dom_thru'])),
            'min_unique_dom_break_c3': lambda f, A, n: min(range(n),
                key=lambda v: (f[v]['unique_dom_of'], -f[v]['c3'])),
            'min_f_break_c3': lambda f, A, n: min(range(n),
                key=lambda v: (f[v]['f'], -f[v]['c3'])),
            'zero_unique_dom_max_c3': lambda f, A, n: max(
                [v for v in range(n) if f[v]['unique_dom_of'] == 0] or list(range(n)),
                key=lambda v: f[v]['c3']),
            'max_dom_thru_minus_free_thru': lambda f, A, n: max(range(n),
                key=lambda v: f[v]['dom_thru'] - f[v]['free_thru']),
            'min_dom_balance': lambda f, A, n: min(range(n),
                key=lambda v: abs(f[v]['dom_above'] - f[v]['dom_below'])),
        }

        results = {}
        for name, selector in rules.items():
            s, t = test_rule(name, selector, tournaments, n)
            results[name] = (s, t)

        # Sort by success rate
        sorted_rules = sorted(results.items(), key=lambda x: -x[1][0]/max(x[1][1],1))
        for name, (s, t) in sorted_rules:
            pct = 100*s/t if t > 0 else 0
            marker = " ***" if pct == 100 else ""
            print(f"  {name}: {s}/{t} ({pct:.1f}%){marker}")

    # PART 2: Look at failures of the best rules at n=7
    print("\n" + "=" * 70)
    print("PART 2: Analyze failures of best rules at n=7")
    print("-" * 50)

    rng = np.random.RandomState(42)
    n = 7
    nsamp = 10000

    # Track failures for top rules
    failure_analysis = defaultdict(list)

    for _ in range(nsamp):
        A = random_tournament(n, rng)
        if compute_b1(A, n) != 0:
            continue

        bad = set()
        for v in range(n):
            B = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad.add(v)

        if len(bad) == 0:
            continue

        cycles = find_3cycles(A, n)
        features = compute_vertex_features(A, n, cycles)

        # Test max_c3
        selected = max(range(n), key=lambda v: features[v]['c3'])
        if selected in bad:
            failure_analysis['max_c3'].append((A, bad, features, selected))

        # Test max_dom_thru
        selected = max(range(n), key=lambda v: features[v]['dom_thru'])
        if selected in bad:
            failure_analysis['max_dom_thru'].append((A, bad, features, selected))

        # Test min_unique_dom_break_c3
        selected = min(range(n), key=lambda v: (features[v]['unique_dom_of'], -features[v]['c3']))
        if selected in bad:
            failure_analysis['min_unique_dom_break_c3'].append((A, bad, features, selected))

    for rule_name, failures in failure_analysis.items():
        print(f"\n{rule_name}: {len(failures)} failures")
        if failures:
            # Analyze first few
            for i, (A, bad, features, sel) in enumerate(failures[:3]):
                print(f"  Failure {i+1}: selected={sel} (BAD), bad={bad}")
                for v in range(n):
                    f = features[v]
                    marker = " [BAD]" if v in bad else " [GOOD]"
                    marker2 = " <-- SELECTED" if v == sel else ""
                    print(f"    v={v}: c3={f['c3']}, score={f['score']}, "
                          f"dom_thru={f['dom_thru']}, unique_dom={f['unique_dom_of']}, "
                          f"f={f['f']}, delta_comp={f['delta_comp']}"
                          f"{marker}{marker2}")

    # PART 3: Can we combine multiple rules to get 100%?
    # Try: if rule A picks good, use it; else try rule B; else rule C.
    print("\n" + "=" * 70)
    print("PART 3: Cascading rules at n=7,8")
    print("-" * 50)

    for n in [7, 8]:
        rng = np.random.RandomState(42)
        nsamp = {7: 10000, 8: 5000}[n]

        cascade_success = 0
        cascade_total = 0

        for _ in range(nsamp):
            A = random_tournament(n, rng)
            if compute_b1(A, n) != 0:
                continue

            bad = set()
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad.add(v)

            if len(bad) == 0:
                continue

            cascade_total += 1
            cycles = find_3cycles(A, n)
            features = compute_vertex_features(A, n, cycles)

            # Cascade: try max_c3, then min_unique_dom, then max_dom_thru
            selected = max(range(n), key=lambda v: features[v]['c3'])
            if selected not in bad:
                cascade_success += 1
                continue

            selected = min(range(n), key=lambda v: features[v]['unique_dom_of'])
            if selected not in bad:
                cascade_success += 1
                continue

            selected = max(range(n), key=lambda v: features[v]['dom_thru'])
            if selected not in bad:
                cascade_success += 1
                continue

            # Try: min delta_comp
            selected = min(range(n), key=lambda v: features[v]['delta_comp'])
            if selected not in bad:
                cascade_success += 1
                continue

        print(f"  n={n}: cascade {cascade_success}/{cascade_total} ({100*cascade_success/max(cascade_total,1):.2f}%)")

    # PART 4: Try "vertex with max dominated-cycle participation that has
    # min unique domination" — focuses on vertices deeply embedded in dominated structure
    print("\n" + "=" * 70)
    print("PART 4: Advanced composite rules")
    print("-" * 50)

    for n in [7, 8, 9]:
        rng = np.random.RandomState(42)
        nsamp = {7: 10000, 8: 5000, 9: 2000}[n]

        rule_stats = defaultdict(lambda: [0, 0])

        for _ in range(nsamp):
            A = random_tournament(n, rng)
            if compute_b1(A, n) != 0:
                continue

            bad = set()
            for v in range(n):
                B = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad.add(v)

            if len(bad) == 0:
                continue

            cycles = find_3cycles(A, n)
            features = compute_vertex_features(A, n, cycles)

            # Rule: vertex with MINIMUM f(v) (least freed cycles)
            # Rationale: good vertex = one whose deletion doesn't free a spanning set
            sel_min_f = min(range(n), key=lambda v: features[v]['f'])
            rule_stats['min_f'][1] += 1
            if sel_min_f not in bad:
                rule_stats['min_f'][0] += 1

            # Rule: vertex with MAXIMUM unique_dom_of (most unique domination)
            # Rationale: if you uniquely dominate many cycles, removing you frees too many
            # So vertex with MAX unique_dom is most likely bad. Pick vertex with 0.
            # But everyone with 0 might also be bad (bridge type).
            # So: pick vertex with MODERATE unique_dom?

            # Rule: vertex in most dominated cycles (deeply embedded in dominated structure)
            sel_max_dom = max(range(n), key=lambda v: features[v]['dom_thru'])
            rule_stats['max_dom_thru'][1] += 1
            if sel_max_dom not in bad:
                rule_stats['max_dom_thru'][0] += 1

            # Rule: vertex with max (dom_thru - unique_dom_of)
            # High dom_thru = in many dominated cycles = removing doesn't free them (they have other doms)
            # Low unique_dom = doesn't uniquely dominate much
            sel = max(range(n), key=lambda v: features[v]['dom_thru'] - features[v]['unique_dom_of'])
            rule_stats['max_dom-unique'][1] += 1
            if sel not in bad:
                rule_stats['max_dom-unique'][0] += 1

            # Rule: vertex with max multi_dom_thru (in many multi-dominated cycles)
            sel = max(range(n), key=lambda v: features[v]['multi_dom_thru'])
            rule_stats['max_multi_dom'][1] += 1
            if sel not in bad:
                rule_stats['max_multi_dom'][0] += 1

            # Rule: vertex with max (c3 - unique_dom_of)
            sel = max(range(n), key=lambda v: features[v]['c3'] - features[v]['unique_dom_of'])
            rule_stats['max_c3-unique'][1] += 1
            if sel not in bad:
                rule_stats['max_c3-unique'][0] += 1

            # Rule: vertex with min (unique_dom_of + delta_comp)
            sel = min(range(n), key=lambda v: features[v]['unique_dom_of'] + features[v]['delta_comp'])
            rule_stats['min_unique+delta'][1] += 1
            if sel not in bad:
                rule_stats['min_unique+delta'][0] += 1

            # Rule: vertex that participates in the most cycles dominated by MULTIPLE external vertices
            # (Not uniquely dominated - so removing this vertex doesn't free those cycles)
            multi_dom_v = [0]*n
            for cyc in cycles:
                doms = get_dom(A, n, cyc)
                if len(doms) >= 2:
                    for u in cyc:
                        multi_dom_v[u] += 1
            sel = max(range(n), key=lambda v: multi_dom_v[v])
            rule_stats['max_in_multidom'][1] += 1
            if sel not in bad:
                rule_stats['max_in_multidom'][0] += 1

        print(f"\nn={n}:")
        for name, (s, t) in sorted(rule_stats.items(), key=lambda x: -x[1][0]/max(x[1][1],1)):
            pct = 100*s/t if t > 0 else 0
            marker = " ***" if pct == 100 else ""
            print(f"  {name}: {s}/{t} ({pct:.1f}%){marker}")

if __name__ == '__main__':
    main()
