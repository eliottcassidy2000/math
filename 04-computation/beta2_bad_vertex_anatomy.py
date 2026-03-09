"""
Anatomy of bad vertices when b1(T)=0.

Goal: Understand WHY at most 3 vertices can be bad, and prove
that good vertices always exist for n >= 4.

A vertex v is "bad" for tournament T if b1(T\v) > b1(T) = 0.
Since b1 <= 1 (THM-107), bad means b1(T\v) = 1.

Key question: What is the structure of the "free component" that
appears in T\v for a bad vertex v?
"""

import numpy as np
from itertools import combinations
import sys
from collections import defaultdict

def bits_to_adj(bits, n):
    """Convert bit-encoding to adjacency matrix."""
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
    """Find all directed 3-cycles."""
    cycles = []
    for a, b, c in combinations(range(n), 3):
        # Check all 2 orientations of the 3-cycle on {a,b,c}
        if A[a][b] and A[b][c] and A[c][a]:
            cycles.append((a, b, c))
        if A[a][c] and A[c][b] and A[b][a]:
            cycles.append((a, c, b))
    return cycles

def cycle_verts(cyc):
    return frozenset(cyc)

def shared_directed_edge(c1, c2):
    """Check if two 3-cycles share a directed edge."""
    edges1 = {(c1[0],c1[1]), (c1[1],c1[2]), (c1[2],c1[0])}
    edges2 = {(c2[0],c2[1]), (c2[1],c2[2]), (c2[2],c2[0])}
    return len(edges1 & edges2) > 0

def cycle_graph_components(cycles):
    """Compute connected components of the 3-cycle adjacency graph."""
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
    """Check if 3-cycle (a,b,c) is dominated by some external vertex."""
    verts = set(cyc)
    ext = [d for d in range(n) if d not in verts]
    for d in ext:
        a, b, c = cyc
        if A[d][a] and A[d][b] and A[d][c]:
            return True  # d beats all
        if A[a][d] and A[b][d] and A[c][d]:
            return True  # all beat d
    return False

def get_dominators(A, n, cyc):
    """Get all external vertices that dominate this cycle."""
    verts = set(cyc)
    ext = [d for d in range(n) if d not in verts]
    doms = []
    for d in ext:
        a, b, c = cyc
        if (A[d][a] and A[d][b] and A[d][c]) or (A[a][d] and A[b][d] and A[c][d]):
            doms.append(d)
    return doms

def compute_b1(A, n):
    """Compute b1 = #{free components}."""
    cycles = find_3cycles(A, n)
    if not cycles:
        return 0
    comps = cycle_graph_components(cycles)
    free_count = 0
    for comp in comps:
        all_free = True
        for ci in comp:
            if is_dominated(A, n, cycles[ci]):
                all_free = False
                break
        if all_free:
            free_count += 1
    return free_count

def delete_vertex(A, n, v):
    """Return adjacency matrix with vertex v deleted."""
    keep = [i for i in range(n) if i != v]
    B = np.zeros((n-1, n-1), dtype=int)
    for i, ki in enumerate(keep):
        for j, kj in enumerate(keep):
            B[i][j] = A[ki][kj]
    return B, keep

def analyze_bad_vertex(A, n, v):
    """Detailed analysis of why vertex v is bad."""
    B, keep = delete_vertex(A, n, v)

    # Cycles of T\v (using original vertex labels for clarity)
    cycles_Tv = find_3cycles(B, n-1)
    # Map back to original labels
    cycles_orig = [(keep[c[0]], keep[c[1]], keep[c[2]]) for c in cycles_Tv]

    comps = cycle_graph_components(cycles_Tv)

    # Find the free component
    free_comp_indices = []
    for comp in comps:
        all_free = True
        for ci in comp:
            if is_dominated(B, n-1, cycles_Tv[ci]):
                all_free = False
                break
        if all_free:
            free_comp_indices.append(comp)

    if not free_comp_indices:
        return None  # Not actually bad

    # Analyze cycles in the free component (in original labeling)
    free_comp = free_comp_indices[0]
    free_cycles_orig = [cycles_orig[ci] for ci in free_comp]
    free_cycles_Tv = [cycles_Tv[ci] for ci in free_comp]

    info = {
        'vertex': v,
        'score': sum(A[v]),
        'free_comp_size': len(free_comp),
        'free_cycles': free_cycles_orig,
    }

    # For each free cycle: what was its dominator set in T?
    dom_info = []
    for cyc_orig in free_cycles_orig:
        doms = get_dominators(A, n, cyc_orig)
        dom_info.append({
            'cycle': cyc_orig,
            'dominators_in_T': doms,
            'uniquely_v_dominated': (doms == [v]),
            'free_in_T': (len(doms) == 0),
            'v_orientation': 'beats_all' if (A[v][cyc_orig[0]] and A[v][cyc_orig[1]] and A[v][cyc_orig[2]]) else
                           'loses_all' if (A[cyc_orig[0]][v] and A[cyc_orig[1]][v] and A[cyc_orig[2]][v]) else
                           'mixed'
        })
    info['cycle_dom_info'] = dom_info

    # Check: all cycles in free comp must have Dom ⊆ {v}
    all_v_only = all(set(d['dominators_in_T']) <= {v} for d in dom_info)
    info['all_dom_subset_v'] = all_v_only

    # Count cycles that were free in T (not dominated at all)
    n_free_in_T = sum(1 for d in dom_info if d['free_in_T'])
    n_uniq_v = sum(1 for d in dom_info if d['uniquely_v_dominated'])
    info['n_free_in_T'] = n_free_in_T
    info['n_uniquely_v_dominated'] = n_uniq_v

    # Vertex set of the free component
    verts_in_free = set()
    for cyc in free_cycles_orig:
        verts_in_free.update(cyc)
    info['verts_in_free_comp'] = verts_in_free
    info['n_verts_in_free_comp'] = len(verts_in_free)

    # v's relationship to free comp vertices
    v_out = sum(1 for u in verts_in_free if A[v][u])
    v_in = sum(1 for u in verts_in_free if A[u][v])
    info['v_beats_in_free'] = v_out
    info['v_loses_in_free'] = v_in

    return info

def main():
    print("=" * 70)
    print("BAD VERTEX ANATOMY: When b1(T)=0, which vertices are bad?")
    print("=" * 70)

    for n in [5, 6, 7]:
        print(f"\n{'='*70}")
        print(f"n = {n}")
        print(f"{'='*70}")

        total = 2 ** (n*(n-1)//2)
        if n <= 6:
            sample_range = range(total)
            is_exhaustive = True
        else:
            rng = np.random.RandomState(42)
            sample_range = [rng.randint(0, total) for _ in range(2000)]
            is_exhaustive = False

        bad_count_dist = defaultdict(int)
        bad_vertex_scores = []
        bad_vertex_details = []
        total_b1_0 = 0
        total_examined = 0

        # Track: which vertex properties correlate with being good?
        source_good = 0
        sink_good = 0
        source_exists = 0

        # Track free component structure
        free_comp_verts_dist = defaultdict(int)
        free_comp_size_dist = defaultdict(int)

        # Track domination type
        dom_type_counts = defaultdict(int)

        for bits in sample_range:
            A = bits_to_adj(bits, n)
            b1 = compute_b1(A, n)
            if b1 != 0:
                continue
            total_b1_0 += 1
            total_examined += 1

            bad_verts = []
            for v in range(n):
                B, _ = delete_vertex(A, n, v)
                if compute_b1(B, n-1) > 0:
                    bad_verts.append(v)

            n_bad = len(bad_verts)
            bad_count_dist[n_bad] += 1

            if n_bad > 0:
                scores = [sum(A[i]) for i in range(n)]
                for v in bad_verts:
                    bad_vertex_scores.append(scores[v])

                # Analyze first bad vertex in detail (limit output)
                if len(bad_vertex_details) < 20:
                    for v in bad_verts[:1]:
                        info = analyze_bad_vertex(A, n, v)
                        if info:
                            bad_vertex_details.append(info)

                # Track free component structure for ALL bad vertices
                for v in bad_verts:
                    info = analyze_bad_vertex(A, n, v)
                    if info:
                        free_comp_verts_dist[info['n_verts_in_free_comp']] += 1
                        free_comp_size_dist[info['free_comp_size']] += 1
                        for ci in info['cycle_dom_info']:
                            if ci['uniquely_v_dominated']:
                                dom_type_counts['uniquely_v'] += 1
                            elif ci['free_in_T']:
                                dom_type_counts['free_in_T'] += 1
                            else:
                                dom_type_counts['other_dom'] += 1

                # Check source/sink
                max_score = max(scores)
                min_score = min(scores)
                source = scores.index(max_score)
                sink = scores.index(min_score)
                if source not in bad_verts:
                    source_good += 1
                if sink not in bad_verts:
                    sink_good += 1
                source_exists += 1

        # Report
        method = "exhaustive" if is_exhaustive else "sampled"
        print(f"\nTotal b1=0 tournaments: {total_b1_0} ({method})")
        print(f"\nBad vertex count distribution:")
        for k in sorted(bad_count_dist.keys()):
            pct = 100*bad_count_dist[k]/total_b1_0 if total_b1_0 > 0 else 0
            print(f"  {k} bad vertices: {bad_count_dist[k]} ({pct:.1f}%)")

        if source_exists > 0:
            print(f"\nWhen bad vertices exist ({source_exists} tournaments):")
            print(f"  Source (max score) is good: {source_good}/{source_exists} ({100*source_good/source_exists:.1f}%)")
            print(f"  Sink (min score) is good: {sink_good}/{source_exists} ({100*sink_good/source_exists:.1f}%)")
            print(f"  Source OR sink good: always? Let me check...")

        if bad_vertex_scores:
            print(f"\nBad vertex score distribution:")
            score_counts = defaultdict(int)
            for s in bad_vertex_scores:
                score_counts[s] += 1
            for s in sorted(score_counts.keys()):
                print(f"  score {s}: {score_counts[s]} times")

        print(f"\nFree component vertex count (in T\\v):")
        for k in sorted(free_comp_verts_dist.keys()):
            print(f"  {k} vertices: {free_comp_verts_dist[k]}")

        print(f"\nFree component cycle count (in T\\v):")
        for k in sorted(free_comp_size_dist.keys()):
            print(f"  {k} cycles: {free_comp_size_dist[k]}")

        print(f"\nDomination type of cycles in free component:")
        for t, c in sorted(dom_type_counts.items()):
            print(f"  {t}: {c}")

        if bad_vertex_details:
            print(f"\nDetailed examples (first few bad vertices):")
            for i, info in enumerate(bad_vertex_details[:5]):
                print(f"\n  Example {i+1}: v={info['vertex']}, score={info['score']}")
                print(f"    Free comp: {info['free_comp_size']} cycles on {info['n_verts_in_free_comp']} verts")
                print(f"    Verts in free comp: {info['verts_in_free_comp']}")
                print(f"    v beats {info['v_beats_in_free']}, loses to {info['v_loses_in_free']} in free comp")
                print(f"    all_dom_subset_v: {info['all_dom_subset_v']}")
                for ci in info['cycle_dom_info'][:4]:
                    print(f"      Cycle {ci['cycle']}: doms={ci['dominators_in_T']}, "
                          f"uniq_v={ci['uniquely_v_dominated']}, free_T={ci['free_in_T']}, "
                          f"v_orient={ci['v_orientation']}")

        # Reset for next n
        bad_vertex_details = []

    # Now: deeper analysis at n=6 to understand the EXACT constraint
    print(f"\n{'='*70}")
    print("DEEP ANALYSIS n=6: Bad vertex triples")
    print(f"{'='*70}")

    n = 6
    total = 2 ** (n*(n-1)//2)

    # Track co-occurrence of bad vertices
    bad_triple_types = defaultdict(int)
    bad_scores_pattern = defaultdict(int)

    # Track: can we always find a vertex with specific property that's good?
    median_good = 0
    extremal_good = 0  # max or min score
    max_dominator_good = 0  # vertex that dominates most cycles
    n_with_bad = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        b1 = compute_b1(A, n)
        if b1 != 0:
            continue

        bad_verts = []
        for v in range(n):
            B, _ = delete_vertex(A, n, v)
            if compute_b1(B, n-1) > 0:
                bad_verts.append(v)

        if not bad_verts:
            continue
        n_with_bad += 1

        scores = [sum(A[i]) for i in range(n)]
        bad_scores = tuple(sorted([scores[v] for v in bad_verts]))
        bad_scores_pattern[bad_scores] += 1

        # Check: is extremal (max or min score) always good?
        max_s = max(scores)
        min_s = min(scores)
        extremal_vertices = [i for i in range(n) if scores[i] == max_s or scores[i] == min_s]
        if any(v not in bad_verts for v in extremal_vertices):
            extremal_good += 1

        # Check: vertex dominating most cycles
        dom_counts = []
        for v in range(n):
            count = 0
            cycles = find_3cycles(A, n)
            for cyc in cycles:
                if v in cyc:
                    continue
                doms = get_dominators(A, n, cyc)
                if v in doms:
                    count += 1
            dom_counts.append(count)
        max_dom_v = max(range(n), key=lambda i: dom_counts[i])
        if max_dom_v not in bad_verts:
            max_dominator_good += 1

    print(f"\nTournaments with b1=0 and bad vertices: {n_with_bad}")
    print(f"Extremal vertex (max/min score) good: {extremal_good}/{n_with_bad} ({100*extremal_good/n_with_bad:.1f}%)")
    print(f"Max-dominator vertex good: {max_dominator_good}/{n_with_bad} ({100*max_dominator_good/n_with_bad:.1f}%)")

    print(f"\nBad vertex score patterns:")
    for pattern, count in sorted(bad_scores_pattern.items(), key=lambda x: -x[1])[:10]:
        print(f"  scores {pattern}: {count}")

    # CRITICAL TEST: For each bad vertex v, what is its score relative to free comp?
    print(f"\n{'='*70}")
    print("CRITICAL: Bad vertex's relation to its freed component")
    print(f"{'='*70}")

    n = 6
    # v beats ALL free comp verts? Loses to all? Mixed?
    v_relation_counts = defaultdict(int)
    v_comp_overlap = defaultdict(int)  # how many free comp verts does v beat?

    for bits in range(total):
        A = bits_to_adj(bits, n)
        b1 = compute_b1(A, n)
        if b1 != 0:
            continue

        for v in range(n):
            B, _ = delete_vertex(A, n, v)
            if compute_b1(B, n-1) == 0:
                continue

            info = analyze_bad_vertex(A, n, v)
            if info:
                fv = info['verts_in_free_comp']
                v_out = info['v_beats_in_free']
                v_in = info['v_loses_in_free']
                nfv = len(fv)
                v_relation_counts[(v_out, v_in)] += 1

    print("\n(v_beats, v_loses) in free component:")
    for key, count in sorted(v_relation_counts.items()):
        print(f"  beats={key[0]}, loses={key[1]}: {count}")

    # ULTIMATE TEST: Is there a simple vertex selection rule that always works?
    print(f"\n{'='*70}")
    print("VERTEX SELECTION RULES (n=6 exhaustive)")
    print(f"{'='*70}")

    n = 6
    rules = {
        'max_score': 0,
        'min_score': 0,
        'max_or_min_score': 0,
        'median_score': 0,
        'max_total_domination': 0,
        'in_most_3cycles': 0,
        'source_of_condensation': 0,
    }
    n_b1_0 = 0

    for bits in range(total):
        A = bits_to_adj(bits, n)
        if compute_b1(A, n) != 0:
            continue
        n_b1_0 += 1

        scores = [sum(A[i]) for i in range(n)]
        cycles = find_3cycles(A, n)

        # Count 3-cycles per vertex
        cyc_count = [0]*n
        for cyc in cycles:
            for v in cyc:
                cyc_count[v] += 1

        # Count domination per vertex
        dom_count = [0]*n
        for cyc in cycles:
            for v in range(n):
                if v in cyc:
                    continue
                doms = get_dominators(A, n, cyc)
                if v in doms:
                    dom_count[v] += 1

        # Compute b1(T\v) for each v
        b1_tv = []
        for v in range(n):
            B, _ = delete_vertex(A, n, v)
            b1_tv.append(compute_b1(B, n-1))

        # Test rules
        max_s = max(scores)
        min_s = min(scores)

        # Max score vertex
        v = scores.index(max_s)
        if b1_tv[v] == 0:
            rules['max_score'] += 1

        # Min score vertex
        v = scores.index(min_s)
        if b1_tv[v] == 0:
            rules['min_score'] += 1

        # Max or min
        candidates = [i for i in range(n) if scores[i] == max_s or scores[i] == min_s]
        if any(b1_tv[v] == 0 for v in candidates):
            rules['max_or_min_score'] += 1

        # Median score
        sorted_scores = sorted(range(n), key=lambda i: scores[i])
        median_v = sorted_scores[n//2]
        if b1_tv[median_v] == 0:
            rules['median_score'] += 1

        # Max total domination
        v = max(range(n), key=lambda i: dom_count[i])
        if b1_tv[v] == 0:
            rules['max_total_domination'] += 1

        # In most 3-cycles
        v = max(range(n), key=lambda i: cyc_count[i])
        if b1_tv[v] == 0:
            rules['in_most_3cycles'] += 1

    print(f"\nTotal b1=0: {n_b1_0}")
    for rule, count in rules.items():
        pct = 100*count/n_b1_0 if n_b1_0 > 0 else 0
        print(f"  {rule}: {count}/{n_b1_0} ({pct:.1f}%)")

if __name__ == '__main__':
    main()
