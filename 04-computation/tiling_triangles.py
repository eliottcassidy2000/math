"""
Triangles on the tiling grid and their relationship to OCF.

Each triple of vertices (a,b,c) with a<b<c determines a triangle on the grid.
The arcs between them determine whether they form a 3-cycle (part of Omega_3).

A 3-cycle exists iff the three arcs form a directed cycle.
In tiling terms, this is a specific bit pattern on the triangle's tiles.

Key insight: the grid tile (r,c) encodes arc between position c-1 and position c-1+r
in the base path. A triangle {a,b,c} involves arcs:
  (a,b): gap b-a, tile at (b-a-1, a)
  (b,c): gap c-b, tile at (c-b-1, b)
  (a,c): gap c-a, tile at (c-a-1, a)

The tile positions form a triangle in the grid!

Question: Is there a local rule on the tiling grid that counts 3-cycles
and thus determines alpha_1 = c_3?
"""
from itertools import permutations, combinations
from collections import defaultdict, Counter

def build_data(n):
    arcs = [(i,j) for i in range(n) for j in range(i+2, n)]
    m = len(arcs)
    arc_idx = {a: k for k, a in enumerate(arcs)}

    # Map each arc to grid coordinates
    arc_grid = {}
    for k, (i,j) in enumerate(arcs):
        r = j - i - 1  # gap - 1
        c = i           # start position
        arc_grid[k] = (r, c)

    return arcs, m, arc_idx, arc_grid

def analyze_triangles(n):
    print(f"\n{'='*70}")
    print(f"TRIANGLE ANALYSIS: n = {n}")
    print(f"{'='*70}")

    arcs, m, arc_idx, arc_grid = build_data(n)

    # For each triple (a,b,c) with a<b<c, identify the three arcs
    # and their grid positions
    print(f"\nTriangles on the grid:")
    triangles = []
    for a in range(n):
        for b in range(a+1, n):
            for c in range(b+1, n):
                # Arcs: (a,b), (b,c), (a,c) — but only non-path arcs
                arc_ab = (a,b) if b > a+1 else None  # path if b=a+1
                arc_bc = (b,c) if c > b+1 else None  # path if c=b+1
                arc_ac = (a,c)  # always non-path since c >= a+2

                # Grid positions
                positions = []
                arc_list = []
                if arc_ab is not None:
                    positions.append(arc_grid[arc_idx[arc_ab]])
                    arc_list.append(arc_ab)
                else:
                    positions.append('PATH')
                    arc_list.append((a,b))
                if arc_bc is not None:
                    positions.append(arc_grid[arc_idx[arc_bc]])
                    arc_list.append(arc_bc)
                else:
                    positions.append('PATH')
                    arc_list.append((b,c))
                positions.append(arc_grid[arc_idx[arc_ac]])
                arc_list.append(arc_ac)

                # How many non-path arcs?
                np_arcs = sum(1 for p in positions if p != 'PATH')

                triangles.append({
                    'vertices': (a,b,c),
                    'arcs': arc_list,
                    'positions': positions,
                    'non_path_count': np_arcs,
                })

                if n <= 6:
                    print(f"  ({a},{b},{c}): arcs={arc_list}, grid={positions}, "
                          f"non-path={np_arcs}")

    # For a 3-cycle to exist on (a,b,c), we need a directed cycle.
    # The arcs are a→b, b→c, a→c (or their reverses).
    # 3-cycle a→b→c→a: need a→b(fwd), b→c(fwd), c→a(bwd on (a,c))
    #   Tiling bits: arc_ab=0(fwd), arc_bc=0(fwd), arc_ac=1(bwd)
    #   But path arcs are always forward! So if (a,b) is path, arc_ab=fwd always.
    # 3-cycle a→c→b→a: need a→c(fwd on (a,c)), c→b(bwd on (b,c)), b→a(bwd on (a,b))
    #   Tiling bits: arc_ac=0(fwd), arc_bc=1(bwd), arc_ab=1(bwd)
    #   But if arc_ab or arc_bc is path, it's always fwd, so this direction can fail.

    # Key: the number of non-path arcs determines the "flexibility" of the triangle
    np_counts = Counter(t['non_path_count'] for t in triangles)
    print(f"\nNon-path arc counts per triangle: {dict(sorted(np_counts.items()))}")

    # Now compute for all tilings: for each triangle, does it form a 3-cycle?
    print(f"\n--- 3-cycle patterns per triangle type ---")

    results = []
    for bits in range(1 << m):
        # Determine arc directions
        # For non-path arcs: bit=0 means forward (a→b), bit=1 means backward (b→a)
        # Path arcs: always forward
        dir_map = {}
        for k, (i,j) in enumerate(arcs):
            dir_map[(i,j)] = 'bwd' if (bits & (1 << k)) else 'fwd'

        # Count 3-cycles
        c3 = 0
        for tri in triangles:
            a, b, c = tri['vertices']
            # Get direction of each arc
            # (a,b): forward=a→b, backward=b→a
            d_ab = dir_map.get((a,b), 'fwd')  # default fwd if path
            d_bc = dir_map.get((b,c), 'fwd')
            d_ac = dir_map.get((a,c), 'fwd')

            # Check for 3-cycle: a→b→c→a or a→c→b→a
            if d_ab == 'fwd' and d_bc == 'fwd' and d_ac == 'bwd':
                c3 += 1  # a→b→c→a
            elif d_ab == 'bwd' and d_bc == 'bwd' and d_ac == 'fwd':
                c3 += 1  # a→c→b→a

        results.append({'bits': bits, 'c3': c3})

    # c3 distribution
    c3_dist = Counter(r['c3'] for r in results)
    print(f"\n3-cycle count distribution:")
    for c3 in sorted(c3_dist.keys()):
        print(f"  c3 = {c3}: {c3_dist[c3]} tilings ({c3_dist[c3]/len(results)*100:.1f}%)")

    # For each triangle, what fraction of tilings make it a 3-cycle?
    print(f"\nPer-triangle 3-cycle probability:")
    for tri in triangles[:20]:  # limit output
        a, b, c = tri['vertices']
        np = tri['non_path_count']
        # Count tilings where this triangle is a 3-cycle
        cycle_count = 0
        for bits in range(1 << m):
            d_ab = 'bwd' if ((a,b) in arc_idx and bits & (1 << arc_idx[(a,b)])) else 'fwd'
            d_bc = 'bwd' if ((b,c) in arc_idx and bits & (1 << arc_idx[(b,c)])) else 'fwd'
            d_ac = 'bwd' if bits & (1 << arc_idx[(a,c)]) else 'fwd'
            if (d_ab == 'fwd' and d_bc == 'fwd' and d_ac == 'bwd') or \
               (d_ab == 'bwd' and d_bc == 'bwd' and d_ac == 'fwd'):
                cycle_count += 1
        prob = cycle_count / (1 << m)
        print(f"  ({a},{b},{c}): non-path={np}, P(3-cycle)={prob:.4f}")

    # KEY: What determines H(T) from the tiling pattern?
    # Hypothesis: H depends primarily on the number of 3-cycles (c3)
    # and the number of "conflicting" 3-cycles (non-disjoint ones)
    print(f"\n--- H vs c3 relationship ---")
    # Need to compute H for each tiling
    for bits in range(1 << m):
        adj = [[False]*n for _ in range(n)]
        for i in range(n-1):
            adj[i][i+1] = True
        for k, (i,j) in enumerate(arcs):
            if bits & (1 << k):
                adj[j][i] = True
            else:
                adj[i][j] = True
        h = sum(1 for p in permutations(range(n))
                if all(adj[p[i]][p[i+1]] for i in range(n-1)))
        results[bits]['h'] = h

    # Group by c3
    c3_h = defaultdict(list)
    for r in results:
        c3_h[r['c3']].append(r['h'])

    print(f"\nH statistics by c3:")
    for c3 in sorted(c3_h.keys()):
        vals = c3_h[c3]
        avg = sum(vals)/len(vals)
        mn, mx = min(vals), max(vals)
        print(f"  c3={c3}: avg_H={avg:.1f}, min={mn}, max={mx}, count={len(vals)}")

    # Correlation
    all_c3 = [r['c3'] for r in results]
    all_h = [r['h'] for r in results]
    mean_c3 = sum(all_c3)/len(all_c3)
    mean_h = sum(all_h)/len(all_h)
    cov = sum((c-mean_c3)*(h-mean_h) for c,h in zip(all_c3, all_h))/len(all_c3)
    std_c3 = (sum((c-mean_c3)**2 for c in all_c3)/len(all_c3))**0.5
    std_h = (sum((h-mean_h)**2 for h in all_h)/len(all_h))**0.5
    corr = cov/(std_c3 * std_h) if std_c3 > 0 and std_h > 0 else 0
    print(f"\n  Correlation(c3, H) = {corr:.4f}")
    print(f"  mean_H = {mean_h:.2f} = (n-1)! = {1}")
    import math
    print(f"  Expected mean_H = (n-1)! = {math.factorial(n-1)}")
    print(f"  Actual mean_H = {mean_h:.2f}")

if __name__ == '__main__':
    for n in range(3, 7):
        analyze_triangles(n)
